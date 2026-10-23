#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Block-wise AGBD model comparisons on *precomputed* Pred_RH/Gs_RH inputs.

12-fold leave-one-site-out validation. The held-out block is never used to fit
AGBD-model parameters or select RH percentiles. This is conditional AGBD-model
validation, NOT leakage-free validation of the upstream GC-WCM RH correction:
precomputed Pred_RH uses ALS/UAV reference RH from each footprint.

Equations, initial guesses, bounds, and optimizer limits reproduce the source
scripts. Intentional changes: both RH input sets are tested for fair comparison;
the original stepwise loop range(20,20) is empty, so its intended greedy RH
selection is implemented *within training folds*. By default it tests 2/3/4
selected RHs; use --stepwise-max-features 10 for the intended 2..10 study.
"""

import argparse
import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import OptimizeWarning, curve_fit
from sklearn.metrics import r2_score

SUFFIXES = tuple(range(5, 101, 5))
RH_INPUTS = ('Pred_RH', 'Gs_RH')
MODEL_NAMES = ('layer', 'integral', 'power', 'stepwise')
POWER_RHS = (20, 50, 70, 90, 100)


def score(y, pred):
    y = np.asarray(y, dtype=float)
    pred = np.asarray(pred, dtype=float)
    residuals = pred - y
    return {
        'n': len(y),
        'R2': float(r2_score(y, pred)) if len(y) > 1 else np.nan,
        'r': float(np.corrcoef(y, pred)[0, 1])
             if len(y) > 1 and np.std(y) > 0 and np.std(pred) > 0 else np.nan,
        'RMSE': float(np.sqrt(np.mean(residuals ** 2))),
        'MAE': float(np.mean(np.abs(residuals))),
        'Bias': float(np.mean(residuals)),   # prediction minus ALS_AGB
    }


def prep_layers(rh):
    padded = np.pad(rh, ((0, 0), (1, 0)), constant_values=0.0)
    return np.maximum(np.diff(padded, axis=1), 0)


def layer_function(delta, gamma, *alphas):
    return np.sum(np.asarray(alphas) * (delta ** gamma), axis=1)


def layer_jac(delta, gamma, *alphas):
    power = delta ** gamma
    logs = np.zeros_like(delta)
    np.log(delta, out=logs, where=delta > 0)
    return np.column_stack((np.sum(np.asarray(alphas) * power * logs, axis=1), power))


def fit_layer(train_rh, train_y, test_rh):
    n = len(SUFFIXES)
    delta_train = prep_layers(train_rh)
    popt, _ = curve_fit(
        layer_function, delta_train, train_y,
        p0=[1.5] + [0.1] * n,
        bounds=([0.5] + [0.0] * n, [3.0] + [5.0] * n),
        jac=layer_jac, maxfev=10000,
    )
    pred = layer_function(prep_layers(test_rh), *popt)
    params = {'gamma': float(popt[0])}
    params.update({f'alpha_RH{rh}': float(x)
                   for rh, x in zip(SUFFIXES, popt[1:])})
    return pred, params


def integral_function(rh, alpha, k):
    padded = np.pad(rh, ((0, 0), (1, 0)), constant_values=0.0)
    delta = np.maximum(np.diff(padded, axis=1), 0)
    heights = np.maximum((padded[:, :-1] + padded[:, 1:]) / 2.0, 0.1)
    return alpha * np.sum(delta * (heights ** k), axis=1)


def fit_integral(train_rh, train_y, test_rh):
    popt, _ = curve_fit(integral_function, train_rh, train_y,
                        p0=[0.01, 1.0], maxfev=2000)
    return integral_function(test_rh, *popt), {
        'alpha': float(popt[0]), 'k': float(popt[1]),
    }


def power_function(rh_transposed, alpha, beta, *raw_weights):
    weights = np.asarray(raw_weights, dtype=float)
    # Same normalized-weight model as original source script.
    weights = weights / np.sum(weights)
    h_eff = np.dot(weights, rh_transposed)
    return alpha * (np.maximum(h_eff, 0.1) ** beta)


def fit_power(train_rh, train_y, test_rh):
    idx = [SUFFIXES.index(rh) for rh in POWER_RHS]
    xtrain = train_rh[:, idx].T
    xtest = test_rh[:, idx].T
    n = len(POWER_RHS)
    p0 = [0.05, 2.0] + [1.0] * n
    bounds = ([0, 0.5] + [0] * n, [10, 4.0] + [np.inf] * n)
    popt, _ = curve_fit(power_function, xtrain, train_y, p0=p0,
                        bounds=bounds, maxfev=5000)
    norm_weights = popt[2:] / np.sum(popt[2:])
    result = {'alpha': float(popt[0]), 'beta': float(popt[1])}
    result.update({f'weight_RH{rh}': float(w)
                   for rh, w in zip(POWER_RHS, norm_weights)})
    return power_function(xtest, *popt), result


def sqrt_function(rh_transposed, *params):
    a, b = params[:2]
    coefficients = np.asarray(params[2:])
    # Original: a * (b + sum_i c_i * sqrt(RH_i+100))**2
    inner = b + np.sum(coefficients[:, None] * np.sqrt(rh_transposed + 100), axis=0)
    return a * inner**2


def stepwise_fits(train_rh, train_y, test_rh, max_features):
    """Perform greedy RH selection ON TRAINING DATA ONLY, return stages 2..K.

    Original source intended greedy steps 2..10, but range(20,20) prevented
    any training. Each feature candidate is fitted with the source's unbounded
    curve_fit, p0, and maxfev. Use training R2 to select each new RH.
    """
    chosen = []
    predictions = {}
    # Train/test arrays have RHs in SUFFIXES order.
    if np.any(train_rh < -100) or np.any(test_rh < -100):
        raise ValueError('Square-root model requires all RH values >= -100.')
    for k in range(1, max_features + 1):
        best_r2 = -np.inf
        best = None
        errors = []
        for j in range(len(SUFFIXES)):
            if j in chosen:
                continue
            trial = chosen + [j]
            xtrain = train_rh[:, trial].T
            p0 = [0.01, 1.0] + [0.1] * len(trial)
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore', category=OptimizeWarning)
                    popt, _ = curve_fit(sqrt_function, xtrain, train_y,
                                        p0=p0, maxfev=2000)
                fitted = sqrt_function(xtrain, *popt)
                if not np.isfinite(popt).all() or not np.isfinite(fitted).all():
                    continue
                r2 = float(r2_score(train_y, fitted))
            except (ValueError, RuntimeError, FloatingPointError,
                    OverflowError, TypeError) as exc:
                errors.append(f'RH{SUFFIXES[j]}: {exc}')
                continue
            if r2 > best_r2:
                best_r2 = r2
                best = (j, popt)
        if best is None:
            raise RuntimeError(f'All candidate fits failed during greedy step {k}; '
                               f'first failures: {errors[:3]}')
        j, params = best
        chosen.append(j)
        if k >= 2:
            stage = f'stepwise_{k}RH'
            preds = sqrt_function(test_rh[:, chosen].T, *params)
            predictions[stage] = (
                preds,
                {'selected_RH': ','.join(f'RH{SUFFIXES[t]}' for t in chosen),
                 'train_selection_R2': best_r2,
                 'a': float(params[0]), 'b': float(params[1]),
                 **{f'c_RH{SUFFIXES[t]}': float(coef)
                    for t, coef in zip(chosen, params[2:])}}
            )
    return predictions


def read_data(path, group_col, target_col, expected_blocks):
    raw = pd.read_csv(path, low_memory=False)
    rh_cols = [f'{prefix}{rh}' for prefix in RH_INPUTS for rh in SUFFIXES]
    required = [group_col, target_col] + rh_cols
    absent = [col for col in required if col not in raw.columns]
    if absent:
        raise ValueError('Missing columns: ' + ', '.join(absent))

    numeric = raw[[target_col] + rh_cols].apply(pd.to_numeric, errors='coerce')
    good_num = np.isfinite(numeric.to_numpy(dtype=float)).all(axis=1)
    good_group = raw[group_col].notna() & raw[group_col].astype(str).str.strip().ne('')
    valid = good_num & good_group.to_numpy()
    n_missing = int((~valid).sum())
    if n_missing:
        print(f'Dropped {n_missing}/{len(raw)} rows with missing/nonfinite data or group IDs.')
    df = raw.loc[valid].copy()
    df[[target_col] + rh_cols] = numeric.loc[valid]
    df.insert(0, '_original_csv_row', np.flatnonzero(valid) + 2)  # header row=1
    df[group_col] = df[group_col].astype(str).str.strip()
    df = df.reset_index(drop=True)
    groups = sorted(df[group_col].unique())
    if len(groups) != expected_blocks:
        raise ValueError(f'Expected {expected_blocks} blocks but found {len(groups)}: {groups}')
    if (df[group_col].value_counts() < 2).any():
        raise ValueError('Every block must contain at least two valid observations.')
    matrices = {prefix: df[[f'{prefix}{rh}' for rh in SUFFIXES]].to_numpy(dtype=float)
                for prefix in RH_INPUTS}
    return df, matrices, groups


def run(args):
    models = tuple(s.strip() for s in args.models.split(',') if s.strip())
    invalid = sorted(set(models) - set(MODEL_NAMES))
    if not models or invalid:
        raise ValueError(f'Unknown models {invalid}; valid: {MODEL_NAMES}')
    if args.stepwise_max_features < 2 or args.stepwise_max_features > 20:
        raise ValueError('--stepwise-max-features must be between 2 and 20.')
    df, matrices, blocks = read_data(args.input, args.block_col,
                                    args.target_col, args.expected_blocks)
    ref = df[args.target_col].to_numpy(dtype=float)
    group_ids = df[args.block_col].to_numpy(dtype=str)
    pred_cols = {}
    params_rows = []
    block_rows = []
    print(f'LOBO: {len(df):,} observations, {len(blocks)} blocks; '
          f'models={models}, inputs={RH_INPUTS}', flush=True)
    if 'stepwise' in models:
        print('Nested greedy selection can take much longer than the other models. '
              f'Stages 2..{args.stepwise_max_features} will be evaluated.', flush=True)

    for fold, held_out in enumerate(blocks, 1):
        test = group_ids == held_out
        train = ~test
        ytrain = ref[train]
        print(f'[{fold:02d}/{len(blocks)}] held-out={held_out}: '
              f'train={train.sum():,}, test={test.sum():,}', flush=True)
        if train.sum() < 22:
            raise ValueError(f'Insufficient training rows for 21-parameter model: {held_out}')
        for input_type, rh in matrices.items():
            xtrain, xtest = rh[train], rh[test]
            to_run = {}
            for model in models:
                try:
                    if model == 'layer':
                        to_run[model] = fit_layer(xtrain, ytrain, xtest)
                    elif model == 'integral':
                        to_run[model] = fit_integral(xtrain, ytrain, xtest)
                    elif model == 'power':
                        to_run[model] = fit_power(xtrain, ytrain, xtest)
                    elif model == 'stepwise':
                        to_run.update(stepwise_fits(xtrain, ytrain, xtest,
                                                    args.stepwise_max_features))
                except Exception as exc:
                    raise RuntimeError(f'Fit failed: held_out={held_out}, input={input_type}, '
                                       f'model={model}. Original error: {exc}') from exc
            for model, (test_pred, parameters) in to_run.items():
                test_pred = np.asarray(test_pred, dtype=float)
                if test_pred.shape != (int(test.sum()),) or not np.isfinite(test_pred).all():
                    raise ValueError(f'Nonfinite or wrong-sized predictions: {model}, '
                                     f'{input_type}, {held_out}')
                key = (model, input_type)
                if key not in pred_cols:
                    pred_cols[key] = np.full(len(df), np.nan)
                pred_cols[key][test] = test_pred
                block_rows.append({
                    'model': model, 'RH_input': input_type, 'held_out_site': held_out,
                    'n_train': int(train.sum()), **score(ref[test], test_pred),
                })
                params_rows.append({
                    'model': model, 'RH_input': input_type, 'held_out_site': held_out,
                    'n_train': int(train.sum()), 'parameters_json': json.dumps(parameters),
                    **({'selected_RH': parameters['selected_RH']}
                       if 'selected_RH' in parameters else {}),
                })

    pooled_rows = []
    for (model, input_type), predictions in pred_cols.items():
        if not np.isfinite(predictions).all():
            raise RuntimeError(f'Incomplete LOBO predictions: {model}, {input_type}')
        name = f'LOBO_{model}_{input_type}'
        if name in df:
            raise ValueError(f'Refusing to overwrite original column: {name}')
        df[name] = predictions
        pooled_rows.append({'model': model, 'RH_input': input_type,
                            **score(ref, predictions)})

    pooled = pd.DataFrame(pooled_rows).sort_values(['RMSE', 'model', 'RH_input'])
    out = Path(args.outdir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / 'LOBO_AGBD_all_predictions.csv', index=False)
    pooled.to_csv(out / 'LOBO_AGBD_all_metrics_pooled.csv', index=False)
    pd.DataFrame(block_rows).to_csv(out / 'LOBO_AGBD_all_metrics_by_block.csv', index=False)
    pd.DataFrame(params_rows).to_csv(out / 'LOBO_AGBD_all_parameters.csv', index=False)
    print('\nPooled scores (all held-out predictions; lower RMSE is better):')
    print(pooled[['model', 'RH_input', 'n', 'R2', 'r', 'RMSE', 'MAE', 'Bias']]
          .to_string(index=False, float_format=lambda v: f'{v:.4f}'))
    print(f'\nSaved 4 files in: {out.resolve()}')
    print('CAUTION: Corrected Pred_RH already uses reference ALS/UAV RH from each '
          'footprint. These are conditional AGBD-model LOBO scores, not '
          'independent end-to-end GC-WCM pipeline scores.')


def main():
    parser = argparse.ArgumentParser(description='LOBO validation: four AGBD model families on Pred_RH and Gs_RH.')
    parser.add_argument('--input', default='data_pred_all_data_v2.csv')
    parser.add_argument('--outdir', default='LOBO_AGBD_comparison')
    parser.add_argument('--block-col', default='site')
    parser.add_argument('--target-col', default='ALS_AGB')
    parser.add_argument('--expected-blocks', type=int, default=12)
    parser.add_argument('--models', default='layer,integral,power,stepwise',
                        help='Comma-separated: layer,integral,power,stepwise')
    parser.add_argument('--stepwise-max-features', type=int, default=4,
                        help='Nested feature-selection steps. Default: 4; use 10 for RH2..RH10.')
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
