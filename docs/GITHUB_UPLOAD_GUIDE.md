# Upload this project to GitHub

## Prepare the final content

1. Unzip the package. The repository root is the folder containing `README.md`, `run_all.py`, `data/`, and `src/`.
2. Read `docs/RELEASE_STATUS.md` and resolve the missing-data and code/manuscript differences before describing the release as complete.
3. Confirm the intended code/data licenses and finalize `CITATION.cff`.
4. Rerun `python run_all.py` from the documented environment and inspect `results/verification.json` plus the manuscript comparison files.

## Browser upload

1. Sign in to GitHub and create a new repository. Suggested name: `gedi-gc-wcm-vlim`.
2. For an initial review, a private repository can be used. The final materials linked in the manuscript need the public access required by the editor.
3. Choose **Add file → Upload files** and upload the extracted project contents. Keep `README.md` at the repository root. Do not upload only the ZIP: readers should be able to browse the numbered scripts directly.
4. Keep the existing folder structure. Upload in multiple batches if necessary. GitHub documents a 25 MiB per-file browser limit and up to 100 files per upload; see [Adding a file to a repository](https://docs.github.com/en/repositories/working-with-files/managing-files/adding-a-file-to-a-repository). Larger raw datasets can be deposited separately and linked with exact retrieval instructions.
5. Commit the files, open `README.md`, and verify that its relative links work.

## Optional Git command-line upload

Create an empty GitHub repository, open a terminal in the extracted project folder, replace `YOUR_USERNAME` with the real account, and run:

```bash
git init
git add .
git commit -m "Add manuscript code and processed analysis data"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/gedi-gc-wcm-vlim.git
git push -u origin main
```

The public repository URL is a placeholder until the repository has actually been created. No remote repository has been created by this package.

## Freeze the manuscript version

Once the materials and manuscript agree, create a release (for example `v1.0.0`) and archive that exact version. GitHub documents how Zenodo can archive a repository release and issue a DOI: [Referencing and citing content](https://docs.github.com/en/repositories/archiving-a-github-repository/referencing-and-citing-content). Add the real repository URL and DOI to the README, citation metadata, and manuscript.

## Data Availability wording

Do not use a claim that all data/code are available until the missing inputs and release gaps have been resolved. At that point, adapt this template using only the materials actually public:

> The data and code used in this study are publicly available at [actual repository URL] and [actual archived-data DOI]. The repository provides the processing workflow, model implementations, environment specifications, analysis inputs, and instructions for reproducing the reported results. The archived release [actual version/DOI] corresponds to the version used in the revised manuscript.

If any data remain unavailable or restricted, state the precise limitation and resolve its acceptability with the editor rather than retaining the complete-availability wording.
