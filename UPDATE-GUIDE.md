# Install this guide update

Extract this ZIP into the repository root with overwrite, or extract to a temporary folder for review and copy its contents into the repo. Paths start at the repo root, with no enclosing folder.

This update changes only guide content, book configuration and the report template. It includes no datasets, analysis scripts, results, or existing CSS. Keep those files in the repo.

## What changed

- Added `00-setup.qmd` and linked it from the index and shared commands.
- Kept the path synthetic practice, DietSwap, Mouse59, then your own study.
- Replaced the old root Mouse59 preparation README with a repository overview.
- Removed outdated dataset command references and corrected Mouse59 folders.
- Kept each analysis command on the shared workflow page.
- Marked independent sample calculations optional while retaining table checks.
- Added short standard Quarto notes and a plain-language glossary.
- Added a full-size Mouse59 composition link and improved the report template.
- Made code visible by default and removed the old fixed book date.

## After extraction

Run `quarto render` from the repository root. If a supplied dataset's output is missing, generate it using `01-run-the-workflow.qmd` first. Review the new setup page, Mouse59 command anchor, and figures in the rendered book.

This content update does not remove unrelated old files from your repository or deployment. Remove stale published pages separately if they remain on the server.

Validation in this environment checks YAML, chapter paths, relative QMD links and heading anchors, code fences, and obsolete references. R analysis and Quarto rendering could not be run here because those programs are unavailable. Figure references depend on your existing result files.
