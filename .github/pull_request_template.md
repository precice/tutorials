TODO: Summarize and motivate the changes, link to issues, remove the checklist entries that are not relevant.

## Checklist

- [ ] I added a summary of any user-facing changes (compared to the last release) in the `changelog-entries/<PRnumber>.md`.
- [ ] If I changed `requirements.txt` files, I regenerated sibling `requirements-reference.txt` files with `python3 tools/releasing/update-requirements-reference.py` (pass a path to update only that directory, or `--all` to refresh everything).

For new tutorials or tutorial cases:

- [ ] I followed the [tutorial folder structure](https://precice.org/community-contribute-to-precice.html#contributing-tutorials)
- [ ] I added/updated the tutorial `README.md`
- [ ] I added/updated the tutorial `metadata.yaml`
- [ ] I added tests in `tools/tests/tests.yaml`
- [ ] I submitted a pull request to the website with:
   - An entry in [`_config.yaml`](https://github.com/precice/precice.github.io/blob/master/_config.yml)
   - A [sidebar entry](https://github.com/precice/precice.github.io/blob/master/_data/sidebars/tutorial_sidebar.yml)
   - An entry in the [tutorials overview](https://github.com/precice/precice.github.io/blob/master/content/tutorials/tutorials.md)
