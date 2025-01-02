# taranys: Contributing Guidelines

## Contribution workflow

If you'd like to write or modify some code for taranys, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [taranys issues](https://github.com/BU-ISCIII/taranys/issues) to avoid duplicating work. **If there isn't one already, please create one so that others know you're working on this**.
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [taranys repository](https://github.com/BU-ISCIII/taranys/) to your GitHub account.
3. Make the necessary changes / additions within your forked repository following the [code style guidelines](#code-style-guidelines).
4. Modify the [`CHANGELOG`](../CHANGELOG.md) file according to your changes in the appropiate section ([X.X.Xdev]), you should register your changes regarding:
   1. Added enhancements
   2. New modules
   3. Fixes
   4. Removed stuff
   5. Requirements added or version update
5. Update any documentation as needed.
6. [Submit a Pull Request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) against the `develop` branch and send the url to the #pipelines-dev channel in slack (if you are not in the slack channel just wait fot the PR to be reviewed and rebased).

If you're not used to this workflow with git, you can start with:

- Some [docs in the bu-isciii wiki](https://github.com/BU-ISCIII/BU-ISCIII/wiki/Github--gitflow).
- [some slides](https://docs.google.com/presentation/d/1PruqGxPQVxtNcuEbOd86mylXorgYIU5a/edit?pli=1#slide=id.p1) (in spanish).
- some github generic docs [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests).
- even their [excellent `git` resources](https://try.github.io/).

### taranys repo branches

taranys repo works with a two branching scheme: `main` and `develop`.

- `main`: stable code only for releases.
- `develop`: new code development for the different modules.

You need to submit your PR always against `develop`. Once approbed, this changes must be **`rebased`** so we do not create empty unwanted merges.

## Tests

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

There are typically two types of tests that run:

### Lint tests

We use black and flake8 linting based on PEP8 guidelines for python coding. You can check more information [here](https://github.com/BU-ISCIII/BU-ISCIII/wiki/Python#linting).

### Code tests

Taranys modules are executed using a test dataset.

Anyhow you should always submit locally tested code!!

### New version bumping and release

In order to create a new release you need to follow the next steps:

1. Set the new version according to [semantic versioning](https://semver.org/), in our particular case, changes in the `hotfix` branch will change the PATCH version (third one), and changes in develop will typicaly change the MINOR version, unless the developing team decides otherwise.
2. Create a PR bumping the new version against `develop`. For bumping a new version just change [this line](https://github.com/BU-ISCIII/taranys/blob/09c00c1ddd11f7489de7757841aff506ef4b7e1d/setup.py#L5) with the new version.
3. Once that PR is merged, create via web another PR against `main` (origin `develop`). This PR would need 2 approvals.
4. [Create a new release](https://docs.github.com/en/repositories/releasing-projects-on-github/managing-releases-in-a-repository) copying the appropiate notes from the `CHANGELOG`.
5. Once the release is approved and merged, you're all set!

PRs from one branch to another, like in a release should be **`merged`** not rebased, so we avoid conflicts and the branch merge is correctly visualize in the commits history.

### Code style guidelines

We follow PEP8 conventions as code style guidelines, please check [here](https://github.com/BU-ISCIII/BU-ISCIII/wiki/Python#pep-8-guidelines-read-the-full-pep-8-documentation) for more detail.

When developing new code, we strongly recommend to implement LogSum() functions from log_summary.py instead of the classic python logging in order to keep track of all the warnings and errors that may appear during any of the processes.

## Getting help

For further information/help, please ask on the  `#pipelines-dev` slack channel or write us an email! ([bionformatica@isciii.es](emailto:bioinformatica@isciii.es)).