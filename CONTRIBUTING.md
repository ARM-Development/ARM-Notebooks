# Contributor's Guide

## Overview

Welcome! This is the main contributors guide for ARM-Notebooks. All
of the content for ARM-Notebooks is hosted on GitHub in a central repository.
From this document you can learn about the many ways that you can contribute to this community project, and where you can find additional information.

## The many ways to contribute

### Starting a conversation

One of the easiest ways to contribute is to communicate with the
ARM Open Science community. You might, for example, suggest changes
to content on the site, ask a clarifying question, or even respond
to another community member’s question. There are a number of ways
to have a dialogue with the ARM community, such as by opening a
GitHub issue, but the easiest method is to simply create (or respond
to) a post on our [discussion
forum](https://discourse.arm.gov/).

### Add a new Jupyter Notebook to ARM-Notebooks

[ARM-Notebooks](https://arm-development.github.io/ARM-Notebooks)
is a collection of material that the ARM Open Science team believes is essential
knowledge for the ARM Community to access and analyze datasets in an open manner. Adding new, or changing
existing ARM-Notebooks content requires contributors to be familiar
with a few technologies and workflows, described in Advanced
Contributions below.

If you'd like to contribute a Jupyter Notebook to these materials, please reference our [template](Templates/notebook-template) viewable on the next page. This template is available to you in `Templates/notebook-template.ipynb` if you've cloned the [source repository](https://github.com/ARM-Development/ARM-Notebooks), or available as a download [directly from GitHub](https://github.com/ARM-Development/ARM-Notebooks/raw/main/Templates/notebook-template.ipynb).

#### Rendering Your Notebook Locally

##### Create a conda environment

The first time you check out this repository, run:

```bash
conda env update -f environment.yml
```

This will create or update the dev environment (`arm-tutorial-dev`).

##### Make sure your file is in the table contents
In order for your new notebook to appear in the book, you will need to make sure it is listed in the table of contents file. This file is labeled "`_toc.yml`".

For example, if we had a new value added product (VAP) notebook, titled `new-vap.ipynb`, in the `VAPs/` directory, we would add a line to the `_toc.yml` file with this path, tagged as `# THIS IS OUR NEW FILE`.

```yaml
format: jb-book
root: README.md
parts:
  - caption: Getting Started
    chapters:
    - file: getting-started.md
    - file: Templates/notebook-template
    - file: Templates/tutorial-schedule
  - caption: Value-added Products
    chapters:
    - file: VAPs/README.md
      sections:
        - file: VAPs/squire/intro-to-squire.ipynb
        - file: VAPs/new-vap/new-vap.ipynb # THIS IS OUR NEW FILE
```

##### Building the book locally

To build the book locally, run the following:

```bash
conda activate arm-tutorial-dev
jupyter-book build .
```

Finally, you can view the book by opening the file `_build/html/index.html` with your favorite web browser. On most platforms you can simply run:

```bash
open _build/html/index.html
```

### Advanced Contributions

Some contributions, such as adding a new Jupyter Notebook, or making
changes to an existing one, require operating directly on the GitHub
repository that maintains the desired source material. In the
parlance of GitHub, contributors are required to “fork a repo”, and
“submit a pull request”. This process, though widely used by open
development projects, is unfortunately complex and somewhat varied
from one project to another. Below, we describe in great detail the
process of setting up GitHub/Git, installing and configuring conda,
and submitting a PR.

## Ready, set, fork!

As noted above contributing directly to a GitHub repository that
is not owned by you is a somewhat complicated process that involves
a number of technologies. In the paragraphs below we describe the
process of "forking" an external repository, making changes, and
submitting your changes back to the owners of the external repository.
But first we discuss how to configure your Python environment with the
conda (anaconda) package manager, and how to configure Git and GitHub. These
are one-time steps that you should not need to perform again.

It's a long journey, but the steps we describe below are common to many, many
open source projects hosted on GitHub.

### Getting started with GitHub and Git

Contributing to the ARM-Notebooks repo requires using GitHub, as
already mentioned, and also Git. The latter, Git, is an open source,
command line tool for collaborative software version control, while
GitHub is an online, web-accessible service that greatly simplifies
using the powerful, yet often complex, Git.

```{Note}
GitHub operates entirely within a web browser. You do not
need to install anything, but you will need to set up a free GitHub
account. Git is a command line tool that is most likely already
installed on your machine, and will need to be run from a “terminal” window, AKA
a “shell”.
```

Using, and even just configuring, Git and GitHub are often the most
daunting aspects of contributing to a GitHub hosted project. Here
are the basic steps for Git/GitHub configuration, all of which must
be performed before the next subsection, forking a repo.

#### GitHub Setup

Create a free GitHub account [here](https://github.com/). Note
GitHub offers free personal use and paid enterprise accounts. The
free account is all that is needed.

#### Git Setup

If not already installed on your machine, download and install the
[latest version of Git](https://git-scm.com/downloads) . Set up
Git with a user name and your email using the steps below. Note,
it is advisable that you use the same user name/email as you did
when setting up your GitHub account, though technically this may
not be necessary. Once git is installed you will need to open a
terminal/shell and type the following commands to configure git:

```
$ git config --global user.name "Your name here"
$ git config --global user.email "your_email@example.com"
```

Don’t type the \$. This simply indicates the command line prompt.

#### Configure your environment to authenticate with GitHub from Git

This is a complicated process and there are two authentication
protocols supported: HTTP or SSH. Either will work fine, but we
find HTTP to be the easiest to set up. Both processes are described
in detail on the GitHub site [here](https://docs.github.com/en/get-started/quickstart/set-up-git#authenticating-with-github-from-git).

For further reading see the [GitHub Getting Started
Guide](https://docs.github.com/en/github/getting-started-with-github).

### Creating a Python environment that will work with ARM-Notebooks

Before starting any development, you’ll need to create an isolated
Python environment that will work for ARM-Notebooks. When we use the term
“Python environment” here we are referring to the Python programming
language plus the myriad packages that go along with it. Because
there are so many Python packages available, maintaining interoperability
between them is a huge challenge. To overcome some of these
difficulties the use of Anaconda or miniconda is required to manage
your Python ecosystem. These package managers allow you to create
a separate, custom Python environment for each specific Python set
of tools. Yes, this unfortunately results in multiple copies of
Python on your system, but it greatly reduces breaking toolchains
whenever a change is made to your Python environment (and is more
reliable than any other solution we’ve encountered). Also, the use
of Anaconda/Miniconda is standard practice for working with Python
in general, not simply for using ARM-Notebooks.

The steps:

1. Install either Anaconda or miniconda
1. Make sure your conda is up to date by running this command from the terminal:

```
conda update conda
```

At this point you have a current version of conda available on your
desktop or laptop. Before using your conda environment to work on
ARM-Notebooks content, you'll need to perform an addtional one-time setup
that is specific to each ARM-Notebooks repo. After the one-time configuration is
complete you will need to "activate" a repo-specific environment whenever
you wish to use it.

More information on installing and using conda may be found
[here](https://foundations.projectARM-Notebooks.org/foundations/conda.html).

### Forking a repo

With Git, GitHub, and conda properly installed and configured, we are ready to
"fork" a repository so
that you can safely make changes to the repository contents without
changing the ARM-Notebooks repos until you are ready. Forking a repository
as described below creates a clone of the ARM-Notebooks repo under your
own account on GitHub. Any changes you make to your repository will
only be seen by you until you are ready to share them with others,
and hopefully “merge” your changes into one of the official ARM-Notebooks
repositories.

Note that the ARM-Notebooks maintainers employ a version of the “Forking Workflow”
to support contributions from the outside world. This workflow is
summarized below, and described in detail
[here](https://www.atlassian.com/git/tutorials/comparing-workflows/forking-workflow).

The steps:

1. Navigate your web browser to the ARM-Notebooks repository. For example, [here](https://github.com/ARM-Development/ARM-Notebooks).
1. Click on the “Fork” icon (upper right). This will create a copy of the ARM-Notebooks repository
   on the GitHub server under **your** account name. You may be prompted
   to sign in. If so, use the GitHub (not Git) account name and password
   that you created when you created your GitHub account above.
1. After successfully forking a ARM-Notebooks repo you should have a copy
   of that repository on the GitHub server under your account name.
   To verify this you can navigate to the GitHub [home page](https://github.com/), sign in if you are not
   already, and click on “your repositories” under the pull down menu
   in the top right corner of the page. You should see the ARM-Notebooks
   repository you just cloned listed. Click on it. This is a remote
   clone of the ARM-Notebooks repo. Changes you make to your copy will not
   impact the contents of the ARM-Notebooks repo. The next step is to make a
   local copy (clone) of the just-cloned GitHub repository on your
   laptop or workstation.

That’s right. After this final step you will now have **two** copies
of the repo, one local and one remote. From a terminal window type:

```
$ git clone https://github.com/YOUR_USER_NAME/ARM-Notebooks.git
$ cd ARM-Notebooks
```

where `YOUR_USER_NAME` is your GitHub user name.

The git command above does two things, and it is important to
understand them. Firstly, it creates a local repository inside of
the hidden directory `ARM-Notebooks/.git`, that is populated with
the contents of the repository you are cloning. You should never
need to operate directly on the contents of the .git directory.
Secondly, it creates a copy of the repo’s assets (e.g. Python files,
documentation, etc.) in the local directory `ARM-Notebooks`. The files in this latter directory are the ones that you will edit.

Remember: you now essentially have two clones of the `ARM-Notebooks` repository, one
on the GitHub server under your account, and one on your local
workstation or laptop.

Next, connect your local copy of the repository to the “upstream” (remote) ARM-Notebooks repository:

```
$ git remote add upstream https://github.com/ARM-Development/ARM-Notebooks.git
```

Finally, create a new branch in your local repository:

```
$ git checkout -b YOUR_BRANCH_NAME
```

Where `YOUR_BRANCH_NAME` is the name that you want to give your local
branch. What name should you choose? If the work that you are doing
is associated with a GitHub issue you should follow the convention:

_issue_XXX_

Where _XXX_ is the GitHub issue number. If it is not associated with a ARM-Notebooks GitHub issue, pick something short and meaningful, e.g. “documentation_cleanup”.

You can now make changes to your local copy of the ARM-Notebooks repo
without having those changes affect either the remote ARM-Notebooks GitHub
repo, your remote, personal GitHub repo, or your local repo in .git,
until you are ready to merge your local changes upstream, first to
your .git local repo, then to your remote GitHub repo, and then
ultimately to the remote ARM-Notebooks GitHub repo. Simple, right? More
on this later.

For further information see the [GitHub docs on forking a repo](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo).

```{Note}
At this point you should have a local copy of the repository
in your current working directory. You can safely make changes to
any of the contents. Once you are ready to contribute your changes
back to the ARM-Notebooks repository you will need to submit a Pull Request
(PR), described later.
```

### Make your local changes

At this point you should be able to make your desired changes to
the local copies of your files. **Always consult the repo-specific contributor’s
guide for information specific to the repo you are working on.**

### Submit a Pull Request (PR)

Once you have completed making changes to your local copy of a
ARM-Notebooks repository and are ready to have your changes merged with a
ARM-Notebooks repository on GitHub, you need to essentially perform the
reverse process used to acquire a copy of the ARM-Notebooks repo, and
submit a PR asking the ARM-Notebooks maintainers to consider your merge
request. The merge will occur between your personal GitHub repository
and the ARM-Notebooks GitHub repository, so you first need to merge any
changes you’ve made in your local copy into your local .git repo.
Next, you need to merge these local changes with your personal
remote repo on GitHub. Finally, you need to submit a request to
merge your personal GitHub repo with the ARM-Notebooks GitHub repo.

Git has lots and lots of commands, each with lots and lots of
options. Here we only cover the very basics. Detailed information
about Git can be found [here](https://git-scm.com/), but your best
friend for figuring out to do things with Git may be Google, and
in particular [StackOverflow](https://stackoverflow.com/).

#### Committing Your Code Locally with Git

Changes you’ve made to your local copy of a repository must be
“committed” (merged) to your local repository (the .git subdirectory)
using Git. You can see any uncommitted changes you’ve made to your
local copy of the repository by running the following command from
anywhere (any directory) within the directory where you ran git
checkout:

```
$ git status
```

If you have added any new files you will need to explicitly add
them to the local repo with:

```
$ git add PATH_TO_NEW_FILE
```

Where `PATH_TO_NEW_FILE` is the path name of the newly created file.

To commit changed files, including new files just added with the above command, run the following command from the root of your local copy:

```
$ git commit PATH_TO_NEW_FILE
```

Which will prompt you for a log message. Please provide something informative. If you make lots of changes, it is best to make multiple commits, broken up into related chunks. E.g. “fixed x”, “added documentation”, “added testing”.

```{Note}
When executing `git commit` after `git add PATH_TO_NEW_FILE`,
specifying the path to the new file isn't stricly necessary. However,
in other instances the file path argument is required. We include it
here to keep things simple.
```

#### Pushing Your Changes to Your Personal GitHub Repository

Once all of your changes have been committed to your local .git
repository you are ready to “push” (merge) them with your personal
GitHub repository. To push your .git repository run the following
command from anywhere within your local copy of the repo:

```
$ git push origin FEATURE_NAME
```

Where `FEATURE_NAME` is the name you gave your branch when you checked
it out before starting to make your changes. Typically, if you are
submitting a PR for a change that addresses an open ARM-Notebooks issue,
the name should be _issue_XXX_ where _XXX_ is the issue number.

After successfully running this command your changes will now be
on GitHub under your personal account, but they are not yet part
of the ARM-Notebooks repo. For that to happen one more step is required:
you must **Make the PR**.

But first:

#### Review your Code

Before you make the actual PR, it is a good idea to review the
changes that you’ve made and to have followed all guidelines in
this document, and any repo-specific guidelines.

To review your changes against the official ARM-Notebooks repository do the following:

1. Navigate your web browser to your GitHub repository. E.g.
   https://github.com/YOUR_USER_NAME/ARM-Notebooks
1. Click on `Compare`
1. Check the `head repository` and `compare` branches are set correctly.
   These should be `YOUR_NAME/ARM-Notebooks_REPO_NAME`, and `BRANCH_NAME`,
   respectively, where `BRANCH_NAME` is the name you gave your branch
   when you pushed your changes to your remote repository on GitHub.

Select the “base repository” and “base”. For “base repository” this
should be the ARM-Notebooks repository, for example `ARM-Development/ARM-Notebooks`.
For “base” this should be the branch on the ARM-Notebooks repository that
you wish to compare against (and subsequently merge with).

At this point you should be able to review changes between your
repositories and the GitHub repository.

#### Make the PR

At long last you are ready to make the actual PR, requesting the
ARM-Notebooks community to review your code, make possible suggestions for
changes, and ultimately merge your repo with ARM-Notebooks. To submit a
pull request:

Navigate your web browser to your GitHub repository. E.g. https://github.com/YOUR_USER_NAME/ARM-Notebooks

Click on the "Pull Request" button

Write a description of your changes in the “Preview Discussion” tab. Give an overview of what this PR does, and be sure to indicate any ARM-Notebooks issues that this PR addresses by number.

Click "Send Pull Request"

This request then goes to the repository maintainers, and they will
review the code. If you need to make more changes, you can make
them in your branch, add them to a new commit, push them to GitHub,
and the pull request will be automatically updated. Pushing them
to GitHub again is done with the command:

```
$ git push origin FEATURE_NAME
```

Congratulations!!! You've submitted a PR!