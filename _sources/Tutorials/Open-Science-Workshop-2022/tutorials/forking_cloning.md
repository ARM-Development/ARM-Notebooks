# Forking a Repository

## Overview:

1. Forking a repository

---

## What's Forking?

Forking is essentially _making a copy_ of a Github Repo! Scenarios where forking a repo is indicated include the following:

1. You wish to collaborate on projects that are hosted on GitHub, but you are not one of that project's _maintainers_ (i.e., you do not have _write permissions_ on it).
1. You wish to experiment with changing or adding new features to a project, and do not immediately intend to _merge_ them into the original project's repo (aka, the _base_ repository).

In a fork, you first create a copy of an existing repository, but store it in your own personal GitHub organization (recall that when you create a GitHub account, the _organization_ name is your GitHub user ID).

Let's say we intend to make some changes to the [Aerosol-Mentors-Processing-Routines Repo](https://github.com/ARM-Development/Aerosol-Mentors-Processing-Routines) repo, that ultimately we'll submit to the original repository as a _Pull request_.

```{note}
Be sure you have logged into GitHub at this time!
```

Notice at the top right of the screen, there is a _Fork_ button:

<img src="https://github.com/ProjectPythiaTutorials/github-arm-2022-02-28/blob/main/tutorial/images/Github_Fork.jpeg?raw=true" alt="Fork">

---

Click on it! You should see your GitHub user id (if you administer any other GitHub organizations, you will see them as well). Click on your user id to complete the _fork_. 

<img src="https://github.com/ProjectPythiaTutorials/github-arm-2022-02-28/blob/main/tutorial/images/Github_ForkDest.jpeg?raw=true" alt="ForkDest">

---

After a few seconds, your browser will be redirected to the forked repo, now residing in your personal GitHub organization!

Notice that the _Fork_ button on the upper right has incremented by one, and there is also is a line relating your fork to the original repo:

<img src="https://github.com/ProjectPythiaTutorials/github-arm-2022-02-28/blob/main/tutorial/images/Github_ForkPost.jpeg?raw=true" alt="ForkPost">

You now have a copy (essentially a clone) of the forked repository, which is now owned by you.

You could, at this point, select one of the files in the repository and use GitHub's built-in editor to make changes to these text-based files. But the typical use case that leverages the collaborative power of GitHub and its command-line cousin, _git_, involves _cloning_ your _forked_ copy of the repo to your local computer, where you can then perform your edits, and (in the case of software) test them on your system.


### What's Next?

In the next lesson, you will make changes to your _fork_ of a repository, and submit a _pull request_ with the changes!

## References

1. [What the Fork?(GitHub Community)](https://github.community/t/what-the-fork/10187)