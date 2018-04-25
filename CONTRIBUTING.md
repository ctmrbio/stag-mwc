# Contributing guidelines for StaG-mwc

## Issue tracker
We use the issue tracker in Github. Submit issues for things such as
bug reports, feature requests, or general improvement discussion topics.

## Submitting changes
The typical procedure to develop new features or fix bugs in StaG-mwc looks
something like this:

1. Fork or clone the repository.
2. Create a branch with a descriptive name based on your intended changes using
   dashes to separate words, e.g. branch-to-add-megahit-assembly-step
3. Insert your code into the respective folders, i.e. scripts, rules and envs.
   Define the entry point of the workflow in the Snakefile and the main
   configuration in the config.yaml file.
4. Commit changes to your fork/clone.
5. Create a pull request (PR) with some motivation behind the work you have
   done and possibly some explanations for tricky bits.
