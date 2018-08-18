Getting Started with Rosetta
============================
Start here: https://www.rosettacommons.org/docs/latest/getting_started/Getting-Started

Rosetta/main (the repository you are looking at) contains the Rosetta source code, database, unit tests and integration tests. The source code is located in source/src can be compiled with SCons using the following commands:

``` sh
$ cd Rosetta/main/source
$ ./scons.py -j<NumOfJobs> mode=[debug/release] [bin]
```

## Notes to Developers
### Don't Fork Rosetta Main
In our current development workflow, all work should be done in branches. Do not fork this repository.  Read [our workflow documentation](https://www.rosettacommons.org/docs/wiki/internal_documentation/GithubWorkflow) to learn how to work with the Rosetta repositories.

We avoid forking, using a branching-only method, to ensure that all the Rosetta code stays accessible.  We would lose access to unmerged forked branches as developers leave the community.  It also ensures we have a single backup of all pushed branches in the event that GitHub were to disappear.

### Interacting with Pull Requests
All changes to master must go through Pull Requests (PR), also described in [our workflow documentation](https://www.rosettacommons.org/docs/wiki/internal_documentation/GithubWorkflow).  Those PRs also get reviewed by the community before merge, [here](https://www.rosettacommons.org/docs/wiki/internal_documentation/GithubWorkflow#workflow-for-using-github_pull-request-pr_pull-request-review-what-to-do-before-review) are tips on structuring your code for easy review.
