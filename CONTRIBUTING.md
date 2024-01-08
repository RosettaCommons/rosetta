Contributing to Rosetta
=======================

Starting Rosetta Developement
-----------------------------

Rosetta uses a Fork-and-PR model for development. To start developing with Rosetta, use the "Fork" button on the RosettaCommons repo to create a new Rosetta fork in your own Github space. (Note that if you wish to alter the content of submodules, you will need to separately fork the submodule.)

You can then clone the repo from your space, start a new branch and edit your code. Feel free to push the code to your forked repo at any time - it will not affect the main Rosetta repository.

When you feel your contributions are ready to be incorporated into the main version of Rosetta, use the Pull Request (PR) feature of Github to open a PR, targeting your branch on your repo to merge into the `main` branch of the RosettaCommons/rosetta repo.

If you haven't already, you will need to sign the Rosetta Contributor Licensing Agreement (See below for details).

Once you've opened your PR, members of the Rosetta development team will take a look at the code and check to make sure it works in the wider framework of Rosetta. They may ask you to make some changes to alter your code to fit with Rosetta conventions, or to fix edge cases. Once your code looks good, they will squash-merge the code into the Rosetta main branch. (Squash-merge means that all of your changes will be re-written to be a single commit.)

See <<<LINK>>> for more information about the recommended workflow.

Contributor Licensing Agreement
-----------------------------

To maintain RosettaCommons's ability to license Rosetta commercially, we ask that people sign a Contributor Licensing Agreement (CLA). (See CLA.md for the current text)

Please read and make sure you understand the agreement prior to signing, as it is a legal agreement. 
You only have to sign it once, as it will apply to all the contributions you submit to the RosettaCommons as PRs.

If you are making contributions as an employee (which includes employees and graduate students of academic institutions), you may have assigned rights for code you develop to your employer. 
In that case, please confirm with your employer that you are authorized to sign the CLA on their behalf.
If you change employers, please re-sign the agreement with your new employer's information.

To sign the CLA, please visit <<CLA WEBSITE>>>.

If you are a member of a RosettaCommons lab, there is a separate Developers Agreement process. 
Please consult an existing member of the lab or your PI for details.

Submodules
----------

Rosetta uses a number of git submodules to hold code and materials which are generally not needed for rountine development and usage, or which would be too large to include directly. Submodules which are needed for building Rosetta should be automatically downloaded by the Rosetta build system, if they haven't already been.

<< QUESTION -- How does submodule locations work with forking? >>

Note that if you wish to alter the content of one of the submodules, you will need to separately fork the submodule and make your modifications in that forked repo. You would then open a PR against the `main` branch of that submodule to incorporate.

Merging the PR on the submodule does not immediately make it availible from the RosettaCommons/rosetta repository -- a separate PR to update the RosettaCommons/rosetta repository will be necessary.

Contributing FAQ
----------------

Q. *Why don't you use the OpenSource-standard "inbound=outbound" model for licensing?*

A. To support RosettaCommons's ability to maintain Rosetta, the RosettaCommons currently charges a fee for commercial usage. 
All profits from the sale of Rosetta go toward supporting the scientific mission of the RosettaCommons.
At the moment, we currently don't see a way to go officially (OSI) Open Source while keeping the same level of support.
As a consequence, the Rosetta license is written in a way which privleges the RosettaCommons and their ability to charge a commercial license fee.
However, licensing inbound contributions under that same license would result in an incompatible mess, particularly if we ever attempted to change the license in the future.

As such, the current code base is licensed under the existing licensing framework, but we ask new contributions to be submitted under more permissive terms, to allow flexibilty for future license changes.
 
Q. *Can I use the Rosetta code in my own software?*

Use of Rosetta code in other projects is governed by the Rosetta license agreement. The combination of Rosetta with another program must maintain the Rosetta licensing terms for those portions which include Rosetta. In particular, use of Rosetta (or PyRosetta) as a library in another program does not negate the requirement for purchasing a separate Rosetta (and/or PyRosetta) license for commercial use of the combined program. Please contact UW CoMotion <license@uw.edu> with any questions about commercial licensing requirements for programs containing Rosetta code.

We note specifically that the Rosetta license is incompatible with the GPL and other such license with "share alike" provisions.

Q. *Can I charge a license fee for my own code?*

Under the Rosetta CLA, you retain rights to the code you develop specifically. You are free to provide that code to others under whatever terms you desire.
However, this only extends to your code in isolation. You may only distribute the combination with Rosetta code (or Rosetta-derived code) under terms compatible with the Rosetta license. Charging a fee for use or distribution of programs containing Rosetta code would count as "commercial use" under the Rosetta license. Please contact UW CoMotion <license@uw.edu> to discuss commercial licensing terms.

Q. *Can I include external libraries into Rosetta?*

Other libraries can be included/linked against Rosetta, provided they have a compatible license. MIT and OSI approved BSD licenses are known to be compatible. GPL and other such license with "share alike" provisions are not. Please open a Github issue to discuss details, as including an additional external library often has complexity issues over and above the licensing issues.

