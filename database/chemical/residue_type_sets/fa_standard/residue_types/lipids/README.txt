I am adding a few lipids so that I can begin modeling glycolipids.  If anyone
after me wishes to add to this database, take note of the following, because
some of the choices I made in making fatty acid and other lipid residues are
not straightforward:

All the lipids here are treated as POLYMER (not LIGAND) types.  This is because
POLYMER in Rosetta does not mean the same thing as it does in the real world.
If the residue connects to another residue, it needs to be treated as a
POLYMER.

Head groups (e.g., glycerol, inositol, sphingosine, etc.) are treated as lower
termini only.  Thus, they do NOT contain LOWER_CONNECT records and do not
(cannot) be patched with a LOWER_TERMINUS patch.  The main chain and all
branches must extend from the head group in the fold tree.  For this reason, I
have chosen to run the main chain from C-omega (e.g., C18 in sphingosine) in
long chain head groups toward C1.  For short chain head groups (e.g., glycerol)
the main chain will begin at C1.  This avoids having to have a dozen or more
chi angles defined.

Some head groups also include CONNECT records to define branches directly in
their .params files.  This is different from how I set up carbohydrates.  The
reason I used patch files for this in carbohydrates is that not all hydroxyls
in the sugars have branches.  If a branch point usually has a branch, one will
need to use a branch to convert it into an unmodified hydroxyl/amine, similar
to what one does for an UPPER_ or LOWER_TERMINUS.  In other words, some head
groups (such as sphingosine) can be thought of as having TWO UPPER_CONNECTs,
but the only way to do this in Rosetta is with CONNECT records.

Fatty acids are treated the opposite way; they are upper termini only.  As
such, they lack UPPER_CONNECTs and do not need to (cannot) be patched with an
UPPER_TERMINUS patch.

All of these decisions are of course up for debate, but please notify me if you
wish to change anything, because I will soon be modeling glycolipids with the
framework as here described.

~Labonte
<JWLabonte@jhu.edu>
