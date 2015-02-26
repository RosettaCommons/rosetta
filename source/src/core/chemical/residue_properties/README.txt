To add a new property or variant type to Rosetta for use in .params
files, patch files, and ResidueTypes, simply add a string to the
general_properties.list or variant_types.list files, respectively, and run
update_ResidueType_enum_files.py.

(You may also consider adding new is_ accessors to ResidueType and Residue.)


Suggested/Strongly Encouraged Use of Properties and Variants:

Properties:
* A property should be used for the following purposes:
  - To designate whether the residue in question is an actual polymeric residue
    or a ligand, e.g., POLYMER, LIGAND
  - To designate to which chemical classification/family the residue belongs,
    e.g., PROTEIN, RNA
  - To designate a global chemical property of the molecule/residue, e.g.,
    POLAR, CHARGED
  - To designate that a residue has been modified at a particular position,
    (but NOT what that modification is,) e.g., OG_MODIFIED, C1_MODIFIED
* Do not combine properties that are distinct: AA and L_STEROCHEMISTRY would
  be preferred properties; L_AA would not be.
* Properties should be named as nouns if a family designation or adjectives
  if a chemical property, e.g., PEPTOID vs. AROMATIC.

Variants:
* A variant should be used to designate a specific variation/modification of a
  base residue, specifying the location if necessary, e.g., LOWER_TERMINUS,
  OG_PHOSPHORYLATED, C1_PHOSPHORYLATED
* Variants are NOT directly synonymous with patches:
  - It is possible for several patch files to be of the same variant.  For
    example, the way to make a LOWER_TERMINUS variant of an AA residue is very
    different from how one makes a LOWER_TERMINUS of an NA residue.
  - It should NOT be possible for a single base residue to have more than one
    patch file specifying the same variant type!

Do not duplicate information: No property should also be a variant and vice
versa.
