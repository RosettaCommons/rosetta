This directory contains files with the sequences of glycans commonly seen in
biology.

Poses can be "glycosylated" with any of the glycans in this folder using:
core::pose::carbohydrates::glycosylate_pose_by_file(
		Pose & pose,
		core::uint const sequence_position,
		std::string const & filename )

Files may be in either the GlycoWorkbench format (.gws) or the IUPAC format
(.iupac).  To add new glycans to this folder, ensure that you have followed one
of the formats correctly.

Brief Explanation of the IUPAC Format:
* General:
  - Prefixes apply to the residue to which they are attached, below indicated
    by residue n.
  - Residues are listed from N to 1, where N is the total number of residues in
    the saccharide.
  - The sequence is parsed by reading to the next hyphen, so hyphens are
    crucial.
* Linkage indication:
  - "(a->x)-" specifies the linkage of residue n, where a is the anomeric
    carbon number of residue (n+1) and x is the oxygen number of residue n.
  - The first residue listed in the annotated sequence (residue N) need not
    have the linkage prefix.  A ->4) ResidueType will automatically be assigned
    by default if not specified.
* Anomer indication:
  - The strings "alpha-" (or "a-") or "beta-" (or "b-") are supplied next,
    which determine the stereochemistry of the anomeric carbon of the residue
    to which it is prefixed.
  - An alpha ResidueType will automatically be assigned by default.
* Stereochemical indication:
  - "L-" or "D-" specifies whether residue n is an L- or D-sugar.
  - The default is "D-".
* 3-Letter code:
  - A three letter code (in sentence case) MUST be supplied next.  This
    specifies the "base sugar name", e.g., Glc is for glucose.
  - A list of all recognized 3-letter codes for sugars can be found in
    database/chemical/carbohydrates/codes_to_roots.map.
* 1-Letter suffix:
  - If no suffix follows, residue n will be linear.
  - If a letter is present, it indicates the ring size, where
    * "f" is furanose,
    * "p" is puranose, and
    * "s" is septanose.
* Branches:
  Branches are indicated using nested brackets and are best explained by
  example:
  beta-D-Galp-(1->4)-[alpha-L-Fucp-(1->3)]-D-GlcpNAc- is:
  
  beta-D-Galp-(1->4)-D-GlcpNAc-
                         |
       alpha-L-Fucp-(1->3)

For more information, see:
www.chem.qmul.ac.uk/iupac/2carb
