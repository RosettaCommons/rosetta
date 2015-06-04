This directory contains lists of alternative 3-letter codes for use in PDB
files and their corresponding Rosetta equivalents.

Files in this directory should be named *.codes for consistency.

Each file must contain at least two columns:
* The first column is an alternate code.
* The second column is the Rosetta 3-letter code.
* The third column (optional) is a HETNAM record designation for that code.

CODES ARE CASE-SENSITIVE!

To add any of these sets of alternative codes to Rosetta, use the
-alternate_3_letter_codes <FILENAMES> option flag.
