To add a new property or variant type to Rosetta for use in .params
files, patch files, and ResidueTypes, simply add a string to the
general_properties.list or variant_types.list files, respectively, and run
update_ResidueType_enum_files.py.

(You may also consider adding new is_ accessors to ResidueType and Residue.)
