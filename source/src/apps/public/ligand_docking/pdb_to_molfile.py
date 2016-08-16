#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import sys, os
from optparse import OptionParser

# Magic spell to make sure Rosetta python libs are on the PYTHONPATH:
def dirup(path, n=1):
    for i in xrange(n): path = os.path.dirname(path)
    return path
sys.path.append(os.path.join( dirup(os.path.abspath(sys.path[0]), 3), "python" ))

from rosetta_py.io.mdl_molfile import *
from rosetta_py.io import pdb


def main(argv):
    '''
    Extracts ligand coordinates from Rosetta PDBs to (re-)generate .mol2 or .sdf files.
    Requires an input .mol2 or .sdf file as a template, with the same atoms in the same order
    as was used to generate the Rosetta parameter files.
    '''
    parser = OptionParser(usage="usage: %prog [FLAGS] TEMPLATE.{mol,sdf,mol2} MODEL1.pdb [MODEL2.pdb ...] > OUTPUT.{sdf,mol2}")
    parser.set_description(main.__doc__)
    # parser.add_option("-short", ["--long"],
    #   action="store|store_true|store_false",
    #   default=True|False|...
    #   type="string|int|float",
    #   dest="opt_name",
    #   help="store value in PLACE",
    #   metavar="PLACE",
    # )
    parser.add_option("-i", "--input-format",
      default=None,
      help="File format for the template, in case it can't be deduced from the file extension.  Legal values are MOL2 and SDF.",
      metavar="FMT",
    )
    parser.add_option("-o", "--output-format",
      default=None,
      help="File format for the output, in case it can't be deduced from the file extension.  Legal values are MOL2 and SDF.  Default is same as input.",
      metavar="FMT",
    )
    parser.add_option("-c", "--comment",
      action="store_true",
      default=False,
      help="Write PDB filename into molfile comment field.",
    )
    parser.add_option("-t", "--title",
      action="store_true",
      default=False,
      help="Write PDB filename into molfile title field.",
    )
    (options, args) = parser.parse_args(args=argv)

    def err(*args, **kwargs):
        sys.stderr.write(*args, **kwargs)
        sys.stderr.write("\n")

    if len(args) >= 2:
        template_file = args[0]
        model_files = args[1:]
    else:
        parser.print_help()
        err("Wrong number of arguments!")
        return 1

    if options.input_format is None:
        f = template_file.lower()
        if f.endswith(".mol2"): options.input_format = "MOL2"
        elif f.endswith(".sdf") or f.endswith(".mol"): options.input_format = "SDF"
    if options.output_format is None:
        options.output_format = options.input_format

    if options.input_format == "MOL2":
        template = next(read_tripos_mol2(template_file, do_find_rings=False))
    elif options.input_format == "SDF":
        template = next(read_mdl_sdf(template_file, do_find_rings=False))
    else:
        parser.print_help()
        err("Unknown input format '%s'" % options.input_format)
        return 2
    # SDF files don't (generally) have unique atom names.
    # Using the same algorithm as molfile_to_params should ensure that names match up.
    uniquify_atom_names(template.atoms)

    outfile = sys.stdout
    if options.output_format == "MOL2":
        out_writer = write_tripos_mol2
    elif options.output_format == "SDF":
        out_writer = write_mdl_sdf
    else:
        parser.print_help()
        err("Unknown output format '%s'" % options.output_format)
        return 3

    # Need the calls to .upper() to treat e.g. "Br" and "BR" the same
    atoms = dict([(pdb_pad_atom_name(a).upper(), a) for a in template.atoms])
    for model_file in model_files:
        if options.comment: template.comment = model_file
        if options.title:   template.title   = model_file
        changed = 0
        #changed_names = set()
        for het in pdb.get_het_atoms( pdb.read_pdb_file(model_file) ):
            name = het.name.upper()
            if name in atoms:
                atom = atoms[name]
                atom.x = het.x
                atom.y = het.y
                atom.z = het.z
                changed += 1
                #changed_names.add(name)
        diff = changed - len(template.atoms)
        if diff < 0:
            err("Missing %i HETATM records in %s" % (-diff, model_file))
        elif diff > 0:
            err("%i duplicate HETATM records in %s" % (diff, model_file))
            # This generally occurs because we have multiple HET residues in the PDB,
            # and some of their atoms have the same names.
            # Nonetheless, this function should Do The Right Thing,
            # because later records overwrite earlier ones, and
            # the dockable ligand always comes last in ligand-docking PDB files.
        out_writer(outfile, [template])

    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
