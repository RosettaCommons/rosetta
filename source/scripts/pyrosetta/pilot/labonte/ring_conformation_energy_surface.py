# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available
# (c) under license.
# (c) The Rosetta software is developed by the contributing members of the
# (c) Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org.
# (c) Questions about this can be addressed to University of Washington UW
# (c) TechTransfer, email: license@u.washington.edu.
"""Brief:   This PyRosetta script selects and applies every ring conformation
to, packs, and scores a monosaccharide pose.  The scores are output as a tab- delimited table of C-P parameters theta vs. phi for use in generating relative
energy surface diagrams for ring conformation moves.

Author:  Jason W. Labonte

"""
# Imports
from os import mkdir
from argparse import ArgumentParser

from rosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    create_score_function, standard_packer_task, Pose, \
                    MoveMap, RotamerTrialsMover, Vector1
from rosetta.core.chemical.carbohydrates import RingConformer, RingConformerSet

# Parse arguments.
parser = ArgumentParser(description=__doc__)
parser.add_argument('pdb_filename',
                    help='the filename of the PDB structure to be evaluated;' +
                    ' must contain a single saccharide residue')
parser.add_argument('--make_movie', action='store_true',
                    help='flag to output a movie frame for each conformer')
parser.add_argument('--mm', action='store_true',
                    help='flag to use the molecular mechanics score function')
args = parser.parse_args()

# Initialize Rosetta.
init(extra_options='-include_sugars -mute all')

# Load pdb file.
initial_pose = pose_from_pdb(args.pdb_filename)

if initial_pose.total_residue() != 1:
    exit('The pdb file must contain a single saccharide residue.')

if not initial_pose.residue(1).is_carbohydrate():
    exit('The pdb file must contain a single saccharide residue.')

info = initial_pose.residue(1).carbohydrate_info()

print 'Generating energy surface for', info.anomer() + "-" + \
                                                       info.full_name() + '...'

# Set up ScoreFunction.
if args.mm:
    sf = create_score_function('mm_std')
else:
    sf = get_fa_scorefxn()

print ' Initial Score:', sf(initial_pose)

# Set up packer.
pt = standard_packer_task(initial_pose)
pt.restrict_to_repacking()
pt.or_include_current(True)

packer = RotamerTrialsMover(sf, pt)

# Prepare data storage.
surface = []
params = []

if args.make_movie:
    dir_name = args.pdb_filename[:-4] + '_movie_frames'
    try:
        mkdir(dir_name)
    except OSError:
        print 'Warning: Directory already exists; files will be overwritten.'

# Prepare temp pose.
pose = Pose()

# Switch based on ring size.
ring_size = info.ring_size()
if ring_size == 3:
    exit(' Functionality for 3-membered rings not coded yet; exiting')
elif ring_size == 4:
    exit(' Functionality for 4-membered rings not coded yet; exiting')
elif ring_size == 5:
    exit(' Functionality for 5-membered rings not coded yet; exiting')
elif ring_size == 6:
    # Generate header for data output.
    header = "\t"
    for phi in range(0, 360 + 1, 30):
        header += str(phi)
        if phi == 360:
            header += "\n"
        else:
            header += "\t"
    surface.append(header)

    # Generate data.
    for theta in range(180, 0 - 1, -45):
        data_row = str(theta) + "\t"
        for phi in range(0, 360 + 1, 30):
            # Reset pose.
            pose.assign(initial_pose)

            # Get and set conformer.
            params = Vector1([0.5, phi, theta])  # q can be any positive real.
            set_ = info.ring_conformer_set()
            conformer = set_.get_ideal_conformer_by_CP_parameters(params)
            pose.set_ring_conformation(1, conformer)

            # Pack conformer.
            packer.apply(pose)

            name = str(conformer)[:str(conformer).find(" ")]  # HACK
            print ' Generated conformer:', name,

            if args.make_movie:
                pose.dump_pdb(dir_name + "/" + name + ".pdb")

            # Output score.
            score = str(sf(pose))
            print ' Score:', score
            data_row += score
            if phi == 360:
                data_row += '\n'
            else:
                data_row += '\t'
        surface.append(data_row)

# Output results.
data_filename = args.pdb_filename[:-4] + '_ring_conformation_energy_surface'
if args.mm:
    data_filename += '_mm'
data_filename += '.data'
with open(data_filename, 'w') as f:
    f.writelines(surface)

print 'Data written to:', data_filename