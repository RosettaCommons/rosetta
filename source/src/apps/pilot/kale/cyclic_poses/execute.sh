#!/usr/bin/env sh

LOAD_CYCLIC_POSE="$HOME/rosetta/source/bin/load_cyclic_pose.linuxgccdebug"
KICK_CYCLIC_POSE="$HOME/rosetta/source/bin/kick_cyclic_pose.linuxgccdebug"

TOPOLOGY='loop'
INPUT_PDB="-in:file:s structures/ideal_$TOPOLOGY.pdb"
OUTPUT_PDB="-loops:output_pdb structures/parsed_$TOPOLOGY.pdb"

if [ "_$2" = "_verbose" ]; then
    VERBOSITY=""
else
    VERBOSITY="-mute all -unmute apps.pilot.kale"
fi

case "_$1" in 

    ("_load")   $LOAD_CYCLIC_POSE $VERBOSITY $INPUT_PDB $OUTPUT_PDB     ;;
    ("_kick")   $KICK_CYCLIC_POSE $VERBOSITY                            ;;
    (*)         echo "Usage: ./execute <load|kick> [verbose]"

esac
