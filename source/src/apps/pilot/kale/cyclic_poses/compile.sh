#!/usr/bin/env sh

case "$1" in

    (load)      ~/rosetta/source/ninja_build.py load_cyclic_pose        ;;
    (kick)      ~/rosetta/source/ninja_build.py kick_cyclic_pose        ;;
    (*)         echo "Usage: ./compile.sh load|kick"                    ;;

esac
