#!/usr/bin/env sh

cd -P .

pose=1srp
fixbb=../../../../../../bin/fixbb

$fixbb                              \
    -in:file:fullatom               \
    -in:file:s $pose.pdb            \
    -out:suffix $pose.min.pdb       \
    -min_pack                       \
    -packing:repack_only            \
    -ignore_unrecognized_res        \
    -overwrite                      \
    -nstruct 1

    #-ex1 -ex2 -extrachi_cutoff 0
