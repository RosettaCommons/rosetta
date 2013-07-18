~/rosetta/rosetta_source/cmake/build_release/general_pair_counting -database ~/rosetta/rosetta_database/ -corrections:score:cenrot -l nonideal_* > ../test/pair/5.cen.pair.seq8.rot.test
~/rosetta/rosetta_source/cmake/build_release/general_pair_counting -database ~/rosetta/rosetta_database/ -corrections:score:cenrot -l nonideal_* -env_dependent_pair 16 40 > ../test/pair/5.cen.pair.seq8.rot.test


#5/17/2013 now we use new general_pair_counting code, which would be a little differ out side the cutoff
