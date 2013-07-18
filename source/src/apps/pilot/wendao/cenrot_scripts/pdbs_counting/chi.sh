for fn in `cat ../all.lst`
do
#	~/rosetta/rosetta_source/cmake/build_release/fa_to_cenrot_score -database ~/rosetta/rosetta_database/ -ignore_unrecognized_res -ignore_zero_occupancy false -mute all -l $fn -unmute pilot.wendao.cenrot -output_cenrot_intcoord > refine/${fn}.chi 2> refine/${fn}.log
#	~/rosetta/rosetta_source/cmake/build_release/fa_to_cenrot_score -database ~/rosetta/rosetta_database/ -ignore_unrecognized_res -ignore_zero_occupancy false -l $fn
done

