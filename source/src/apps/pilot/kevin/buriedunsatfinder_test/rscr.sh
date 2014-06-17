/Users/khouli/Rosetta/main/source/cmake/build_release/rosetta_scripts \
	-database /Users/khouli/Rosetta/main/database \
	-score:weights talaris2013 \
	-l pdbs.list \
	-parser:protocol bunsat2filt.xml \
	-unmute devel.buns.BuriedUnsatisfiedPolarsCalculator2 \
	-unmute devel.vardist_solaccess \
	-overwrite \
	-out:nooutput \
	-out:pdb_gz \
	-ignore_unrecognized_res
#-unmute devel.buns.BuriedUnsatHbondFilter2 \
#-out:nooutput \
