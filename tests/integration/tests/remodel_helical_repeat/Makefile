default:
	~/src/Rosetta/main/source/bin/rosetta_scripts.hdf5.linuxgccrelease -nstruct 250 -out:file:silent 1.out -parser:protocol layeredDesign.xml -num_trajectory 1 -save_top 1 @flags -remodel:staged_sampling true -out:level 400 
clean:
	rm 1.out 1.pdb tmp
