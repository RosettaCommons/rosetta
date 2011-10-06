
########################
How to get this sample_source setup:

1) put directory of pdbs into the input/ directory, e.g. 
   ln -s /path/to/datasets/top4400/input/top4400pdbs .

2) run setup_resfiles.py to generate rosetta_inputs.db3:
   python setup_resfiles.py --data_dir top4400pdbs

3) run convert_htam_names.py to convert Reduce names to hydrogen names:
   python convert hatm_names.py --data_dir top4400pdbs --output_dir top4400pdbs_rosetta_named_hydrogens

4) setup the list of structures to use:
   ls top4400pdbs_rosetta_named_hydrogens/ > all_pdbs.list' 

5) if you didn't use the path 'top4400pdbs_rosetta_named_hydrogens', edit the '-in:path' option in in flags.TEMPLATE 

6) add this sample_source to features/sample_source/benchmark.list

##########################

Note: You should use the top8000 set instead. It is newer and some
problems have been cleaned up.
