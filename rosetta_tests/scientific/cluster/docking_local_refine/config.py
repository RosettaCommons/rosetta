{
	"nstruct.description": "Number of predictions to make for each target",
	"nstruct": 180,

	"extra_chi_rotamers.description": "Flags used to control how much rotamer sampling should be done",
	"extra_chi_rotamers": "-ex1\n-ex2\n-extrachi_cutoff 0",

        "energy_function_flags.description": "Indicate energy function parameters",
        "energy_function_flags": "",

	"svn_benchmark_path.description" : "svn path to the docking benchmark set and where it should get checked out",
	"svn_benchmark_path": "https://svn.rosettacommons.org/source/trunk/mini.data/tests/scientific/cluster/docking",

	"benchmark_data_path.description": "Where the input data should be prepared",
	"benchmark_data_path": "inputs",

	"benchmark_data_set.description" : "Look in '<svn_benchmark_path>/<benchmark_data_set>.tar.gz for the input data",
	"benchmark_data_set": "zdock_benchmark_v4",

	"input_target_path.description" : "where the fully preprocessed input targets should end up",
	"input_target_path": "inputs/bound_targets",

	"receptor_target_extension.description" : "r <- 'receptor': the larger partner",
	"receptor_target_extension": "_r_b.pdb",

	"ligand_target_extension.description" : "l <- 'ligand': the smaller partner",
	"ligand_target_extension": "_l_b.pdb",

	"output_run_log_path.description": "Where the log files should be generated.",
	"output_run_log_path": ".",

	"output_decoy_path.description": "Where the decoys should be generated.",
	"output_decoy_path": "decoys",
	"output_silentfile_extension.description": "The output silentfiles will be <output_decoy_path>/<target><input_target_extension><output_silentfile_extension>.",
	"output_silentfile_extension": "_rl_b.silent.gz",

	"output_score_path.description": "Where the score files should be generated.",
	"output_score_path": ".",
	"output_scorefile_extension.description": "The output scorefiles will be <output_score_path>/<target><input_target_extension><output_scorefile_extension>.",
	"output_scorefile_extension": "_rl_b.sc",

        "output_analysis_path.description": "Where the analysis plots should be generated.",
        "output_analysis_path": "files",

	"test_results_log": ".results.log",
	"test_results_yaml": ".results.yaml",

	"condor_priority.description": "The priority of the condor job",
	"condor_priority": -10,

	"condor_queue.description": "Number of cores to request with condor job",
	"condor_queue": 50
}
