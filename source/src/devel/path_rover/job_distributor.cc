// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///@author sheffler

// MPI
// jk Note: if MPI is used, it must be included BEFORE some other std libraries
// jk We'll keep the include here for safety since some Rosetta headers include other std libs
#ifdef USEMPI
#include "mpi.h"
#include "random_numbers.h"
#include <ObjexxFCL/Time_Date.hh>
#define NUM_MPI_ARGS 4 // number of ints sent in each MPI communication
#endif // MPI

// Rosetta Headers
#include "job_distributor.h"
#include "after_opts.h"
#include "abrelax.h"
#include "analyze_interface_ddg.h"
#include "antibody_modeling.h"
#include "assemble_domains.h"
#include "barcode_stats.h"
#include "constraints.h"
#include "counters.h"
#include "crankshaft.h"
#include "csa_main.h"
#include "design.h"
#include "decoystats.h"
#include "decoy_features.h"
#include "design_structure.h"
#include "dipolar.h"
#include "disulfides.h"
#include "dock_structure.h"
#include "docking.h"
#include "docking_ns.h"
#include "evolve.h"
#include "files_paths.h"
#include "filters.h"
#include "flexpep_dock.h"
#include "fold_abinitio.h"
#include "fold_loops.h"
#include "fold_membrane.h"
#include "force_barcode.h"
#include "fragments.h"
#include "fullatom.h"
#include "hbonds.h"
#include "idealize.h"
#include "initialize.h"
#include "input_pdb.h"
#include "jumping_minimize.h"
#include "jumping_pairings.h"
#include "knots.h"
#include "ligand.h"
#include "ligand_ns.h"
#include "loops.h"
#include "loop_relax.h"
#include "maps.h"
#include "maps_ns.h"
#include "map_sequence.h"
#include "minimize.h"
#include "misc.h"
#include "monte_carlo.h"
#include "namespace_cold.h"
#include "namespace_options.h"
#include "native.h"
#include "options.h"
#include "orient_rms.h"
#include "output_decoy.h"
#include "pathways.h" // Barak & Angela 29/11/2007
#include "pack_fwd.h"
#include "packing_measures.h"
#include "param.h"
#include "pdbstats.h"
#include "pH_main.h"
#include "pKa_mode.h"
#include "pose.h"
#include "pose_io.h"
#include "pose_looping.h"
#include "pose_main.h"
#include "prof.h"
#include "DomainInsertionMode.h"
#include "ramachandran.h"
#include "read_aa_ss.h"
#include "read_aaproperties.h"
#include "read_paths.h"
#include "recover.h"
#include "refine_structure.h"
#include "refold.h"
#include "relax_structure.h"
#include "repeat.h"
#include "runlevel.h" //For benchmark.
#include "score.h"
#include "start.h"
#include "status.h"
#include "structure.h"
#include "taboo_search.h"
#include "trajectory.h"
#include "vdw.h"

//Utility Headers
#include <utility/basic_sys_util.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/gzip_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

// C++ Headers
#include "time.h"

//Graphics
#ifdef GL_GRAPHICS
#include "gl_graphics.h"
#include <pthread.h>
#endif

// BOINC
#ifdef BOINC
#include "boinc_rosetta_util.h"

#ifdef _WIN32
#include "boinc_win.h"
#endif
#include "diagnostics.h"
#include "boinc_api.h"
#include "util.h"
#include "boinc_rosetta_util.h"
#include "boinc_rosetta_graphics.h"
#include "watchdog.h"

#ifdef BOINC_GRAPHICS
#include "graphics_api.h"
#include "graphics_lib.h"
#include "graphics_data.h"
#endif

#define BOINC_MAX_NO_PROGRESS_INIT_CNT 5 // maximum allowed restarts w/ no progress
#endif

namespace main_job_distributor {
	job_distributor * jd;
}


/*

// jk Sample code demonstrating usage of "loop_over_function", which can be put into any file

#include "job_distributor.h"

// Definition of DATA_CLASS, which holds data needed for evaluation of "eval"
class DATA_CLASS {
public:
	int a;
	DATA_CLASS();
	inline DATA_CLASS(int const ina) { a=ina; return; };
	inline ~DATA_CLASS() { return; };
};

// Function which will be evaluated with various values of i, j, k
// Note: instance f of DATA_CLASS is available, but only as const
void eval( const void * in_dc, int const i, int const j, int const k ) {
	const DATA_CLASS * dc = (DATA_CLASS*) in_dc;
	std::cout << "value of a, i inside eval are: " << dc->a << ' ' << i << std::endl;
}

// Instantiation of a DATA_CLASS object, invokation of "loop_over_function"
DATA_CLASS f( 5 );
main_job_distributor::jd->loop_over_function( eval, &f, 5 );

*/


////////////////////////////////
///  job_distributor classes
////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::type
///
/// @brief:  Return the job_distributor type
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 11/12/06
/////////////////////////////////////////////////////////////////////////////////
std::string job_distributor::type() {

	return "base";

}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::initialize
///
/// @brief: setup the job distributor
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::initialize() {
	using namespace param;

	startnm.dimension( MAX_START );
	outnm.dimension( MAX_START );
	startch.dimension( MAX_START);
	startch = ' ';

	initialize_rosetta();

	num_decoys=0;
	attempted_decoys=0;
	starting_pdbs_skipped=0;

	curr_startnum=0;
	curr_outnum=number_of_output;
	reinitialize=false;

	// for pose Monte_carlo checkpointing
	mc_checkpoint_last_count = 0;
	mc_checkpoint_current_count = 0;
	mc_checkpoint_file = "mc_checkpoint";
	update_mc_checkpoint_time();

	return;

}


////////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::initialize_rosetta
///
/// @brief
///
/// @detailed
///
/// @global_read
/// atom_vdw_set files_paths.h
///
/// @global_write
/// atom_weights static
///
/// @remarks
/// many indirect changes to global variables
///
/// @references
///
/// @authors car 8/18/2003
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::initialize_rosetta() {

	using namespace cold;
	using namespace files_paths;
	using namespace param;

	output_command();
	add_args_from_file();
	read_paths();

	get_rosetta_options( mode, number_of_output );
	setup_start_list( nstartnum, startnm, outnm, startch );

	read_rama(); // Ramachandran table
	structure::SecondaryStructureEval::load_and_cache_phi_theta_bins_from_file(); // secondary structure scoring terms
	init_fold(); // precompute refold variables
	setup_atomvdw(); // select default atom radii
	initialize_aaproperties(); // Set amino acid properties
	read_residue_paircutoffs();
	select_rotamer_set( "default" ); // initialize rotamer set parameters

	minimize_reset();
	initialize_torsion_logicals();
	setup_decoystats();
	packing_ns::initialize_packing_measures();
	decoy_features_ns::decoy_features_initialize();
	initialize_taboo();
	initialize_evolve();
	prof::reset(); // for profiling on the fly

}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::check_pose1
///
/// @brief: used when mode == pose1
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
bool job_distributor::check_pose1() {

	if ( mode == "pose1" ) return true;
	return false;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::run_main_pose1
///
/// @brief: used when mode == pose1
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::run_main_pose1() {

	main_pose1( number_of_output, nstartnum );
	if (runlevel_ns::benchmark) finish(); //Print out the silly "DONE" statement.
	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::check_csa
///
/// @brief: used for csa mode
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
bool job_distributor::check_csa() {

	if (truefalseoption("clnt" )) return true;  // Rosetta client mode
	return false;

}



//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::run_csa
///
/// @brief: used for csa mode
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::run_csa() {

	using namespace files_paths;

	protein_chain = startch(1);
	rosetta_client_mode_main(mode);
	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::skip_starting_pdb
///
/// @brief: called when starting coors cannot be obtained
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::skip_starting_pdb() {

	std::cout << "starting coordinates were not obtained" << std::endl;
	std::cout << "skipping this structure" << std::endl;
	++starting_pdbs_skipped;
	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::get_next_job_num
///
/// @brief: set the start_inx and output_inx for next job. Returns true if
///          more jobs remain, returns false if we've already done the last job
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
bool job_distributor::get_next_job_num() {

	++curr_outnum;
	if ( curr_outnum > number_of_output ) {
		++curr_startnum;
		if ( curr_startnum > nstartnum ) return false;
		curr_outnum=1;
	}

	return true;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::prepare_for_next_startnum
///
/// @brief: set the start_inx and output_inx such that
///      next time we'll move on to a new start_num
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::prepare_for_next_startnum() {

	curr_outnum = number_of_output;
	return;

}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::made_output_decoy
///
/// @brief: called when an output decoy has been made
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::made_output_decoy() {

	++num_decoys;
	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::done_attempt
///
/// @brief: called when an attempt at an output decoy has been completed,
///           but the output decoy filters have not yet been applied
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::done_attempt() {

	++attempted_decoys;
	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::call_abrelax
///
/// @brief: essentially a stub, required due to a derived class
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::call_abrelax() {

	abrelax(0);
	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::setup_job
///
/// @brief: setup for the next job. Return true if a job has been set up
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
bool job_distributor::setup_job() {

	using namespace files_paths;

	int prev_startnum(curr_startnum);

	if ( get_next_job_num() ) {

		std::string status("okay");
		if ( curr_startnum != prev_startnum ) {
			// jk we've moved on to the next input structure
			start_file = startnm(curr_startnum);
			output_file = outnm(curr_startnum);
			protein_chain = startch(curr_startnum);
			delete_START();
			initialize_start( mode, curr_startnum, number_of_output, status );
#ifdef GL_GRAPHICS
			open_windows();
#endif
			if ( status == "fail" ) {
				skip_starting_pdb();
			}
			prev_startnum = curr_startnum;
		}

		if ( status == "okay" ) {

			if (reinitialize) { // This is an abrelax thing, if we're switching between a homolog and a query.
				reinitialize_query(); // Go back to homolog sequence, fragments, etc.
				initialize_start( mode, curr_startnum, number_of_output, status );
				reinitialize = false;
			}
			return true;

		} else {  // status of reading input was not "okay"
			prepare_for_next_startnum();
			return setup_job(); // note: recursive!
		}

	} else {
		return false;
	}

	return true;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::run_job
///
/// @brief: call the appropriate function to run this job, given the mode
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::run_job() {

	using namespace files_paths;

	if (curr_outnum > number_of_output) return; // useful in score mode ... skip run_job.

	barcode_initialize_decoy(); // should happen before checking for previous decoy ("file_exists")

	reset_decoy_attempt_counter();

	bool file_exists(false);
	open_decoy_pdb(mode, curr_startnum, curr_outnum, file_exists);

	//rh Initialize pKa mode and fill dGprotonation table
	pKa_mode::initialize_pKa_start(curr_outnum);

	bool accepted(false);
	bool failed(false);
	while ( !accepted && !file_exists && !failed &&
					( ( !accept_all || no_prev_decoy_attempts() ) ) ) {

		increment_decoy_attempt_counter();
		initialize_decoy();

		repeat_begin_iter(curr_outnum); // added by sheffler 11/5/04

		if ( mode == "loops" ) {
			mode_title = "Fold loops";
			fold_loops();
		} else if ( mode == "assemble" ) {
			mode_title = "Domain assembly";
			assemble_domains();
		} else if ( mode == "refine" ) {
			mode_title = "Refinement";
			refine_structure();
		} else if ( mode == "design" ) {
			mode_title = "Design";
			design_structure();
		} else if ( mode == "dock" ) {
			mode_title = "Docking";
			dock_structure( failed );
		} else if ( mode == "flexpepdock") { // barak 12/Feb/2007
			mode_title = "Flexible Peptide Docking";
			flexpep_dock::invoke_flexpep_dock_from_misc(failed);
		} else if (mode == "pathways") { // barak & angela 29/Nov/2007
			mode_title = "Conformation Pathways Generator";
			pathways::pathways_generator_main(failed);
		} else if ( mode == "relax" ) {
			mode_title = "Relax";
			relax_structure(); // bqian: should we add attempted decoys too?
		} else if ( mode == "idealize" ) {
			mode_title = "Idealize";
			idealize(failed);
		} else if ( mode == "membrane" ) {
			mode_title = "Membrane";
			fold_membrane();
		} else if ( mode == "abinitio" ) {
			mode_title = "Ab initio";
			fold_abinitio();
		} else if ( mode == "pdbstats" ) {
			mode_title = "PDB stats";
			get_pdbstats();
		} else if ( mode == "interface" ) {
			mode_title = "Interface";
			analyze_interface_ddg();
		} else if ( mode == "bc_stats" ) {
			mode_title = "BC stats";
			get_barcode_stats();
		} else if ( mode == "pKa" ) {
			mode_title = "pKa";
			pKa_mode::main_protocol();
		} else if ( mode == "abrelax" ) {
			mode_title = "Ab initio + relax";
			call_abrelax();
		} else if (mode == "pose_looping" ) {
			mode_title = "Pose Looping";
			looping_main();
		} else if (mode == "domain_insertion" ) {
			mode_title = "Domain Insertion";
			domain_insertion_main();
		}
		else if (mode == "antibody_modeler" ) {
			mode_title = "Antibody Modeler";
			antibody_modeling();
		}

		if ( ! failed ) {
			done_attempt();
			prof::show();
			store_low_info();
			output_decoy(accepted);
			if (homolog_to_query_mapping) reinitialize = true;
			repeat_end_iter(curr_outnum,accepted); // added by sheffler 11/5/2004
		}

	}

	if ( !file_exists && !failed ) made_output_decoy();
	report_progress();

	return;

}

void job_distributor::report_progress()
{
	return;
}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::finish
///
/// @brief: end Rosetta (cleanly)
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::finish() {

	delete_START();
	prof::show();

	finish_writing_files();
	finish_summary( std::cout );

	repeat_final(); // added by sheffler 11/5/04
	decoy_features_ns::decoy_features_final();

#ifdef GL_GRAPHICS
	set_worker_done( true );
#endif

	return;

}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::finish_writing_files
///
/// @brief: write the last set of files
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::finish_writing_files() {

	using namespace files_paths;

	if ( mode == "pdbstats" ) output_pdbstats();
	if ( get_decoystats_flag() && ! get_ds_outpdbonly_flag() ) output_decoystats();
	if ( get_output_silent_gz_flag() ) {
		std::vector< std::string > file_list;
		std::string const basename( pdb_out_path + code + protein_name + ".out" );
		file_list.push_back( basename );
		file_list.push_back( basename + ".bonds" );
		file_list.push_back( basename + ".rot_templates");
		file_list.push_back( basename + ".full_torsions");
		for( int i = 0; i < int( file_list.size() ); ++i ) {
			std::string const filename( file_list[i] );
			utility::io::izstream in_stream ( filename );
			if ( in_stream ) {
				in_stream.close();
				in_stream.clear();
				std::cout << "GZIP SILENT FILE: " << filename << std::endl;
				utility::file::gzip( filename, true );
			} // file exists
		} // loop over each file
	}
	if ( get_output_scorefile_gz_flag() ) {
		std::string scorefile_name = get_scorefile_name();
		std::cout << "GZIP SCORE FILE: " << scorefile_name << std::endl;
		utility::file::gzip( scorefile_name, true );
	}

	return;
}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::finish_summary
///
/// @brief: write the final summary to stdout
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::finish_summary( std::ostream & iunit ) {

	using namespace files_paths;

	iunit << "======================================================" << std::endl;
	iunit << "DONE ::" << I( 6, std::max(curr_startnum-1,0) ) <<
		" starting structures built " << I( 9, std::max(curr_outnum-1,0) ) << " (nstruct) times" << std::endl;
	iunit << "This process generated " <<
		I( 6, num_decoys ) << " decoys from " << I( 7, attempted_decoys ) <<
		" attempts" << std::endl;
	if ( require_start ) iunit << space( 23 ) <<
												 I( 6, starting_pdbs_skipped ) << " starting pdbs were skipped" << std::endl;
	iunit << "======================================================" << std::endl;

	return;
}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::dump_status
///
/// @brief: dump contents of member data, eg. for debugging
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::dump_status( std::ostream & iunit ) {

	using namespace files_paths;

	iunit << "=== reporting job_distributor contents ===" << std::endl;
	iunit << "mode is: " << mode << std::endl;
	iunit << "nstartnum is: " << nstartnum << std::endl;
	iunit << "number_of_output is: " << number_of_output << std::endl;
	iunit << "curr_startnum is: " << curr_startnum << std::endl;
	iunit << "curr_outnum is: " << curr_outnum << std::endl;
	iunit << "num_decoys is: " << num_decoys << std::endl;
	iunit << "attempted_decoys is: " << attempted_decoys << std::endl;
	iunit << "starting_pdbs_skipped is: " << starting_pdbs_skipped << std::endl;
	iunit << "reinitialize is: " << reinitialize << std::endl;
	iunit << "=== done listing job_distributor contents ===" << std::endl;

	return;
}


//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::loop_over_function
///
/// @brief:  Loop over a function (passed in as first argument),
///          passing to it a ptr to a data class, then three loop counters
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 11/12/06
/////////////////////////////////////////////////////////////////////////////////
void job_distributor::loop_over_function(
	void(*function)(const void * data_class, int const i, int const j, int const k),
	void * data_class,
	int const imax,
	int const jmax,
	int const kmax
) {

	for ( int i=1; i<=imax; ++i ) {
		for ( int j=1; j<=jmax; ++j ) {
			for ( int k=1; k<=kmax; ++k ) {
				(*function)(data_class,i,j,k);
			}
		}
	}

	return;

}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::get_curr_outnum
///
/// @brief
/// read access to private data
///
/// @detailed
/// For MPI jobs, output PDBs must be named by either their
/// processor of origin, or which of the nstruct simulations they represent;
/// naming is handled automatically by the job_distributor if the default
/// output-pdb-creation path is followed; most design-mode protocols do not
/// folow this path (by setting files_paths::output_coord to false)
/// This function gives read access to the jd's which-of-the-nstruct private
/// data so that ambiguity over output names can be avoided for protocols that
/// control their own output behavior.
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: apl
///
/// @last_modified: 01/09/07
/////////////////////////////////////////////////////////////////////////////////
int job_distributor::get_curr_outnum() const {return curr_outnum;}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::get_mc_checkpoint_last_count
///
/// @brief
/// read access to private data
///
/// @detailed
/// The last pose Monte_carlo checkpoint step for the current nstruct has to be
/// saved and persist outside of the Monte_carlo object. This returns the last
/// Monte_carlo step that was checkpointed.
///
/// @authors: dekim
///
/// @last_modified: 04/11/07
//////////////////////////////////////////////////////////////////////////////
int job_distributor::get_mc_checkpoint_last_count(){
	return mc_checkpoint_last_count;
}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::set_mc_checkpoint_last_count
///
/// @brief
/// write access to private data
///
/// @detailed
/// Sets the last Pose Monte_carlo checkpoint step.
///
/// @authors: dekim
///
/// @last_modified: 04/11/07
//////////////////////////////////////////////////////////////////////////////
void job_distributor::set_mc_checkpoint_last_count(int last_cnt){
	mc_checkpoint_last_count = last_cnt;
}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::get_mc_checkpoint_current_count
///
/// @brief
/// read access to private data
///
/// @detailed
/// Gets the current Pose Monte_carlo checkpoint step.
///
/// @authors: dekim
///
/// @last_modified: 04/11/07
//////////////////////////////////////////////////////////////////////////////
int job_distributor::get_mc_checkpoint_current_count(){
	return mc_checkpoint_current_count;
}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::iterate_mc_checkpoint_current_count
///
/// @brief
/// iterate private data
///
/// @detailed
/// Iterates the current Pose Monte_carlo checkpoint step.
///
/// @authors: dekim
///
/// @last_modified: 04/11/07
//////////////////////////////////////////////////////////////////////////////
void job_distributor::iterate_mc_checkpoint_current_count(){
	mc_checkpoint_current_count++;
}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::get_mc_checkpoint_time
///
/// @brief
/// read access to private data
///
/// @detailed
/// Gets the time of the last Pose Monte_carlo checkpoint. This is used to
/// determine when to checkpoint given a time interval. The interval can be
/// set in seconds using the -checkpointing_interval argument.
///
/// @authors: dekim
///
/// @last_modified: 04/11/07
//////////////////////////////////////////////////////////////////////////////
time_t job_distributor::get_mc_checkpoint_time() const {
	return mc_checkpoint_time;
}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::update_mc_checkpoint_time
///
/// @brief
/// update checkpoint time private data
///
/// @detailed
/// Updates the time of the last Pose Monte_carlo checkpoint.
///
/// @authors: dekim
///
/// @last_modified: 04/11/07
//////////////////////////////////////////////////////////////////////////////
void job_distributor::update_mc_checkpoint_time() {
	time(&mc_checkpoint_time);
}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::get_mc_checkpoint_file
///
/// @brief
/// read access to private data
///
/// @detailed
/// Provides the checkpoint filename to the pose Monte_carlo object.
///
/// @authors: dekim
///
/// @last_modified: 04/11/07
//////////////////////////////////////////////////////////////////////////////
std::string job_distributor::get_mc_checkpoint_file() const {
	return mc_checkpoint_file;
}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::reset_mc_checkpoint
///
/// @brief
/// clear checkpoint information
///
/// @detailed
/// Clears the pose Monte_carlo checkpoint information. This should be called
/// after each nstruct if checkpointing is used.
///
/// @authors: dekim
///
/// @last_modified: 04/11/07
//////////////////////////////////////////////////////////////////////////////
void job_distributor::reset_mc_checkpoint() {

	if (!get_do_pose_checkpointing()) return;

#ifdef BOINC
	boinc_begin_critical_section();
#endif

	utility::io::ozstream out_stream( mc_checkpoint_file );
	if (!out_stream.good() ) {
		std::cout << "STOP: cant open file: " << mc_checkpoint_file << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}
	out_stream << "0 0 0";
	out_stream.close();
	out_stream.clear();

	mc_checkpoint_last_count = 0;
	mc_checkpoint_current_count = 0;
	update_mc_checkpoint_time();

#ifdef BOINC
  boinc_end_critical_section();
#endif
}

//////////////////////////////////////////////////////////////////////////////
/// @begin job_distributor::skip_to_mc_checkpoint
///
/// @brief
/// skip checkpoint step
///
/// @detailed
/// Determine whether to skip to get to the last step that was checkpointed.
///
/// @authors: dekim
///
/// @last_modified: 04/11/07
//////////////////////////////////////////////////////////////////////////////
bool job_distributor::skip_to_mc_checkpoint() {
	return ( mc_checkpoint_current_count < mc_checkpoint_last_count );
}

#ifdef BOINC

////////////////////////////////
///  BOINC_job_distributor classes
////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::type
///
/// @brief:  Return the job_distributor type
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 11/12/06
/////////////////////////////////////////////////////////////////////////////////
std::string BOINC_job_distributor::type() {

	return "BOINC";

}

//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::initialize
///
/// @brief: setup the BOINC job distributor
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void BOINC_job_distributor::initialize()
{

	// Reset namespace vars as appropriate for BOINC
	set_BOINC_defaults();

	// read BOINC params
	read_BOINC_params();

	// do the initializations from the base class
	job_distributor::initialize();

	// member data of the BOINC_job_distributor class
	farlx_stage=0;

	// set the copies in the boinc namespace
	boinc_params::orig_number_of_output = number_of_output;
	boinc_params::orig_nstartnm = nstartnum;

	// restore decoy progress for BOINC
	restoreDecoyInfo( attempted_decoys, num_decoys, farlx_stage );

	if (num_decoys > 0) {
		update_boinc_params( num_decoys, number_of_output );
		//chu bugfix 2006-10-05
		//not only update number_of_output, but also temporarily for curr_outnum
		//otherwise, jd.get_next_job_num() will get confused. The correct
		//curr_outnum will be set there.
		curr_outnum = number_of_output;

		if (num_decoys >= number_of_output || boinc_params::pct_complete >= 1) {
			curr_startnum = nstartnum+1;
			curr_outnum++;
			finish();
		}
	} else {
		// boinc_params::pct_complete = 0.01;
		boinc_params::pct_complete = 0.00; //Start counting at 0.0!
	}
	boinc_fraction_done(boinc_params::pct_complete);
	std::cout << "BOINC :: " << utility::timestamp() << " :: mode: " << mode << " :: nstartnum: " <<
		nstartnum << " :: number_of_output: " << number_of_output <<
		" :: num_decoys: " << num_decoys << " :: pct_complete: " << boinc_params::pct_complete << std::endl;

	if (!boinc_is_standalone()) {
		// check number of restarts with no progress (if app is switching before decoys are made)
		checkNoProgressInitCount( attempted_decoys, farlx_stage );
		if ( boinc_params::no_progress_init_cnt > BOINC_MAX_NO_PROGRESS_INIT_CNT ) {
			std::cerr << "Too many restarts with no progress. Keep application in memory while preempted." << std::endl;
			curr_startnum = nstartnum+1;
			curr_outnum++;
			finish();
		}
	}
	rosetta_has_started();
	watchdog_start();

	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::checkNoProgressInitCount
///
/// @brief: setup the BOINC job distributor
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
bool BOINC_job_distributor::checkNoProgressInitCount( int& attempts, int& farlxstage )
{

	using namespace boinc_params;
  bool saved = false;
	int prev_initcnt = 0;
	int prev_attemptcnt = 0;
	int prev_farlxstage = 0;
	int prev_mc_checkpoint_count = 0;

	// read last attempts count from init_cnt file
	utility::io::izstream initcnt_istream( init_cnt_filename );
	if ( initcnt_istream ) {
		initcnt_istream >> prev_initcnt >> prev_attemptcnt >> prev_farlxstage >>
				prev_mc_checkpoint_count;
		initcnt_istream.close();
		initcnt_istream.clear();
	}
	// get last mc_checkpoint
	if (get_do_pose_checkpointing()) {
		utility::io::izstream in_stream( mc_checkpoint_file );
		if (in_stream) {
			int cnt, isfullatom, ntrials;
			in_stream >> cnt >> isfullatom >> ntrials;
			in_stream.close();
			in_stream.clear();
			mc_checkpoint_last_count = cnt;
		}
	}
	// iterate no_progress_init_cnt if the following have not increased
	no_progress_init_cnt = (
			prev_attemptcnt < attempts ||													// nstuct attempts
			prev_farlxstage < farlxstage ||												// farlx_stage checkpoints
			prev_mc_checkpoint_count < mc_checkpoint_last_count		// mc_checkpoints
	) ? prev_initcnt : prev_initcnt + 1;

	// update counts
  utility::io::ozstream initcnt_ostream( init_cnt_filename );
  if ( initcnt_ostream ) {
    initcnt_ostream << no_progress_init_cnt << " " << attempts << " " <<
				farlxstage << " " << mc_checkpoint_last_count;
    initcnt_ostream.close();
    initcnt_ostream.clear();
    saved = true;
  }

  return saved;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::call_abrelax
///
/// @brief: essentially a stub, required due to a derived class
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void BOINC_job_distributor::call_abrelax() {

	abrelax(farlx_stage);
	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::run_main_pose1
///
/// @brief: setup when mode == pose1
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void BOINC_job_distributor::run_main_pose1() {

	main_pose1( number_of_output, nstartnum );
	restoreDecoyInfo( attempted_decoys, num_decoys, farlx_stage );
	update_boinc_params( num_decoys, number_of_output );
	curr_startnum=nstartnum+1;
	curr_outnum++;
	attempted_decoys = num_decoys;

	finish();

	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::set_BOINC_defaults
///
/// @brief: reset namespace vars as appropriate for BOINC
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void BOINC_job_distributor::set_BOINC_defaults() {

	files_paths::default_nstruct = 10;

	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::read_BOINC_params
///
/// @brief: read BOINC params from stdout.txt
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void BOINC_job_distributor::read_BOINC_params() {

	// allow upload of stdout.txt
	char filename[256];
	boinc_resolve_filename("stdout.txt", filename, 256);
	freopen(filename, "a", stdout);
	realafteroption( "cpu_frac", boinc_project_prefs::default_max_cpu , boinc_project_prefs::default_max_cpu );
	realafteroption( "frame_rate", boinc_project_prefs::default_max_fps , boinc_project_prefs::default_max_fps );
	intafteroption( "cpu_run_time", boinc_project_prefs::default_cpu_run_time, boinc_project_prefs::default_cpu_run_time );
	files_paths::mode_title = "Initializing";

	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::made_output_decoy
///
/// @brief: called when an output decoy has been made
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void BOINC_job_distributor::made_output_decoy() {

	job_distributor::made_output_decoy();

	// From WCG
	// Don't checkpoint if a prior PDB exists (i.e. a checkpoint file.)
	// otherwise we'll be overwriting the previous checkpoint with invalid data.
	// Checkpointing and timer update has been moved to boinc/boinc_rosetta_util.cc.
	// This allows jumping to also call the same routines.
	if ( boinc_checkpoint_in_main_loop(
				 attempted_decoys,num_decoys,number_of_output,farlx_stage) ) {
		curr_startnum = nstartnum+1;
		curr_outnum++;
		finish();
	}

	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::report_progress
///
/// @brief: report how far we've come
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void BOINC_job_distributor::report_progress() {

	std::cout << "BOINC :: " << utility::timestamp()
						<< " :: report_progress() :: num_decoys: " << num_decoys
						<< " :: number_of_output: " << number_of_output
						<< " :: pct_complete: " << boinc_params::pct_complete
						<< std::endl;
	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::finish
///
/// @brief: end Rosetta (cleanly)
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void BOINC_job_distributor::finish() {

	using namespace files_paths;

	delete_START();
	prof::show();

	boinc_begin_critical_section();
	boinc_fraction_done(1);
	if ( runlevel_ns::runlevel == runlevel_ns::silent &&
		!utility::file::file_exists( pdb_out_path + code + protein_name + ".out" ) &&
		!utility::file::file_exists( pdb_out_path + code + protein_name + ".out.gz" )  ) {
			utility::file::create_blank_file( pdb_out_path + code + protein_name + ".out" );
	}

	finish_writing_files();

	boinc_end_critical_section();

	finish_summary( std::cerr );
	finish_summary( std::cout );

	repeat_final(); // added by sheffler 11/5/04
	decoy_features_ns::decoy_features_final();

#ifdef GL_GRAPHICS
	set_worker_done( true );
#endif

	std::cerr << std::endl << std::endl;
	std::cerr << "BOINC :: Watchdog shutting down..." << std::endl;
	rosetta_has_finished();
	watchdog_finish();
	std::cerr << "BOINC :: BOINC support services shutting down..." << std::endl;
	boinc_finish(0); // After this function is called nothing executes beyond it.

	return;
}

//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::finish_summary
///
/// @brief: write the final summary to stdout
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void BOINC_job_distributor::finish_summary( std::ostream & iunit ) {

	using namespace files_paths;

	double cputime = 0.0;
	boinc_wu_cpu_time(cputime); // get wu cpu run time

	iunit << "======================================================" << std::endl;
	iunit << "DONE ::" << I( 6, std::max(curr_startnum-1,0) ) <<
		" starting structures " << I( 8, cputime ) << " cpu seconds" << std::endl;
	iunit << "This process generated " <<
		I( 6, num_decoys ) << " decoys from " << I( 7, attempted_decoys ) <<
		" attempts" << std::endl;
	if ( require_start ) iunit << space( 23 ) <<
												I( 6, starting_pdbs_skipped ) << " starting pdbs were skipped" << std::endl;
	iunit << "======================================================" << std::endl;

	return;
}

//////////////////////////////////////////////////////////////////////////////
/// @begin BOINC_job_distributor::get_next_job_num
///
/// @brief: set the start_inx and output_inx for next job. Returns true if
///          more jobs remain, returns false if we've already done the last job
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks: output_inx for next job in BOINC is determined by num_decoys,
///           which is either read from rosetta_decoy_cnt.txt or ZERO
///
/// @references
///
/// @authors: chu
///
/// @last_modified: 10/05/06
/////////////////////////////////////////////////////////////////////////////////
bool BOINC_job_distributor::get_next_job_num() {

	++curr_outnum;
	if ( curr_outnum > number_of_output ) {
		++curr_startnum;
		if ( curr_startnum > nstartnum ) return false;
		curr_outnum=num_decoys+1;
	}

	return true;

}

#endif // BOINC






#ifdef USEMPI

////////////////////////////////
///  MPI_generic_job_distributor classes
////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_generic_job_distributor::type
///
/// @brief:  Return the job_distributor type
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 11/12/06
/////////////////////////////////////////////////////////////////////////////////
std::string MPI_generic_job_distributor::type() {

	return "MPI_generic";

}

//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_generic_job_distributor::initialize
///
/// @brief: setup the job distributor
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void MPI_generic_job_distributor::initialize() {

	using namespace basic::options;

	cpu_num=0;
	is_server=false;

	// Initialize MPI
	std::cout << "about to initialize MPI" << std::endl;
	MPI::Init(arg_count,arg_vector);

	cpu_num = MPI::COMM_WORLD.Get_rank();

	// redirect output to a separate file for each process
	char outFile[100];
	sprintf(outFile, "rosetta.mpi.out%i",cpu_num);
	freopen(outFile, "w",stdout);

	// Reset namespace vars as appropriate for MPI
	set_MPI_defaults();

	if ( cpu_num == 0 ) {
		is_server=true;
		std::cout << "this node is the server" << std::endl;
	} else {
		std::cout << "client number is: " << cpu_num << std::endl;
	}

	// Invoke the base class initialization
	job_distributor::initialize();

	sync_random_seeds();

	return;

}

//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_generic_job_distributor::sync_random_seed
///
/// @brief: reset the random seed on all cpus to match
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void MPI_generic_job_distributor::sync_random_seeds() {

	// jk distribute the random seed from the main_node
	// for MPI_generic_job_distributor, clients take the seed of the main node
	int jran = -get_server_seed(); // ran3 takes a negative seed
	ran3_seed( jran );

	return;

}



//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_generic_job_distributor::get_server_seed
///
/// @brief: use a requested seed or get the random seed from the clock time
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
int MPI_generic_job_distributor::get_server_seed() {

	int server_seed;

	if ( is_server ) {
		FArray1D_int timearray( 3 );
		if ( truefalseoption("constant_seed") ) {
			intafteroption("jran",1111111,server_seed);
		} else {
			itime(timearray);
			server_seed = 1 + 20*(3600*timearray(1)+60*timearray(2)+timearray(3));
		}
	}

	MPI::COMM_WORLD.Bcast(&server_seed,1,MPI::INT,0);

	return server_seed;

}


////////////////////////////////////////////////////////////////////////////////
/// @begin MPI_generic_job_distributor::set_MPI_defaults
///
/// @brief: reset namespace vars as appropriate for MPI
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors jk Aug. 11/06
///
/// @last_modified
/////////////////////////////////////////////////////////////////////////////////
void MPI_generic_job_distributor::set_MPI_defaults() {

	// Since the server will control running jobs,
	// thus preventing multiple cpus from running the same thing,
	// only create files AFTER running the desired simulation rather than before
	// (this prevents "missing" output files upon restarting)
	files_paths::touch_before_starting = false;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_generic_job_distributor::run_server
///
/// @brief: run a server to distribute jobs
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void MPI_generic_job_distributor::run_server() {

	// server waits for requests, then answers them using "get_next_task"
	int const nprocs = MPI::COMM_WORLD.Get_size();
	int num_available_clients = nprocs-1;
	int outgoing_task[NUM_MPI_ARGS];
	int tag=1;
	std::cout << "server running with " << num_available_clients << " clients" << std::endl;
	while ( num_available_clients > 0 ) {

		// wait for a client to request a job
		MPI::Status status;
		int source=MPI::ANY_SOURCE;
		int msg;
		MPI::COMM_WORLD.Recv(&msg, 1, MPI::INT, source, tag, status);

		// reply to this client, using results from a call to get_next_task
		int dest=status.Get_source();
		// jk Note: calling "get_next_task()" sets jd_i, jd_j, jd_k
		if ( get_next_task() ) {
			outgoing_task[0]=1;
			outgoing_task[1]=jd_i;
			outgoing_task[2]=jd_j;
			outgoing_task[3]=jd_k;
			MPI::COMM_WORLD.Send(outgoing_task, NUM_MPI_ARGS, MPI::INT, dest, tag);
		} else {
			std::cout << "done with client " << dest << std::endl;
			--num_available_clients; // this client won't run anything else
		}

	}

	// tell each client to finish
	// Note: wait until all clients are done before doing this, so that
	// all clients perform a syncronized exit from control of this server
	for ( int dest=1; dest<=(nprocs-1); ++dest ) {
		outgoing_task[0]=0;
		outgoing_task[1]=0;
		outgoing_task[2]=0;
		outgoing_task[3]=0;
		MPI::COMM_WORLD.Send(outgoing_task, NUM_MPI_ARGS, MPI::INT, dest, tag);
	}

	std::cout << "ending server" << std::endl;
	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_generic_job_distributor::get_next_job_num
///
/// @brief: set the start_inx and output_inx for next job. Returns true if
///          more jobs remain, returns false if we've already done the last job
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
bool MPI_generic_job_distributor::get_next_job_num() {

	sync_random_seeds();
	return job_distributor::get_next_job_num();

}


//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_generic_job_distributor::get_next_task
///
/// @brief: set the start_inx and output_inx for next job. Returns true if
///          more jobs remain, returns false if we've already done the last job
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
bool MPI_generic_job_distributor::get_next_task() {

	if ( is_server ) {
		++jd_k;
		if ( jd_k > jd_kmax ) {
			jd_k=1;
			++jd_j;
			if ( jd_j > jd_jmax ) {
				jd_j=1;
				++jd_i;
			}
		}
		if ( jd_i > jd_imax ) return false;
		return true;

	} else {

		int tag=1;

		// tell server (node 0) that we're ready for a task by sending "1"
		int dest=0;
		int count=0;
		int msg=1;
		std::cout << "requesting task from server" << std::endl;
		MPI::COMM_WORLD.Send(&msg, count, MPI::INT, dest, tag);

		// wait for a server response, use this to set curr_startnum and curr_outnum
		MPI::Status status;
		int incoming_task[NUM_MPI_ARGS];
		int source=0;
		MPI::COMM_WORLD.Recv(incoming_task, NUM_MPI_ARGS, MPI::INT, source, tag, status);

		// the first incoming int marks "success", the others are jd_i, jd_j, jd_k
		if ( incoming_task[0] == 1 ) {
			jd_i=incoming_task[1];
			jd_j=incoming_task[2];
			jd_k=incoming_task[3];
			return true;
		} else {
			return false;
		}
	}

	return false; // never used

}

//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_generic_job_distributor::finish
///
/// @brief: end Rosetta (cleanly)
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void MPI_generic_job_distributor::finish() {

	std::cout << "closing MPI" << std::endl;
	MPI::Finalize();
	std::cout << "exiting Rosetta" << std::endl;
	job_distributor::finish();
	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_generic_job_distributor::loop_over_function
///
/// @brief:  Use MPI to loop over a function (passed in as first argument),
///          passing to it a ptr to a data class, then three loop counters
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 11/12/06
/////////////////////////////////////////////////////////////////////////////////
void MPI_generic_job_distributor::loop_over_function(
	void(*function)(const void * data_class, int const i, int const j, int const k),
	void * data_class,
	int const imax,
	int const jmax,
	int const kmax
) {

	// Use MPI to distribute small "tasks"

	// Note: i is outermost loop, k is innermost loop
	jd_i = 0;
	jd_j = jmax;
	jd_k = kmax;
	jd_imax = imax;
	jd_jmax = jmax;
	jd_kmax = kmax;

	if ( is_server ) {
		run_server();
	} else {
		while ( get_next_task() ) {
			(*function)(data_class,jd_i,jd_j,jd_k);
		}
	}

	sync_random_seeds();

	return;

}



////////////////////////////////
///  MPI_full_job_distributor classes
////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_full_job_distributor::type
///
/// @brief:  Return the job_distributor type
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 11/12/06
/////////////////////////////////////////////////////////////////////////////////
std::string MPI_full_job_distributor::type() {

	return "MPI_full";

}


//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_full_job_distributor::initialize
///
/// @brief: setup the job distributor
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void MPI_full_job_distributor::initialize() {

	// Invoke the base class initialization
	MPI_generic_job_distributor::initialize();

	// jk distribute the random seed from the main_node
	// for MPI_full_job_distributor, clients offset the random number with
	// their rank, so that each client has a unique seed.
	int jran = -(get_server_seed() + cpu_num); // ran3 takes a negative seed
	ran3_seed( jran );

	// clients will leave this initialization function and continue running Rosetta
	// server will enter the "run_server" function
	if ( is_server ) {
		run_server();
		finish();
	}

	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_full_job_distributor::prepare_for_next_startnum
///
/// @brief: in the base class, sets the start_inx and output_inx
///      such that next time we'll move on to a new start_num. When
///      using MPI, however, we can't really do this (because more jobs
///      may have been started. Instead, just do nothing.
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
void MPI_full_job_distributor::prepare_for_next_startnum() {

	return;

}


//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_full_job_distributor::get_next_task
///
/// @brief: set the start_inx and output_inx for next job. Returns true if
///          more jobs remain, returns false if we've already done the last job
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
bool MPI_full_job_distributor::get_next_task() {

	if ( is_server ) {
		bool const return_val = get_next_job_num();
		jd_i = curr_startnum;
		jd_j = curr_outnum;
		return return_val;
	} else {
		std::cout << "ERROR: MPI_full_job_distributor called get_next_task()." << std::endl;
		std::cout << "This should never happen...." << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		return false;
	}

	return false; // never used

}


//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_full_job_distributor::get_next_job_num
///
/// @brief: set the start_inx and output_inx for next job. Returns true if
///          more jobs remain, returns false if we've already done the last job
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 06/24/06
/////////////////////////////////////////////////////////////////////////////////
bool MPI_full_job_distributor::get_next_job_num() {

	if ( is_server ) {
		return job_distributor::get_next_job_num();
	} else {
		bool const return_val = MPI_generic_job_distributor::get_next_task();
		curr_startnum=jd_i;
		curr_outnum=jd_j;
		return return_val;
	}

	return false; // never used

}


//////////////////////////////////////////////////////////////////////////////
/// @begin MPI_full_job_distributor::loop_over_function
///
/// @brief:  Loop over a function (passed in as first argument),
///          passing to it a ptr to a data class, then three loop counters
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors: jk
///
/// @last_modified: 11/12/06
/////////////////////////////////////////////////////////////////////////////////
void MPI_full_job_distributor::loop_over_function(
	void(*function)(const void * data_class, int const i, int const j, int const k),
	void * data_class,
	int const imax,
	int const jmax,
	int const kmax
) {

	// Call the NON-MPI version of this function, since we're already using
	// MPI to distribute starting jobs (so we can't use it to distribute small "tasks")

	job_distributor::loop_over_function( function, data_class, imax, jmax, kmax);
	return;

}


#endif // MPI

