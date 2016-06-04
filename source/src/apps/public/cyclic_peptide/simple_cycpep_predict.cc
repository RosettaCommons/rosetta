// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/public/cyclic_peptide/simple_cycpep_predict.cc
/// @brief Predicts structure of N-to-C cyclic peptides from amino acid sequence.
/// @details
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifdef BOINC
#include <utility/boinc/boinc_util.hh>
#include <protocols/boinc/boinc.hh>
#include "boinc_zip.h"
#endif // BOINC -- keep this first.

// MPI headers
#ifdef USEMPI
#include <mpi.h> //Message Passing Interface for parallelization -- keep this second
#endif

//General includes
#include <basic/options/option.hh>
#include <basic/options/keys/cyclic_peptide.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>
#include <numeric/conversions.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <numeric/constants.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/protein_interface_design/filters/HbondsToResidueFilter.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/Energies.hh>

//Application-specific includes
#ifdef USEMPI
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI.hh>
#endif
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh>

#include <utility/vector1.hh>
#include <stdio.h>

//Tracer:
static THREAD_LOCAL basic::Tracer TR( "apps.public.cyclic_peptide.simple_cycpep_predict" );

/*************************/
/* FUNCTION PROTOTYPES:  */
/*************************/

#ifdef USEMPI
/// @brief If rank > 0, wait for a message from proc 0 to continue.  If rank
/// is zero, send the message.
/// #details Only used in MPI mode.
void wait_for_proc_zero();
#endif

#ifdef USEMPI
/// @brief In MPI mode, set the MPI-specific variables that the parallel version
/// of the protocol will need.
void
set_MPI_vars(
	int &MPI_rank,
	int &MPI_n_procs,
	core::Size &total_hierarchy_levels,
	utility::vector1 < core::Size > &procs_per_hierarchy_level,
	utility::vector1< core::Size > &batchsize_by_level,
	std::string &sort_by,
	bool &select_highest,
	core::Real &output_fraction,
	std::string &output_filename
);
#endif

/*************************/
/* FUNCTION DEFINITIONS: */
/*************************/

#ifdef USEMPI
/// @brief If rank > 0, wait for a message from proc 0 to continue.  If rank
/// is zero, send the message.
/// #details Only used in MPI mode.
void
wait_for_proc_zero() {
	char mybuf('A'); //A dummy piece of information to send.
	MPI_Bcast( &mybuf, 1, MPI_CHAR, 0, MPI_COMM_WORLD );
}
#endif

#ifdef USEMPI
/// @brief In MPI mode, set the MPI-specific variables that the parallel version
/// of the protocol will need.
void
set_MPI_vars(
	int &MPI_rank,
	int &MPI_n_procs,
	core::Size &total_hierarchy_levels,
	utility::vector1 < core::Size > &procs_per_hierarchy_level,
	utility::vector1< core::Size > &batchsize_by_level,
	std::string &sort_by,
	bool &select_highest,
	core::Real &output_fraction,
	std::string &output_filename
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string const errormsg( "Error in simple_cycpep_predict application: ");

	//Get the rank and number of processes.
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &MPI_n_procs );

	if (TR.Debug.visible()) {
		TR.Debug << "MPI rank " << MPI_rank << " of " << MPI_n_procs << " reporting for duty." << std::endl;
	}

	if(MPI_rank == 0 ) {
		runtime_assert_string_msg( option[cyclic_peptide::MPI_processes_by_level].user(), errormsg + "In MPI mode, the user MUST specify the number of communication levels, and the number of processes in each level, with the \"-cyclic_peptide:MPI_processes_by_level\" flag." );
	}
	utility::vector1< long int > user_specified_hierarchy_levels( option[cyclic_peptide::MPI_processes_by_level]() );
	if(MPI_rank == 0 ) {
		runtime_assert_string_msg( user_specified_hierarchy_levels.size() > 1, errormsg + "The number of communication levels specified with the \"-cyclic_peptide:MPI_processes_by_level\" flag must be greater than one." );
		runtime_assert_string_msg( user_specified_hierarchy_levels.size() <= static_cast<core::Size>(MPI_n_procs), errormsg + "The number of communication levels specified with the \"-cyclic_peptide:MPI_processes_by_level\" flag must be less than or equal to the number of processes launched." );
		runtime_assert_string_msg( user_specified_hierarchy_levels[1] == 1, errormsg + "The first communication level must have a single, master process.");
		core::Size proccount(0);
		for(core::Size i=1; i<=user_specified_hierarchy_levels.size(); ++i) {
			runtime_assert_string_msg( user_specified_hierarchy_levels[i] > 0, errormsg + "Each level of communication specified with the \"-cyclic_peptide:MPI_processes_by_level\" flag must have at least one process associated with it." );
			proccount += static_cast<core::Size>( user_specified_hierarchy_levels[i] );
			if(i>1) {
				runtime_assert_string_msg( user_specified_hierarchy_levels[i] >= user_specified_hierarchy_levels[i-1], errormsg + "Each level of communication specified with the \"-cyclic_peptide:MPI_processes_by_level\" flag must have an equal or greater number of processes as compared to its parent." );
			}
		}
		runtime_assert_string_msg( proccount == static_cast<core::Size>(MPI_n_procs), errormsg + "The total number of processes specified with the \"-cyclic_peptide:MPI_processes_by_level\" flag must match the total number of processes launched." );
	}
	total_hierarchy_levels = static_cast< core::Size >( user_specified_hierarchy_levels.size() );
	procs_per_hierarchy_level.clear();
	procs_per_hierarchy_level.reserve( total_hierarchy_levels );
	for(core::Size i=1; i<=total_hierarchy_levels; ++i) {
		procs_per_hierarchy_level.push_back( static_cast< core::Size >( user_specified_hierarchy_levels[i] ) );
	}

	if(MPI_rank == 0 ) {
		runtime_assert_string_msg( option[cyclic_peptide::MPI_batchsize_by_level].user(), errormsg + "Batch sizes for each level must be specified." );
		runtime_assert_string_msg( option[cyclic_peptide::MPI_batchsize_by_level]().size() == total_hierarchy_levels - 1, errormsg + "The number of values provided with -cyclic_peptide::MPI_batchsize_by_level must be one less than the number of communication levels." );
		for(core::Size i=1; i<total_hierarchy_levels; ++i) {
			runtime_assert_string_msg( option[cyclic_peptide::MPI_batchsize_by_level]()[i] > 0, errormsg + "The batch size must be greater than zero." );
			if(i>1) {
				runtime_assert_string_msg( option[cyclic_peptide::MPI_batchsize_by_level]()[i] <= option[cyclic_peptide::MPI_batchsize_by_level]()[i-1], errormsg + "The lower level batch sizes must be smaller than the upper level batch sizes." );
			}
		}
	}
	batchsize_by_level.clear();
	batchsize_by_level.reserve( total_hierarchy_levels - 1 );
	for(core::Size i=1; i<total_hierarchy_levels; ++i) {
		batchsize_by_level.push_back ( static_cast< core::Size >( option[cyclic_peptide::MPI_batchsize_by_level]()[i] ) );
	}

	sort_by = option[cyclic_peptide::MPI_sort_by]();
	select_highest = option[cyclic_peptide::MPI_choose_highest]();
	output_fraction = static_cast<core::Real>( option[cyclic_peptide::MPI_output_fraction]() );

	if( MPI_rank == 0) {
		runtime_assert_string_msg( !option[out::file::o].user(), errormsg + "The -out:file:o option cannot be used with the simple_cycpep_predict app in MPI mode.  Only silent file output is permitted." );
		runtime_assert_string_msg( option[out::file::silent].user(), errormsg + "A silent file for output must be specified for the simple_cycpep_predict app in MPI mode." );
	}
	output_filename = option[out::file::silent]();

	wait_for_proc_zero();
}
#endif


/*************************/
/* MAIN: */
/*************************/
int
main( int argc, char * argv [] )
{
	try {
		protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::register_options(); //Has to be before init().
		devel::init(argc, argv);

#ifdef USEMPI
		//Variables for MPI -- only used in MPI mode:
		int MPI_rank(0); //Note -- must be int for MPI
		int MPI_n_procs(0);
		utility::vector1< core::Size > batchsize_per_level;
		//core::Size hierarchy_level(0); //0 for emperor process, 1 for first-level masters, 2 for second-level masters ... total_hierarchy_levels - 1 for slaves/workers.
		core::Size total_hierarchy_levels(0); //Note that for a 3-level hierarchy, this will be 3, but the highest hierarchy_level index will be 2.
		utility::vector1< core::Size > procs_per_hierarchy_level; //The number of processes at each level of the communications hierarchy.  Should be 1, n>=1, m>=n, etc.
		std::string sort_by("");
		bool select_highest(false);
		core::Real output_fraction(1.0);
		std::string output_filename("");
		set_MPI_vars( MPI_rank, MPI_n_procs, total_hierarchy_levels, procs_per_hierarchy_level, batchsize_per_level, sort_by, select_highest, output_fraction, output_filename ); //Get the values of these vars (only used in MPI mode).
#endif

		if ( TR.visible() ) {
#ifdef USEMPI
			if( MPI_rank == 0 ) { //If this is the emperor:
				TR << "Starting simple_cycpep_predict.cc, MPI version." << std::endl;
				TR << "Application created 16 September 2015 by Vikram K. Mulligan, Ph.D., Baker laboratory." << std::endl;
				TR << "For questions, contact vmullig@uw.edu." << std::endl << std::endl;
				TR << "Launching " << MPI_n_procs << " parallel processes with " << total_hierarchy_levels + 1 << " levels of communication." << std::endl;
				TR << "Batch size is ";
				for(core::Size i=1, imax=batchsize_per_level.size(); i<=imax; ++i) {
					TR << batchsize_per_level << " (level " << i << ")";
					if(i<imax) TR << ", ";
				}
				TR << "." << std::endl;
				TR.flush();
			}
#else
			TR << "Starting simple_cycpep_predict.cc" << std::endl;
			TR << "Application created 16 September 2015 by Vikram K. Mulligan, Ph.D., Baker laboratory." << std::endl;
			TR << "For questions, contact vmullig@uw.edu." << std::endl << std::endl;
			TR.flush();
#endif
		}

#ifdef USEMPI
		protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication_MPI this_app( MPI_rank, MPI_n_procs, total_hierarchy_levels, procs_per_hierarchy_level, batchsize_per_level, sort_by, select_highest, output_fraction, output_filename );
#else
		protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication this_app;
#endif
		this_app.run();

	} catch ( utility::excn::EXCN_Base& excn ) {
		TR.Error << "Exception caught: " << std::endl;
		excn.show( TR.Error );
		TR.Error.flush();
		return -1;
	}

	if ( TR.visible() ) {
		TR << "Finished simple_cycpep_predict.cc.  Exiting." << std::endl;
		TR.flush();
	}

#ifdef USEMPI
	MPI_Barrier(MPI_COMM_WORLD); //Wait until all procs reach this point.
	MPI_Finalize(); //MPI finalization and cleanup.
#endif

	return 0;
}
