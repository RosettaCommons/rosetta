// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/helical_bundle/helical_bundle_predict.cc
/// @brief An application that attempts to predict helical bundle structures, without using PDB fragments.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// MPI headers
#ifdef USEMPI
#include <mpi.h> //Message Passing Interface for parallelization -- keep this first
#endif

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/helical_bundle_predict/HelicalBundlePredictApplication.hh>
#ifdef USEMPI
#include <protocols/helical_bundle_predict/HelicalBundlePredictApplication_MPI.hh>
#include <protocols/cyclic_peptide_predict/util.hh>
#endif //USEMPI


// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/pointer/memory.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/helical_bundle_predict.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR("apps.public.helical_bundle.helical_bundle_predict");


/// @brief Indicate what options are relevant for this application.
void register_options() {
	protocols::helical_bundle_predict::HelicalBundlePredictApplicationOptions::register_options(); //Static function call.
}


/// @brief Set relevant options.
/// @details Triggers reads from disk!
void
get_options_from_options_collection(
	protocols::helical_bundle_predict::HelicalBundlePredictApplicationOptions & options,
	utility::options::OptionCollection const & option_collection
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	static std::string const errmsg( "Error in helical_bundle_predict application: " );
	runtime_assert_string_msg( option_collection[in::file::fasta].user() || option_collection[helical_bundle_predict::sequence_file].user(), errmsg + "The user must provide a sequence in FASTA format or as a whitespace-separated series of residue type names, using the \"-in:file:fasta\" commandline flag for FASTA input or the \"helical_bundle_predict:sequence_file\" commandline for full name input." );
	runtime_assert_string_msg( !( option_collection[in::file::fasta].user() && option_collection[helical_bundle_predict::sequence_file].user() ), errmsg + "The user must provide EITHER the \"-in:file:fasta\" flag OR the \"-helical_bundle_predict:sequence_file\" flag, not both." );
	if ( option_collection[in::file::fasta].user() ) {
		runtime_assert_string_msg( option_collection[in::file::fasta]().size() == 1, errmsg + "The user must provide one and only one sequence in FASTA format using the \"-in:file:fasta\" commandline flag." );
	}
	runtime_assert_string_msg( option_collection[helical_bundle_predict::helix_assignment_file].user(), errmsg + "The user must provide helix assignments using the \"-helical_bundle_predict:helix_assignment_file\" commandline flag." );

	if ( option_collection[in::file::fasta].user() ) {
		options.set_fasta_file( option_collection[in::file::fasta]()[1] );
	} else {
		runtime_assert( option_collection[helical_bundle_predict::sequence_file].user() ); //Should be true at this point; checked above.
		options.set_sequence_file( option_collection[helical_bundle_predict::sequence_file]() );
	}
	options.set_helix_assignment_file( option_collection[helical_bundle_predict::helix_assignment_file]() );

	// The following triggers read from disk.  Modify tihs for multi-threaded or multi-process processing:
	options.read_inputs();
}

#ifndef USEMPI
/// @brief In the non-MPI build, set the file output options.
void
set_file_output(
	protocols::helical_bundle_predict::HelicalBundlePredictApplication & application
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[out::file::silent].user() ) {
		application.set_output_prefix_and_suffix( option[out::prefix].user() ? option[ out::prefix ]() : "result" , option[out::suffix]() );
		application.set_silent_output();
		TR << "Configuring for Rosetta silent file output." << std::endl;
	} else {
		application.set_output_prefix_and_suffix( option[out::prefix].user() ? option[ out::prefix ]() : "result_" , option[out::suffix]() );
		if ( option[out::mmCIF].value() ) {
			application.set_output_format( core::import_pose::CIF_file );
			TR << "Configuring for CIF output." << std::endl;
		} else if ( option[out::mmtf].value() ) {
			application.set_output_format( core::import_pose::MMTF_file );
			TR << "Configuring for MMTF output." << std::endl;
		} else {
			application.set_output_format( core::import_pose::PDB_file );
			TR << "Configuring for PDB output." << std::endl;
		}
	}
}
#endif

/// @brief Entry point for program execution.
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::helical_bundle_predict;

		register_options();
		devel::init( argc, argv );

		HelicalBundlePredictApplicationOptionsOP options( utility::pointer::make_shared< HelicalBundlePredictApplicationOptions >() );
		get_options_from_options_collection( *options, basic::options::option ); //Read from global options, here.  (This is in the main() function, so global options aren't that risky here.)


#ifdef USEMPI
		//Variables for MPI -- only used in MPI mode:
		int MPI_rank(0); //Note -- must be int for MPI
		int MPI_n_procs(0);
		utility::vector1< core::Size > batchsize_per_level;
		core::Size total_hierarchy_levels(0); //Note that for a 3-level hierarchy, this will be 3, but the highest hierarchy_level index will be 2.
		utility::vector1< core::Size > procs_per_hierarchy_level; //The number of processes at each level of the communications hierarchy.  Should be 1, n>=1, m>=n, etc.
		std::string sort_by("");
		bool select_highest(false);
		core::Real output_fraction(1.0);
		std::string output_filename("");
		core::Real lambda( 0.5 );
		core::Real kbt( 1.0 );
		core::Size threads_per_worker(1);
		bool compute_rmsd_to_lowest( false );
		core::Real compute_pnear_to_lowest_fract(0.0);
		bool compute_sasa_metrics( false );
		protocols::cyclic_peptide_predict::set_MPI_vars(
			MPI_rank, MPI_n_procs, total_hierarchy_levels, procs_per_hierarchy_level, batchsize_per_level,
			sort_by, select_highest, output_fraction, output_filename, lambda, kbt, compute_rmsd_to_lowest,
			compute_pnear_to_lowest_fract, compute_sasa_metrics, threads_per_worker, "helical_bundle_predict"
		); //Get the values of these vars (only used in MPI mode).
#endif

#ifdef USEMPI
		HelicalBundlePredictApplication_MPI application(
			MPI_rank, MPI_n_procs, HelicalBundlePredictApplication::create_centroid_scorefunction(),
			HelicalBundlePredictApplication::create_fullatom_scorefunction(), total_hierarchy_levels,
			procs_per_hierarchy_level, batchsize_per_level, sort_by, select_highest, output_fraction,
			output_filename, lambda, kbt, compute_rmsd_to_lowest, compute_pnear_to_lowest_fract, compute_sasa_metrics,
			threads_per_worker
		);
		application.set_options( options );
#else // !USEMPI
		HelicalBundlePredictApplication application( options );
		set_file_output( application );
#endif //USEMPI

		if (
#ifdef USEMPI
			MPI_rank == 0 &&
#endif //USEMPI
				TR.visible()
				) {
			TR << "Starting helical_bundle_predict application." << std::endl;
			TR << "Application created 15 October 2018 by Vikram K. Mulligan, Flatiron Institute." << std::endl;
			TR << "Note that this application is currently unpublished.  If you use it in a publication, please include the original author." << std::endl;
			TR << "For questions, please contact vmulligan@flatironinstitute.org." << std::endl;
		}

		//Register with citation manager:
		basic::citation_manager::CitationManager::get_instance()->add_citation(
			utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
			"helical_bundle_predict",
			basic::citation_manager::CitedModuleType::Application,
			"Vikram K. Mulligan",
			"Center for Computational Biology, Flatiron Institute",
			"vmulligan@flatironinstitute.org",
			"Wrote the application and developed the protocol."
			)
		);

#ifdef USEMPI
	if( MPI_rank == 0 ) {
		TR << "The MPI mode of helical_bundle_predict is using " << MPI_n_procs << " MPI processes distributed into " << total_hierarchy_levels << " layers." << std::endl;
#ifdef MULTI_THREADED
		if( threads_per_worker > 1 ) {
			TR << "Each process will launch " << threads_per_worker << " threads." << std::endl;
		}
#endif //MULTI_THREADED
	}
#endif //USEMPI

		application.run();

#ifdef USEMPI
		if( MPI_rank == 0 )
#endif
		{
			//Output from CitationManager:
			basic::citation_manager::CitationManager::get_instance()->write_all_citations_and_unpublished_author_info();
			if ( TR.visible() ) {
				TR << "Finished execution of helical_bundle_predict application." << std::endl;
				TR.flush();
			}
		}

	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

#ifdef USEMPI
	MPI_Barrier(MPI_COMM_WORLD); //Wait until all procs reach this point.
	MPI_Finalize(); //MPI finalization and cleanup.
#endif

	return 0;
}
