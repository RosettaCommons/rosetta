// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/util.cc
/// @brief Utility code for the simple_cycpep_predict app and its MPI version.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)

#ifdef USEMPI
#include <mpi.h> //Message Passing Interface for parallelization -- keep this first
#endif //USEMPI

// Unit Headers
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_JobResultsSummary.hh>
#include <protocols/cyclic_peptide_predict/util.hh>

// Package Headers
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/cyclic_peptide.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

// C++ headers
#include <cstdio>

static basic::Tracer TR( "protocols.cyclic_peptide_predict.util" );

namespace protocols {
namespace cyclic_peptide_predict {

#ifdef USEMPI
/// @brief If rank > 0, wait for a message from proc 0 to continue.  If rank
/// is zero, send the message.
/// @details Only used in MPI mode.
void
wait_for_proc_zero() {
	char mybuf('A'); //A dummy piece of information to send.
	MPI_Bcast( &mybuf, 1, MPI_CHAR, 0, MPI_COMM_WORLD );
}

/// @brief In MPI mode, get the number of hierarchy levels from the options system.
utility::vector1< long int >
get_user_specified_hierarchy_levels(
	int const MPI_rank,
	int const MPI_n_procs,
	std::string const & errormsg
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(MPI_rank == 0 ) {
		runtime_assert_string_msg( option[cyclic_peptide::MPI_processes_by_level].user() != /*XOR*/ (option[cyclic_peptide::MPI_auto_2level_distribution].user() && option[cyclic_peptide::MPI_auto_2level_distribution].user() ),
		errormsg + "In MPI mode, the user must EITHER use the \"-cyclic_peptide:MPI_auto_2level_distribution\" flag, or specify the number of communication levels and the number of processes in each level, with the \"-cyclic_peptide:MPI_processes_by_level\" flag." );
	}
	MPI_Barrier(MPI_COMM_WORLD); //Wait until all procs reach this point -- ensures that processes 0 does this check first.

	utility::vector1< long int > levelvect;
	if( option[cyclic_peptide::MPI_processes_by_level].user() ) {
		levelvect = option[cyclic_peptide::MPI_processes_by_level]();
	} else if ( option[cyclic_peptide::MPI_auto_2level_distribution].user() && option[cyclic_peptide::MPI_auto_2level_distribution].value() ) {
		levelvect.resize(2);
		levelvect[1] = 1;
		runtime_assert( MPI_n_procs > 1 );
		levelvect[2] = MPI_n_procs - 1;
	}
	debug_assert( levelvect.size() > 0 );
	return levelvect;
}

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
	std::string &output_filename,
	core::Real &lambda,
	core::Real &kbt,
	bool &compute_rmsd_to_lowest,
	core::Real &compute_pnear_to_this_fract,
	bool &compute_sasa_metrics,
	core::Size &threads_per_worker,
	std::string const & app_name
) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string const errormsg( "Error in " + app_name + " application: ");

	//Get the rank and number of processes.
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &MPI_n_procs );

	if (TR.Debug.visible()) {
		TR.Debug << "MPI rank " << MPI_rank << " of " << MPI_n_procs << " reporting for duty." << std::endl;
	}

	utility::vector1< long int > user_specified_hierarchy_levels( get_user_specified_hierarchy_levels( MPI_rank, MPI_n_procs, errormsg ) );
	if(MPI_rank == 0 ) {
		runtime_assert_string_msg( user_specified_hierarchy_levels.size() > 1, errormsg + "The number of communication levels specified with the \"-cyclic_peptide:MPI_processes_by_level\" flag must be greater than one." );
		runtime_assert_string_msg( user_specified_hierarchy_levels.size() <= static_cast<core::Size>(MPI_n_procs), errormsg + "The number of communication levels specified with the \"-cyclic_peptide:MPI_processes_by_level\" flag must be less than or equal to the number of processes launched." );
		runtime_assert_string_msg( user_specified_hierarchy_levels[1] == 1, errormsg + "The first communication level must have a single director process.  This will control intermediate managers and front-line workers.");
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
		runtime_assert_string_msg( !option[out::file::o].user(), errormsg + "The -out:file:o option cannot be used with the " + app_name + " app in MPI mode.  Only silent file output is permitted." );
		runtime_assert_string_msg( option[out::file::silent].user(), errormsg + "A silent file for output must be specified for the " + app_name + " app in MPI mode." );
	}
	output_filename = option[out::file::silent]();

	lambda = option[cyclic_peptide::MPI_pnear_lambda]();
	kbt = option[cyclic_peptide::MPI_pnear_kbt]();
	if( MPI_rank == 0 ) {
		runtime_assert_string_msg( lambda > 0, errormsg + "The -cyclic_peptide:MPI_pnear_lambda option must be greater than zero." );
		runtime_assert_string_msg( kbt > 0, errormsg + "The -cyclic_peptide:MPI_pnear_kbt option must be greater than zero." );
	}
	compute_rmsd_to_lowest = option[cyclic_peptide::compute_rmsd_to_lowest].value();
	compute_pnear_to_this_fract = option[ cyclic_peptide::compute_pnear_to_this_fract ].value();
	runtime_assert_string_msg( compute_pnear_to_this_fract >= 0.0 && compute_pnear_to_this_fract <= 1.0, errormsg + "The -cyclic_peptide:compute_pnear_to_this_fract option must be accompanied by a floating-point number between 0 and 1 inclusive.  Value " + std::to_string(compute_pnear_to_this_fract) + " is invalid." );
	compute_sasa_metrics = option[cyclic_peptide::compute_ensemble_sasa_metrics].value();

#ifdef MULTI_THREADED
	runtime_assert_string_msg( option[cyclic_peptide::threads_per_worker]() > 0, errormsg + "The -cyclic_peptide:threads_per_worker option's value must be greater than zero." );
	threads_per_worker = static_cast<core::Size>( option[cyclic_peptide::threads_per_worker]() );
#else
	threads_per_worker = 1;
#endif
	wait_for_proc_zero();
}
#endif //USEMPI

/// @brief Given a list of job summaries, sort these by some criterion (e.g. energies, rmsd, etc.) from lowest to highest.
/// @brief Uses the quicksort algorithm (recursive).
/// @param[in,out] jobsummaries The list of job summaries.  At the end of this operation, this is sorted from lowest to highest by the criterion specified.
/// @param[in] sort_type The criterion based on which we're sorting the list.
void
sort_jobsummaries_list(
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > &jobsummaries,
	HIERARCHICAL_HYBRID_JD_MPI_SORT_TYPE const sort_type
) {
	if ( jobsummaries.size() < 2 ) return; //Do nothing for zero- or one-length arrays.
	HierarchicalHybridJD_JobResultsSummaryOP pivot( jobsummaries[1] );

	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > lesser_summaries;
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > greater_summaries;

	for ( core::Size i=2, imax=jobsummaries.size(); i<=imax; ++i ) {
		if ( sort_type == SORT_BY_ENERGIES ) {
			if ( jobsummaries[i]->pose_energy() < pivot->pose_energy() ) {
				lesser_summaries.push_back( HierarchicalHybridJD_JobResultsSummaryOP(jobsummaries[i]) );
			} else {
				greater_summaries.push_back( HierarchicalHybridJD_JobResultsSummaryOP(jobsummaries[i]) );
			}
		} else if ( sort_type == SORT_BY_RMSD ) {
			if ( jobsummaries[i]->rmsd() < pivot->rmsd() ) {
				lesser_summaries.push_back( HierarchicalHybridJD_JobResultsSummaryOP(jobsummaries[i]) );
			} else {
				greater_summaries.push_back( HierarchicalHybridJD_JobResultsSummaryOP(jobsummaries[i]) );
			}
		} else if ( sort_type == SORT_BY_HBONDS ) {
			if ( jobsummaries[i]->hbonds() < pivot->hbonds() ) {
				lesser_summaries.push_back( HierarchicalHybridJD_JobResultsSummaryOP(jobsummaries[i]) );
			} else {
				greater_summaries.push_back( HierarchicalHybridJD_JobResultsSummaryOP(jobsummaries[i]) );
			}
		} else if ( sort_type == SORT_BY_CIS_PEPTIDE_BONDS ) {
			if ( jobsummaries[i]->cis_peptide_bonds() < pivot->cis_peptide_bonds() ) {
				lesser_summaries.push_back( HierarchicalHybridJD_JobResultsSummaryOP(jobsummaries[i]) );
			} else {
				greater_summaries.push_back( HierarchicalHybridJD_JobResultsSummaryOP(jobsummaries[i]) );
			}
		}
	}

	sort_jobsummaries_list( lesser_summaries, sort_type );
	sort_jobsummaries_list( greater_summaries, sort_type );

	debug_assert( lesser_summaries.size() + greater_summaries.size() + 1 == jobsummaries.size() );

	core::Size i(1), j(1);
	core::Size const lessersize(lesser_summaries.size());
	do {
		if ( i <= lessersize ) {
			jobsummaries[i] = lesser_summaries[i];
		} else if ( i == lessersize + 1 ) {
			jobsummaries[i] = pivot;
		} else {
			jobsummaries[i] = greater_summaries[j];
			++j;
		}
		++i;
	} while ( i<=jobsummaries.size() );

	debug_assert( j == greater_summaries.size() + 1 );
}

/// @brief Given an existing, sorted list of job summaries and a new, sorted list of job summaries, merge the two lists
/// into a new, sorted list and replace the existing list with the result.
/// @param[in,out] list_to_sort_into This is an existing, sorted list of job summaries, which will have elements added to it.
/// @param[in] additional_list These are the new job summaries (presorted), which will be merged with the existing list.
/// @param[in] sort_type The criterion based on which we're sorting the list.
void mergesort_jobsummaries_list (
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > &list_to_sort_into,
	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > const &additional_list,
	HIERARCHICAL_HYBRID_JD_MPI_SORT_TYPE const sort_type
) {
	//Quick and easy if one list is empty:
	if ( list_to_sort_into.size() == 0 ) {
		list_to_sort_into = additional_list;
		return;
	}
	if ( additional_list.size() == 0 ) {
		return; //Do nothing in this case.
	}

	utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > const origlist_copy( list_to_sort_into );  //Make a copy of the original list.
	core::Size const origlist_size(origlist_copy.size());
	core::Size const additional_size(additional_list.size());
	core::Size const combined_size( origlist_size + additional_size );


	list_to_sort_into.resize( combined_size );

	core::Size original_index(1);
	core::Size additional_index(1);
	for ( core::Size i=1; i<=combined_size; ++i ) {
		if ( sort_type == SORT_BY_ENERGIES ) {
			if ( additional_index > additional_size || ( !(original_index>origlist_size) && origlist_copy[original_index]->pose_energy() < additional_list[additional_index]->pose_energy() ) ) {
				list_to_sort_into[i] = origlist_copy[original_index];
				++original_index;
			} else {
				list_to_sort_into[i] = additional_list[additional_index];
				++additional_index;
			}
		} else if ( sort_type == SORT_BY_RMSD ) {
			if ( additional_index > additional_size || ( !(original_index>origlist_size) && origlist_copy[original_index]->rmsd() < additional_list[additional_index]->rmsd() ) ) {
				list_to_sort_into[i] = origlist_copy[original_index];
				++original_index;
			} else {
				list_to_sort_into[i] = additional_list[additional_index];
				++additional_index;
			}
		} else if ( sort_type == SORT_BY_HBONDS ) {
			if ( additional_index > additional_size || ( !(original_index>origlist_size) && origlist_copy[original_index]->hbonds() < additional_list[additional_index]->hbonds() ) ) {
				list_to_sort_into[i] = origlist_copy[original_index];
				++original_index;
			} else {
				list_to_sort_into[i] = additional_list[additional_index];
				++additional_index;
			}
		} else if ( sort_type == SORT_BY_CIS_PEPTIDE_BONDS ) {
			if ( additional_index > additional_size || ( !(original_index>origlist_size) && origlist_copy[original_index]->cis_peptide_bonds() < additional_list[additional_index]->cis_peptide_bonds() ) ) {
				list_to_sort_into[i] = origlist_copy[original_index];
				++original_index;
			} else {
				list_to_sort_into[i] = additional_list[additional_index];
				++additional_index;
			}
		}
	}

	debug_assert( original_index == origlist_size + 1);
	debug_assert( additional_index == additional_size + 1);
}

/// @brief Given a filename, read and parse the file, returning a list of canonical residues allowed at each position
/// and a list of noncanonicals allowed at each position.
/// @details This does the actual file read.  It is NOT THREADSAFE.  The file format is a series of lines with the pattern:
/// residue_index residuetype_1_fullname residuetype_2_fullname residuetype_3_fullname ...
/// Anything after a pound sign should be ignored.  A line with DEFAULT in place of the resiude index should be interpreted
/// as providing default values, which should be stored as map key 0.
/// @param[in] filename The file name from which we'll read.
/// @param[out] allowed_canonicals_by_position A map of [position->vector of strings of full names] listing the allowed canonical
/// residue types at each position.  Reset and populated by this function.  Key 0 indicates default settings applied anywhere
/// that lacks a map key.
/// @param[out] allowed_noncanonicals_by_position A map of [position->vector of strings of full names] listing the allowed noncanonical
/// residue types at each position.  Reset and populated by this function.  Key 0 indicates default settings applied anywhere
/// that lacks a map key.
void
read_peptide_design_file(
	std::string const &filename,
	std::map < core::Size, utility::vector1 < std::string > > &allowed_canonicals_by_position,
	std::map < core::Size, utility::vector1 < std::string > > &allowed_noncanonicals_by_position
) {
	//Clear these maps:
	allowed_canonicals_by_position.clear();
	allowed_noncanonicals_by_position.clear();

	//Open the input file:
	utility::io::izstream infile;
	infile.open(filename);
	runtime_assert_string_msg( infile.good(), "Error in protocols::cyclic_peptide_predict::read_peptide_design_file():  Unable to open the specified file (\"" + filename + "\") for read!" );

	TR << "Opened " << filename << " for read." << std::endl;

	std::string curline(""); //Buffer for current line.
	std::string curline2("");  //Second buffer.

	//Read the file:
	while ( getline(infile, curline) ) {
		if ( curline.size() < 1 ) continue; //Ignore blank lines.

		//Terminate line at comment (pound) symbol, and trim terminal whitespace:
		std::string::size_type pound = curline.find('#', 0);
		if ( pound != std::string::npos ) {
			curline2 = curline.substr(0, pound);
		} else {
			curline2 = curline;
		}
		utility::trim( curline2, " \n\t" );

		if ( curline2.size() < 1 ) continue; //Ignore lines with nothing but comments and/or whitespace.

		std::istringstream curline3( curline2 );

		//I'll recycle the curline buffer for parsing the line.
		//First, get the residue index.
		curline3 >> curline;
		core::Size curres(0);
		if ( curline != "DEFAULT" ) {
			std::istringstream curline4( curline );
			curline4 >> curres;
			runtime_assert_string_msg( !curline4.fail() && !curline4.bad(), "Error in protocols::cyclic_peptide_predict::read_peptide_design_file():  Unable to parse \"" + curline + "\" as a residue index." );
		}

		runtime_assert_string_msg( allowed_canonicals_by_position.count(curres) == 0 && allowed_noncanonicals_by_position.count(curres) == 0,
			"Error in protocols::cyclic_peptide_predict::read_peptide_design_file(): Residue " + curline + " has already been defined." );

		utility::vector1 < std::string > canonicals;
		utility::vector1 < std::string > noncanonicals;

		while ( !curline3.eof() ) {
			curline3 >> curline;
			if ( is_canonical(curline) ) {
				canonicals.push_back(curline);
			} else {
				noncanonicals.push_back(curline);
			}
		}

		allowed_canonicals_by_position[curres] = canonicals;
		allowed_noncanonicals_by_position[curres] = noncanonicals;
	} //while getline

	infile.close();

	TR << "Read " << filename << " and closed file." << std::endl;
}

/// @brief Given a residue name, return true if this is one of the 20
/// canonical amino acids, false otherwise.
bool
is_canonical(
	std::string const &resname
) {
	if ( resname == "ALA" || resname == "CYS" || resname == "ASP" || resname == "GLU" || resname == "PHE" || resname == "GLY" ||
			resname == "HIS" || resname == "HIS_D" || resname == "ILE" || resname == "LYS" || resname == "LEU" || resname == "MET" ||
			resname == "ASN" || resname == "PRO"   || resname == "GLN" || resname == "ARG" || resname == "SER" || resname == "THR" ||
			resname == "VAL" || resname == "TRP"   || resname == "TYR" ) {
		return true;
	}
	return false;
}

/// @brief Given an ASCII file name, read the contents into a string.  If from_database is true, the read is from the database.
/// @param[out] output_string The string that will be filled with the file contents.  Overwritten by this operation.
/// @param[in] filename The name of the file to read.  If from_database is true, this is a relative database path.
/// @param[in] from_database If true, the file is assumed to be in the database.  If false, the path is relative the execution
/// directory, or is absolute.
void
read_file_into_string(
	std::string &output_string,
	std::string const &filename,
	bool const from_database
) {
	using namespace utility::io;

	izstream infile;
	if ( from_database ) {
		basic::database::open( infile, filename );
	} else {
		infile.open( filename );
	}
	runtime_assert_string_msg( infile.good(), "Error in protocols::cyclic_peptide_predict::read_file_into_string():  Unable to open file " + filename + " for read!" );
	output_string.clear();
	utility::slurp( infile, output_string);
	infile.close();
}

} //cyclic_peptide_predict
} //protocols
