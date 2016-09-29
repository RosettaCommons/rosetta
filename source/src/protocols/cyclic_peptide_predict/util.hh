// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/util.hh
/// @brief Utility code for the simple_cycpep_predict app and its MPI version.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.  This
/// version uses MPI for cross-communication with parallel processes.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_predict_util_hh
#define INCLUDED_protocols_cyclic_peptide_predict_util_hh

// Unit Headers
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI_JobResultsSummary.fwd.hh>

// Package Headers
#include <core/types.hh>

// Project Headers

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <cstdio>
#include <string>
#include <map>

namespace protocols {
namespace cyclic_peptide_predict {

enum SIMPLE_CYCPEP_PREDICT_MPI_SORT_TYPE {
	SORT_BY_ENERGIES=1,
	SORT_BY_RMSD,
	SORT_BY_HBONDS,
	SORT_BY_CIS_PEPTIDE_BONDS
};

/// @brief Given a list of job summaries, sort these by some criterion (e.g. energies, rmsd, etc.) from lowest to highest.
/// @brief Uses the quicksort algorithm (recursive).
/// @param[in,out] jobsummaries The list of job summaries.  At the end of this operation, this is sorted from lowest to highest by the criterion specified.
/// @param[in] sort_type The criterion based on which we're sorting the list.
void sort_jobsummaries_list(
	utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > &jobsummaries,
	SIMPLE_CYCPEP_PREDICT_MPI_SORT_TYPE const sort_type
);

/// @brief Given an existing, sorted list of job summaries and a new, sorted list of job summaries, merge the two lists
/// into a new, sorted list and replace the existing list with the result.
/// @param[in,out] list_to_sort_into This is an existing, sorted list of job summaries, which will have elements added to it.
/// @param[in] additional_list These are the new job summaries (presorted), which will be merged with the existing list.
/// @param[in] sort_type The criterion based on which we're sorting the list.
void mergesort_jobsummaries_list (
	utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > &list_to_sort_into,
	utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > const &additional_list,
	SIMPLE_CYCPEP_PREDICT_MPI_SORT_TYPE const sort_type
);

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
void read_peptide_design_file(
	std::string const &filename,
	std::map < core::Size, utility::vector1 < std::string > > &allowed_canonicals_by_position,
	std::map < core::Size, utility::vector1 < std::string > > &allowed_noncanonicals_by_position
);

/// @brief Given a residue name, return true if this is one of the 20
/// canonical amino acids, false otherwise.
bool is_canonical( std::string const &resname );

/// @brief Given an ASCII file name, read the contents into a string.  If from_database is true, the read is from the database.
/// @param[out] output_string The string that will be filled with the file contents.  Overwritten by this operation.
/// @param[in] filename The name of the file to read.  If from_database is true, this is a relative database path.
/// @param[in] from_database If true, the file is assumed to be in the database.  If false, the path is relative the execution
/// directory, or is absolute.
void read_file_into_string( std::string &output_string, std::string const &filename, bool const from_database);

} //cyclic_peptide
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_predict_util_hh
