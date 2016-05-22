// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/cyclic_peptide_predict/util.cc
/// @brief Utility code for the simple_cycpep_predict app and its MPI version.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)

// Unit Headers
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI_JobResultsSummary.hh>
#include <protocols/cyclic_peptide_predict/util.hh>

// Package Headers
#include <basic/options/option.hh>
#include <core/types.hh>

// option key includes

//numeric headers

// Utility headers
#include <basic/Tracer.hh>

// C++ headers
#include <stdio.h>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide_predict.util" );

namespace protocols {
namespace cyclic_peptide_predict {

/// @brief Given a list of job summaries, sort these by some criterion (e.g. energies, rmsd, etc.) from lowest to highest.
/// @brief Uses the quicksort algorithm (recursive).
/// @param[in,out] jobsummaries The list of job summaries.  At the end of this operation, this is sorted from lowest to highest by the criterion specified.
/// @param[in] sort_type The criterion based on which we're sorting the list.
void
sort_jobsummaries_list(
	utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > &jobsummaries,
	SIMPLE_CYCPEP_PREDICT_MPI_SORT_TYPE const sort_type
) {
	if ( jobsummaries.size() < 2 ) return; //Do nothing for zero- or one-length arrays.
	SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP pivot( jobsummaries[1] );

	utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > lesser_summaries;
	utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > greater_summaries;

	for ( core::Size i=2, imax=jobsummaries.size(); i<=imax; ++i ) {
		if ( sort_type == SORT_BY_ENERGIES ) {
			if ( jobsummaries[i]->pose_energy() < pivot->pose_energy() ) {
				lesser_summaries.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP(jobsummaries[i]) );
			} else {
				greater_summaries.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP(jobsummaries[i]) );
			}
		} else if ( sort_type == SORT_BY_RMSD ) {
			if ( jobsummaries[i]->rmsd() < pivot->rmsd() ) {
				lesser_summaries.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP(jobsummaries[i]) );
			} else {
				greater_summaries.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP(jobsummaries[i]) );
			}
		} else if ( sort_type == SORT_BY_HBONDS ) {
			if ( jobsummaries[i]->hbonds() < pivot->hbonds() ) {
				lesser_summaries.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP(jobsummaries[i]) );
			} else {
				greater_summaries.push_back( SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP(jobsummaries[i]) );
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
	utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > &list_to_sort_into,
	utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > const &additional_list,
	SIMPLE_CYCPEP_PREDICT_MPI_SORT_TYPE const sort_type
) {
	//Quick and easy if one list is empty:
	if ( list_to_sort_into.size() == 0 ) {
		list_to_sort_into = additional_list;
		return;
	}
	if ( additional_list.size() == 0 ) {
		return; //Do nothing in this case.
	}

	utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > const origlist_copy( list_to_sort_into );  //Make a copy of the original list.
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
		}
	}

	debug_assert( original_index == origlist_size + 1);
	debug_assert( additional_index == additional_size + 1);
}

} //cyclic_peptide_predict
} //protocols
