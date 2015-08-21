// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/optimize_weights/DGBindOptEData.cc
/// @brief Implementation of the OptEPositionData class that handles interface ddG optimization
/// @author Ron Jacak

#ifdef USEMPI
#include <mpi.h>
#endif

// Unit headers
#include <protocols/optimize_weights/OptEData.hh>
#include <protocols/optimize_weights/DDGBindOptEData.hh>
#include <ObjexxFCL/format.hh>
#include <utility/vector1.functions.hh>  // to get arg_min
#include <basic/Tracer.hh>

// option key includes

#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>


using namespace core;
using namespace core::scoring;
using namespace ObjexxFCL::format;

namespace protocols {
namespace optimize_weights {

static thread_local basic::Tracer TR( "DDGBindOptEData" );

typedef core::chemical::AA AA;

///
/// @brief
/// Initialize all of the member variables to 0.
///
DDGBindOptEData::DDGBindOptEData():
	experimental_ddG_bind_(0.0)
{}

///
DDGBindOptEData::~DDGBindOptEData() {}

///
Real
DDGBindOptEData::get_score(
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	Size const num_energy_dofs,
	int const num_ref_dofs,
	int const num_total_dofs,
	EnergyMap const & fixed_terms,
	ScoreTypes const & free_score_list,
	ScoreTypes const & fixed_score_list
) const
{
	return process_score( TR, false, component_weights, vars, dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs, fixed_terms, free_score_list, fixed_score_list );
}

///
void
DDGBindOptEData::print_score(
	std::ostream & ostr,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	Size const num_energy_dofs,
	int const num_ref_dofs,
	int const num_total_dofs,
	EnergyMap const & fixed_terms,
	ScoreTypes const & free_score_list,
	ScoreTypes const & fixed_score_list
) const
{
	process_score( ostr, true, component_weights, vars, dE_dvars, num_energy_dofs, num_ref_dofs, num_total_dofs, fixed_terms, free_score_list, fixed_score_list );
}

///
/// @brief
/// One method to do the score processing which takes a boolean dictating whether to print to an ostream or not. With this function, changes
/// to how scoring works only need to be made in one place as opposed to two (when get_score() and print_score() both had scoring logic in them).
/// Avoiding code duplication is good.
///
Real
DDGBindOptEData::process_score(
	std::ostream & ostr,
	bool print,
	optimization::Multivec const & component_weights,
	optimization::Multivec const & vars,
	optimization::Multivec & dE_dvars,
	/// Basically, turn over all the private data from OptEMultiFunc
	Size const num_energy_dofs,
	int const num_ref_dofs,
	int const,
	EnergyMap const & fixed_terms,
	ScoreTypes const & /* score_list */,
	ScoreTypes const & fixed_score_list
) const
{
	using namespace core;
	using namespace core::optimization;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility;

	// if there are no structures to go through, return immediately
	if ( wt_complexes_.size() == 0 || mutant_complexes_.size() == 0 ||
			wt_unbounds_.size()  == 0 || mutant_unbounds_.size()  == 0 ) return 0.0;


	// these vectors are sized to the number of structures there are for each of the four types of structures
	// they'll be used to determine which structure has the best energy
	utility::vector1< Real > wt_complexes_energies( wt_complexes_.size(), 0.0 );
	utility::vector1< Real > mutant_complexes_energies( mutant_complexes_.size(), 0.0 );
	utility::vector1< Real > wt_unbounds_energies( wt_unbounds_.size(), 0.0 );
	utility::vector1< Real > mutant_unbounds_energies( mutant_unbounds_.size(), 0.0 );

	// go through and come up with a total score for each structure in the 4 types of structues
	//
	// wt_complexes_ is a vector1 of SingleStructureData (SSD) objects. this for loop iterates over each free weight and
	// takes the unweighted energy for the current free term in wt_complexes_[jj] and multiplies it by the weight in vars.
	// The SSD objects have a vector1 of Reals accessible by free_data() and fixed_data() member functions.  These store
	// the total unweighted energies by score type for the entire pose.

	for ( Size ii = 1; ii <= num_energy_dofs; ++ii ) {
		for ( Size jj = 1; jj <= wt_complexes_.size(); ++jj ) {
			// cap the fa_rep term at some value - this at least keeps it around for most of the mutants; this is how it would be done for wt_complexes
			//#ifdef CAP_FA_REP
			// if ( ( score_list[ ii ] == fa_rep ) && ( vars[ ii ] * wts_[ jj ]->free_data()[ ii ] > 10 ) ) { wt_energies[ jj ] += 10; }
			// else
			//#endif
			wt_complexes_energies[ jj ] += vars[ ii ] * wt_complexes_[ jj ]->free_data()[ ii ];
		}
		for ( Size jj = 1; jj <= mutant_complexes_.size(); ++jj ) {
			mutant_complexes_energies[ jj ] += vars[ ii ] * mutant_complexes_[ jj ]->free_data()[ ii ];
		}
		for ( Size jj = 1; jj <= wt_unbounds_.size(); ++jj ) {
			wt_unbounds_energies[ jj ] += vars[ ii ] * wt_unbounds_[ jj ]->free_data()[ ii ];
		}
		for ( Size jj = 1; jj <= mutant_unbounds_.size(); ++jj ) {
			mutant_unbounds_energies[ jj ] += vars[ ii ] * mutant_unbounds_[ jj ]->free_data()[ ii ];
		}
	}

	for ( Size ii = 1; ii <= fixed_score_list.size(); ++ii ) {
		for ( Size jj = 1; jj <= wt_complexes_.size(); ++jj ) {
			//#ifdef CAP_FA_REP
			// if ( ( fixed_score_list[ ii ] == fa_rep ) && ( fixed_terms[ fixed_score_list[ ii ] ] * wts_[ jj ]->fixed_data()[ ii ] > 10 ) ) { wt_energies[ jj ] += 10; }
			// else
			//#endif
			wt_complexes_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * wt_complexes_[ jj ]->fixed_data()[ ii ];
		}
		for ( Size jj = 1; jj <= mutant_complexes_.size(); ++jj ) {
			mutant_complexes_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * mutant_complexes_[ jj ]->fixed_data()[ ii ];
		}
		for ( Size jj = 1; jj <= wt_unbounds_.size(); ++jj ) {
			wt_unbounds_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * wt_unbounds_[ jj ]->fixed_data()[ ii ];
		}
		for ( Size jj = 1; jj <= mutant_unbounds_.size(); ++jj ) {
			mutant_unbounds_energies[ jj ] += fixed_terms[ fixed_score_list[ ii ] ] * mutant_unbounds_[ jj ]->fixed_data()[ ii ];
		}
	}

	// num_energy_dofs is the number of free, non-reference energy parameters in the run; num_ref_dofs are the reference energies
	// but only add in the reference energy if they're being optimized in this optE run (not always the case)
	if ( num_ref_dofs != 0 ) {
		for ( Size ii = 1; ii <= mutations_.size(); ++ii ) {
			for ( Size jj = 1; jj <= wt_complexes_.size(); ++jj ) {
				wt_complexes_energies[ jj ] += vars[ num_energy_dofs + mutations_[ ii ].second.first ];
			}
			for ( Size jj = 1; jj <= mutant_complexes_.size(); ++jj ) {
				mutant_complexes_energies[ jj ] += vars[ num_energy_dofs + mutations_[ ii ].second.second ];
			}
			for ( Size jj = 1; jj <= wt_unbounds_.size(); ++jj ) {
				wt_unbounds_energies[ jj ] += vars[ num_energy_dofs + mutations_[ ii ].second.first ];
			}
			for ( Size jj = 1; jj <= mutant_unbounds_.size(); ++jj ) {
				mutant_unbounds_energies[ jj ] += vars[ num_energy_dofs + mutations_[ ii ].second.second ];
			}
		}
	}

	//TR << "process_score(): weighted structure energies, wt complexes: [ ";
	//for ( Size jj = 1; jj <= wt_complexes_.size(); ++jj ) { TR << F(6,1,wt_complexes_energies[ jj ]) << ", "; }
	//TR << "]" << std::endl;

	// This is where we branch on how the score is calculated.  The simplest approach is to take the minimum energy
	// of all the wts and all the muts and subtract them to get the ddG.  The mean-based approach uses the difference
	// of the average of all muts and average of all wts. Finally, the boltzmann approach calculates a boltzmann
	// probability for the wts and muts to get a score.
	// The two latter ways have never been implemented for this.

	// Do things the old-fashioned way: best energy mut - best energy wt

	// the next four lines just get the index to the best energy in each vector
	Size const best_wt_complex_index      = arg_min( wt_complexes_energies );
	Size const best_mutant_complex_index  = arg_min( mutant_complexes_energies );
	Size const best_wt_unbounds_index     = arg_min( wt_unbounds_energies );
	Size const best_mutant_unbounds_index = arg_min( mutant_unbounds_energies );

	Real const best_wt_complex_energy      = wt_complexes_energies[     best_wt_complex_index ];
	Real const best_mutant_complex_energy  = mutant_complexes_energies[ best_mutant_complex_index ];
	Real const best_wt_unbounds_energy     = wt_unbounds_energies[      best_wt_unbounds_index ];
	Real const best_mutant_unbounds_energy = mutant_unbounds_energies[  best_mutant_unbounds_index ];

	Real dG_bind_wt     = best_wt_complex_energy     - best_wt_unbounds_energy;
	Real dG_bind_mutant = best_mutant_complex_energy - best_mutant_unbounds_energy;

	Real predicted_ddG_bind = dG_bind_mutant     - dG_bind_wt;
	Real ddG_bind_diff      = predicted_ddG_bind - experimental_ddG_bind_;
	Real ddG_bind_diff_sq   = ddG_bind_diff * ddG_bind_diff;

	// It might be good to adjust the slope here in the same way that Ian Davis was doing for dG bind calculations.
	// PdbBind 2007 core set
	// Linear regression of components of interface_delta gave an intercept between -3 and -4;
	// including this allows a better overall fit to the real binding energy.
	//Real const TdS = 3.5; // kcal/mol
	//Real const pred_dG = (boundE - unboundE) - TdS;
	//Real const diff_dG = pred_dG - deltaG_bind_;
	//Real const sq_err = diff_dG * diff_dG;

	if ( print ) {
		ostr << "DDGBind " << A( 20, tag() ) << X(1) << "pred: " << F(6,2,predicted_ddG_bind) << " exp: " << F(6,2,experimental_ddG_bind_)
			<< " diff^2: " << F(7,2,ddG_bind_diff_sq )
			<< " cmptwt_diff^2: " << F(7,2,component_weights[ ddG_bind_correlation ]*ddG_bind_diff_sq)
			<< std::endl;

		TR << "process_score(): "
			<< "wt bound[ " << best_wt_complex_energy << " ] - wt unbound[ " << best_wt_unbounds_energy << " ] = " << dG_bind_wt
			<< ", mut bound[ " << best_mutant_complex_energy << " ] - mut unbound[ " << best_mutant_unbounds_energy << " ] = " << dG_bind_mutant
			<< ", predicted[ " << predicted_ddG_bind << " ] - experimental[ " << experimental_ddG_bind_ << " ] = error[ " << ddG_bind_diff << " ], "
			<< "error^2[ " << ddG_bind_diff_sq << " ], "
			<< "weighted error^2[ " << component_weights[ ddG_bind_correlation ] * ddG_bind_diff_sq << " ]" << ", tag: " << this->tag()
			<< std::endl;

	} else {

		for ( core::Size e_dof = 1; e_dof <= num_energy_dofs; ++e_dof ) {
			// may want to deal with really bad repulsive energy cases here somehow
			dE_dvars[ e_dof ] += 2 * component_weights[ ddG_bind_correlation ] * ddG_bind_diff *
				( mutant_complexes_[ best_mutant_complex_index ]->free_data()[ e_dof ]
				- mutant_unbounds_[ best_mutant_unbounds_index ]->free_data()[ e_dof ] ) -
				( wt_complexes_[ best_wt_complex_index ]->free_data()[ e_dof ]
				- wt_unbounds_[ best_wt_unbounds_index ]->free_data()[ e_dof ] );
		}
		if ( num_ref_dofs != 0 ) {
			for ( Size ii = 1; ii <= mutations_.size(); ++ii ) {
				dE_dvars[ num_energy_dofs + mutations_[ ii ].second.second ] += 2 * component_weights[ ddG_bind_correlation ] * ddG_bind_diff;
				dE_dvars[ num_energy_dofs + mutations_[ ii ].second.first  ] -= 2 * component_weights[ ddG_bind_correlation ] * ddG_bind_diff;
			}
		}
	}

	return component_weights[ ddG_bind_correlation ] * ddG_bind_diff_sq;

}

///
OptEPositionDataType
DDGBindOptEData::type() const {
	return ddG_bind_correlation;
}

///
/// @details
/// Determine the upper and lower bounds on the unweighted component energy terms at this "position" and deposit them in the passed-in
/// EnergyMap objects.  Called by the IterativeOptE driver class to print extra information to the minimization data file.
///
void
DDGBindOptEData::range( ScoreTypes const & free_score_list, ScoreTypes const & fixed_score_list, EnergyMap & lower_bound, EnergyMap & upper_bound ) const {

	for ( Size ii = 1; ii <= wt_complexes_.size(); ++ii ) {
		update_range( wt_complexes_[ ii ],  free_score_list, fixed_score_list, lower_bound, upper_bound );
	}
	for ( Size ii = 1; ii <= mutant_complexes_.size(); ++ii ) {
		update_range( mutant_complexes_[ ii ], free_score_list, fixed_score_list, lower_bound, upper_bound );
	}
	for ( Size ii = 1; ii <= wt_unbounds_.size(); ++ii ) {
		update_range( wt_unbounds_[ ii ],  free_score_list, fixed_score_list, lower_bound, upper_bound );
	}
	for ( Size ii = 1; ii <= mutant_unbounds_.size(); ++ii ) {
		update_range( mutant_unbounds_[ ii ], free_score_list, fixed_score_list, lower_bound, upper_bound );
	}
}

///
core::Size
DDGBindOptEData::size() const {
	return wt_complexes_.size() + mutant_complexes_.size() + wt_unbounds_.size() + mutant_unbounds_.size();
}

///
/// Only used for user feedback. Nothing in the code uses the result from this to allocate memory.
///
core::Size
DDGBindOptEData::memory_use() const {

	Size total = sizeof( DDGBindOptEData ) +
		sizeof( SingleStructureData ) * wt_complexes_.size() +
		sizeof( SingleStructureData ) * mutant_complexes_.size() +
		sizeof( SingleStructureData ) * wt_unbounds_.size() +
		sizeof( SingleStructureData ) * mutant_unbounds_.size();

	if ( wt_complexes_.size() > 0 ) {
		total += sizeof( Real ) * ( wt_complexes_[ 1 ]->free_data().size() + wt_complexes_[ 1 ]->fixed_data().size() ) * wt_complexes_.size();
	}
	if ( mutant_complexes_.size() > 0 ) {
		total += sizeof( Real ) * ( mutant_complexes_[ 1 ]->free_data().size() + mutant_complexes_[ 1 ]->fixed_data().size() ) * mutant_complexes_.size();
	}
	if ( wt_unbounds_.size() > 0 ) {
		total += sizeof( Real ) * ( wt_unbounds_[ 1 ]->free_data().size() + wt_unbounds_[ 1 ]->fixed_data().size() ) * wt_unbounds_.size();
	}
	if ( mutant_unbounds_.size() > 0 ) {
		total += sizeof( Real ) * ( mutant_unbounds_[ 1 ]->free_data().size() + mutant_unbounds_[ 1 ]->fixed_data().size() ) * mutant_unbounds_.size();
	}

	return total;
}


#ifdef USEMPI
///
///
void
DDGBindOptEData::send_to_node( int const destination_node, int const tag ) const {

	/// 1. Experimental DDG, wt_aa, mut_aa
	Real experimental_ddG_bind = experimental_ddG_bind_; // stupid const pointer
	MPI_Send( & experimental_ddG_bind, 1, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	Size n_mutations = mutations_.size();
	MPI_Send( & n_mutations, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	for ( Size ii = 1; ii <= mutations_.size(); ++ii ) {
		Size position = mutations_[ ii ].first;
		MPI_Send( & position, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );
		int wt_aa( mutations_[ ii ].second.first ), mut_aa( mutations_[ ii ].second.second );
		MPI_Send( & wt_aa, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );
		MPI_Send( & mut_aa, 1, MPI_INT, destination_node, tag, MPI_COMM_WORLD );
	}

	// 2a. n wt complexes
	Size n_wt_complexes = wt_complexes_.size();
	//TR << "sending n wt complexes to node " << destination_node << ": " << n_wt_complexes <<  std::endl;
	MPI_Send( & n_wt_complexes, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// 2b. n mut complexes
	Size n_mut_complexes = mutant_complexes_.size();
	//TR << "sending n mut complexes to node " << destination_node << ": " << n_mut_complexes <<  std::endl;
	MPI_Send( & n_mut_complexes, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	// 2c. n wt unbounds
	Size n_wt_unbounds = wt_unbounds_.size();
	//TR << "sending n wt unbounds to node " << destination_node << ": " << n_wt_unbounds <<  std::endl;
	MPI_Send( & n_wt_unbounds, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// 2d. n mut unbounds
	Size n_mut_unbounds = mutant_unbounds_.size();
	//TR << "sending n mut unbounds to node " << destination_node << ": " << n_mut_unbounds <<  std::endl;
	MPI_Send( & n_mut_unbounds, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	if ( n_wt_complexes == 0 || n_mut_complexes == 0 || n_wt_unbounds == 0 || n_mut_unbounds == 0 )
		return;

	/// 3. n free weights
	Size n_free = wt_complexes_[ 1 ]->free_data().size();
	//TR << "sending n_free to node " << destination_node << ": " << n_free << std::endl;
	MPI_Send( & n_free, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );

	/// 4. n fixed weights
	Size n_fixed = wt_complexes_[ 1 ]->fixed_data().size();
	//TR << "sending n_fixed to node " << destination_node  << ": " << n_fixed << std::endl;
	MPI_Send( & n_fixed, 1, MPI_UNSIGNED_LONG, destination_node, tag, MPI_COMM_WORLD );


	/// --------- wt complexes ---------
	Real * free_data = new Real[ n_free * n_wt_complexes ];
	Real * fixed_data = new Real[ n_fixed * n_wt_complexes ];
	for ( Size ii = 1; ii <= n_wt_complexes; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ] = wt_complexes_[ ii ]->free_data()[ jj ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ] = wt_complexes_[ ii ]->fixed_data()[ jj ];
		}
	}
	//TR << "sending wt complexes free_data to node " << destination_node << ": " << free_data <<  std::endl;
	/// 5. wt complexes free data
	MPI_Send( free_data, n_wt_complexes * n_free, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	//TR << "sending wt complexes fixed_data to node " << destination_node << ": " << fixed_data <<  std::endl;
	/// 6. wt complexes fixed data
	MPI_Send( fixed_data, n_wt_complexes * n_fixed, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	delete [] free_data;
	delete [] fixed_data;

	/// --------- mut complexes ---------
	free_data = new Real[ n_free * n_mut_complexes ];
	fixed_data = new Real[ n_fixed * n_mut_complexes ];
	for ( Size ii = 1; ii <= n_mut_complexes; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ] = mutant_complexes_[ ii ]->free_data()[ jj ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ] = mutant_complexes_[ ii ]->fixed_data()[ jj ];
		}
	}
	/// 7. mut complexes free data
	//TR << "sending mut complexes free_data to node " << destination_node << ": " << free_data <<  std::endl;
	MPI_Send( free_data, n_mut_complexes * n_free, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 8. mut complexes fixed data
	//TR << "sending mut complexes fixed_data to node " << destination_node << ": " << fixed_data <<  std::endl;
	MPI_Send( fixed_data, n_mut_complexes * n_fixed, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	delete [] free_data;
	delete [] fixed_data;

	/// --------- wt unbounds ---------
	free_data = new Real[ n_free * n_wt_unbounds ];
	fixed_data = new Real[ n_fixed * n_wt_unbounds ];
	for ( Size ii = 1; ii <= n_wt_unbounds; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ] = wt_unbounds_[ ii ]->free_data()[ jj ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ] = wt_unbounds_[ ii ]->fixed_data()[ jj ];
		}
	}
	//TR << "sending wt unbounds free_data to node " << destination_node << ": " << free_data <<  std::endl;
	/// 9. wt unbounds free data
	MPI_Send( free_data, n_wt_unbounds * n_free, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	//TR << "sending wt unbounds fixed_data to node " << destination_node << ": " << fixed_data <<  std::endl;
	/// 10. fixed data
	MPI_Send( fixed_data, n_wt_unbounds * n_fixed, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	delete [] free_data;
	delete [] fixed_data;

	/// --------- mut unbounds ---------
	free_data = new Real[ n_free * n_mut_unbounds ];
	fixed_data = new Real[ n_fixed * n_mut_unbounds ];
	for ( Size ii = 1; ii <= n_mut_unbounds; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ] = mutant_unbounds_[ ii ]->free_data()[ jj ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ] = mutant_unbounds_[ ii ]->fixed_data()[ jj ];
		}
	}
	/// 11. mut unbounds free data
	//TR << "sending mut unbounds free_data to node " << destination_node << ": " << free_data <<  std::endl;
	MPI_Send( free_data, n_mut_unbounds * n_free, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	/// 12. mut unbounds fixed data
	//TR << "sending mut unbounds fixed_data to node " << destination_node << ": " << fixed_data <<  std::endl;
	MPI_Send( fixed_data, n_mut_unbounds * n_fixed, MPI_DOUBLE, destination_node, tag, MPI_COMM_WORLD );

	delete [] free_data;
	delete [] fixed_data;

	OptEPositionData::send_to_node( destination_node, tag );

}

///
void
DDGBindOptEData::receive_from_node( int const source_node, int const tag )
{
	MPI_Status stat;
	//TR << "receive_from_node(): receiving data from node... " << source_node << std::endl;

	/// 1. Experimental DDG, wt_aa, mut_aa
	MPI_Recv( &experimental_ddG_bind_, 1, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );
	//TR << "receive_from_node(): rec'd exp ddGbind: " << experimental_ddG_bind_ << std::endl;

	Size n_mutations;
	MPI_Recv( & n_mutations, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	mutations_.resize( n_mutations );
	//TR << "receive_from_node(): rec'd n_mutations: " << n_mutations << std::endl;

	for ( Size ii = 1; ii <= mutations_.size(); ++ii ) {
		Size position( 0 );
		MPI_Recv( & position, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
		//TR << "receive_from_node(): rec'd position: " << position << std::endl;
		int wt_aa(0), mut_aa(0);
		MPI_Recv( & wt_aa, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
		MPI_Recv( & mut_aa, 1, MPI_INT, source_node, tag, MPI_COMM_WORLD, &stat );
		//TR << "receive_from_node(): rec'd wt_aa and mut_aa: " << wt_aa << ", " << mut_aa << std::endl;
		core::chemical::AA wt_aa_  = static_cast< core::chemical::AA > ( wt_aa );
		core::chemical::AA mut_aa_ = static_cast< core::chemical::AA > ( mut_aa );
		mutations_[ ii ] = std::make_pair( position, std::make_pair( wt_aa_, mut_aa_ ) );
	}

	/// 2a. n wt complexes
	Size n_wt_complexes( 0 );
	MPI_Recv( & n_wt_complexes, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	//TR << "receive_from_node(): rec'd n_wt_complexes: " << n_wt_complexes << std::endl;

	/// 2b. n mut complexes
	Size n_mut_complexes( 0 );
	MPI_Recv( & n_mut_complexes, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	//TR << "receive_from_node(): rec'd n_mut_complexes: " << n_mut_complexes << std::endl;

	/// 2c. n wt unbounds
	Size n_wt_unbounds( 0 );
	MPI_Recv( & n_wt_unbounds, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	//TR << "receive_from_node(): rec'd n_wt_unbounds: " << n_wt_unbounds << std::endl;

	/// 2d. n mut unbounds
	Size n_mut_unbounds( 0 );
	MPI_Recv( & n_mut_unbounds, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	//TR << "receive_from_node(): rec'd n_mut_unbounds: " << n_mut_unbounds << std::endl;

	if ( n_wt_complexes == 0 || n_mut_complexes == 0 || n_wt_unbounds == 0 || n_mut_unbounds == 0 ) return;
	wt_complexes_.reserve( n_wt_complexes );
	mutant_complexes_.reserve( n_mut_complexes );
	wt_unbounds_.reserve( n_wt_unbounds );
	mutant_unbounds_.reserve( n_mut_unbounds );

	/// 3. n free weights
	Size n_free( 0 );
	MPI_Recv( & n_free, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	//TR << "receive_from_node(): rec'd n_free: " << n_free << std::endl;

	/// 4. n fixed weights
	Size n_fixed( 0 );
	MPI_Recv( & n_fixed, 1, MPI_UNSIGNED_LONG, source_node, tag, MPI_COMM_WORLD, &stat );
	//TR << "receive_from_node(): rec'd n_fixed: " << n_fixed << std::endl;


	/// --------- wt complexes ---------
	Real * free_data = new Real[ n_free * n_wt_complexes ];
	Real * fixed_data = new Real[ n_fixed * n_wt_complexes ];

	/// 5. wt complexes free data
	MPI_Recv( free_data, n_wt_complexes * n_free, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 6. wt complexes fixed data
	MPI_Recv( fixed_data, n_wt_complexes * n_fixed, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	utility::vector1< Real > free_data_v( n_free );
	utility::vector1< Real > fixed_data_v( n_fixed );
	for ( Size ii = 1; ii <= n_wt_complexes; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data_v[ jj ] = free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data_v[ jj ] = fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ];
		}
		wt_complexes_.push_back( SingleStructureDataOP( new SingleStructureData( free_data_v, fixed_data_v ) ) );
	}
	delete [] free_data; /*free_data = 0; */free_data_v.resize(0);
	delete [] fixed_data; /*fixed_data = 0; */fixed_data_v.resize(0);

	/// --------- mut complexes ---------
	free_data = new Real[ n_free * n_mut_complexes ];
	fixed_data = new Real[ n_fixed * n_mut_complexes ];

	/// 7. free data
	MPI_Recv( free_data, n_mut_complexes * n_free, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 8. fixed data
	MPI_Recv( fixed_data, n_mut_complexes * n_fixed, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	free_data_v.resize( n_free );
	fixed_data_v.resize( n_fixed );
	for ( Size ii = 1; ii <= n_mut_complexes; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data_v[ jj ] = free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data_v[ jj ] = fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ];
		}
		mutant_complexes_.push_back( SingleStructureDataOP( new SingleStructureData( free_data_v, fixed_data_v ) ) );
	}
	delete [] free_data; /*free_data = 0; */free_data_v.resize(0);
	delete [] fixed_data; /*fixed_data = 0; */fixed_data_v.resize(0);

	/// --------- wt unbounds ---------
	free_data = new Real[ n_free * n_wt_unbounds ];
	fixed_data = new Real[ n_fixed * n_wt_unbounds ];

	/// 5. wt unbounds free data
	MPI_Recv( free_data, n_wt_unbounds * n_free, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 6. wt complexes fixed data
	MPI_Recv( fixed_data, n_wt_unbounds * n_fixed, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	free_data_v.resize( n_free );
	fixed_data_v.resize( n_fixed );
	for ( Size ii = 1; ii <= n_wt_unbounds; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data_v[ jj ] = free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data_v[ jj ] = fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ];
		}
		wt_unbounds_.push_back( SingleStructureDataOP( new SingleStructureData( free_data_v, fixed_data_v ) ) );
	}
	delete [] free_data; /*free_data = 0; */free_data_v.resize(0);
	delete [] fixed_data; /*fixed_data = 0; */fixed_data_v.resize(0);

	/// --------- mut unbounds ---------
	free_data = new Real[ n_free * n_mut_unbounds ];
	fixed_data = new Real[ n_fixed * n_mut_unbounds ];

	/// 7. free data
	MPI_Recv( free_data, n_mut_unbounds * n_free, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	/// 8. fixed data
	MPI_Recv( fixed_data, n_mut_unbounds * n_fixed, MPI_DOUBLE, source_node, tag, MPI_COMM_WORLD, &stat );

	free_data_v.resize( n_free );
	fixed_data_v.resize( n_fixed );
	for ( Size ii = 1; ii <= n_mut_unbounds; ++ ii ) {
		for ( Size jj = 1; jj <= n_free; ++jj ) {
			free_data_v[ jj ] = free_data[ ( ii - 1 ) * n_free + ( jj - 1 ) ];
		}
		for ( Size jj = 1; jj <= n_fixed; ++jj ) {
			fixed_data_v[ jj ] = fixed_data[ ( ii - 1 ) * n_fixed + ( jj - 1 ) ];
		}
		mutant_unbounds_.push_back( SingleStructureDataOP( new SingleStructureData( free_data_v, fixed_data_v ) ) );
	}
	delete [] free_data; free_data = 0; free_data_v.resize(0);
	delete [] fixed_data; fixed_data = 0; fixed_data_v.resize(0);

	//std::cout << "receive_from_node(): receiving tag information from slave node." << std::endl;
	OptEPositionData::receive_from_node( source_node, tag );

}
#endif


void
DDGBindOptEData::set_experimental_ddg_bind( Real exp_ddg_bind ) {
	experimental_ddG_bind_ = exp_ddg_bind;
}


void
DDGBindOptEData::add_mutation( std::pair< Size, std::pair < AA, AA > > mutation ) {
	mutations_.push_back( mutation );
}


void
DDGBindOptEData::add_wt_complex( SingleStructureDataOP wt ) {
	wt_complexes_.push_back( wt );
}


void
DDGBindOptEData::add_mutant_complex( SingleStructureDataOP mut ) {
	mutant_complexes_.push_back( mut );
}


void
DDGBindOptEData::add_wt_unbounds( SingleStructureDataOP wt ) {
	wt_unbounds_.push_back( wt );
}


void
DDGBindOptEData::add_mutant_unbounds( SingleStructureDataOP mut ) {
	mutant_unbounds_.push_back( mut );
}


} // namespace optimize_weights
} // namespace protocols
