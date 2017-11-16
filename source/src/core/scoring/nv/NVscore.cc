// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/NV/NVscore.cc
/// @details  Implementation of Neighbor Vector estimation ( from Durham EA, Et. Al. "Solvent Accessible Surface Area Approximations for Protein Structure Prediction"
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#include <core/scoring/nv/NVscore.hh>
#include <core/scoring/nv/NVscoreCreator.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/TwelveANeighborGraph.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AA.hh>

#include <core/conformation/Residue.hh>

#include <core/kinematics/Jump.hh>

#include <utility/vector1.hh>

#include <numeric/constants.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


namespace core {
namespace scoring {
namespace nv {


/// @details This must return a fresh instance of the NVscore class,
/// never an instance already in use
methods::EnergyMethodOP
NVscoreCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new NVscore );
}

ScoreTypes
NVscoreCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( neigh_count );
	sts.push_back( neigh_vect );
	sts.push_back( neigh_vect_raw );
	return sts;
}


static basic::Tracer TR( "core.scoring.NVscore" );


NVscore::NVscore() :
	parent( methods::EnergyMethodCreatorOP( new NVscoreCreator ) ),
	lookup_table_(ScoringManager::get_instance()->get_NVLookupTable() )
{
	//lbound defaults to 3.3 and ubound defaults to 11.1.  If you change these values the lookup table may no longer be accurate
	lower_bound_ = basic::options::option[ basic::options::OptionKeys::score::NV_lbound]();
	upper_bound_ = basic::options::option[ basic::options::OptionKeys::score::NV_ubound]();

	lower_bound_squared_ = lower_bound_*lower_bound_;
	upper_bound_squared_ = upper_bound_*upper_bound_;
}


methods::EnergyMethodOP NVscore::clone() const
{
	return methods::EnergyMethodOP( new NVscore(*this) );
}


void NVscore::setup_for_scoring(pose::Pose &pose, const ScoreFunction &) const
{
	pose.update_residue_neighbors();
}

void NVscore::setup_for_packing(pose::Pose &pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const
{
	pose.update_residue_neighbors();
}

void NVscore::setup_for_derivatives(pose::Pose & /*pose*/, const ScoreFunction &) const
{
	//pose.update_residue_neighbors();
}

void NVscore::setup_for_minimizing(pose::Pose & /*pose*/, ScoreFunction const & ,kinematics::MinimizerMapBase const &) const
{
	//pose.update_residue_neighbors();
}


void NVscore::indicate_required_context_graphs(utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[twelve_A_neighbor_graph] = true;
}

///Calculate the weighted neighbor count given an upper and lower bound
Real NVscore::neighbor_weight(Vector::Value const & distance2) const
{
	if ( distance2 <= lower_bound_squared_ ) {
		//neighbor count score is 1 if less than the lower bound
		return 1;
	} else if ( distance2 >= upper_bound_squared_ ) {
		//neighbor count score is 0 if greater than upper bound
		return 0 ;
	} else if ( (lower_bound_squared_ < distance2) && (upper_bound_squared_ > distance2) ) {
		//if between upper and lower bound, score follows a smooth function
		core::Real distance(sqrt(distance2));
		Real weight = ( cos( ( (distance-lower_bound_) / (upper_bound_-lower_bound_) ) * numeric::constants::r::pi ) + 1 )/2.0;
		return(weight);
	}
	return(0);
}

void NVscore::residue_energy( conformation::Residue const &current_residue,  pose::Pose const & pose, EnergyMap & emap) const
{

	Real neighbor_count(0);
	Vector neighbor_vector_sum(0,0,0);

	conformation::ResidueOPs::iterator poseIT;
	//use the coordinates of residue neighbor atom for all calcuations
	Vector current_vector(current_residue.nbr_atom_xyz());

	TwelveANeighborGraph const & graph(pose.energies().twelveA_neighbor_graph());

	for ( utility::graph::Graph::EdgeListConstIter
			node_index  = graph.get_node( current_residue.seqpos() )->const_edge_list_begin(),
			node_index_end = graph.get_node( current_residue.seqpos() )->const_edge_list_end();
			node_index != node_index_end; ++node_index ) {

		//get the residue to compare to the current residue
		core::Size comparison_residue_index((*node_index)->get_other_ind(current_residue.seqpos()));
		conformation::Residue const & comparison_residue = pose.residue(comparison_residue_index);
		//you don't want to compare a residue to itself
		if ( current_residue.seqpos() == comparison_residue.seqpos() ) continue;
		Vector const & comparison_vector(comparison_residue.nbr_atom_xyz());
		//calculate the distance between the two residues
		Vector::Value distance2 = current_vector.distance_squared(comparison_vector);
		//get the weighted neighbor count
		Real weight = neighbor_weight(distance2);

		//calculate the weighted neighbor vector for this pair and sum
		if ( weight != 0 ) {
			Vector weighted_vector = ( (comparison_vector-current_vector) / sqrt(distance2)) * weight;
			neighbor_count += weight;
			neighbor_vector_sum += weighted_vector;
		}
	}


	if ( neighbor_count == 0.0 ) return; // do not try to divide by zero

	Vector average_sum = neighbor_vector_sum/neighbor_count;
	//neighbor vector score is the norm of the average sum of all neighbor vectors
	Real neighbor_vector = average_sum.norm();

	core::chemical::AA aa_type = current_residue.aa();
	//use the neighbor vector score to look up the potential from the knowledge base
	Real nv_potential = lookup_table_.get_potentials(aa_type,neighbor_vector);

	emap[ neigh_vect ] += nv_potential;
	emap[ neigh_vect_raw ] += neighbor_vector;
	emap[ neigh_count ] += neighbor_count;
}

core::Size
NVscore::version() const
{
	return 1; // Initial versioning
}

} //NV
} //scoring
} //core
