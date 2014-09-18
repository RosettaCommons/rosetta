// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/IntermolEnergy.cc
/// @brief  Cost of bringing two chains together.
/// @author Rhiju Das


// Unit headers
#include <core/scoring/methods/IntermolEnergy.hh>
#include <core/scoring/methods/IntermolEnergyCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/scoring/EnergyMap.hh>

// Utility headers
#include <utility/vector1.hh>

// C++
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "core.scoring.methods.IntermolEnergy" );

/////////////////////////////////////////////////////////////////////////////////////
//
// Created in attempt to fit 'Turner rules' for METHODS.
//  assumes 1 M standard state -- later will allow
//  specification of strand concentration(s) from
//  command line.
//
/////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the IntermolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
IntermolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new IntermolEnergy;
}

ScoreTypes
IntermolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( intermol );
	return sts;
}


/// c-tor
IntermolEnergy::IntermolEnergy() :
	parent( new IntermolEnergyCreator ),
	penalty_at_1M_( 2.30 ), // calibrated outside.
	log_conc_( std::log( basic::options::option[ basic::options::OptionKeys::score::conc ]() ) ) // in kT
{}

/// clone
methods::EnergyMethodOP
IntermolEnergy::clone() const
{
	return new IntermolEnergy;
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
void
IntermolEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {

	using namespace core::pose::full_model_info;
	make_sure_full_model_info_is_setup( pose );
	Size const num_chains_frozen = get_num_chains_frozen( pose );
	totals[ intermol ] = num_chains_frozen * ( penalty_at_1M_ - log_conc_ );
} // finalize_total_energy


///////////////////////////////////////////////////////////////////////////////
// quick graph traversal lets us figure out number of connected components.
// this may repeat work in LoopGraph and LoopClose Energy, in which case unify!
Size
IntermolEnergy::get_num_chains_frozen( pose::Pose const & pose ) const {
	using namespace core::pose::full_model_info;
	Size const num_chains = const_full_model_info( pose ).cutpoint_open_in_full_model().size() + 1;

	utility::vector1< utility::vector1< Size > > chains_connected;
	utility::vector1< Size > blank_vector;
	for ( Size n = 1; n <= num_chains; n++ ) blank_vector.push_back( false );
	for ( Size n = 1; n <= num_chains; n++ ) chains_connected.push_back( blank_vector );

	get_chains_connected( pose, chains_connected );

	// alternatively, there must be a graph/ sublibrary somewhere that does this.
	Size const num_subgraphs = get_number_of_connected_subgraphs( chains_connected );
	runtime_assert( num_subgraphs <= num_chains );
	return ( num_chains - num_subgraphs );
}

///////////////////////////////////////////////////////////////////////////////
Size
IntermolEnergy::get_number_of_connected_subgraphs( utility::vector1< utility::vector1< Size > > const & chains_connected ) const {
	Size const num_chains = chains_connected.size();
	utility::vector1< Size > colors( num_chains, 0 );
	Size current_color = 0;
	for ( Size n = 1; n <= num_chains; n++ ){
		if ( colors[ n ] ) continue;
		current_color++;
		colors[ n ] =  current_color;
		color_connected( colors, chains_connected, n, current_color );
	}
	return current_color;
}

///////////////////////////////////////////////////////////////////////////////
void
IntermolEnergy::color_connected( utility::vector1< Size > & colors,
																 utility::vector1< utility::vector1< Size > > const & chains_connected,
																 Size const n, Size const current_color ) const {
	Size const num_chains = chains_connected.size();
	for ( Size m = (n+1); m <= num_chains; m++ ){
		if ( chains_connected[ n ][ m ] && !colors[ m ] ){
			colors[ m ] = current_color;
			color_connected( colors, chains_connected, m, current_color );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void
IntermolEnergy::get_chains_connected( pose::Pose const & pose, utility::vector1< utility::vector1< Size > > &  chains_connected ) const {
	using namespace core::pose::full_model_info;
	utility::vector1< Size > const chains = figure_out_chain_numbers_from_full_model_info_const( pose );
	utility::vector1< Size > frozen_chains;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		if ( !frozen_chains.has_value( chains[ n ] ) ) frozen_chains.push_back( chains[ n ] );
	}
	for ( Size i = 1; i <= frozen_chains.size(); i++ ){
		for ( Size j = 1; j <= frozen_chains.size(); j++ ){
			chains_connected[ i ][ j ] = true;
			chains_connected[ j ][ i ] = true;
		}
	}
	// recurse through any daughter poses.
	utility::vector1< pose::PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();
	for ( Size k = 1; k <= other_pose_list.size(); k++ )	get_chains_connected( *other_pose_list[ k ], chains_connected );
}

///////////////////////////////////////////////////////////////////////////////
void
IntermolEnergy::eval_atom_derivative(
	id::AtomID const &,
	pose::Pose const &,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const &,
	Vector &,
	Vector &
 	) const
{
	// no op.
} // eval atom derivative

core::Size
IntermolEnergy::version() const
{
	return 1; // Initial versioning
}



} // methods
} // scoring
} // core
