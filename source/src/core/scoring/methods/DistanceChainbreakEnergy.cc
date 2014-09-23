// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.cc
/// @brief  Atom pair energy functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)


// Unit headers
#include <core/scoring/methods/DistanceChainbreakEnergy.hh>
#include <core/scoring/methods/DistanceChainbreakEnergyCreator.hh>

// Package headers


// Project headers
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/EnergyMap.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>

// Numeric headers

// utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "core.scoring.DistanceChainbreak", basic::t_info );

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the DistanceChainbreakEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
DistanceChainbreakEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new DistanceChainbreakEnergy;
}

ScoreTypes
DistanceChainbreakEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( distance_chainbreak );
	return sts;
}



DistanceChainbreakEnergy::DistanceChainbreakEnergy() :
	parent( methods::EnergyMethodCreatorOP( new DistanceChainbreakEnergyCreator ) )
{}

/// called at the end of energy evaluation
void
DistanceChainbreakEnergy::finalize_total_energy(
  pose::Pose & pose,
  ScoreFunction const &,
  EnergyMap & totals
) const
{
  using conformation::Residue;
  Real total_dev(0.0);
	Real const dist_target( 1.32 ); //square root of dist2_target from r++ jumping_loops.cc
	tr.Trace << "called! cuts: " << pose.fold_tree().num_cutpoint() << std::endl;
  for ( int n=1; n<= pose.fold_tree().num_cutpoint(); ++n ) {
    int const cutpoint( pose.fold_tree().cutpoint( n ) );
    Residue const & lower_rsd( pose.residue( cutpoint ) );
    if ( !lower_rsd.has_variant_type( chemical::CUTPOINT_LOWER ) ) continue;
		tr.Trace << "cutpoint " << n << "has CUTPOINT variant" << std::endl;
    Residue const & upper_rsd( pose.residue( cutpoint+1 ) );
    assert( upper_rsd.has_variant_type( chemical::CUTPOINT_UPPER ) );
		//    Size const nbb( lower_rsd.mainchain_atoms().size() );

		total_dev +=
			std::abs( dist_target - ( upper_rsd.xyz( upper_rsd.lower_connect_atom() ).distance( lower_rsd.xyz( lower_rsd.upper_connect_atom() ) ) ) );

		//r++ does it a little differently
		//r++ distancechainbreak is calculated as sqrt( dist_target**2 - dist**2 )
		// not sure if that is that is the right way to do it.
		//Distance chain break should only be the distance between 'C' of cutpoint and
		//'N' of cutpoint + 1
	}

  assert( std::abs( totals[ distance_chainbreak ] ) < 1e-3 );
  totals[ distance_chainbreak ] = total_dev;
}


/// @brief DistanceChainbreak Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
DistanceChainbreakEnergy::indicate_required_context_graphs(
   utility::vector1< bool > & /*context_graphs_required*/
) const
{}
core::Size
DistanceChainbreakEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace methods
} // namespace scoring
} // namespace core
