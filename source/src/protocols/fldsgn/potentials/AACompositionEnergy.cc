// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fldsgn/potentials/AACompositionEnergy.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit headers
#include <protocols/fldsgn/potentials/AACompositionEnergy.hh>
#include <protocols/fldsgn/potentials/AACompositionEnergyCreator.hh>

// Package headers
#include <core/conformation/Residue.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/EnergyMethod.hh>

#include <utility>
#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace potentials {

AACompositionEnergyCreator::EnergyMethodOP
AACompositionEnergyCreator::create_energy_method( EnergyMethodOptions const & ) const
{
	return AACompositionEnergyCreator::EnergyMethodOP( new AACompositionEnergy );
}

AACompositionEnergyCreator::ScoreTypes
AACompositionEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back( core::scoring::aa_cmp );
	return sts;
}


/// @brief default constructor
AACompositionEnergy::AACompositionEnergy() :
	parent( core::scoring::methods::EnergyMethodCreatorOP( new AACompositionEnergyCreator ) )
{}


/// @brief default constructor
AACompositionEnergy::AACompositionEnergy( std::map< AA, std::pair< Real, Real > > const & comp_constraint_aas ) :
	parent( core::scoring::methods::EnergyMethodCreatorOP( new AACompositionEnergyCreator ) ),
	comp_constraint_aas_( comp_constraint_aas )
{
	initialize();
}


/// @brief copy constructor
AACompositionEnergy::AACompositionEnergy( AACompositionEnergy const & /*src*/ ) = default;


/// @brief
AACompositionEnergy::~AACompositionEnergy() = default;


/// @brief clone
AACompositionEnergy::EnergyMethodOP
AACompositionEnergy::clone() const
{
	return AACompositionEnergy::EnergyMethodOP( new AACompositionEnergy( *this ) );
}


/// @brief use GoPotential
void
AACompositionEnergy::set_comp_constraint_aa( std::map< AA, std::pair< Real, Real > > const & comp_constraint_aas )
{
	comp_constraint_aas_ = comp_constraint_aas;
	initialize();
}


/// @brief
void
AACompositionEnergy::initialize()
{}


/// @brief
void
AACompositionEnergy::residue_energy(
	Residue const & rsd,
	Pose const & pose,
	EnergyMap & emap
) const
{
	std::map< AA, Size > hist;
	Size total_types( core::chemical::num_aa_types );
	for ( Size i=1; i<=total_types; i++ ) {
		auto aa = AA( i );
		hist.insert( std::map< AA, Size >::value_type( aa, 0 ) );
	}
	for ( Size i=1; i<=pose.size(); i++ ) {
		hist[ pose.aa( i ) ] ++;
	}


	auto itr( comp_constraint_aas_.find( rsd.aa() ) );
	if ( itr != comp_constraint_aas_.end() ) {
		// histgoram if the residue at rsd.seqpos() is mutated to rsd.aa()
		hist[ itr->first ] ++;
		hist[ pose.aa( rsd.seqpos() ) ] --;
	}

	// basic energy
	Real e( 0.0 );
	itr = comp_constraint_aas_.begin();
	while ( itr != comp_constraint_aas_.end() ) {

		std::pair< Real, Real > thresholds( itr->second );
#ifndef WIN32
		auto const lower = (Size)round( thresholds.first * pose.size() );
		auto const upper = (Size)round( thresholds.second * pose.size() );
#else
		Size const lower = (Size)double( thresholds.first * pose.size() + 0.5 );
		Size const upper = (Size)double( thresholds.second * pose.size() + 0.5 );
#endif

		if ( hist[ itr->first ] < lower ) {
			e += numeric::square( hist[ itr->first ] - lower );
		} else if ( hist[ itr->first ] > upper ) {
			e += numeric::square( hist[ itr->first ] - upper );
		}

		//std::cout << core::chemical::name_from_aa( itr->first ) << " " << hist[ itr->first ] <<
		//  " " << thresholds.first << " " << thresholds.second << " " << lower << " " << upper << std::endl;
		++itr;

	}
	emap[ core::scoring::aa_cmp ] = e;

	//std::cout << rsd.seqpos() << " " << core::chemical::name_from_aa( pose.aa( rsd.seqpos() ) ) << " "
	//<< core::chemical::name_from_aa( rsd.aa() ) << " " << hist[ rsd.aa() ] << " " << e << std::endl;

}


/// @brief ReferenceEnergy is context independent; indicates that no
/// context graphs are required
void
AACompositionEnergy::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}


} // potentials
} // fldsgn
} // protocols
