// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ReferenceEnergyNoncanonical.cc
/// @brief  Reference energy method implementation
/// @author Fang-Chieh Chou (fcchou@stanford.edu)

// Unit headers
#include <core/scoring/methods/ReferenceEnergyNoncanonical.hh>
#include <core/scoring/methods/ReferenceEnergyNoncanonicalCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the ReferenceEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ReferenceEnergyNoncanonicalCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	if ( options.has_method_weights( ref_nc ) ) {
		return methods::EnergyMethodOP( new ReferenceEnergyNoncanonical( options.method_weights( ref_nc ) ) );
	} else {
		return methods::EnergyMethodOP( new ReferenceEnergyNoncanonical );
	}
}

ScoreTypes
ReferenceEnergyNoncanonicalCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( ref_nc );
	return sts;
}


ReferenceEnergyNoncanonical::ReferenceEnergyNoncanonical() :
	parent( EnergyMethodCreatorOP( new ReferenceEnergyNoncanonicalCreator ) )
{
	init_res_list();
}

ReferenceEnergyNoncanonical::ReferenceEnergyNoncanonical( utility::vector1< Real > const & weight_list ):
	parent( EnergyMethodCreatorOP( new ReferenceEnergyNoncanonicalCreator ) ),
	weights_( weight_list )
{
	init_res_list();
}


ReferenceEnergyNoncanonical::~ReferenceEnergyNoncanonical() {}

EnergyMethodOP
ReferenceEnergyNoncanonical::clone() const
{
	return EnergyMethodOP( new ReferenceEnergyNoncanonical( weights_ ) );
}

/// This is a terrible terrible terrible hack that will do for now.
void
ReferenceEnergyNoncanonical::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{
	// ignore scoring residues which have been marked as "REPLONLY" residues (only the repulsive energy will be calculated)
	if ( rsd.has_variant_type( core::chemical::REPLONLY ) ) return;

	std::string const res_name = rsd.name3();
	Size index_find = 0;
	for ( Size i = 1; i <= res_list_.size(); ++i ) {
		if ( res_name == res_list_[i] ) {
			index_find = i;
			break;
		}
	}

	if ( index_find == 0 ) return;

	if ( ( ! weights_.empty() ) && index_find <= weights_.size() ) {
		emap[ ref_nc ] += weights_[ index_find ];
		return;
	}

	//AMW TODO: this assigns ref style reference energies to beta AAs.
	switch ( index_find ) {
	case 1 :  emap[ ref_nc ] +=   0.16; break;
	case 2 :  emap[ ref_nc ] +=   1.70; break;
	case 3 :  emap[ ref_nc ] +=  -0.67; break;
	case 4 :  emap[ ref_nc ] +=  -0.81; break;
	case 5 :  emap[ ref_nc ] +=   0.63; break;
	case 6 :  emap[ ref_nc ] +=  -0.17; break;
	case 7 :  emap[ ref_nc ] +=   0.56; break;
	case 8 :  emap[ ref_nc ] +=   0.24; break;
	case 9 :  emap[ ref_nc ] +=  -0.65; break;
	case 10 : emap[ ref_nc ] +=  -0.10; break;
	case 11 : emap[ ref_nc ] +=  -0.34; break;
	case 12 : emap[ ref_nc ] +=  -0.89; break;
	case 13 : emap[ ref_nc ] +=   0.02; break;
	case 14 : emap[ ref_nc ] +=  -0.97; break;
	case 15 : emap[ ref_nc ] +=  -0.98; break;
	case 16 : emap[ ref_nc ] +=  -0.37; break;
	case 17 : emap[ ref_nc ] +=  -0.27; break;
	case 18 : emap[ ref_nc ] +=   0.29; break;
	case 19 : emap[ ref_nc ] +=   0.91; break;
	case 20 : emap[ ref_nc ] +=   0.51; break;
	default :
		return;
	}
}


Real
ReferenceEnergyNoncanonical::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap const &
) const
{
	return 0.0;
}

void
ReferenceEnergyNoncanonical::init_res_list()
{
	if ( !res_list_.empty() ) return;
	//Must all be 3-character string!
	std::string const res_list [] = {"B3A", "B3C", "B3D", "B3E", "B3F", "B3G", "B3H", "B3I", "B3K", "B3L", "B3M", "B3N",
		"B3P", "B3Q", "B3R", "B3S", "B3T", "B3V", "B3W", "B3Y"};
	Size const res_list_size = sizeof( res_list ) / sizeof( res_list [0] );
	for ( Size i = 0; i != res_list_size; ++i ) {
		res_list_.push_back(res_list[i]);
	}
}


/// @brief ReferenceEnergy is context independent; indicates that no
/// context graphs are required
void
ReferenceEnergyNoncanonical::indicate_required_context_graphs( utility::vector1< bool > & ) const
{}
core::Size
ReferenceEnergyNoncanonical::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

