// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/energy_methods/ReferenceEnergyNoncanonical.cc
/// @brief  Reference energy method implementation
/// @author Fang-Chieh Chou (fcchou@stanford.edu)

// Unit headers
#include <core/energy_methods/ReferenceEnergyNoncanonical.hh>
#include <core/energy_methods/ReferenceEnergyNoncanonicalCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/conformation/Residue.hh>

#include <utility>
#include <utility/vector1.hh>


namespace core {
namespace energy_methods {



/// @details This must return a fresh instance of the ReferenceEnergy class,
/// never an instance already in use
core::scoring::methods::EnergyMethodOP
ReferenceEnergyNoncanonicalCreator::create_energy_method(
	core::scoring::methods::EnergyMethodOptions const & options
) const {
	if ( options.has_method_weights( core::scoring::ref_nc ) ) {
		return utility::pointer::make_shared< ReferenceEnergyNoncanonical >( options.method_weights( core::scoring::ref_nc ) );
	} else {
		return utility::pointer::make_shared< ReferenceEnergyNoncanonical >();
	}
}

core::scoring::ScoreTypes
ReferenceEnergyNoncanonicalCreator::score_types_for_method() const {
	using namespace core::scoring;
	ScoreTypes sts;
	sts.push_back( ref_nc );
	return sts;
}


ReferenceEnergyNoncanonical::ReferenceEnergyNoncanonical() :
	parent( utility::pointer::make_shared< ReferenceEnergyNoncanonicalCreator >() )
{
	init_res_list();
}

ReferenceEnergyNoncanonical::ReferenceEnergyNoncanonical( utility::vector1< Real > const & weight_list ):
	parent( utility::pointer::make_shared< ReferenceEnergyNoncanonicalCreator >() ),
	weights_( weight_list )
{
	init_res_list();
}


ReferenceEnergyNoncanonical::~ReferenceEnergyNoncanonical() = default;

core::scoring::methods::EnergyMethodOP
ReferenceEnergyNoncanonical::clone() const
{
	return utility::pointer::make_shared< ReferenceEnergyNoncanonical >( weights_ );
}

/// This is a terrible terrible terrible hack that will do for now.
void
ReferenceEnergyNoncanonical::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	core::scoring::EnergyMap & emap
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

	if ( index_find == 0 ) {
		// try if restype contains reference weight info as a numeric property
		// grab directly from full list; get_numeric_property("REFERENCE") function doesn't work
		std::map< std::string, core::Real > const &
			numeric_properties = rsd.type().properties().numeric_properties();

		auto property_it( numeric_properties.find( "REFERENCE" ) );
		if ( property_it == numeric_properties.end() ) {
			return;
		} else {
			emap[ core::scoring::ref_nc ] += property_it->second;
		}
		return; // skip below
	}

	if ( ( ! weights_.empty() ) && index_find <= weights_.size() ) {
		emap[ core::scoring::ref_nc ] += weights_[ index_find ];
		return;
	}

	//AMW TODO: this assigns ref style reference energies to beta AAs.
	switch ( index_find ) {
	case 1 :  emap[ core::scoring::ref_nc ] +=   0.16; break;
	case 2 :  emap[ core::scoring::ref_nc ] +=   1.70; break;
	case 3 :  emap[ core::scoring::ref_nc ] +=  -0.67; break;
	case 4 :  emap[ core::scoring::ref_nc ] +=  -0.81; break;
	case 5 :  emap[ core::scoring::ref_nc ] +=   0.63; break;
	case 6 :  emap[ core::scoring::ref_nc ] +=  -0.17; break;
	case 7 :  emap[ core::scoring::ref_nc ] +=   0.56; break;
	case 8 :  emap[ core::scoring::ref_nc ] +=   0.24; break;
	case 9 :  emap[ core::scoring::ref_nc ] +=  -0.65; break;
	case 10 : emap[ core::scoring::ref_nc ] +=  -0.10; break;
	case 11 : emap[ core::scoring::ref_nc ] +=  -0.34; break;
	case 12 : emap[ core::scoring::ref_nc ] +=  -0.89; break;
	case 13 : emap[ core::scoring::ref_nc ] +=   0.02; break;
	case 14 : emap[ core::scoring::ref_nc ] +=  -0.97; break;
	case 15 : emap[ core::scoring::ref_nc ] +=  -0.98; break;
	case 16 : emap[ core::scoring::ref_nc ] +=  -0.37; break;
	case 17 : emap[ core::scoring::ref_nc ] +=  -0.27; break;
	case 18 : emap[ core::scoring::ref_nc ] +=   0.29; break;
	case 19 : emap[ core::scoring::ref_nc ] +=   0.91; break;
	case 20 : emap[ core::scoring::ref_nc ] +=   0.51; break;
		// "OMA", "OMC", "OMG", "OMU", "PSU", "H2U", "PUR", "1MA"
		// GACU 4.14  3.58  2.82  3.76
	case 21 : emap[ core::scoring::ref_nc ] +=   3.58; break;
	case 22 : emap[ core::scoring::ref_nc ] +=   2.82; break;
	case 23 : emap[ core::scoring::ref_nc ] +=   4.14; break;
	case 24 : emap[ core::scoring::ref_nc ] +=   3.76; break;
	case 25 : emap[ core::scoring::ref_nc ] +=   3.76; break;
	case 26 : emap[ core::scoring::ref_nc ] +=   3.76; break;
	case 27 : emap[ core::scoring::ref_nc ] +=   3.58; break;
	case 28 : emap[ core::scoring::ref_nc ] +=   3.58; break;
	default :
		return;
	}
}


Real
ReferenceEnergyNoncanonical::eval_dof_derivative(
	id::DOF_ID const &,
	id::TorsionID const &,
	pose::Pose const &,
	core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap const &
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
		"B3P", "B3Q", "B3R", "B3S", "B3T", "B3V", "B3W", "B3Y", "OMA", "OMC", "OMG", "OMU", "PSU", "H2U", "PUR", "1MA" };
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


} // scoring
} // core

