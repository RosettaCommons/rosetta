// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/rna/RNA_EnergyMethodOptions.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/rna/RNA_EnergyMethodOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.scoring.rna.RNA_EnergyMethodOptions" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace rna {

//Constructor
RNA_EnergyMethodOptions::RNA_EnergyMethodOptions():
	syn_G_potential_bonus_( 0.0 ),
	torsion_potential_( "" ),
	suiteness_bonus_( "" )
{
	initialize_from_options();
}

//Destructor
RNA_EnergyMethodOptions::~RNA_EnergyMethodOptions()
{}

/////////////////////////////
void RNA_EnergyMethodOptions::initialize_from_options() {
	syn_G_potential_bonus_ = basic::options::option[ basic::options::OptionKeys::score::syn_G_potential_bonus ]();
	torsion_potential_ = basic::options::option[ basic::options::OptionKeys::score::rna_torsion_potential ]();
	suiteness_bonus_ = basic::options::option[ basic::options::OptionKeys::score::suiteness_bonus ]();
}

////////////////////////////
bool
operator==( RNA_EnergyMethodOptions const & a, RNA_EnergyMethodOptions const & b ){

	return ( ( a.syn_G_potential_bonus_ == b.syn_G_potential_bonus_ ) &&
		( a.torsion_potential_ == b.torsion_potential_  ) &&
		( a.suiteness_bonus_ == b.suiteness_bonus_ ) );
}

////////////////////////////
std::ostream &
operator<< ( std::ostream & out, const RNA_EnergyMethodOptions & options ){
	options.show( out );
	return out;
}


void
RNA_EnergyMethodOptions::show( std::ostream & out ) const
{
	out <<"RNA_EnergyMethodOptions::show: syn_G_potential_bonus: " << syn_G_potential_bonus_ << std::endl;
	out <<"RNA_EnergyMethodOptions::show: torsion_potential: " << torsion_potential_ << std::endl;
	out <<"RNA_EnergyMethodOptions::show: suiteness_bonus: " << suiteness_bonus_ << std::endl;
}

} //rna
} //scoring
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::rna::RNA_EnergyMethodOptions::save( Archive & arc ) const {
	arc( CEREAL_NVP( syn_G_potential_bonus_ ) ); // core::Real
	arc( CEREAL_NVP( torsion_potential_ ) ); // std::string
	arc( CEREAL_NVP( suiteness_bonus_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::rna::RNA_EnergyMethodOptions::load( Archive & arc ) {
	arc( syn_G_potential_bonus_ ); // core::Real
	arc( torsion_potential_ ); // std::string
	arc( suiteness_bonus_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::rna::RNA_EnergyMethodOptions );
CEREAL_REGISTER_TYPE( core::scoring::rna::RNA_EnergyMethodOptions )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_rna_RNA_EnergyMethodOptions )
#endif // SERIALIZATION
