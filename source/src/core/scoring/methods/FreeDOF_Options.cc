// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/methods/FreeDOF_Options.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/scoring/methods/FreeDOF_Options.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.scoring.methods.FreeDOF_Options" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace methods {

//Constructor
FreeDOF_Options::FreeDOF_Options()
{
	initialize_from_options();
}

//Destructor
FreeDOF_Options::~FreeDOF_Options()
{}

/////////////////////////////
void FreeDOF_Options::initialize_from_options() {
	free_suite_bonus_ = basic::options::option[ basic::options::OptionKeys::score::free_suite_bonus ]();
	free_2HOprime_bonus_ = basic::options::option[ basic::options::OptionKeys::score::free_2HOprime_bonus ]();
	free_sugar_bonus_ = basic::options::option[ basic::options::OptionKeys::score::free_sugar_bonus ](); // this is -1.0 by default (also ad ho
	pack_phosphate_penalty_ = basic::options::option[ basic::options::OptionKeys::score::pack_phosphate_penalty ]();
	free_side_chain_bonus_ = basic::options::option[ basic::options::OptionKeys::score::free_side_chain_bonus ]();
}

////////////////////////////
bool
operator==( FreeDOF_Options const & a, FreeDOF_Options const & b ){

	return ( ( a.free_suite_bonus_ == b.free_suite_bonus_ ) &&
		( a.free_2HOprime_bonus_ == b.free_2HOprime_bonus_ ) &&
		( a.free_sugar_bonus_ == b.free_sugar_bonus_ ) &&
		( a.pack_phosphate_penalty_ == b.pack_phosphate_penalty_ ) &&
		( a.free_side_chain_bonus_ == b.free_side_chain_bonus_ ) );
}

////////////////////////////
std::ostream &
operator<< ( std::ostream & out, const FreeDOF_Options & options ){
	options.show( out );
	return out;
}


void
FreeDOF_Options::show( std::ostream & out ) const
{
	out <<"FreeDOF_Options::show: free_suite_bonus: " << free_suite_bonus_ << std::endl;
	out <<"FreeDOF_Options::show: free_2HOprime_bonus: " << free_2HOprime_bonus_ << std::endl;
	out <<"FreeDOF_Options::show: free_sugar_bonus: " << free_sugar_bonus_ << std::endl;
	out <<"FreeDOF_Options::show: pack_phosphate_penalty: " << pack_phosphate_penalty_ << std::endl;
	out <<"FreeDOF_Options::show: free_side_chain_bonus: " << free_side_chain_bonus_ << std::endl;
}

} //methods
} //scoring
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::methods::FreeDOF_Options::save( Archive & arc ) const {
	arc( CEREAL_NVP( free_suite_bonus_ ) ); // core::Real
	arc( CEREAL_NVP( free_2HOprime_bonus_ ) ); // core::Real
	arc( CEREAL_NVP( free_sugar_bonus_ ) ); // core::Real
	arc( CEREAL_NVP( pack_phosphate_penalty_ ) ); // core::Real
	arc( CEREAL_NVP( free_side_chain_bonus_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::methods::FreeDOF_Options::load( Archive & arc ) {
	arc( free_suite_bonus_ ); // core::Real
	arc( free_2HOprime_bonus_ ); // core::Real
	arc( free_sugar_bonus_ ); // core::Real
	arc( pack_phosphate_penalty_ ); // core::Real
	arc( free_side_chain_bonus_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::methods::FreeDOF_Options );
CEREAL_REGISTER_TYPE( core::scoring::methods::FreeDOF_Options )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_methods_FreeDOF_Options )
#endif // SERIALIZATION
