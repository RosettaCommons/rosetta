// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/constraint_generator/ConstraintGeneratorFactory.cc
/// @brief  Implementation of the class for instantiating arbitrary ConstraintGenerators
///         from a string --> ConstraintGeneratorCreator map
/// @author Tom Linsky ( tlinsky at uw dot edu )

// Unit headers
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>

// Package headers
#include <protocols/constraint_generator/ConstraintGenerator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

namespace protocols {
namespace constraint_generator {

ConstraintGeneratorFactory::ConstraintGeneratorFactory():
	utility::SingletonBase< ConstraintGeneratorFactory >(),
	creator_map_()
{}

ConstraintGeneratorFactory *
ConstraintGeneratorFactory::create_singleton_instance()
{
	return new ConstraintGeneratorFactory;
}

void
ConstraintGeneratorFactory::factory_register( ConstraintGeneratorCreatorOP creator )
{
	if ( creator_map_.find( creator->keyname() ) != creator_map_.end() ) {
		std::string const err_msg = "Factory Name Conflict: Two or more ConstraintGeneratorCreators registered with the name " + creator->keyname();
		utility_exit_with_message(  err_msg );
	}
	creator_map_[ creator->keyname() ] = creator;
}

bool ConstraintGeneratorFactory::has_type( std::string const & selector_type ) const
{
	return creator_map_.find( selector_type ) != creator_map_.end();
}

ConstraintGeneratorOP
ConstraintGeneratorFactory::new_constraint_generator(
	std::string const & constraint_generator_name,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const
{
	if ( ! has_type( constraint_generator_name ) ) {
		std::string err_msg =  "No ConstraintGeneratorCreator with the name '" + constraint_generator_name + "' has been registered with the ConstraintGeneratorFactory";
		throw utility::excn::EXCN_Msg_Exception( err_msg );
	}
	auto iter = creator_map_.find( constraint_generator_name );
	ConstraintGeneratorOP new_constraint_generator = iter->second->create_constraint_generator();
	new_constraint_generator->parse_my_tag( tag, datamap );
	return new_constraint_generator;
}

} //namespace constraint_generator
} //namespace protocols

// Singleton instance and mutex static data members
namespace utility {

using protocols::constraint_generator::ConstraintGeneratorFactory;

#ifdef MULTI_THREADED
template<> std::mutex utility::SingletonBase< ConstraintGeneratorFactory >::singleton_mutex_{};
template<> std::atomic< ConstraintGeneratorFactory * > utility::SingletonBase< ConstraintGeneratorFactory >::instance_( NULL );
#else
template<> ConstraintGeneratorFactory * utility::SingletonBase< ConstraintGeneratorFactory >::instance_( nullptr );
#endif

} // namespace utility
