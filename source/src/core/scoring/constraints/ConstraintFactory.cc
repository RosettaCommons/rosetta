// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/constraints/ConstraintFactory.cc
/// @brief Factory for creating various types of constraints.
/// @author Greg Taylor <gktaylor@u.washington.edu>

// Unit headers
#include <core/scoring/constraints/ConstraintFactory.hh>

// Package headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintCreator.hh>

#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// Singleton instance and mutex static data members
namespace utility {

using core::scoring::constraints::ConstraintFactory;

#if defined MULTI_THREADED
template <> std::mutex utility::SingletonBase< ConstraintFactory >::singleton_mutex_{};
template <> std::atomic< ConstraintFactory * > utility::SingletonBase< ConstraintFactory >::instance_( 0 );
#else
template <> ConstraintFactory * utility::SingletonBase< ConstraintFactory >::instance_( 0 );
#endif

}

namespace core {
namespace scoring {
namespace constraints {

ConstraintFactory *
ConstraintFactory::create_singleton_instance()
{
	return new ConstraintFactory;
}

/// @details Private constructor insures correctness of singleton.
ConstraintFactory::ConstraintFactory() {}

/*
void
ConstraintFactory::add_type( ConstraintOP new_func ) {
std::string type_name = new_func->type();
if ( type_name == "UNKNOWN_TYPE" ) {
utility_exit_with_message(
"failed to register Constraint... define method type() for your constraint!"
);
}
cst_types_[ type_name ] = new_func;
}

void
ConstraintFactory::add_type( std::string const & type_name, ConstraintOP new_func ) {
cst_types_[ type_name ] = new_func;
}*/

void
ConstraintFactory::factory_register( ConstraintCreatorCOP creator )
{
	cst_types_[ creator->keyname() ] = creator;
}

ConstraintOP
ConstraintFactory::newConstraint( std::string const & type_name )
{
	ConstraintCreatorMap::const_iterator iter = cst_types_.find( type_name );
	if ( iter != cst_types_.end() ) {
		return iter->second->create_constraint();
	} else {

		using std::string;
		using utility::vector1;
		string msg("ConstraintFactory::newConstraint:  ");
		msg += type_name + " does not name a known ConstraintType --> " +
			"check spelling or register new Constraint type in ConstraintFactory!";

		msg += "known types are:\n";
		vector1< string > types = get_cst_names();
		for ( vector1< string >::const_iterator it = types.begin(), end = types.end(); it != end; ++it ) {
			msg += *it + "\n";
		}

		utility_exit_with_message( msg );
	}
	return 0;
}

utility::vector1< std::string >
ConstraintFactory::get_cst_names() const {
	using std::string;
	using utility::vector1;

	vector1< string > cst_names;
	for ( ConstraintCreatorMap::const_iterator
			it = cst_types_.begin(), end = cst_types_.end(); it != end; ++it ) {
		cst_names.push_back( it->first );
	}

	return cst_names;
}

/// @details WARNING WARNING WARNING NOT THREADSAFE!
void ConstraintFactory::replace_creator( ConstraintCreatorCOP creator )
{
	cst_types_[ creator->keyname() ] = creator;
}

ConstraintCreatorCOP
ConstraintFactory::get_creator( std::string const & type_name )
{
	ConstraintCreatorMap::const_iterator iter = cst_types_.find( type_name );
	if ( iter != cst_types_.end() ) {
		return iter->second;
	} else {

		using std::string;
		using utility::vector1;
		string msg("ConstraintFactory::get_creator:  ");
		msg += type_name + " does not name a known ConstraintType --> " +
			"check spelling or register new Constraint type in ConstraintFactory!";

		msg += "known types are:\n";
		vector1< string > types = get_cst_names();
		for ( vector1< string >::const_iterator it = types.begin(), end = types.end(); it != end; ++it ) {
			msg += *it + "\n";
		}

		utility_exit_with_message( msg );
	}
	return 0;
}


/*ConstraintFactory::ConstraintFactory(void) {
// initialization of functions which this factory knows how to instantiate
add_type( new AtomPairConstraint( id::AtomID(), id::AtomID(), NULL ) );
add_type( new AngleConstraint( id::AtomID(), id::AtomID(), id::AtomID(), NULL ) );
add_type( new DihedralConstraint( id::AtomID(), id::AtomID(), id::AtomID(), id::AtomID(), NULL ) );
add_type( new InterfaceConstraint( id::AtomID(), NULL ) );
add_type( new BindingSiteConstraint() );
add_type( new BindingSiteConstraintResidues() );
add_type( new BigBinConstraint() );
add_type( new MultiConstraint() );
add_type( new AmbiguousConstraint() );
add_type( new KofNConstraint() );
add_type( new CoordinateConstraint() );
add_type( new LocalCoordinateConstraint() );
add_type( new DunbrackConstraint() );
add_type( new AmbiguousNMRDistanceConstraint() );
add_type( new AmbiguousNMRConstraint() );
add_type( new ResidueTypeConstraint() );

add_type( new AmbiguousConstraint() );
add_type( new AmbiguousNMRConstraint() );
add_type( new AmbiguousNMRDistanceConstraint() );
add_type( new AngleConstraint( id::AtomID(), id::AtomID(), id::AtomID(), NULL ) );
add_type( new AtomPairConstraint( id::AtomID(), id::AtomID(), NULL ) );
add_type( new BigBinConstraint() );
add_type( new BindingSiteConstraint() );
add_type( new BindingSiteConstraintResidues() );
add_type( new CoordinateConstraint() );
add_type( new DihedralConstraint( id::AtomID(), id::AtomID(), id::AtomID(), id::AtomID(), NULL ) );
add_type( new DunbrackConstraint() );
add_type( new KofNConstraint() );
add_type( new LocalCoordinateConstraint() );
add_type( new MultiConstraint() );


}*/

} // constraints
} // scoring
} // core
