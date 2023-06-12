// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/PatchOperationFactory.hh
/// @brief Construct PatchOperations with a factory creator system.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_PatchOperationFactory_hh
#define INCLUDED_core_chemical_PatchOperationFactory_hh

// Unit Headers
#include <core/chemical/PatchOperationFactory.fwd.hh>
#include <core/chemical/PatchOperation.fwd.hh>

#include <core/types.hh>

// Utility Headers
#include <utility/SingletonBase.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/VirtualBase.hh>

// c++ headers
#include <map>
#include <set>
#include <string>

namespace core {
namespace chemical {

/// @brief The PatchOperationCreator is responsible for creating a PatchOperation from input line(s)
class PatchOperationCreator : public utility::VirtualBase
{
public:
	PatchOperationCreator() = default;

	/// @brief Return a new mover.
	virtual PatchOperationOP create_operation(std::string const & line, std::istream & input, std::map< std::string, core::Real > const & atomic_charge_reassignments) const = 0;

	/// @brief Return the tag name associated with this factory.
	virtual std::string keyname() const = 0;

};

/// @brief This templated class will register an instance of an
/// PatchOperationCreator (class T) with the PatchOperationFactory.  It will ensure
/// that no PatchOperationCreator is registered twice, and, centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class PatchOperationRegistrator : public utility::factory::WidgetRegistrator< PatchOperationFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< PatchOperationFactory, T > parent;
public:
	PatchOperationRegistrator() : parent() {}
};


class PatchOperationFactory : public utility::SingletonBase< PatchOperationFactory >
{
public:
	friend class utility::SingletonBase< PatchOperationFactory >;

	typedef std::map< std::string, PatchOperationCreatorOP > CreatorMap;

public:

	void factory_register( PatchOperationCreatorOP creator );

	/// @brief Create a mover given its identifying string, the full line, and the istream advanced to after the tag.
	PatchOperationOP newPatchOperation( std::string const & tag, std::string const & line, std::istream & input, std::map< std::string, core::Real > const & atomic_charge_reassignments ) const;

private:
	PatchOperationFactory() = default;

	// Unimplemented -- uncopyable
	PatchOperationFactory( PatchOperationFactory const & ) = delete;
	PatchOperationFactory const & operator = ( PatchOperationFactory const & ) = delete;

private:

	CreatorMap creator_map_;
};

} //namespace moves
} //namespace protocols

#endif
