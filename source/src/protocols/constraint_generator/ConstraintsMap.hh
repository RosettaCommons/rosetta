// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/constraint_generator/ConstraintsMap.hh
/// @brief Cacheable data map to cache constraint pointers in the pose.
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_constraint_generator_ConstraintsMap_hh
#define INCLUDED_protocols_constraint_generator_ConstraintsMap_hh

#include <protocols/constraint_generator/ConstraintsMap.fwd.hh>

// Core headers
#include <core/scoring/constraints/Constraint.fwd.hh>

// Basic headers
#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>

#ifdef    SERIALIZATION // 1
#include <cereal/access.fwd.hpp> // 5
#include <cereal/types/polymorphic.fwd.hpp> // 9
#endif // SERIALIZATION

namespace protocols {
namespace constraint_generator {

/// @brief Cacheable data map to cache constraint pointers in the pose.
class ConstraintsMap : public basic::datacache::CacheableData {
private:
	// private types
	typedef core::scoring::constraints::ConstraintCOPs ConstraintCOPs;
	typedef std::map< std::string, ConstraintCOPs > NameToConstraintsMap;

public:
	typedef NameToConstraintsMap::iterator iterator;
	typedef NameToConstraintsMap::const_iterator const_iterator;

public:
	ConstraintsMap();

	virtual ~ConstraintsMap();

	virtual basic::datacache::CacheableDataOP
	clone() const;

public:
	/// @brief Insert csts into the ConstraintsMap under the name given.
	/// @param[in] name Map key name under which constraints will be stored
	/// @param[in] csts Constraints to store
	/// @returns ConstraintsMap::iterator to new map item.
	iterator
	insert( std::string const & name, ConstraintCOPs const & csts );

	void
	erase( iterator const & erase_me );

	iterator
	find( std::string const & name );

	const_iterator
	find( std::string const & name ) const;

	iterator
	begin();

	const_iterator
	begin() const;

	iterator
	end();

	const_iterator
	end() const;

	/// @brief prints out comma-separated list of constraint set names
	std::string
	valid_names_string() const;

private:
	/// @brief returns vector of constraint set names
	utility::vector1< std::string >
	valid_names() const;

private:
	NameToConstraintsMap cst_map_;

#ifdef   SERIALIZATION
public:
	template < class Archive >
	void save( Archive & arc ) const; // 6

	template < class Archive >
	void load( Archive & arc ); // 7
#endif // SERIALIZATION
};

} //protocols
} //constraint_generator

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_constraint_generator_ConstraintsMap ) // 8
#endif // SERIALIZATION

#endif //INCLUDED_protocols_constraint_generator_ConstraintsMap_hh

