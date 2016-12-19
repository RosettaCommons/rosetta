// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/ResidueTypeSetCache.hh
/// @brief
/// @details Internal implementation classs for ResidueTypeSet. Do not expose to the outside world.
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_chemical_ResidueTypeSetCache_HH
#define INCLUDED_core_chemical_ResidueTypeSetCache_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/chemical/ResidueTypeSetCache.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AA.hh>
#include <map>
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

typedef std::pair< AA, std::pair< utility::vector1< std::string >, utility::vector1< VariantType > > > AA_VariantsExceptions;

class ResidueTypeSetCache: public utility::pointer::ReferenceCount {

public:

	//constructor
	ResidueTypeSetCache( ResidueTypeSet const & rsd_type_set );

	//destructor
	~ResidueTypeSetCache() override;

	ResidueTypeSetCacheOP
	clone( ResidueTypeSet const & rsd_type_set ) const;

public:

	/// @details Main accessor function into ResidueTypeSetCache
	ResidueTypeCOP
	name_map( std::string const & name_in ) const;

	void
	add_residue_type( ResidueTypeCOP residue_type );

	void
	remove_residue_type( std::string const & name );

	void
	update_residue_type( ResidueTypeCOP residue_type_original, ResidueTypeCOP residue_type_new );

	bool
	has_generated_residue_type( ResidueTypeCOP rsd_type ) const;

	bool
	has_generated_residue_type( std::string const & rsd_name ) const;

	ResidueTypeCOPs
	generated_residue_types();

	void
	add_prohibited( std::string const & rsd_name );

	bool
	is_prohibited( std::string const & rsd_name ) const;

	ResidueTypeCOPs
	get_all_types_with_variants_aa( AA aa,
		utility::vector1< std::string > const & variants,
		utility::vector1< VariantType > const & exceptions );

	void clear_cached_maps();

	/// @brief information on residue types whose name3's can be changed by patches.
	std::map< std::string, std::set< std::string > > const &
	name3_generated_by_base_residue_name() {
		if ( ! cache_up_to_date_ ) { regenerate_cached_maps(); }
		return name3_generated_by_base_residue_name_;
	}

	/// @brief interchangeability groups that appear upon patch application.
	std::map< std::string, std::set< std::string > > const &
	interchangeability_group_generated_by_base_residue_name() {
		if ( ! cache_up_to_date_ ) { regenerate_cached_maps(); }
		return interchangeability_group_generated_by_base_residue_name_;
	}

private:

	void regenerate_cached_maps();

private:
	// Default constructor and copy constructor don't make sense, due to the need for a ResidueTypeSet reference: use clone() instead.
	ResidueTypeSetCache() = delete;
	ResidueTypeSetCache( ResidueTypeSetCache const & ) = delete;

private:
	ResidueTypeSet const & rsd_type_set_;

	////////////////////////////////////////////////////////////////////////////
	// Following must always be up to date.
	////////////////////////////////////////////////////////////////////////////
	/// @brief map to ResidueType pointers by unique residue id
	std::map< std::string, ResidueTypeCOP > name_map_;

	/// @brief annotation about types which theoretically may exist but we don't want
	/// @details For example, for PDB components which duplicate standard types.
	std::set< std::string > prohibited_types_;

	////////////////////////////////////////////////////////////////////////////
	// These pieces of cached data can get recomputed on the fly.
	////////////////////////////////////////////////////////////////////////////

	/// @brief map to ResidueType pointers by AA enum -- may deprecate soon.
	//std::map< AA, ResidueTypeCOPs > aa_map_;

	/// @brief map to ResidueType pointers by 3-letter string name -- may deprecate soon.
	//std::map< std::string, ResidueTypeCOPs > name3_map_;

	/// @brief caching queries based on aa & variants to avoid recomputation with ResidueTypeFinder
	///
	std::map< AA_VariantsExceptions, ResidueTypeCOPs > cached_aa_variants_map_;

	///////////////////////////////////////////////////////////////////////////
	// These need to be recomputed in a batchwise fashion,
	// as signaled by cache_up_to_date_
	//////////////////////////////////////////////////////////////////////////

	/// @brief True if the (batch) computed cached data is up to date.
	bool cache_up_to_date_ = false;

	/// @brief information on residue types whose name3's can be changed by patches.
	std::map< std::string, std::set< std::string > > name3_generated_by_base_residue_name_;

	/// @brief interchangeability groups that appear upon patch application.
	std::map< std::string, std::set< std::string > > interchangeability_group_generated_by_base_residue_name_;

#ifdef    SERIALIZATION
public:
	friend class cereal::access;
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > static void load_and_construct( Archive & arc, cereal::construct< ResidueTypeSetCache > & construct );
#endif // SERIALIZATION

};

} //chemical
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_ResidueTypeSetCache )
#endif // SERIALIZATION

#endif
