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

// Package headers
#include <core/chemical/ResidueTypeSetCache.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AA.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>
#include <set>
#include <tuple>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.fwd.hh>
#endif

namespace core {
namespace chemical {

typedef std::pair< AA, std::pair< utility::vector1< std::string >, utility::vector1< VariantType > > > AA_VariantsExceptions;

/// @brief A type for key lookup.  This consists of a BASE TYPE ResidueTypeCOP, a list of variants (by enum), a list of custom variants (by string),
/// a list of variant exceptions (by enum), and a bool for whether to turn off metapatches.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
typedef std::tuple<
	ResidueTypeCOP,
	utility::vector1< VariantType >, //Included VariantTypes
	utility::vector1< std::string >, //Included special VariantTypes without enums
	utility::vector1< VariantType >, //Excluded VariantTypes
	bool //Disable metapatches
	> Restype_Variants_VariantStrings_Exceptions;

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

	/// @details Like name_map, but returns a nullptr rather than raising an error if the entry can't be found
	ResidueTypeCOP
	name_map_or_null( std::string const & name_in ) const;

	void
	add_residue_type( ResidueTypeCOP residue_type );

	/// @brief Add a ResidueType to the cache which isn't strictly in the associated ResidueTypeSet,
	/// but is included here for efficiency/convenience
	void
	add_pass_through( ResidueTypeCOP residue_type );

	/// @brief Is the ResidueType one of the pass-through convenience types?
	bool
	is_pass_through( std::string const & name_in );

	void
	remove_residue_type( std::string const & name );

	void
	update_residue_type( ResidueTypeCOP residue_type_original, ResidueTypeCOP residue_type_new );

	bool
	has_generated_residue_type( ResidueTypeCOP rsd_type ) const;

	bool
	has_generated_residue_type( std::string const & rsd_name ) const;

	//ResidueTypeCOPs
	//generated_residue_types();

	void
	add_prohibited( std::string const & rsd_name );

	bool
	is_prohibited( std::string const & rsd_name ) const;

	/// @brief Returns whether or not the all_types_with_variants_aa map already
	/// has an entry for the given combination of aa, variants, and exceptions.
	/// If so, then the cached data may be directly retrieved.
	bool
	all_types_with_variants_aa_already_cached(
		AA aa,
		utility::vector1< std::string > const & variants,
		utility::vector1< VariantType > const & exceptions ) const;

	void
	cache_all_types_with_variants_aa(
		AA aa,
		utility::vector1< std::string > const & variants,
		utility::vector1< VariantType > const & exceptions,
		ResidueTypeCOPs cached_types
	);

	ResidueTypeCOPs
	retrieve_all_types_with_variants_aa(
		AA aa,
		utility::vector1< std::string > const & variants,
		utility::vector1< VariantType > const & exceptions ) const;

	/// @brief Have we already cached this combination of base type, variants, and exceptions?
	/// @note This function is *not* threadsafe.  It is assumed that thread locking will be handled by whatever is calling this.
	/// @param[in] base_type A ResidueTypeCOP to a base residue type, used for looking up the variant.
	/// @param[in] variants A list of VariantTypes that the returned ResidueTypes *must* have, used for looking up the variant.
	/// @param[in] variant_strings A list of custom VariantTypes (that don't have enums) that the returned ResidueTypes *must*
	/// have, used for looking up the variant.
	/// @param[in] exceptions A list of VariantTypes that are ignored in matching.
	/// @param[in] no_metapatches If true, metapatches are ignored.
	/// @returns TRUE if we have cached this combination, FALSE otherwise.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	bool
	all_types_with_variants_residuetypecops_already_cached(
		ResidueTypeCOP base_type,
		utility::vector1< VariantType > const & variants,
		utility::vector1< std::string > const & variant_strings,
		utility::vector1< VariantType > const & exceptions,
		bool const no_metapatches
	) const;

	/// @brief Cache this combination of base type, variants, and exceptions.
	/// @note This function is *not* threadsafe.  It is assumed that thread locking will be handled by whatever is calling this.
	/// @param[in] base_type A ResidueTypeCOP to a base residue type, used for looking up the variant.
	/// @param[in] variants A list of VariantTypes that the returned ResidueTypes *must* have, used for looking up the variant.
	/// @param[in] variant_strings A list of custom VariantTypes (that don't have enums) that the returned ResidueTypes *must*
	/// have, used for looking up the variant.
	/// @param[in] exceptions A list of VariantTypes that are ignored in matching.
	/// @param[in] no_metapatches If true, metapatches are ignored.
	/// @param[in] types_to_cache A vector of ResidueTypeCOPs to associate with the key defined by the five variables above.
	/// @returns Nothing.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	void
	cache_all_types_with_variants_residuetypecops(
		ResidueTypeCOP base_type,
		utility::vector1< VariantType > const & variants,
		utility::vector1< std::string > const & variant_strings,
		utility::vector1< VariantType > const & exceptions,
		bool const no_metapatches,
		ResidueTypeCOPs const & types_to_cache
	);

	/// @brief Retrieve this combination of base type, variants, and exceptions.  NOTE THAT THIS ASSUMES THIS INFORMATION WAS
	/// ALREADY CACHED!
	/// @note This function is *not* threadsafe.  It is assumed that thread locking will be handled by whatever is calling this.
	/// @param[in] base_type A ResidueTypeCOP to a base residue type, used for looking up the variant.
	/// @param[in] variants A list of VariantTypes that the returned ResidueTypes *must* have, used for looking up the variant.
	/// @param[in] variant_strings A list of custom VariantTypes (that don't have enums) that the returned ResidueTypes *must*
	/// have, used for looking up the variant.
	/// @param[in] exceptions A list of VariantTypes that are ignored in matching.
	/// @param[in] no_metapatches If true, metapatches are ignored.
	/// @returns A vector of ResidueTypeCOPs with the given base type and the associated variant types.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	ResidueTypeCOPs
	retrieve_all_types_with_variants_residuetypecops(
		ResidueTypeCOP base_type,
		utility::vector1< VariantType > const & variants,
		utility::vector1< std::string > const & variant_strings,
		utility::vector1< VariantType > const & exceptions,
		bool const no_metapatches
	) const;

	//ResidueTypeCOPs
	//get_all_types_with_variants_aa( AA aa,
	// utility::vector1< std::string > const & variants,
	// utility::vector1< VariantType > const & exceptions );

	void clear_cached_maps();

	/// @brief The RTSC performs just-in-time updates on data members that are accessed
	/// through two of its methods -- before calling these methods, the ResidueTypeSet may
	/// need to obtain a write lock on the RTSC. Thse are:
	/// - name3_generated_by_base_residue_name, and
	/// - interchangeability_group_generated_by_base_residue_name
	bool
	maps_up_to_date() const {
		return cache_up_to_date_;
	}

	/// @brief information on residue types whose name3's can be changed by patches.
	/// @details To call this function, the ResidueTypeSet will need to obtain a
	/// write lock if the maps_up_to_date function returns false. Otherwise, a read lock
	/// will suffice.
	std::map< std::string, std::set< std::string > > const &
	name3_generated_by_base_residue_name();

	/// @brief interchangeability groups that appear upon patch application.
	/// @details To call this function, the ResidueTypeSet will need to obtain a
	/// write lock if the maps_up_to_date function returns false. Otherwise, a read lock
	/// will suffice.
	std::map< std::string, std::set< std::string > > const &
	interchangeability_group_generated_by_base_residue_name();

	/// @brief Does this cache already have a residue with the given name3
	/// This is likely to be slow, as it looks at all the generated residue types.
	bool
	has_restype_with_name3( std::string const & name3 ) const;

#ifdef MULTI_THREADED
	/// @brief The %ResidueTypeSet will lock this instance for read and write opperations
	/// by accessing the mutex that this RTSC stores.
	utility::thread::ReadWriteMutex & read_write_mutex();
#endif

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

	/// @brief ResidueTypes which don't directly belong to the associated ResidueTypeSet,
	/// but are included in the cache for efficiency.
	std::set< std::string > pass_throughs_;

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

	/// @brief Caching queries based on ResidueTypeCOP and variants to avoid recomputation with ResidueTypeFinder.
	/// @details The key is a tuple of < BASE residue type COP, included variant list (enum and string), excluded variant list (enum and string) >.
	/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
	std::map< Restype_Variants_VariantStrings_Exceptions, ResidueTypeCOPs  > cached_restype_variants_map_;

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

#ifdef MULTI_THREADED
	utility::thread::ReadWriteMutexOP read_write_mutex_;
	//std::recursive_mutex cache_mutex_;
#endif

};

} //chemical
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_ResidueTypeSetCache )
#endif // SERIALIZATION

#endif
