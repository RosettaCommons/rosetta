// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/ResidueTypeSetCache.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_chemical_ResidueTypeSetCache_HH
#define INCLUDED_core_chemical_ResidueTypeSetCache_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/chemical/ResidueTypeSetCache.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/types.hh>
#include <map>

namespace core {
namespace chemical {

	typedef std::pair< AA, std::pair< utility::vector1< std::string >, utility::vector1< VariantType > > > AA_VariantsExceptions;

	class ResidueTypeSetCache: public utility::pointer::ReferenceCount {

	public:

		//constructor
		ResidueTypeSetCache( ResidueTypeSet const & rsd_type_set );

		//destructor
		~ResidueTypeSetCache();

	public:

		///@details Main accessor function into ResidueTypeSetCache
		ResidueType const &
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

		ResidueTypeCOPs
		get_all_types_with_variants_aa( AA aa,
																		utility::vector1< std::string > const & variants,
																		utility::vector1< VariantType > const & exceptions );

		/// @brief query ResidueTypes by their AA enum type. Does not handle aa_unk and does not handle most new patches. Legacy function.
		ResidueTypeCOPs
		aa_map_DO_NOT_USE( AA const & aa );

		/// @brief query ResidueTypes by their 3-letter name.  Does not handle aa_unk and does not handle most new patches. Legacy function.
		ResidueTypeCOPs
		name3_map_DO_NOT_USE( std::string const & name );

		void clear_cached_maps();

	private:

		ResidueTypeSet const & rsd_type_set_;

		////////////////////////////////////////////////////////////////////////////
		// Following must always be up to date.
		////////////////////////////////////////////////////////////////////////////
		/// @brief map to ResidueType pointers by unique residue id
		std::map< std::string, ResidueTypeCOP > name_map_;

		////////////////////////////////////////////////////////////////////////////
		// Following can get recomputed if custom_residue_types are
		// added to the ResidueTypeSet -- flag such an event through 'up_to_date'
		////////////////////////////////////////////////////////////////////////////
		/// @brief map to ResidueType pointers by AA enum -- may deprecate soon.
		std::map< AA, ResidueTypeCOPs > aa_map_;

		/// @brief map to ResidueType pointers by 3-letter string name -- may deprecate soon.
		std::map< std::string, ResidueTypeCOPs > name3_map_;

		/// @brief caching queries based on aa & variants to avoid recomputation with ResidueTypeFinder
		std::map< AA_VariantsExceptions, ResidueTypeCOPs > cached_aa_variants_map_;

	};

} //chemical
} //core

#endif
