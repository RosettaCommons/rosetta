// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/CacheableResidueTypeSets.hh
/// @brief A (Pose-cacheable) container for ResidueTypeSets
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_CacheableResidueTypeSets_hh
#define INCLUDED_core_chemical_CacheableResidueTypeSets_hh

#include <core/chemical/CacheableResidueTypeSets.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/PoseResidueTypeSet.fwd.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

///@brief A (Pose-cacheable) container for ResidueTypeSets
class CacheableResidueTypeSets : public basic::datacache::CacheableData {

public:

	CacheableResidueTypeSets();
	CacheableResidueTypeSets(CacheableResidueTypeSets const & src);
	CacheableResidueTypeSets & operator=(CacheableResidueTypeSets const & src);

	virtual ~CacheableResidueTypeSets();

	virtual
	basic::datacache::CacheableDataOP
	clone() const;

	void
	clear();

	/// @brief Do we have a 'mode' ResidueTypeSet already instantiated?
	bool
	has_res_type_set( TypeSetMode mode ) const;

	/// @brief Return a ResidueTypeSet of the appropriate type,
	/// If one doesn't already exist, return a null pointer.
	PoseResidueTypeSetCOP
	get_res_type_set( TypeSetMode mode = FULL_ATOM_t ) const;

	/// @brief Replace the current ResidueTypeSet of the given mode with the given RTS.
	/// If mode is INVALID_t (the recommended default) the mode will be autodetermined from the rts.
	void
	set_res_type_set( PoseResidueTypeSetOP rts, TypeSetMode mode = INVALID_t  );

private:

	utility::vector1< PoseResidueTypeSetOP > res_type_sets_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //core
} //chemical


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_CacheableResidueTypeSets )
#endif // SERIALIZATION

#endif //INCLUDED_core_chemical_CacheableResidueTypeSets_hh





