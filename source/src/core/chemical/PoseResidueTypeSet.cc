// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core/chemical/PoseResidueTypeSet.cc
/// @brief A ResidueTypeSet which can be cached in the Pose
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/PoseResidueTypeSet.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/GlobalResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSetCache.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/Metapatch.hh>

#include <basic/Tracer.hh>

#include <utility/thread/ReadWriteMutex.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.chemical.PoseResidueTypeSet" );


namespace core {
namespace chemical {

PoseResidueTypeSet::PoseResidueTypeSet():
	ResidueTypeSet()
{}

PoseResidueTypeSet::~PoseResidueTypeSet()= default;

// We want a semi-shallow copy: that is, modifiable datamembers are cloned,
// but all const-pointers (including those to ResidueTypes) are shared.
// The base class copy constructor should give us that.
PoseResidueTypeSet::PoseResidueTypeSet( PoseResidueTypeSet const & /*src*/ ) = default;

PoseResidueTypeSet::PoseResidueTypeSet(core::chemical::ResidueTypeSetCOP deflt_rts):
	ResidueTypeSet()
{
	default_rts( deflt_rts ); // Need function call to set mode & type sets.
}

PoseResidueTypeSetOP
PoseResidueTypeSet::clone() const {
	return PoseResidueTypeSetOP( new PoseResidueTypeSet( *this ) );
}

void
PoseResidueTypeSet::default_rts(core::chemical::ResidueTypeSetCOP setting) {
	default_rts_ = setting;
	atom_type_set( setting->atom_type_set() );
	element_set( setting->element_set() );
	mm_atom_type_set( setting->mm_atom_type_set() );
	orbital_type_set( setting->orbital_type_set() );
	// Reset the all_patch_maps, and re-add this RTS's item
	all_patch_map_ = default_rts_->patch_map();
	for ( auto const & p : ResidueTypeSet::patches() ) {
		all_patch_map_[ p->name() ].push_back( p );
	}
	all_metapatch_map_ = default_rts_->metapatch_map();
	for ( auto const & p : ResidueTypeSet::metapatches() ) {
		all_metapatch_map_[ p->name() ] = p;
	}
#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard write_lock( cache_object()->read_write_mutex() );
#endif
	cache_object()->clear_cached_maps();
}

core::chemical::ResidueTypeSetCOP
PoseResidueTypeSet::default_rts() const {
	return default_rts_;
}

ResidueTypeCOP
PoseResidueTypeSet::generate_residue_type_write_locked( std::string const & name ) const {
	if ( cache_object()->has_generated_residue_type( name ) ) {
		return cache_object()->name_map( name );
	}

	std::string rsd_name_base, patch_name;
	figure_out_last_patch_from_name( name, rsd_name_base, patch_name );

	if ( patch_name.size() == 0 ) { // If this is the non-patched base type
		// The cache should already have the base types known exclusively to this PoseRTS,
		// so if we reach this point we'll just defer to the default rts
		if ( default_rts_ ) {
			ResidueTypeCOP retval( default_rts_->name_mapOP( name ) );
			cache_object()->add_pass_through( retval ); // Make grabbing from cache more efficient
			return retval;
		} else {
			return nullptr; // No default RTS means we're now stuck.
		}
	}

	ResidueTypeCOP rsd_base_ptr = generate_residue_type_write_locked( rsd_name_base );
	if ( ! rsd_base_ptr ) {
		return nullptr;
	}

	// Note: We're deliberately invoking the base class method here, as we only want the (meta)patches from this RTS itself,
	// and not any of the (meta)patches from the default_rts_ object.
	ResidueTypeCOP patched_type( apply_patch( rsd_base_ptr, patch_name, ResidueTypeSet::patch_map(), ResidueTypeSet::metapatch_map() ) );

	if ( patched_type == nullptr ) {
		// Okay, try patching with the full set of patches, but only if the base RT is a local one.
		if ( cache_object()->has_generated_residue_type( rsd_name_base ) && ! cache_object()->is_pass_through( rsd_name_base ) ) {
			patched_type = apply_patch( rsd_base_ptr, patch_name, patch_map(), metapatch_map() );
		} else if ( default_rts_ ) {
			// Otherwise, this is a purely default_rts type: return it.
			ResidueTypeCOP retval( default_rts_->name_mapOP( name ) );
			cache_object()->add_pass_through( retval ); // Make grabbing from cache more efficient
			return retval;
		}
	}

	if ( patched_type == nullptr ) {
		return nullptr;
	} else {
		cache_object()->add_residue_type( patched_type );
		return patched_type;
	}
}

ResidueTypeCOP
PoseResidueTypeSet::get_d_equivalent( ResidueTypeCOP l_rsd ) const {
	// This really should look into the local ResidueTypes and patches to see if a non-standard
	// Residue has an equivalent ... but the way we currently load the correspondences doesn't permit that.
	if ( default_rts_ ) {
		return default_rts_->get_d_equivalent( l_rsd );
	} else {
		return nullptr;
	}
}

ResidueTypeCOP
PoseResidueTypeSet::get_l_equivalent( ResidueTypeCOP d_rsd ) const {
	// This really should look into the local ResidueTypes and patches to see if a non-standard
	// Residue has an equivalent ... but the way we currently load the correspondences doesn't permit that.
	if ( default_rts_ ) {
		return default_rts_->get_l_equivalent( d_rsd );
	} else {
		return nullptr;
	}
}

ResidueTypeCOP
PoseResidueTypeSet::get_mirrored_type( ResidueTypeCOP original_rsd ) const {
	// This really should look into the local ResidueTypes and patches to see if a non-standard
	// Residue has an equivalent ... but the way we currently load the correspondences doesn't permit that.
	if ( default_rts_ ) {
		return default_rts_->get_mirrored_type( original_rsd );
	} else {
		return nullptr;
	}
}

bool
PoseResidueTypeSet::generates_patched_residue_type_with_name3( std::string const & base_residue_name, std::string const & name3 ) const {
	bool this_generates( ResidueTypeSet::generates_patched_residue_type_with_name3( base_residue_name, name3 ) );
	if ( !this_generates && default_rts_ ) {
		this_generates = default_rts_->generates_patched_residue_type_with_name3( base_residue_name, name3 );
	}
	return this_generates;
}

bool
PoseResidueTypeSet::generates_patched_residue_type_with_interchangeability_group( std::string const & base_residue_name, std::string const & interchangeability_group ) const {
	bool this_generates( ResidueTypeSet::generates_patched_residue_type_with_interchangeability_group( base_residue_name, interchangeability_group ) );
	if ( !this_generates && default_rts_ ) {
		this_generates = default_rts_->generates_patched_residue_type_with_interchangeability_group( base_residue_name, interchangeability_group );
	}
	return this_generates;
}

ResidueTypeCOPs
PoseResidueTypeSet::base_residue_types() const {
	ResidueTypeCOPs base_types( ResidueTypeSet::base_residue_types() ); // Make a copy, as we can't modify the current one.
	if ( default_rts_ ) {
		base_types.append( default_rts_->base_residue_types() );
	}
	return base_types;
}

ResidueTypeCOPs
PoseResidueTypeSet::unpatchable_residue_types() const {
	ResidueTypeCOPs unpatchable_types( ResidueTypeSet::unpatchable_residue_types() );
	if ( default_rts_ ) {
		unpatchable_types.append( default_rts_->unpatchable_residue_types() );
	}
	return unpatchable_types;
}

utility::vector1< PatchCOP >
PoseResidueTypeSet::patches() const {
	utility::vector1< PatchCOP > patch_list( ResidueTypeSet::patches() );
	if ( default_rts_ ) {
		patch_list.append( default_rts_->patches() );
	}
	return patch_list;
}

utility::vector1< MetapatchCOP >
PoseResidueTypeSet::metapatches() const {
	utility::vector1< MetapatchCOP > metapatch_list( ResidueTypeSet::metapatches() );
	if ( default_rts_ ) {
		metapatch_list.append( default_rts_->metapatches() );
	}
	return metapatch_list;
}

/// @brief Add a patch object to the RTS
void
PoseResidueTypeSet::add_patch(PatchCOP p) {
	ResidueTypeSet::add_patch(p);
	all_patch_map_[ p->name() ].push_back( p );
}

/// @brief Add a metapatch object to the RTS
void
PoseResidueTypeSet::add_metapatch(MetapatchCOP p) {
	ResidueTypeSet::add_metapatch(p);
	all_metapatch_map_[ p->name() ] = p;
}

bool
PoseResidueTypeSet::has_metapatch( std::string const & name ) const {
	return all_metapatch_map_.count( name ) != 0;
}

MetapatchCOP
PoseResidueTypeSet::metapatch( std::string const & name ) const {
	if ( all_metapatch_map_.count( name ) == 0 ) {
		utility_exit_with_message(  "Metapatch " + name + " not in the metapatch map!" );
	}
	return all_metapatch_map_.find( name )->second;
}


void
PoseResidueTypeSet::atom_type_set(AtomTypeSetCOP atom_types) {
	if ( atom_type_set() ) {
		TR.Warning << "Resetting atom type set." << std::endl;
	}
	ResidueTypeSet::atom_type_set( atom_types );
}

void
PoseResidueTypeSet::element_set(ElementSetCOP elements) {
	if ( element_set() ) {
		TR.Warning << "Resetting element type set." << std::endl;
	}
	ResidueTypeSet::element_set( elements );
}

void
PoseResidueTypeSet::mm_atom_type_set(MMAtomTypeSetCOP mm_atom_types) {
	if ( mm_atom_type_set() ) {
		TR.Warning << "Resetting MM atom type set." << std::endl;
	}
	ResidueTypeSet::mm_atom_type_set( mm_atom_types );
}

void
PoseResidueTypeSet::orbital_type_set(orbitals::OrbitalTypeSetCOP orbital_types) {
	if ( orbital_type_set() ) {
		TR.Warning << "Resetting orbital type set." << std::endl;
	}
	ResidueTypeSet::orbital_type_set( orbital_types );
}

ResidueTypeCOPs
PoseResidueTypeSet::get_all_types_with_variants_aa( AA aa, utility::vector1< std::string > const & variants ) const {
	// Note: Doesn't quite properly handle types in this set itself which mask types in the default_rts_ object with the same name.
	ResidueTypeCOPs the_types( ResidueTypeSet::get_all_types_with_variants_aa(aa, variants) );
	if ( default_rts_ ) {
		the_types.append( default_rts_->get_all_types_with_variants_aa( aa, variants ) );
	}
	return the_types;
}

ResidueTypeCOPs
PoseResidueTypeSet::get_all_types_with_variants_aa( AA aa,
	utility::vector1< std::string > const & variants,
	utility::vector1< VariantType > const & exceptions ) const
{
	// Note: Doesn't quite properly handle types in this set itself which mask types in the default_rts_ object with the same name.
	ResidueTypeCOPs the_types( ResidueTypeSet::get_all_types_with_variants_aa(aa, variants, exceptions) );
	if ( default_rts_ ) {
		the_types.append( default_rts_->get_all_types_with_variants_aa( aa, variants, exceptions ) );
	}
	return the_types;
}


// These could be hoisted by using statements, but PyRosetta currently has issues with that.
void
PoseResidueTypeSet::add_base_residue_type( ResidueTypeOP new_type ) {
	ResidueTypeSet::add_base_residue_type(new_type);
}
void
PoseResidueTypeSet::add_base_residue_type( std::string const &  filename ) {
	ResidueTypeSet::add_base_residue_type( filename );
}
void
PoseResidueTypeSet::read_files_for_base_residue_types( utility::vector1< std::string > const & filenames ) {
	ResidueTypeSet::read_files_for_base_residue_types( filenames );
}
void
PoseResidueTypeSet::add_unpatchable_residue_type( ResidueTypeOP new_type ) {
	ResidueTypeSet::add_unpatchable_residue_type(new_type);
}
void
PoseResidueTypeSet::add_unpatchable_residue_type( std::string const &  filename ) {
	ResidueTypeSet::add_unpatchable_residue_type( filename );
}
void
PoseResidueTypeSet::read_files_for_unpatchable_residue_types( utility::vector1< std::string > const & filenames ) {
	ResidueTypeSet::read_files_for_unpatchable_residue_types( filenames );
}
void
PoseResidueTypeSet::remove_base_residue_type( std::string const & name ) {
	ResidueTypeSet::remove_base_residue_type( name );
}
void
PoseResidueTypeSet::remove_unpatchable_residue_type( std::string const & name ) {
	ResidueTypeSet::remove_unpatchable_residue_type( name );
}
void
PoseResidueTypeSet::add_patches(
	utility::vector1< std::string > const & patch_filenames,
	utility::vector1< std::string > const & metapatch_filenames
) {
	ResidueTypeSet::add_patches( patch_filenames, metapatch_filenames );
}
void
PoseResidueTypeSet::set_merge_behavior_manager( MergeBehaviorManagerCOP mbm) {
	ResidueTypeSet::set_merge_behavior_manager( mbm );
}
///////////////// end hoist work-around

/// @brief Attempt to lazily load the given residue type from data.
bool
PoseResidueTypeSet::lazy_load_base_type( std::string const & rsd_base_name ) const
{
	// We don't have any special lazy loading for the PoseResidueTypeSet
	// (The default_rts will lazy load when we ask it for it's ResidueTypes.)
	return cache_object()->has_generated_residue_type( rsd_base_name );
}


} //core
} //chemical

#ifdef    SERIALIZATION

template< class Archive >
void
core::chemical::PoseResidueTypeSet::save( Archive & arc ) const {
	arc( cereal::base_class< ResidueTypeSet >( this ) );
	arc( CEREAL_NVP( default_rts_ ) ); // core::chemical::ResidueTypeSetCOP

	// Will be regenerated:
	// EXEMPT all_patch_map_
	// EXEMPT all_metapatch_map_
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::PoseResidueTypeSet::load( Archive & arc ) {
	arc( cereal::base_class< ResidueTypeSet >( this ) );

	ResidueTypeSetCOP default_rts_in;
	arc( default_rts_in ); // EXEMPT default_rts_
	default_rts( default_rts_in ); // Will regenerate all_patch_map_ and all_metapatch_map
	// We do it this way, as we don't want to make copies of (meta)patches in the GlobalRTS unnecessarily.
	// EXEMPT all_patch_map_
	// EXEMPT all_metapatch_map_
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::PoseResidueTypeSet );
CEREAL_REGISTER_TYPE( core::chemical::PoseResidueTypeSet )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_PoseResidueTypeSet )
#endif // SERIALIZATION
