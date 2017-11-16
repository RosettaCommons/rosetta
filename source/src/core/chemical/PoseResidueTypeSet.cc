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

#include <basic/Tracer.hh>

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

PoseResidueTypeSet::~PoseResidueTypeSet(){}

// We want a semi-shallow copy: that is, modifiable datamembers are cloned,
// but all const-pointers (including those to ResidueTypes) are shared.
// The base class copy constructor should give us that.
PoseResidueTypeSet::PoseResidueTypeSet( PoseResidueTypeSet const & src ) :
	ResidueTypeSet( src ),
	default_rts_( src.default_rts_ )
{}

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
}

core::chemical::ResidueTypeSetCOP
PoseResidueTypeSet::default_rts() const {
	return default_rts_;
}

ResidueTypeCOP
PoseResidueTypeSet::name_mapOP( std::string const & name ) const {
	ResidueTypeCOP restype( ResidueTypeSet::name_mapOP( name ) );
	if ( restype ) {
		return restype;
	} else if ( default_rts_ ) {
		return default_rts_->name_mapOP( name );
	} else {
		return nullptr;
	}
}

ResidueTypeCOP
PoseResidueTypeSet::name_mapOP_write_locked( std::string const & name ) const
{
	// Begin by attempting the recursion on myself, using the parental version of this function
	ResidueTypeCOP restype( ResidueTypeSet::name_mapOP_write_locked( name ) );
	if ( restype ) {
		return restype;
	} else if ( default_rts_ ) {
		// I may be write locked, but the default_rts that I follow is not yet write locked; invoke
		// its non-write-locked version of name_mapOP
		return default_rts_->name_mapOP( name );
	} else {
		return nullptr;
	}
}


bool
PoseResidueTypeSet::has_name( std::string const & name ) const {
	if ( generate_residue_type( name ) ) {
		return true;
	} else if ( default_rts_ ) {
		return default_rts_->has_name( name );
	} else {
		return false;
	}
}

bool
PoseResidueTypeSet::has_name_write_locked( std::string const & name ) const {
	// Begin by attempting the recursion on myself, using the parental version of this function
	if ( generate_residue_type_write_locked( name ) ) {
		return true;
	} else if ( default_rts_ ) {
		// I may be write locked, but the default_rts that I follow is not yet write locked; invoke
		// its non-write-locked version of name_mapOP
		return default_rts_->has_name( name );
	} else {
		return false;
	}
}

bool
PoseResidueTypeSet::generates_patched_residue_type_with_name3( std::string const & base_residue_name, std::string const & name3 ) const {
	debug_assert( default_rts_ );
	return default_rts_->generates_patched_residue_type_with_name3( base_residue_name, name3 );
}

bool
PoseResidueTypeSet::generates_patched_residue_type_with_interchangeability_group( std::string const & base_residue_name, std::string const & interchangeability_group ) const {
	debug_assert( default_rts_ );
	return default_rts_->generates_patched_residue_type_with_interchangeability_group( base_residue_name, interchangeability_group );
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

std::map< std::string, utility::vector1< PatchCOP > >
PoseResidueTypeSet::patch_map() const {
	// We want the local patches to override the default patches, so we add in that order
	std::map< std::string, utility::vector1< PatchCOP > > patch_map;
	if ( default_rts_ ) {
		std::map< std::string, utility::vector1< PatchCOP > > dflt_patch_map( default_rts_->patch_map() );
		patch_map.insert( dflt_patch_map.begin(), dflt_patch_map.end() );
	}
	std::map< std::string, utility::vector1< PatchCOP > > parent_patch_map( ResidueTypeSet::patch_map() );
	patch_map.insert( parent_patch_map.begin(), parent_patch_map.end() );
	return patch_map;
}

bool
PoseResidueTypeSet::has_metapatch( std::string const & name ) const {
	if ( ResidueTypeSet::has_metapatch(name) ) {
		return true;
	} else if ( default_rts_ ) {
		return default_rts_->has_metapatch(name);
	} else {
		return false;
	}
}

MetapatchCOP
PoseResidueTypeSet::metapatch( std::string const & name ) const {
	if ( ResidueTypeSet::has_metapatch( name ) ) {
		return ResidueTypeSet::metapatch( name );
	} else if ( default_rts_ && default_rts_->has_metapatch( name ) ) {
		return default_rts_->metapatch( name );
	} else {
		utility_exit_with_message(  "Metapatch " + name + " not in the metapatch map!" );
	}
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
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::PoseResidueTypeSet::load( Archive & arc ) {
	arc( cereal::base_class< ResidueTypeSet >( this ) );
	arc( default_rts_ ); // core::chemical::ResidueTypeSetCOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::PoseResidueTypeSet );
CEREAL_REGISTER_TYPE( core::chemical::PoseResidueTypeSet )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_PoseResidueTypeSet )
#endif // SERIALIZATION
