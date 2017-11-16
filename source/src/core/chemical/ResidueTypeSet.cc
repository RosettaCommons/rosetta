// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ResidueTypeSet.cc
/// @brief ResidueTypeSet class
///
/// @details
/// This class is an abstract base class for other ResidueTypeSets.
///  ( Currently this is just the GlobalResidueTypeSet and the PoseResidueTypeSet.)
///
/// @author Rocco Moretti (rmorettiase@gmail.com)
/// Phil Bradley
/// Steven Combs
/// Rhiju Das


// Rosetta headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSetCache.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/MergeBehaviorManager.hh>
#include <core/chemical/Metapatch.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/adduct_util.hh>
#include <core/chemical/util.hh>
#include <core/chemical/Orbital.hh> /* for copying ResidueType */
#include <core/chemical/ResidueConnection.hh> /* for copying ResidueType */

#include <core/chemical/gasteiger/GasteigerAtomTyper.hh>

// Basic headers
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <fstream>
#include <string>
//#include <sstream>
#include <set>
#include <algorithm>

// option key includes
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/chemical/orbitals/AssignOrbitals.hh>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>

#ifdef MULTI_THREADED
#include <utility/thread/ReadWriteMutex.hh>
#endif

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

using namespace basic::options;

namespace core {
namespace chemical {

static basic::Tracer TR( "core.chemical.ResidueTypeSet" );

ResidueTypeSet::ResidueTypeSet( TypeSetMode mode /*= INVALID_t*/ ) :
	mode_( mode ),
	merge_behavior_manager_( new MergeBehaviorManager ),
	cache_( ResidueTypeSetCacheOP( new ResidueTypeSetCache( *this ) ) )
{}

ResidueTypeSet::~ResidueTypeSet() = default;

///////////////////////////////////////////////////////////////////////////////

void
ResidueTypeSet::atom_type_set(AtomTypeSetCOP atom_types) {
	if ( ! atom_types ) {
		TR.Warning << "Setting a ResidueTypeSet's AtomTypeSet to a null pointer!" << std::endl;
	} else if ( mode() != INVALID_t && atom_types->mode() != mode() ) {
		TR.Warning << "Using an AtomTypeSet of mode " << atom_types->mode()
			<< " with a ResidueTypeSet of mode " << mode() << std::endl;
	} else if ( mode() == INVALID_t ) {
		mode_ = atom_types->mode();
	}
	atom_types_ = atom_types;
}

//////////////////////////////////////////////////////////////////
/// @details since residue id is unique, it only returns
/// one residue type or exit without match.
///
//////////////////////////////////////////////////////////////////
ResidueType const &
ResidueTypeSet::name_map( std::string const & name_in ) const
{
	std::string name = name_in;
	if ( name_in == "CYD" ) {
		name = "CYS:disulfide";
	}
	ResidueTypeCOP restype( name_mapOP( name ) );
	runtime_assert_string_msg( restype != nullptr, "The residue " + name + " could not be generated.  Has a suitable params file been loaded?  (Note that custom params files not in the Rosetta database can be loaded with the -extra_res or -extra_res_fa command-line flags.)"  );
	return *restype;
}

//////////////////////////////////////////////////////////////////
///
/// MOST RESIDUE TYPES WILL BE GENERATED THROUGH THIS FUNCTION.
///
//////////////////////////////////////////////////////////////////
ResidueTypeCOP
ResidueTypeSet::name_mapOP( std::string const & name_in ) const
{
	std::string const name = fixup_patches( name_in );
	{ // scope the read lock
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard read_lock( cache_->read_write_mutex() );
#endif
		if ( cache_object()->has_generated_residue_type( name) ) {
			return cache_object()->name_map( name );
		}
	}

#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard write_lock( cache_->read_write_mutex() );
#endif
	return name_mapOP_write_locked( name );
}

bool
ResidueTypeSet::generate_residue_type( std::string const & rsd_name ) const
{
	{ // scope the read lock
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard read_lock( cache_->read_write_mutex() );
#endif
		if ( cache_object()->has_generated_residue_type( rsd_name ) ) {
			return true;
		}
	}
#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard write_lock( cache_->read_write_mutex() );
#endif
	return generate_residue_type_write_locked( rsd_name );
}

/// @details The calling function must first obtain a write lock on the ResidueTypeSetCache
ResidueTypeCOP
ResidueTypeSet::name_mapOP_write_locked( std::string const & name_in ) const
{
	std::string const name = fixup_patches( name_in );
	if ( generate_residue_type_write_locked( name ) ) {
		return cache_object()->name_map( name );
	} else {
		return ResidueTypeCOP( nullptr );
	}
}

/// @details Instantiates ResidueType recursively, peeling off the last-most patch and generating
/// the patched ResidueType from the base ResidueType and the corresponding Patch, recursively
/// generating the base ResidueType if necessary. The function that calls this function must first
/// obtain a write lock
bool
ResidueTypeSet::generate_residue_type_write_locked( std::string const & rsd_name ) const
{
	if ( cache_object()->has_generated_residue_type( rsd_name ) ) return true; // already generated

	// get name (which holds patch information)
	std::string rsd_name_base, patch_name;
	figure_out_last_patch_from_name( rsd_name, rsd_name_base, patch_name );

	if ( patch_name.size() == 0 ) { // If this is the non-patched base type
		if ( ! cache_object()->has_generated_residue_type( rsd_name_base ) ) {
			lazy_load_base_type( rsd_name_base );
		}
		return cache_object()->has_generated_residue_type( rsd_name ); // Is generated?
	}

	// now apply patches.
	ResidueTypeCOP rsd_base_ptr = name_mapOP_write_locked( rsd_name_base );
	if ( ! rsd_base_ptr ) { return false; }

	ResidueType const & rsd_base( *rsd_base_ptr );
	runtime_assert( rsd_base.finalized() );

	// I may have to create this patch, if it is metapatch-derived!
	// This version preserves this function's constness by not
	// adding to patches_ or the patch_map.
	// You are looking for the name of the base residue type
	if ( patch_name.find( "MP-" ) != std::string::npos ) {

		std::string metapatch_name = utility::string_split( patch_name, '-' )[3];

		// This patch needs to be generated by this metapatch.
		// PHE-CD1-methylated
		// Atom name: second element if you split the patch name by '-'
		std::string atom_name = utility::string_split( patch_name, '-' )[2];

		if ( ! has_metapatch( metapatch_name ) ) {
			TR.Debug << "Metapatch " << metapatch_name << " not in the metapatch map!" << std::endl;
			return false; // Can't generate this with this ResidueTypeSet
		}

		// Add buffering spaces just in case the resultant patch is PDB-naming sensitive
		// we need enough whitespace -- will trim later
		PatchCOP needed_patch = metapatch( metapatch_name )->get_one_patch( /*rsd_base, */"  " + atom_name + "  " );

		ResidueTypeOP rsd_instantiated ( needed_patch->apply( rsd_base ) );

		if ( rsd_instantiated == nullptr ) {
			return false; // utility_exit_with_message(  "Failed to apply: " + p->name() + " to " + rsd_base.name() );
		}
		if ( option[ OptionKeys::in::file::assign_gasteiger_atom_types ] ) {
			gasteiger::GasteigerAtomTypeSetCOP gasteiger_set(
				ChemicalManager::get_instance()->gasteiger_atom_type_set() );
			gasteiger::assign_gasteiger_atom_types( *rsd_instantiated, gasteiger_set, false );
		}
		if ( option[ OptionKeys::in::add_orbitals] ) {
			orbitals::AssignOrbitals( rsd_instantiated ).assign_orbitals();
		}

		//Set the pointer to the base type:
		if ( rsd_base.get_base_type_cop() ) {
			rsd_instantiated->set_base_type_cop( rsd_base.get_base_type_cop() );
			if ( rsd_base.is_base_type() && rsd_base.mainchain_potentials_match( *rsd_instantiated ) ) { //If we're making a copy of a base type and the mainchain potentials match...
				rsd_instantiated->reset_mainchain_torsion_potential_names();
			}
		} else {
			rsd_instantiated->set_base_type_cop( rsd_base_ptr );
		}

		cache_object()->add_residue_type( rsd_instantiated );
		return true;

	} else {
		// patch_map() returns by value - need to keep a reference to the object around while we're using it.
		std::map< std::string, utility::vector1< PatchCOP > > patch_mapping( patch_map() );
		if ( patch_mapping.find( patch_name ) == patch_mapping.end() ) return false;
		utility::vector1< PatchCOP > const & patches = patch_mapping.find( patch_name )->second;
		bool patch_applied( false );

		// sometimes patch cases are split between several patches -- look through all:
		for ( PatchCOP p : patches ) {
			if ( !p->applies_to( rsd_base ) ) continue;

			runtime_assert( !patch_applied ); // patch cannot be applied twice.

			ResidueTypeOP rsd_instantiated = p->apply( rsd_base );

			if ( rsd_instantiated == nullptr ) {
				return false; // utility_exit_with_message(  "Failed to apply: " + p->name() + " to " + rsd_base.name() );
			}
			if ( option[ OptionKeys::in::file::assign_gasteiger_atom_types ] ) {
				gasteiger::GasteigerAtomTypeSetCOP gasteiger_set(
					ChemicalManager::get_instance()->gasteiger_atom_type_set() );
				gasteiger::assign_gasteiger_atom_types( *rsd_instantiated, gasteiger_set, false );
			}
			if ( option[ OptionKeys::in::add_orbitals] ) {
				orbitals::AssignOrbitals( rsd_instantiated ).assign_orbitals();
			}

			//Set the pointer to the base type:
			if ( rsd_base.get_base_type_cop() ) {
				rsd_instantiated->set_base_type_cop( rsd_base.get_base_type_cop() );
				if ( rsd_base.is_base_type() && rsd_base.mainchain_potentials_match( *rsd_instantiated ) ) { //If we're making a copy of a base type...
					rsd_instantiated->reset_mainchain_torsion_potential_names();
				}
			} else {
				rsd_instantiated->set_base_type_cop( rsd_base_ptr );
			}

			cache_object()->add_residue_type( rsd_instantiated );
			patch_applied = true;
		}
		return patch_applied;
	}
}


void
ResidueTypeSet::add_patches(
	utility::vector1< std::string > const & patch_filenames,
	utility::vector1< std::string > const & metapatch_filenames
)
{
	for ( std::string const & filename : patch_filenames ) {
		PatchOP p( new Patch(mode_) );
		p->read_file( filename );
		patches_.push_back( p );
		patch_map_[ p->name() ].push_back( p );
	}

	for ( std::string const & filename : metapatch_filenames ) {
		MetapatchOP p( new Metapatch );
		p->read_file( filename );
		metapatches_.push_back( p );
		metapatch_map_[ p->name() ] = p;
	}

#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard write_lock( cache_object()->read_write_mutex() );
#endif
	cache_object()->clear_cached_maps();
}



/// @brief Check if a base type (like "SER") generates any types with another name3 (like "SEP")
bool
ResidueTypeSet::generates_patched_residue_type_with_name3( std::string const & base_residue_name,
	std::string const & name3 ) const
{
#ifdef MULTI_THREADED
	{ // scope the read lock
		utility::thread::ReadLockGuard read_lock( cache_object()->read_write_mutex() );
		if ( cache_object()->maps_up_to_date() ) {
			std::map< std::string, std::set< std::string > > const & name3_generated_by_base_residue_name(
				cache_object()->name3_generated_by_base_residue_name() );

			if ( name3_generated_by_base_residue_name.find( base_residue_name ) ==
					name3_generated_by_base_residue_name.end() ) return false;
			std::set< std::string> const & name3_set = name3_generated_by_base_residue_name.find( base_residue_name )->second;
			return ( name3_set.count( name3 ) );
		}
	}

	//otherwise, obtain a write lock and let the cache update its internal maps
	utility::thread::WriteLockGuard write_lock( cache_object()->read_write_mutex() );
#endif
	std::map< std::string, std::set< std::string > > const & name3_generated_by_base_residue_name(
		cache_object()->name3_generated_by_base_residue_name() );

	if ( name3_generated_by_base_residue_name.find( base_residue_name ) ==
			name3_generated_by_base_residue_name.end() ) return false;
	std::set< std::string> const & name3_set = name3_generated_by_base_residue_name.find( base_residue_name )->second;
	return ( name3_set.count( name3 ) );
}

/// @brief Check if a base type (like "CYS") generates any types with a new
/// interchangeability group (like "SCY" (via cys_acetylated))
bool
ResidueTypeSet::generates_patched_residue_type_with_interchangeability_group( std::string const & base_residue_name,
	std::string const & interchangeability_group ) const
{
#ifdef MULTI_THREADED
	{ // scope the read lock
		utility::thread::ReadLockGuard read_lock( cache_object()->read_write_mutex() );
		if ( cache_object()->maps_up_to_date() ) {
			std::map< std::string, std::set< std::string > > const & interchangeability_group_generated_by_base_residue_name(
				cache_object()->interchangeability_group_generated_by_base_residue_name() );

			if ( interchangeability_group_generated_by_base_residue_name.find( base_residue_name ) ==
					interchangeability_group_generated_by_base_residue_name.end() ) {
				return false;
			}
			std::set< std::string> const & interchangeability_group_set =
				interchangeability_group_generated_by_base_residue_name.find( base_residue_name )->second;
			return interchangeability_group_set.count( interchangeability_group );
		}
	}

	// otherwise, obtain a write lock and let the cache update its internal maps
	utility::thread::WriteLockGuard write_lock( cache_object()->read_write_mutex() );
#endif
	std::map< std::string, std::set< std::string > > const & interchangeability_group_generated_by_base_residue_name(
		cache_object()->interchangeability_group_generated_by_base_residue_name() );

	if ( interchangeability_group_generated_by_base_residue_name.find( base_residue_name ) ==
			interchangeability_group_generated_by_base_residue_name.end() ) {
		return false;
	}
	std::set< std::string> const & interchangeability_group_set =
		interchangeability_group_generated_by_base_residue_name.find( base_residue_name )->second;
	return interchangeability_group_set.count( interchangeability_group );

}

void
ResidueTypeSet::prep_restype( ResidueTypeOP new_type ) {

	if ( option[ OptionKeys::in::file::assign_gasteiger_atom_types ] ) {
		gasteiger::GasteigerAtomTypeSetCOP gasteiger_set(
			ChemicalManager::get_instance()->gasteiger_atom_type_set() );
		gasteiger::assign_gasteiger_atom_types( *new_type, gasteiger_set, false );
	}

	if ( option[ OptionKeys::in::add_orbitals] ) {
		orbitals::AssignOrbitals add_orbitals_to_residue(new_type);
		add_orbitals_to_residue.assign_orbitals();
	}

}


void
ResidueTypeSet::add_base_residue_type( std::string const & filename )
{
	ResidueTypeOP rsd_type( read_topology_file( filename, atom_type_set(), element_set(), mm_atom_type_set(), orbital_type_set() ) );
	add_base_residue_type( rsd_type );
}

void
ResidueTypeSet::read_files_for_base_residue_types(
	utility::vector1< std::string > const & filenames
)
{
	for ( std::string const & filename : filenames ) {
		add_base_residue_type( filename );
	}
}

void
ResidueTypeSet::add_base_residue_type( ResidueTypeOP new_type )
{
	debug_assert( new_type );
	if ( mode() != new_type->mode() ) {
		TR.Warning << "ResidueType " << new_type->name() << " of mode " << new_type->mode()
			<< " is being added to a ResidueTypeSet of mode " << mode() << std::endl;
		// But we're doing it anyway (though we probably shouldn't).
	}

	prep_restype( new_type );
#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard write_lock( cache_object()->read_write_mutex() );
#endif
	cache_object()->add_residue_type( new_type );
	cache_object()->clear_cached_maps();
	base_residue_types_.push_back( new_type );
}

void
ResidueTypeSet::add_unpatchable_residue_type( std::string const & filename )
{
	ResidueTypeOP rsd_type( read_topology_file( filename, atom_type_set(), element_set(), mm_atom_type_set(), orbital_type_set() ) );
	add_unpatchable_residue_type( rsd_type );
}

void
ResidueTypeSet::read_files_for_unpatchable_residue_types(
	utility::vector1< std::string > const & filenames
)
{
	for ( std::string const & filename : filenames ) {
		add_unpatchable_residue_type( filename );
	}
}

void
ResidueTypeSet::add_unpatchable_residue_type( ResidueTypeOP new_type )
{
	debug_assert( new_type );
	if ( mode() != new_type->mode() ) {
		TR.Warning << "ResidueType " << new_type->name() << " of mode " << new_type->mode()
			<< " is being added to a ResidueTypeSet of mode " << mode() << std::endl;
		// But we're doing it anyway (though we probably shouldn't).
	}

	prep_restype( new_type );

#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard write_lock( cache_object()->read_write_mutex() );
#endif
	cache_object()->add_residue_type( new_type );
	cache_object()->clear_cached_maps();
	unpatchable_residue_types_.push_back( new_type );
}

void
ResidueTypeSet::remove_base_residue_type( std::string const & name )
{
	if ( ! has_name( name ) ) {
		utility_exit_with_message( "Attempting to remove ResidueType " + name + " from a ResidueTypeSet which doesn't contain it.");
	}

#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard write_lock( cache_object()->read_write_mutex() );
#endif

	ResidueTypeCOP rsd_type( cache_object()->name_map( name ) );
	ResidueTypeCOPs::iterator res_it = std::find( base_residue_types_.begin(), base_residue_types_.end(), rsd_type );
	runtime_assert( res_it != base_residue_types_.end() );
	base_residue_types_.erase( res_it );

	// TODO: Do we need to also remove any patched residue types from this one?
	cache_object()->remove_residue_type( name );
}

void
ResidueTypeSet::remove_unpatchable_residue_type( std::string const & name )
{
	if ( ! has_name( name ) ) {
		utility_exit_with_message( "Attempting to remove ResidueType " + name + " from a ResidueTypeSet which doesn't contain it.");
	}

#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard write_lock( cache_object()->read_write_mutex() );
#endif

	ResidueTypeCOP rsd_type( cache_object()->name_map( name ) );
	ResidueTypeCOPs::iterator res_it = std::find( unpatchable_residue_types_.begin(), unpatchable_residue_types_.end(), rsd_type );
	runtime_assert( res_it != unpatchable_residue_types_.end() );
	unpatchable_residue_types_.erase( res_it );

	cache_object()->remove_residue_type( name );
}

//////////////////////////////////////////////////////////////////////////////
/// @details helper function used during replacing residue types after, e.g., orbitals. Could possibly expand to update all maps.
bool
ResidueTypeSet::update_base_residue_types_if_replaced( ResidueTypeCOP rsd_type, ResidueTypeCOP rsd_type_new )
{
	if ( !base_residue_types_.has_value( rsd_type ) ) return false;
	base_residue_types_[ base_residue_types_.index( rsd_type ) ] = rsd_type_new;
	// TODO: Do we need to remove the old one from the cache?
	return true;
}

//////////////////////////////////////////////////////////////////////////////
/// @details Invoked only during residue-type construction of get_residue_type_write_locked,
/// and so it must invoke the version of has_name that presumes a write lock has been
/// obtained already (and will thus not invoke another function that tries to establish
/// either a read or a write lock on the ResidueTypeSetCache)
void
ResidueTypeSet::figure_out_last_patch_from_name( std::string const & rsd_name,
	std::string & rsd_name_base,
	std::string & patch_name ) const
{
	Size pos = rsd_name.find_last_of( PATCH_LINKER );
	if ( pos != std::string::npos ) {
		rsd_name_base = rsd_name.substr( 0, pos );
		patch_name    = rsd_name.substr( pos + 1 );
	} else { // Patch linker not found
		rsd_name_base = rsd_name;
		patch_name    = "";
	}

	// For D patch, it's the first letter.
	if ( patch_name.size() == 0 && rsd_name[ 0 ] == 'D' && has_name_write_locked( rsd_name.substr( 1 ) ) ) {
		rsd_name_base = rsd_name.substr( 1 );
		patch_name    = "D";
	}

	// For chiral-flip nucleic acid patch, it's the first letter.
	if ( patch_name.size() == 0 && rsd_name[ 0 ] == 'L' && has_name_write_locked( rsd_name.substr( 1 ) ) ) {
		rsd_name_base = rsd_name.substr( 1 );
		patch_name    = "L";
	}
}

//////////////////////////////////////////////////////////////////////////////
// Helpful ResidueTypeFinder back-hooks
//////////////////////////////////////////////////////////////////////////////

/// @brief Get the base ResidueType with the given aa type and variants
/// @details Returns 0 if one does not exist.
ResidueTypeCOP
ResidueTypeSet::get_representative_type_aa( AA aa, utility::vector1< std::string > const & variants ) const {
	return ResidueTypeFinder( *this ).aa( aa ).variants( variants ).get_representative_type();
}

ResidueTypeCOP
ResidueTypeSet::get_representative_type_aa( AA aa ) const {
	utility::vector1< std::string > const variants; // blank
	return get_representative_type_aa( aa, variants );
}

/// @brief Get the base ResidueType with the given name1 and variants
/// @details Returns 0 if one does not exist.
ResidueTypeCOP
ResidueTypeSet::get_representative_type_name1( char name1, utility::vector1< std::string > const & variants  ) const {
	return ResidueTypeFinder( *this ).name1( name1 ).variants( variants ).get_representative_type();
}

ResidueTypeCOP
ResidueTypeSet::get_representative_type_name1( char name1  ) const {
	utility::vector1< std::string > const variants; // blank
	return get_representative_type_name1( name1, variants );
}

/// @brief Get the base ResidueType with the given name3 and variants
/// @details Returns 0 if one does not exist.
ResidueTypeCOP
ResidueTypeSet::get_representative_type_name3( std::string const &  name3, utility::vector1< std::string > const & variants  ) const {
	return ResidueTypeFinder( *this ).name3( name3 ).variants( variants ).get_representative_type();
}

ResidueTypeCOP
ResidueTypeSet::get_representative_type_name3( std::string const &  name3  ) const {
	utility::vector1< std::string > const variants; // blank
	return get_representative_type_name3( name3, variants );
}

ResidueTypeCOP
ResidueTypeSet::get_representative_type_base_name( std::string const & base_name ) const {
	return ResidueTypeFinder( *this ).residue_base_name( base_name ).get_representative_type();
}

/// @brief Gets all non-patched types with the given aa type
ResidueTypeCOPs
ResidueTypeSet::get_base_types_aa( AA aa ) const {
	return ResidueTypeFinder( *this ).aa( aa ).get_possible_base_residue_types();
}

/// @brief Get all non-patched ResidueTypes with the given name1
ResidueTypeCOPs
ResidueTypeSet::get_base_types_name1( char name1 ) const {
	return ResidueTypeFinder( *this ).name1( name1 ).get_possible_base_residue_types();
}

/// @brief Get all non-patched ResidueTypes with the given name3
ResidueTypeCOPs
ResidueTypeSet::get_base_types_name3( std::string const & name3 ) const {
	return ResidueTypeFinder( *this ).name3( name3 ).get_possible_base_residue_types();
}

/// @brief Given a D-residue, get its L-equivalent.
/// @details Returns NULL if there is no equivalent, true otherwise.  Throws an error if this is not a D-residue.
/// Preserves variant types.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
ResidueTypeCOP
ResidueTypeSet::get_d_equivalent(
	ResidueTypeCOP l_rsd
) const {
	runtime_assert_string_msg( l_rsd, "Error in core::chemical::ResidueTypeSet::get_d_equivalent(): A null pointer was passed to this function!" );
	runtime_assert_string_msg( l_rsd->is_l_aa(), "Error in core::chemical::ResidueTypeSet::get_d_equivalent(): The residue passed to this function is not an L_AA!" );

	ResidueTypeCOP l_basetype( l_rsd->get_base_type_cop() );
	if ( !l_to_d_mapping_.count(l_basetype) ) return ResidueTypeCOP(); //Returns NULL pointer if there's no D-equivalent.

	ResidueTypeCOP d_basetype( l_to_d_mapping_.at(l_basetype) );

	//Add back the variants, as efficiently as possible:
	utility::vector1<VariantType> const variant_type_list( l_rsd->variant_type_enums() );
	utility::vector1<std::string> const & custom_variant_type_list( l_rsd->custom_variant_types() );

	ResidueTypeCOP d_rsd ( ResidueTypeFinder( *this ).base_type(d_basetype).variants( variant_type_list, custom_variant_type_list ).get_representative_type() );

	return d_rsd;
}

/// @brief Given an L-residue, get its D-equivalent.
/// @details Returns NULL if there is no equivalent, true otherwise.  Throws an error if this is not an L-residue.
/// Preserves variant types.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
ResidueTypeCOP
ResidueTypeSet::get_l_equivalent(
	ResidueTypeCOP d_rsd
) const {
	runtime_assert_string_msg( d_rsd, "Error in core::chemical::ResidueTypeSet::get_l_equivalent(): A null pointer was passed to this function!" );
	runtime_assert_string_msg( d_rsd->is_d_aa(), "Error in core::chemical::ResidueTypeSet::get_l_equivalent(): The residue passed to this function is not a D_AA!" );

	ResidueTypeCOP d_basetype( d_rsd->get_base_type_cop() );
	if ( !d_to_l_mapping_.count(d_basetype) ) return ResidueTypeCOP(); //Returns NULL pointer if there's no D-equivalent.

	ResidueTypeCOP l_basetype( d_to_l_mapping_.at(d_basetype) );

	//Add back the variants, as efficiently as possible:
	utility::vector1<VariantType> const variant_type_list( d_rsd->variant_type_enums() );
	utility::vector1<std::string> const & custom_variant_type_list( d_rsd->custom_variant_types() );

	ResidueTypeCOP l_rsd ( ResidueTypeFinder( *this ).base_type(l_basetype).variants( variant_type_list, custom_variant_type_list ).get_representative_type() );

	return l_rsd;
}

/// @brief Given a residue, get its mirror-image type.
/// @details Returns the same residue if this is an ACHIRAL type (e.g. gly), the D-equivalent for an L-residue, the L-equivalent of a D-residue,
/// or NULL if this is an L-residue with no D-equivalent (or a D- with no L-equivalent).  Preserves variant types.
ResidueTypeCOP
ResidueTypeSet::get_mirrored_type(
	ResidueTypeCOP original_rsd
) const {
	if ( original_rsd->is_achiral_backbone() ) return original_rsd;

	runtime_assert_string_msg( original_rsd->is_d_aa() || original_rsd->is_l_aa(), "Error in core::chemical::ResidueTypeSet::get_mirror_type(): The residue type must be achiral, or must have the D_AA or L_AA property." );

	if ( original_rsd->is_d_aa() ) return get_l_equivalent(original_rsd);
	if ( original_rsd->is_l_aa() ) return get_d_equivalent(original_rsd);

	return ResidueTypeCOP();
}

/// @brief Gets all types with the given aa type and variants
/// @details The number of variants must match exactly.
/// (It's assumed that the passed VariantTypeList contains no duplicates.)
ResidueTypeCOPs
ResidueTypeSet::get_all_types_with_variants_aa( AA aa, utility::vector1< std::string > const & variants ) const
{
	utility::vector1< VariantType > exceptions;
	return get_all_types_with_variants_aa( aa, variants, exceptions );
}

ResidueTypeCOPs
ResidueTypeSet::get_all_types_with_variants_aa( AA aa,
	utility::vector1< std::string > const & variants,
	utility::vector1< VariantType > const & exceptions ) const
{
	{ // scope the read lock
#ifdef MULTI_THREADED
		utility::thread::ReadLockGuard read_lock( cache_object()->read_write_mutex() );
#endif
		if ( cache_object()->all_types_with_variants_aa_already_cached( aa, variants, exceptions ) ) {
			return cache_object()->retrieve_all_types_with_variants_aa( aa, variants, exceptions );
		}
	}

	// Without creating a write lock, go collect the appropriate set of RTs for the
	// given query. Wait until this statement completes before creating the write lock,
	// or the locks that the RTF creates will create a deadlock
	ResidueTypeCOPs rts = ResidueTypeFinder( *this ).aa( aa ).variants( variants ).variant_exceptions( exceptions ).get_all_possible_residue_types();

#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard write_lock( cache_object()->read_write_mutex() );
#endif
	cache_object()->cache_all_types_with_variants_aa( aa, variants, exceptions, rts );
	return rts;
}

/// @brief Gets all types with the given name1 and variants
/// @brief Get all non-patched ResidueTypes with the given name1
/// @details The number of variants must match exactly.
/// (It's assumed that the passed VariantTypeList contains no duplicates.)
ResidueTypeCOPs
ResidueTypeSet::get_all_types_with_variants_name1( char name1, utility::vector1< std::string > const & variants ) const {
	return ResidueTypeFinder( *this ).name1( name1 ).variants( variants ).get_all_possible_residue_types();
}

/// @brief Gets all types with the given name3 and variants
/// @details The number of variants must match exactly.
/// (It's assumed that the passed VariantTypeList contains no duplicates.)
ResidueTypeCOPs
ResidueTypeSet::get_all_types_with_variants_name3( std::string const &  name3, utility::vector1< std::string > const & variants ) const {
	return ResidueTypeFinder( *this ).name3( name3 ).variants( variants ).get_all_possible_residue_types();
}

/// @brief any residue types with this name3?
bool
ResidueTypeSet::has_name3( std::string const & name3 ) const
{
	return ( get_representative_type_name3( name3 ) != nullptr );
}

bool
ResidueTypeSet::has_interchangeability_group( std::string const & interchangeability_group ) const
{
	ResidueTypeCOP rsd_type = ResidueTypeFinder( *this ).interchangeability_group( interchangeability_group ).get_representative_type();
	return ( rsd_type != nullptr );
}

///////////////////////////////////////////////////////////////////////////////
/// @details Return the first match with both base ResidueType id and variant_type name.  Abort if there is no match.
/// @note    Currently, this will not work for variant types defined as alternate base residues (i.e., different .params
///          files).
/// @remark  TODO: This should be refactored to make better use of the new ResidueProperties system. ~Labonte
ResidueType const &
ResidueTypeSet::get_residue_type_with_variant_added(
	ResidueType const & init_rsd,
	VariantType const new_type ) const
{
	if ( init_rsd.has_variant_type( new_type ) ) return init_rsd;

	// Find all residues with the same base name as init_rsd.
	std::string const base_name( residue_type_base_name( init_rsd ) );

	// the desired set of variant types:
	utility::vector1< std::string > target_variants( init_rsd.properties().get_list_of_variants() );
	if ( !init_rsd.has_variant_type(new_type) ) {
		target_variants.push_back( ResidueProperties::get_string_from_variant( new_type ) );
	}

	ResidueTypeCOP rsd_type = ResidueTypeFinder( *this ).residue_base_name( base_name ).variants( target_variants ).get_representative_type();

	if ( rsd_type == nullptr ) {
		utility_exit_with_message( "unable to find desired variant residue: " + init_rsd.name() + " " + base_name + " " +
			ResidueProperties::get_string_from_variant( new_type ) );
	}

	return *rsd_type;
}

MergeBehaviorManager const &
ResidueTypeSet::merge_behavior_manager() const
{
	return *merge_behavior_manager_;
}

///////////////////////////////////////////////////////////////////////////////
/// @note    Currently, this will not work for variant types defined as alternate base residues (i.e., different .params
///          files).
/// @remark  TODO: This should be refactored to make better use of the new ResidueProperties system. ~Labonte
ResidueType const &
ResidueTypeSet::get_residue_type_with_variant_removed(
	ResidueType const & init_rsd,
	VariantType const old_type) const
{
	if ( !init_rsd.has_variant_type( old_type ) ) return init_rsd;

	// find all residues with the same base name as init_rsd
	std::string const base_name( residue_type_base_name( init_rsd ) );

	// the desired set of variant types:
	utility::vector1< std::string > target_variants( init_rsd.properties().get_list_of_variants() );
	target_variants.erase( std::find( target_variants.begin(), target_variants.end(),
		ResidueProperties::get_string_from_variant( old_type ) ) );

	ResidueTypeCOP rsd_type = ResidueTypeFinder( *this ).residue_base_name( base_name ).variants( target_variants ).get_representative_type();

	if ( rsd_type == nullptr ) {
		utility_exit_with_message( "unable to find desired non-variant residue: " + init_rsd.name() + " " + base_name +
			" " + ResidueProperties::get_string_from_variant( old_type ) );
	}

	return *rsd_type;
}

void
ResidueTypeSet::set_merge_behavior_manager( MergeBehaviorManagerCOP mbm) {
	merge_behavior_manager_ = mbm;
}

// @details Copy constructor -- Not exposed by the Base class itself, or by the GlobalRTS, but it is used by the PoseRTS.
// The semantics are that we want a semi-shallow copy: all const data is shared (including ResidueTypes), but non-const
// data gets copied. This means we need to clone the ResidueTypeSetCache, so that two RTS don't share the same ResidueTypeSetCache
ResidueTypeSet::ResidueTypeSet( ResidueTypeSet const & src ) :
	utility::pointer::enable_shared_from_this< ResidueTypeSet >( src ),
	atom_types_( src.atom_types_), // const, so can share
	elements_( src.elements_ ), // const, so can share
	mm_atom_types_( src.mm_atom_types_ ), // const, so can share
	orbital_types_( src.orbital_types_ ), // const, so can share
	mode_( src.mode_ ),
	merge_behavior_manager_( src.merge_behavior_manager_ ), // const, so can share
	//cache_( src.cache_->clone(*this) ), // DEEP-ish COPY NEEDED!
	base_residue_types_( src.base_residue_types_ ), // individual entries are const, so can share
	unpatchable_residue_types_( src.unpatchable_residue_types_ ), // individual entries are const, so can share
	patches_( src.patches_ ), // individual entries are const, so can share
	metapatches_( src.metapatches_ ), // individual entries are const, so can share
	patch_map_( src.patch_map_ ), // individual entries are const, so can share
	metapatch_map_( src.metapatch_map_ ) // individual entries are const, so can share
	//l_to_d_mapping_;
	//d_to_l_mapping_;
{
#ifdef MULTI_THREADED
	utility::thread::ReadLockGuard read_lock( src.cache_->read_write_mutex() );
#endif
	cache_ = src.cache_->clone(*this); // DEEP-ish COPY NEEDED!
}

} // chemical
} // core

#ifdef    SERIALIZATION
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
bool
core::chemical::ResidueTypeSet::has( core::chemical::ResidueTypeCOP restype ) const {
#ifdef MULTI_THREADED
	utility::thread::ReadLockGuard read_lock( cache_->read_write_mutex() );
#endif
	return cache_->has_generated_residue_type( restype );
}

template< class Archive >
void
core::chemical::ResidueTypeSet::save( Archive & arc ) const {

	arc( CEREAL_NVP( atom_types_ ) );
	arc( CEREAL_NVP( elements_ ) );
	arc( CEREAL_NVP( mm_atom_types_ ) );
	arc( CEREAL_NVP( orbital_types_ ) );

	arc( CEREAL_NVP( mode_ ) ); // TypeSetMode

	arc( CEREAL_NVP( merge_behavior_manager_ ) ); // MergeBehaviorManagerCOP
	arc( CEREAL_NVP( cache_ ) ); // ResidueTypeSetCacheOP

	// Cereal should be smart enough to know that the COPs are duplicated from the Cache object,
	// and will only store them once.
	arc( CEREAL_NVP( base_residue_types_ ) ); // ResidueTypeCOPs
	arc( CEREAL_NVP( unpatchable_residue_types_ ) ); // ResidueTypeCOPs

	arc( CEREAL_NVP( patches_ ) ); // utility::vector1< PatchCOP >
	arc( CEREAL_NVP( metapatches_ ) ); // utility::vector1< MetapatchCOP >
	arc( CEREAL_NVP( patch_map_ ) ); // std::map< std::string, utility::vector1< PatchCOP > >
	arc( CEREAL_NVP( metapatch_map_ ) ); // std::map< std::string, MetapatchCOP >

	arc( CEREAL_NVP( l_to_d_mapping_ ) ); // std::map < ResidueTypeCOP, ResidueTypeCOP >
	arc( CEREAL_NVP( d_to_l_mapping_ ) ); // std::map < ResidueTypeCOP, ResidueTypeCOP >
}


/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::ResidueTypeSet::load( Archive & arc ) {
	using namespace core::chemical;

	arc( atom_types_ );
	arc( elements_ );
	arc( mm_atom_types_ );
	arc( orbital_types_ );

	arc( mode_ ); // TypeSetMode

	MergeBehaviorManagerOP merge_behavior_manager;
	arc( merge_behavior_manager );
	merge_behavior_manager_ = merge_behavior_manager;

	arc( cache_ ); // ResidueTypeSetCacheOP

	arc( base_residue_types_ );
	arc( unpatchable_residue_types_ );

	utility::vector1< PatchOP > patches;
	arc( patches );
	patches_ = patches;

	utility::vector1< MetapatchOP > metapatches;
	arc( metapatches );
	metapatches_ = metapatches;

	std::map< std::string, utility::vector1< PatchOP > > patch_map;
	arc( patch_map );
	patch_map_.clear(); patch_map_.insert( patch_map.begin(), patch_map.end() );

	std::map< std::string, MetapatchOP > metapatch_map;
	arc( metapatch_map );
	metapatch_map_.clear(); metapatch_map_.insert( metapatch_map.begin(), metapatch_map.end() );

	arc( l_to_d_mapping_ );
	arc( d_to_l_mapping_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::ResidueTypeSet );
CEREAL_REGISTER_TYPE( core::chemical::ResidueTypeSet )

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_ResidueTypeSet )
#endif // SERIALIZATION
