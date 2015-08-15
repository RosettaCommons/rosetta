// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////////////////////////////
/// @file ResidueTypeSet.cc
///
/// @brief
/// ResidueTypeSet class
///
/// @details
/// This class is responsible for iterating through the sets of residue types, including, but not limited to, amino
/// acids, nucleic acids, peptoid residues, and monosaccharides.  It first reads through a file that contains the
/// location of residue types in the database.  At the beginning of that file are the atom types, mm atom types,
/// element sets, and orbital types that will be used.  The sets are all for fa_standard.  If a new type of atom are
/// added for residues, this is where they would be added.  Once it assigns the types, it then reads in extra residue
/// params that are passed through the command line.  Finally, patches are applied to all residues added.
///
/// @author
/// Phil Bradley
/// Steven Combs - these comments
/// Rocco Moretti - helper functions for more efficient ResidueTypeSet use
/// Rhiju Das - 'on-the-fly' residue type instantiation/caching.
/////////////////////////////////////////////////////////////////////////


// Rosetta headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/adduct_util.hh>
#include <core/chemical/util.hh>
#include <core/chemical/Orbital.hh> /* for copying ResidueType */
#include <core/chemical/ResidueConnection.hh> /* for copying ResidueType */

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
#include <sstream>
#include <set>
#include <algorithm>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/chemical/orbitals/AssignOrbitals.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>

using namespace basic::options;

namespace core {
namespace chemical {

static thread_local basic::Tracer tr( "core.chemical.ResidueTypeSet" );

///////////////////////////////////////////////////////////////////////////////
/// @brief c-tor from directory
ResidueTypeSet::ResidueTypeSet(
		std::string const & name,
		std::string const & directory
) :
	name_( name ),
	database_directory_(directory)
{}

void ResidueTypeSet::init(
		std::vector< std::string > const & extra_res_param_files, // defaults to empty
		std::vector< std::string > const & extra_patch_files // defaults to empty
)
{
	using namespace basic::options;

	clock_t const time_start( clock() );

	on_the_fly_ = option[ OptionKeys::chemical::on_the_fly ]();
	if ( on_the_fly_ && ( option[ OptionKeys::packing::adducts ].user() ) ) {
		tr << tr.Red << "WARNING! Adducts are not compatible with on_the_fly residue_type_set. Setting on_the_fly to be false. You may have a long load-up time for this ResidueTypeSet. See notes on fixing this in place_adducts() in ResidueTypeSet.cc." << tr.Reset << std::endl;
		on_the_fly_ = false;
	}

	//XRW_B_T1
	//coarse::RuleSetOP coarsify_rule_set;
	//XRW_E_T1
	//ResidueTypeSetCAP fine_res_set;

	// read ResidueTypes
	{
		std::string const list_filename( database_directory_ + "residue_types.txt" );
		utility::io::izstream data( list_filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + list_filename + '\n' );
		}
		std::string line, tag;
		while ( getline( data, line ) ) {
			// Skip empty lines and comments.
			if ( line.size() < 1 || line[0] == '#' ) continue;

			// kp don't consider files for protonation versions of the residues if flag pH_mode is not used
			// to make sure even applications that use ResidueTypeSet directly never run into problems
			bool no_proton_states = false;
			if ( line.size() > 20 ){
				if ( ( !option[ OptionKeys::pH::pH_mode ].user() ) &&
						( line.substr (14,6) == "proton" ) ) {
					no_proton_states = true;
				}
			}
			if ( no_proton_states ) continue;

			// Skip carbohydrate ResidueTypes unless included with include_sugars flag.
			if ( ( ! option[ OptionKeys::in::include_sugars ] ) &&
					( line.substr( 0, 27 ) == "residue_types/carbohydrates" ) ) {
				continue;
			}

			// Skip lipid ResidueTypes unless included with include_lipids flag.
			if ( ( ! option[ OptionKeys::in::include_lipids ] ) &&
					( line.substr( 0, 20 ) == "residue_types/lipids" ) ) {
				continue;
			}

			//Skip mineral surface ResidueTypes unless included with surface_mode flag.
			if ((!option[OptionKeys::in::include_surfaces]) &&
				 (line.substr(0, 29) == "residue_types/mineral_surface"))
			{
				continue;
			}

			// Parse lines.
			std::istringstream l( line );
			l >> tag;
			if ( tag == "ATOM_TYPE_SET" ) {
				l >> tag;
				atom_types_ = ChemicalManager::get_instance()->atom_type_set( tag );
			} else if ( tag == "ELEMENT_SET" ) {
				l >> tag;
				elements_ = ChemicalManager::get_instance()->element_set( tag );
			} else if ( tag == "MM_ATOM_TYPE_SET" ) {
				l >> tag;
				mm_atom_types_ = ChemicalManager::get_instance()->mm_atom_type_set( tag );
			} else if(tag == "ORBITAL_TYPE_SET"){
				l >> tag;
				orbital_types_ = ChemicalManager::get_instance()->orbital_type_set(tag);
			} else {
				std::string const filename( database_directory_ + line );

				ResidueTypeOP rsd_type( read_topology_file(
						filename, atom_types_, elements_, mm_atom_types_, orbital_types_, get_self_weak_ptr() ) );

				residue_types_.push_back( rsd_type );
			}
		}

		BOOST_FOREACH(std::string filename, extra_res_param_files){
			ResidueTypeOP rsd_type( read_topology_file(
					filename, atom_types_, elements_, mm_atom_types_, orbital_types_, get_self_weak_ptr() ) );
			residue_types_.push_back( rsd_type );
		}

		update_residue_maps();

		base_residue_types_ = residue_types_;

	}  // ResidueTypes read

	// now apply patches
	{
		std::string const list_filename( database_directory_+"/patches.txt" );
		utility::io::izstream data( list_filename.c_str() );

		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open patch list file: "+list_filename );
		}

		// Read the command line and avoid applying patches that the user has requested be
		// skipped.  The flag allows the user to specify a patch by its name or by its file.
		// E.g. "SpecialRotamer" or "SpecialRotamer.txt".  Directory names will be ignored if given.
		std::set< std::string > patches_to_avoid;
		if ( option[ OptionKeys::chemical::exclude_patches ].user() ) {
			utility::vector1< std::string > avoidlist =
					option[ OptionKeys::chemical::exclude_patches ];
			for ( Size ii = 1; ii <= avoidlist.size(); ++ii ) {
				utility::file::FileName fname( avoidlist[ ii ] );
				patches_to_avoid.insert( fname.base() );
			}
		}

		utility::vector1< std::string > patch_filenames(extra_patch_files);
		// Unconditional loading of listed patches is deliberate --
		// if you specified it explicitly, you probably want it to load.
		std::string line;
		while ( getline( data,line) ) {
			if ( line.size() < 1 || line[0] == '#' ) continue;

			// get rid of any comment lines.
			line = utility::string_split( line, '#' )[1];
			line = utility::string_split( line, ' ' )[1];

			// Skip branching/conjugation patches unless included with read_pdb_link_records flag.
			if ( ( ! option[ OptionKeys::in::file::read_pdb_link_records ] ) &&
					( line.substr( 0, 17 ) == "patches/branching" ) ) {
				continue;
			}

			// Skip carbohydrate patches unless included with include_sugars flag.
			if ( ( ! option[ OptionKeys::in::include_sugars ] ) &&
					( line.substr( 0, 21 ) == "patches/carbohydrates" ) ) {
				continue;
			}

			// Skip this patch if the "patches_to_avoid" set contains the named patch.
			utility::file::FileName fname( line );
			if ( patches_to_avoid.find( fname.base() ) != patches_to_avoid.end() ) {
				tr << "While generating ResidueTypeSet " << name_ <<
						": Skipping patch " << fname.base() << " as requested" << std::endl;
				continue;
			}

			patch_filenames.push_back( database_directory_ + line );
		}

		// kdrew: include list allows patches to be included while being commented out in patches.txt,
		// useful for testing non-canonical patches.
		//tr << "include_patches activated? " <<
		//		option[ OptionKeys::chemical::include_patches ].active() << std::endl;
		if ( option[ OptionKeys::chemical::include_patches ].active() ) {
			utility::vector1< std::string > includelist =
					option[ OptionKeys::chemical::include_patches ];
			for ( Size ii = 1; ii <= includelist.size(); ++ii ) {
				utility::file::FileName fname( includelist[ ii ] );
				if ( !utility::file::file_exists( database_directory_ + includelist[ ii ] ) ) {
					tr.Warning << "Could not find: " << database_directory_+includelist[ii]  << std::endl;
					continue;
				}
				patch_filenames.push_back( database_directory_ + includelist[ ii ]);
				tr << "While generating ResidueTypeSet " << name_ <<
						": Including patch " << fname << " as requested" << std::endl;
			}
		}

		//fpd  if missing density is to be read correctly, we will have to also load the terminal truncation variants
		if ( option[ OptionKeys::in::missing_density_to_jump ]()
				|| option[ OptionKeys::in::use_truncated_termini ]() ) {
			if ( std::find( patch_filenames.begin(), patch_filenames.end(), database_directory_ + "patches/NtermTruncation.txt" )
					== patch_filenames.end())
				patch_filenames.push_back( database_directory_ + "patches/NtermTruncation.txt" );
			if ( std::find( patch_filenames.begin(), patch_filenames.end(), database_directory_ + "patches/CtermTruncation.txt" )
					== patch_filenames.end())
				patch_filenames.push_back( database_directory_ + "patches/CtermTruncation.txt" );

				// we'll assume user wants this
				option[ OptionKeys::chemical::override_rsd_type_limit ].value(true);

		}

		apply_patches( patch_filenames );
	}  // patches applied

	// Generate combinations of adducts as specified by the user
	place_adducts();

	if(option[ OptionKeys::in::add_orbitals]){
		for( Size ii = 1 ; ii <= residue_types_.size() ; ++ii ) {
			if ( residue_types_[ii]->finalized() ) {
				ResidueTypeOP rsd_type_clone = residue_types_[ii]->clone();
				orbitals::AssignOrbitals( rsd_type_clone ).assign_orbitals();
				update_base_residue_types_if_replaced( residue_types_[ ii ], rsd_type_clone );
				residue_types_[ ii ] = rsd_type_clone;
			}
		}
		update_residue_maps();
	}

	tr << "Finished initializing " << name_ << " residue type set.  ";
	tr << "Created " << residue_types_.size() << " residue types" << std::endl;
	tr << "Total time to initialize " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;

	// Sanity check. Might be good to tighten this limit -- fa_standard appears under 1500 residue types in most use cases.
	Size const MAX_RESIDUE_TYPES = 3000;
	if ( residue_types_.size() > MAX_RESIDUE_TYPES && !option[ OptionKeys::chemical::override_rsd_type_limit ]() ){
		if ( !option[ OptionKeys::mistakes::restore_pre_talaris_2013_behavior ]() )  {
			for ( Size n = 1; n <= residue_types_.size(); n++ )	 tr << residue_types_[n]->name() << std::endl;
			tr << "Created " << residue_types_.size() << " residue types" << std::endl;
		}
		std::string const error_message = "Number of residue types is greater than MAX_RESIDUE_TYPES. Rerun with -override_rsd_type_limit. Or if you have introduced a bunch of patches, consider declaring only the ones you want to use at the top of your app (with the options) with the command option[ chemical::include_patches ].push_back( ... ). ";
		std::cerr << error_message << std::endl;
		if ( !option[ OptionKeys::mistakes::restore_pre_talaris_2013_behavior ]() ) utility_exit_with_message( error_message );
	}
}


ResidueTypeSet::ResidueTypeSet() {}


ResidueTypeSet::~ResidueTypeSet() {}

///////////////////////////////////////////////////////////////////////////////
/// @details the file contains a list of names of residue type parameter files
/// stored in the database path
void
ResidueTypeSet::read_list_of_residues(
	std::string const & list_filename
)
{
	// read the files
	utility::vector1< std::string > filenames;
	{
		utility::io::izstream data( list_filename.c_str() );
		std::string line;
		while ( getline( data, line ) ) {
			// add full database path to the AA.params filename
			filenames.push_back( basic::database::full_name( line ) );
		}
		data.close();
	}

	read_files( filenames );
}

///////////////////////////////////////////////////////////////////////////////
void
ResidueTypeSet::read_files(
	utility::vector1< std::string > const & filenames
)
{
	for ( Size ii=1; ii<= filenames.size(); ++ii ) {
		ResidueTypeCOP rsd_type( read_topology_file( filenames[ii], atom_types_, elements_, mm_atom_types_,orbital_types_, get_self_weak_ptr() ) );
		residue_types_.push_back( rsd_type );
		base_residue_types_.push_back( rsd_type );
	}

	update_residue_maps();
}

///////////////////////////////////////////////////////////////////////////////
void
ResidueTypeSet::apply_patches(
	utility::vector1< std::string > const & filenames
)
{
	for ( Size ii=1; ii<= filenames.size(); ++ii ) {
		PatchOP p( new Patch );
		p->read_file( filenames[ii] );
		patches_.push_back( p );
		patch_map_[ p->name() ].push_back( p );
	}

	// These "replace_residue" patches are a special case, and barely in use anymore.
	// In their current implementation, they actually do *not* change the name of the ResidueType.
	// But we probably should keep track of their application
	//  by updating name() of residue; and hold copies of the replaced residues without patch applied in, e.g., replaced_name_map_.
	// That's going to require a careful refactoring of how residue types are accessed (e.g., can't just use name3_map() anymore),
	//  which I might do later. For example, there's code in SwitchResidueTypeSet that will look for "MET" when it really should
	//  look for "MET:protein_centroid_with_HA" and know about this patch.
	// For now, apply them first, and force application/instantiation later.
	//	-- rhiju
	for ( Size ii=1; ii<= patches_.size(); ++ii ) {
		PatchCOP p( patches_[ ii ] );
		Size current_n_residue_types( residue_types_.size() );
		for ( Size i=1; i<= current_n_residue_types; ++i ) {
			ResidueType const & rsd_type( *residue_types_[ i ] );
			if ( p->applies_to( rsd_type ) && p->replaces( rsd_type ) ) {
				runtime_assert( rsd_type.finalized() );
				ResidueTypeCOP rsd_type_new = p->apply( rsd_type );
				update_base_residue_types_if_replaced( residue_types_[ i ], rsd_type_new );
				residue_types_[ i ] = rsd_type_new;
			}
		}
	}

	// This is the main loop.
	// Note that this spawns a number of residue types that is exponentially large in the number of patches!
	// Use of placeholder residue types (-on_the_fly) significantly decreases the memory & time footprint of this operation,
	//  but still this probably should not happen -- rhiju.
	for ( Size ii=1; ii<= patches_.size(); ++ii ) {
		PatchCOP p( patches_[ ii ] );
		Size current_n_residue_types( residue_types_.size() );
		for ( Size i=1; i<= current_n_residue_types; ++i ) {
			ResidueType const & rsd_type( *residue_types_[ i ] );
			if ( p->applies_to( rsd_type ) ) {
				if ( p->replaces( rsd_type ) ) {
					continue; // Now replace_residue patches are applied in a 'first pass' above.
				}
				bool const instantiate = (!on_the_fly_) || rsd_type.is_carbohydrate() /*hack for now -- carbohydrate info is complicated.*/;
				if ( instantiate ) {
					runtime_assert( rsd_type.finalized() );
				}
				ResidueTypeCOP new_rsd_type( p->apply( rsd_type, instantiate ) );
				if ( new_rsd_type != 0 ) {
					residue_types_.push_back( new_rsd_type );
					if ( new_rsd_type->name3() != rsd_type.name3() ) {
						// name3 is often used to figure out base_residue_type on which to apply patches; sometimes it switches, e.g., CYS-->CYD,
						// and we need to know that.
						name3_generated_by_base_residue_name_[ residue_type_base_name( rsd_type ) ].insert( new_rsd_type->name3() );
					}
				} else {
					utility_exit_with_message( "Could not apply patch " + p->name() + " to " + rsd_type.name() );
				}
			}
		}
	}

	update_residue_maps();
}


ResidueTypeCOPs const &
ResidueTypeSet::interchangeability_group_map_DO_NOT_USE( std::string const & name ) const
{
	std::map< std::string, ResidueTypeCOPs >::const_iterator iter = interchangeability_group_map_.find( name );
	if ( iter == interchangeability_group_map_.end() ) {
		return empty_residue_list_;
	}
	make_sure_instantiated( iter->second );
	return iter->second;
}

bool
ResidueTypeSet::has_interchangeability_group( std::string const & name ) const
{
	std::map< std::string, ResidueTypeCOPs >::const_iterator iter = interchangeability_group_map_.find( name );
	if ( iter == interchangeability_group_map_.end() ) {
		return false;
	}
	return true;
}

/// @details 3-letter name is not unique to each ResidueType
/// for example, 3-letter name "HIS" matches both his tautomers,
/// HIS and HIS_D. Return an empty list if no match is found.
ResidueTypeCOPs const &
ResidueTypeSet::name3_map_DO_NOT_USE( std::string const & name ) const
{
	debug_assert( name.size() == 3 );
	std::map< std::string, ResidueTypeCOPs >::const_iterator iter = name3_map_.find( name );
	if ( iter == name3_map_.end() ) {
		return empty_residue_list_;
	}
	/*bool did_instantiation = */ make_sure_instantiated( iter->second );
	//	if ( did_instantiation ) tr << tr.Red << "Instantiating name3_map for: " << name << tr.Reset << std::endl;
	return iter->second;
}

/// @details since residue id is unique, it only returns
/// one residue type or exit without match.
ResidueType const &
ResidueTypeSet::name_map( std::string const & name_in ) const
{
	std::string const name = fixup_patches( name_in );
	ResidueTypeCOP rsd_type( 0 );
	if ( name_map_.find( name ) != name_map_.end() ) {
		rsd_type = name_map_.find( name )->second;
		//  Do not delete yet -- may revive later to 'properly' handle replace_residue_types -- rhiju.
		//	} else if ( replaced_name_map_.find( name ) != replaced_name_map_.end() ) {
		//		rsd_type = replaced_name_map_.find( name )->second;
	} else {
		for ( std::map< std::string, ResidueTypeCOP >::const_iterator it = name_map_.begin();
					it != name_map_.end(); ++it ) tr << it->first << std::endl;
		//  Do not delete yet -- may revive later to 'properly' handle replace_residue_types -- rhiju.
		//		for ( std::map< std::string, ResidueTypeCOP >::const_iterator it = replaced_name_map_.begin();
		//					it != replaced_name_map_.end(); ++it ) tr << it->first << std::endl;
		utility_exit_with_message( "unrecognized residue name '"+name+"'" );
	}
	make_sure_instantiated( rsd_type );
	return *( rsd_type );
}

/// @brief Check if a base type (like "SER") generates any types with another name3 (like "SEP")
bool
ResidueTypeSet::generates_patched_residue_type_with_name3( std::string const base_residue_name, std::string const name3 ) const
{
	if ( name3_generated_by_base_residue_name_.find( base_residue_name ) ==
			 name3_generated_by_base_residue_name_.end() ) return false;
	std::set< std::string> const & name3_set = name3_generated_by_base_residue_name_.find( base_residue_name )->second;
	return ( name3_set.count( name3 ) );
}

/// @details Instantiates ResidueType on-the-fly if it isn't finalized.
///  Watch out: Does not obey constness! Not thread-safe if on_the_fly_ is on!
bool
ResidueTypeSet::make_sure_instantiated( ResidueTypeCOP const & rsd_type ) const
{
	if ( rsd_type->finalized() ) return false;

	runtime_assert( on_the_fly_ )
	// get name (which holds patch information)
	std::string rsd_name_base, patch_name;
	figure_out_last_patch_from_name( rsd_type->name(), rsd_name_base, patch_name );
	if ( patch_name.size() == 0 ) utility_exit_with_message( "Problem: " + rsd_type->name() + " not instantiated." );

	// now apply patch.
	ResidueType const & rsd_base = name_map( rsd_name_base );
	runtime_assert( rsd_base.finalized() );

	runtime_assert( patch_map_.find( patch_name ) != patch_map_.end() );
	utility::vector1< PatchCOP > const & patches = patch_map_.find( patch_name )->second;
	bool patch_applied( false );

	// sometimes patch cases are split between several patches -- look through all:
	for ( Size n = 1; n <= patches.size(); n++ ) {
		PatchCOP p = patches[ n ];

		if ( !p->applies_to( rsd_base ) ) continue;
		runtime_assert( !patch_applied ); // patch cannot be applied twice.

		ResidueTypeOP rsd_instantiated = p->apply( rsd_base );

		if ( rsd_instantiated == 0 ) {
			utility_exit_with_message(  "Failed to apply: " + p->name() + " to " + rsd_base.name() );
		} else {

			if (option[ OptionKeys::in::add_orbitals]) orbitals::AssignOrbitals( rsd_instantiated ).assign_orbitals();

			replace_residue_type_in_set_defying_constness( rsd_type, *rsd_instantiated );
		}
		patch_applied = true;
	}

	if ( !patch_applied ) {
		utility_exit_with_message( "What? Patch " + patch_name + " does not apply to " + rsd_base.name() + "  but user is asking for " + rsd_type->name() );
	}

	if ( !patch_applied ) {
		utility_exit_with_message( "What? Patch " + patch_name + " does not apply to " + rsd_base.name() + "  but user is asking for " + rsd_type->name() );
	}

	return true;
}

/// @details Instantiates ResidueType on-the-fly if it isn't finalized.
bool
ResidueTypeSet::make_sure_instantiated( utility::vector1< ResidueTypeCOP > const & rsd_types ) const
{
	bool did_instantiation( false );
	for ( Size n = 1; n <= rsd_types.size(); n++ ) {
		if ( make_sure_instantiated( rsd_types[ n ] ) ) did_instantiation = true;
	}
	return did_instantiation;
}

//////////////////////////////////////////////////////////////////////////////
void
ResidueTypeSet::figure_out_last_patch_from_name( std::string const rsd_name,
																 std::string & rsd_name_base,
																 std::string & patch_name ) const
{
	Size pos = rsd_name.find_last_of( PATCH_LINKER );
	runtime_assert( pos != std::string::npos );
	rsd_name_base = rsd_name.substr( 0, pos );
	patch_name    = rsd_name.substr( pos + 1 );
}

/// @details helper function used during replacing residue types after, e.g., orbitals. Could possibly expand to update all maps.
void
ResidueTypeSet::update_base_residue_types_if_replaced( ResidueTypeCOP rsd_type, ResidueTypeCOP rsd_type_new )
{
	if ( !base_residue_types_.has_value( rsd_type ) ) return;
	base_residue_types_[ base_residue_types_.index( rsd_type ) ] = rsd_type_new;
}


/// @details OH MY GAWRSH! oh my gosh oh my gosh oh my gosh. Gasp!
///  If we get rid of aa_map_DO_NOT_USE, etc. we actually might be able to get rid of this call.
///  Could instead simply change out the ResidueTypeCOP held in residue_types_ and base_residue_types_ (see ,e.g., update_base_residue_types_if_replaced), and use a
///  mutex to ensure thread safety.
void
ResidueTypeSet::replace_residue_type_in_set_defying_constness( ResidueTypeCOP rsd_type, ResidueType const & rsd_new ) const
{
	const_cast< ResidueType & >( *rsd_type ) = rsd_new;
}


ResidueTypeCOP
ResidueTypeSet::get_residue_type_with_patch( PatchCOP patch, ResidueTypeCOP rsd_type ) const
{
	runtime_assert( patch->applies_to( *rsd_type ) );
	runtime_assert( patches_.has_value( patch ) );
	runtime_assert( residue_types_.has_value( rsd_type ) );
	runtime_assert( rsd_type->finalized() );
	std::string name = rsd_type->name() + PATCH_LINKER + patch->name();
	return name_map( name ).get_self_ptr();
}

/// @brief Get the base ResidueType with the given aa type
/// @details Returns 0 if one does not exist.
ResidueTypeCOP
ResidueTypeSet::get_representative_type_aa( AA aa ) const {
	return ResidueTypeFinder( *this ).aa( aa ).get_representative_type();
}

/// @brief Get the base ResidueType with the given name1
/// @details Returns 0 if one does not exist.
ResidueTypeCOP
ResidueTypeSet::get_representative_type_name1( char name1 ) const {
	return ResidueTypeFinder( *this ).name1( name1 ).get_representative_type();
}

/// @brief Get the base ResidueType with the given name3
/// @details Returns 0 if one does not exist.
ResidueTypeCOP
ResidueTypeSet::get_representative_type_name3( std::string const &  name3 ) const {
	return ResidueTypeFinder( *this ).name3( name3 ).get_representative_type();
}

/// @brief Get the base ResidueType with the given aa type and variants
/// @details Returns 0 if one does not exist.
ResidueTypeCOP
ResidueTypeSet::get_representative_type_aa( AA aa, utility::vector1< std::string > const & variants ) const {
	return ResidueTypeFinder( *this ).aa( aa ).variants( variants ).get_representative_type();
}

/// @brief Get the base ResidueType with the given name1 and variants
/// @details Returns 0 if one does not exist.
ResidueTypeCOP
ResidueTypeSet::get_representative_type_name1( char name1, utility::vector1< std::string > const & variants ) const {
	return ResidueTypeFinder( *this ).name1( name1 ).variants( variants ).get_representative_type();
}

/// @brief Get the base ResidueType with the given name3 and variants
/// @details Returns 0 if one does not exist.
ResidueTypeCOP
ResidueTypeSet::get_representative_type_name3( std::string const &  name3, utility::vector1< std::string > const & variants ) const {
	return ResidueTypeFinder( *this ).name3( name3 ).variants( variants ).get_representative_type();
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
ResidueTypeSet::get_base_types_name3( std::string const &  name3 ) const {
	return ResidueTypeFinder( *this ).name3( name3 ).get_possible_base_residue_types();
}

/// @brief Gets all types with the given aa type and variants
/// @details The number of variants must match exactly.
/// (It's assumed that the passed VariantTypeList contains no duplicates.)
ResidueTypeCOPs
ResidueTypeSet::get_all_types_with_variants_aa( AA aa, utility::vector1< std::string > const & variants ) const {
	std::pair< AA, utility::vector1< std::string > > query( std::make_pair( aa, variants ) );
	if ( cached_aa_variants_map_.find( query ) == cached_aa_variants_map_.end() ) {
		cached_aa_variants_map_[ query ] = ResidueTypeFinder( *this ).aa( aa ).variants( variants ).get_all_possible_residue_types();
	}
	return cached_aa_variants_map_[ query ];
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

bool
ResidueTypeSet::has_name( std::string const & name ) const
{
	return ( name_map_.find( name ) != name_map_.end() );
}

bool
ResidueTypeSet::has_name3( std::string const & name3 ) const
{
	return ( name3_map_.find( name3 ) != name3_map_.end() );
}

/// @details similar to name3_map, return all matched residue types
/// or an empty list.
ResidueTypeCOPs const &
ResidueTypeSet::aa_map_DO_NOT_USE( AA const & aa ) const
{
	std::map< AA, ResidueTypeCOPs >::const_iterator iter =	aa_map_.find( aa );
	if ( iter == aa_map_.end() ) {
		return empty_residue_list_;
	}
	/*bool did_instantiation =*/ make_sure_instantiated( iter->second );
	//	if ( did_instantiation ) tr << tr.Red << "Instantiating aa_map for: " << aa << tr.Reset << std::endl;
	return aa_map_.find( aa )->second;
}


///////////////////////////////////////////////////////////////////////////////
//private
/// @brief clear residue  maps
void
ResidueTypeSet::clear_residue_maps()
{
	aa_map_.clear();
	name3_map_.clear();
	name_map_.clear();
	aas_defined_.clear();
}


///////////////////////////////////////////////////////////////////////////////
//private
/// @details check that residue id map should be unique,
/// sort the aas_defined list and make its members unique
void
ResidueTypeSet::update_residue_maps()
{
	clear_residue_maps();

	for ( ResidueTypeCOPs::iterator rsdtype_it(residue_types_.begin() ), rsdtype_end( residue_types_.end() );
			rsdtype_it != rsdtype_end; ++rsdtype_it ) {
		ResidueTypeCOP rsd( *rsdtype_it );
		add_residue_type_to_maps( rsd );
	}
	aas_defined_.sort();
	aas_defined_.unique();
}

void
ResidueTypeSet::add_residue_type_to_maps( ResidueTypeCOP rsd_ptr )
{

	// name should be unique!
	if ( name_map_.count( rsd_ptr->name() ) ) {
		std::stringstream err_msg;
		err_msg << "Attempting to add a residue type with name '" + rsd_ptr->name() + "' ";
		err_msg << "but this name is already taken." << std::endl;
		err_msg << "Please check your residue parameter files." << std::endl;
		utility_exit_with_message(err_msg.str());
	}
	name_map_[ rsd_ptr->name() ] = rsd_ptr;

	// map by AA
	if ( rsd_ptr->aa() != aa_unk ) {
		aa_map_[ rsd_ptr->aa() ].push_back( rsd_ptr );
	}

	// add aa type
	aas_defined_.push_back( rsd_ptr->aa() );

	// map by pdb string
	//	tr << "ADDING name3 to map " << rsd_ptr->name3() << "  based on rsd_type " << rsd_ptr->name() << std::endl;
	name3_map_[ rsd_ptr->name3() ].push_back( rsd_ptr );

	// map by interchangeability group
	interchangeability_group_map_[ rsd_ptr->interchangeability_group() ].push_back( rsd_ptr );

	// For specialty amino acids, add them to the name three maps both with their PDB strings and
	// with their specialty string -- the first three letters of the residue name.
	// E.g., CYD will appear in both lists for name3_map_[ "CYS" ] and name3_map_[ "CYD" ]
	if ( rsd_ptr->name3() != rsd_ptr->name().substr(0,3) ) {
		name3_map_[ rsd_ptr->name().substr(0,3) ].push_back( rsd_ptr );
	}

}

void
ResidueTypeSet::remove_residue_type_from_maps( ResidueTypeCOP rsd )
{
	//Assert rather than utility exit, because this should have been called by remove residue, which utility exits
	debug_assert(has_name(rsd->name()));

	ResidueTypeCOPs::iterator const_res_it(std::find(residue_types_.begin(),residue_types_.end(),rsd));
	residue_types_.erase(const_res_it);

	name_map_.erase(rsd->name());

	// clear this residue type from the aa_map_
	if(rsd->aa() != aa_unk) {
		ResidueTypeCOPs::iterator aa_it(std::find(aa_map_[rsd->aa()].begin(),aa_map_[rsd->aa()].end(),rsd));
		aa_map_[rsd->aa()].erase(aa_it);
		if(aa_map_[rsd->aa()].size() == 0)
		{
			aa_map_.erase(rsd->aa());

			// only clear the aa type from the aas_defined_ list if it was the last aa of this type in the aa_map_
			// special case for aa_unk ahead
			std::list<AA>::iterator aa_it(std::find(aas_defined_.begin(),aas_defined_.end(),rsd->aa()));
			aas_defined_.erase(aa_it);
		}
	} else {
		// if this was the only residue type in the ResidueTypeSet with an aa() of aa_unk, then
		// remove aa_unk from the aas_defined_ list.
		bool only_aa_unk = true;
		for ( ResidueTypeCOPs::const_iterator iter = residue_types_.begin(); iter != residue_types_.end(); ++iter ) {
			if ( (*iter)->aa() == aa_unk ) {
				only_aa_unk = false;
				break;
			}
		}
		if ( only_aa_unk ) {
			aas_defined_.erase( std::find( aas_defined_.begin(), aas_defined_.end(), aa_unk ) );
		}
	}

	// clear this residue type from the name3_map_
	ResidueTypeCOPs::iterator name3_it(std::find(name3_map_[rsd->name3()].begin(),name3_map_[rsd->name3()].end(),rsd));
	name3_map_[rsd->name3()].erase(name3_it);
	if(name3_map_[rsd->name3()].size() == 0)
	{
		name3_map_.erase(rsd->name3());
	}

	// clear this residue type from the interchangeability_group_map_
	ResidueTypeCOPs::iterator intgrp_it( std::find(
		interchangeability_group_map_[ rsd->interchangeability_group() ].begin(),
		interchangeability_group_map_[ rsd->interchangeability_group() ].end(),
		rsd ));
	interchangeability_group_map_[ rsd->interchangeability_group() ].erase( intgrp_it );
	if ( interchangeability_group_map_[ rsd->interchangeability_group() ].size() ) {
		interchangeability_group_map_.erase( rsd->interchangeability_group() );
	}

	if(rsd->name3() != rsd->name().substr(0,3))
	{
		name3_it = std::find(name3_map_[rsd->name3().substr(0,3)].begin(),name3_map_[rsd->name3().substr(0,3)].end(),rsd);
		name3_map_[rsd->name3().substr(0,3)].erase(name3_it);
		if(name3_map_[rsd->name3().substr(0,3)].size() == 0)
		{
			name3_map_.erase(rsd->name3().substr(0,3));
		}
	}

}

/// @brief beginning of aas_defined_ list
std::list< AA >::const_iterator
ResidueTypeSet::aas_defined_begin() const
{
	return aas_defined_.begin();
}

/// @brief end of aas_defined_ list
std::list< AA >::const_iterator
ResidueTypeSet::aas_defined_end() const
{
	return aas_defined_.end();
}

///////////////////////////////////////////////////////////////////////////////
/// @details selection done by ResidueTypeSelector class
void
ResidueTypeSet::select_residues_DO_NOT_USE(
	ResidueTypeSelector const & selector,
	ResidueTypeCOPs & matches
) const
{
	for ( ResidueTypeCOPs::const_iterator iter=residue_types_.begin(), iter_end = residue_types_.end(); iter!= iter_end; ++iter ) {
		make_sure_instantiated( *iter );
		if ( selector[ **iter ] ) {
			matches.push_back( *iter );
		}
	}
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

	if ( rsd_type == 0 ) {
		utility_exit_with_message( "unable to find desired variant residue: " + init_rsd.name() + " " + base_name + " " +
															 ResidueProperties::get_string_from_variant( new_type ) );
	}

	return *rsd_type;
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

	if ( rsd_type == 0 ) {
		utility_exit_with_message( "unable to find desired non-variant residue: " + init_rsd.name() + " " + base_name +
															 " " + ResidueProperties::get_string_from_variant( old_type ) );
	}

	return *rsd_type;
}

///////////////////////////////////////////////////////////////////////////////
/// @details Generation of new residue types augmented by adduct atoms
/// @note    This is almost perfectly consistent with on_the_fly residue type sets.
///            Why not just make adducts a patch? Then it will all work together nicely.
///            Just need to have the -adduct string (add_map) turned into a PatchSelector
///                 -- rhiju
void
ResidueTypeSet::place_adducts()
{
	// First parse the command line for requested adducts
	utility::options::StringVectorOption & add_set
		= option[ OptionKeys::packing::adducts ];

	// No adducts, skip out
	if( add_set.size() == 0 ) return;

	// Convert to a map that takes a string descriptor of an adduct and
	// gives the max number of adducts of that class to apply, so
	// a command line option of -adducts <adduct_type> 2 for a type that has
	// 5 entries in the rsd param file will create all combinations with up to
	// 2 total adducts.

	AdductMap add_map = parse_adduct_string( add_set );

	// Error check each requested adduct from the command line, and
	// complain if there are no examples in any residues.  This function
	// will not return if
	error_check_requested_adducts( add_map, residue_types_ );

	// Set up a starting point map where the int value is the number
	// of adducts of a given type placed
		AdductMap blank_map( add_map );
		for( AdductMap::iterator add_iter = blank_map.begin(),
						end_iter = blank_map.end() ;
					add_iter != end_iter ; ++add_iter ) {
			add_iter->second = 0;
		}

	// Process the residues in turn
	for ( ResidueTypeCOPs::const_iterator iter= residue_types_.begin(), iter_end = residue_types_.end();
				iter != iter_end; ++iter ) {
		ResidueType const & rsd( **iter );
		if ( !rsd.finalized() ) continue;
		AdductMap count_map( blank_map );
		utility::vector1< bool > add_mask( rsd.defined_adducts().size(), false  );
		create_adduct_combinations( rsd, add_map, count_map, add_mask, rsd.defined_adducts().begin() );
	}

//	utility_exit_with_message( "Debug stop point \n" );

	update_residue_maps();

}

void
ResidueTypeSet::add_residue_type( ResidueTypeOP new_type )
{
	if(option[ OptionKeys::in::add_orbitals]){
		orbitals::AssignOrbitals add_orbitals_to_residue(new_type);
		add_orbitals_to_residue.assign_orbitals();
	}
	residue_types_.push_back( new_type );
 	if ( residue_type_base_name( *new_type ) == new_type->name() ) base_residue_types_.push_back( new_type );
	add_residue_type_to_maps( new_type );
	aas_defined_.sort();
	aas_defined_.unique();
}

void
ResidueTypeSet::add_residue_type(std::string const & filename)
{
	ResidueTypeOP rsd_type( read_topology_file( filename, atom_types_, elements_, mm_atom_types_, orbital_types_, get_self_weak_ptr() ) );
	add_residue_type(rsd_type);
}

void
ResidueTypeSet::remove_residue_type(std::string const & name)
{

	if(!has_name(name))
	{
		utility_exit_with_message("ResidueTypeSet does not have a residue called "+name+ " so it cannot be deleted.");
	}

	ResidueTypeCOP type_to_remove( name_map_[ name ] );

	remove_residue_type_from_maps(type_to_remove);
}

/// @brief Create correct combinations of adducts for a residue type
void
ResidueTypeSet:: create_adduct_combinations(
	ResidueType const & rsd,
	AdductMap ref_map,
	AdductMap count_map,
	utility::vector1< bool > add_mask,
	utility::vector1< Adduct >::const_iterator work_iter
)
{

	if( work_iter == rsd.defined_adducts().end() ) {
		// Skip the 'no adduct' case - that has already been
		// made when reading in files
		if( std::find( add_mask.begin(), add_mask.end(), true ) == add_mask.end() ) {
			return;
		}
		// Make this combo and return;
		//		std::cout << "Making an adduct" << std::endl;

		utility::vector1< Adduct >::const_iterator add_iter = rsd.defined_adducts().begin() ;
		BOOST_FOREACH(bool make, add_mask){
			std::cout << "Adduct " << add_iter->adduct_name() << " make is " << make << std::endl;
			++add_iter;
		}

		// Farm this out to a helper function
		residue_types_.push_back( apply_adducts_to_residue( rsd, add_mask ) );

		return;
	}

	// Traverse the 'make' branch for this adduct if:
	// 1. The adduct is in the map of requested adducts
	// 2. we haven't exceeded the count limit for this adduct
	AdductMap::iterator test_iter =
			ref_map.find( work_iter->adduct_name() );

	if ( test_iter != ref_map.end() &&
				count_map[ test_iter->first ] < ref_map[ test_iter->first ]   ) {
		AdductMap new_count_map( count_map );
		new_count_map[ work_iter->adduct_name() ]++;
		utility::vector1< bool > new_add_mask( add_mask );
		// This following line may not work if the Adducts are no longer
		// stored in a vector
		new_add_mask[ work_iter - rsd.defined_adducts().begin() + 1 ] = true;
		create_adduct_combinations( rsd, ref_map, new_count_map, new_add_mask, work_iter+1 );
	}

	// Always traverse the 'do not make' for this adduct
	// The count is not incremented, and the mask is left at the default (false)
	AdductMap new_count_map( count_map );
	utility::vector1< bool > new_add_mask( add_mask );
	create_adduct_combinations( rsd, ref_map, new_count_map, new_add_mask, work_iter+1 );

}

} // pose
} // core
