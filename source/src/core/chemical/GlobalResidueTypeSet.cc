// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////////////////////////////
/// @file GlobalResidueTypeSet.cc
///
/// @brief
/// GlobalResidueTypeSet class
///
/// @details
/// This class is responsible for iterating through the sets of residue types, including, but not limited to, amino
///  acids, nucleic acids, peptoid residues, and monosaccharides.  It first reads through a file that contains the
///  location of residue types in the database.  At the beginning of that file are the atom types, mm atom types,
///  element sets, and orbital types that will be used.  The sets are all for fa_standard.  If a new type of atom is
///  added for residues, this is where they would be added.  Once it assigns the types, it then reads in extra residue
///  params that are passed through the command line.
/// Then, the class reads in patches.
/// There can be an exponentially large number of possible residue types that can be built from base_residue_types and
///  patch applications -- they are created on the fly when requested by accessor functions like has_name() and
///  get_all_types_with_variants_aa().
/// Later, the class can accept 'unpatchable residue types' (to which patches will not be added).
///
/// @author
/// Phil Bradley
/// Steven Combs - these comments
/// Rocco Moretti - helper functions for more efficient GlobalResidueTypeSet use
/// Rhiju Das - 'on-the-fly' residue type generation; allowing unpatchable residue types; ResidueTypeSetCache.
/////////////////////////////////////////////////////////////////////////


// Rosetta headers
#include <core/chemical/GlobalResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSetCache.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/Metapatch.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MergeBehaviorManager.hh>
#include <core/chemical/sdf/MolFileIOReader.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/adduct_util.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueDatabaseIO.hh>
#include <core/chemical/Orbital.hh> /* for copying ResidueType */
#include <core/chemical/ResidueConnection.hh> /* for copying ResidueType */

#include <core/chemical/gasteiger/GasteigerAtomTyper.hh>
#include <core/chemical/mmCIF/mmCIFParser.hh>

// Basic headers
#include <basic/database/open.hh>
#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/sql_database/types.hh>

// C++ headers
#include <fstream>
#include <string>
#include <set>
#include <algorithm>

// option key includes
#include <basic/options/keys/pH.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <core/chemical/orbitals/AssignOrbitals.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/file/file_sys_util.hh>

using namespace basic::options;

namespace core {
namespace chemical {

static THREAD_LOCAL basic::Tracer TR( "core.chemical.GlobalResidueTypeSet" );

///////////////////////////////////////////////////////////////////////////////
/// @brief c-tor from directory
GlobalResidueTypeSet::GlobalResidueTypeSet(
	std::string const & cm_name,
	std::string const & directory
) :
	ResidueTypeSet(),
	database_directory_(directory)
{
	clock_t const time_start( clock() );

	name( cm_name );
	set_merge_behavior_manager( MergeBehaviorManagerCOP( new MergeBehaviorManager( directory ) ) );
	load_shadowed_ids( directory );

	init_restypes_from_database(); // Will also load the sub-typesets
	init_restypes_from_commandline();

	init_patches_from_commandline(); // Patch resolution is order dependent - allow commandline to overrule
	init_patches_from_database();
	deal_with_patch_special_cases();

	// Generate combinations of adducts as specified by the user
	place_adducts();

	// Components file? (Will be empty if we're not doing this.)
	if ( option[ OptionKeys::in::file::load_PDB_components ] || option[ OptionKeys::in::file::PDB_components_file ].user() ) {
		pdb_components_filename( option[ OptionKeys::in::file::PDB_components_file ].value() );
	}

	TR << "Finished initializing " << name() << " residue type set.  ";
	TR << "Created " << base_residue_types().size() + unpatchable_residue_types().size() << " residue types" << std::endl;
	TR << "Total time to initialize " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;
}

bool sort_patchop_by_name( PatchOP p, PatchOP q ) {
	return ( p->name() < q->name() );
}

void
GlobalResidueTypeSet::init_restypes_from_database() {

	using namespace basic::options;

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
		// to make sure even applications that use GlobalResidueTypeSet directly never run into problems
		// AMW: I have to exclude these purely for the Matcher--could a Matcher expert weigh in?
		// At least they all work with mm_params etc. now.
		// This just means that the "enumerate RTs" app should run with pH mode
		bool no_proton_states = false;
		if ( line.size() > 20 ) {
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
		if ( (!option[OptionKeys::in::include_surfaces]) &&
				(line.substr(0, 29) == "residue_types/mineral_surface") ) {
			continue;
		}

		// Parse lines.
		std::istringstream l( line );
		l >> tag;
		// Note on thread safety. This method is (probably) being called during RTS initialization within
		// the ChemicalManager. As long as we don't attempt to get a RTS from the CM we should avoid deadlock issues,
		// even as we access other type sets.
		if ( tag == "TYPE_SET_MODE" || tag == "TYPE_SET_CATEGORY" ) {
			l >> tag;
			mode( type_set_mode_from_string( tag ) );
		} else if ( tag == "ATOM_TYPE_SET" ) {
			l >> tag;
			atom_type_set( ChemicalManager::get_instance()->atom_type_set( tag ) );
		} else if ( tag == "ELEMENT_SET" ) {
			l >> tag;
			element_set( ChemicalManager::get_instance()->element_set( tag ) );
		} else if ( tag == "MM_ATOM_TYPE_SET" ) {
			l >> tag;
			mm_atom_type_set( ChemicalManager::get_instance()->mm_atom_type_set( tag ) );
		} else if ( tag == "ORBITAL_TYPE_SET" ) {
			l >> tag;
			orbital_type_set( ChemicalManager::get_instance()->orbital_type_set(tag) );
		} else {
			std::string const filename( database_directory_ + line );

			ResidueTypeOP rsd_type( read_topology_file(
				filename, atom_type_set(), element_set(), mm_atom_type_set(), orbital_type_set() ) );

			add_base_residue_type( rsd_type );
			TR.Debug << "Loading ResidueType " << rsd_type->name() << " from file " << filename << std::endl;
		}
	}

}

void GlobalResidueTypeSet::init_restypes_from_commandline() {

	utility::vector1< std::string > extra_params_files( params_files_from_commandline() );

	for ( std::string const & filename : extra_params_files ) {
		ResidueTypeOP rsd_type( read_topology_file(
			filename, atom_type_set(), element_set(), mm_atom_type_set(), orbital_type_set() ) );
		add_base_residue_type( rsd_type );
		TR.Debug << "Loading ResidueType " << rsd_type->name() << " from file " << filename << std::endl;
	}

	if ( mode () == FULL_ATOM_t ) {
		utility::vector1< ResidueTypeOP > extra_residues( extra_nonparam_restypes_from_commandline() );
		for ( ResidueTypeOP const & rsd_type : extra_residues ) {
			add_base_residue_type( rsd_type );
		}
	} else if ( mode() == CENTROID_t ) {
		// Add code to load as fa, convert to centroid, then add to the set
		// See GitHub PR # 1311
	} else {
		// Add code to load as fa, convert to appropriate type set, then add to the set.
	}

	if ( mode() == FULL_ATOM_t ) {
		load_residue_types_from_sql_database();
	} // Sad that we can't load non-fullatom types from the database.

}

utility::vector1< std::string >
GlobalResidueTypeSet::params_files_from_commandline() const {
	using namespace basic::options;

	utility::vector1< std::string > extra_params_files;

	if ( mode() == FULL_ATOM_t ) {

		if ( option[OptionKeys::in::file::extra_res_fa].active() ) {
			for ( utility::file::FileName const & fname : option[OptionKeys::in::file::extra_res_fa] ) {
				extra_params_files.push_back(fname.name());
			}
		}

		if ( option[OptionKeys::in::file::extra_res_path].active() ) {
			for ( utility::file::PathName const & pname : option[OptionKeys::in::file::extra_res_path] ) {
				std::string directory=pname.name();

				utility::vector1<std::string> files;
				utility::file::list_dir(directory, files);
				TR.Debug << "Loading params files from " << directory << ": " << std::endl;
				for ( std::string const & file : files ) {
					if ( file.find("param") != std::string::npos ) {
						TR.Debug << file<< ", ";
						std::string path = directory+'/'+file;
						extra_params_files.push_back(path);
					}
				}
				TR.Debug<< std::endl;
			}
		}

		if ( option[OptionKeys::in::file::extra_res_batch_path].active() ) {
			for ( utility::file::PathName const & pname : option[OptionKeys::in::file::extra_res_batch_path] ) {
				std::string directory=pname.name();

				utility::vector1<std::string> subdirs;
				utility::file::list_dir(directory, subdirs);
				TR.Debug << "Loading batch params files from " << directory << ": " << std::endl;
				for ( std::string const & subdir : subdirs ) {
					if ( subdir == "." || subdir == ".." ) {
						continue;
					}
					utility::vector1<std::string> files;
					utility::file::list_dir(directory+"/"+subdir,files);
					for ( std::string const & file : files ) {
						if ( file.find("param")!=std::string::npos ) {
							TR.Debug << file << ", ";
							std::string path = directory+'/'+subdir+'/'+file;
							extra_params_files.push_back(path);
						}
					}
				}
				TR.Debug << std::endl;
			}
		}

	} else if ( mode() == CENTROID_t ) {

		if ( option[OptionKeys::in::file::extra_res_cen].active() ) {
			for ( utility::file::FileName const & fname : option[OptionKeys::in::file::extra_res_cen] ) {
				extra_params_files.push_back(fname.name());
			}
		}

	}

	// generically specify extra res (not necessarily part of fa_standard) -- will get added to
	//  any and every residue_type_set instantiated.
	if ( option[OptionKeys::in::file::extra_res].active() ) {
		for ( utility::file::FileName const & fname : option[OptionKeys::in::file::extra_res] ) {
			extra_params_files.push_back(fname.name());
		}
	}

	return extra_params_files;
}

utility::vector1< ResidueTypeOP >
GlobalResidueTypeSet::extra_nonparam_restypes_from_commandline() const {
	using namespace  basic::options;

	utility::vector1< ResidueTypeOP > extra_residues;

	// Regardless of what the current ResidueTypeSet is, we want to load these residues in full atom mode
	// It's the responsibility of the caller to convert, if necessary.
	//
	// Note on thread safety. This method is (probably) being called during RTS initialization within
	// the ChemicalManager. As long as we don't attempt to get a RTS from the CM we should avoid deadlock issues,
	// even as we access other type sets.
	core::chemical::AtomTypeSetCOP atom_types = ChemicalManager::get_instance()->atom_type_set("fa_standard");
	core::chemical::ElementSetCOP elements = ChemicalManager::get_instance()->element_set("default");
	core::chemical::MMAtomTypeSetCOP mm_atom_types = ChemicalManager::get_instance()->mm_atom_type_set("fa_standard");
	//core::chemical::orbitals::OrbitalTypeSetCOP orbital_types = ChemicalManager::get_instance()->orbital_type_set("fa_standard");

	// Read extra ResidueTypes from SDF files
	if ( option[OptionKeys::in::file::extra_res_mol].active() ) {
		sdf::MolFileIOReader molfile_reader;
		for ( utility::file::FileName const & filename : option[OptionKeys::in::file::extra_res_mol] ) {
			utility::vector1< sdf::MolFileIOMoleculeOP > data( molfile_reader.parse_file( filename ) );
			utility::vector1< ResidueTypeOP > rtvec( sdf::convert_to_ResidueTypes( data, /* load_rotamers= */ true, atom_types, elements, mm_atom_types ) );
			TR << "Reading " << rtvec.size() << " residue types from the " << data.size() << " models in " << filename << std::endl;
			extra_residues.append( rtvec );
		}
	}

	// Read extra ResidueTypes from mmCIF files
	if ( option[OptionKeys::in::file::extra_res_mmCIF].active() ) {
		mmCIF::mmCIFParser mmCIF_parser;
		for ( utility::file::FileName const & filename : option[OptionKeys::in::file::extra_res_mmCIF] ) {
			utility::vector1< sdf::MolFileIOMoleculeOP> molecules( mmCIF_parser.parse( filename ) );
			utility::vector1< ResidueTypeOP > rtvec( sdf::convert_to_ResidueTypes( molecules, true, atom_types, elements, mm_atom_types ) );
			extra_residues.append( rtvec );
		}
	}

	return extra_residues;
}

void
GlobalResidueTypeSet::load_residue_types_from_sql_database() {

	debug_assert( mode() == FULL_ATOM_t );

	if ( basic::options::option[basic::options::OptionKeys::in::file::extra_res_database].user() ) {
		utility::sql_database::DatabaseMode::e database_mode(
			utility::sql_database::database_mode_from_name(
			basic::options::option[basic::options::OptionKeys::in::file::extra_res_database_mode]));
		std::string database_name(basic::options::option[basic::options::OptionKeys::in::file::extra_res_database]);
		std::string database_pq_schema(basic::options::option[basic::options::OptionKeys::in::file::extra_res_pq_schema]);


		utility::sql_database::sessionOP db_session(
			basic::database::get_db_session(database_mode, database_name, database_pq_schema));

		ResidueDatabaseIO residue_database_interface;

		if ( basic::options::option[basic::options::OptionKeys::in::file::extra_res_database_resname_list].user() ) {
			utility::file::FileName residue_list = basic::options::option[basic::options::OptionKeys::in::file::extra_res_database_resname_list];
			utility::io::izstream residue_name_file(residue_list);
			std::string residue_name;
			while ( residue_name_file >> residue_name )
					{
				//residue_name_file >> residue_name;
				TR <<residue_name <<std::endl;
				ResidueTypeOP new_residue(
					residue_database_interface.read_residuetype_from_database(
					atom_type_set(),
					element_set(),
					mm_atom_type_set(),
					orbital_type_set(),
					"fa_standard",  // NOTE: We're only running this function when the mode is full atom
					residue_name,
					db_session));

				add_base_residue_type( new_residue );
			}

		} else {
			utility::vector1<std::string> residue_names_in_database( residue_database_interface.get_all_residues_in_database(db_session));
			for ( Size index =1; index <= residue_names_in_database.size(); ++index ) {
				ResidueTypeOP new_residue(
					residue_database_interface.read_residuetype_from_database(
					atom_type_set(),
					element_set(),
					mm_atom_type_set(),
					orbital_type_set(),
					"fa_standard",
					residue_names_in_database[index],
					db_session));

				add_base_residue_type( new_residue );
			}
		}
	}
}

void GlobalResidueTypeSet::init_patches_from_database() {
	using namespace basic::options;

	utility::vector1< std::string > patch_filenames;
	utility::vector1< std::string > metapatch_filenames;

	// Load patches

	// Read the command line and avoid applying patches that the user has requested be
	// skipped.  The flag allows the user to specify a patch by its name or by its file.
	// E.g. "SpecialRotamer" or "SpecialRotamer.txt".  Directory names will be ignored if given.
	std::set< std::string > patches_to_avoid;
	if ( option[ OptionKeys::chemical::exclude_patches ].user() ) {
		for ( utility::file::FileName const & fname : option[ OptionKeys::chemical::exclude_patches ] ) {
			patches_to_avoid.insert( fname.base() );
		}
	}

	std::string const list_filename( database_directory_+"/patches.txt" );
	utility::io::izstream data( list_filename.c_str() );

	if ( !data.good() ) {
		utility_exit_with_message( "Unable to open patch list file: "+list_filename );
	}

	// Unconditional loading of listed patches is deliberate --
	// if you specified it explicitly, you probably want it to load.
	std::string line;
	while ( getline( data,line) ) {
		if ( line.size() < 1 || line[0] == '#' ) continue;

		// get rid of any comment lines.
		line = utility::string_split( line, '#' )[1];
		line = utility::string_split( line, ' ' )[1];

		// Skip carbohydrate patches unless included with include_sugars flag.
		if ( ( ! option[ OptionKeys::in::include_sugars ] ) &&
				( line.substr( 0, 21 ) == "patches/carbohydrates" ) ) {
			continue;
		}

		// Skip this patch if the "patches_to_avoid" set contains the named patch.
		// AMW: keeping this because "patches_to_avoid" is explicitly asked for
		// on the command line.
		utility::file::FileName fname( line );
		if ( patches_to_avoid.find( fname.base() ) != patches_to_avoid.end() ) {
			TR << "While generating GlobalResidueTypeSet " << name() <<
				": Skipping patch " << fname.base() << " as requested" << std::endl;
			continue;
		}

		patch_filenames.push_back( database_directory_ + line );
	}

	// kdrew: include list allows patches to be included while being commented out in patches.txt,
	// useful for testing non-canonical patches.
	// Retaining this just because maybe you haven't put your patch in the list yet
	if ( option[ OptionKeys::chemical::include_patches ].active() ) {
		for ( utility::file::FileName const & fname : option[ OptionKeys::chemical::include_patches ] ) {
			if ( !utility::file::file_exists( database_directory_ + fname.name() ) ) {
				TR.Warning << "Could not find: " << database_directory_ + fname.name()  << std::endl;
				continue;
			}
			patch_filenames.push_back( database_directory_ + fname.name() );
			TR << "While generating GlobalResidueTypeSet " << name() <<
				": Including patch " << fname << " as requested" << std::endl;
		}
	}

	//fpd  if missing density is to be read correctly, we will have to also load the terminal truncation variants
	if ( option[ OptionKeys::in::missing_density_to_jump ]()
			|| option[ OptionKeys::in::use_truncated_termini ]() ) {
		if ( std::find( patch_filenames.begin(), patch_filenames.end(), database_directory_ + "patches/NtermTruncation.txt" )
				== patch_filenames.end() ) {
			patch_filenames.push_back( database_directory_ + "patches/NtermTruncation.txt" );
		}
		if ( std::find( patch_filenames.begin(), patch_filenames.end(), database_directory_ + "patches/CtermTruncation.txt" )
				== patch_filenames.end() ) {
			patch_filenames.push_back( database_directory_ + "patches/CtermTruncation.txt" );
		}

	}

	// Also obtain metapatch filenames
	std::string const meta_filename( database_directory_+"/metapatches.txt" );
	utility::io::izstream data2( meta_filename.c_str() );

	if ( data2.good() ) { // Not all typeset directories have metapatches
		std::string mpline;
		while ( getline( data2, mpline ) ) {
			if ( mpline.size() < 1 || mpline[0] == '#' ) continue;
			metapatch_filenames.push_back( database_directory_ + mpline );
		}
	} else {
		TR.Debug << "Skipping metapatch loading for " << name() << " as metapatches.txt can't be found. " << std::endl;
	}

	add_patches( patch_filenames, metapatch_filenames );
}

void GlobalResidueTypeSet::init_patches_from_commandline() {
	using namespace basic::options;

	utility::vector1< std::string > extra_patch_files;

	if ( mode() == FULL_ATOM_t ) {
		if ( option[OptionKeys::in::file::extra_patch_fa].active() ) {
			for ( utility::file::FileName const & fname : option[OptionKeys::in::file::extra_patch_fa] ) {
				extra_patch_files.push_back( fname.name());
			}
		}
	} else if ( mode() == CENTROID_t ) {
		if ( option[OptionKeys::in::file::extra_patch_cen].active() ) {
			for ( utility::file::FileName const & fname : option[OptionKeys::in::file::extra_patch_cen] ) {
				extra_patch_files.push_back( fname.name());
			}
		}
	}

	utility::vector1< std::string > extra_metapatch_files;

	add_patches( extra_patch_files, extra_metapatch_files );
}

void
GlobalResidueTypeSet::deal_with_patch_special_cases()
{
	// separate this to handle a set of base residue types and a set of patches...
	// this would allow addition of patches and/or base_residue_types at stages after initialization.

	// The "replace_residue" patches are a special case, and barely in use anymore.
	// In their current implementation, they actually do *not* change the name of the ResidueType.
	// But we probably should keep track of their application
	//  by updating name() of residue; and hold copies of the replaced residues without patch applied in, e.g., replaced_name_map_.
	// That's going to require a careful refactoring of how residue types are accessed (e.g., can't just use name3_map() anymore),
	//  which I might do later. For example, there's code in SwitchGlobalResidueTypeSet that will look for "MET" when it really should
	//  look for "MET:protein_centroid_with_HA" and know about this patch.
	// For now, apply them first, and force application/instantiation later.
	// -- rhiju
	for ( PatchCOP p : patches() ) {
		for ( ResidueTypeCOP rsd_type : base_residue_types() ) {
			if ( p->applies_to( *rsd_type ) ) {
				if ( p->replaces( *rsd_type ) ) {
					runtime_assert( rsd_type->finalized() );
					ResidueTypeCOP rsd_type_new = p->apply( *rsd_type );
					cache_object()->update_residue_type( rsd_type, rsd_type_new );
					runtime_assert( update_base_residue_types_if_replaced( rsd_type, rsd_type_new ) );
				}
			}
		}
	}

	// separate this to handle a set of base residue types and a set of patches.
	// this would allow addition of patches and/or base_residue_types at stages after initialization.
	for ( PatchCOP p : patches() ) {
		for ( ResidueTypeCOP rsd_type : base_residue_types() ) {
			if ( p->applies_to( *rsd_type ) && p->adds_properties( *rsd_type ).has_value( "D_AA" ) ) {
				ResidueTypeOP new_rsd_type = p->apply( *rsd_type );
				new_rsd_type->base_name( new_rsd_type->name() ); //D-residues have their own base names.
				new_rsd_type->reset_base_type_cop(); //This is now a base type, so its base type pointer must be NULL.

				add_base_residue_type( new_rsd_type );

				// Store the D-to-L and L-to-D mappings:
				runtime_assert_string_msg(
					l_to_d_mapping().count( rsd_type ) == 0,
					"Error in core::chemical::ResidueTypeSet::apply_patches: A D-equivalent for " + rsd_type->name() + " has already been defined."
				);
				l_to_d_mapping()[ rsd_type ] = new_rsd_type;
				runtime_assert_string_msg(
					d_to_l_mapping().count( new_rsd_type ) == 0,
					"Error in core::chemical::ResidueTypeSet::apply_patches: An L-equivalent for " + new_rsd_type->name() + " has already been defined."
				);
				d_to_l_mapping()[ new_rsd_type ] = rsd_type;
			}

			if ( p->applies_to( *rsd_type ) && p->adds_properties( *rsd_type ).has_value( "L_RNA" ) ) {
				ResidueTypeOP new_rsd_type = p->apply( *rsd_type );
				new_rsd_type->base_name( new_rsd_type->name() ); //L-RNA residues have their own base names.
				new_rsd_type->reset_base_type_cop(); //This is now a base type, so its base type pointer must be NULL.

				add_base_residue_type( new_rsd_type );

				// Store the D-to-L and L-to-D mappings -- in reverse, of course!
				runtime_assert_string_msg(
					l_to_d_mapping().count( new_rsd_type ) == 0,
					"Error in core::chemical::ResidueTypeSet::apply_patches: A D-equivalent for " + rsd_type->name() + " has already been defined." );
				l_to_d_mapping()[ new_rsd_type ] = rsd_type;
				runtime_assert_string_msg(
					d_to_l_mapping().count( rsd_type ) == 0,
					"Error in core::chemical::ResidueTypeSet::apply_patches: An L-equivalent for " + new_rsd_type->name() + " has already been defined." );
				d_to_l_mapping()[ rsd_type ] = new_rsd_type;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
/// @details Generation of new residue types augmented by adduct atoms
/// @note    This is almost perfectly consistent with on_the_fly residue type sets.
///            Why not just make adducts a patch? Then it will all work together nicely.
///            Just need to have the -adduct string (add_map) turned into a PatchSelector
///                 -- rhiju
void
GlobalResidueTypeSet::place_adducts() {

	// First parse the command line for requested adducts
	utility::options::StringVectorOption & add_set
		= basic::options::option[ basic::options::OptionKeys::packing::adducts ];

	// No adducts, skip out
	if ( add_set.size() == 0 ) return;

	// Convert to a map that takes a string descriptor of an adduct and
	// gives the max number of adducts of that class to apply, so
	// a command line option of -adducts <adduct_type> 2 for a type that has
	// 5 entries in the rsd param file will create all combinations with up to
	// 2 total adducts.

	AdductMap add_map = parse_adduct_string( add_set );

	// Error check each requested adduct from the command line, and
	// complain if there are no examples in any residues.  This function
	// will not return if
	error_check_requested_adducts( add_map, base_residue_types() );

	// Set up a starting point map where the int value is the number
	// of adducts of a given type placed
	AdductMap blank_map( add_map );
	for ( auto & add_iter : blank_map ) {
		add_iter.second = 0;
	}

	// Process the residues in turn
	ResidueTypeCOPs residue_types = ResidueTypeFinder( *this ).base_property( DNA ).variant_exceptions( utility::tools::make_vector1( LOWER_TERMINUS_VARIANT, UPPER_TERMINUS_VARIANT, ADDUCT_VARIANT ) ).get_all_possible_residue_types();

	for ( ResidueTypeCOP const & rsd : residue_types ) {
		if ( !rsd || !rsd->finalized() ) continue;
		AdductMap count_map( blank_map );
		utility::vector1< bool > add_mask( rsd->defined_adducts().size(), false  );
		for ( ResidueTypeOP const & newtype : create_adduct_combinations( *rsd, add_map, count_map, add_mask, rsd->defined_adducts().begin() ) ) {
			add_unpatchable_residue_type( newtype );
		}
	}
}

GlobalResidueTypeSet::~GlobalResidueTypeSet() {}

//////////////////////////////////////////////////////////////////////////////
// The useful stuff:  Accessor functions
//////////////////////////////////////////////////////////////////////////////

/// @brief   checks if name exists.
/// @details actually instantiates the residue type if it does not exist.
bool
GlobalResidueTypeSet::has_name( std::string const & name ) const
{
	return generate_residue_type( name );
}

/// @brief Gets all types with the given aa type and variants
/// @details The number of variants must match exactly.
/// (It's assumed that the passed VariantTypeList contains no duplicates.)
ResidueTypeCOPs
GlobalResidueTypeSet::get_all_types_with_variants_aa( AA aa, utility::vector1< std::string > const & variants ) const
{
	utility::vector1< VariantType > exceptions;
	return cache_object()->get_all_types_with_variants_aa( aa, variants, exceptions );
}

ResidueTypeCOPs
GlobalResidueTypeSet::get_all_types_with_variants_aa( AA aa,
	utility::vector1< std::string > const & variants,
	utility::vector1< VariantType > const & exceptions ) const
{
	return cache_object()->get_all_types_with_variants_aa( aa, variants, exceptions );
}

void
GlobalResidueTypeSet::remove_base_residue_type( std::string const & )
{
	utility_exit_with_message("Shouldn't remove types from a GlobalResidueTypeSet!!");
}

void
GlobalResidueTypeSet::remove_unpatchable_residue_type( std::string const & )
{
	utility_exit_with_message("Shouldn't remove types from a GlobalResidueTypeSet!!");
}

/// @brief From a file, read which IDs shouldn't be loaded from the components.
void
GlobalResidueTypeSet::load_shadowed_ids( std::string const & directory, std::string const & filename /* = "shadow_list.txt" */ ) {

	TR.Debug << "Loading shadowed PDB IDs from " << directory + filename << std::endl;

	shadowed_ids_.clear();

	utility::io::izstream file( directory + filename );
	if ( ! file.good()  ) {
		TR << "For ResidueTypeSet " << name() << " there is no " << filename << " file to list known PDB ids." << std::endl;
		TR << "    This will turn off PDB component loading for ResidueTypeSet " << name() << std::endl;
		TR << "    Expected file: " << directory + filename << std::endl;
		return;
	}
	std::string line;
	getline( file, line );
	while ( file.good() ) {
		utility::trim( line ); // inplace;
		if ( line[0] != '#' ) {
			shadowed_ids_.insert( line );
		}
		getline( file, line );
	}
	if ( shadowed_ids_.size() == 0 ) {
		TR.Warning << "For ResidueTypeSet " << name() << ", " << filename << " doesn't have any entries." << std::endl;
		TR.Warning << "    This will turn off PDB component loading for ResidueTypeSet " << name() << std::endl;
	}
}

/// @brief Attempt to lazily load the given residue type from data.
bool
GlobalResidueTypeSet::lazy_load_base_type( std::string const & rsd_base_name ) const
{
	if ( cache_object()->has_generated_residue_type( rsd_base_name ) ) { return true; }
	if ( cache_object()->is_prohibited( rsd_base_name ) ) { return false; }

	core::chemical::ResidueTypeOP new_rsd_type;

	// These are heuristics to figure out where to load the data from.

	// Heuristic: if the ResidueTypeName begins with 'pdb_', then it's loaded from the chemical components directory
	if ( rsd_base_name.find("pdb_") == 0 ) {
		std::string short_name( utility::strip( rsd_base_name.substr( 4, rsd_base_name.size() ) ) );
		if ( shadowed_ids_.size() > 0 && shadowed_ids_.count( short_name ) == 0 ) {
			new_rsd_type = load_pdb_component( short_name );
			if ( new_rsd_type ) {
				// Duplicate detection is handled by the shadowed file -- if it's not shadowed, we load the component
				new_rsd_type->name( "pdb_" + short_name );
				TR << "Loading '" << short_name << "' from the PDB components dictionary for residue type '" << rsd_base_name << "'" << std::endl;
			}
		} else {
			if ( shadowed_ids_.size() == 0 ) {
				TR.Debug << "Not loading '" << short_name << "' from PDB components dictionary because components are turned off for this ResidueTypeSet." << std::endl;
			} else {
				TR.Debug << "Not loading '" << short_name << "' from PDB components dictionary because it is shadowed in the ResidueTypeSet." << std::endl;
			}
			cache_object()->add_prohibited( rsd_base_name );
			return false;
		}
	}

	// Finish up with the new residue type.
	if ( new_rsd_type ) {
		cache_object()->add_residue_type( new_rsd_type );
	}
	return cache_object()->has_generated_residue_type( rsd_base_name );

}

/// @brief Load a residue type from the components dictionary.
ResidueTypeOP
GlobalResidueTypeSet::load_pdb_component( std::string const & pdb_id ) const {
	static THREAD_LOCAL bool warned_about_missing_file( false );
	if ( pdb_components_filename_.size() ) {
		utility::io::izstream filestream( pdb_components_filename_ );
		if ( !filestream.good() ) {
			std::string db_filename( basic::database::full_name( pdb_components_filename_, false ) );
			filestream.open( db_filename );

			if ( !filestream.good() ) {
				if ( ! warned_about_missing_file ) {
					warned_about_missing_file = true;
					TR.Warning << "PDB component dictionary file not found at (./)" << pdb_components_filename_ << std::endl;
					TR.Warning << "   or in the Rosetta database at " << db_filename << std::endl;
					TR.Warning << "   For information on how to obtain the file and set it for use with Rosetta, visit: \n\n";
					TR.Warning << "  https://www.rosettacommons.org/docs/latest/build_documentation/Build-Documentation#setting-up-rosetta-3_obtaining-additional-files_pdb-chemical-components-dictionary  \n" << std::endl;
				}
				return ResidueTypeOP( nullptr );
			}
		}

		std::string entry( "data_" + pdb_id );
		std::string line;
		std::string lines;
		//std::cout << "Finding '" << entry <<"' " << entry.size() << std::endl;
		mmCIF::mmCIFParser mmCIF_parser;
		while ( filestream.good() ) {
			getline( filestream, line );
			if ( line.size() == entry.size() ) {
				//std::cout << line << std::endl;
				if ( line == entry ) {
					lines += line;
					//lines.push_back( line);
					getline( filestream, line);
					while ( line.substr(0, 5) != "data_" && filestream.good() ) {
						lines += line + '\n';
						//lines.push_back( line);
						getline( filestream, line);
					}
					break;
				}
			}
		}

		if ( lines.size() == 0 ) {
			TR.Warning << "Could not find: '" << pdb_id << "' in pdb components file '" << pdb_components_filename_
				<< "'! Skipping residue..." << std::endl;
			return ResidueTypeOP( nullptr );
		}

		utility::vector1< core::chemical::sdf::MolFileIOMoleculeOP> molecules;
		molecules.push_back( mmCIF_parser.parse( lines, pdb_id) );
		core::chemical::ResidueTypeOP new_rsd_type( core::chemical::sdf::convert_to_ResidueType( molecules ) );

		return new_rsd_type;
	}
	return ResidueTypeOP( nullptr );
}

} // pose
} // core
