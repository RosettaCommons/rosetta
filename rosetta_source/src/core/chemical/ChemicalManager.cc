// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin ChemicalManager
///
/// @brief
/// Chemical manager class
///
/// @detailed
/// The Chemical Manager is a singleton class, which means that it can only been initialized once (exist once in memory). Once initialized,
/// you can call it by simply access it via:
///
/// core::chemical::AtomTypeSetCAP atom_types =
/// core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
///
/// You can substitute AtomTypeSet, with whatever is seen below (residue_type_set, mm_atom_type_set, orbital_type_set).
/// In the below functions, the "tag_in" refers to fullatom, centroid, which basically tells what type of set to load in.
/// The chemical manager will call functions within the AtomTypeSet, MMAtomTypeSet, ResidueTypeSet, etc etc. The classes type set
/// reads in files from the database to create atom types, residue types, and mmatom types. The information from those files are stored
/// in the type class.
///
///
///
/// @authors
/// Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// Steven Combs - comments
///
///
/// @last_modified December 6 2010
/////////////////////////////////////////////////////////////////////////


#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/sdf/mol_parser.hh>
#include <core/chemical/ResidueDatabaseIO.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
// Project headers
#include <basic/database/sql_utils.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <core/chemical/ResidueType.hh>



namespace core {
namespace chemical {

static basic::Tracer tr("core.chemical.ChemicalManager");

/// @brief set initial value as no instance
ChemicalManager* ChemicalManager::instance_( 0 );

/// @brief static function to get the instance of ( pointer to) this singleton class
ChemicalManager * ChemicalManager::get_instance()
{
	if ( instance_ == 0 )
	{
		 instance_ = new ChemicalManager();
	}
	return instance_;
}

/// @brief private constructor to guarantee the singleton
ChemicalManager::ChemicalManager(){}

/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
AtomTypeSetCAP
ChemicalManager::atom_type_set( std::string const & tag )
{
	AtomTypeSets::const_iterator iter( atom_type_sets_.find( tag ) );
	if ( iter == atom_type_sets_.end() ) {
		// read from file
		std::string const directory( basic::database::full_name( "chemical/atom_type_sets/"+tag+"/" ) );
		AtomTypeSetOP new_set( new AtomTypeSet( directory ) );
		iter = atom_type_sets_.insert( std::make_pair( tag, new_set ) ).first;
	}
	return iter->second();
}

/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
ElementSetCAP
ChemicalManager::element_set( std::string const & tag )
{
	ElementSets::const_iterator iter( element_sets_.find( tag ) );
	if ( iter == element_sets_.end() ) {
		// read from file
		std::string const filename( basic::database::full_name( "chemical/element_sets/"+tag+"/element_properties.txt" ));
		ElementSetOP new_set( new ElementSet() );
		new_set->read_file( filename );
		iter = element_sets_.insert( std::make_pair( tag, new_set ) ).first;
	}
	return iter->second();
}

/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
orbitals::OrbitalTypeSetCAP
ChemicalManager::orbital_type_set( std::string const & tag )
{
	OrbitalTypeSets::const_iterator iter( orbital_type_sets_.find( tag ) );
	if ( iter == orbital_type_sets_.end() ) {
		// read from file
		std::string const directory( basic::database::full_name( "chemical/orbital_type_sets/"+tag+"/" ) );
		orbitals::OrbitalTypeSetOP new_set( new orbitals::OrbitalTypeSet( directory ) );
		iter = orbital_type_sets_.insert( std::make_pair( tag, new_set ) ).first;
	}
	return iter->second();
}

/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
MMAtomTypeSetCAP
ChemicalManager::mm_atom_type_set( std::string const & tag )
{
	MMAtomTypeSets::const_iterator iter( mm_atom_type_sets_.find( tag ) );
	if ( iter == mm_atom_type_sets_.end() ) {
		// read from file
		std::string const filename( basic::database::full_name( "chemical/mm_atom_type_sets/"+tag+"/mm_atom_properties.txt" ));
		MMAtomTypeSetOP new_set( new MMAtomTypeSet() );
		new_set->read_file( filename );
		iter = mm_atom_type_sets_.insert( std::make_pair( tag, new_set ) ).first;
	}
	return iter->second();
}


// /// @details if the tag is not in the map, input it from a database file and add it
// /// to the map for future look-up.
//CSDAtomTypeSetCAP
//ChemicalManager::csd_atom_type_set( std::string const & tag )
//{
//	CSDAtomTypeSets::const_iterator iter( csd_atom_type_sets_.find( tag ) );
//	if ( iter == csd_atom_type_sets_.end() ) {
		// read from file
//		std::string const filename( basic::database::full_name( "chemical/csd_atom_type_sets/"+tag+"/csd_atom_properties.txt" ));
//		CSDAtomTypeSetOP new_set( new CSDAtomTypeSet() );
//		new_set->read_file( filename );
//		iter = csd_atom_type_sets_.insert( std::make_pair( tag, new_set ) ).first;
//	}
//	return iter->second();
//}




///@ details if the tag is not in the map, input it from a database file and add it
///to the map for future look-up.
ResidueTypeSetCAP
ChemicalManager::residue_type_set( std::string const & tag )
{

	using namespace basic;

	ResidueTypeSets::const_iterator iter( residue_type_sets_.find( tag ) );
	if ( iter == residue_type_sets_.end() ) {
		// Look for additional residue .params files specified on the cmd line
		std::vector< std::string > extra_params_files;
		std::vector<core::chemical::ResidueTypeOP> extra_residues;
		if(tag == FA_STANDARD) {
			utility::options::FileVectorOption & fvec
				= basic::options::option[ basic::options::OptionKeys::in::file::extra_res_fa ];
			for(Size i = 1, e = fvec.size(); i <= e; ++i) {
				utility::file::FileName fname = fvec[i];
				extra_params_files.push_back(fname.name());
			}

			utility::options::PathVectorOption & pvec
				= basic::options::option[basic::options::OptionKeys::in::file::extra_res_path];
			// convert Pathname->string->char*, glob it, convert char*->string
			for(Size i=1, e= pvec.size(); i<=e; i++){
				utility::vector1<std::string> files;
				std::string directory=pvec[i].name();

				utility::file::list_dir(directory, files);
				tr.Debug<< std::endl;
				for(size_t j=1; j<= files.size(); j++){
					if (files[j].find("param")!=std::string::npos){
						tr.Debug << files[j]<< ", ";
						std::string path= directory+'/'+files[j];
						extra_params_files.push_back(path);
					}
				}
				tr.Debug<< std::endl;
			}

			utility::options::FileVectorOption & mdlvec
				= basic::options::option[basic::options::OptionKeys::in::file::extra_res_mol];
			core::chemical::AtomTypeSetCAP atom_types = atom_type_set("fa_standard");
			core::chemical::ElementSetCAP elements = element_set("fa_standard");
			core::chemical::MMAtomTypeSetCAP mm_atom_types = mm_atom_type_set("fa_standard");
			core::chemical::orbitals::OrbitalTypeSetCAP orbital_types = orbital_type_set("fa_standard");
			for(Size i=1, e = mdlvec.size(); i <= e;++i)
			{
				utility::file::FileName filename = mdlvec[i];
				core::chemical::sdf::MolFileParser parser(filename.name());
				parser.parse_mol_file(atom_types, elements, mm_atom_types, orbital_types);
				extra_residues.push_back(parser.GetResidueTypeOP());
			}

			if(basic::options::option[basic::options::OptionKeys::in::file::extra_res_database].user())
			{
				std::string database_name = basic::options::option[basic::options::OptionKeys::in::file::extra_res_database];
				std::string database_mode = basic::options::option[basic::options::OptionKeys::in::file::extra_res_database_mode];

				utility::sql_database::sessionOP db_session;

				if(database_mode=="sqlite3")
				{
					db_session = basic::database::get_db_session(database_name,database_mode,true,false);
				}else
				{
					db_session = basic::database::get_db_session(database_name,database_mode,false,false);
				}
				ResidueDatabaseIO residue_database_interface;

				if(basic::options::option[basic::options::OptionKeys::in::file::extra_res_database_resname_list].user())
				{
					utility::file::FileName residue_list = basic::options::option[basic::options::OptionKeys::in::file::extra_res_database_resname_list];
					utility::io::izstream residue_name_file(residue_list);
					std::string residue_name;
					while(residue_name_file >> residue_name)
					{
						//residue_name_file >> residue_name;
						//tr <<residue_name <<std::endl;
						ResidueTypeOP new_residue(
							residue_database_interface.read_residuetype_from_database(
								atom_types,
								elements,
								mm_atom_types,
								orbital_types,
								"fa_standard",
								residue_name,
								db_session));
						extra_residues.push_back(new_residue);
					}

				}else
				{
					utility::vector1<std::string> residue_names_in_database( residue_database_interface.get_all_residues_in_database(db_session));
					for(Size index =1; index <= residue_names_in_database.size();++index)
					{
						ResidueTypeOP new_residue(
							residue_database_interface.read_residuetype_from_database(
								atom_types,
								elements,
								mm_atom_types,
								orbital_types,
								"fa_standard",
								residue_names_in_database[index],
								db_session));
						extra_residues.push_back(new_residue);
					}
				}
			}



		} else if(tag == CENTROID) {
			utility::options::FileVectorOption & fvec
				= basic::options::option[ basic::options::OptionKeys::in::file::extra_res_cen ];
			for(Size i = 1, e = fvec.size(); i <= e; ++i) {
				utility::file::FileName fname = fvec[i];
				extra_params_files.push_back(fname.name());
			}
		}
		// read from file
		tr.Debug << "CHEMICAL_MANAGER: read residue types: " << tag << std::endl;
		// redirecting to new icoor folder
		std::string temp_str( basic::database::full_name( "chemical/residue_type_sets/"+tag ) );
		if(tag == FA_STANDARD) {
			if (basic::options::option[basic::options::OptionKeys::corrections::chemical::icoor_05_2009]) {
				temp_str += "_05.2009_icoor";
			}
		}
		temp_str += "/";



		std::string const directory( temp_str );
		ResidueTypeSetOP new_set( new ResidueTypeSet( tag, directory, extra_params_files ) );
		ResidueTypeSetCAP new_setCAP(*new_set);

		for(core::Size index =0 ;index < extra_residues.size();++index)
		{
			//tr << extra_residues[index]->name3() <<std::endl;
			new_set->add_residue_type(extra_residues[index]);
			extra_residues[index]->residue_type_set(new_setCAP);
			ResidueTypeSetCAP new_set_cap(new_set.get());
			extra_residues[index]->residue_type_set(new_set_cap);
		}

		iter = residue_type_sets_.insert( std::make_pair( tag, new_set ) ).first;

	}
	return iter->second();
}


///@ details if the tag is not in the map, input it from a database file and add it
///to the map for future look-up.
ResidueTypeSet &
ChemicalManager::nonconst_residue_type_set( std::string const & tag )
{
	// trigger initialization if necessary:
	residue_type_set( tag );

	return *( residue_type_sets_.find( tag )->second );
}



// global data
/// @brief tag name for querying fullatom chemical type set.
std::string const FA_STANDARD( "fa_standard" );
/// @brief tag name for querying centroid chemical type set.
std::string const CENTROID( "centroid" );
/// @brief tag name for querying coarse-grained chemical type set.
std::string const COARSE_TWO_BEAD( "coarse_two_bead" );
/// @brief tag name for querying hybrid fullatom+centroid chemical type set.
std::string const HYBRID_FA_STANDARD_CENTROID( "hybrid_fa_standard_centroid" );
/// @brief tag name for querying RNA chemical type set.
std::string const RNA( "rna" );
/// @brief tag name for querying COARSE_RNA chemical type set.
std::string const COARSE_RNA( "coarse_rna" );

} // namespace core
} // namespace chemical


/// THIS TURNED OUT TO BE MORE TROUBLE THAN IT WAS WORTH: but maybe some day...


/**
/// @details  Duplicate a ResidueTypeSet, preparatory to modifying it in some way DOES NOT DUPLICATE ITS ATOMTYPESETS
/// Uses ResidueTypeSet::clone()

void
ChemicalManager::copy_residue_type_set(
																			 std::string const & old_name,
																			 std::string const & new_name
																			 )
{
	residue_type_set( old_name ); // triggers initialization if necessary
	if ( residue_type_sets_.find( new_name ) ) {
		utility_exit_with_message( "new name is already being used!" );
	}
	residue_type_sets_.insert( std::make_pair( new_name, residue_type_sets_.find( old_name )->second->clone() ) );
}

/// @details  Duplicate an AtomTypeSet, preparatory to modifying it in some way
/// Uses AtomTypeSet::clone()

void
ChemicalManager::copy_atom_type_set(
																		std::string const & old_name,
																		std::string const & new_name
																		)
{
	atom_type_set( old_name ); // triggers initialization if necessary
	if ( atom_type_sets_.find( new_name ) ) {
		utility_exit_with_message( "new name is already being used!" );
	}
	atom_type_sets_.insert( std::make_pair( new_name, atom_type_sets_.find( old_name )->second->clone() ) );
}
**/

