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
#include <core/chemical/IdealBondLengthSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/sdf/MolFileIOReader.hh>
#include <core/chemical/sdf/MolFileIOData.hh>
#include <core/chemical/ResidueDatabaseIO.hh>
#include <core/chemical/util.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
// Project headers
#include <basic/database/sql_utils.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/sql_database/types.hh>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <core/chemical/ResidueType.hh>

#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>


// Singleton instance and mutex static data members
namespace utility {

using core::chemical::ChemicalManager;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< ChemicalManager > ::singleton_mutex_;
template <> std::atomic< ChemicalManager * > utility::SingletonBase< ChemicalManager >::instance_( 0 );
#else
template <> ChemicalManager * utility::SingletonBase< ChemicalManager >::instance_( 0 );
#endif

}


namespace core {
namespace chemical {

static thread_local basic::Tracer TR( "core.chemical.ChemicalManager" );

/// @brief private constructor to guarantee the singleton
ChemicalManager::ChemicalManager(){}

ChemicalManager *
ChemicalManager::create_singleton_instance()
{
	return new ChemicalManager;
}

/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
AtomTypeSetCOP
ChemicalManager::atom_type_set( std::string const & tag )
{
	AtomTypeSets::const_iterator iter;
	{ // scope for the read lock guard
#if defined MULTI_THREADED && CXX11
		utility::thread::ReadLockGuard lock( atomtype_mutex_ );
#endif
		iter = atom_type_sets_.find( tag );
	}
	if ( iter == atom_type_sets_.end() ) {

		// bind the relevant atom-type-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< AtomTypeSet > () > func =
				boost::bind( &ChemicalManager::create_atom_type_set, this, boost::cref(tag) );

#if defined MULTI_THREADED && defined CXX11
		iter = utility::thread::create_and_insert( func, atomtype_mutex_, tag, atom_type_sets_ );
#else
		AtomTypeSetOP newset = func();
		iter = atom_type_sets_.insert( std::make_pair( tag, newset ) ).first;
#endif
	}
	return iter->second;
}

///@details Actually go and create an AtomTypeSet
AtomTypeSetOP
ChemicalManager::create_atom_type_set(
		std::string const & tag
) const
{
	std::string const directory( basic::database::full_name( "chemical/atom_type_sets/"+tag+"/" ) );
	AtomTypeSetOP new_set( new AtomTypeSet( directory ) );
	// optionally add extra parameters from files given on the command line (see util.hh)
	modify_atom_properties_from_command_line( tag, *new_set );
	// optionally add extra parameters from files given on the command line (see util.hh)
	add_atom_type_set_parameters_from_command_line( tag, *new_set );
	return new_set;
}


/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
ElementSetCOP
ChemicalManager::element_set( std::string const & tag )
{
	ElementSets::const_iterator iter;
	{ // scope for the ReadLockGuard
#if defined MULTI_THREADED && defined CXX11
		utility::thread::ReadLockGuard lock( elem_mutex_ );
#endif
		iter = element_sets_.find( tag );
	}
	if ( iter == element_sets_.end() ) {

		// bind the relevant atom-type-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< ElementSet > () > func =
				boost::bind( &ChemicalManager::create_element_set, this, boost::cref(tag) );

#if defined MULTI_THREADED && defined CXX11
		iter = utility::thread::create_and_insert( func, elem_mutex_, tag, element_sets_ );
#else
		ElementSetOP newset = func();
		iter = element_sets_.insert( std::make_pair( tag, newset ) ).first;
#endif
	}
	return iter->second;
}
ElementSetOP
ChemicalManager::create_element_set( std::string const & tag ) const
{
	// read from file
	std::string const filename( basic::database::full_name( "chemical/element_sets/" + tag + "/element_properties.txt" )); //there are only one type of elements!!!!!!
	ElementSetOP new_set( new ElementSet() );
	new_set->read_file( filename );
	return new_set;
}


/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
orbitals::OrbitalTypeSetCOP
ChemicalManager::orbital_type_set( std::string const & tag )
{
	OrbitalTypeSets::const_iterator iter;
	{
#if defined MULTI_THREADED && CXX11
		utility::thread::ReadLockGuard lock( orbtype_mutex_ );
#endif
		iter = orbital_type_sets_.find( tag );
	}
	if ( iter == orbital_type_sets_.end() ) {

		// bind the relevant atom-type-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< orbitals::OrbitalTypeSet > () > func =
				boost::bind( &ChemicalManager::create_orbital_type_set, this, boost::cref(tag) );

#if defined MULTI_THREADED && defined CXX11
		iter = utility::thread::create_and_insert( func, orbtype_mutex_, tag, orbital_type_sets_ );
#else
		orbitals::OrbitalTypeSetOP newset = func();
		iter = orbital_type_sets_.insert( std::make_pair( tag, newset ) ).first;
#endif
	}
	return iter->second;
}

orbitals::OrbitalTypeSetOP
ChemicalManager::create_orbital_type_set( std::string const & tag ) const
{
	// read from file
	std::string const directory( basic::database::full_name( "chemical/orbital_type_sets/"+tag+"/" ) );
	return orbitals::OrbitalTypeSetOP( new orbitals::OrbitalTypeSet( directory ) );
}

/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
MMAtomTypeSetCOP
ChemicalManager::mm_atom_type_set( std::string const & tag )
{
	MMAtomTypeSets::const_iterator iter;
	{
#if defined MULTI_THREADED && CXX11
		utility::thread::ReadLockGuard lock( mmatomtype_mutex_ );
#endif
		iter = mm_atom_type_sets_.find( tag );
	}
	if ( iter == mm_atom_type_sets_.end() ) {

		// bind the relevant mm-atom-type-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< MMAtomTypeSet > () > func =
				boost::bind( &ChemicalManager::create_mm_atom_type_set, this, boost::cref(tag) );

#if defined MULTI_THREADED && defined CXX11
		iter = utility::thread::create_and_insert( func, mmatomtype_mutex_, tag, mm_atom_type_sets_ );
#else
		MMAtomTypeSetOP newset = func();
		iter = mm_atom_type_sets_.insert( std::make_pair( tag, newset ) ).first;
#endif
	}
	return iter->second;
}

/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
gasteiger::GasteigerAtomTypeSetCOP
ChemicalManager::gasteiger_atom_type_set( std::string const & tag /* = "default" */ )
{
	GasteigerAtomTypeSets::const_iterator iter( gasteiger_atom_type_sets_.find( tag ) );
	if ( iter == gasteiger_atom_type_sets_.end() ) {
		// read from file
		std::string filename( basic::database::full_name( "chemical/gasteiger/"+tag+"/atom_type_data.txt" ));
		ElementSetCAP elements( element_set(tag) );
		gasteiger::GasteigerAtomTypeSetOP new_set( new gasteiger::GasteigerAtomTypeSet( elements ) );
		new_set->read_file( filename );
		iter = gasteiger_atom_type_sets_.insert( std::make_pair( tag, new_set ) ).first;
		filename =  basic::database::full_name( "chemical/gasteiger/"+tag+"/bond_lengths.txt" );
		new_set->read_bond_file(filename);
	}
	return iter->second;
}
MMAtomTypeSetOP
ChemicalManager::create_mm_atom_type_set( std::string const & tag ) const
{
	// read from file
	std::string const filename( basic::database::full_name( "chemical/mm_atom_type_sets/"+tag+"/mm_atom_properties.txt" ));
	MMAtomTypeSetOP new_set( new MMAtomTypeSet() );
	new_set->read_file( filename );
	return new_set;
}



///@ details if the tag is not in the map, input it from a database file and add it
///to the map for future look-up.
ResidueTypeSetCOP
ChemicalManager::residue_type_set( std::string tag )
{

	using namespace basic;

	//There are a few experimental RNA residue type sets under
	//consideration.  The long term plan is to fold these into
	//fa_standard.  In the short term, override in ChemicalManager to
	//detect their use (via cmdline flag) and override base RNA to
	//something fancier.  Since these are covered by an override, they
	//are NOT covered by proper global-data-string tags RNA_PHENIX and
	//RNA_PROT_ERRASER.  SML & FC
	if ( tag == "rna" ) {
		using basic::options::OptionKeys::rna::corrected_geo;
		using basic::options::OptionKeys::rna::rna_prot_erraser;
		bool const use_corrected_rna_geo = basic::options::option[ corrected_geo ].value();
		bool const use_RNA_and_protein = basic::options::option[ rna_prot_erraser ].value();
		if ( use_corrected_rna_geo && !use_RNA_and_protein ) {
			tag = "rna_phenix";
		} else if ( use_corrected_rna_geo && use_RNA_and_protein ) {
			tag = "rna_prot_erraser";
		} else if ( !use_corrected_rna_geo && use_RNA_and_protein ) {
			utility_exit_with_message("cannot use -rna:rna_prot_erraser without -rna:corrected_geo");
		} //final case - if both are false - leave tag as-is
	} //if ( tag == "rna" ) {

	ResidueTypeSets::const_iterator iter;
	{ // scope for the ReadLockGuard
#if defined MULTI_THREADED && defined CXX11
		utility::thread::ReadLockGuard lock( restype_mutex_ );
#endif
		iter = residue_type_sets_.find( tag );
	}
	if ( iter == residue_type_sets_.end() ) {

		// bind the relevant residue-type-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< ResidueTypeSet > () > func =
				boost::bind( &ChemicalManager::create_residue_type_set, this, boost::cref(tag) );

#if defined MULTI_THREADED && defined CXX11
		iter = utility::thread::create_and_insert( func, restype_mutex_, tag, residue_type_sets_ );
#else
		ResidueTypeSetOP newset = func();
		iter = residue_type_sets_.insert( std::make_pair( tag, newset ) ).first;
#endif
	}
	return iter->second;
}

ResidueTypeSetOP
ChemicalManager::create_residue_type_set( std::string const & tag ) const {
	// Look for additional residue .params files specified on the cmd line
	std::vector< std::string > extra_params_files;
	std::vector< std::string > extra_patch_files;
	std::vector<core::chemical::ResidueTypeOP> extra_residues;

	if(tag == FA_STANDARD) {

		//this whole thing is desperately in need of some method extraction -- holy cow it does!
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
			TR.Debug<< std::endl;
			for(size_t j=1; j<= files.size(); j++){
				if (files[j].find("param")!=std::string::npos){
					TR << files[j]<< ", ";
					std::string path= directory+'/'+files[j];
					extra_params_files.push_back(path);
				}
			}
			TR.Debug<< std::endl;
		}

		utility::options::PathVectorOption & pvec_batch
		= basic::options::option[basic::options::OptionKeys::in::file::extra_res_batch_path];
		for(Size i=1, e= pvec_batch.size(); i<=e; i++){
			utility::vector1<std::string> subdirs;
			std::string directory=pvec_batch[i].name();

			utility::file::list_dir(directory, subdirs);
			TR.Debug<< std::endl;
			for ( size_t j=1; j<= subdirs.size();++j) {
				if ( subdirs[j] == "." || subdirs[j] == "..") {
					continue;
				}
				utility::vector1<std::string> files;
				utility::file::list_dir(directory+"/"+subdirs[j],files);
				for(size_t k=1; k<= files.size(); k++){
					if (files[k].find("param")!=std::string::npos){
						TR.Debug << files[k]<< ", ";
						std::string path= directory+'/'+subdirs[j]+'/'+files[k];
						extra_params_files.push_back(path);
					}
				}
			}

		}

		utility::options::FileVectorOption & molfilevec
		= basic::options::option[basic::options::OptionKeys::in::file::extra_res_mol];

		// this function itself does not (directly) modify any member data of class ChemicalManager,
		// but it is allowed to (indirectly) modify the singleton instance (which, to be fair,
		// is this instance) through singleton accessor functions.  In particular,
		// it will lead to the creation of the following *Sets for fa_standard.
		core::chemical::AtomTypeSetCOP atom_types = ChemicalManager::get_instance()->atom_type_set("fa_standard");
		core::chemical::ElementSetCOP elements = ChemicalManager::get_instance()->element_set("default");
		core::chemical::MMAtomTypeSetCOP mm_atom_types = ChemicalManager::get_instance()->mm_atom_type_set("fa_standard");
		core::chemical::orbitals::OrbitalTypeSetCOP orbital_types = ChemicalManager::get_instance()->orbital_type_set("fa_standard");

		sdf::MolFileIOReader molfile_reader;
		for(Size i=1, e = molfilevec.size(); i <= e;++i)
		{
			utility::file::FileName filename = molfilevec[i];
			utility::vector1< sdf::MolFileIOMoleculeOP > data( molfile_reader.parse_file( filename ) );
			utility::vector1< ResidueTypeOP > rtvec( sdf::convert_to_ResidueType( data ) );
			if( rtvec.size() > 1 ) {
				TR.Warning << "molfile " << filename << " has more than one model -- behavior towards this file may change in the future." << std::endl;
			}
			for( core::Size ii(1); ii <= rtvec.size(); ++ii ) {
				extra_residues.push_back( rtvec[ii] );
			}
		}

		if(basic::options::option[basic::options::OptionKeys::in::file::extra_res_database].user())
		{
			utility::sql_database::DatabaseMode::e database_mode(
					utility::sql_database::database_mode_from_name(
							basic::options::option[basic::options::OptionKeys::in::file::extra_res_database_mode]));
			std::string database_name(basic::options::option[basic::options::OptionKeys::in::file::extra_res_database]);
			std::string database_pq_schema(basic::options::option[basic::options::OptionKeys::in::file::extra_res_pq_schema]);


			utility::sql_database::sessionOP db_session(
					basic::database::get_db_session(database_mode, database_name, database_pq_schema));

			ResidueDatabaseIO residue_database_interface;

			if(basic::options::option[basic::options::OptionKeys::in::file::extra_res_database_resname_list].user())
			{
				utility::file::FileName residue_list = basic::options::option[basic::options::OptionKeys::in::file::extra_res_database_resname_list];
				utility::io::izstream residue_name_file(residue_list);
				std::string residue_name;
				while(residue_name_file >> residue_name)
				{
					//residue_name_file >> residue_name;
					TR <<residue_name <<std::endl;
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

		// Patches
		utility::options::FileVectorOption & pfvec
		= basic::options::option[ basic::options::OptionKeys::in::file::extra_patch_fa ];
		for(Size i = 1, e = pfvec.size(); i <= e; ++i) {
			extra_patch_files.push_back( pfvec[i].name());
		}

	} else if(tag == CENTROID) {
		utility::options::FileVectorOption & fvec
		= basic::options::option[ basic::options::OptionKeys::in::file::extra_res_cen ];
		for(Size i = 1, e = fvec.size(); i <= e; ++i) {
			utility::file::FileName fname = fvec[i];
			extra_params_files.push_back(fname.name());
		}
		// Patches
		utility::options::FileVectorOption & pfvec
		= basic::options::option[ basic::options::OptionKeys::in::file::extra_patch_cen ];
		for(Size i = 1, e = pfvec.size(); i <= e; ++i) {
			extra_patch_files.push_back( pfvec[i].name());
		}
	}

	// generically specify extra res (not necessarily part of fa_standard) -- will get added to
	//  any and every residue_type_set instantiated.
	utility::options::FileVectorOption & fvec
	= basic::options::option[ basic::options::OptionKeys::in::file::extra_res ];
	for(Size i = 1, e = fvec.size(); i <= e; ++i) {
		utility::file::FileName fname = fvec[i];
		extra_params_files.push_back(fname.name());
	}

	// read from file
	TR.Debug << "CHEMICAL_MANAGER: read residue types: " << tag << std::endl;
	// redirecting to new icoor folder
	std::string temp_str( basic::database::full_name( "chemical/residue_type_sets/"+tag ) );
	if ( tag == FA_STANDARD ) {
		if ( basic::options::option[basic::options::OptionKeys::corrections::chemical::icoor_05_2009]) {
			temp_str += "_05.2009_icoor";
		} else if ( basic::options::option[basic::options::OptionKeys::mistakes::chemical::pre_talaris2013_geometries ]) {
			temp_str += "_pre_talaris2013";
		}
	} else if ( tag == CENTROID ) {
		if ( basic::options::option[basic::options::OptionKeys::mistakes::chemical::pre_talaris2013_geometries ]) {
			temp_str += "_pre_talaris2013";
		}
	}
	temp_str += "/";


	std::string const directory( temp_str );
	ResidueTypeSetOP new_set( new ResidueTypeSet( tag, directory ) );
	new_set->init( extra_params_files, extra_patch_files );

	for(core::Size index =0 ;index < extra_residues.size();++index)
	{
		//TR << extra_residues[index]->name3() <<std::endl;
		new_set->add_residue_type(extra_residues[index]);
		extra_residues[index]->residue_type_set(ResidueTypeSetCAP(new_set));
	}

	return new_set;
}

///@ details if the tag is not in the map, input it from a database file and add it
///to the map for future look-up.
/// THIS FUNCTION IS DECIDEDLY NOT THREADSAFE!
ResidueTypeSet &
ChemicalManager::nonconst_residue_type_set( std::string const & tag )
{
	// trigger initialization if necessary:
	residue_type_set( tag );

	return *( residue_type_sets_.find( tag )->second );
}

ResidueTypeSetOP
ChemicalManager::nonconst_residue_type_set_op( std::string const & tag )
{
	// trigger initialization if necessary:
	residue_type_set( tag );

	return residue_type_sets_.find( tag )->second;
}

/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
IdealBondLengthSetCOP
ChemicalManager::ideal_bond_length_set( std::string const & tag )
{
	IdealBondLengthSets::const_iterator iter;
	{ // scope for the ReadLockGuard
#if defined MULTI_THREADED && defined CXX11
		utility::thread::ReadLockGuard lock( idealbondlength_mutex_ );
#endif
		iter = ideal_bond_length_sets_.find( tag );
	}
	if ( iter == ideal_bond_length_sets_.end() ) {

		// bind the relevant ideal-bond-length-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< IdealBondLengthSet > () > func =
				boost::bind( &ChemicalManager::create_ideal_bond_length_set, this, boost::cref(tag) );

#if defined MULTI_THREADED && defined CXX11
		iter = utility::thread::create_and_insert( func, idealbondlength_mutex_, tag, ideal_bond_length_sets_ );
#else
		IdealBondLengthSetOP newset = func();
		iter = ideal_bond_length_sets_.insert( std::make_pair( tag, newset ) ).first;
#endif
	}
	return iter->second;
}

IdealBondLengthSetOP
ChemicalManager::create_ideal_bond_length_set( std::string const & tag ) const
{
	// read from file
	std::string const filename( basic::database::full_name( "chemical/atom_type_sets/"+tag+"/ideal_bond_lengths.txt" ));
	IdealBondLengthSetOP new_set( new IdealBondLengthSet() );
	new_set->read_file( filename );
	return new_set;
}

// global data
/// @brief tag name for querying fullatom chemical type set.
std::string const FA_STANDARD( "fa_standard" );
/// @brief tag name for querying centroid chemical type set.
std::string const CENTROID( "centroid" );
/// @brief tag name for querying centroid_rot chemical type set.
std::string const CENTROID_ROT( "centroid_rot" );
/// @brief tag name for querying coarse-grained chemical type set.
std::string const COARSE_TWO_BEAD( "coarse_two_bead" );
/// @brief tag name for querying hybrid fullatom+centroid chemical type set.
std::string const HYBRID_FA_STANDARD_CENTROID( "hybrid_fa_standard_centroid" );
/// @brief tag name for querying RNA chemical type set.
std::string const FA_RNA = "rna";
/// @brief tag name for querying COARSE_RNA chemical type set.
std::string const COARSE_RNA( "coarse_rna" );

} // namespace core
} // namespace chemical
