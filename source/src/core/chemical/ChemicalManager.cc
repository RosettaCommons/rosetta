// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// Chemical manager class
///
/// @details
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
/// @author
/// Andrew Leaver-Fay (leaverfa@email.unc.edu)
/// Steven Combs - comments
///
///
/////////////////////////////////////////////////////////////////////////


#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/IdealBondLengthSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/GlobalResidueTypeSet.hh>
#include <core/chemical/sdf/MolFileIOReader.hh>
#include <core/chemical/mmCIF/mmCIFParser.hh>
//#include <core/chemical/sdf/MolFileIOData.hh>
#include <core/chemical/util.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
// Project headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <utility/file/file_sys_util.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <core/chemical/ResidueType.hh>

#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
//#include <boost/function.hpp>


namespace core {
namespace chemical {

static THREAD_LOCAL basic::Tracer TR( "core.chemical.ChemicalManager" );

/// @brief private constructor to guarantee the singleton
ChemicalManager::ChemicalManager(){}

/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
AtomTypeSetCOP
ChemicalManager::atom_type_set( std::string const & tag )
{
	AtomTypeSets::const_iterator iter;
	{ // scope for the read lock guard
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( atomtype_mutex_ );
#endif
		iter = atom_type_sets_.find( tag );
	}
	if ( iter == atom_type_sets_.end() ) {

		// bind the relevant atom-type-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< AtomTypeSet > () > func =
			boost::bind( &ChemicalManager::create_atom_type_set, this, boost::cref(tag) );

#if defined MULTI_THREADED
		iter = utility::thread::create_and_insert( func, atomtype_mutex_, tag, atom_type_sets_ );
#else
		AtomTypeSetOP newset = func();
		iter = atom_type_sets_.insert( std::make_pair( tag, newset ) ).first;
#endif
	}
	return iter->second;
}

/// @details Actually go and create an AtomTypeSet
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
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( elem_mutex_ );
#endif
		iter = element_sets_.find( tag );
	}
	if ( iter == element_sets_.end() ) {

		// bind the relevant atom-type-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< ElementSet > () > func =
			boost::bind( &ChemicalManager::create_element_set, this, boost::cref(tag) );

#if defined MULTI_THREADED
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
	ElementSetOP new_set( new ElementSet( tag ) );
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
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( orbtype_mutex_ );
#endif
		iter = orbital_type_sets_.find( tag );
	}
	if ( iter == orbital_type_sets_.end() ) {

		// bind the relevant atom-type-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< orbitals::OrbitalTypeSet > () > func =
			boost::bind( &ChemicalManager::create_orbital_type_set, this, boost::cref(tag) );

#if defined MULTI_THREADED
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
	return orbitals::OrbitalTypeSetOP( new orbitals::OrbitalTypeSet( directory, tag ) );
}

/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
MMAtomTypeSetCOP
ChemicalManager::mm_atom_type_set( std::string const & tag )
{
	MMAtomTypeSets::const_iterator iter;
	{
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( mmatomtype_mutex_ );
#endif
		iter = mm_atom_type_sets_.find( tag );
	}
	if ( iter == mm_atom_type_sets_.end() ) {

		// bind the relevant mm-atom-type-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< MMAtomTypeSet > () > func =
			boost::bind( &ChemicalManager::create_mm_atom_type_set, this, boost::cref(tag) );

#if defined MULTI_THREADED
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
		gasteiger::GasteigerAtomTypeSetOP new_set( new gasteiger::GasteigerAtomTypeSet( elements, tag ) );
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
	MMAtomTypeSetOP new_set( new MMAtomTypeSet( tag ) );
	new_set->read_file( filename );
	return new_set;
}


/// @brief query residue_type_set by a type
ResidueTypeSetCOP
ChemicalManager::residue_type_set( TypeSetMode type_set_mode ) {
	std::string standard_name( string_from_type_set_mode( type_set_mode ) );
	return residue_type_set( standard_name );
}

/// @ details if the tag is not in the map, input it from a database file and add it
///to the map for future look-up.
ResidueTypeSetCOP
ChemicalManager::residue_type_set( std::string const & tag )
{

	using namespace basic;

	// no longer supporting legacy (uncorrected) rna bond params -- probably should deprecate this option entirely.
	runtime_assert ( basic::options::option[ basic::options::OptionKeys::rna::corrected_geo ]() );

	GlobalResidueTypeSets::const_iterator iter;
	{ // scope for the ReadLockGuard
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( restype_mutex_ );
#endif
		iter = residue_type_sets_.find( tag );
	}
	if ( iter == residue_type_sets_.end() ) {

		// bind the relevant residue-type-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< GlobalResidueTypeSet > () > func =
			boost::bind( &ChemicalManager::create_residue_type_set, this, boost::cref(tag) );

#if defined MULTI_THREADED
		iter = utility::thread::create_and_insert( func, restype_mutex_, tag, residue_type_sets_ );
#else
		GlobalResidueTypeSetOP newset = func();
		iter = residue_type_sets_.insert( std::make_pair( tag, newset ) ).first;
#endif
	}
	return iter->second;
}

bool
ChemicalManager::has_residue_type_set( std::string const & tag ) {
#if defined MULTI_THREADED && defined CXX11
	utility::thread::ReadLockGuard lock( restype_mutex_ );
#endif
	return( residue_type_sets_.find( tag )  != residue_type_sets_.end() );
}

GlobalResidueTypeSetOP
ChemicalManager::create_residue_type_set( std::string const & tag ) const {

	// read from file
	TR.Debug << "CHEMICAL_MANAGER: read residue types: " << tag << std::endl;
	// redirecting to new icoor folder
	std::string temp_str( basic::database::full_name( "chemical/residue_type_sets/"+tag ) );
	if ( tag == FA_STANDARD ) {
		if ( basic::options::option[basic::options::OptionKeys::corrections::chemical::icoor_05_2009] ) {
			temp_str += "_05.2009_icoor";
		} else if ( basic::options::option[basic::options::OptionKeys::mistakes::chemical::pre_talaris2013_geometries ] ) {
			temp_str += "_pre_talaris2013";
		}
	} else if ( tag == CENTROID ) {
		if ( basic::options::option[basic::options::OptionKeys::mistakes::chemical::pre_talaris2013_geometries ] ) {
			temp_str += "_pre_talaris2013";
		}
	}
	temp_str += "/";


	std::string const directory( temp_str );

	// Note on thread safety: The GlobalResidueTypeSet constructor can (and will) call back into the ChemicalManager
	// to initialize things like AtomTypeSets, etc. So long as the RTS/ATS/etc. are protected under different mutexes,
	// there shouldn't be any issues with a deadlock.
	GlobalResidueTypeSetOP new_set( new GlobalResidueTypeSet( tag, directory ) );
	return new_set;
}

/// @details if the tag is not in the map, input it from a database file and add it
/// to the map for future look-up.
IdealBondLengthSetCOP
ChemicalManager::ideal_bond_length_set( std::string const & tag )
{
	IdealBondLengthSets::const_iterator iter;
	{ // scope for the ReadLockGuard
#if defined MULTI_THREADED
		utility::thread::ReadLockGuard lock( idealbondlength_mutex_ );
#endif
		iter = ideal_bond_length_sets_.find( tag );
	}
	if ( iter == ideal_bond_length_sets_.end() ) {

		// bind the relevant ideal-bond-length-set creation function and its arguments so that
		// it can be passed to the utility::thread::create_and_insert function.
		boost::function< utility::pointer::shared_ptr< IdealBondLengthSet > () > func =
			boost::bind( &ChemicalManager::create_ideal_bond_length_set, this, boost::cref(tag) );

#if defined MULTI_THREADED
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
/// @brief tag name for querying hybrid fullatom+centroid chemical type set.
std::string const HYBRID_FA_STANDARD_CENTROID( "hybrid_fa_standard_centroid" );
/// @brief tag name for querying COARSE_RNA chemical type set.
std::string const COARSE_RNA( "coarse_rna" );

TypeSetMode
type_set_mode_from_string( std::string const & mode ) {
	if ( mode == FA_STANDARD ) return FULL_ATOM_t;
	if ( mode == "full_atom" ) return FULL_ATOM_t;
	if ( mode == "default" ) return DEFAULT_t;
	if ( mode == CENTROID ) return CENTROID_t;
	if ( mode == CENTROID_ROT ) return CENTROID_ROT_t;
	if ( mode == HYBRID_FA_STANDARD_CENTROID ) return HYBRID_FA_STANDARD_CENTROID_t;
	if ( mode == COARSE_RNA ) return COARSE_RNA_t;
	utility_exit_with_message("String '"+mode+"' not recognized as a TypeSetMode.");
}

std::string
string_from_type_set_mode( TypeSetMode mode ) {
	switch ( mode ) {
	case FULL_ATOM_t :
		return FA_STANDARD;
	case DEFAULT_t :
		return "default";
	case CENTROID_t :
		return CENTROID;
	case CENTROID_ROT_t :
		return CENTROID_ROT;
	case HYBRID_FA_STANDARD_CENTROID_t :
		return HYBRID_FA_STANDARD_CENTROID;
	case COARSE_RNA_t :
		return COARSE_RNA;
	case INVALID_t :
		return "INVALID_CATEGORY";
	default :
		TR.Error << "Value " << mode << " is not a valid TypeSetMode." << std::endl;
		utility_exit_with_message("Can't convert TypeSetMode to string.");
	}
}

std::ostream &
operator <<( std::ostream & out, TypeSetMode mode ) {
	out << string_from_type_set_mode( mode );
	return out;
}

} // namespace core
} // namespace chemical
