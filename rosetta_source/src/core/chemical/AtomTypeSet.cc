// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin AtomTypeSet
///
/// @brief
/// A class for reading in the atom type properties
///
/// @detailed
/// This class reads in the atom_properties.txt file which contains the "chemical" information for atoms.
/// This does not contain the actual properties, but sets the properties through the AtomType class.
/// This class is called by the ChemicalManager
///
///
///
/// @authors
/// Phil Bradley
/// Steven Combs - comments
///
///
/// @last_modified December 6 2010
/////////////////////////////////////////////////////////////////////////
// Unit headers
#include <core/chemical/AtomTypeSet.hh>

// Project headers
#include <basic/Tracer.hh>

// Utility headers
// Commented by inclean daemon #include <utility/exit.hh>

// C++ headers
// Commented by inclean daemon #include <string>
// Commented by inclean daemon #include <sstream>
#include <fstream>
#include <iostream>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <utility/io/izstream.hh>

namespace core {
namespace chemical {

static basic::Tracer tr("core.chemical");


////////////////////////////////////////////////////////////////////////////////

AtomTypeSet::AtomTypeSet( std::string const & directory )
{
	directory_ = directory;

	read_file( directory + "/atom_properties.txt" );

	utility::io::izstream data( ( directory+"/extras.txt" ).c_str() );
	if ( data.good() ) { // add extra data
		std::string line;
		while( getline( data, line ) ) {
			if ( line.size() && line[0] == '#' ) continue;
			add_parameters_from_file( directory+"/"+line );
		}
	}
	data.close();

}



/// @brief file I/O
///
/// @details initialize an AtomTypeSet from an external file "filename",
/// and set parameters and properties for each AtomType.
/// Refer to minirosetta_database_stock/chemical/atom_type_sets/fa_standard/atom_properties.txt
/// for file format
///
void
AtomTypeSet::read_file( std::string const & filename )
{
	utility::io::izstream data( filename.c_str() ); 

	if ( !data.good() ) utility_exit_with_message( "Unable to open atomset file: "+filename );

	// parse the header line
	utility::vector1< std::string > tags;
	{ // scope
		std::string line, tag, tag2;
		getline( data, line );
		std::istringstream l( line );
		l >> tag >> tag2;
		if ( tag != "NAME" || tag2 != "ATOM" ) {
			utility_exit_with_message("AtomTypeSet::read_file: bad first line: "+	line );
		}
		l >> tag;
		while ( !l.fail() ) {
			tags.push_back( tag );
			l >> tag;
		}
	}

	// now parse the rest of the file
	Size const ntags( tags.size() );
	{
		using namespace basic;

		std::string line, tag, name_wo_whitespace;
		while ( getline( data,line ) ) {
			std::istringstream l( line );
			l >> name_wo_whitespace;
			if ( l.fail() || name_wo_whitespace.find("#",0) == 0 ) continue; // skip comment,blank lines
			l >> tag;
			if ( l.fail() || tag.size() < 1 ) {
				utility_exit_with_message("bad line: "+line);
			}

			//			std::string const name( line.substr(0,4) );
			std::string const element( tag );
			AtomType* atom_type_ptr( new AtomType( name_wo_whitespace, element ) );

			// now parse the parameters
			for ( Size i=1; i<= ntags; ++i ) {
				Real setting;
				l >> setting;
				atom_type_ptr->set_parameter( tags[i], setting );
			}
			if ( l.fail() ) {
				utility_exit_with_message("bad line: "+line);
			}

			// now parse the properties
			l >> tag;
			while ( !l.fail() && tag.find("#",0) != 0) {
				atom_type_ptr->set_property( tag, true );
				l >> tag;
			}

			// add this to the list
			atoms_.push_back( atom_type_ptr );
			//		atom_type_index_[ name ] = atoms_.size();
			if ( atom_type_index_.count( name_wo_whitespace ) ) {
				utility_exit_with_message("AtomTypeSet:: duplicate atom name "+name_wo_whitespace);
			}
			atom_type_index_[ name_wo_whitespace ] = atoms_.size();
			tr.Debug << "New atom type: " << name_wo_whitespace << ' ' << element << std::endl; //std::endl;
		}
	} // scope


}


///////////////////////////////////////////////////////////////////////////////
/// @details  Private helper function for filling in default values in the fxn add_parameters_from_file
/// Enables the user to specify a default parameter set to be used and then provide a few modifications.
///
/// eg in the dna_interface lj-radii parameter set we shrink the hydrogens but leave the rest unchanged.
///
/// @note This function is very SLOW, but that should be OK we only use it a bit right at the start.
Real
AtomTypeSet::get_default_parameter( std::string const & param_name, std::string const & atm_name ) const
{
	AtomType const & atom_type( *( atoms_ [ atom_type_index( atm_name ) ] ) );
	if ( has_extra_parameter( param_name ) ) {
		return atom_type.extra_parameter( extra_parameter_index( param_name ) );
	} else {
		// get hardcoded params from atomtype

		if ( param_name == "LJ_RADIUS" ) {
			return atom_type.lj_radius();

		} /// etc etc add LJ_WDEPTH, LK_***

	}
	utility_exit_with_message( "unrecognized parameter type: "+param_name );
	return 0.0; // appease compiler
}


///////////////////////////////////////////////////////////////////////////////
void
AtomTypeSet::add_parameters_from_file( std::string const & filename )
{

	// parse the header line
	utility::vector1< std::string > tags;

	Size const index_offset( extra_parameter_indices_.size() );

	utility::vector1< std::string > default_parameter_names;
	utility::vector1< std::string > lines;
	{ // read all the lines from the file
		utility::io::izstream data( filename.c_str() ); 

		if ( !data.good() ) utility_exit_with_message( "Unable to open atomset parameter file: "+filename );
		std::string line, tag;
		while ( getline( data,line ) ) {
			std::istringstream l( line );
			l >> tag;
			if ( l.fail() || tag.size() < 1 || tag[0] == '#' ) continue; // skip blank lines or comments

			if ( tag == "NAME" ) {
				l >> tag;
				while ( !l.fail() ) {
					if ( tag[0] == '#' ) break;
					tags.push_back( tag );
					extra_parameter_indices_[ tag ] = tags.size() + index_offset;
					l >> tag;
				}
			} else if ( tag == "DEFAULT" ) {
				l >> tag;
				while ( !l.fail() ) {
					if ( tag[0] == '#' ) break;
					default_parameter_names.push_back( tag );
					l >> tag;
				}
			} else {
				lines.push_back( line );
			}
		}
		data.close();
	}

	if ( tags.empty() ) utility_exit_with_message("AtomTypeSet::read_file: missing NAME line");
	if ( !default_parameter_names.empty() && default_parameter_names.size() != tags.size() ) {
		std::cout << "AtomTypeSet:: number of params doesnt match number of defaults " <<
			default_parameter_names.size() << ' ' << tags.size() << std::endl;
		utility_exit();
	}

	// now parse the rest of the file
	Size const ntags( tags.size() );
	std::map< std::string, utility::vector1< Real > > all_parameters;
	{
		std::string tag, name_wo_whitespace;
		for ( Size ii=1; ii<= lines.size(); ++ii ) {
			std::string const & line( lines[ii] );
			std::istringstream l( line );
			l >> tag; // name_wo_whitespace
			if ( tag.find("#",0) == 0 ) continue; // skip comment lines

			//			std::string const name( line.substr(0,4) );

			// now parse the parameters
			utility::vector1< Real > parameters;
			for ( Size i=1; i<= ntags; ++i ) {
				Real setting;
				l >> setting;
				parameters.push_back( setting );
			}
			if ( l.fail() || parameters.size() != tags.size() ) {
				utility_exit_with_message("bad line: "+line);
			}

			all_parameters[ tag ] = parameters;
		} // loop over lines from file
	}

	// now fill in the data
	for ( std::map< std::string, int >::const_iterator
					iter = atom_type_index_.begin(), iter_end = atom_type_index_.end(); iter != iter_end; ++iter ) {
		std::string const & name( iter->first );
		int const atom_index( iter->second );

		//if ( name.size() < 4 ) continue; // ignore stripped ws versions of the atom names
		std::map< std::string, utility::vector1< Real > >::const_iterator iter2( all_parameters.find( name ) );
		utility::vector1< Real > params;

		if ( iter2 == all_parameters.end() ) {
			if ( !default_parameter_names.empty() ) {
				for ( Size i=1; i<= default_parameter_names.size(); ++i ) {
					Real const default_param( get_default_parameter( default_parameter_names[i], name ) );
					basic::T("core.chemical.AtomTypeSet") << "Using default parameter " << default_parameter_names[i] << " = " <<
						default_param << " in place of " << tags[i] << " for atomtype " << name << '\n';
					params.push_back( default_param );
				}
			} else {
				iter2 = all_parameters.find( "****" );
				if ( iter2 == all_parameters.end() ) {
					utility_exit_with_message( "no parameters specified for atom type: "+name+" in file "+filename );
				}
				params = iter2->second;
			}
		} else {
			params = iter2->second;
      //pbadebug
      //std::cout << "params " << params << std::endl;
		}
		//utility::vector1< Real > const & params( iter2->second );
		assert( params.size() == tags.size() );
		for ( Size i=1; i<= tags.size(); ++i ) {
			atoms_[ atom_index ]->set_extra_parameter( i + index_offset, params[i] );
		}
	} // loop over atom names in this AtomTypeSet

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


} // chemical
} // core
