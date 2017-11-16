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
/// A class for reading in the atom type properties
///
/// @details
/// This class reads in the atom_properties.txt file which contains the "chemical" information for atoms.
/// This does not contain the actual properties, but sets the properties through the AtomType class.
/// This class is called by the ChemicalManager
///
///
///
/// @author
/// Phil Bradley
/// Steven Combs - comments
///
/////////////////////////////////////////////////////////////////////////
// Unit headers
#include <core/chemical/AtomTypeSet.hh>

// Project headers
#include <basic/Tracer.hh>

#include <iostream>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/util.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <basic/database/open.hh>
#include <basic/database/sql_utils.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

using namespace basic::options;

namespace core {
namespace chemical {

static basic::Tracer tr( "core.chemical.AtomTypeSet" );


////////////////////////////////////////////////////////////////////////////////

AtomTypeSet::AtomTypeSet( std::string const & directory ):
	mode_( INVALID_t ),
	directory_( directory )
{

	read_meta_file( directory+"/meta.txt" );

	read_file( directory + "/atom_properties.txt" );

	utility::io::izstream data( ( directory+"/extras.txt" ).c_str() );
	if ( data.good() ) { // add extra data
		std::string line;
		while ( getline( data, line ) ) {
			if ( line.size() && line[0] == '#' ) continue;
			add_parameters_from_file( directory+"/"+line );
		}
	}
	data.close();

	legacy_command_line_post_processing();

	clone_atom_types_from_commandline();
}

AtomTypeSet::AtomTypeSet(
	std::string const & name,
	utility::sql_database::sessionOP db_session):
	mode_( INVALID_t )
{

	directory_ = basic::database::full_name( "chemical/atom_type_sets/" + name);

	// TODO: There probably should be a loading of the mode here, though where that would be set
	// in the database is a good question.

	{ // add atom type to atom type set
		std::string stmt_string =
			"SELECT name FROM atom_types WHERE atom_type_set_name = ?;";
		cppdb::statement stmt(
			basic::database::safely_prepare_statement(stmt_string, db_session));
		stmt.bind(1, name);
		cppdb::result res(basic::database::safely_read_from_database(stmt));

		std::string atom_type_name;
		while ( res.next() ) {
			res >> atom_type_name;
			AtomType & atom_type(
				create_atom_type_from_database(
				name, atom_type_name, db_session));
			read_atom_type_properties_table(
				name, atom_type, db_session);
			read_atom_type_extra_parameters_table(
				name, atom_type, db_session);
		}
	}

	{ // set the extra parameter indices
		std::string stmt_string =
			"SELECT DISTINCT\n"
			"\tparameter\n"
			"FROM\n"
			"\tatom_type_extra_parameters\n"
			"WHERE\n"
			"\tatom_type_set_name = ?\n"
			"ORDER BY\n"
			"\tparameter\n";
		cppdb::statement stmt(
			basic::database::safely_prepare_statement(stmt_string, db_session));
		stmt.bind(1, name);
		cppdb::result res(
			basic::database::safely_read_from_database(stmt));

		Size extra_parameter_index(1);
		std::string extra_parameter_name;
		while ( res.next() ) {
			res >> extra_parameter_name;
			extra_parameter_indices_[extra_parameter_name] = extra_parameter_index;
			++extra_parameter_index;
		}
	}
}


AtomTypeSet::~AtomTypeSet() {
	// The atoms in the atom type set are kept with raw pointers so they
	// must be deleted to prevent memory leaks
	for ( Size i=1; i <= atom_type_index_.size(); ++i ) {
		delete atoms_[i];
	}
}

/// @detail  The directory is like '$ROSETTA3_DB/rosetta_database/chemical/atom_type_sets/<atom_type_set_name>/'
/// Return 'atom_type_set_name'
/// Note: strip off the trailing slash, if it exists
std::string
AtomTypeSet::name() const {
	Size const last_char_pos(directory_.find_last_not_of('/'));
	if ( last_char_pos == std::string::npos || last_char_pos == 0 ) return directory_;

	Size first_char_pos(directory_.find_last_of('/', last_char_pos - 1) + 1);
	if ( first_char_pos == std::string::npos ) first_char_pos = 0;

	return directory_.substr(first_char_pos, last_char_pos - first_char_pos + 1);
}

/// @brief lookup the atom_type by the atom_type_name string
int
AtomTypeSet::atom_type_index( std::string const & atom_type_name ) const
{
	auto iter( atom_type_index_.find( atom_type_name ) );
	if ( iter == atom_type_index_.end() ) {
		std::string trimmed( utility::trim(atom_type_name) ); // Try with stripped whitespace
		iter = atom_type_index_.find( trimmed );
		if ( iter == atom_type_index_.end() ) {
			utility_exit_with_message("unrecognized atom_type_name '"+atom_type_name+"'");
		}
	}
	return iter->second;
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
			utility_exit_with_message("AtomTypeSet::read_file: bad first line: "+ line );
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

			//   std::string const name( line.substr(0,4) );
			std::string const element( tag );
			auto* atom_type_ptr( new AtomType( name_wo_whitespace, element ) );

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
			// AMW: fixing up cppcheck errors
			// string::find here is unnecessary because you are just checking the first character
			while ( !l.fail() && tag.length() != 0 && tag[0] != '#' ) {//tag.find("#",0) != 0 ) {
				atom_type_ptr->set_property( tag, true );
				l >> tag;
			}

			// add this to the list
			atoms_.push_back( atom_type_ptr );
			//  atom_type_index_[ name ] = atoms_.size();
			if ( atom_type_index_.count( name_wo_whitespace ) ) {
				utility_exit_with_message("AtomTypeSet:: duplicate atom name "+name_wo_whitespace);
			}
			atom_type_index_[ name_wo_whitespace ] = atoms_.size();
			tr.Debug << "New atom type: " << name_wo_whitespace << ' ' << element << std::endl; //std::endl;
		}
	} // scope


}

/// @brief Read in meta information from the given file
void
AtomTypeSet::read_meta_file( std::string const & filename ) {

	utility::io::izstream data( filename.c_str() );
	if ( !data.good() ) {
		return; // Be silent if the directory doesn't have a meta file.
	}

	std::string line, tag;
	while ( getline( data, line ) ) {
		// Skip empty lines and comments.
		if ( line.size() < 1 || line[0] == '#' ) continue;
		std::stringstream l( line );
		l >> tag;
		if ( tag == "TYPE_SET_MODE" || tag == "TYPE_SET_CATEGORY" ) {
			l >> tag;
			mode_ = type_set_mode_from_string( tag );
		}
	}

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
		} else if ( param_name == "LJ_WDEPTH" ) {
			return atom_type.lj_wdepth();
		} else if ( param_name == "LK_VOLUME" ) {
			return atom_type.lk_volume();
		} else if ( param_name == "LK_DGFREE" ) {
			return atom_type.lk_dgfree();
		} else if ( param_name == "LK_LAMBDA" ) {
			return atom_type.lk_lambda();
		}
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
		std::cerr << "AtomTypeSet:: number of params doesnt match number of defaults " <<
			default_parameter_names.size() << ' ' << tags.size() << std::endl;
		utility_exit();
	}

	// now parse the rest of the file
	Size const ntags( tags.size() );
	std::map< std::string, utility::vector1< Real > > all_parameters;
	{
		std::string tag;//, name_wo_whitespace;
		for ( Size ii=1; ii<= lines.size(); ++ii ) {
			std::string const & line( lines[ii] );
			std::istringstream l( line );
			l >> tag; // name_wo_whitespace
			if ( tag.find("#",0) == 0 ) continue; // skip comment lines

			//   std::string const name( line.substr(0,4) );

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

		std::string paramsrc("UNK");
		if ( iter2 == all_parameters.end() ) {
			if ( !default_parameter_names.empty() ) {
				paramsrc = "from_defaults";
				for ( Size i=1; i<= default_parameter_names.size(); ++i ) {
					Real const default_param( get_default_parameter( default_parameter_names[i], name ) );
					// this output doesnt seem to get written out!? so I added more verbose output below (PB)
					tr << "Using default parameter " << default_parameter_names[i] << " = " <<
						default_param << " in place of " << tags[i] << " for atomtype " << name << '\n';
					params.push_back( default_param );
				}
			} else {
				paramsrc = "from_****";
				iter2 = all_parameters.find( "****" );
				if ( iter2 == all_parameters.end() ) {
					utility_exit_with_message( "no parameters specified for atom type: "+name+" in file "+filename );
				}
				params = iter2->second;
			}
		} else {
			paramsrc = "from_file";
			params = iter2->second;
			//pbadebug
			//std::cout << "params " << params << std::endl;
		}
		//utility::vector1< Real > const & params( iter2->second );
		runtime_assert( params.size() == tags.size() );
		for ( Size i=1; i<= tags.size(); ++i ) {
			tr.Trace << "setting extra parameter: " << tags[i] << ' ' <<atoms_[ atom_index ]->name() << ' ' <<
				paramsrc << ' ' << params[i] << std::endl;
			atoms_[ atom_index ]->set_extra_parameter( i + index_offset, params[i] );
		}
	} // loop over atom names in this AtomTypeSet

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


AtomType &
AtomTypeSet::create_atom_type_from_database(
	std::string const & atom_type_set_name,
	std::string const & atom_type_name,
	utility::sql_database::sessionOP db_session
) {
	std::string stmt_string =
		"SELECT\n"
		"\telement,\n"
		"\tlennard_jones_radius REAL,\n"
		"\tlennard_jones_well_depth REAL,\n"
		"\tlazaridis_karplus_lambda REAL,\n"
		"\tlazaridis_karplus_degrees_of_freedom REAL,\n"
		"\tlazaridis_karplus_volume REAL\n"
		"FROM\n"
		"\tatom_types\n"
		"WHERE\n"
		"\tatom_type_set_name = ? AND name = ?;";

	cppdb::statement stmt(
		basic::database::safely_prepare_statement(stmt_string, db_session));
	stmt.bind(1, atom_type_set_name);
	stmt.bind(2, atom_type_name);
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	if ( !res.next() ) {
		utility_exit_with_message(
			"could not find atom '" + atom_type_name + "' in '" +
			atom_type_set_name + "'.");
	}

	std::string element;
	Real lennard_jones_radius;
	Real lennard_jones_well_depth;
	Real lazaridis_karplus_lambda;
	Real lazaridis_karplus_degrees_of_freedom;
	Real lazaridis_karplus_volume;

	res
		>> element
		>> lennard_jones_radius
		>> lennard_jones_well_depth
		>> lazaridis_karplus_lambda
		>> lazaridis_karplus_degrees_of_freedom
		>> lazaridis_karplus_volume;

	auto * atom_type_ptr(new AtomType(atom_type_name, element));
	atom_type_ptr->set_parameter("LJ_RADIUS", lennard_jones_radius);
	atom_type_ptr->set_parameter("LJ_WDEPTH", lennard_jones_well_depth);
	atom_type_ptr->set_parameter("LK_LAMBDA", lazaridis_karplus_lambda);
	atom_type_ptr->set_parameter("LK_DGFREE", lazaridis_karplus_degrees_of_freedom);
	atom_type_ptr->set_parameter("LK_VOLUME", lazaridis_karplus_volume);

	atoms_.push_back(atom_type_ptr);
	atom_type_index_[atom_type_ptr->name()] = atoms_.size();
	return *atom_type_ptr;
}


void
AtomTypeSet::read_atom_type_properties_table(
	std::string const & atom_type_set_name,
	AtomType & atom_type,
	utility::sql_database::sessionOP db_session
) {
	std::string stmt_string =
		"SELECT\n"
		"\tproperty\n"
		"FROM\n"
		"\tatom_type_properties\n"
		"WHERE\n"
		"\tatom_type_set_name = ? AND name = ?;";

	cppdb::statement stmt(basic::database::safely_prepare_statement(stmt_string, db_session));
	stmt.bind(1, atom_type_set_name);
	stmt.bind(2, atom_type.name());
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	std::string property;
	while ( res.next() ) {
		res >> property;
		atom_type.add_property(property);
	}
}

void
AtomTypeSet::read_atom_type_extra_parameters_table(
	std::string const & atom_type_set_name,
	AtomType & atom_type,
	utility::sql_database::sessionOP db_session
) {
	std::string stmt_string =
		"SELECT\n"
		"\tvalue\n"
		"FROM\n"
		"\tatom_type_extra_parameters\n"
		"WHERE\n"
		"\tatom_type_set_name = ? AND name = ?\n"
		"ORDER BY\n"
		"\tparameter;";

	cppdb::statement stmt(basic::database::safely_prepare_statement(stmt_string, db_session));
	stmt.bind(1, atom_type_set_name);
	stmt.bind(2, atom_type.name());
	cppdb::result res(basic::database::safely_read_from_database(stmt));

	//std::string parameter;
	Real value;
	Size parameter_index(1);
	while ( res.next() ) {
		res >> value;
		atom_type.set_extra_parameter(parameter_index, value);
		++parameter_index;
	}
}


/// @details Do not add anything here.
///     These options are largely deprecated in favor of options that fine tune
///       score terms (EnergyMethodOptions), which are better controlled than
///       changing AtomTypeSet, which is almost a
///       global setting in Rosetta.
void
AtomTypeSet::legacy_command_line_post_processing()
{
	if ( option[ OptionKeys::chemical::enlarge_H_lj ] ) enlarge_h_lj_wdepth( *this ); // better set
	if ( option[ OptionKeys::chemical::no_hbonds_to_ether_oxygens ] ) {
		utility_exit_with_message( "-no_hbonds_to_ether_oxygens is deprecated. Instead use -score:hb_exclude_ether_oxygens or, less preferred, -chemical:unset_acceptor_ether_oxygens." );
	}
	if ( option[ OptionKeys::chemical::unset_acceptor_ether_oxygens ] ) unset_acceptor_ether_oxygens( *this );
}


void
AtomTypeSet::clone_atom_types_from_commandline()
{
	using std::string;
	using utility::vector1;

	if ( !basic::options::option[ basic::options::OptionKeys::chemical::clone_atom_types ].user() ) return;


	vector1< string > const clone_strings
		( basic::options::option[ basic::options::OptionKeys::chemical::clone_atom_types ]() );

	std::string const errmsg( "-clone_atom_types format should be:: -clone_atom_types <set1>:<atomname1>:<cloned-atomname1> <set2>:<atomname2>:<cloned-atomname2> ...; for example: '-chemical:clone_atom_types fa_standard:OOC:OOC2' ");


	for ( Size i=1; i<= clone_strings.size(); ++i ) {

		std::string const & mod( clone_strings[i] );

		Size const pos1( mod.find(":") );
		if ( pos1 == std::string::npos ) utility_exit_with_message(errmsg);
		std::string const atomset_tag( mod.substr(0,pos1) );
		if ( atomset_tag != name() ) continue;

		Size const pos2( mod.substr(pos1+1).find(":") );
		if ( pos2 == std::string::npos ) utility_exit_with_message(errmsg);
		std::string const atom_name( mod.substr(pos1+1,pos2) );
		if ( !has_atom( atom_name ) ) utility_exit_with_message(errmsg+". Nonexistent atomname: "+atom_name);

		std::string const new_atom_name( mod.substr(pos1+1+pos2+1) );
		if ( has_atom( new_atom_name ) ) {
			utility_exit_with_message(errmsg+". Duplicate atomname: "+new_atom_name);
		}

		tr.Trace << "clone_atom_types_from_commandline:: cloning " << atomset_tag << ' ' << atom_name <<
			" to " << new_atom_name << std::endl;

		Size const atom_index( atom_type_index( atom_name ) );

		auto* new_atom_type_ptr( new AtomType( *atoms_[ atom_index ] ) );
		new_atom_type_ptr->name( new_atom_name );
		atoms_.push_back( new_atom_type_ptr );
		atom_type_index_[ new_atom_name ] = atoms_.size();
	}


}

} // chemical
} // core

#ifdef SERIALIZATION

template < class Archive >
void core::chemical::serialize_atom_type_set( Archive & arc, core::chemical::AtomTypeSetCOP ptr )
{
	if ( ! ptr ) {
		bool ptr_is_nonnull( false );
		arc( CEREAL_NVP( ptr_is_nonnull ) );
	} else {
		bool ptr_is_nonnull( true );
		arc( CEREAL_NVP( ptr_is_nonnull ) );
		std::string typeset_name( ptr->name() ); // Assumes that the name can be used to extract it from the ChemicalManager
		arc( CEREAL_NVP( typeset_name ) );
	}
}
INSTANTIATE_FOR_OUTPUT_ARCHIVES( void, core::chemical::serialize_atom_type_set, core::chemical::AtomTypeSetCOP );

template < class Archive >
void core::chemical::deserialize_atom_type_set( Archive & arc, core::chemical::AtomTypeSetCOP & ptr )
{
	bool ptr_is_nonnull( true ); arc( ptr_is_nonnull );
	if ( ptr_is_nonnull ) {
		std::string typeset_name;
		arc( typeset_name );
		ptr = core::chemical::ChemicalManager::get_instance()->atom_type_set( typeset_name );
	} else {
		ptr = nullptr;
	}
}
INSTANTIATE_FOR_INPUT_ARCHIVES( void, core::chemical::deserialize_atom_type_set, core::chemical::AtomTypeSetCOP & );

CEREAL_REGISTER_DYNAMIC_INIT( core_chemical_AtomTypeSet )
#endif // SERIALIZATION
