// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @file Metapatch.cc
///
/// @brief A structure for applying multiple Patches procedurally
///
/// @details
/// Suppose you want to chlorinate all 'chlorinatable positions' and have a bajillion RTs result. NOW. YOU. CAN.
/////////////////////////////////////////////////////////////////////////


// Unit headers
#include <core/chemical/Patch.hh>
#include <core/chemical/Metapatch.hh>
#include <core/chemical/AtomPropertiesManager.hh>

// Basic header
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <fstream>


namespace core {
namespace chemical {

/// @details Auto-generated virtual destructor
Metapatch::~Metapatch() = default;

static THREAD_LOCAL basic::Tracer tr( "core.chemical" );

Metapatch::Metapatch() :
	name_( std::string() ),
	types_( utility::vector1< std::string >() ),
	//cases_( utility::vector1< PatchCase >(1,PatchCase()) ),
	pertinent_property_( NO_ATOM_PROPERTY ),
	case_lines_( utility::vector1< std::string >() ),
	selector_( ResidueTypeSelector() )
{}

/// @details - first read in all lines from the file, discarding # comment lines
/// - parse input lines for Metapatch name and variant types (NAME, TYPES)
/// - parse input lines for general ResidueTypeSelector defined for this Patch (BEGIN_SELECTOR, END_SELECTOR)
/// - parse input lines to create each case accordingly (BEGIN_CASE, END_CASE)
/// @note keep the order to avoid triggering parsing errors
void
Metapatch::read_file( std::string const & filename )
{
	// clear old data
	tr.Debug << "Reading metapatch file: " << filename << std::endl;

	name_ = "";
	types_.clear();
	selector_.clear();

	utility::vector1< std::string > lines;
	{ // read the lines file
		utility::io::izstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message("Cannot find metapatch file: "+filename);
		}
		std::string line;
		while ( getline( data,line ) ) {
			std::string const tag( tag_from_line( line ) );
			if ( tag.size() && tag[0] != '#' ) lines.push_back( line );
		}
	}

	// misc parsing
	for ( uint i=1; i<= lines.size(); ++i ) {
		std::istringstream l(lines[i]);
		std::string tag;
		l >> tag;
		if ( tag == "NAME" ) {
			l >> name_;
		} else if ( tag == "TYPES" ) {
			std::string t;
			l >> t;
			while ( !l.fail() ) {
				types_.push_back( t );
				l >> t;
			}
		} else if ( tag == "PERTINENT_PROPERTY" ) {
			std::string t;
			l >> t;
			pertinent_property_ = AtomPropertiesManager::property_from_string( t );
		}
	}

	// build the residue selector
	{
		bool in_selector( false );
		for ( uint i=1; i<= lines.size(); ++i ) {
			std::string tag( tag_from_line( lines[i] ) );
			if ( tag == "BEGIN_CASE" ) {
				break;
			} else if ( tag == "BEGIN_SELECTOR" ) {
				debug_assert( !in_selector );
				in_selector = true;
			} else if ( tag == "END_SELECTOR" ) {
				debug_assert( in_selector );
				in_selector = false;
			} else if ( in_selector ) {
				selector_.add_line( lines[i] );
			}
		}
	}

	// get the cases
	utility::vector1< std::string > case_lines;
	bool in_case( false );
	while ( !lines.empty() ) {
		// look for a case
		std::string tag( tag_from_line( lines[1] ) );
		if ( tag == "BEGIN_CASE" ) {
			debug_assert( case_lines.empty() );
			debug_assert( !in_case );
			in_case = true;
		} else if ( tag == "END_CASE" ) {
			PatchCaseOP new_case( case_from_lines( case_lines ) );
			// There must only be one case in a metapatch, for now...
			case_lines_ = case_lines;
			if ( new_case ) cases_.push_back( new_case );
			case_lines.clear();
			in_case = false;
		} else if ( in_case ) case_lines.push_back( lines[1] );

		lines.erase( lines.begin() );
	}
}

PatchCOP
Metapatch::get_one_patch( /*ResidueType const & rsd_type, */std::string const & atom_name ) const
{
	PatchOP p( new Patch );

	Size first = atom_name.find_first_not_of(' ');
	Size last = atom_name.find_last_not_of(' ');

	std::string trimmed_atom = atom_name.substr( first, last-first+1 );

	// Have to explicitly set this, since we aren't calling Patch::read_file :-(
	// If ever we want a replacing metapatch (why??) make metapatch::read_file
	// grab a tag and pass that value on here.
	p->replaces_residue_type( false );

	// Loop through cases until we find one that applies
	utility::vector1< std::string > substituted_lines = case_lines_;
	// Replace mentions of "blank" with the atom name.
	for ( Size l = 1; l <= substituted_lines.size(); ++l ) {
		if ( substituted_lines[l].find("blank") != std::string::npos ) {
			substituted_lines[l] = substituted_lines[l].replace(
				substituted_lines[l].find("blank"), 5, atom_name );
		}
	}

	utility::vector1< std::string > types = types_;
	for ( Size i = 1; i <= types_.size(); ++i ) {
		if ( types[i].find( "blank" ) != std::string::npos ) {
			types[i] = types[i].replace( types[i].find("blank"), 5, trimmed_atom );
		}
	}

	PatchCaseOP new_case( case_from_lines( substituted_lines ) );
	p->add_case( new_case );
	p->set_name( "MP-" + trimmed_atom + "-" + name_ );
	p->set_selector( selector_ );
	p->types( types );

	return p;
}

utility::vector1< std::string >
Metapatch::atoms( ResidueType const & rsd_type ) const
{
	utility::vector1< std::string > good_atoms;
	utility::vector1< std::string > patch_names = get_patch_names( rsd_type );
	for ( Size i = 1; i <= rsd_type.natoms(); ++i ) {

		if ( !meets_requirements( rsd_type, i ) ) continue;

		Size first = rsd_type.atom_name( i ).find_first_not_of(' ');
		Size last = rsd_type.atom_name( i ).find_last_not_of(' ');
		std::string trimmed_atom = rsd_type.atom_name( i ).substr( first, last-first+1 );

		bool cont = false;
		for ( Size pn = 1; pn <= patch_names.size(); ++pn ) {
			if ( patch_names[ pn ].find( rsd_type.name3() ) == std::string::npos ) continue;
			utility::vector1< std::string > elems = utility::string_split( patch_names[pn], '-' );
			if ( trimmed_atom < elems[2] ) cont = true;
		}
		if ( cont ) continue;

		good_atoms.push_back( trimmed_atom );
	}
	return good_atoms;
}

} // chemical
} // core
