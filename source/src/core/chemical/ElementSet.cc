// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/ElementSet.cc
/// @brief
/// @author P. Douglas Renfrew (renfrew@unc.edu)

// Unit headers
#include <core/chemical/ElementSet.hh>
#include <core/chemical/Element.hh>

// Project headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>  

// C++ headers
#include <fstream>
#include <iostream>

#include <utility/vector1.hh>


namespace core {
namespace chemical {

static basic::Tracer tr("core.chemical");

ElementSet::ElementSet():
	element_index_(),
	elements_()
{
}

ElementSet::~ElementSet() {}

/// @details Initialize an ElementSet from an external file "filename",
/// and set parameters and properties for each Element.
/// Refer to minirosetta_database_stock/chemical/element_properties.txt
/// for file format
///
void
ElementSet::read_file( std::string const & filename )
{
	 utility::io::izstream data( filename.c_str() ); 

	if ( !data.good() ) utility_exit_with_message( "Unable to open element file: "+filename );

	// now parse the rest of the file
	{
		using namespace basic;

		std::string line;
		// parse the header line
		getline( data, line ); // throw out the header line (currently it is just for file readability)
		while ( getline( data,line ) ) {
			utility::trim(line, " \t\n"); // remove leading and trailing spaces
			if ( line.empty() > 0 ) continue; //skip blank lines
			if ( line.find("#",0) == 0 ) continue; // skip comment lines

			std::istringstream l( line );
			std::string symbol, name;
			Real weight;
			Size z, mass;

			l >> z;

			if ( l.fail() ) {
				utility_exit_with_message("bad line: "+line);
			}

			l >> symbol >> name >> weight >> mass;

			if ( l.fail() ) {
				utility_exit_with_message("bad line: "+line);
			}

			Element* e( new Element( z, symbol, name, weight, mass) );

			// add this to the list
			elements_.push_back( e );
			if ( element_index_.count( symbol ) ) {
				utility_exit_with_message("ElementSet:: duplicate element symbol "+symbol);
			}
			element_index_[ symbol ] = elements_.size();
			tr.Debug << "New element: " << symbol << std::endl;
		}
	}
}

/// @details This function iterates over each element in the element_index_ map and
/// prints both keys. It is only used for debugging.
void
ElementSet::print_all_types()
{
	for( std::map< std::string, int >::const_iterator i = element_index_.begin(), e = element_index_.end(); i != e; ++i )
		{
			std::cout << (*i).first << " " << (*i).second << std::endl;
		}
}

} // chemical
} // core
