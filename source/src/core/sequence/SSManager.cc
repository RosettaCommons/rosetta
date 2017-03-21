// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/sequence/SSManager.cc
/// @brief class for converting SS to int. For use with the ss hashed fragment database
/// @author TJ Brunette ( tjbrunette@gmail.com )

#include <core/chemical/ResidueType.hh>
#include <core/sequence/SSManager.hh>
#include <core/pose/Pose.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <cmath>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "core.sequence.SSManager" );

namespace core {
namespace sequence {


/// @brief default constructor
SSManager::SSManager() : utility::pointer::ReferenceCount()
{
}
// @brief Auto-generated virtual destructor
SSManager::~SSManager() {}


/// @brief transform ss index to symbol
char
SSManager::index2symbol( Size const & idx )
{
	switch( idx ) {
	case 1 :
		return 'H';
	case 2 :
		return 'L';
	case 3 :
		return 'E';
	default :
		TR.Error << "Unrecognized ss index: " << idx << std::endl;
		utility_exit_with_message("Unrecognized ss index");
		return 0;
	}
}


/// @brief transform abego symbol to index
Size
SSManager::symbol2index( char const & symbol )
{
	switch( symbol ) {
	case 'H' :
	case 'h' :
		return 1;
	case 'L' :
	case 'l' :
		return 2;
	case 'E' :
	case 'e' :
		return 3;
	default :
		TR.Error << "Unrecognized ss index: " << symbol << std::endl;
		utility_exit_with_message("Unrecognized ss index");
		return 0;
	}
}

/// @brief transform abego symbol string to base5 index. This is used to quickly pool the abego from Alex's hd5 database
Size SSManager::symbolString2index( std::string symbolString){
	std::string allowedChar = "HLE";
	for ( Size ii=0; ii<symbolString.size(); ++ii ) { //check for only allowed characters
		if ( allowedChar.find(symbolString.substr(ii,1)) == std::string::npos ) {
			utility_exit_with_message("Looking for " + symbolString.substr(ii,1) + " which shouldn't exist");
		}
	}
	Size index = 0;
	for ( Size ii=0; ii<symbolString.size(); ++ii ) {
		Size symbolValue = symbol2index( symbolString[ii]);
		Size ss_index_radix = pow(3,ii);
		index += ss_index_radix*(symbolValue-1);
	}
	return(index);
}

/// @brief transform abego symbol string to base5 index. This is used to quickly pool the abego from Alex's hd5 database
std::string SSManager::index2symbolString( Size index,Size length){
	Size tmp_index=index;
	std::string symbolString ="";
	for ( Size ii=0; ii<length; ++ii ) {
		Size index = tmp_index % 3;
		tmp_index = tmp_index/3;
		symbolString+=index2symbol(index+1);
	}
	return(symbolString);
}


} // namespace sequence
} // namespace core

