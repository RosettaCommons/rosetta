// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/sequence/SSManager.hh
/// @details header file for class of to convert SS to numeric object. Used by the hdf5 data structure to shring DB storage requirements
/// @author TJ Brunette( tjbrunette@gmail.com )

#ifndef INCLUDED_core_sequence_SSManager_hh
#define INCLUDED_core_sequence_SSManager_hh

#include <core/sequence/SSManager.fwd.hh>

// package headers
#include <core/pose/Pose.fwd.hh>

// utility headers
#include <core/types.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace sequence {

/// @brief struct



/// @brief manager for abego
class SSManager : public utility::pointer::ReferenceCount {
public: // typedef
	typedef core::Size Size;
	typedef std::string String;

public:

	/// @brief default constructor
	SSManager();

	/// @brief value constructor
	virtual ~SSManager() ;

public:
	/// @brief transform abego index to symbol
	char index2symbol( Size const & idx );

	/// @brief transform abego symbol to index
	Size symbol2index( char const & symbol );

	/// @brief transform abego symbol string to base5 index
	Size symbolString2index( std::string symbolString);

	/// @brief transform abego string to abego base5 index
	std::string index2symbolString( Size base5index,Size length);



};

} // namespace util
} // namespace core

#endif /* INCLUDED_core_sequence_hh */
