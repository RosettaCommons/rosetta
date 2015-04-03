// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/Exceptions.hh
/// @details responsibilities:
/// @author Ken Jung
/// @author Christopher Miles (cmiles@uw.edu)
#ifndef INCLUDED_protocols_loophash_Exceptions_HH
#define INCLUDED_protocols_loophash_Exceptions_HH

// Project headers
#include <core/types.hh>

// Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

// C/C++
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace loophash {

class EXCN_DB_IO_Failed: public utility::excn::EXCN_Msg_Exception {
 public:
  EXCN_DB_IO_Failed( std::string filename, std::string action ) :
  utility::excn::EXCN_Msg_Exception( "Error in opening File='" + filename + "' for operation='" + action + "'" ),
	filename_( filename ), action_( action ) {};

 private:
  std::string filename_;
  std::string action_;
};

class EXCN_Invalid_Hashmap: public utility::excn::EXCN_Msg_Exception {
 public:
  EXCN_Invalid_Hashmap( core::Size size ) :
  utility::excn::EXCN_Msg_Exception( "Invalid hashmap of size " + utility::to_string( size ) + " requested." )
	// size_( size )
{};

 private:
  // KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
  // core::Size size_;
};

class EXCN_No_Extra_Data_To_Write: public utility::excn::EXCN_Msg_Exception {
 public:
  EXCN_No_Extra_Data_To_Write():
  utility::excn::EXCN_Msg_Exception( "Cannot write DB with extras if extra data is nonexistent." ) {};
};

class EXCN_Wrong_DB_Format: public utility::excn::EXCN_Msg_Exception {
 public:
  EXCN_Wrong_DB_Format( std::string filename ):
  utility::excn::EXCN_Msg_Exception( filename + " does not have the correct format." ),
	filename_( filename ) {};

 private:
  std::string filename_;
};

class EXCN_bbdb_Merge_Failed: public utility::excn::EXCN_Msg_Exception {
 public:
  EXCN_bbdb_Merge_Failed( std::string reason ):
  utility::excn::EXCN_Msg_Exception( reason ) {};

  EXCN_bbdb_Merge_Failed( bool /*masterlib_extra*/, bool /*secondlib_extra*/ ):
  utility::excn::EXCN_Msg_Exception( "LoopHashLibrary.extra_'s not equal" )
	// masterlib_extra_( masterlib_extra ),
	// secondlib_extra_( secondlib_extra )
{};

 private:
  // KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
  // bool masterlib_extra_;
  // KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
  // bool secondlib_extra_;
};

}  //namespace loophash
}  //namespace protocols

#endif
