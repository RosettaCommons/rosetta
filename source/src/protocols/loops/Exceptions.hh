// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// what's the point of having annotations if they're wrong? Serendipitously,
// this entire file is an Exception to the coding guidelines!
/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange


#ifndef INCLUDED_protocols_loops_Exceptions_HH
#define INCLUDED_protocols_loops_Exceptions_HH

// Utility headers
#include <utility/excn/Exceptions.hh>

namespace protocols {
namespace loops {

class EXCN_Loop_not_closed : public utility::excn::EXCN_Msg_Exception {
  typedef EXCN_Msg_Exception Parent;
public:
  EXCN_Loop_not_closed( std::string msg = "") : EXCN_Msg_Exception( "failed to close loop " + msg ){};
};

} // namespace loops
} // namespace protocols

#endif //INCLUDED_protocols_loops_Exceptions_HH
