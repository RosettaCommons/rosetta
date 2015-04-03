// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange


#ifndef INCLUDED_core_id_Exceptions_hh
#define INCLUDED_core_id_Exceptions_hh


// Unit Headers
//#include <protocols/topology_broker/Exceptions.fwd.hh>

// Utility Headers
#include <core/id/NamedAtomID.hh>
#include <utility/excn/Exceptions.hh>


namespace core {
namespace id {

class EXCN_AtomNotFound: public utility::excn::EXCN_Msg_Exception {
public:
  EXCN_AtomNotFound( NamedAtomID const& );
	NamedAtomID const& atom() { return id_; }
private:
  NamedAtomID id_;
};


}
}

#endif
