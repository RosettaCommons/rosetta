// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_Exceptions_hh
#define INCLUDED_protocols_topology_broker_Exceptions_hh


// Unit Headers
//#include <protocols/topology_broker/Exceptions.fwd.hh>

// Utility Headers

#include <utility/excn/Exceptions.hh>


namespace protocols {
namespace topology_broker {

class EXCN_TopologyBroker : public utility::excn::Exception {
public:
	EXCN_TopologyBroker(char const *file, int line, std::string const & m) : utility::excn::Exception(file, line, "\n[TopologyBroker Exception]: " + m) {}
};

class EXCN_Input :  public utility::excn::BadInput {
public:
	EXCN_Input(char const *file, int line, std::string const& msg ) : utility::excn::BadInput(file, line, "") {
		add_msg("*************** Error in Broker Setup: *************** \n");
		add_msg(msg);
		add_msg("\n\n********** Check your (inconsistent) input *************** \n");
	}
};

// this is more like a range error --- asking for something which isn't there
class EXCN_Unknown : public EXCN_TopologyBroker {
public:
	using EXCN_TopologyBroker::EXCN_TopologyBroker;
};

class EXCN_FailedBroking : public EXCN_TopologyBroker {
public:
	EXCN_FailedBroking(char const *file, int line, std::string const& msg) : EXCN_TopologyBroker(file, line, msg) {
		add_msg("\nFailed to mediate between different Claimers...");
	}
};

class EXCN_FilterFailed : public EXCN_TopologyBroker {
public:
	EXCN_FilterFailed(char const *file, int line, std::string const& msg) : EXCN_TopologyBroker(file, line, msg) {
		add_msg("\n[FILTER] failed... ");
	}

};


}
}

#endif
