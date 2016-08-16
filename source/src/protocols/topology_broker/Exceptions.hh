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

class EXCN_TopologyBroker : public virtual utility::excn::EXCN_Msg_Exception {
	typedef EXCN_Msg_Exception Parent;
protected:
	EXCN_TopologyBroker() : EXCN_Msg_Exception( "" ){};
	virtual void show( std::ostream& os ) const {
		os << "\n[TopologyBroker Exception]: ";
		Parent::show( os );
	}
};

class EXCN_Input : public EXCN_TopologyBroker, public utility::excn::EXCN_BadInput {
	typedef EXCN_TopologyBroker Parent;
public:
	EXCN_Input( std::string const& msg ) : utility::excn::EXCN_Msg_Exception( msg ) {};
	virtual void show( std::ostream& os ) const {
		os << "*************** Error in Broker Setup: *************** \n";
		Parent::show( os );
		os << "\n\n********** Check your (inconsistent) input *************** \n";
	}
};

// this is more like a range error --- asking for something which isn't there
class EXCN_Unknown : public EXCN_TopologyBroker {
public:
	EXCN_Unknown( std::string const& msg ) : utility::excn::EXCN_Msg_Exception( msg ) {};
};

class EXCN_FailedBroking : public EXCN_TopologyBroker {
	typedef EXCN_TopologyBroker Parent;
public:
	EXCN_FailedBroking( std::string const& msg ) : utility::excn::EXCN_Msg_Exception( msg ) {};
	virtual void show( std::ostream& os ) const {
		Parent::show( os );
		os << "Failed to mediate between different Claimers. " << std::endl;
	}
};

class EXCN_FilterFailed : public EXCN_TopologyBroker {
	typedef EXCN_TopologyBroker Parent;
public:
	EXCN_FilterFailed( std::string const& msg ) : utility::excn::EXCN_Msg_Exception( msg ) {};

	virtual void show( std::ostream& os ) const {
		Parent::show( os );
		os << "[FILTER] failed... " << std::endl;
	}

};


}
}

#endif
