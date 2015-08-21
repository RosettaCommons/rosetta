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


#ifndef INCLUDED_protocols_noesy_assign_Exceptions_hh
#define INCLUDED_protocols_noesy_assign_Exceptions_hh


// Unit Headers
//#include <devel/topology_broker/Exceptions.fwd.hh>

// Utility Headers
#include <core/id/NamedAtomID.hh>
#include <utility/excn/Exceptions.hh>
#include <protocols/noesy_assign/PeakAssignment.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace noesy_assign {

class EXCN_NoesyAssign : public virtual utility::excn::EXCN_Msg_Exception {
	typedef EXCN_Msg_Exception Parent;
protected:
	EXCN_NoesyAssign() : EXCN_Msg_Exception( "" ){};
public:
	virtual void show( std::ostream& os ) const {
		os << "\n[NOE Exception]: ";
		Parent::show( os );
	}
};

class EXCN_UnknownAtomname : public EXCN_NoesyAssign {
public:
	EXCN_UnknownAtomname( std::string const& msg )
	: utility::excn::EXCN_Msg_Exception( msg ) {};
};


class EXCN_UnknownResonance : public EXCN_NoesyAssign {
public:
	EXCN_UnknownResonance( core::id::NamedAtomID atom, std::string const& msg )
	: utility::excn::EXCN_Msg_Exception( msg ), atom_( atom ) {};

	core::id::NamedAtomID const& atom() { return atom_; }
	virtual void show( std::ostream& os ) const {
		utility::excn::EXCN_Msg_Exception::show( os );
		os << "Resonance for atom " << atom_ << " not found ";
	}
private:
	core::id::NamedAtomID atom_;
};

class EXCN_AssignmentNotFound : public EXCN_NoesyAssign {
public:
	EXCN_AssignmentNotFound( PeakAssignment const& assignment, std::string const& msg )
	: utility::excn::EXCN_Msg_Exception( msg ), assignment_( assignment ) {};
	PeakAssignment assignment_;
};

class EXCN_FileFormat : public EXCN_NoesyAssign {
public:
	EXCN_FileFormat( std::string msg ) : utility::excn::EXCN_Msg_Exception( msg ) {};
};

} //namespace
}
#endif
