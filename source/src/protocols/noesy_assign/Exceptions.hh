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


#ifndef INCLUDED_protocols_noesy_assign_Exceptions_hh
#define INCLUDED_protocols_noesy_assign_Exceptions_hh


// Utility Headers
#include <core/id/NamedAtomID.hh>
#include <utility/excn/Exceptions.hh>
#include <protocols/noesy_assign/PeakAssignment.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace noesy_assign {

class EXCN_NoesyAssign : public utility::excn::Exception {
public:
	EXCN_NoesyAssign(char const *file, int line, std::string const& msg) : utility::excn::Exception(file, line, "\n[NOE Exception]: " + msg) {}
};

class EXCN_UnknownAtomname : public EXCN_NoesyAssign {
public:
	using EXCN_NoesyAssign::EXCN_NoesyAssign;
};


class EXCN_UnknownResonance : public EXCN_NoesyAssign {
public:
	EXCN_UnknownResonance(char const *file, int line, core::id::NamedAtomID atom, std::string const& msg ) : EXCN_NoesyAssign(file, line, msg), atom_( atom ) {
		utility::excn::Exception::add_msg("Resonance for atom " + atom_.to_string() +" not found!");
	}

	core::id::NamedAtomID const& atom() { return atom_; }
private:
	core::id::NamedAtomID atom_;
};

class EXCN_AssignmentNotFound : public EXCN_NoesyAssign {
public:
	EXCN_AssignmentNotFound(char const *file, int line, PeakAssignment const& assignment, std::string const& msg ) : EXCN_NoesyAssign(file, line, msg), assignment_(assignment) {};
	PeakAssignment assignment_;
};

class EXCN_FileFormat : public EXCN_NoesyAssign {
public:
	using EXCN_NoesyAssign::EXCN_NoesyAssign;
};

} //namespace
}
#endif
