// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/conformation/membrane/Exceptions.hh
///
/// @brief      Exception Hierarchy for Membrane framework Code
/// @details    Exception hierarchy extended from utility exceptions - specific
///    to when things go wring in the membrane code. Should be used exclusively in
///    the membrane framework and within JD2
///    Last Updated: 7/23/14
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_Exceptions_hh
#define INCLUDED_core_conformation_membrane_Exceptions_hh

// Utility Headers
#include <utility/excn/Exceptions.hh>

namespace core {
namespace conformation {
namespace membrane {

/// @brief Parent Exception - Exception Membrane
class EXCN_Membrane : public virtual utility::excn::EXCN_Msg_Exception {

public:

	// Constructor
	EXCN_Membrane() : utility::excn::EXCN_Msg_Exception( "" ) {};

	// Show Exception
	virtual void show( std::ostream & os ) const {
		os << "\n[Membrane Exception]:";
		EXCN_Msg_Exception::show( os );
	}

}; // class EXCN_Membrane

/// @brief Resource Manager Exception
class EXCN_Resource_Definition : public EXCN_Membrane {

public:

	// Constructor
	EXCN_Resource_Definition( std::string const & msg ) : utility::excn::EXCN_Msg_Exception( msg ){};

	// Show Exception
	virtual void show( std::ostream & os ) const {
		os << "========================== Resource Loader Exception ==========================\n";
		os << "Please refine your Resource Definition file and try again. For instructions on how to \n";
		os << "properly construct your resource definition file for the membrane code, please visit the \n";
		os << "Rosetta documentation";
		EXCN_Msg_Exception::show(os);
	}

}; // class EXCN_Resource_Definition

/// @brief Illegal Arguments Exception
class EXCN_Illegal_Arguments : public EXCN_Membrane {

public:

	// Constructor
	EXCN_Illegal_Arguments( std::string const & msg ) : utility::excn::EXCN_Msg_Exception( msg ){};

	// Show Exception
	virtual void show( std::ostream & os ) const {
		os << "Illegal Arguments Exception!";
		EXCN_Msg_Exception::show(os);

	}

}; // class EXCN_Membrane_Bounds

/// @brief Membrane Out of Bounds Exception
class EXCN_Membrane_Bounds : public EXCN_Membrane {

public:

	// Constructor
	EXCN_Membrane_Bounds( std::string const & msg ) : utility::excn::EXCN_Msg_Exception( msg ){};

	// Show Exception
	virtual void show( std::ostream & os ) const {
		os << "Membrane out of Bounds Exception!";
		EXCN_Msg_Exception::show(os);
	}

}; // class EXCN_Membrane_Bounds

/// @brief Virtual Residue Definition Exception
class EXCN_VirtualRsd : public EXCN_Membrane {

public:

	// Constructor
	EXCN_VirtualRsd( std::string const & msg ) : utility::excn::EXCN_Msg_Exception( msg ){};

	// Show Exception
	virtual void show( std::ostream & os ) const {
		os << "Virtual Residue Definition Exception";
		EXCN_Msg_Exception::show(os);
	}

}; // class EXCN_VirtualRsd

/// @brief Non Membrane Pose Exception
class EXCN_NonMembrane : public EXCN_Membrane {

public:

	// Contructor
	EXCN_NonMembrane( std::string const & msg ) : utility::excn::EXCN_Msg_Exception( msg ){};

	// Show Exception
	virtual void show( std::ostream & os ) const {
		os << "Pose is not initialized as a membrane pose!";
		EXCN_Msg_Exception::show(os);
	}

}; // class EXCN_NonMembrane

/// @brief Fold tree Exception for membrane proteins
class EXCN_MembraneFoldTree : public EXCN_Membrane {

public:

	// Constructor
	EXCN_MembraneFoldTree( std::string const & msg ) : utility::excn::EXCN_Msg_Exception( msg ){};

	// Show exception
	virtual void show( std::ostream & os ) const {
		os << "Membrane FoldTree Error!";
		EXCN_Msg_Exception::show(os);
	}

}; // class EXCN_MembraneFoldTree

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_Exceptions_hh

