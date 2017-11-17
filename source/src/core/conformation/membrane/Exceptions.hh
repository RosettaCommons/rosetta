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
class EXCN_Membrane : public utility::excn::Exception {

public:

	EXCN_Membrane(char const *file, int line, std::string const &m) : utility::excn::Exception(file, line, "\n[Membrane Exception]: ") {
		add_msg(m);
	};

}; // class EXCN_Membrane

/// @brief Resource Manager Exception
class EXCN_Resource_Definition : public EXCN_Membrane {

public:

	// Constructor
	EXCN_Resource_Definition(char const *file, int line, std::string const & m) : EXCN_Membrane(file, line, "") {
		add_msg("========================== Resource Loader Exception ==========================\n"
				"Please refine your Resource Definition file and try again. For instructions on how to \n"
				"properly construct your resource definition file for the membrane code, please visit the \n"
				"Rosetta documentation");
		add_msg(m);
	};

};


class EXCN_Illegal_Arguments : public EXCN_Membrane {

public:

	// Constructor
	EXCN_Illegal_Arguments(char const *file, int line, std::string const &m) : EXCN_Membrane(file, line, "Illegal Arguments Exception!\n")
	{
		add_msg(m);
	};

}; // class EXCN_Membrane_Bounds

/// @brief Membrane Out of Bounds Exception
class EXCN_Membrane_Bounds : public EXCN_Membrane {

public:

	// Constructor
	EXCN_Membrane_Bounds(char const *file, int line, std::string const & m) : EXCN_Membrane(file, line, "Membrane out of Bounds Exception!")
	{
		add_msg(m);
	};

}; // class EXCN_Membrane_Bounds

/// @brief Virtual Residue Definition Exception
class EXCN_VirtualRsd : public EXCN_Membrane {

public:

	// Constructor
	EXCN_VirtualRsd(char const *file, int line, std::string const & m) : EXCN_Membrane(file, line, "Virtual Residue Definition Exception")
	{
		add_msg(m);
	};

}; // class EXCN_VirtualRsd

/// @brief Non Membrane Pose Exception
class EXCN_NonMembrane : public EXCN_Membrane {

public:

	// Contructor
	EXCN_NonMembrane(char const *file, int line, std::string const & m) : EXCN_Membrane(file, line, "Pose is not initialized as a membrane pose!")
	{
		add_msg(m);
	};

}; // class EXCN_NonMembrane

/// @brief Fold tree Exception for membrane proteins
class EXCN_MembraneFoldTree : public EXCN_Membrane {

public:

	// Constructor
	EXCN_MembraneFoldTree(char const *file, int line, std::string const & m) : EXCN_Membrane(file, line, "Membrane FoldTree Error!")
	{
		add_msg(m);
	};

}; // class EXCN_MembraneFoldTree

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_Exceptions_hh
