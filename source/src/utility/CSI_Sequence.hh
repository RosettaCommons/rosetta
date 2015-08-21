// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/utility/CSI_Sequence.hh
/// @brief  Terminal ASCII codes
/// @author Sergey Lyskov, @modified by Caleb Geniesse


#ifndef INCLUDED_utility_CSI_Sequence_hh
#define INCLUDED_utility_CSI_Sequence_hh

#include <utility/CSI_Sequence.fwd.hh>

#include <ostream>
#include <string>


namespace utility {

/// @brief Class to hold all Terminal ASCII codes as static data for CSI_Sequence.
///        Note: that on non-tty terminals all codes will initialized as empty to avoid polution of Rosetta logs.
class CSI_Sequence
{
public:
	/// @brief constructor
	CSI_Sequence(std::string sequence_);

	/// @brief operator to output our sequence so we can write: std::cout << CSI_SequenceObject
	friend std::ostream & operator << (std::ostream & os, CSI_Sequence const &sq) { os << sq.sequence; return os; }

private:
	std::string sequence;
};


/// @details Constant static string objects to hold various ASCII CSI codes
///          Codes below is all Hogwarts-approved magic numbers, so do not modify them.
///          For reference see: http://en.wikipedia.org/wiki/ANSI_escape_code#CSI_codes
static std::string const  CSI_Reset("\x1b[0m"),
CSI_Bold("\x1b[1m"),
CSI_Underline("\x1b[4m"),
CSI_Black("\x1b[30m"),
CSI_Red("\x1b[31m"),
CSI_Green("\x1b[32m"),
CSI_Yellow("\x1b[33m"),
CSI_Blue("\x1b[34m"),
CSI_Magenta("\x1b[35m"),
CSI_Cyan("\x1b[36m"),
CSI_White("\x1b[37m"),
CSI_bgBlack("\x1b[40m"),
CSI_bgRed("\x1b[41m"),
CSI_bgGreen("\x1b[42m"),
CSI_bgYellow("\x1b[43m"),
CSI_bgBlue("\x1b[44m"),
CSI_bgMagenta("\x1b[45m"),
CSI_bgCyan("\x1b[46m"),
CSI_bgWhite("\x1b[47m");


} // utility

#endif // INCLUDED_utility_CSI_Sequence_HH
