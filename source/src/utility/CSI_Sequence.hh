// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/utility/CSI_Sequence.hh
/// @brief  Terminal ASCII codes
/// @author Sergey Lyskov, @modified by Caleb Geniesse


#ifndef INCLUDED_utility_CSI_Sequence_hh
#define INCLUDED_utility_CSI_Sequence_hh

#include <utility/CSI_Sequence.fwd.hh>

#include <string>
#include <ostream>  // for ostream, operator<<


namespace utility {

/// @brief Class to hold all Terminal ASCII codes as static data for CSI_Sequence.
///        Note: that on non-tty terminals the codes will not print to avoid polution of Rosetta logs.
class CSI_Sequence
{
public:
	/// @brief constructor
	CSI_Sequence(std::string const & sequence = "") : sequence_( sequence ) {}

	// @details This is explicit such that we don't get odd int conversions based on
	explicit CSI_Sequence(CSI::CSI_Enum sq );

	// No! No implicit conversion to string : we need the special output overloading.
	// operator std::string() const { return sequence_; }

	CSI_Sequence operator+ ( CSI_Sequence const & right ) const {
		return CSI_Sequence( sequence_ + right.sequence_ );
	}

	/// @brief operator to output our sequence so we can write: std::cout << CSI_SequenceObject
	friend std::ostream & operator << (std::ostream & os, CSI_Sequence const &sq) {
		if ( ! CSI_Sequence::suppress_CSI_seq() ) { os << sq.sequence_; }
		return os;
	}

	/// @brief If called, suppress all future printing of CSI codes
	static void suppress_CSI_codes() {
		suppress_CSI_seq() = true;
	}

private:
	std::string sequence_;

	// Should we suppress the use of CSI Sequences?
	static bool & suppress_CSI_seq();
};

/// @details Initialize from a CSI_Enum
///          Codes below is all Hogwarts-approved magic numbers, so do not modify them.
///          For reference see: http://en.wikipedia.org/wiki/ANSI_escape_code#CSI_codes
inline
CSI_Sequence::CSI_Sequence( CSI::CSI_Enum sq ) {
	switch( sq ) {
	case CSI::Nothing :   sequence_ = ""; break;
	case CSI::Reset :     sequence_ = "\x1b[0m"; break;
	case CSI::Bold :      sequence_ = "\x1b[1m"; break;
	case CSI::Underline : sequence_ = "\x1b[4m"; break;
	case CSI::Black :     sequence_ = "\x1b[30m"; break;
	case CSI::Red :       sequence_ = "\x1b[31m"; break;
	case CSI::Green :     sequence_ = "\x1b[32m"; break;
	case CSI::Yellow :    sequence_ = "\x1b[33m"; break;
	case CSI::Blue :      sequence_ = "\x1b[34m"; break;
	case CSI::Magenta :   sequence_ = "\x1b[35m"; break;
	case CSI::Cyan :      sequence_ = "\x1b[36m"; break;
	case CSI::White :     sequence_ = "\x1b[37m"; break;
	case CSI::Default :   sequence_ = "\x1b[39m"; break;
	case CSI::bgBlack :   sequence_ = "\x1b[40m"; break;
	case CSI::bgRed :     sequence_ = "\x1b[41m"; break;
	case CSI::bgGreen :   sequence_ = "\x1b[42m"; break;
	case CSI::bgYellow :  sequence_ = "\x1b[43m"; break;
	case CSI::bgBlue :    sequence_ = "\x1b[44m"; break;
	case CSI::bgMagenta : sequence_ = "\x1b[45m"; break;
	case CSI::bgCyan :    sequence_ = "\x1b[46m"; break;
	case CSI::bgWhite :   sequence_ = "\x1b[47m"; break;
	case CSI::bgDefault : sequence_ = "\x1b[49m"; break;
	default :             sequence_ = ""; break;
	};
}

/// @details Functions to return CSI_Seqeunce objects.
///          These are not static constants, due to the static initilization order fiasco.
inline CSI_Sequence CSI_Nothing()   { return CSI_Sequence(CSI::Nothing)  ; }
inline CSI_Sequence CSI_Reset()     { return CSI_Sequence(CSI::Reset)    ; }
inline CSI_Sequence CSI_Bold()      { return CSI_Sequence(CSI::Bold)     ; }
inline CSI_Sequence CSI_Underline() { return CSI_Sequence(CSI::Underline); }
inline CSI_Sequence CSI_Black()     { return CSI_Sequence(CSI::Black)    ; }
inline CSI_Sequence CSI_Red()       { return CSI_Sequence(CSI::Red)      ; }
inline CSI_Sequence CSI_Green()     { return CSI_Sequence(CSI::Green)    ; }
inline CSI_Sequence CSI_Yellow()    { return CSI_Sequence(CSI::Yellow)   ; }
inline CSI_Sequence CSI_Blue()      { return CSI_Sequence(CSI::Blue)     ; }
inline CSI_Sequence CSI_Magenta()   { return CSI_Sequence(CSI::Magenta)  ; }
inline CSI_Sequence CSI_Cyan()      { return CSI_Sequence(CSI::Cyan)     ; }
inline CSI_Sequence CSI_White()     { return CSI_Sequence(CSI::White)    ; }
inline CSI_Sequence CSI_Default()   { return CSI_Sequence(CSI::Default)  ; }
inline CSI_Sequence CSI_bgBlack()   { return CSI_Sequence(CSI::bgBlack)  ; }
inline CSI_Sequence CSI_bgRed()     { return CSI_Sequence(CSI::bgRed)    ; }
inline CSI_Sequence CSI_bgGreen()   { return CSI_Sequence(CSI::bgGreen)  ; }
inline CSI_Sequence CSI_bgYellow()  { return CSI_Sequence(CSI::bgYellow) ; }
inline CSI_Sequence CSI_bgBlue()    { return CSI_Sequence(CSI::bgBlue)   ; }
inline CSI_Sequence CSI_bgMagenta() { return CSI_Sequence(CSI::bgMagenta); }
inline CSI_Sequence CSI_bgCyan()    { return CSI_Sequence(CSI::bgCyan)   ; }
inline CSI_Sequence CSI_bgWhite()   { return CSI_Sequence(CSI::bgWhite)  ; }
inline CSI_Sequence CSI_bgDefault() { return CSI_Sequence(CSI::bgDefault); }

namespace CSI {

inline
std::ostream & operator << (std::ostream & os, CSI_Enum const &sq) {
	os << CSI_Sequence( sq );
	return os;
}

} // CSI

} // utility

#endif // INCLUDED_utility_CSI_Sequence_HH
