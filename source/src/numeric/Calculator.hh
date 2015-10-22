// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  numeric/Calculator.hh
/// @brief a string-input calculator
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#ifndef INCLUDED_numeric_Calculator_hh
#define INCLUDED_numeric_Calculator_hh

#include <numeric/Calculator.fwd.hh>

#include <numeric/types.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <iosfwd>
#include <map>

// IMPORTANT - Keep all the boost crazyness confined to the .cpp,
//    so we limit the associated template madness to a single compilation unit.

namespace numeric {

class Calculator : public utility::pointer::ReferenceCount {

public:
	virtual ~Calculator();

	Calculator(std::string const & equation);

	/// @brief Calculate the value of the equation, putting it in output
	/// Return true if the computation failed.
	bool compute( std::map<std::string, Real> & values, Real & output ) const;

private:

	std::string equation_;
};

} // numeric

#endif // INCLUDED_numeric_Calculator_HH


