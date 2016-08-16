// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/constraints/AlignerFactory.hh
/// @brief
/// @author James Thompson <tex@u.washington.edu>

#ifndef INCLUDED_core_sequence_AlignerFactory_hh
#define INCLUDED_core_sequence_AlignerFactory_hh

// Package headers
#include <core/sequence/Aligner.fwd.hh>

// Project headers

// C++ Headers

#include <iostream>


namespace core {
namespace sequence {

class AlignerFactory {
public:
	/// @brief returns an AlignerOP
	static AlignerOP get_aligner( std::string const & type );
};

} // sequence
} // core

#endif
