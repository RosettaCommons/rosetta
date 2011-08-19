// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/init.hh
/// @brief  Core Init functions
/// @author Sergey Lyskov
///


#ifndef INCLUDED_core_init_hh
#define INCLUDED_core_init_hh

// STL
#include <string>

#include <utility/vector1.hh>
#include <numeric/random/random.hh>

namespace core {

	/// @brief Init basic core systems: options system, random system.
	void init(int argc, char * argv []);

	/// @brief wrapper for core system Init
	void init( utility::vector1<std::string> const & args );

	/// @brief Initialize random generator systems (and send debug io to tracer with seed/mode info).
	void init_random_generators(int const start_seed, numeric::random::RND_RunType run_type, std::string const & RGtype);
} // namespace core

#endif // INCLUDED_core_init_HH
