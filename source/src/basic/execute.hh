// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/basic/execute.hh
/// @brief  Various functions to run external executables
/// @author Sergey Lyskov
///
/// Note: this code is placed in to 'basic' level instead of utility because we want to have access to Tracers so verbosity of all external calls could be controlled

#ifndef INCLUDED_basic_execute_hh
#define INCLUDED_basic_execute_hh

#include <string>

namespace basic {


/// @brief Struct that hold result code + output of external process
struct ExecutionResult {
public:
	int result; // result code, zero if no error and otherwise will contain error code
	std::string output;

	#ifdef HAS_MOVE_SEMANTICS
		ExecutionResult(ExecutionResult &&r) =default; //{ result=r.result; output=std::move(r.output); }
    #endif
};

/// @brief excute provided command_line though shell and return exit_code and output
ExecutionResult execute(std::string const & message, std::string const & command_line, bool terminate_on_failure=true, bool silent=false);


} // namespace basic



#endif // INCLUDED_basic_execute_hh
