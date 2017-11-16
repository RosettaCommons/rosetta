// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/basic/TracerImpl.fwd.hh
/// @brief  Tracer IO system
/// @author Sergey Lyskov
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_basic_TracerImpl_fwd_hh
#define INCLUDED_basic_TracerImpl_fwd_hh

#include <ostream>
#include <utility/pointer/owning_ptr.hh>

namespace basic {

/// @brief
/// Priority levels for T() and Tracer object, modeled on the log4j project and its offspring.
/// Priorities in Tracer are still ints so users can pass other arbitrary integer values (for now).
enum TracerPriority {
	t_fatal   = 0,   //< The FATAL level designates very severe error events that will presumably lead the application to abort.
	t_error   = 100, //< The ERROR level designates error events that might still allow the application to continue running.
	t_warning = 200, //< The WARN level designates potentially harmful situations.
	t_info    = 300, //< The INFO level designates informational messages that highlight the progress of the application at coarse-grained level.
	t_debug   = 400, //< The DEBUG level designates fine-grained informational events that are most useful to debug an application.
	t_trace   = 500  //< The TRACE level designates finer-grained informational events than the DEBUG level.
};

class TracerImpl;

typedef utility::pointer::shared_ptr< TracerImpl > TracerImplOP;
typedef utility::pointer::shared_ptr< TracerImpl const > TracerImplCOP;

}

#endif // INCLUDED_basic_tracerImpl_FWD_HH


