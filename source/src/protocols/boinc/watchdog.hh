// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_boinc_watchdog_hh
#define INCLUDED_protocols_boinc_watchdog_hh

#include <string>

#ifndef _WIN32
#include "pthread.h"
#endif

namespace protocols {
namespace boinc {
namespace watchdog {

// protocols can set this pose as the global bailout - if the watchdog kicks in it will write out *this*
// pose and give it a special label to be identified as the Bailout ( W_xxx )
// Currently this is only set by the CheckPointer.
#ifdef WIN32
#else
extern pthread_mutex_t bailout_mutex;
#endif
extern std::string bailout_silent_structure;
extern std::string bailout_silent_structure_header;

void
watchdog_start();

void
watchdog_finish();

void*
main_watchdog( void* );

} // namespace watchdog
} // namespace boinc
} // namespace protocols


#endif // INCLUDED_protocols_boinc_watchdog_HH
