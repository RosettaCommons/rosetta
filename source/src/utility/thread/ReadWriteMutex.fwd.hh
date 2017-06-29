// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/thread/ReadWriteMutex.hh
/// @brief  Declaration of the ReadWriteMutex class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_utility_thread_ReadWriteMutex_FWD_HH
#define INCLUDED_utility_thread_ReadWriteMutex_FWD_HH

#ifdef MULTI_THREADED

#include <utility/pointer/owning_ptr.hh>

namespace utility {
namespace thread {

class ReadWriteMutex;

typedef utility::pointer::shared_ptr< ReadWriteMutex > ReadWriteMutexOP;
typedef utility::pointer::shared_ptr< ReadWriteMutex const > ReadWriteMutexCOP;


}
}

#endif // MULTI_THREADED

#endif
