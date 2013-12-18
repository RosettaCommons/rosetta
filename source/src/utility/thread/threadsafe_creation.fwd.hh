// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

#ifndef INCLUDED_utility_thread_threadsafe_singleton_FWD_HH
#define INCLUDED_utility_thread_threadsafe_singleton_FWD_HH

namespace utility {
namespace thread {

template < class T >
inline
void
safely_create_singleton( T * & instance );

}
}

#endif
