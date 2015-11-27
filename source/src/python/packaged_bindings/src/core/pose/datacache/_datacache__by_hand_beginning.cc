// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
///
/// @author Sergey Lyskov


#include "boost/python.hpp"

#include <core/pose/datacache/CacheableObserver.hh>
#include <basic/datacache/DataCache.hh>

namespace bp = boost::python;

void __datacache_by_hand_beginning__()
{ 
  // Add templated abstract base class wrapper for DataCache objects used in core.pose.datacache and core.pose
  boost::python::class_< basic::datacache::DataCache< core::pose::datacache::CacheableObserver >, utility::pointer::shared_ptr< basic::datacache::DataCache< core::pose::datacache::CacheableObserver > > >( "DataCache_CacheableObserver" );
}
