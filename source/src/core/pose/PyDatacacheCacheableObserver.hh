// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/PyDatacacheCacheableObserver.hh
/// @brief  Concreat version of template class for PyRosetta
/// @author Sergey Lyskov

#ifndef INCLUDED_core_pose_PyDatacacheCacheableObserver_hh
#define INCLUDED_core_pose_PyDatacacheCacheableObserver_hh

// unit headers
#include <core/pose/datacache/CacheableObserver.hh>
#include <basic/datacache/DataCache.hh>


namespace core {
namespace pose {

// PyRosetta concreate class
class Py_basic_datacache_DataCache_CacheableObserver : public basic::datacache::DataCache< datacache::CacheableObserver > {};


} // namespace pose
} // namespace core


#endif /* INCLUDED_core_pose_PyDatacacheCacheableObserver_hh */
