// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/residue_datacache.hh
/// @brief  enum for storing data in the Residue's datacache
/// @author pbradley@fhcrc.org

#ifndef INCLUDED_core_conformation_residue_datacache_hh
#define INCLUDED_core_conformation_residue_datacache_hh

// unit headers

// #include <core/conformation/residue_datacache.fwd.hh>

// // package headers
// #include <core/conformation/Residue.fwd.hh>

// #include <basic/datacache/CacheableData.hh>
// #include <basic/datacache/DataCache.hh>


namespace core {
namespace conformation {
namespace residue_datacache {

enum ResidueDataCacheIndex {
	LK_BALL_INFO = 1,
	n_cacheable_types = LK_BALL_INFO
};


} // namespace residue_datacache
} // namespace conformation
} // namespace core


#endif /* INCLUDED_basic_datacache_ResidueDataCache_HH */
