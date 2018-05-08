// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/sewing/hashing/EdgeMapGenerator.fwd.hh
/// @brief Class used to generate edge files for use with SMART_SEWING devel.
/// @author guffysl (guffy@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_hashing_EdgeMapGenerator_fwd_hh
#define INCLUDED_protocols_sewing_hashing_EdgeMapGenerator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace sewing {
namespace hashing {

class EdgeMapGenerator;

typedef utility::pointer::shared_ptr< EdgeMapGenerator > EdgeMapGeneratorOP;
typedef utility::pointer::shared_ptr< EdgeMapGenerator const > EdgeMapGeneratorCOP;


} //hashing
} //sewing
} //protocols


#endif //INCLUDED_protocols_sewing_hashing_EdgeMapGenerator_fwd_hh





