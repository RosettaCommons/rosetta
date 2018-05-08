// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/hashing/Hasher.fwd.hh
/// @brief A geometric hashing class used by the SEWING protocol
/// @author Minnie Langlois (minnie@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_hashing_Hasher_fwd_hh
#define INCLUDED_protocols_sewing_hashing_Hasher_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>



// Forward
namespace protocols {
namespace sewing {
namespace hashing {

class Hasher;
struct HasherSettings;
typedef utility::pointer::shared_ptr< Hasher > HasherOP;
typedef utility::pointer::shared_ptr< Hasher const > HasherCOP;



} //hashing
} //sewing
} //protocols


#endif //INCLUDED_protocols_sewing_hashing_Hasher_fwd_hh





