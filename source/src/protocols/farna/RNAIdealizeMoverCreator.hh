// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/farna/RNAIdealizeMoverCreator.hh
/// @brief Slowly accomodate movement from non-ideal to ideal bond lengths and angles by repeated minimization
/// @author Andy Watkins (amw579@nyu.edu)

#ifndef INCLUDED_protocols_farna_RNAIdealizeMoverCreator_hh
#define INCLUDED_protocols_farna_RNAIdealizeMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace farna {

class RNAIdealizeMoverCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP
	create_mover() const;

	virtual std::string
	keyname() const;

};

} //protocols
} //farna

#endif //INCLUDED_protocols/farna_RNAIdealizeMover_fwd_hh
