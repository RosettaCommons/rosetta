// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/GlycanRelaxMoverCreator.hh
/// @brief Main mover for Glycan Relax, which optimizes glycans in a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_carbohydrates_GlycanRelaxMoverCreator_hh
#define INCLUDED_protocols_carbohydrates_GlycanRelaxMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>


namespace protocols {
namespace carbohydrates {

class GlycanRelaxMoverCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();



};

} //protocols
} //carbohydrates

#endif //INCLUDED_protocols/carbohydrates_GlycanRelaxMover_fwd_hh



