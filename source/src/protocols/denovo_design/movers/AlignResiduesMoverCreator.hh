// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/AlignResiduesMoverCreator.hh
/// @brief Aligns one residue onto another
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_movers_AlignResiduesMoverCreator_hh
#define INCLUDED_protocols_denovo_design_movers_AlignResiduesMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace denovo_design {
namespace movers {

class AlignResiduesMoverCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP
	create_mover() const;

	virtual std::string
	keyname() const;

};

} //protocols
} //denovo_design
} //movers

#endif //INCLUDED_protocols/denovo_design/movers_AlignResiduesMover_fwd_hh
