// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/free_peptide_movers/LongAxisRotationFreePeptideMover.hh
/// @brief Rotate the free peptide along the axis formed by its first and last CA atom
/// @author xingjiepan (xingjiepan@gmail.com)


#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_LongAxisRotationFreePeptideMover_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_LongAxisRotationFreePeptideMover_hh

#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/LongAxisRotationFreePeptideMover.fwd.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/FreePeptideMover.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {
namespace free_peptide_movers {

///@brief Rotate the free peptide along the axis formed by its first and last CA atom
class LongAxisRotationFreePeptideMover : public FreePeptideMover {

public:

	// If random is set to false, the rotation redian will
	// be the given value. If random is set to true, the
	// rotation redian will be a random number between
	// radian and -radian.
	LongAxisRotationFreePeptideMover(Real radian, bool random=false);

	virtual ~LongAxisRotationFreePeptideMover();

	virtual void apply(FreePeptide &free_peptide);

private:

	Real radian_ = 0;
	bool random_ = false;
};


} //protocols
} //backbone_moves
} //local_backbone_mover
} //free_peptide_movers



#endif //INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_LongAxisRotationFreePeptideMover_hh





