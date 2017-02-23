// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/free_peptide_movers/CircularPermuteFreePeptideMover.hh
/// @brief Circularly permute the torsions on the free peptide.
/// @author xingjiepan (xingjiepan@gmail.com)


#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_CircularPermuteFreePeptideMover_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_CircularPermuteFreePeptideMover_hh

#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/CircularPermuteFreePeptideMover.fwd.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/FreePeptideMover.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {
namespace free_peptide_movers {

///@brief Circularly permute the torsions on the free peptide.
class CircularPermuteFreePeptideMover : public FreePeptideMover {

public:

	CircularPermuteFreePeptideMover(Size off_set, bool direction = true);

	virtual ~CircularPermuteFreePeptideMover();

	virtual void apply(FreePeptide &free_peptide);

private:

	Size off_set_ = 0;

	// If direction is true, the torsions of i + 1 th residue will
	// be given to the i th residue. Otherwise the opposite direction.
	bool direction_ = true;
};


} //protocols
} //backbone_moves
} //local_backbone_mover
} //free_peptide_movers



#endif //INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_CircularPermuteFreePeptideMover_hh





