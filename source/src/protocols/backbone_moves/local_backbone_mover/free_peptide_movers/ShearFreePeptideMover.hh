// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/free_peptide_movers/ShearFreePeptideMover.hh
/// @brief Change the psi torsion of a residue and phi torsion of its subsequent residue by a opposite value.
/// @author xingjiepan (xingjiepan@gmail.com)


#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_ShearFreePeptideMover_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_ShearFreePeptideMover_hh

#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/ShearFreePeptideMover.fwd.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/FreePeptideMover.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {
namespace free_peptide_movers {

///@brief Change the psi torsion of a residue and phi torsion of its subsequent residue by a opposite value.
class ShearFreePeptideMover : public FreePeptideMover {

public:

	ShearFreePeptideMover(Size seqpos, Real torsion_degree, bool random = false);

	virtual ~ShearFreePeptideMover();

	virtual void apply(FreePeptide &free_peptide);

private:

	bool random_ = false;
	Size seqpos_ = 0;
	Real torsion_degree_ = 0;
};


} //protocols
} //backbone_moves
} //local_backbone_mover
} //free_peptide_movers



#endif //INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_ShearFreePeptideMover_hh





