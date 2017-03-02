// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/free_peptide_movers/TranslationFreePeptideMover.hh
/// @brief Translate the free peptide by a cartesian vector.
/// @author xingjiepan (xingjiepan@gmail.com)


#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_TranslationFreePeptideMover_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_TranslationFreePeptideMover_hh

#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/TranslationFreePeptideMover.fwd.hh>
#include <protocols/backbone_moves/local_backbone_mover/free_peptide_movers/FreePeptideMover.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>

namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {
namespace free_peptide_movers {

///@brief Translate the free peptide by a cartesian vector.
class TranslationFreePeptideMover : public FreePeptideMover {

public:

	///@brief Translate the free peptide by a given vector
	TranslationFreePeptideMover(xyzVector <Real> v_translate);

	///@brief Translate the free peptide by generating a random vector
	/// whose amplitude is smaller than a given value
	TranslationFreePeptideMover(Real max_amplitude);

	virtual ~TranslationFreePeptideMover();

	virtual void apply(FreePeptide &free_peptide);

private:

	bool random_ = false;
	Real max_amplitude_ = 0;
	xyzVector <Real> v_translate_;

};


} //protocols
} //backbone_moves
} //local_backbone_mover
} //free_peptide_movers



#endif //INCLUDED_protocols_backbone_moves_local_backbone_mover_free_peptide_movers_TranslationFreePeptideMover_hh





