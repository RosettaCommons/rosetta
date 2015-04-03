// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /devel/DenovoProteinDesign/SequenceFitnessMover.hh
/// @brief
/// @author

#ifndef INCLUDED_devel_denovo_protein_design_SequenceFitnessMover_hh
#define INCLUDED_devel_denovo_protein_design_SequenceFitnessMover_hh

// Unit Headers
#include <devel/denovo_protein_design/SequenceFitnessMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <utility/vector1.hh>


// Utility Headers

namespace devel {
namespace denovo_protein_design {

/// @details
class SequenceFitnessMover : public protocols::moves::Mover {

public:

	/// @brief
	SequenceFitnessMover(
	);

	virtual ~SequenceFitnessMover();

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

private:

};//end SequenceFitnessMover

}//namespace denovo_protein_design
}//namespace devel

#endif // INCLUDED_devel_DenovoProteinDesign_SequenceFitnessMover_HH
