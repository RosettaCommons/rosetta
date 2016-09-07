// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /devel/CreateStartingStructureMover/CreateStartingStructureMover.hh
/// @brief
/// @author

#ifndef INCLUDED_devel_denovo_protein_design_CreateStartingStructureMover_hh
#define INCLUDED_devel_denovo_protein_design_CreateStartingStructureMover_hh

// Unit Headers
#include <devel/denovo_protein_design/CreateStartingStructureMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


// Utility Headers

namespace devel {
namespace denovo_protein_design {

/// @details
class CreateStartingStructureMover : public protocols::moves::Mover {

public:

	/// @brief
	CreateStartingStructureMover();

	~CreateStartingStructureMover() override;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override;

private:

};//end CreateStartingStructureMover

}//namespace denovo_protein_design
}//namespace devel

#endif // INCLUDED_devel_CreateStartingStructureMover_CreateStartingStructureMover_HH
