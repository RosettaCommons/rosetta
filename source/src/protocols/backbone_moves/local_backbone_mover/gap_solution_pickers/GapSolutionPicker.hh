// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backbone_moves/local_backbone_mover/gap_solution_pickers/GapSolutionPicker.hh
/// @brief Base class for picking solution for a gap.
/// @author xingjiepan (xingjiepan@gmail.com)


#ifndef INCLUDED_protocols_backbone_moves_local_backbone_mover_gap_solution_pickers_GapSolutionPicker_hh
#define INCLUDED_protocols_backbone_moves_local_backbone_mover_gap_solution_pickers_GapSolutionPicker_hh

#include <protocols/backbone_moves/local_backbone_mover/types.hh>
#include <protocols/backbone_moves/local_backbone_mover/FreePeptide.fwd.hh>
#include <protocols/backbone_moves/local_backbone_mover/gap_solution_pickers/GapSolutionPicker.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// BOOST
#include <boost/noncopyable.hpp>

namespace protocols {
namespace backbone_moves {
namespace local_backbone_mover {
namespace gap_solution_pickers {

///@brief Base class for picking solution for a gap.
class GapSolutionPicker :
	public utility::pointer::ReferenceCount, protected boost::noncopyable {

public:

	virtual Size pick(core::pose::Pose const & pose, FreePeptide const & free_peptide,
		vector1<vector1<Real> > const & pivot_torsions, Size const pivot) const = 0;

private:

};


} //protocols
} //backbone_moves
} //local_backbone_mover
} //gap_solution_pickers



#endif //INCLUDED_protocols_backbone_moves_local_backbone_mover_gap_solution_pickers_GapSolutionPicker_hh





