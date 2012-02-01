// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   
/// @brief  
/// @author Javier Castellanos ( javiercv@uw.edu )

#ifndef INCLUDED_devel_constrained_sequence_design_SequenceConstraintLoader_hh
#define INCLUDED_devel_constrained_sequence_design_SequenceConstraintLoader_hh

// Package Headers
#include <devel/constrained_sequence_design/SequenceConstraint.hh>
#include <devel/constrained_sequence_design/SequenceConstraintSet.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/jd2/parser/DataLoader.hh>
#include <protocols/moves/DataMap.fwd.hh>


// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>


namespace devel {
namespace constrained_sequence_design {

class SequenceConstraintLoader : public protocols::jd2::parser::DataLoader {

public:
	SequenceConstraintLoader();
	virtual ~SequenceConstraintLoader();

	virtual
	void load_data(
		core::pose::Pose const & pose,
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data
	) const;

};

} // namespace devel
} // namespace constrained_sequence_design

#endif
