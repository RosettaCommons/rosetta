// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/swa/rna/RNA_TorsionScreener.hh
/// @brief Screener checking whether the rna torsions are resonable
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_sampler_screener_RNA_TorsionScreener_HH
#define INCLUDED_protocols_sampler_screener_RNA_TorsionScreener_HH

// Unit headers
#include <protocols/stepwise/sampler/screener/RNA_TorsionScreener.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/rna/RNA_SuiteName.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace stepwise {
namespace sampler {
namespace screener {

class RNA_TorsionScreener : public utility::pointer::ReferenceCount {
public:
	RNA_TorsionScreener();

	/// @brief screen the pose
	bool screen( core::pose::Pose const & pose, core::Size const suite ) ;

private:
	core::Real angle_conv( core::Real const input );

	core::pose::rna::RNA_SuiteNameOP suitename_;
};

} //screener
} //sampler
} //stepwise
} //protocols
#endif
