// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AlignEvaluator
/// @brief base class for evaluators that evaluate based on some sort of sequence
/// alignment.
/// @author James Thompson

#ifndef INCLUDED_protocols_evaluation_AlignEvaluator_hh
#define INCLUDED_protocols_evaluation_AlignEvaluator_hh

#include <core/sequence/SequenceAlignment.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace evaluation {

class AlignEvaluator : public SingleValuePoseEvaluator< core::Real > {

public:
	AlignEvaluator(
		core::pose::PoseCOP native_pose,
		std::string tag = "",
		bool report_aln_components = true,
		core::sequence::SequenceAlignmentOP aln = nullptr
	);

	~AlignEvaluator() override;

	void apply(
		core::pose::Pose & pose,
		std::string tag,
		core::io::silent::SilentStruct & ss
	) const override = 0;

	/// @brief outdated method - don't use!
	core::Real apply(
		core::pose::Pose & /*pose*/
	) const override {
		utility_exit_with_message(
			"Called AlignEvaluator::apply( Pose & pose ). Don't do that!\n"
		);
		return 0.0;
	}

	core::sequence::SequenceAlignmentOP get_alignment(
		core::pose::Pose const & pose
	) const;

	void report_aln_components( bool const setting );

	bool report_aln_components() const;

	std::string tag() const;

	core::pose::PoseCOP native_pose() const;

private:
	core::pose::PoseCOP native_pose_;
	std::string tag_;
	bool report_aln_components_;
	core::sequence::SequenceAlignmentOP aln_;
}; // AlignEvaluator

} // evaluation
} // protocols

#endif
