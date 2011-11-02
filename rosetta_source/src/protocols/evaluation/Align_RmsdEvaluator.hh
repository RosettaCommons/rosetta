// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Align_RmsdEvaluator
/// @author James Thompson

#ifndef INCLUDED_protocols_evaluation_Align_RmsdEvaluator_hh
#define INCLUDED_protocols_evaluation_Align_RmsdEvaluator_hh

#include <protocols/evaluation/AlignEvaluator.hh>
// AUTO-REMOVED #include <protocols/evaluation/util.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/types.hh>

#include <ObjexxFCL/FArray2D.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace evaluation {

class Align_RmsdEvaluator : public AlignEvaluator {

public:
	Align_RmsdEvaluator(
		core::pose::PoseCOP native_pose,
		std::string tag = "",
		bool calc_gdt = true,
		core::sequence::SequenceAlignmentOP aln = 0
	);

	~Align_RmsdEvaluator();

	virtual void apply(
		core::pose::Pose & pose,
		std::string tag,
		core::io::silent::SilentStruct & ss
	) const;

	virtual core::Real apply(
		core::pose::Pose & /*pose*/
	) const {
		return 0;
	}

	void report_gdt_components( bool const setting ) {
		report_gdt_components_ = setting;
		if ( setting ) calc_gdt(true);
	}

	bool report_gdt_components() const {
		return report_gdt_components_;
	}

	void calc_gdt( bool const setting ) {
		calc_gdt_ = setting;
	}

	bool calc_gdt() const {
		return calc_gdt_;
	}

private:
	bool calc_gdt_;
	bool report_gdt_components_;
}; // Align_RmsdEvaluator

} // evaluation
} // protocols

#endif
