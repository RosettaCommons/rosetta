// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AlignEvaluator
/// @author James Thompson

#include <protocols/evaluation/AlignEvaluator.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>

// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>

#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/model_quality/rms.hh>

#include <core/chemical/ResidueType.hh>
#include <core/id/types.hh>
#include <utility/vector1.hh>


// C++ headers

static thread_local basic::Tracer tr( "protocols.evaluation.AlignEvaluator" );

namespace protocols {
namespace evaluation {

AlignEvaluator::~AlignEvaluator() {}

AlignEvaluator::AlignEvaluator(
	core::pose::PoseCOP native_pose,
	std::string tag,
	bool report_aln_components,
	core::sequence::SequenceAlignmentOP aln
) :
	SingleValuePoseEvaluator< core::Real >("align_rms" + tag),
	native_pose_(native_pose),
	tag_(tag),
	report_aln_components_(report_aln_components),
	aln_(aln)
{}

core::sequence::SequenceAlignmentOP AlignEvaluator::get_alignment(
	core::pose::Pose const & pose
) const {
	if ( aln_ ) return aln_;
	// calculate a naive alignment if we don't already have an alignment.
	using namespace core::sequence;
	SequenceOP model_seq( new Sequence( pose.sequence(),  "model",  1 ) );
	SequenceOP native_seq( new Sequence( native_pose_->sequence(), "native", 1 ) );
	tr.Debug << "aligning: " << std::endl;
	tr.Debug << model_seq ->to_string() << std::endl;
	tr.Debug << native_seq->to_string() << std::endl;
	SequenceAlignmentOP aln_op( new SequenceAlignment );
	*aln_op = align_naive(model_seq,native_seq);
	tr.Debug << "alignment: " << std::endl;
	tr.Debug << *aln_op;
	tr.flush_all_channels();
	return aln_op;
}


void AlignEvaluator::report_aln_components( bool const setting ) {
	report_aln_components_ = setting;
}

bool AlignEvaluator::report_aln_components() const {
	return report_aln_components_;
}

std::string AlignEvaluator::tag() const {
	return tag_;
}

core::pose::PoseCOP AlignEvaluator::native_pose() const {
	return native_pose_;
}

} // evaluation
} // protocols
