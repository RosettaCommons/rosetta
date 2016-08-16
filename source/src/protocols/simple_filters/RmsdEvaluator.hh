// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_simple_filters_RmsdEvaluator_hh
#define INCLUDED_protocols_simple_filters_RmsdEvaluator_hh


// Unit Headers
#include <protocols/evaluation/PoseEvaluator.hh>

// Package Headers

// Project Headers
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <list>

#include <core/scoring/rms_util.hh>
#include <protocols/loops/Loops.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {


class RmsdEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	RmsdEvaluator( core::pose::PoseCOP, core::Size start, core::Size end, std::string tag = "", bool bGDT = true );
	RmsdEvaluator( core::pose::PoseCOP, std::string tag = "", bool bGDT = true);
	~RmsdEvaluator();
	virtual void
	apply ( core::pose::Pose& pose, std::string tag, core::io::silent::SilentStruct &pss) const;

	/// @brief evaluate pose
	virtual core::Real apply( core::pose::Pose& ) const;

	void report_gdt_components( bool const setting ){ report_gdt_components_ = setting; }

private:
	core::pose::PoseCOP rmsd_pose_;
	core::Size start_;
	core::Size end_;
	bool bGDT_;
	std::string tag_;
	bool report_gdt_components_;
};

class SelectRmsdEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	SelectRmsdEvaluator( core::pose::PoseCOP, core::scoring::ResidueSelection const& selection, std::string tag = "", bool CAonly=true  );
	SelectRmsdEvaluator( core::pose::PoseCOP, utility::vector1< core::Size> const& selection, std::string tag = "", bool CAonly=true  );

	//work it out by yourself from missing density == whacky random coords
	SelectRmsdEvaluator( core::pose::PoseCOP, std::string tag = "", bool CAonly=true   );

	//work it out by yourself from missing density == whacky random coords
	SelectRmsdEvaluator( core::pose::Pose const&, std::string tag = "", bool CAonly=true );

	/// @brief evaluate pose
	virtual core::Real apply( core::pose::Pose& ) const;

private:
	core::pose::PoseCOP rmsd_pose_;
	core::scoring::ResidueSelection selection_;
	std::string tag_;
	bool CAonly_;
};

class SelectGdtEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	SelectGdtEvaluator( core::pose::PoseCOP, core::scoring::ResidueSelection const& selection, std::string tag = "" );
	SelectGdtEvaluator( core::pose::PoseCOP, utility::vector1< core::Size> const& selection, std::string tag = "" );

	//work it out by yourself from missing density == whacky random coords
	SelectGdtEvaluator( core::pose::PoseCOP, std::string tag = "" );
	~SelectGdtEvaluator();
	//work it out by yourself from missing density == whacky random coords
	SelectGdtEvaluator( core::pose::Pose const&, std::string tag = "" );

	/// @brief evaluate pose
	virtual core::Real apply( core::pose::Pose& ) const;

private:
	core::pose::PoseCOP rmsd_pose_;
	core::scoring::ResidueSelection selection_;
	std::string tag_;
};

class SelectMaxsubEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	SelectMaxsubEvaluator( core::pose::PoseCOP, core::scoring::ResidueSelection const& selection, std::string tag = "", core::Real rmsd_threshold = 4.0  );
	SelectMaxsubEvaluator( core::pose::PoseCOP, utility::vector1< core::Size> const& selection, std::string tag = "", core::Real rmsd_threshold = 4.0  );

	//work it out by yourself from missing density == whacky random coords
	SelectMaxsubEvaluator( core::pose::PoseCOP, std::string tag = "", core::Real rmsd_threshold = 4.0  );

	//work it out by yourself from missing density == whacky random coords
	SelectMaxsubEvaluator( core::pose::Pose const&, std::string tag = "", core::Real rmsd_threshold = 4.0 );

	/// @brief evaluate pose
	virtual core::Real apply( core::pose::Pose& ) const;

private:
	core::pose::PoseCOP rmsd_pose_;
	core::scoring::ResidueSelection selection_;
	std::string tag_;
	core::Real rmsd_threshold_;
};

class SymmetricRmsdEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	SymmetricRmsdEvaluator( core::pose::PoseCOP, std::string tag );
	~SymmetricRmsdEvaluator();
	virtual core::Real apply ( core::pose::Pose& ) const;

private:
	core::pose::PoseCOP rmsd_pose_;

};

class LoopRmsdEvaluator : public evaluation::SingleValuePoseEvaluator< core::Real > {
public:
	LoopRmsdEvaluator( core::pose::PoseCOP, protocols::loops::Loops, std::string tag, bool CA_only, bool superimpose );
	LoopRmsdEvaluator( core::pose::PoseCOP, protocols::loops::Loops, protocols::loops::Loops core, std::string tag, bool CA_only, bool superimpose );
	virtual core::Real apply ( core::pose::Pose& ) const;

private:
	core::pose::PoseCOP rmsd_pose_;
	protocols::loops::Loops loops_;
	protocols::loops::Loops core_;
	bool CAonly_;
	bool superimpose_;
};

}
}

#endif
