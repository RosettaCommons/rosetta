// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/RNA_ChainClosableGeometryStubBasedScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_RNA_ChainClosableGeometryStubBasedScreener_HH
#define INCLUDED_protocols_stepwise_screener_RNA_ChainClosableGeometryStubBasedScreener_HH

#include <protocols/stepwise/screener/StepWiseResiduePairScreener.hh>
#include <protocols/stepwise/screener/RNA_ChainClosableGeometryStubBasedScreener.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace stepwise {
namespace screener {

class RNA_ChainClosableGeometryStubBasedScreener: public StepWiseResiduePairScreener {

public:

	//constructor
	RNA_ChainClosableGeometryStubBasedScreener( modeler::rna::checker::RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_checker,
		utility::vector1< core::pose::PoseOP > screening_pose_list,
		core::kinematics::Stub const & moving_res_base_stub,
		Size const reference_res,
		utility::vector1< core::conformation::ResidueOP > moving_rsd_at_origin_list,
		bool const using_predefined_moving_rsd_list = true );


	//constructor
	RNA_ChainClosableGeometryStubBasedScreener( modeler::rna::checker::RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_checker,
		utility::vector1< core::pose::PoseOP > screening_pose_list,
		core::kinematics::Stub const & moving_res_base_stub,
		Size const chain_break_reference_res );

	//destructor
	~RNA_ChainClosableGeometryStubBasedScreener();

public:

	void
	get_update( sampler::StepWiseSamplerBaseOP sampler );

	bool
	check_screen();

	std::string
	name() const { return "RNA_ChainClosableGeometryStubBasedScreener"; }

	StepWiseScreenerType
	type() const { return RNA_CHAIN_CLOSABLE_GEOMETRY_STUB_BASED; }

private:

	modeler::rna::checker::RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_checker_;

	utility::vector1< core::pose::PoseOP > screening_pose_list_;
	core::conformation::ResidueOP moving_rsd_at_origin_;
	utility::vector1< core::conformation::ResidueOP > moving_rsd_at_origin_list_;
	core::kinematics::Stub const & moving_res_base_stub_;
	Size const reference_res_;

	bool const using_predefined_moving_rsd_list_;

};

} //screener
} //stepwise
} //protocols

#endif
