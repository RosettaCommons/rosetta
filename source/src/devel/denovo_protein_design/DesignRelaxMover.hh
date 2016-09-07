// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /protocols/moves/MoverStub.hh
/// @brief
/// @author

#ifndef INCLUDED_devel_denovo_protein_design_DesignRelaxMover_hh
#define INCLUDED_devel_denovo_protein_design_DesignRelaxMover_hh

// Unit Headers
#include <devel/denovo_protein_design/DesignRelaxMover.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility Headers

#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>


namespace devel {
namespace denovo_protein_design {

/// @details
class DesignRelaxMover : public protocols::moves::Mover {

public:

	/// @brief

	// default constructor will set designtask from commanline and resfile
	// with design using softrep_design_wts+patch and relax using score12+patch
	DesignRelaxMover();

	// ctor with hand made designtask
	DesignRelaxMover(
		core::pack::task::TaskFactoryOP designtaskfactory
	);

	// ctor with hand made designtask and control of scorefunctions
	DesignRelaxMover(
		core::pack::task::TaskFactoryOP designtaskfactory,
		core::scoring::ScoreFunctionOP designfxn,
		core::scoring::ScoreFunctionOP relaxfxn
	);

	~DesignRelaxMover() override;

	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override;

private:
	core::pack::task::TaskFactoryOP designtaskfactory_;
	core::scoring::ScoreFunctionOP designfxn_;
	core::scoring::ScoreFunctionOP relaxfxn_;

};//end DesignRelaxMover

}//namespace denovo_protein_design
}//namespace devel

#endif // INCLUDED_devel_DenovoProteinDesign_DesignRelaxMover_HH
