// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_CCD_hh
#define INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_CCD_hh


#include <protocols/loops/loop_mover/perturb/LoopMover_CCD.fwd.hh>
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>

///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {
namespace loop_mover {
namespace perturb {

////////////////////////////////////////////////////////////////////////////////////////
/// @details the main function to model one single loop in centroid mode. The
/// modeling algorithm is fragment_ccd_min_trial, which consists of perturbing
/// the loop conformation by fragment insertion , then close the loop by CCD loop
/// closure, then minimize the loop conformation and finally subject it to Monte
/// Carlo acceptance or rejection. The loop conformation will be initialized as
/// extended conformation if it is specified in the loop definition, resembling
/// ab initio loop modeling in real practice. The loop has to be long enough for
/// inserting certain length of fragments.
/////////////////////////////////////////////////////////////////////////////////////////
class LoopMover_Perturb_CCD: public loop_mover::IndependentLoopMover {
public:
	//constructor
	LoopMover_Perturb_CCD();

	LoopMover_Perturb_CCD(
		protocols::loops::LoopsOP  loops_in
	);

	LoopMover_Perturb_CCD(
		protocols::loops::LoopsOP  loops_in,
		core::scoring::ScoreFunctionOP  scorefxn
	);

	LoopMover_Perturb_CCD(
		protocols::loops::LoopsOP  loops_in,
		core::scoring::ScoreFunctionOP  scorefxn,
		core::fragment::FragSetOP fragset
	);

	//destructor
	~LoopMover_Perturb_CCD();

	// XRW TEMP  virtual std::string get_name() const;

	void show(std::ostream & output=std::cout) const override;

	/// @brief Clone this object
	protocols::moves::MoverOP clone() const override;

	void set_default_settings(){}
	virtual protocols::loops::LoopsCOP get_loops() const;
	virtual core::scoring::ScoreFunctionOP get_scorefxn() const;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



protected:
	std::vector< core::fragment::FragSetOP > frag_libs_;

	loop_mover::LoopResult model_loop( core::pose::Pose & pose,
		protocols::loops::Loop const & loop ) override;

	basic::Tracer & tr() const override;
};

std::ostream &operator<< ( std::ostream &os, LoopMover_Perturb_CCD const &mover );

} //namespace perturb
} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_CCD_HH
