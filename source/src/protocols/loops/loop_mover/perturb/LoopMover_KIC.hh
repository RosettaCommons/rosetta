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
/// @author Daniel J. Mandell
/// @author Amelie Stein

#ifndef INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_KIC_hh
#define INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_KIC_hh

#include <protocols/loops/loop_mover/perturb/LoopMover_KIC.fwd.hh>
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>


// C++ Headers


///////////////////////////////////////////////////////////////////////////////
namespace protocols {
namespace loops {
namespace loop_mover {
namespace perturb {

class LoopMover_Perturb_KIC: public loop_mover::IndependentLoopMover {
public:

	LoopMover_Perturb_KIC();

	LoopMover_Perturb_KIC(
		protocols::loops::LoopsOP loops_in
	);

	LoopMover_Perturb_KIC(
		protocols::loops::LoopsOP  loops_in,
		core::scoring::ScoreFunctionOP  scorefxn
	);

	//destructor
	~LoopMover_Perturb_KIC();

	// XRW TEMP  virtual std::string get_name() const;

	void show(std::ostream & output=std::cout) const override;

	void set_extended_torsions(
		core::pose::Pose & pose,
		Loop const & loop
	) override;

	void set_default_settings();


	core::Size get_max_kic_build_attempts() const
	{
		return max_kic_build_attempts_;
	}

	void set_max_kic_build_attempts( core::Size max_kic_build_attempts )
	{
		max_kic_build_attempts_ = max_kic_build_attempts;
	}


	/// @brief Clone this object
	protocols::moves::MoverOP clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );



private:

	core::Size max_seglen_; // maximum KIC segment length
	bool recover_low_;
	core::Size max_kic_build_attempts_;
	core::Size remodel_kic_attempts_;

protected:

	loop_mover::LoopResult model_loop(
		core::pose::Pose & pose,
		protocols::loops::Loop const & loop
	) override;

	virtual basic::Tracer & tr() const override;
};

std::ostream &operator<< ( std::ostream &os, LoopMover_Perturb_KIC const &mover );

} //namespace perturb
} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_KIC_HH
