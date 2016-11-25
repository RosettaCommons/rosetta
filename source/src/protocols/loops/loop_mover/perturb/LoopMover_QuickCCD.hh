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

#ifndef INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_QuickCCD_hh
#define INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_QuickCCD_hh

// package headers
#include <protocols/loops/loop_mover/perturb/LoopMover_QuickCCD.fwd.hh>
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>
#include <protocols/moves/Mover.hh>

// project headers
#include <core/types.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


// C++ Headers


namespace protocols {
namespace loops {
namespace loop_mover {
namespace perturb {

/// @brief LoopMover utilizing fragment insertion, ccd loop closure, and
/// minimization
class LoopMover_Perturb_QuickCCD: public IndependentLoopMover {

public: // construct/destruct

	/// @brief Loops constructor
	/// @remarks Will be initialized with centroid level score function 'score4L'.
	LoopMover_Perturb_QuickCCD();


	/// @brief Loops constructor
	/// @param[in] loops_in the set of loops to model
	/// @remarks Will be initialized with centroid level score function 'score4L'.
	LoopMover_Perturb_QuickCCD(
		protocols::loops::LoopsOP loops_in
	);

	/// @brief Loops & ScoreFunction constructor
	/// @param[in] loops_in the set of loops to model
	/// @param[in] scorefxn desired ScoreFunction
	LoopMover_Perturb_QuickCCD(
		protocols::loops::LoopsOP loops_in,
		core::scoring::ScoreFunctionOP scorefxn
	);

	/// @brief Loops & ScoreFunction constructor
	/// @param[in] loops_in the set of loops to model
	/// @param[in] scorefxn desired ScoreFunction
	/// @param[in] fragset is the FragSet to be used
	/// line?
	LoopMover_Perturb_QuickCCD(
		protocols::loops::LoopsOP loops_in,
		core::scoring::ScoreFunctionOP scorefxn,
		core::fragment::FragSetOP fragset
	);

	//destructor
	~LoopMover_Perturb_QuickCCD();

	// XRW TEMP  virtual std::string get_name() const;

public: // virtual constructors

	/// @brief clone this object
	protocols::moves::MoverOP clone() const override;

public: // accessors

	/// @brief randomize loops prior to loop modeling?
	inline
	bool randomize_loop() const {
		return randomize_loop_;
	}


public: // mutators


	/// @brief indicate whether loops should be randomized prior to modeling
	inline
	void randomize_loop( bool const flag ) {
		randomize_loop_ = flag;
	}


	/// @brief set default settings
	/// @details default settings are as follows:
	///  <ul>
	///      <li> randomize_loop() = true
	///  </ul>
	void set_default_settings() {
		randomize_loop_ = true;
	}

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected: // virtual loop operations

	/// @brief loop modeling protocol implementation
	LoopResult model_loop(
		core::pose::Pose & pose,
		protocols::loops::Loop const & loop
	) override;

	basic::Tracer & tr() const override;

protected: // data. should be private!

	/// @brief randomize loops prior to performing loop modeling?
	bool randomize_loop_;
}; // class LoopMover_Perturb_QuickCCD


void fast_ccd_close_loops(
	core::pose::Pose & pose,
	protocols::loops::Loop const & loop,
	core::kinematics::MoveMap & mm
);

} //namespace perturb
} //namespace loop_mover
} //namespace loops
} //namespace protocols

#endif //INCLUDED_protocols_loops_loop_mover_perturb_LoopMover_QuickCCD_hh
