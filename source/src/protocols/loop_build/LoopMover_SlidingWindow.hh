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

#ifndef INCLUDED_protocols_loop_build_LoopMover_SlidingWindow_hh
#define INCLUDED_protocols_loop_build_LoopMover_SlidingWindow_hh

// package headers
#include <protocols/loops/loop_mover/IndependentLoopMover.hh>
#include <protocols/moves/Mover.hh>

// project headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <utility/vector1.hh>


// C++ Headers


namespace protocols {
namespace loop_build {


/// @brief LoopMover utilizing fragment insertion, ccd loop closure, and minimization
class LoopMover_SlidingWindow: public loops::loop_mover::IndependentLoopMover {


public: // construct/destruct


	/// @brief Loops constructor
	/// @remarks Will be initialized with centroid level score function 'score4L'.
	LoopMover_SlidingWindow();

	/// @brief Loops constructor
	/// @param[in] loops_in the set of loops to model
	/// @param[in] frags_from_file read fragments from files specified on command line?
	/// @remarks Will be initialized with centroid level score function 'score4L'.
	LoopMover_SlidingWindow(
		protocols::loops::LoopsOP loops_in
	);


	/// @brief Loops & ScoreFunction constructor
	/// @param[in] loops_in the set of loops to model
	/// @param[in] scorefxn desired ScoreFunction
	/// @param[in] frags_from_file read fragments from files specified on command line?
	LoopMover_SlidingWindow(
		protocols::loops::LoopsOP loops_in,
		core::scoring::ScoreFunctionOP scorefxn
	);


public: // virtual constructors


	/// @brief clone this object

	protocols::moves::MoverOP clone() const override {
		return protocols::moves::MoverOP( new LoopMover_SlidingWindow( *this ) );
	}


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
	loops::loop_mover::LoopResult model_loop(
		core::pose::Pose & pose,
		protocols::loops::Loop const & loop
	) override;

	basic::Tracer & tr() const override;

protected: // data


	/// @brief list of fragment libraries to use
	std::vector< core::fragment::ConstantLengthFragSetOP > frag_libs_;


	/// @brief randomize loops prior to performing loop modeling?
	bool randomize_loop_;


};

/*  Undefined, commenting out to fix PyRosetta build
void fast_ccd_close_loops(
core::pose::Pose & pose,
protocols::loops::Loop const & loop,
core::kinematics::MoveMap & mm
); */


} //namespace loop_build
} //namespace protocols

#endif //INCLUDED_protocols_loop_build_LoopMover_QuickCCD_HH
