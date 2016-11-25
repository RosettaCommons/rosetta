// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_loop_modeling_LoopBuilder_HH
#define INCLUDED_protocols_loop_modeling_LoopBuilder_HH

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopBuilder.fwd.hh>
#include <protocols/loop_modeling/LoopModelerTests.fwd.hh> //For friendship.

// Core headers
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/kinematic_closure/KicMover.fwd.hh>
#include <protocols/loop_modeling/refiners/MinimizationRefiner.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/kinematic_closure/perturbers/Perturber.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace loop_modeling {

/// @brief Build loops from scratch.
/// @author Kale Kundert
/// @author Roland A. Pache, PhD
///
/// @details Building a loop from scratch is useful in two scenarios.  The
/// first is when there's missing density that needs to be modeled, and the
/// second is when the whole loop modeling algorithm needs to be benchmarked.
/// This mover uses kinematic closure (KIC) to build loops.  By default, the
/// loop are built by picking phi and psi values from a Ramachandran
/// distribution and setting all other DOFs to ideal values.  Phi and psi
/// values can also be picked using fragment libraries.  Loop building succeeds
/// when a model is found that passes a more-lenient-than-usual bump check.  If
/// no such model is found after 1000 iterations, the mover gives up and
/// reports failure.
///
/// This process can be very slow for long loops, because there's nothing
/// guiding the algorithm towards the right solution.  By default, torsions are
/// just being randomly picked, and they often won't fit in the relatively
/// narrow space that's available.  The problem is worse for interior loops
/// than it is for surface loops, of course.  This algorithm seems to work well
/// enough on 12 residue loops, but beyond that it may be necessary to develop
/// a smarter algorithm that preferentially builds into free space.

class LoopBuilder : public LoopMover {

	friend class ::LoopModelerTests;

public:

	/// @brief Default constructor.
	LoopBuilder();

	/// @brief Default destructor.
	~LoopBuilder() override;

	/// @copydoc LoopMover::get_name
	// XRW TEMP  string get_name() const override { return "LoopBuilder"; }

public:

	/// @brief Configure from a RosettaScripts tag.
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		Pose const & pose) override;

public:

	/// @brief Use the given fragment libraries when building the loop.
	void use_fragments(
		utility::vector1<core::fragment::FragSetCOP> const & frag_libs);

	/// @brief Return the number of times KIC will be invoked before the
	/// LoopBuilder gives up.
	Size get_max_attempts() const;

	/// @brief Specify how many time to invoke KIC before giving up.
	void set_max_attempts(Size attempts);

	/// @brief Get the score function to be used on the next call to apply().
	core::scoring::ScoreFunctionOP get_score_function();

	/// @brief Set the score function to be used on the next call to apply().
	void set_score_function(core::scoring::ScoreFunctionOP score_function);

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static void
	get_score_function_attributes( utility::tag::AttributeList & attlist );

protected:

	/// @brief Attempt to find a reasonable loop conformation without using any
	/// information from the original coordinates.
	bool do_apply(Pose & pose, Loop const & loop) override;

private:

	/// @brief Idealize the loop with a give fold tree
	/// @detail The fold tree of the pose is restored after the idealization
	void idealize_loop(Pose & pose, Loop const & loop, core::kinematics::FoldTree &ft) const;

private:

	protocols::kinematic_closure::KicMoverOP kic_mover_;
	refiners::MinimizationRefinerOP minimizer_;
	Size max_attempts_;

};

}
}

#endif

