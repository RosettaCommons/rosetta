// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/backrub/BackrubMover.hh
/// @brief definition of BackrubMover class and functions
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_protocols_backrub_BackrubMover_hh
#define INCLUDED_protocols_backrub_BackrubMover_hh

#include <protocols/backrub/BackrubMover.fwd.hh>


// Core Headers
#include <core/id/DOF_ID_Range.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/tree/Atom.fwd.hh>

// Protocols Headers
#include <protocols/branch_angle/BranchAngleOptimizer.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/backrub/BackrubSegment.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <utility/keys/Key3Vector.hh>

// Numeric Headers
#include <numeric/IntervalSet.fwd.hh>

// C++ Headers
#include <map>

#include <utility/vector0_bool.hh>
#include <utility/vector1.hh>
#include <numeric/NumericTraits.hh>


namespace protocols {
namespace backrub {


/// @brief class for applying backrub moves to arbitrary protein segments
class BackrubMover : public protocols::canonical_sampling::ThermodynamicMover {

public:

	BackrubMover();

	BackrubMover(
		BackrubMover const & mover
	);

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	void
	init_with_options();

	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	virtual
	void
	initialize_simulation(
		core::pose::Pose & pose,
		protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
		core::Size cycle   //non-zero if trajectory is restarted
	);

	// @brief apply the backrub move to a Pose object
	void
	apply(
		Pose & pose
	);
	virtual std::string get_name() const;

	/// @brief get the segment ID for a given starting and ending atom
	/// @details
	/// if the segment ID does not exist, 0 is returned
	Size
	segment_id(
		core::id::AtomID start_atomid,
		core::id::AtomID end_atomid
	)
	{
		for (Size i = 1; i <= segments_.size(); ++i) {
			if (segments_[i].start_atomid() == start_atomid && segments_[i].end_atomid() == end_atomid) return i;
		}

		return 0;
	}

	/// @brief determine whether a given segment exists
	bool
	has_segment(
		core::id::AtomID start_atomid,
		core::id::AtomID end_atomid
	)
	{
		return segment_id(start_atomid, end_atomid);
	}

	/// @brief get a particular segment
	BackrubSegment const &
	segment(
		core::Size segment_id
	)
	{
		return segments_[segment_id];
	}

	/// @brief remove all segements
	void
	clear_segments()
	{
		segments_.clear();
	}

	/// @brief get the number of segments
	core::Size
	num_segments()
	{
		return segments_.size();
	}

	/// @brief add a segment to the mover
	Size
	add_segment(
		core::id::AtomID start_atomid,
		core::id::AtomID end_atomid,
		core::Real max_angle_disp = 0
	);

	/// @brief add segments between a list of mainchain atomids
	core::Size
	add_mainchain_segments(
		utility::vector1<core::id::AtomID> atomids,
		core::Size min_atoms,
		core::Size max_atoms
	);

	/// @brief add segments between mainchain atoms in contiguous stretches of residues
	core::Size
	add_mainchain_segments(
		utility::vector1<core::Size> resnums,
		utility::vector1<std::string> atomnames,
		core::Size min_atoms,
		core::Size max_atoms
	);

	/// @brief add segments between mainchain atoms using stored parameters
	core::Size
	add_mainchain_segments();

	/// @brief add segments between mainchain atoms using command line options
	core::Size
	add_mainchain_segments_from_options();

	/// @brief optimize branching atoms for all segment pivot atoms
	void
	optimize_branch_angles(
		Pose & pose
	);

	/// @brief get residues to pivot if no segments manually defined
	utility::vector1<core::Size> const &
	pivot_residues() const;

	/// @brief set residues to pivot if no segments manually defined
	void
	set_pivot_residues(
		utility::vector1<core::Size> const & pivot_residues
	);

	/// @brief get atom names to pivot if no segments manually defined
	utility::vector1<std::string> const &
	pivot_atoms() const;

	/// @brief set atom names to pivot if no segments manually defined
	void
	set_pivot_atoms(
		utility::vector1<std::string> const & pivot_atoms
	);

	/// @brief get minimum segment length if no segments manually defined
	core::Size
	min_atoms() const;

	/// @brief set minimum segment length if no segments manually defined
	void
	set_min_atoms(
		core::Size min_atoms
	);

	/// @brief get maximum segment length if no segments manually defined
	core::Size
	max_atoms() const;

	/// @brief set maximum segment length if no segments manually defined
	void
	set_max_atoms(
		core::Size max_atoms
	);

	/// @brief get maximum angular displacement for 4 atom segments
	core::Real
	max_angle_disp_4() const;

	/// @brief set maximum angular displacement for 4 atom segments
	void
	set_max_angle_disp_4(
		core::Real max_angle_disp_4
	);

	/// @brief get maximum angular displacement for 7 atom segments
	core::Real
	max_angle_disp_7() const;

	/// @brief set maximum angular displacement for 7 atom segments
	void
	set_max_angle_disp_7(
		core::Real max_angle_disp_7
	);

	/// @brief get maximum angular displacement slope for other atom segments
	core::Real
	max_angle_disp_slope() const;

	/// @brief set maximum angular displacement slope for other atom segments
	void
	set_max_angle_disp_slope(
		core::Real max_angle_disp_slope
	);
	
	/// @brief get whether rotation angle is customized
	virtual
	bool
	custom_angle() const;

	/// @brief set whether rotation angle is customized
	virtual
	void
	set_custom_angle(
		bool custom_angle
	);
	
	/// @brief get whether detailed balance is preserved (i.e. no branch angle optimization during moves)
	virtual
	bool
	preserve_detailed_balance() const;

	/// @brief set whether detailed balance is preserved (i.e. no branch angle optimization during moves)
	virtual
	void
	set_preserve_detailed_balance(
		bool preserve_detailed_balance
	);

	/// @brief get whether to exit during initialize_simulation if the mm_bend term isn't turned on
	bool
	require_mm_bend() const;

	/// @brief set whether to exit during initialize_simulation if the mm_bend term isn't turned on
	void
	set_require_mm_bend(
		bool require_mm_bend
	);

	/// @brief get the TorsionIDs perturbed by the mover during moves, along with their ranges
	virtual
	utility::vector1<core::id::TorsionID_Range>
	torsion_id_ranges(
		core::pose::Pose & pose
	);

	/// @brief get the DOF_IDs perturbed by the mover during moves, along with their ranges
	virtual
	utility::vector1<core::id::DOF_ID_Range>
	dof_id_ranges(
		core::pose::Pose & pose
	);

	/// @brief do a manual rotation about the given segment
	void
	rotate_segment(
		Pose & pose,
		Size segment_id,
		core::Real angle,
		utility::vector0<core::Real> & start_constants,
		utility::vector0<core::Real> & end_constants
	);

	/// @brief get a random angle for a segment subject to bond angle and rotation constraints
	core::Real
	random_angle(
		Pose & pose,
		Size segment_id,
		utility::vector0<core::Real> & start_constants,
		utility::vector0<core::Real> & end_constants
	);

	/// @brief get the branch angle optimizer stored in the mover
	protocols::branch_angle::BranchAngleOptimizer &
	branchopt()
	{
		return branchopt_;
	}

	/// @brief set the id for the next move (i.e. nonrandom)
	core::Size
	next_segment_id() const;

	/// @brief set the id for the next move (i.e. nonrandom)
	void
	set_next_segment_id(
		core::Size next_segment_id
	);

	/// @brief get id the last backrub segment moved
	core::Size
	last_segment_id() const;

	/// @brief get the name of the atom at the start of the last segment
	std::string
	last_start_atom_name() const;

	/// @brief get the name of the atom at the end of the last segment
	std::string
	last_end_atom_name() const;

	/// @brief set the rotation angle for the next move
	void
	set_next_angle(
		core::Real next_angle
	);
	
	/// @brief get the rotation angle for the next move
	core::Real
	next_angle() const;
	
	/// @brief get the last rotation angle
	core::Real
	last_angle() const;

	/// @brief update string describing the move type
	void
	update_type();

private:

	utility::vector1<protocols::backrub::BackrubSegment> segments_;
	protocols::branch_angle::BranchAngleOptimizer branchopt_;
	std::map<protocols::backrub::BackrubSegment::BondAngleKey, core::Size> bond_angle_map_;
	utility::vector1<core::Size> pivot_residues_;
	utility::vector1<std::string> pivot_atoms_;
	core::Size min_atoms_;
	core::Size max_atoms_;
	core::Real max_angle_disp_4_;
	core::Real max_angle_disp_7_;
	core::Real max_angle_disp_slope_;
	core::Size next_segment_id_;
	core::Size last_segment_id_;
	std::string last_start_atom_name_;
	std::string last_end_atom_name_;
	core::Real next_angle_;
	core::Real last_angle_;
	bool preserve_detailed_balance_;
	bool require_mm_bend_;
	bool custom_angle_;
};

/// @brief calculate constants necessary for calculating internal angles/derivatives
void
backrub_rotation_constants(
	core::kinematics::tree::AtomCOP PM2_atom,
	core::kinematics::tree::AtomCOP PM1_atom,
	core::kinematics::tree::AtomCOP P_atom,
	core::kinematics::tree::AtomCOP PP1_atom,
	core::kinematics::tree::AtomCOP PP2_atom,
	core::kinematics::tree::AtomCOP REF_atom,
	utility::vector0<double> & constants,
	core::Real const alpha_min = 0,
	core::Real const alpha_max = numeric::NumericTraits<core::Real>::pi(),
	numeric::IntervalSet<core::Real> *tau_intervals = NULL
);

/// @brief calculate internal coordinate values for any tau value
void
backrub_rotation_angles(
	utility::vector0<core::Real> const & constants,
	core::Real const tau,
	core::Real & bondange,
	core::Real & torsion1,
	core::Real & torsion2
);

} // moves
} // protocols

#endif
