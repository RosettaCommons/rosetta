// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/ExtendedPoseBuilder.hh
/// @brief Builds a pose using a blueprint from a structure architect
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_components_ExtendedPoseBuilder_hh
#define INCLUDED_protocols_denovo_design_components_ExtendedPoseBuilder_hh

// Unit headers
#include <protocols/denovo_design/components/ExtendedPoseBuilder.fwd.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.fwd.hh>
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/types.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/VirtualBase.hh>

// C++ headers

#include <core/kinematics/FoldTree.fwd.hh> // AUTO IWYU For FoldTree

namespace protocols {
namespace denovo_design {
namespace components {

/// @brief Builds a pose from StructureData
class PoseBuilder : public utility::VirtualBase {
public:
	PoseBuilder();

	virtual core::pose::PoseOP
	apply( StructureData const & sd ) const = 0;

	virtual PoseBuilderOP
	clone() const = 0;

	/// @brief Returns whether or not to output debugging PDBs will be outputted which could be useful for fixing problems (true=output extra PDBs; false=don't output extra PDBs)
	bool
	debug() const { return debug_; }

	/// @brief Sets whether or not to output debugging PDBs will be outputted which could be useful for fixing problems (true=output extra PDBs; false=don't output extra PDBs)
	void
	set_debug( bool const debug ) { debug_ = debug; }

	/// @brief Writes a debugging pdb (as described above). Does nothing if "debug_" is false and only outputs a PDB if "debug_" is true
	bool
	write_debug_pdb( core::pose::Pose const & pose, std::string const & filename ) const;

private:
	/// @brief If true, some debugging PDBs will be outputted which could be useful for fixing problems
	bool debug_;
};

class ExtendedPoseBuilder : public PoseBuilder {
public:
	typedef utility::vector1< core::Size > Resids;

public:
	ExtendedPoseBuilder();
	~ExtendedPoseBuilder() override;

	core::pose::PoseOP
	apply( StructureData const & sd ) const override;

	PoseBuilderOP
	clone() const override;

private:
	/// @brief Finds the segments that have template conformations
	SegmentNames
	find_template_segments( StructureData const & sd ) const;

protected:
	/// @brief Creates a pose that contains only the template segments
	// this is protected so it can be unit tested --- should probably simplify this more object-oriented eventually
	core::pose::PoseOP
	create_template_pose( StructureData const & sd, SegmentNames const & template_segments ) const;

	/// @brief Extends a pose containing only template segments so that it also includes de novo segments
	void
	extend_pose(
		core::pose::Pose & pose,
		StructureData const & sd,
		SegmentNames const & template_segments ) const;

private:
	Resids
	calc_chain_endings( StructureData const & sd ) const;


private:
};

void
prepend_new_residues(
	core::pose::Pose & pose,
	core::Size const num_residues,
	core::Size insert_pos,
	core::Size const anchor,
	std::string const & type,
	ResidueDihedrals const & lower_dihedrals );

void
append_new_residues(
	core::pose::Pose & pose,
	core::Size const num_residues,
	core::Size insert_pos,
	core::Size const anchor,
	std::string const & type,
	ResidueDihedrals const & upper_dihedrals );

void
append_new_chain_from_template_segment(
	core::pose::Pose & pose,
	Segment const & segment );

void
append_residues_from_template_segment(
	core::pose::Pose & pose,
	Segment const & prev_segment,
	Segment const & segment );

/// @brief modifies teh ft in the pose, returns the original
core::kinematics::FoldTree
modify_ft_for_residue_insertion( core::pose::Pose & pose, core::Size const safe_res );

} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_ExtendedPoseBuilder_hh
