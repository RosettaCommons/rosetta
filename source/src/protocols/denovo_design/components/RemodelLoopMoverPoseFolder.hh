// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/RemodelLoopMoverPoseFolder.hh
/// @brief Folds residues in a pose using RemodelLoopMover
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_components_RemodelLoopMoverPoseFolder_hh
#define INCLUDED_protocols_denovo_design_components_RemodelLoopMoverPoseFolder_hh

// Unit headers
#include <protocols/denovo_design/components/RemodelLoopMoverPoseFolder.fwd.hh>
#include <protocols/denovo_design/components/PoseFolder.hh>

// Protocol headers
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// Core headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace denovo_design {
namespace components {

///@brief Folds residues in a pose using RemodelLoopMover
class RemodelLoopMoverPoseFolder : public protocols::denovo_design::components::PoseFolder {
public:
	typedef protocols::denovo_design::components::PoseFolder PoseFolder;
	typedef protocols::denovo_design::components::PoseFolderOP PoseFolderOP;

public:
	RemodelLoopMoverPoseFolder();

	virtual ~RemodelLoopMoverPoseFolder();

	static std::string
	class_name();

	PoseFolderOP
	clone() const;

	/// @brief performs folding
	/// @param pose    - The pose to be folded, with all residues added.  The pose should be prepared with
	///                  any necessary cutpoints added before giving to the PoseFolder. Torsions in the pose
	///                  should be adjusted, and no residues should be added or removed.
	/// @param movable - Subset of residues for which new backbone conformations will be sampled. Residues
	///                  specified as 'True' in movable must also be present in one or more Loops in order
	///                  to be folded. Movable's size must match pose.size()
	/// @param loops   - Loops to be folded.  Cutpoints specified here must be match the cutpoints found in
	///                  the pose. Residues not within any loop should not be folded. Residues contained
	///                  in a loop but not in the movable set should not be folded.
	/// @throws EXCN_Fold if anything goes wrong in folding. Derived classes should throw this.
	virtual void
	apply(
		core::pose::Pose & pose,
		core::select::residue_selector::ResidueSubset const & movable,
		protocols::loops::Loops const & loops ) const;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

public:
	void
	set_scorefxn( core::scoring::ScoreFunction const & sfxn );

public:
	static core::scoring::ScoreFunctionOP
	default_score_function();

private:
	void
	remove_cutpoints( StructureData & sd, protocols::loops::Loops const & loops ) const;

	protocols::moves::MoverOP
	create_remodel_loop_mover(
		core::pose::Pose const & pose,
		StructureData const & sd,
		core::select::residue_selector::ResidueSubset const & movable,
		protocols::loops::Loops const & loops ) const;

	core::kinematics::MoveMapOP
	create_false_movemap(
		core::pose::Pose const & pose,
		core::select::residue_selector::ResidueSubset const & movable ) const;

	core::scoring::ScoreFunctionOP
	create_scorefxn( core::pose::Pose const & pose ) const;

private:
	core::scoring::ScoreFunctionCOP scorefxn_;
};

} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_RemodelLoopMoverPoseFolder_hh
