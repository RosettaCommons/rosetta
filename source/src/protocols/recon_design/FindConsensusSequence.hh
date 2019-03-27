// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/recon_design/FindConsensusSequence.hh
/// @brief Takes in multiple poses from the VectorPoseJobDistributor and finds the consensus
/// sequence that optimizes energy of all input poses. Used in conjuction with MSDMover
/// at the end of a protocol to make sure that you end up with one multistate solution.
/// Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_protocols_recon_design_FindConsensusSequence_hh
#define INCLUDED_protocols_recon_design_FindConsensusSequence_hh

#include <protocols/recon_design/FindConsensusSequence.fwd.hh>
#include <protocols/moves/VectorPoseMover.hh>
#include <protocols/moves/Mover.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.fwd.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace recon_design {

/// @brief Takes in multiple poses from the VectorPoseJobDistributor and finds the consensus
/// sequence that optimizes energy of all input poses. Used in conjuction with MSDMover
/// at the end of a protocol to make sure that you end up with one multistate solution.
/// Only accessible through recon application.
class FindConsensusSequence : public moves::VectorPoseMover {

public:

	/// @brief empty constructor
	FindConsensusSequence();

	~FindConsensusSequence();

	moves::MoverOP clone() const override;
	moves::MoverOP fresh_instance() const override;
	std::string get_name() const override;

	/// @brief Read options from RosettaScripts
	void parse_my_tag(
		TagCOP,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		Pose const & ) override;

	/// @brief Specify XML schema
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Movers derived from VectorPoseMover must define two apply methods:
	/// apply and apply_mpi, depending on whether the protocol is run in MPI or not
	/// Running in MPI will distribute each pose to a separate node, so the mover will
	/// operate on one pose
	void apply( Pose & pose ) override;
	void apply_mpi( Pose & pose ) override;

	/// @brief Populates designable_residues_ with a list of designable residues corresponding
	/// to the poses in poses_, corresponding element-wise (designables_residues[0] matches
	/// to poses_[0], etc)
	void parse_resfiles();

	/// @brief Get the resfile corresponding to the pose at index. If
	/// only one resfile is present then it will be returned regardless
	/// of the value of index
	std::string resfile_at ( core::Size index );

	/// @brief Based on all the input poses, find the optimal AA at position
	/// res_link_index. candidate_AAs specifies all of the AAs present in any
	/// of poses_ at position res_link_index. Fitness is defined implicitly
	/// as the sum of energy of all poses.
	void pick_consensus_AA( core::Size res_link_index,
		utility::vector1<std::string> candidate_AAs );

	void pick_consensus_AA_mpi( core::pose::Pose & pose,
		utility::vector1<std::string> candidate_AAs,
		core::Size pose_position );

	/// @brief Set up the packer to be used when substituting different
	/// candidate AAs and repacking before evaluating energy
	void initialize_packer();

	/// @brief Place a candidate AA (trial_AA) onto a pose (pose_copy)
	/// at position pose_position and repack and measure the energy.
	/// Function returns the energy of that pose with the trial_AA
	core::Real test_AA( core::pose::PoseOP pose_copy,
		std::string trial_AA,
		core::Size pose_position );

	/// Getters and setters
	core::scoring::ScoreFunctionOP score_function() const;
	void score_function( core::scoring::ScoreFunctionOP );

	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP );

	void resfiles ( utility::vector1< std::string > & resfiles );

	utility::vector1< utility::vector1< core::Size > > designable_residues();

private:
	/// Scorefunction to be used
	core::scoring::ScoreFunctionOP sfxn_;

	/// TaskFactory to specify tasks in repacking when substituting AAs
	core::pack::task::TaskFactoryOP task_factory_;

	/// List of resfile file names, must be either 1 resfile, or same length as poses_
	utility::vector1< std::string > resfiles_;

	/// Vector of designable residues for each pose in poses_.
	/// Must be same length as poses_.
	/// Also each pose must have same number of designable residues
	/// (i.e. designable_residues[1].size() == designable_residues[2].size())
	utility::vector1< utility::vector1< core::Size > > designable_residues_;

	/// Packer to be used when substituting AAs
	minimization_packing::PackRotamersMoverOP packer_;

	/// Output extra debug messages
	bool debug_ = false;

	/// Variables used in MPI
	core::Size rank_ = 0;
	core::Size n_procs_ = 1;
	bool master_ = false;

};

} //recon_design
} //protocols

#endif
