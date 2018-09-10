// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/AntibodySeqDesignTFCreator.hh
/// @brief Class for creating a TaskFactory from SeqDesignOptions for use in various protocols
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_AntibodySeqDesignTFCreator_hh
#define INCLUDED_protocols_antibody_design_AntibodySeqDesignTFCreator_hh

#include <protocols/antibody/design/AntibodySeqDesignTFCreator.fwd.hh>


#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.fwd.hh>
#include <protocols/antibody/task_operations/AddCDRProfilesOperation.fwd.hh>
#include <protocols/task_operations/ConservativeDesignOperation.fwd.hh>
#include <protocols/task_operations/ResidueProbDesignOperation.fwd.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/simple_task_operations/RestrictToLoopsAndNeighbors.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <map>

namespace protocols {
namespace antibody {
namespace design {



/// @brief Create a TaskFactory or individual TaskOperations based on a set of options
/// These options are mainly for per-CDR and per-framework control of sequence design
///
class AntibodySeqDesignTFCreator : public utility::pointer::ReferenceCount {
public:

	/// @brief Constructor with default set of options.  You probably do not want this!
	AntibodySeqDesignTFCreator(AntibodyInfoCOP ab_info, bool force_north_paper_db =  false);

	/// @brief Constructor with a CDRSeqDesignOptionsOP for each CDR
	AntibodySeqDesignTFCreator(
		AntibodyInfoCOP ab_info,
		utility::vector1<CDRSeqDesignOptionsOP> const design_options,
		bool force_north_paper_db = false,
		core::Size stem_size = 2);


	//AntibodySeqDesignTFCreator(AntibodySeqDesignTFCreator const & src);

	virtual ~AntibodySeqDesignTFCreator();

	AntibodySeqDesignTFCreator( AntibodySeqDesignTFCreator const & src );

	AntibodySeqDesignTFCreatorOP
	clone () const;

	/// @brief Set design options for single CDR
	void
	set_cdr_design_options(CDRNameEnum cdr, CDRSeqDesignOptionsCOP design_options);

	/// @brief Set design options for All CDRs
	void
	set_cdr_design_options(utility::vector1<CDRSeqDesignOptionsOP> const design_options);

	/// @brief Get modifiable options
	utility::vector1<CDRSeqDesignOptionsOP>
	get_cdr_design_options();

	/// @brief Get modifiable options
	CDRSeqDesignOptionsOP
	get_cdr_design_options(CDRNameEnum cdr);

	// Undefined, commenting out to fix PyRosetta build  void set_antibody_info(AntibodyInfoCOP ab_info);


public:
	//////// TaskFactory and TaskOperation Generation ///////////////////

	/// @brief Normally, we do not want to design the antigen residues - but perhaps someday we will.
	/// If False, turn design off for the antigen for any TF creation.
	void
	design_antigen( bool antigen_design );

	/// @brief Design any framework residues included in the task. Default True.
	/// If False, turn design off for the framework for any TF creation.
	void
	design_framework( bool framework_design );

	/// @brief Create the FULL TaskFactory for antibody sequence design.
	/// There are no limits to repacking or design - AKA NO RestrictToLoops.
	/// Optionally disable framework residues for design
	/// Optionally disable antigen residues for design.
	/// Optionally disable framework conservative residues
	/// -> Combine with other Tasks for general design requirements.
	core::pack::task::TaskFactoryOP
	generate_tf_seq_design(core::pose::Pose const & pose, bool disable_non_designing_cdrs = false);

	/// @brief Create the TaskFactory used for sequence design during the GraftDesign stage.
	/// This limits design to the passed CDR, optionally designing neighbor CDRs specified in options and included in min.
	///
	/// @details
	///  Basically will create a TF including CDRs you are minimizing and any neighbors.
	///  It will then use settings in your Options classes and class settings such as framework or antigen design to create the TF,
	///  disabling CDRs that are not allowed to design as well as framework or antigen regions.
	///
	/// Optionally design any neighbor framework residues
	/// Optionally design any neighbor antigen residues
	/// Optionally disable framework conservative residues
	///
	core::pack::task::TaskFactoryOP
	generate_tf_seq_design_graft_design(
		core::pose::Pose const & pose,
		CDRNameEnum cdr,
		utility::vector1<bool>const & neighbor_cdr_min);

	/// @brief Explicitly Generate a TF for framework-optimization only based design.
	/// This will only design neighboring framework residues of the grafted CDR.
	/// It will ignore the CDRs set to design in the SeqDesign options.
	//core::pack::task::TaskFactoryOP
	//generate_tf_seq_design_graft_design_framework_only(
	// core::pose::Pose const & pose,
	// CDRNameEnum cdr,
	// utility::vector1<bool> neighbor_cdr_min);

	/// @brief Setup the CDRProfilesOperation and conservative/basic design
	///  according to primary and fallback strategies and which CDRs are set to design.
	///  Pre-load the data
	///
	task_operations::AddCDRProfilesOperationOP
	generate_task_op_cdr_profiles(core::pose::Pose const & pose);

	/// @brief Create a TaskOp to limit Packing and Design to only CDR loops with design on.  Use neighbor distance.
	protocols::simple_task_operations::RestrictToLoopsAndNeighborsOP
	generate_task_op_cdr_design(core::pose::Pose const & pose, bool design_neighbors = true) const;

	/// @brief Create a TaskOp for only CDR loops set to True in the boolean vector.
	protocols::simple_task_operations::RestrictToLoopsAndNeighborsOP
	generate_task_op_cdr_design(core::pose::Pose const & pose, utility::vector1<bool> cdrs, bool design_neighbors = true) const;

	/// @brief Create a TaskOp to limit Packing and Design to CDR loops and neighbors.
	protocols::simple_task_operations::RestrictToLoopsAndNeighborsOP
	generate_task_op_all_cdr_design( core::pose::Pose const & pose, bool design_neighbors = true ) const;


	///Some helper functions

	/// @brief Turns off CDRs for design that are set to off.
	void
	disable_design_for_non_designing_cdrs(core::pack::task::TaskFactoryOP tf, const core::pose::Pose & pose);

	void
	disable_proline_design(core::pack::task::TaskFactoryOP tf, const core::pose::Pose & pose) const;

	void
	disable_disallowed_aa( core::pack::task::TaskFactoryOP tf, const core::pose::Pose & pose ) const;

	/// @brief Create a TaskOp for profile-based design of CDRs according to SeqDesign options.
	//task_operations::ResidueProbDesignOperationOP
	//generate_task_op_cdr_profile(core::pose::Pose const & pose);

	/// @brief Create a TaskOp for conservative-based design of CDRs according to SeqDesign options.
	//task_operations::ConservativeDesignOperationOP
	//generate_task_op_cdr_conservative(core::pose::Pose const & pose);



public:

	///////// General Options to create TF //////////////

	/// @brief Repack neighbors of CDR's being designed within this distance.
	void
	neighbor_detection_dis( core::Real const neighbor_distance );

	/// @brief Keep proline fixed for design. If using profile design, this should not really come into play.
	void
	design_proline( bool const setting );


	/// @brief Use the Conservative Design TaskOP if designing Framework residues.
	/// Default true. Recommended.
	void
	set_design_framework_conservative( bool design_framework_conservative );

	/// @brief Enable design of 100% conserved framework positions during TF generation.
	/// Default false.  Will be expanded.
	void
	set_design_framework_conserved_res( bool design_framework_conserved_res );


	/// @brief Use these weights during probabilistic design for data that is normally zero.
	void
	set_zero_prob_weight_at(core::Real const weight);

	/// @brief Use conservative mutations (or alternative method) instead of using cluster sequence probabilities for design
	/// if the number of sequences in the particular CDR's cluster probability data is lower than this cutoff. Default is 10.  This is why we try and stay in type 1 lengths during graft.
	void
	set_probability_data_cutoff(core::Size const cutoff);

	/// @brief Enable design of the first 2 and last 3 residues of the H3 loop.  These are off by default as to help hinder the transition from extended
	/// to kinked and vice versa during sequence design.
	void
	set_design_H3_stem(bool design_H3_stem);

	///@brief Set to sample all available AAs per position instead of sampling based on weights
	void
	set_no_probability( bool no_probability);

private:

	void
	setup_default_options();

	void
	read_command_line_options();


private:
	/// @brief Get a LoopsOP for CDRs set to design
	protocols::loops::LoopsOP
	get_design_cdr_loops( core::pose::Pose const & pose, core::Size stem_size = 0) const;

	/// @brief Explicitly get design cdr loops with stem.
	protocols::loops::LoopsOP
	get_design_cdr_loops_with_stem( core::pose::Pose const & pose ) const;

private:

	protocols::task_operations::ConservativeDesignOperationOP
	get_framework_conservative_op(core::pose::Pose const & pose);

	protocols::simple_task_operations::RestrictToLoopsAndNeighborsOP
	get_general_loop_task_op(protocols::loops::LoopsOP loops, bool design_neighbors = true ) const;

	/// @brief Add restrictions for non-CDR positions and residue types according to options.
	void
	add_extra_restrict_operations(core::pack::task::TaskFactoryOP tf, const core::pose::Pose & pose) const;


private:

	AntibodyInfoCOP ab_info_;

	utility::vector1<CDRSeqDesignOptionsOP> cdr_design_options_;

	core::Real zero_prob_weight_;
	core::Real neighbor_dis_;

	core::Size prob_cutoff_;
	core::Size profile_picking_rounds_;
	core::Size stem_size_;

	bool design_proline_;

	bool design_antigen_;
	bool design_framework_;

	bool design_framework_conservative_;
	bool design_framework_conserved_res_;

	bool design_h3_stem_;

	bool force_north_paper_db_;

	utility::vector1<bool> no_data_cdrs_; //CDRs that were unable to load profile-based data.
	bool use_outliers_;
	bool no_probability_ = false;

};


}
}
}


#endif //INCLUDED_ AntibodySeqDesignTFCreator.hh
