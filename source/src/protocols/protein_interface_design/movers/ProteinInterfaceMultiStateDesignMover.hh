// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ProteinInterfaceMultiStateDesignMover.hh
/// @brief
/// @author Sarel Fleishman (sarelf@uw.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_ProteinInterfaceMultiStateDesignMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_ProteinInterfaceMultiStateDesignMover_hh
#include <protocols/protein_interface_design/movers/ProteinInterfaceMultiStateDesignMover.fwd.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <protocols/multistate_design/MultiStatePacker.fwd.hh>
// AUTO-REMOVED #include <protocols/genetic_algorithm/GeneticAlgorithm.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <protocols/genetic_algorithm/GeneticAlgorithm.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

///@brief wraps protein-interface specific considerations around the general multistate design / genetic algorithm framework
class ProteinInterfaceMultiStateDesignMover : public protocols::simple_moves::PackRotamersMover {
public:
	typedef multistate_design::PosType PosType;
	typedef core::pack::task::PackerTaskCOP PackerTaskCOP;
	typedef utility::pointer::owning_ptr< genetic_algorithm::GeneticAlgorithm >
		GeneticAlgorithmOP;

public:
	ProteinInterfaceMultiStateDesignMover();
	virtual ~ProteinInterfaceMultiStateDesignMover();
	virtual void apply( Pose & );
	virtual std::string get_name() const;
	void output_results( Pose & );

	virtual void parse_my_tag(
		TagCOP const,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		moves::Movers_map const &,
		Pose const & );
	virtual moves::MoverOP fresh_instance() const;
	virtual moves::MoverOP clone() const;
	void restrict_sequence_profile( core::pose::Pose const & pose, core::pack::task::PackerTaskOP const ptask ) const;
	unsigned long sequence_space( core::pack::task::PackerTaskCOP ptask ) const;
	void output_alternative_states( core::pose::Pose const & output_pose ) const;

private:
	void initialize( Pose & );
	
	using protocols::simple_moves::PackRotamersMover::run;
	void run();
/// @brief add target and competitor states
	void add_states( Pose const & );

private:
	GeneticAlgorithmOP gen_alg_;
	// direct use of MultiStatePacker is only for outputting results
	// (GeneticAlgorithm also holds pointer to it and uses it heavily)
	multistate_design::MultiStatePackerOP multistate_packer_;
	// option flags/parameters: constructor defaults to command line options
	// parse_my_tag method may change them
	core::Size rb_jump_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Size generations_, pop_size_, num_packs_, pop_from_ss_, numresults_;
	core::Real fraction_by_recombination_, mutate_rate_, boltz_temp_, anchor_offset_;
	// checkpointing options
	std::string checkpoint_prefix_;
	core::Size checkpoint_interval_;
	bool checkpoint_gz_, checkpoint_rename_;

	bool unbound_, unfolded_, input_is_positive_;
	bool use_unbound_for_sequence_profile_;// use a poly-ala unbound state to decide which residues to allow at each position
	core::Real bump_threshold_; //for residues to be allowed in seq_profile
	bool compare_energy_to_ground_state_; // set internally by the mover
///states
	utility::vector1< core::pose::PoseOP > state_poses_;
	utility::vector1< core::pose::PoseOP > saved_state_poses_; /// for dumping out pdbs of negative states at the end
	utility::vector1< bool > state_positive_;
	utility::vector1< bool > state_unfolded_;
	utility::vector1< bool > state_unbound_;
	utility::vector1< core::pack::task::TaskFactoryCOP > state_task_factory_; // a state-specific task factory that will be applied to each state in addition to the task factory that is implied by the design process
	std::string fname_prefix_;
};

} //namespace movers
} // namespace protein_interface_design
} // namespace protocols

#endif
