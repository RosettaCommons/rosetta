// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#ifndef INCLUDED_protocols_hybridization_HybridizeProtocol_hh
#define INCLUDED_protocols_hybridization_HybridizeProtocol_hh

#include <protocols/hybridization/HybridizeProtocol.fwd.hh>
#include <protocols/hybridization/FoldTreeHybridize.hh>

#include <protocols/moves/Mover.hh>

#include <protocols/loops/Loops.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <utility/file/FileName.hh>

namespace protocols {
namespace hybridization {

class HybridizeProtocol : public protocols::moves::Mover {

public:
	HybridizeProtocol();

	void init();

	void update_last_template();

	void add_template(
		std::string template_fn,
		std::string cst_fn,
		std::string symmdef_file = "NULL",
		core::Real weight = 1.,
		utility::vector1<char> rand_chains = utility::vector1<core::Size>(0) );

	void add_null_template(
		core::pose::PoseOP template_pose,
		std::string cst_fn,
		std::string symmdef_file = "NULL",
		core::Real weight = 1. );

	void add_template(
		core::pose::PoseOP template_pose,
		std::string cst_fn,
		std::string symmdef_file = "NULL",
		core::Real weight = 1.,
		utility::vector1<char> rand_chains = utility::vector1<core::Size>(0),
		std::string filename="default" );

	void validate_template(
		std::string filename,
		std::string fasta,
		core::pose::PoseOP template_pose,
		bool & align_pdb_info
	);

	// pick a starting template using assigned weights
	core::Size pick_starting_template();

	// optional rb docking after stage 1
	void do_intrastage_docking(core::pose::Pose & pose);

	// align all templates to a reference model
	void
	align_templates_by_domain(core::pose::PoseOP & ref_pose, utility::vector1 <Loops> domains=utility::vector1 <Loops>());

	// align 1 templates to a reference model
	void
	align_by_domain(core::pose::Pose & pose, core::pose::Pose const & ref_pose, utility::vector1 <Loops> domains);

	void domain_parse_templates(core::Size nres);

	//fpd  alternate version of stage 1 with no hybridization
	void
	initialize_and_sample_loops(
		core::pose::Pose &pose,
		core::pose::PoseOP chosen_templ,
		protocols::loops::Loops template_contigs_icluster,
		core::scoring::ScoreFunctionOP scorefxn);

	//fpd   add fragment-derived constraints to the pose
	void
	add_fragment_csts( core::pose::Pose &pose );

	// check fragments ... if they do not exist dynamically allocate them
	void check_and_create_fragments( Pose & );

	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag( TagCOP, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, Pose const & );

	// set options
	void set_batch_relax( core::Size newval ) { batch_relax_ = newval; }
	void add_big_fragments( core::fragment::FragSetOP newval ) { fragments_big_.push_back(newval); }
	void add_small_fragments( core::fragment::FragSetOP newval ) { fragments_small_.push_back(newval); }
	void set_stage1_scorefxn( core::scoring::ScoreFunctionOP newval ) { stage1_scorefxn_ = newval; }
	void set_stage2_scorefxn( core::scoring::ScoreFunctionOP newval ) { stage2_scorefxn_ = newval; }
	void set_stage2pack_scorefxn( core::scoring::ScoreFunctionOP newval ) { stage2pack_scorefxn_ = newval; }
	void set_stage2min_scorefxn( core::scoring::ScoreFunctionOP newval ) { stage2min_scorefxn_ = newval; }
	void set_stage1_increase_cycles( core::Real newval ) { stage1_increase_cycles_ = newval; }
	void set_stage2_increase_cycles( core::Real newval ) { stage2_increase_cycles_ = newval; }
	void set_fullatom_scorefxn( core::scoring::ScoreFunctionOP newval ) { fa_scorefxn_ = newval; }

private:
	// parsible options
	core::Real stage1_probability_, stage1_increase_cycles_, stage2_increase_cycles_, stage25_increase_cycles_;
	core::Size stage1_1_cycles_, stage1_2_cycles_, stage1_3_cycles_, stage1_4_cycles_;
	core::Real stage2_temperature_;

	core::Real frag_1mer_insertion_weight_;   // 1mer insertion weights
	core::Real small_frag_insertion_weight_;  // 3mer insertion weights
	core::Real big_frag_insertion_weight_;    // 9mer insertion weights
	core::Real chunk_insertion_weight_;       // chunk insertion weights
	core::Real frag_weight_aligned_;          // fragment insertion rate in templates
	bool auto_frag_insertion_weight_;         // automatically set fragment insertion weights
	core::Size max_registry_shift_;

	bool add_hetatm_, realign_domains_, realign_domains_stage2_, no_global_frame_, linmin_only_;
	core::Real add_non_init_chunks_;
	bool seqfrags_only_, nofragbias_, skip_long_min_;
	core::Real hetatm_self_cst_weight_, hetatm_prot_cst_weight_;
	core::scoring::ScoreFunctionOP stage1_scorefxn_, stage2_scorefxn_, fa_scorefxn_;
	core::scoring::ScoreFunctionOP stage2pack_scorefxn_, stage2min_scorefxn_;  // cenrot
	std::string fa_cst_fn_;
	std::string disulf_file_;
	core::Size cartfrag_overlap_;

	// more options
	bool csts_from_frags_;  // generate dihedral constraints from fragments
	int max_contig_insertion_;  // don't insert contigs larger than this size
	bool min_after_stage1_; // tors min after stage1
	core::Real fragprob_stage2_;
	core::Real randfragprob_stage2_;

	// use cenrot mode
	bool cenrot_;

	// ddomain parameters + domain definitions
	core::Real pcut_,hcut_;
	core::Size length_;
	utility::vector1< utility::vector1< protocols::loops::Loops > > domains_all_templ_;

	// relax parameters
	core::Size batch_relax_, relax_repeats_;

	// abinitio fragment sets
	utility::vector1 <core::fragment::FragSetOP> fragments_big_;   // 9mers/fragA equivalent in AbrelaxApplication
	utility::vector1 <core::fragment::FragSetOP> fragments_small_; // 3mers/fragB equivalent in AbrelaxApplication

	// native pose, aln
	core::pose::PoseOP native_;
	core::sequence::SequenceAlignmentOP aln_;

	// template information
	utility::vector1 < core::Size > non_null_template_indices_;
	utility::vector1 < core::pose::PoseOP > templates_;  // template poses
	utility::vector1 < core::pose::PoseOP > templates_aln_;  // aligned template PDBs (deep copy for multidom, shallow for single dom)
	utility::vector1 < std::string > template_fn_;       // template file tags
	utility::vector1 < std::string > template_cst_fn_;   // template constraint file tags
	utility::vector1 < std::string > symmdef_files_;     // template symmdef files
	utility::vector1 < core::Real > template_weights_;   // template weights
	utility::vector1 < protocols::loops::Loops > template_chunks_;    // template secstruct definitions
	utility::vector1 < protocols::loops::Loops > template_contigs_;   // template continuous pieces
	utility::vector1 < utility::vector1<char> > randomize_chains_;  // per-template chain randomization

	// strand pairings
	std::string pairings_file_;
	utility::vector1<core::Size> sheets_;
	utility::vector1<core::Size> random_sheets_;
	bool filter_templates_;
	utility::vector1< std::pair< core::Size, core::Size > > strand_pairs_;

	// per-residue control
	utility::vector1<bool> residue_sample_template_; // using template fragments
	utility::vector1<bool> residue_sample_abinitio_; // using torsion-based ab initio fragments
	utility::vector1<core::Size> residue_max_registry_shift_; // restraints between chains

	// local docking ater stage 1
	//    this should be incorporated in stage 1 rather than after for fold and dock style sampling!
	bool jump_move_;
	core::Size jump_move_repeat_;

	// constraint
	bool keep_pose_constraint_;
	utility::vector1 < core::Size > user_csts_;
};

} // hybridization
} // protocols

#endif
