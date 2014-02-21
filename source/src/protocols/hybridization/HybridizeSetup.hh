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

#ifndef INCLUDED_protocols_hybridization_HybridizeSetup_hh
#define INCLUDED_protocols_hybridization_HybridizeSetup_hh

#include <protocols/hybridization/HybridizeSetup.fwd.hh>
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

class HybridizeSetup : public utility::pointer::ReferenceCount {
public:
	HybridizeSetup();

	void init();
	
	void add_template(
					  std::string template_fn,
					  std::string cst_fn,
					  std::string symmdef_file = "NULL",
					  core::Real weight = 1.,
					  core::Real domain_assembly_weight = 0.,
					  utility::vector1<core::Size> cst_reses = utility::vector1<core::Size>(0) );
	
	void add_template(
					  core::pose::PoseOP template_pose,
					  std::string cst_fn,
					  std::string symmdef_file = "NULL",
					  core::Real weight = 1.,
					  core::Real domain_assembly_weight = 0.,
					  utility::vector1<core::Size> cst_reses = utility::vector1<core::Size>(0),
					  std::string name="default" );
	
	void add_starting_templates(core::Size t) {
		starting_templates_.push_back(t);
	}
	
	void pick_starting_template();

	void read_template_structures(utility::file::FileName template_list);
	void read_template_structures(utility::vector1 < utility::file::FileName > const & template_filenames);
	
	
	utility::vector1 <Loops>
	expand_domains_to_full_length(utility::vector1 < utility::vector1 < Loops > > all_domains, Size ref_domains_index, Size n_residues);
	
	void
	align_by_domain(utility::vector1<core::pose::PoseOP> & poses, utility::vector1 < Loops > domains, core::pose::PoseCOP & ref_pose);
	
	void
	align_by_domain(core::pose::Pose & pose, core::pose::Pose const & ref_pose, utility::vector1 <Loops> domains);
	
	void
	realign_templates(core::pose::PoseCOP ref_pose);

	void run_domain_assembly();
	
	// check fragments ... if they do not exist dynamically allocate them
	void check_and_create_fragments( core::pose::Pose const & );

	void add_big_fragments( core::fragment::FragSetOP newval ) { fragments_big_.push_back(newval); }
	void add_small_fragments( core::fragment::FragSetOP newval ) { fragments_small_.push_back(newval); }
	
	core::Size initial_template_index() const { return initial_template_index_; }
	utility::vector1 < std::string > template_file_names() {return template_fn_;}
	utility::vector1 < core::pose::PoseOP > template_poses() const { return templates_; }
	utility::vector1 < std::string > template_cst_fn() const {return template_cst_fn_;}
	utility::vector1 < core::Real > template_wts() const {return template_weights_;}
	utility::vector1 < core::Real > domain_assembly_weights() const {return domain_assembly_weights_;}
	utility::vector1 < std::string > symmdef_files() const {return symmdef_files_;}

	utility::vector1 < protocols::loops::Loops > template_chunks() const {return template_chunks_;}
	utility::vector1 < protocols::loops::Loops > template_contigs() const {return template_contigs_;}
	
	utility::vector1 <core::fragment::FragSetOP> fragments_big() const { return fragments_big_; } // 9mers/fragA equivalent in AbrelaxApplication
	utility::vector1 <core::fragment::FragSetOP> fragments_small() const { return fragments_small_; } // 3mers/fragB equivalent in AbrelaxApplication
	
	bool domain_assembly() const {return domain_assembly_;}
	void set_domain_assembly(bool v) {domain_assembly_ = v;}
	void set_align_domains_to_template(bool v) {align_domains_to_template_ = v;}
	void set_align_domains_to_pose(bool v) {align_domains_to_pose_ = v;}
	void add_fragments_small(core::fragment::FragSetOP v) {
		fragments_small_.push_back(v);
	}
	void add_fragments_big(core::fragment::FragSetOP v) {
		fragments_big_.push_back(v);
	}
	
	//DDomain
	void set_domain_hcut(core::Real v) {
		hcut_ = v;
	}
	void set_domain_pcut(core::Real v) {
		pcut_ = v;
	}
	void set_domain_length(core::Size v) {
		length_ = v;
	}
	
	core::Size nres_tgt_asu() const {return nres_tgt_;}
	core::Size nres_protein_asu() const {return nres_protein_tgt_;}
	void set_nres_tgt(core::Size nres_tgt_in) {
		nres_tgt_ = nres_tgt_in;
	}
	void set_nres_protein( core::Size nres_protein_tgt_in ) {
		nres_protein_tgt_ = nres_protein_tgt_in;
	}

	// hetatm
	void set_add_hetatm(bool v) {
		add_hetatm_ = v;
	}
	void set_hetatm_self_cst_weight(core::Real v) {
		hetatm_self_cst_weight_ = v;
	}
	void set_hetatm_prot_cst_weight(core::Real v) {
		hetatm_prot_cst_weight_ = v;
	}
	bool add_hetatm() const {return add_hetatm_;}
	core::Real hetatm_self_cst_weight() const {return hetatm_self_cst_weight_;}
	core::Real hetatm_prot_cst_weight() const {return hetatm_prot_cst_weight_;}
	
	virtual HybridizeSetupOP clone() const;

private:
	// parsible options
	utility::vector1 < core::Size > starting_templates_;
	core::Size initial_template_index_;
	
	// 1mer fragment insertion weight where fragments are not allowed (across anchors) , vs. chunk insertion + big and small frags
	core::Real frag_1mer_insertion_weight_;
	// small fragment insertion weight where big fragments are not allowed (across anchors) , vs. chunk insertion + big frags
	core::Real small_frag_insertion_weight_;
	core::Real big_frag_insertion_weight_; // fragment insertion weight, vs. chunk insertion + small gap frags
	core::Real frag_weight_aligned_; // fragment insertion to the aligned region, vs. unaligned region
	bool auto_frag_insertion_weight_; // automatically set fragment insertion weights
	core::Size max_registry_shift_;
	bool domain_assembly_, add_hetatm_;
	bool align_domains_to_template_, align_domains_to_pose_, add_non_init_chunks_;
	core::Real hetatm_self_cst_weight_, hetatm_prot_cst_weight_;
	std::string disulf_file_;
	core::Size cartfrag_overlap_;
	
	// ddomain options
	core::Real pcut_,hcut_;
	core::Size length_;
	
	// abinitio frag9,frag3 flags
	utility::vector1 <core::fragment::FragSetOP> fragments_big_;  // 9mers/fragA equivalent in AbrelaxApplication
	utility::vector1 <core::fragment::FragSetOP> fragments_small_; // 3mers/fragB equivalent in AbrelaxApplication
	
	// native pose, aln
	core::pose::PoseOP native_;
	core::sequence::SequenceAlignmentOP aln_;
	
	// template information (all from hybrid.config)
	utility::vector1 < core::pose::PoseOP > templates_;
	utility::vector1 < std::string > template_fn_;
	utility::vector1 < std::string > template_cst_fn_;
	utility::vector1 < std::string > symmdef_files_;
	utility::vector1 < core::Real > template_weights_;
	utility::vector1 < core::Real > domain_assembly_weights_;
	utility::vector1 < protocols::loops::Loops > template_chunks_;
	utility::vector1 < protocols::loops::Loops > template_contigs_;
	utility::vector1 < utility::vector1<core::Size> > template_cst_reses_;
	core::Real template_weights_sum_;
	
	utility::vector1< protocols::loops::Loops > domains_;
	
	// strand pairings
	std::string pairings_file_;
	utility::vector1<core::Size> sheets_;
	utility::vector1<core::Size> random_sheets_;
	bool filter_templates_;
	utility::vector1< std::pair< core::Size, core::Size > > strand_pairs_;
	
	// number of residues in asu without VRTs
	core::Size nres_tgt_;
	core::Size nres_protein_tgt_;
};
	
class HybridizeSetupMover : public protocols::moves::Mover {
public:
	HybridizeSetupMover();

	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag( utility::tag::TagCOP const, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, Pose const & );
	
private:
	HybridizeSetupOP hybridize_setup_;
};

} // hybridization
} // protocols

#endif
