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
/// @author Yifan Song
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_hybridization_CartesianSampler_hh
#define INCLUDED_protocols_hybridization_CartesianSampler_hh

#include <protocols/hybridization/InsertChunkMover.hh>
#include <protocols/hybridization/CartesianSampler.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>

#include <core/scoring/ScoreFunction.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>
#include <numeric/random/WeightedSampler.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

class CartesianSampler: public protocols::moves::Mover {
public:
	CartesianSampler();
	CartesianSampler( utility::vector1<core::fragment::FragSetOP> fragments_in );

	// initialize options to defaults
	void init();

	// run the protocol
	void apply(core::pose::Pose & pose) override;

	// set the centroid scorefunction
	void set_scorefunction(core::scoring::ScoreFunctionOP scorefxn_in) { scorefxn_=scorefxn_in; }

	void set_mc_scorefunction(core::scoring::ScoreFunctionOP scorefxn_in) { mc_scorefxn_=scorefxn_in; }

	// set the fullatom scorefunction (only used for some option sets)
	void set_fa_scorefunction(core::scoring::ScoreFunctionOP scorefxn_in) { fa_scorefxn_=scorefxn_in; }

	std::string get_name() const override { return "CartesianSampler"; }

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const & ,
		moves::Movers_map const & ,
		core::pose::Pose const & pose ) override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	// accessors
	void set_userpos(std::set<core::Size> const &user_pos_in) { user_pos_ = user_pos_in; }
	void set_strategy( std::string strategy_in ) { fragment_bias_strategies_.push_back( strategy_in ); }
	void set_rms_cutoff(core::Real rms_cutoff_in) { rms_cutoff_ = rms_cutoff_in; }
	void set_fullatom(bool fullatom_in) {fullatom_ = fullatom_in; }
	void set_bbmove(bool bbmove_in) { bbmove_ = bbmove_in; }
	void set_restore_csts(bool restore_csts_in) { restore_csts_ = restore_csts_in; }
	void set_frag_sizes(utility::vector1<core::Size> const &frag_sizes_in) { frag_sizes_ = frag_sizes_in; }
	void set_nminsteps(core::Size nminsteps_in) { nminsteps_ = nminsteps_in; }

	void set_ncycles(core::Size ncycles_in) { ncycles_=ncycles_in; }
	void set_overlap(core::Size overlap_in) { overlap_=overlap_in; }
	void set_nfrags(core::Size nfrags_in) { nfrags_=nfrags_in; }
	void set_temp(core::Real temp_in) { temp_=temp_in; }


protected:
	// apply a sequence fragment
	bool apply_frame( core::pose::Pose & pose, core::fragment::Frame &frame );


	void compute_fragment_bias( core::pose::Pose & pose );


	void update_fragment_library_pointers( );


	void apply_constraints( core::pose::Pose & pose );

	// get frag->pose transform, return RMS
	core::Real get_transform(
		core::pose::Pose const &pose, core::pose::Pose const &frag, core::Size startpos,
		core::Vector &preT, core::Vector &postT, numeric::xyzMatrix< core::Real > &R);

	// transform fragment
	void apply_transform( core::pose::Pose &frag, core::Vector const &preT, core::Vector const &postT, numeric::xyzMatrix< core::Real > const &R);

	// apply endpoint constraints to fragment
	void apply_fragcsts( core::pose::Pose &working_frag, core::pose::Pose const &pose, core::Size start );

private:
	// parameters
	core::Size ncycles_, overlap_, nminsteps_;
	core::Real rms_cutoff_;

	// fragments
	utility::vector1<core::Size> frag_sizes_;
	core::Size nfrags_;
	utility::vector1<core::fragment::FragSetOP> fragments_;
	utility::vector1<boost::unordered_map<core::Size, core::fragment::Frame> > library_; // [ fragsize: { pos: frag_frame }, ]

	// fragment bias
	//std::string fragment_bias_strategy_;
	utility::vector1<std::string> fragment_bias_strategies_;
	utility::vector1<numeric::random::WeightedSampler> frag_bias_;
	core::Real temp_;

	// FragmentBiasAssigner
	bool cumulate_prob_, exclude_residues_, include_residues_;
	std::set<core::Size> user_pos_;
	std::set<core::Size> residues_to_exclude_;
	std::set<core::Size> residues_to_include_;

	// dump pdb right before changing anything
	bool dump_pdb_;
	std::string dump_pdb_tag_;

	// for auto strategy
	core::Real automode_scorecut_;
	// for rama, geometry strategy
	core::Real score_threshold_;

	// control residue window size to refine
	// for modes: user, auto, rama, geometry, and "include_residues"
	core::Size rsd_wdw_to_refine_;

	// freeze floppy tails: when see jump at the end, how many residues are you chewing back to disallow fragment stitching moves, default=3
	int wdw_to_freeze_;
	bool freeze_endpoints_;

	// selection bias
	std::string selection_bias_;

	// reference model
	core::pose::Pose ref_model_;
	core::Real ref_cst_weight_;
	core::Real ref_cst_maxdist_;
	bool input_as_ref_;
	bool restore_csts_;
	bool debug_;
	bool fullatom_, bbmove_, recover_low_;
	char force_ss_;
	protocols::loops::LoopsOP loops_;

	// scorefunctions
	core::scoring::ScoreFunctionOP scorefxn_, fa_scorefxn_, mc_scorefxn_;  // mc_scorefxn allows us to minimize and eval with different scorefxns

	// ligands
	bool add_hetatm_;
	core::Real hetatm_self_cst_weight_;
	core::Real hetatm_prot_cst_weight_;

}; //class CartesianSampler

} // hybridize
//} // comparative_modeling
} // protocols

#endif
