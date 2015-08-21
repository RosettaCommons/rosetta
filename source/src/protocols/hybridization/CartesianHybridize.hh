// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Yifan Song
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_hybridization_CartesianHybridize_hh
#define INCLUDED_protocols_hybridization_CartesianHybridize_hh

#include <protocols/hybridization/InsertChunkMover.hh>
#include <protocols/hybridization/CartesianHybridize.fwd.hh>

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
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

using namespace core;
using namespace protocols::moves;
using namespace protocols::loops;

class CartesianHybridize: public protocols::moves::Mover {
public:
	CartesianHybridize();

	CartesianHybridize(
		utility::vector1 < core::pose::PoseOP > const & templates_in,
		utility::vector1 < core::Real > const & template_wts_in,
		utility::vector1 < protocols::loops::Loops > const & template_chunks_in,
		utility::vector1 < protocols::loops::Loops > const & template_contigs_in,
		core::fragment::FragSetOP fragments9_in );

	// initialize options to defaults
	void init();

	// run the protocol
	void apply(core::pose::Pose & pose);

	// set the centroid scorefunctions
	void set_scorefunction(core::scoring::ScoreFunctionOP scorefxn_in);
	void set_min_scorefunction( core::scoring::ScoreFunctionOP scorefxn_in );
	void set_pack_scorefunction( core::scoring::ScoreFunctionOP scorefxn_in );

	// set options
	void set_increase_cycles(core::Real increase_cycles_in) { increase_cycles_=increase_cycles_in; }
	void set_no_global_frame(bool no_global_frame_in) { no_global_frame_=no_global_frame_in; }
	void set_linmin_only(bool linmin_only_in) { linmin_only_=linmin_only_in; }
	void set_cartfrag_overlap(core::Size cartfrag_overlap_in) { cartfrag_overlap_=cartfrag_overlap_in; }
	void set_seqfrags_only(bool seqfrags_only_in) { seqfrags_only_=seqfrags_only_in; }
	void set_skip_long_min(bool skip_long_min_in) { skip_long_min_=skip_long_min_in; }
	void set_cenrot(bool cenrot_in) { cenrot_=cenrot_in; }
	void set_temperature(core::Real temp_in) { temperature_ = temp_in; }
	void set_fragment_probs(core::Real prob, core::Real randprob) {
		fragprob_ = prob;
		randfragprob_ = randprob;
	}

	void set_max_insertion(int max_in) { max_contig_insertion_ = max_in; }

	void set_per_residue_controls(
		utility::vector1<bool> const &residue_sample_template_in,
		utility::vector1<bool> const &residue_sample_abinitio_in) {
		residue_sample_template_ = residue_sample_template_in;
		residue_sample_abinitio_ = residue_sample_abinitio_in;
	}

	std::string get_name() const { return "CartesianHybridize"; }

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

protected:
	// apply a homologue fragment
	void apply_frag( core::pose::Pose &pose, core::pose::Pose &templ, protocols::loops::Loop &frag, bool superpose=true);

	// apply a sequence fragment
	void apply_frame( core::pose::Pose & pose, core::fragment::Frame &frame );

private:
	// parser
	bool align_templates_to_pose_;

	// parameters
	core::Real increase_cycles_;
	core::Size ncycles_, cartfrag_overlap_;
	bool no_global_frame_, linmin_only_;
	bool seqfrags_only_;
	bool skip_long_min_;
	bool cenrot_;
	core::Real temperature_;
	int max_contig_insertion_;  // don't insert contigs larger than this size
	core::Real fragprob_;
	core::Real randfragprob_;

	// fragments
	utility::vector1 < core::pose::PoseOP > templates_;
	utility::vector1 < core::Real > template_wts_;
	utility::vector1 < protocols::loops::Loops > template_contigs_;
	core::fragment::FragSetOP fragments9_;
	boost::unordered_map<core::Size, core::fragment::Frame> library_;

	// per-residue controls
	utility::vector1<bool> residue_sample_template_; // using template fragments
	utility::vector1<bool> residue_sample_abinitio_; // using torsion-based ab initio fragments

	// scorefunctions
	core::scoring::ScoreFunctionOP lowres_scorefxn_, min_scorefxn_, cenrot_repack_scorefxn_;

}; //class CartesianHybridize

} // hybridize
//} // comparative_modeling
} // protocols

#endif
