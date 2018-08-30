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


#ifndef INCLUDED_protocols_loop_grower_FragmentExtension_hh
#define INCLUDED_protocols_loop_grower_FragmentExtension_hh

#include <protocols/loop_grower/FragmentExtension.fwd.hh>

#include <iostream>
#include <fstream>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <map>

#include <core/scoring/ScoreFunction.hh>
#include <core/sequence/Sequence.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loop_grower/LoopGrower.hh>

#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>
#include <numeric/random/WeightedSampler.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>
#include <core/id/SequenceMapping.hh>


#include <utility/tag/Tag.fwd.hh>



namespace protocols {
namespace loop_grower {

class FragmentExtension: public protocols::moves::Mover {
public:
	FragmentExtension();

	// run the protocol
	void apply(core::pose::Pose & pose);

	std::string get_name() const { return "FragmentExtension"; }

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const & ,
		moves::Movers_map const & ,
		core::pose::Pose const & pose );

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	void cluster_fragments( utility::vector1<core::fragment::FragSetOP> & fragments, core::Real fragfilter );

	utility::vector1<protocols::loops::Loop>
	get_unaligned( core::id::SequenceMapping const & seqmap ) const;

	static
	std::string
	mover_name();


private:
	// build parameters
	core::Size beamwidth_, fragtrials_, pack_min_cycles_, direction_, maxloopsize_, master_beam_width_, rmswindow_, steps_, parallelcount_, montecarlorounds_, sheetcriteria_;
	bool debug_, dumpbeam_, dumpfinalbeam_, dumpfragments_, dumprms_, dumperrors_, minimize_, nativegrow_,
		greedy_, read_from_file_, assign_incomplete_, parametercheck_, cenrot_, writebeams_, readbeams_, clustercheck_, rescorebeams_, writelps_, fafilter_, samplesheets_,
		trackfragments_, filterprevious_, famin_, rephasemap_, checksymm_, asymmdump_, continuous_sheets_, auto_stop_;
	core::Real chainbreak_, continuous_weight_, rmscutoff_, fragcluster_, beamscorecutoff_, fa_bonus_, window_dens_weight_, master_beam_cutoff_, sheetbonus_, sheet_tolerance_, sc_scale_;
	core::scoring::ScoreFunctionOP sf_, cen_sf_, cenrot_sf_;
	int looporder_;

	//native
	core::pose::PoseOP native_;

	// sequences
	core::sequence::SequenceOP fullength_seq_;
	utility::vector1< core::Size > cbreaks_;
	std::string fastaname_, storedbeams_, filterbeams_, coordfile_, skeleton_file_;
	core::sequence::SequenceOP fragseq;
	// fragments
	utility::vector1<core::Size> frag_sizes_;
	core::Size nfrags_ ,fragmelt_, minmelt_;
	utility::vector1<core::fragment::FragSetOP> fragments_;
};

class LoopComparator {
public:
	LoopComparator(): one_body_score_(0){}

	void
	apply(core::pose::Pose & pose);

	void
	push(LoopPartialSolutionStore const solutionset){
		solutionsets_.push_back(solutionset);
	}
	void
	fill_pose(core::pose::Pose & pose);

	void set_scores( core::pose::Pose & pose );

	void set_density_sf( core::scoring::ScoreFunctionOP sf_dens ) {
		sf_dens_ = sf_dens;
	}
	void set_vdw_sf( core::scoring::ScoreFunctionOP sf_vdw ){
		sf_vdw_ = sf_vdw;
	}
	void set_sequence( core::sequence::SequenceOP fullength_seq){
		fullength_seq_ = fullength_seq;
	}
	void set_loops( utility::vector1< protocols::loops::Loop > loops ){
		loops_ = loops;
	}

	void set_read_from_file(bool read_from_file) {read_from_file_ = read_from_file;}

	void set_assign_incomplete(bool assign_incomplete) {assign_incomplete_ = assign_incomplete;}

	core::Real get_one_body_score(){ return one_body_score_; }

	//Read the "lpsfile.txt" file to load all the beams
	void read_from_file();

private:
	utility::vector1<LoopPartialSolutionStore> solutionsets_;
	core::scoring::ScoreFunctionOP sf_dens_, sf_vdw_;
	utility::vector1< utility::vector1<core::Real> > one_body_energies_;
	std::map< std::pair< core::Size,core::Size >,utility::vector1< utility::vector1< core::Real> > > two_body_energies_map_;
	core::sequence::SequenceOP fullength_seq_;
	utility::vector1< protocols::loops::Loop > loops_;
	core::Real one_body_score_;
	bool read_from_file_, assign_incomplete_;


};
} // hybridize
} // protocols

#endif
