// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/movers/LoopHashMoverWrapper.hh
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_loophash_LoopHashMoverWrapper_hh
#define INCLUDED_protocols_loophash_LoopHashMoverWrapper_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/silent/SilentStruct.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/FastRelax.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
// AUTO-REMOVED #include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace loophash {

using core::Real;
using core::Size;
using std::string;

class LoopHashMoverWrapper : public protocols::moves::Mover
{
public:
	LoopHashMoverWrapper();

	void apply( core::pose::Pose & pose );
	core::pose::PoseOP get_additional_output();

	virtual std::string get_name() const;
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	protocols::moves::MoverOP clone() const { return( protocols::moves::MoverOP( new LoopHashMoverWrapper( *this ) ) ); }
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new LoopHashMoverWrapper ); }
	virtual ~LoopHashMoverWrapper();

	Real min_bbrms() const;
	void min_bbrms( Real const min_bbrms );
	Real max_bbrms() const;
	void max_bbrms( Real const max_bbrms );

	Real min_rms() const;
	void min_rms( Real const min_rms );
	Real max_rms() const;
	void max_rms( Real const max_rms );

	Size max_nstruct() const;
	void max_nstruct( Size const max_nstruct );

	utility::vector1< Size > loop_sizes() const;
	void add_loop_size( Size const loop_size );

	void cenfilter( protocols::filters::FilterOP cenfilter );
	void fafilter( protocols::filters::FilterOP fafilter );
	void ranking_fafilter( protocols::filters::FilterOP ranking_fafilter );
	void ranking_cenfilter( protocols::filters::FilterOP filter ){ ranking_cenfilter_ = filter; }
	protocols::filters::FilterOP ranking_cenfilter() const{ return ranking_cenfilter_; }
	void relax_mover( protocols::relax::FastRelaxOP relax_mover );

private:
	// loopsampler
	utility::vector1< Size > loop_sizes_;
	LoopHashLibraryOP library_;
	protocols::relax::FastRelaxOP fastrelax_;         // the relax mover (will be applied in batches)

	// params
	Real min_bbrms_, max_bbrms_, min_rms_, max_rms_;  // loophash centroid generation params
	Size start_res_, stop_res_, max_struct_, max_radius_, max_struct_per_radius_;                       // residues to loophash over
	Size max_nstruct_;                                // maximum number of residues per position to generate
	Size ncentroid_, nfullatom_;                      // number of structures to carry over to subsequent stages
	Size batch_size_;                                 // batch relax batch size
	protocols::filters::FilterOP cenfilter_, ranking_cenfilter_, fafilter_, ranking_fafilter_; // cen/fafilter actually prune decoys based on the filter's apply function. ranking_fa_filter uses the filter's report_sm method to rank and report the best decoys

	// prefiltering
	core::Size nprefilter_;
	core::scoring::ScoreFunctionOP prefilter_scorefxn_;

	// structure store
	bool ideal_, filter_by_phipsi_;  // should we save space and assume structure is ideal?
	std::vector< std::pair< Real, core::io::silent::SilentStructOP > > all_structs_;
	core::Real sample_weight_const_; //dflt 1.0 ; sets the same sample weight throughout the pose
};


} // loophash
} // protocols


#endif /*INCLUDED_protocols_loophash_movers_LoopHashMoverWrapper_HH*/
