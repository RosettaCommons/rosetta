// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/loophash/movers/LoopHashDiversifier.hh
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_loophash_LoopHashDiversifier_hh
#define INCLUDED_protocols_loophash_LoopHashDiversifier_hh

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

#include <utility/vector1.hh>


namespace protocols {
namespace loophash {

using core::Real;
using core::Size;
using std::string;

class LoopHashDiversifier : public protocols::moves::Mover
{
public:
	LoopHashDiversifier();

	LoopHashDiversifier(
		LoopHashLibraryOP library,
		core::Real min_inter_ss_bbrms,
		core::Real max_inter_ss_bbrms,
		core::Real min_intra_ss_bbrms,
		core::Real max_intra_ss_bbrms,
		core::Real min_rms,
		core::Real max_rms,
		core::Size start_res,
		core::Size stop_res,
		core::Size window_size,
		core::Size max_radius,
		core::Size max_struct,
		core::Size num_iterations,
		core::Size num_try_div,
		bool diversify_loop_only,
		bool ideal,
		bool filter_by_phipsi,
		protocols::filters::FilterOP cenfilter,
		protocols::filters::FilterOP ranking_cenfilter,
		core::scoring::ScoreFunctionOP scorefxn_cen_cst,
		core::scoring::ScoreFunctionOP scorefxn_rama_cst
	);

	void apply( core::pose::Pose & pose );
	// core::pose::PoseOP get_additional_output();

	virtual std::string get_name() const;
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	);
	protocols::moves::MoverOP clone() const { return( protocols::moves::MoverOP( new LoopHashDiversifier( *this ) ) ); }
	protocols::moves::MoverOP fresh_instance() const { return protocols::moves::MoverOP( new LoopHashDiversifier ); }
	virtual ~LoopHashDiversifier();

	Real min_inter_ss_bbrms() const;
	void min_inter_ss_bbrms( Real const min_inter_ss_bbrms );
	Real max_inter_ss_bbrms() const;
	void max_inter_ss_bbrms( Real const max_inter_ss_bbrms );

	Real min_intra_ss_bbrms() const;
	void min_intra_ss_bbrms( Real const min_intra_ss_bbrms );
	Real max_intra_ss_bbrms() const;
	void max_intra_ss_bbrms( Real const max_intra_ss_bbrms );

	Real min_rms() const;
	void min_rms( Real const min_rms );
	Real max_rms() const;
	void max_rms( Real const max_rms );

	Size num_iterations() const;
	void num_iterations( Size const num_iterations );

	Size num_try_div() const;
	void num_try_div( Size const num_try_div );

	utility::vector1< Size > loop_sizes() const;
	void add_loop_size( Size const loop_size );

	void cenfilter( protocols::filters::FilterOP cenfilter );
	void ranking_cenfilter( protocols::filters::FilterOP filter ){ ranking_cenfilter_ = filter; }
	protocols::filters::FilterOP ranking_cenfilter() const{ return ranking_cenfilter_; }

private:
	// loopsampler
	utility::vector1< Size > loop_sizes_;
	LoopHashLibraryOP library_;

	//torsion-angle rmsd for residue windows that span two pieces of secondary structure
	Real min_inter_ss_bbrms_, max_inter_ss_bbrms_;

	//torsion-angle rmsd for residue windows that are a single pieces of secondary structure
	Real min_intra_ss_bbrms_, max_intra_ss_bbrms_;

	//Anstrom Rmsd cuttofs for loophash generated structures
	Real min_rms_, max_rms_;

	//Residues to loophash over
	Size start_res_, stop_res_;

	//loophash window size & loophash fragment size
	Size window_size_;

	//max loophash radius
	Size max_radius_;

	//Max models to create in loophash (number of fragment tried is 200x this)
	Size max_struct_;

	//Number of loophash runs to execute
	Size num_iterations_;

	//Number of loophash runs to try in each execution
	Size num_try_div_;

	bool diversify_loop_only_;

	// should we save space and assume structure is ideal?
	bool ideal_, filter_by_phipsi_;

	//cen actually prune decoys based on the filter's apply function.
	protocols::filters::FilterOP cenfilter_, ranking_cenfilter_;

	core::scoring::ScoreFunctionOP scorefxn_cen_cst_, scorefxn_rama_cst_;

	// structure store
	std::vector< std::pair< Real, core::io::silent::SilentStructOP > > all_structs_;
};


} // loophash
} // protocols


#endif /*INCLUDED_protocols_loophash_movers_LoopHashDiversifier_HH*/
