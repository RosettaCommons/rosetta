// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

/// @brief Use LoopHash to generate low-resolution alternate conformations for 
/// the given pose.
/// @details It's possible to limit the amount of movement based on the 
/// secondary structure, such that secondary structure elements either move 
/// less than loops or don't move at all.  The algorithm in pseudo-code:
///
/// - Convert the pose to centroid mode.
///
/// - For num_iterations():
///
/// 	- Randomly pick a window to diversify with LoopHash.
///
/// 	- If \p diversify_loop_only_ is \t true and there are no loops in the 
/// 	   window, skip the rest of the iteration.
///
/// 	- Set the amount of backbone movement to allow (i.e. minimum and maximum 
/// 	   allowed RMSD of backbone segments) based on whether or not the 
/// 	   secondary structure changes within the window.
///
/// 	- Keep all the insertions that pass a filter.
///
/// - Convert the best-scoring diversified pose to full-atom mode and apply it.
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
		std::string const & start_res,
		std::string const & stop_res,
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

	void apply( core::pose::Pose & pose ) override;
	// core::pose::PoseOP get_additional_output();

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;
	protocols::moves::MoverOP clone() const override { return( protocols::moves::MoverOP( new LoopHashDiversifier( *this ) ) ); }
	protocols::moves::MoverOP fresh_instance() const override { return protocols::moves::MoverOP( new LoopHashDiversifier ); }
	~LoopHashDiversifier() override;

	core::Real min_inter_ss_bbrms() const;
	void min_inter_ss_bbrms( core::Real const min_inter_ss_bbrms );
	core::Real max_inter_ss_bbrms() const;
	void max_inter_ss_bbrms( core::Real const max_inter_ss_bbrms );

	core::Real min_intra_ss_bbrms() const;
	void min_intra_ss_bbrms( core::Real const min_intra_ss_bbrms );
	core::Real max_intra_ss_bbrms() const;
	void max_intra_ss_bbrms( core::Real const max_intra_ss_bbrms );

	core::Real min_rms() const;
	void min_rms( core::Real const min_rms );
	core::Real max_rms() const;
	void max_rms( core::Real const max_rms );

	core::Size num_iterations() const;
	void num_iterations( core::Size const num_iterations );

	core::Size num_try_div() const;
	void num_try_div( core::Size const num_try_div );

	utility::vector1< core::Size > loop_sizes() const;
	void add_loop_size( core::Size const loop_size );

	void cenfilter( protocols::filters::FilterOP cenfilter );
	void ranking_cenfilter( protocols::filters::FilterOP filter ){ ranking_cenfilter_ = filter; }
	protocols::filters::FilterOP ranking_cenfilter() const{ return ranking_cenfilter_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:
	// loopsampler
	utility::vector1< core::Size > loop_sizes_;
	LoopHashLibraryOP library_;

	//torsion-angle rmsd for residue windows that span two pieces of secondary structure
	core::Real min_inter_ss_bbrms_, max_inter_ss_bbrms_;

	//torsion-angle rmsd for residue windows that are a single pieces of secondary structure
	core::Real min_intra_ss_bbrms_, max_intra_ss_bbrms_;

	//Anstrom Rmsd cuttofs for loophash generated structures
	core::Real min_rms_, max_rms_;

	//Residues to loophash over
	std::string start_res_, stop_res_;

	//loophash window size & loophash fragment size
	core::Size window_size_;

	//max loophash radius
	core::Size max_radius_;

	//Max models to create in loophash (number of fragment tried is 200x this)
	core::Size max_struct_;

	//Number of loophash runs to execute
	core::Size num_iterations_;

	//Number of loophash runs to try in each execution
	core::Size num_try_div_;

	bool diversify_loop_only_;

	// should we save space and assume structure is ideal?
	bool ideal_, filter_by_phipsi_;

	//cen actually prune decoys based on the filter's apply function.
	protocols::filters::FilterOP cenfilter_, ranking_cenfilter_;

	core::scoring::ScoreFunctionOP scorefxn_cen_cst_, scorefxn_rama_cst_;

	// structure store
	std::vector< std::pair< core::Real, core::io::silent::SilentStructOP > > all_structs_;
};


} // loophash
} // protocols


#endif /*INCLUDED_protocols_loophash_movers_LoopHashDiversifier_HH*/
