// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   SplitAndMixPoseMover.hh
/// @brief  Splits a Pose according to a ResidueSelector. A se gment defined between multiple chains will
///         end up splitted too.
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

#ifndef INCLUDED_protocols_fold_from_loops_SplitAndMixPoseMover_hh
#define INCLUDED_protocols_fold_from_loops_SplitAndMixPoseMover_hh

#include <basic/Tracer.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/hybridization/util.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>
#include <protocols/fold_from_loops/SplitAndMixPoseMover.fwd.hh>

#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

// XSD XRW Includes
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <string>
#include <boost/tuple/tuple.hpp>

namespace protocols {
namespace fold_from_loops {

class SplitAndMixPoseMover : public moves::Mover {

	typedef boost::tuple<Size, Size> SubPoseIDPair;

public:
	/// @brief Empty Constructor
	SplitAndMixPoseMover();
	/// @brief Destructor
	inline ~SplitAndMixPoseMover(){};

	/// @brief ResidueSelector Setter
	inline void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector ){ selector_ = selector; };
	/// @brief ResidueSelector Getter
	inline core::select::residue_selector::ResidueSelectorCOP residue_selector() const { return selector_; };
	/// @brief ResidueRanges Getter
	inline core::select::residue_selector::ResidueRangesOP ranges() const { return ranges_; }
	/// @brief Re-merge order reSet to empty
	inline void reset_order(){ if ( order_.size() > 0 ) order_ = utility::vector1< core::Size >(); };
	/// @brief Re-merge order Setter
	inline void set_order( utility::vector1< core::Size > order ){ order_ = order; };
	/// @brief Re-merge order Setter
	void set_order( std::string order );
	/// @brief Re-merge order Setter from 1 to max
	inline void set_order( core::Size max ){ reset_order(); for ( core::Size i=1; i<=max; ++i ) { order_.push_back( i );} };
	/// @brief Re-merge order Getter
	inline utility::vector1< core::Size > order() const { return order_; };
	/// @brief Merge Chains Setter
	/// If True, this will only work as long as the segments of the same chain are ordered together.
	inline void set_merge_chains( bool merge ){ merge_chains_ = merge; };
	/// @brief Merge Chains Getter
	inline bool merge_chains() const { return merge_chains_; };
	/// @brief Merge Chains Default
	inline static bool default_merge_chains() { return false; };
	// @brief Try All Pairs Setter
	/// If True, at each loop of -nstruct it will generate a new pose combining two different subposes.
	inline void set_try_all_pairs( bool try_all_pairs ){ try_all_pairs_ = try_all_pairs; };
	/// @brief Try All Pairs Getter
	inline bool try_all_pairs() const { return try_all_pairs_; };
	/// @brief Try All Pairs Default
	inline static bool default_try_all_pairs() { return false; };
	// @brief Exclude Consecutive Setter
	/// If True, and try_all_pairs_ True, consecutive pairs are excluded from the list
	inline void set_exclude_consecutive( bool exclude_consecutive ){ exclude_consecutive_ = exclude_consecutive; };
	/// @brief Exclude Consecutive Getter
	inline bool exclude_consecutive() const { return exclude_consecutive_; };
	/// @brief Exclude Consecutive Default
	inline static bool default_exclude_consecutive() { return true; };
	// @brief MaxDistance Setter
	/// If Set to a positive value, that is the max distance between SubPose1-Cterm and SubPose2-Nterm allowed to concatenate them
	inline void set_max_distance( core::Real max_distance ){ max_distance_ = max_distance; };
	/// @brief MaxDistance Getter
	inline core::Real max_distance() const { return max_distance_; };
	/// @brief MaxDistance Default
	inline static core::Real default_max_distance() { return -1; };
	/// @brief GlobalI Getter
	inline core::Size get_global_i() const { return global_i_; };
	/// @brief Get Number of Segments
	inline core::Size count_segments() const { return ranges_->size(); };
	/// @brief Calculate Number of Segments
	core::Size count_segments( core::pose::Pose const & pose );

	/// @brief Apply Mover
	void apply( core::pose::Pose & pose );

	inline std::string get_name() const { return mover_name(); };
	inline static std::string mover_name() { return "SplitAndMixPoseMover"; };
	inline moves::MoverOP clone() const { return moves::MoverOP( new SplitAndMixPoseMover( *this ) ); };
	inline moves::MoverOP fresh_instance() const{ return moves::MoverOP( new SplitAndMixPoseMover() ); };
	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &data,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	void fix_ranges( core::pose::Pose const & pose, core::select::residue_selector::ResidueRanges const & ranges );
	utility::vector1< core::pose::PoseOP > split_pose( core::pose::Pose const & pose );
	core::pose::PoseOP merge_poses( utility::vector1< core::pose::PoseOP > const & pose_list );
	void transfer_conformation( core::pose::Pose & pose, core::pose::Pose const & merged );
	inline void increase_global_count() { if ( global_i_ == pairs_.size() ) { global_i_ = 1;} else { ++global_i_; } };
	void fill_pair_list( core::Size total_subposes );
	void adapt_order_to_pair( utility::vector1< core::pose::PoseOP > const & pose_list );
	bool NC_distance_filter( core::pose::Pose const & npose, core::pose::Pose const & cpose );

private:
	core::select::residue_selector::ResidueSelectorCOP selector_;
	core::select::residue_selector::ResidueRangesOP    ranges_;
	utility::vector1< core::Size >                     order_;
	bool                                               merge_chains_;
	utility::vector1< SubPoseIDPair >                  pairs_;
	bool                                               try_all_pairs_;
	bool                                               exclude_consecutive_;
	core::Real                                         max_distance_;
	core::Size                                         global_i_;
};

// boost::get<0>(*it);
// Size length = boost::get<1>(*it);
}
}

#endif
