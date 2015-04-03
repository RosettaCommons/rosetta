// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/PlaceStubMover.hh
/// @brief definition of classes for grafting hotspots into a pose
/// @author Sarel Fleishman (sarelf@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_PlaceSimultaneouslyMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_PlaceSimultaneouslyMover_hh

// Project Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// C++ headers

// Unit headers
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/protein_interface_design/movers/PlacementAuctionMover.fwd.hh>
#include <protocols/protein_interface_design/movers/PlacementMinimizationMover.fwd.hh>
#include <protocols/protein_interface_design/filters/StubScoreFilter.fwd.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief choose a stub based on mc sampling, and place it on the pose.
/// Iterates over stubs until one matches criteria.
class PlaceSimultaneouslyMover : public simple_moves::DesignRepackMover
{
public:
	typedef std::pair< protocols::hotspot_hashing::HotspotStubSetOP, std::pair< protocols::hotspot_hashing::HotspotStubOP, core::Size > > StubSetStubPos;
	typedef std::pair< protocols::simple_moves::DesignRepackMoverOP, core::Real > MoverRealPair;
public:
	// used to define pairs of design/repack movers with an associated bool that
	// determines whether to use the stub-based foldtree with cuts etc., or to
	// use the default foldtree
	typedef std::pair< simple_moves::DesignRepackMoverOP, bool > DesignMoverFoldTreePair;
	typedef std::pair< simple_moves::DesignRepackMoverOP, core::Real > DesignMoverRealPair;
	/// ResidueAuction is keyed by energy => we select the residue,stub,stubset combination with the best energy for each stubset,stub combination
	typedef std::pair< protocols::hotspot_hashing::HotspotStubSetOP, protocols::hotspot_hashing::HotspotStubOP > StubsetStubPair;
	typedef std::pair< core::Real, std::pair< core::Size, StubsetStubPair > > ResidueAuctionItem;
	typedef std::multimap< core::Real, std::pair< core::Size, StubsetStubPair > > ResidueAuction;
public:
	PlaceSimultaneouslyMover();
	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const {
		 return protocols::moves::MoverOP( new PlaceSimultaneouslyMover );
	}

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	/// @brief minimize towards clouds of stubs made up of all the stub sets
	/// @brief if bb_cst score is 0 return false, o/w true
	bool minimize_no_bb( core::pose::Pose & pose ) const;
	/// @brief minimize simultaneously towards the placed stubs from each of the stub sets
	/// returns false if one of the filters fails
	void minimize_all( core::pose::Pose & pose, core::Size const minimization_steps ) const;
	/// @brief pair each stub set with a position on the scaffold
	/// @brief if no mutually exclusive matches are found for all sets, return false
	bool pair_sets_with_positions( core::pose::Pose & pose );
/// @brief conducts user-specified design movers. Returns true if the energy per residue filter passes for each of the placed hotspots
	void design( core::pose::Pose & pose );
	/// @brief will be removed
	bool place_stubs( core::pose::Pose & pose ) const;
	virtual void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );
	void final_cleanup( core::pose::Pose & pose );
	/// @brief removes and reinstates coordinate constraints for all placed hotspots according to coord_sdev
	void refresh_coordinate_constraints( core::pose::Pose & pose, core::Real const coord_sdev );
	core::pack::task::PackerTaskOP create_task_for_hotspot_packing( core::pose::Pose const & );
	core::pack::task::PackerTaskOP create_task_for_allhotspot_packing( core::pose::Pose const & );
  void add_coordinatecst_for_hotspot_packing( core::pose::Pose & );
	void stub_sets( utility::vector1< StubSetStubPos > const sets ) { stub_sets_ = sets; }
	void host_chain( core::Size const host_chain ) { host_chain_ = host_chain; }
	virtual ~PlaceSimultaneouslyMover();

private:
	///where is the stub to be placed
	core::Size host_chain_;
	core::pose::Pose saved_pose_;//used to restore the pose if placement fails
	/// stub sets and placed positions
	utility::vector1< StubSetStubPos > stub_sets_;
	utility::vector1< StubSetStubPos > saved_stub_sets_;
	std::map< protocols::hotspot_hashing::HotspotStubSetOP, protocols::filters::FilterOP > stub_set_filters_;

	/// the maximum distance for a stub to be considered a a neighbour to a host residue
	core::Real max_cb_cb_dist_;
	core::Real coor_cst_cutoff_;
	/// Maximum per-residue energy the placed stub can have to be saved
	core::Real stub_energy_threshold_;
	utility::vector1< MoverRealPair  > minimization_movers_;
	core::Size minimization_repeats_before_placement_;//how many times to repeat the minimization procedure after hs placement
	core::Size minimization_repeats_after_placement_;//how many times to repeat the minimization procedure after hs placement
	core::scoring::constraints::ConstraintCOPs saved_bb_constraints_;
	protocols::filters::FilterOP after_placement_filter_, final_filter_;
	utility::vector1< MoverRealPair > design_movers_;
	core::Size explosion_;
	/// This task factory is used by placesimultaneously in two ways: 1) internally to create a packer task for
	/// hotspot placement and user-defined design movers. 2) For design movers downstream of PlaceSimultaneously
	/// so that they are aware of the choices made by PlaceSimultaneously. Note that item 1 is not implemented
	/// in PlaceStub. This means that PlaceSimultaneously *always* uses the task factory, whereas PlaceStub only
	/// uses it if the user specifies TaskAware design movers.
	core::pack::task::TaskFactoryOP residue_level_tasks_for_placed_hotspots_;
	PlacementAuctionMoverOP auction_;
	/// @brief does the pose have any active constraints?
	protocols::protein_interface_design::filters::StubScoreFilterOP stub_score_filter_;
	PlacementMinimizationMoverOP rbstub_minimization_;
//	bool user_defined_auction_; // in which case, auction will not be called in PlaceSim
//	bool user_defined_stub_score_filter_; // in which case, preliminary minimization will not be called in PlaceSim
//	bool user_defined_bbstub_minimization_;
	core::scoring::constraints::ConstraintCOPs saved_coord_constraints_;
};

} //movers
} //protein_interface_design
} //protocols

#endif /*INCLUDED_protocols_protein_interface_design_movers_PlaceSimultaneously_HH*/
