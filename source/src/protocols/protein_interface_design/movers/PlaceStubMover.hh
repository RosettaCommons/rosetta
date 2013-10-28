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
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_protein_interface_design_movers_PlaceStubMover_hh
#define INCLUDED_protocols_protein_interface_design_movers_PlaceStubMover_hh
#include <protocols/protein_interface_design/movers/PlaceStubMover.fwd.hh>

// Project Headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <basic/datacache/DataMapObj.hh>
// C++ headers

// Unit headers
#include <basic/datacache/DataMap.fwd.hh>

#include <core/kinematics/FoldTree.fwd.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>



namespace protocols {
namespace protein_interface_design {
namespace movers {

/// @brief choose a stub based on mc sampling, and place it on the pose.
/// Iterates over stubs until one matches criteria.
class PlaceStubMover : public simple_moves::DesignRepackMover
{
public:
	// used to define pairs of design/repack movers with an associated bool that
	// determines whether to use the stub-based foldtree with cuts etc., or to
	// use the default foldtree
	typedef std::pair< simple_moves::DesignRepackMoverOP, bool > DesignMoverFoldTreePair;
	typedef std::pair< simple_moves::DesignRepackMoverOP, core::Real > DesignMoverRealPair;
public:
	PlaceStubMover();

	PlaceStubMover(
		protocols::hotspot_hashing::HotspotStubSetOP stub_set,
		core::Real score_threshold,
		core::Size const host_chain,
		protocols::filters::FilterOP final_filter,
		bool const hurry=false,
		bool const triage_positions=true,
		core::Real stub_energy_threshold = 1.0
	);
	virtual ~PlaceStubMover();

	protocols::moves::MoverOP clone() const;
	protocols::moves::MoverOP fresh_instance() const {
		 return protocols::moves::MoverOP( new PlaceStubMover );
	}

	void stub_minimize_movers( utility::vector1< DesignMoverRealPair > const & dmrp ) {
		 stub_minimize_movers_ = dmrp;
	}

	utility::vector1< DesignMoverRealPair > const&
	stub_minimize_movers() const
	{ return stub_minimize_movers_; }

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual void parse_my_tag( utility::tag::TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

private: // member functions
	void
	place_stub( core::pose::Pose & pose, core::conformation::Residue const res_stub, core::Size const res_num );

	void stub_based_atom_tree( core::pose::Pose & pose, core::conformation::Residue const res_stub, core::Real const cst_sdev );

	bool SelectStubIteratively( protocols::hotspot_hashing::HotspotStubSet::Hs_vec::const_iterator stub_it );

	bool StubMinimize( core::pose::Pose & pose,
		protocols::hotspot_hashing::HotspotStubCOP = NULL,
		core::Size const host_res = 0,
		bool const hurry = false );
	void refresh_coordinate_constraints( core::pose::Pose & pose, core::Real const sdev );

	/// @brief resets pose's constraints upon stub failure
	/// remove constraints if they exist. To be used on failure
	void cst_cleanup( core::pose::Pose & pose );
	/// clean everything before exiting
	void final_cleanup( core::pose::Pose & pose );

private: // data members
	/// maximum bonus_value for accepting a stub
	core::Real score_threshold_; // dflt 0.0
	///where is the stub to be placed
	core::Size host_chain_; // dflt 2
	/// stub set we're choosing from
	protocols::hotspot_hashing::HotspotStubSetOP stub_set_;
	/// or use fold tree rb jump.
	bool add_constraints_; /// dflt false
	/// what std to use for coordinate cst in each design mover
	utility::vector1< core::Real > coord_cst_std_;
	/// Locations where stubs are not allowed to be placed
	/// @note This is similar to DesignRepackMover::prevent_repacking() but
	/// we want to allow repacking at steps after stub placement.
	utility::vector1< core::Real > disallowed_host_pos_;
	/// movers for stub minimization, vector of pairs of movers and whether or not
	/// to apply bb constraints during the mover
	utility::vector1< DesignMoverRealPair > stub_minimize_movers_;
	///utility::vector1< core::Real > coord_cst_std_stub_minimize_;
	utility::vector1< DesignMoverFoldTreePair > design_movers_;
	///immediately after placement + minimization
	protocols::filters::FilterOP after_placement_filter_;
	/// a filter at the last stage of placement. Defaults to TrueFilter
	protocols::filters::FilterOP final_filter_;
	core::scoring::constraints::HarmonicFuncOP coord_cst_func_;
	/// where stubs were placed and whether they use constraints. vector is
	/// necessary to maintain the order of the placed stubs
	utility::vector1< std::pair< core::Size, bool > > placed_stubs_;

	/// this is the foldtree with which we come into the first PlaceStubMover.
	/// It will be passed between movers for the user to decide that it should be used.
	core::kinematics::FoldTreeOP default_fold_tree_;
	/// saves the coordinate constraints that this mover has associated with the
	/// pose. At the end of the run, they should be removed
	core::scoring::constraints::ConstraintCOPs curr_coordinate_constraints_;
	bool leave_coord_csts_after_placement_;//defaults to false
	core::Real post_placement_sdev_;//If coord_csts are left on, what should the sdev be? defaults to 1.0

	core::scoring::constraints::ConstraintCOPs previous_coordinate_constraints_;
	core::scoring::constraints::ConstraintCOPs saved_bb_constraints_;
	utility::vector1< core::Size > saved_prevent_repacking_;
	utility::vector1< std::pair< core::Size, bool > > saved_placed_stubs_;

	/// use inverse rotamers to place the scaffold
	bool place_scaffold_; //dflt false
	/// the maximum distance for a stub to be considered a a neighbour to a host residue
	core::Real max_cb_cb_dist_; // dflt 4.0
	/// Should we speed up StubMinimize at the expense of accuracy?
	bool hurry_; //dflt true
	/// Should we triage host positions based on smart criteria, such as current position?
	/// Leave false if partners have not already been docked.
	bool triage_positions_; //dflt true
	/// Maximum per-residue energy the placed stub can have to be saved
	core::Real stub_energy_threshold_;
	/// This task factory is used by placed stubs to add prevent repacking instructions to design movers that
	/// would be invoked *after* PlaceStub finishes successfully. It's a way to communicate the placements that
	/// should not be changed to design movers down stream. Note that this is a different implementation than
	/// the one used in PlaceSimultaneously.
	core::pack::task::TaskFactoryOP residue_level_tasks_for_placed_hotspots_;
  utility::pointer::owning_ptr< basic::datacache::DataMapObj< utility::vector1< core::Size > > > residue_numbers_; /// dflt NULL; a vector of residue numbers placed on the basic::datacache::DataMap which specifies all the placed residues. Useful to communicate between movers and filters, without the pesky NotifyMovers strategy
  std::string user_defined_name_; // reserved for keeping the name of the current stub placement mover
};

} //movers
} //protein_interface_design
} //protocols

#endif /*INCLUDED_protocols_protein_interface_design_movers_PlaceStubMover_HH*/

