// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file WaterBoxMover.hh

#ifndef INCLUDED_protocols_moves_WaterBoxMover_hh
#define INCLUDED_protocols_moves_WaterBoxMover_hh

#include <protocols/simple_moves/WaterBoxMover.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/annealer/FixbbPwatSimAnnealer.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/thread/ReadWriteMutex.hh>

#include <numeric/xyzVector.hh>

// C++ Headers
#include <string>
#include <unordered_map>

namespace protocols {
namespace simple_moves {

// typedef for convenience
typedef utility::vector1< utility::vector1< core::Vector > > PWatRotamerCloud;

// store info on rotamer generation from backbone residues
struct WaterRot {
	WaterRot( numeric::xyzVector< core::Real > const &coords, std::string aa, std::string aatm, std::string abase1, std::string abase2) :
		coords_(coords), aa_(aa), aatm_(aatm), abase1_(abase1), abase2_ (abase2)
	{ }

	numeric::xyzVector< core::Real > coords_;
	std::string aa_, aatm_, abase1_, abase2_;
};


//
// a database storing info on backbone generation
class WaterRotsDB {
private:
	utility::vector1< WaterRot > protein_rots_,ligand_rots_;
	bool init_;

public:
	WaterRotsDB(): init_(false) {};

	bool
	is_initialized() {
		return init_;
	}

	core::Size
	n_protein_rots() {
		return protein_rots_.size();
	}

	core::Size
	n_ligand_rots() {
		return ligand_rots_.size();
	}

	void
	initialize();

	WaterRot const &
	protein_rot(core::Size i) {
		return protein_rots_[i];
	}

	WaterRot const &
	ligand_rot(core::Size i) {
		return ligand_rots_[i];
	}
};

//
// a struct/class pair for performing fast clustering of water points
struct AtomHashNode {
	core::Size idx, residx;
	core::Vector xyz;
};

class AtomHash {
public:
	AtomHash() {
		bindis_ = 4.0;
	}

	AtomHash( core::Real bindis ) {
		bindis_ = bindis;
	}

	// adds a single water to the DB
	// give it a unique index
	void
	add_point( core::Size residx, core::Vector xyz );

	// adds a single water to the DB
	// give it a unique index
	core::Size
	get_neighborcount( core::Vector xyz, core::Real dlim );

	// removes all points and returns clusters
	// 3 different distance criteria checks:
	//    1) identify i,j pairs with different resids and d<cut1 (1.5A)
	//    2) throw out pairs w.i cut2 of each other (0.35A)
	//    3) group into "rotamer clouds" w/i cut3 of each other (3A)
	// this function carries out 1
	void
	trim_to_heterogeneous_clusters(core::Real dlim);

	// this function carries out 2&3 (and clears the array)
	void
	reset_and_get_clusters(
		utility::vector1< utility::vector1<core::Vector> > &clusters,
		core::Real dredundant,
		core::Real dcluster
	);

	int
	npoints() { return nodes_.size(); }

private:
	std::unordered_multimap< int, AtomHashNode > nodes_;
	core::Real bindis_;
};


//
// A class that solvates a pose
class WaterBoxMover : public protocols::moves::Mover {
public:
	WaterBoxMover() : Mover("WaterBoxMover") {
		init();
	}

	WaterBoxMover(core::scoring::ScoreFunctionOP sfin): Mover("WaterBoxMover"), sf_ (sfin) {
		init();
	}

	void init();

	static std::string mover_name();

	virtual void apply( Pose & pose );

	// This function deletes waters from a pose (all or virtualized-only)
	void
	delete_waters( core::pose::Pose & pose, bool virt_only );

	// Build one-sided water rotamer clouds
	void
	build_backbone_rotamer_clouds(
		core::pose::Pose const & pose,
		PWatRotamerCloud & new_rotamers
	);

	// Build rotamer-overlap water rotamer clouds
	void
	build_lkboverlap_rotamer_clouds(
		core::pose::Pose const & pose,
		PWatRotamerCloud & new_rotamers
	);

	void
	attach_rotamer_clouds_to_pose_and_rotset(
		core::pose::Pose &pose,
		PWatRotamerCloud &watercloud
	);

	/// setup the initial packing step
	void
	setup_pack( Pose & pose );

	// run packing step
	void
	run_pack( Pose & pose, utility::vector0< int > rot_to_pack, utility::vector1< core::pack::annealer::PointDwell > & all_rot );

	// update packer task to allow packing of waters
	core::pack::task::PackerTaskOP
	update_packer_task(
		core::pose::Pose const & pose,
		core::pack::task::PackerTaskCOP & packer_task
	);

	/// cluster water clouds
	void
	cluster_rotset( utility::vector1< core::pack::annealer::PointDwell > &rotset );

	/// find closest protein residue
	Size
	find_closest( Pose const & pose, core::Vector Ocoord );

	// setters
	void set_lkb_overlap_dist( core::Real val ) { lkb_overlap_dist_=val; }
	void set_lkb_clust_rad( core::Real val ) { lkb_clust_rad_=val; }
	void set_lkb_rotset_radius( core::Real val ) { lkb_rotset_radius_=val; }
	void set_dwell_cutoff( core::Real val ) { dwell_cutoff_=val; }
	void set_clust_radius( core::Real val ) { clust_radius_=val; }
	void set_clust_cutoff( core::Real val ) { clust_cutoff_=val; }
	void set_watlim_scale( core::Real val ) { watlim_scale_=val; }
	void set_gen_fixed( bool val ) { gen_fixed_=val; }
	void set_taskop( core::pack::task::PackerTaskCOP val ) { task_=val; }

	// RS stuff
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		filters::Filters_map const &,
		moves::Movers_map const &,
		core::pose::Pose const & pose );

	static
	utility::tag::XMLSchemaComplexTypeGeneratorOP
	define_water_box_mover_schema();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	// helper function to add a point water to a pose anchored to 'resnum'
	void add_point_water( Pose & pose, core::Vector O, core::Size resnum );

	// helper function: report native recovery stats
	void get_water_recovery( core::pose::Pose const & pose, bool incl_vrt=false );

private:
	// scorefunction
	core::scoring::ScoreFunctionOP sf_;

	// packertask
	core::pack::task::PackerTaskCOP task_;
	core::pack::task::TaskFactoryCOP task_factory_;

	// rotamer data
	core::pack::rotamer_set::RotamerSetsOP rotamer_sets_;
	core::pack::interaction_graph::AnnealableGraphBaseOP ig_;

	// options
	std::string mode_;  // main option: the mode we run in
	core::Real lkb_overlap_dist_, lkb_clust_rad_, lkb_rotset_radius_; // lkb sidechain water overlap, redunant, and rotamer cutoffs
	core::Real dwell_cutoff_; // dwell cutoff after packing, before clustering
	core::Real clust_radius_, clust_cutoff_; // clustering parameters after packing
	core::Real watlim_scale_; // limit number of waters to an absolute factor (scaled by #protein res)
	bool gen_fixed_;

	core::pose::PoseOP native_; // native (if given)

	// the single DB instance and read lock
#ifdef MULTI_THREADED
	static utility::thread::ReadWriteMutex db_mutex_;
#endif
	static WaterRotsDB water_rots_db_;
};

} // moves
} // protocols

#endif
