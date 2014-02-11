// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzdes/EnzdesFlexBBProtocol.hh
///
/// @brief
/// @author Florian Richter




#ifndef INCLUDED_protocols_enzdes_EnzdesFlexBBProtocol_hh
#define INCLUDED_protocols_enzdes_EnzdesFlexBBProtocol_hh


#include <protocols/enzdes/EnzdesBaseProtocol.hh>
// AUTO-REMOVED #include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesLoopsFile.hh>

#include <protocols/backrub/BackrubMover.fwd.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.fwd.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>

#include <core/pose/Pose.hh> /// Replace pack_region_ala_pose_ with a PoseOP to remove this header

#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.fwd.hh>

#include <utility/vector1.hh>


namespace protocols{
namespace enzdes{

class EnzdesFlexBBProtocol;
typedef utility::pointer::owning_ptr< EnzdesFlexBBProtocol > EnzdesFlexBBProtocolOP;
typedef utility::pointer::access_ptr< EnzdesFlexBBProtocol const > EnzdesFlexBBProtocolCAP;
typedef utility::pointer::access_ptr< EnzdesFlexBBProtocol > EnzdesFlexBBProtocolAP;

class EnzdesFlexibleRegion;
typedef utility::pointer::owning_ptr< EnzdesFlexibleRegion > EnzdesFlexibleRegionOP;
typedef utility::pointer::owning_ptr< EnzdesFlexibleRegion const > EnzdesFlexibleRegionCOP;

class EnzdesFlexBBProtocol : public protocols::enzdes::EnzdesBaseProtocol
{

public:

	EnzdesFlexBBProtocol();
	~EnzdesFlexBBProtocol();

	void apply( core::pose::Pose & pose);
	virtual std::string get_name() const;

	static void register_options();

	bool
	is_flexible( core::Size seqpos ) const;

	bool
	is_remodelable( core::Size seqpos ) const;

	void
	get_tenA_neighbor_residues(
  	core::pose::Pose const & pose,
  	utility::vector1<bool> & residue_positions
	) const;

	core::pack::task::PackerTaskOP
	modified_task(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & orig_task
	);

	void
	test_flexbb_rotamer_sets(
		core::pose::Pose const & pose,
		core::pack::task::PackerTaskCOP task
	);

	void
	remap_resid(
		core::pose::Pose const & pose,
		core::id::SequenceMapping const & smap
	);

	void
	add_flexible_region(
		core::Size start,
		core::Size end,
		core::pose::Pose const & pose,
		bool clear_existing
	);

	EnzdesFlexibleRegionCOP
	enz_flexible_region(
		core::Size region ) const{
		return flex_regions_[ region ]; }

	toolbox::match_enzdes_util::EnzdesLoopsFileCOP
	enz_loops_file() const {
		return enz_loops_file_; }

	void
	generate_ensemble_for_region(
		core::pose::Pose & pose,
		core::Size region
	);

	bool
	minimize_flexible_region(
		core::pose::Pose & pose,
		core::Size region,
		core::scoring::ScoreFunctionCOP scorefxn,
		std::set< core::Size > const & chi_to_move,
		bool const including_CA_angles,
		core::Real min_tolerance
	);

	void
	generate_alc_ensemble_for_region(
		core::pose::Pose & pose,
		core::Size region
	);

	void
	generate_backrub_ensemble_for_region(
		core::pose::Pose & pose,
		core::Size region
	);

	void
	determine_flexible_regions(
		core::pose::Pose const & pose,
		core::pack::task::PackerTaskCOP task
	);

protected:


	bool
	assemble_next_best_loop_combination(
		core::pose::Pose & pose );

	bool
	hack_assemble_next_loop_combination(
		core::pose::Pose & pose
	);

	bool
	recover_loops_from_file(
		core::pose::Pose const & pose );

	void
	setup_catalytic_residue_minimization_for_region(
		core::pose::Pose const & pose,
		core::Size region );

	utility::vector1< EnzdesFlexibleRegionOP > flex_regions_;

private:

	toolbox::match_enzdes_util::EnzdesLoopsFileCOP enz_loops_file_;

	utility::vector1< core::Size > fragment_counters_;

	core::pose::Pose pack_region_ala_pose_;
	utility::vector1< core::Real > native_fragment_bb_energies_;

	protocols::backrub::BackrubMoverOP brub_mover_;
	protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP kinematic_mover_;

	core::Real mc_kt_low_, mc_kt_high_;
	core::Size brub_min_atoms_;
	core::Size brub_max_atoms_;
	core::Size loop_ensemble_size_;
	core::Size loopgen_trials_;
	utility::vector1< std::string > brub_pivot_atoms_;

	//std::map< core::Size, core::conformation::ResidueCOP > initial_ligand_positions_;

	//stuff for minimizing catalytic residue chis only during loop generation
	bool minimize_cats_;
	core::kinematics::MoveMapOP catmin_movemap_;
	protocols::simple_moves::MinMoverOP catmin_mover_;
	core::scoring::ScoreFunctionOP catmin_sfxn_;

}; //class EnzdesFlexBBProtocol


class EnzdesFlexibleRegion : public core::fragment::Frame {

	typedef std::pair< core::Size, core::Real > SizeRealPair;
	typedef core::fragment::Frame Super;

public:

	EnzdesFlexibleRegion(
		core::Size index_in,
		core::Size start,
		core::Size end,
		core::Size nr_res,
		core::pose::Pose const & pose,
		EnzdesFlexBBProtocolCAP enz_prot
	);

	~EnzdesFlexibleRegion();

	core::Size
	no_ranked_frags() const{
		return frag_designabilities_.size(); }

	bool
	contains_catalytic_res() const;

	bool remodelable() const{
		return remodelable_; }

	core::Size
	remodel_min_length() const{
		return remodel_min_length_; }

	core::Size remodel_max_length() const{
		return remodel_max_length_; }

	toolbox::match_enzdes_util::EnzdesLoopInfoCOP
	enz_loop_info() const;

	core::Size
	index() const {
		return index_; }

	void
	declare_remodelable(
		core::Size min_length,
		core::Size max_length
	);

	core::Real
	deltaE_best( core::Size const frag_rank ) const;


	core::fragment::FragDataOP
	assemble_enzdes_fragdata(
		core::pose::Pose const & pose
	);

	static bool
	compare_SizeRealPairs( SizeRealPair const & first, SizeRealPair const & second ){
		if( first.second < second.second ) return true;
		else return false;
	}

	utility::vector1< core::Size > const &
	positions() const{
		return positions_;
	}

	core::Real
	extract_lig_designability_score(
		core::pose::Pose const & pose,
		core::pack::task::PackerTaskCOP task,
		core::Real & backgroundE
	);

	core::Real
	get_region_mm_bend_score( core::pose::Pose const & pose ) const;

	void
	apply_ranked_fragment(
		core::pose::Pose & pose,
		core::Size frag_rank
	);

	void
	hack_fillup_frag_designabilities();

	void
	sort_ensemble_by_designability(
		core::pose::Pose const & ref_pose,
		core::scoring::ScoreFunctionOP scorefxn,
		core::pack::task::PackerTaskCOP task
	);

	core::Real
	calculate_rotamer_set_design_targets_partition_sum(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunctionCOP scorefxn,
		core::pack::task::PackerTaskCOP task
	) const;


	bool
	examine_new_loopconf(
		core::pose::Pose const & pose,
		core::pose::Pose & template_pose,
		utility::vector1< core::pose::PoseOP > & compare_poses,
		utility::vector1< core::Real > &  rmsd_to_native
	);

	bool
	minimize_region(
		core::pose::Pose & pose,
		core::scoring::ScoreFunctionCOP scorefxn,
		std::set< core::Size > const & chi_to_move,
		bool const including_CA_angles,
		core::Real min_tolerance
	);

	core::Real
	get_region_total_score( core::pose::Pose const & pose ) const;

	bool remap_resid(
		core::pose::Pose const & pose,
		core::id::SequenceMapping const& smap
	);

	void
	scale_target_proximity_to_starting_conformation( core::Real factor ) {
		target_proximity_to_native_conformation_ *= factor;
	}

	void
	set_target_proximity_to_starting_conformation( core::Real proximity ) {
		target_proximity_to_native_conformation_ = proximity;
	}


	void
	scale_target_proximity_to_other_conformations( core::Real factor ) {
		target_proximity_to_other_conformations_ *= factor;
	}


	std::set< core::Size >
	get_10A_neighbors(
		core::pose::Pose const & pose
	) const;

private:

	//index of this region in the overall scheme of things
	core::Size index_;

	utility::vector1< core::Size > positions_;

	//conformation of the native
	core::fragment::FragDataOP native_conf_;

	//the protocol that this class belongs to
	EnzdesFlexBBProtocolCAP enzdes_protocol_;

	std::set< core::Size > const & design_targets_;

	//a list that holds a designability score for each fragment
	//(fragments accessed by their index in the base class FragData list)
	utility::vector1< SizeRealPair > frag_designabilities_;

	core::Real target_proximity_to_native_conformation_;
	core::Real target_proximity_to_other_conformations_;


	bool remodelable_;
	core::Size remodel_min_length_, remodel_max_length_;
	utility::vector1< std::string > desired_remodel_ss_strings_;

}; //class EnzdesFlexibleRegion


} //namespace enzdes
} //namespace protocols




#endif // INCLUDED_protocols_enzdes_EnzdesFixBBProtocol_HH
