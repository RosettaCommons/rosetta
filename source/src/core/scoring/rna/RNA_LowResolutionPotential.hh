// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_LowResolutionPotential.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_rna_RNA_LowResolutionPotential_hh
#define INCLUDED_core_scoring_rna_RNA_LowResolutionPotential_hh

#include <core/types.hh>

// Unit headers
#include <core/scoring/rna/RNA_LowResolutionPotential.fwd.hh>
#include <core/scoring/rna/RNA_RawBaseBaseInfo.hh>

// Package headers
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>

// Project headers
#include <core/id/AtomID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>
#include <ObjexxFCL/FArray5D.hh>

namespace core {
namespace scoring {
namespace rna {

////////////////////////////////////////////////////////////////////////////////////////////////////

class RNA_LowResolutionPotential : public utility::pointer::ReferenceCount {

public:
	RNA_LowResolutionPotential();


	void
	update_rna_centroid_info(
		pose::Pose & pose
	) const;

	void
	update_rna_base_base_interactions(
		pose::Pose & pose
	) const;

	void
	update_rna_base_pair_list(
		pose::Pose & pose
	) const;

	void
 	finalize( pose::Pose & pose ) const;

	Real
	rna_base_backbone_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Vector const & centroid1,
		Vector const & centroid2,
		core::kinematics::Stub const & stub1,
		core::kinematics::Stub const & stub2 ) const;

	Real
	rna_backbone_backbone_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2 ) const;

	Real
	rna_repulsive_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2 ) const;

	// void initialize_atom_numbers_for_backbone_score_calculations( pose::Pose & pose ) const; // Deprecated.Commenting out to make Python bindings compile.

	Real
	get_rna_stack_score(
		Distance const x, // used in fading
		Distance const y, // used in fading
		Distance const z, // used in fading
		Real & deriv_x = dummy_deriv,
		Real & deriv_y = dummy_deriv,
		Real & deriv_z = dummy_deriv ) const;

	void
	eval_atom_derivative_base_base(
		id::AtomID const & atom_id,
 		pose::Pose const & pose,
		scoring::EnergyMap const & weights,
 		Vector & F1,
 		Vector & F2 ) const;

	void
	eval_atom_derivative_rna_base_backbone(
		id::AtomID const & atom_id,
 		pose::Pose const & pose,
 		Vector & F1,
 		Vector & F2 ) const;

	void
	eval_atom_derivative_rna_backbone_backbone(
		id::AtomID const & atom_id,
 		pose::Pose const & pose,
 		Vector & F1,
 		Vector & F2 ) const;

	void
	eval_atom_derivative_rna_repulsive(
		id::AtomID const & atom_id,
 		pose::Pose const & pose,
 		Vector & F1,
 		Vector & F2 ) const;

	bool
	check_clear_for_stacking(
	 pose::Pose & pose,
	 Size const & i,
	 int const & sign	) const;

	bool
	check_forming_base_pair(
		pose::Pose & pose,
		Size const & i,
		Size const & j ) const;

	void
	eval_rna_base_pair_energy(
		rna::RNA_RawBaseBaseInfo & rna_raw_base_base_info,
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		Vector const & centroid1,
		Vector const & centroid2,
		core::kinematics::Stub const & stub1,
		core::kinematics::Stub const & stub2 ) const;

	void
	more_precise_base_pair_classification( bool const & value ){ more_precise_base_pair_classification_ = value; }

	Distance
	base_backbone_distance_cutoff() const { return base_backbone_distance_cutoff_;}

	Distance
	base_backbone_z_cutoff() const { return base_backbone_z_cutoff_;}

	Distance
	base_backbone_rho_cutoff() const { return base_backbone_rho_cutoff_;}

	void
	get_zeta_cutoff(
									conformation::Residue const & res_i,
									Real & zeta_hoogsteen_cutoff,
									Real & zeta_sugar_cutoff
									) const;

	bool
	check_for_base_neighbor(
    conformation::Residue const & rsd1,
		Vector const & heavy_atom_j,
		Real & atom_cutoff_weight ) const;

private:

	void
	eval_rna_base_pair_energy_one_way(
		rna::RNA_RawBaseBaseInfo & rna_raw_base_base_info,
		conformation::Residue const & res_i,
		conformation::Residue const & res_j,
		pose::Pose const & pose,
		Vector const & centroid1,
		Vector const & centroid2,
		core::kinematics::Stub const & stub1,
		core::kinematics::Stub const & stub2 ) const;

	Real
	rna_base_backbone_pair_energy_one_way(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		Vector const & centroid1,
		core::kinematics::Stub const & stub1 ) const;

	Size
	find_backbone_oxygen_atom(
     ObjexxFCL::FArray2D < Size > const & atom_number_for_backbone_score_calculations,
		 Size const & i,
		 Size const & atom_num_i ) const;


	Size
	find_backbone_oxygen_atom(
		 Size const & atom_num_i ) const;

private:

	void initialize_rna_basepair_xy();
	void initialize_rna_axis();
	void initialize_rna_stagger();
	void initialize_RNA_backbone_oxygen_atoms();
	void initialize_atom_numbers_for_backbone_score_calculations();
	void initialize_rna_base_backbone_xy();
	void initialize_rna_backbone_backbone_weights();
	void initialize_rna_backbone_backbone();
	void initialize_rna_repulsive_weights();
	void initialize_more_precise_base_pair_cutoffs();

	void
	fill_atom_numbers_for_backbone_oxygens( core::chemical::ResidueTypeSetCAP & rsd_set, core::chemical::AA const & aa );

	bool
	check_atom_numbers_for_backbone_oxygens( core::chemical::ResidueTypeSetCAP & rsd_set, core::chemical::AA const & aa ) const;


	Real get_rna_axis_score(
   Real const cos_theta,
	 Real & deriv = dummy_deriv ) const;

	Real get_rna_stagger_score( Distance const height, Real & deriv = dummy_deriv ) const;

	Real get_rna_basepair_xy(
    Distance const x,
		Distance const y,
		Distance const z, // z is used for fading.
		Real const cos_theta,
		conformation::Residue const & res_i,
		conformation::Residue const & res_j,
		bool const deriv_check = true,
		Real & deriv_x = dummy_deriv,
		Real & deriv_y = dummy_deriv,
		Real & deriv_z = dummy_deriv ) const;

	Real
	get_rna_backbone_backbone_score(
	  Distance const & r,
		Size const & atom_num_j_bin,
		Real & deriv = dummy_deriv ) const;

	Real
	get_rna_repulsive_score(
   Distance const & r,
	 Size const & atom_num_j_bin,
	 Real & deriv = dummy_deriv ) const;


	Real
	rna_backbone_backbone_pair_energy_one_way(
			conformation::Residue const & rsd1,
			conformation::Residue const & rsd2
 ) const;

	Real
	rna_repulsive_pair_energy_one_way(
			conformation::Residue const & rsd1,
			conformation::Residue const & rsd2 ) const;

	Real
	get_rna_base_backbone_xy(
    Distance const x,
		Distance const y,
		Distance const z,
		conformation::Residue const & res_i,
		Size const & atom_num_j_bin,
		bool const deriv_check = false,
		Real & deriv_x = dummy_deriv,
		Real & deriv_y = dummy_deriv,
		Real & deriv_z = dummy_deriv	) const;

	Real
	get_rna_backbone_backbone_xy( Distance const x,
													 Distance const y,
													 conformation::Residue const & res_i,
													 Size const & atom_num_j_bin ) const;

	void
	setup_precise_zeta_cutoffs( chemical::AA const & na_rad,
															std::string const & hoogsteen_cutoff_atom,
															std::string const & sugar_cutoff_atom
															);

private: // data

	// Knowledge-based statistics -- read in from files?
	ObjexxFCL::FArray5D < Real > rna_basepair_xy_;
	ObjexxFCL::FArray1D < Real > rna_stagger_;
	ObjexxFCL::FArray1D < Real > rna_axis_;

	utility::vector1< std::string > RNA_backbone_oxygen_atoms_;
	utility::vector1< Size > atom_numbers_for_backbone_score_calculations_;

	ObjexxFCL::FArray4D < Real > rna_base_backbone_xy_;

	ObjexxFCL::FArray1D < Real > rna_backbone_backbone_weight_;
	ObjexxFCL::FArray1D < Real > rna_backbone_backbone_potential_;

	ObjexxFCL::FArray1D < Real > rna_repulsive_weight_;

	ObjexxFCL::FArray1D < Real > zeta_hoogsteen_cutoff_precise_;
	ObjexxFCL::FArray1D < Real > zeta_sugar_cutoff_precise_;

	//Some parameters for scoring.
	Distance const rna_basepair_radius_cutoff_;
	Distance const rna_basepair_stagger_cutoff_;
  Real const rna_basepair_radius_cutoff2_;
	Distance const basepair_xy_bin_width_;
	Size const basepair_xy_num_bins_;
	Size const basepair_xy_table_size_;
	Distance const basepair_xy_z_fade_zone_;

	Distance const base_stack_min_height_;
	Distance const base_stack_max_height_;
  Distance const base_stack_radius_;
	Real const base_stack_radius2_;
	Distance const base_stack_z_fade_zone_;
	Distance const base_stack_rho_fade_zone_;

	Real const axis_bin_width_;
	Size const axis_num_bins_;

	Size const stagger_num_bins_;
	Distance const stagger_bin_width_;
	Distance const stagger_distance_cutoff_;

	Distance const base_backbone_bin_width_;
	Size const base_backbone_num_bins_;
	Size const base_backbone_table_size_;

	Distance const base_backbone_distance_cutoff_;
	Distance const base_backbone_z_cutoff_;
	Distance const base_backbone_rho_cutoff_;
	Distance const base_backbone_atom_dist_cutoff_;
	Distance const base_backbone_z_fade_zone_;
	Distance const base_backbone_rho_fade_zone_;
	bool const base_backbone_check_atom_neighbor_;

	Real const backbone_backbone_bin_width_;
	Real const backbone_backbone_distance_cutoff_;
	Size const backbone_backbone_num_bins_;

	Real const rna_repulsive_max_penalty_;
	Real const rna_repulsive_screen_scale_;
	Real const rna_repulsive_distance_cutoff_;
	Real const rna_repulse_all_;

	Size const num_RNA_base_pair_orientations_;
	Size const num_RNA_backbone_oxygen_atoms_;
	Size const num_RNA_res_types_;

	Size o2prime_index_within_special_backbone_atoms_;
	Size o2p_index_within_special_backbone_atoms_;

	bool const interpolate_;
	bool const fade_;

	bool const rna_verbose_;

	bool more_precise_base_pair_classification_;

public:
	// This must be public for PyRosetta.  It's quite silly.
	static Real dummy_deriv;

};

} //rna
} //scoring
} //core

#endif
