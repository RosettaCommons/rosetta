// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/GridScorer.hh
///
/// @brief  Grid-based scoring for GA ligand docking.  Full grid-based implementation of beta16.
/// @author Hahnbeom Park and Frank DiMaio

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_GridScorer_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_GridScorer_hh

#include <protocols/ligand_docking/GALigandDock/GridScorer.fwd.hh>

#include <protocols/ligand_docking/GALigandDock/LigandConformer.hh>
#include <protocols/ligand_docking/GALigandDock/GridHash3D.hh>
#include <protocols/ligand_docking/GALigandDock/RotamerData.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#include <core/energy_methods/CartesianBondedEnergy.fwd.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/conformation/Atom.hh>
#include <core/optimization/MinimizerMap.fwd.hh>
#include <core/optimization/MinimizerOptions.fwd.hh>
#include <core/scoring/lkball/LK_BallEnergy.fwd.hh>
#include <core/scoring/etable/coulomb/Coulomb.fwd.hh>


#include <core/scoring/ScoreFunction.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <numeric/xyzVector.hh>
#include <numeric/cubic_polynomial.hh>
#include <utility/vector1.hh>


#include <cmath>
#include <chrono>

#include <core/scoring/hbonds/HBondDatabase.fwd.hh> // AUTO IWYU For HBondDatabaseCOP
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

/// @brief Grid representation of scorefunction
class GridScorer {
public:

	GridScorer( core::scoring::ScoreFunctionOP sfxn );
	~GridScorer();

	/// @brief  calulate the bounding box
	void
	prepare_grid(
		core::pose::Pose const &pose,
		utility::vector1< core::Size > const &lig_resnos );

	/// @brief  pass list of residue (acutually rsdtypes) to construct grids for
	void
	get_grid_atomtypes(
		utility::vector1< core::conformation::Residue > const & rsds,
		utility::vector1< bool > const & sconly
	);

	void
	get_grid_all_atomtypes();

	void
	calculate_grid(
		core::pose::Pose const &pose,
		utility::vector1< core::Size > const & lig_resnos,
		utility::vector1< core::Size > const &movingSCs
	);

	/// @brief  calculate energies on the grid
	core::Real
	score( LigandConformer const &lig, bool soft=false);

	/// @brief calculate energies on the grid for initial pool generation
	core::Real
	score_init( LigandConformer const &lig, bool soft=false, core::Real init_dens_weight=100 );

	core::Real
	score_init( core::pose::Pose &pose, LigandConformer const &lig, bool soft, core::Real init_dens_weight=100 );

	/// @brief calculate electron density on the grid
	core::Real
	density_score( LigandConformer const &lig );

	core::Real
	calculate_ligand_density_correlation( int resid, core::conformation::Residue const &rsd, core::pose::Pose const &pose);//, core::Size gene_number);

	//calculates penalty for ligand density correlation for identification
	core::Real
	calculate_density_penalty( core::Size nres, core::conformation::Residue const &lig, core::pose::Pose const &pose, core::Real unpenalized_density );

	core::Real
	calculate_pose_density_correlation( core::pose::Pose const &pose );

	core::Real
	calculate_pocket_density_correlation( core::pose::Pose const &pose );

	/// @brief  gets the clash energy of a point
	core::Real
	point_clash_energy(
		numeric::xyzVector< core::Real > X,
		std::string atype_in="",
		bool soft=false
	);

	/// @brief  gets the solvation energy of a point
	core::Real
	point_energy(
		numeric::xyzVector< core::Real > X,
		std::string atype_in="",
		core::Real weightelec_in=0.0
	);

	/// @brief  get a residue:background energy
	ReweightableRepEnergy
	get_1b_energy(
		core::conformation::Residue const &res_i,
		core::scoring::lkball::LKB_ResidueInfoOP lkbrinfo,
		bool include_bb,
		bool soft=false,
		bool has_density_map=false
	);

	/// @brief  get a residue:residue energy
	ReweightableRepEnergy
	get_2b_energy(
		core::pose::Pose &pose,
		core::conformation::Residue const &res_i,
		core::scoring::lkball::LKB_ResidueInfoOP lkbrinfo_i,
		bool incl_bb_i,
		core::conformation::Residue const &res_j,
		core::scoring::lkball::LKB_ResidueInfoOP lkbrinfo_j,
		bool incl_bb_j,
		bool soft=false
	);

	//core::Size
	//get_hbond_map( core::conformation::Residue const &res_i, core::conformation::Residue const &res_j );
	std::map< std::pair < core::Size, core::Size >, std::vector < core::Size > >
	get_hbond_map( core::pose::Pose pose, core::Size const lig_resno );

	/// @brief  calculate energies on the grid
	core::Real
	score( core::pose::Pose &pose, LigandConformer const &lig, bool soft=false ) ;

	/// @brief  calculate derivatives on the grid
	void
	derivatives(
		core::pose::Pose &pose,
		LigandConformer const &lig,
		core::optimization::MinimizerMap &min_map // pass by non-const ref so we can sum derivs inside
	);

	/// @brief  (1b) dof derivatives
	core::Real
	dof_derivative(
		core::pose::Pose &pose,
		core::optimization::MinimizerMap &min_map,
		core::id::DOF_ID const & dof_id,
		core::id::TorsionID const & torsion_id
	);

	void
	debug_deriv(
		core::pose::Pose &pose,
		LigandConformer const &lig,
		core::optimization::MinimizerMap &min_map
	);

	/// @brief  minimize a ligand conformer
	core::Real
	optimize(
		LigandConformer &lig,
		utility::vector1<core::Real> ramp_schedule,
		utility::vector1< PlaceableRotamers > &placeable_rotdb, // by non-const ref since we update "in place"
		RotamerPairEnergies &rot_energies                       // by non-const ref since we update "in place"
	);

	/// @brief  subroutine: run packing
	core::Real
	packer_loop(
		LigandConformer &lig,
		utility::vector1< PlaceableRotamers > &placeable_rotdb, // by non-const ref since we update "in place"
		RotamerPairEnergies &rotamer_energies                   // by non-const ref since we update "in place"
	);

	/// @brief  subroutine: min loop
	core::Real
	minimizer_loop(
		LigandConformer &lig,
		core::optimization::MinimizerOptions const &minopt
	);

	/// @brief  check to see if a point falls w/i the grid boundaries
	bool
	is_point_in_grid(
		numeric::xyzVector< core::Real > x
	) const;

	/// @detail  check to see if sidechain lies within the grid
	///   "angle buffer" makes sure residues near the edge point toward middle (higer values == stricter, 0 = no filter)
	///   "padding buffer" adds a region on the outside where residues may not move
	bool
	is_residue_in_grid(
		core::conformation::Residue const &res,
		core::Real angle_buffer,
		core::Real padding_buffer
	) const;

	/// @detail check to see if sidechain lies within the grid (ALTERNATE version)
	///   similar to the above, but pass in eigenvalues/vectors that define pocket shape
	bool
	is_residue_in_grid(
		core::conformation::Residue const &res,
		core::Real padding_buffer,
		numeric::xyzVector< core::Real > const &eigval,
		numeric::xyzMatrix< core::Real > const &eigvec
	) const;

	core::Real get_padding() { return bbox_padding_; }
	core::Real get_voxel_spacing() { return voxel_spacing_; }

	void set_voxel_spacing( core::Real voxel_spacing_in ) { voxel_spacing_ = voxel_spacing_in; }
	void set_bbox_padding( core::Real bbox_padding_in ) { bbox_padding_ = bbox_padding_in; }
	void set_hash_gridding( core::Real hash_gridding_in ) { hash_grid_ = hash_gridding_in; }
	void set_hash_subgridding( core::Size hash_subgridding_in ) { hash_subgrid_ = hash_subgridding_in; }

	void set_exact( bool exactin ) { exact_ = exactin; }
	bool get_exact( ) { return exact_; }

	void set_debug( bool debugin ) { debug_ = debugin; }
	bool get_debug( ) { return debug_; }

	void set_grid_dim_with_maxRad( core::Real inval);
	core::Real get_maxRad() {return maxRad_;}

	numeric::xyzVector< core::Real > get_ligand_com() {return lig_com_;}

	void set_force_exact_min( bool setting ) { force_exact_min_ = setting; }

	// protocol stuff
	void set_maxiter_minimize( core::Size setting ) { maxiter_minimize_ = setting; }
	bool get_maxiterminimize( ) { return maxiter_minimize_; }
	void set_packer_cycles( core::Size setting ) { packer_cycles_ = setting; }
	core::Size get_packer_cycles( ) { return packer_cycles_; }
	// smoothing factor on gridded scores only : returns true if the grids were modified
	bool set_smoothing( core::Real setting );
	core::Real get_smoothing( ) { return smoothing_; }

	bool set_elec_scale( core::Real setting );
	core::Real get_elec_scale( ) { return elec_scale_; }

	// scalefactor on fa_rep
	void set_w_rep( core::Real setting ) { w_rep_ = setting; }
	core::Real get_w_rep( ) { return w_rep_; }

	// out of boundary scale
	void set_out_of_bound_e( core::Real setting ) { out_of_bound_e_ = setting; }

	core::scoring::ScoreFunctionOP get_sfxn() const { return sfxn_->clone(); }

	void
	report_and_reset_timers(
		std::chrono::duration<core::Real> &pack_time,
		std::chrono::duration<core::Real> &min_time
	) {
		pack_time = pack_time_;
		min_time = min_time_;
		pack_time_ = min_time_ = std::chrono::duration<core::Real>{};
	}

	bool has_atom_type( int atype ){ return !(raw_faatr_.find( atype ) == raw_faatr_.end()); }

private:
	/// @brief apply the convolution specified in "smoothing_"
	void
	do_convolution_and_compute_coeffs(
		ObjexxFCL::FArray3D< float > const &rawdata,
		ObjexxFCL::FArray3D< float > &smoothed_coeffs,
		core::Real smoothing,
		bool inverted = false // take min instead of max
	);

	// helper functions to move interpolated point to grid and compute an "out of bounds" penalty
	core::Real
	move_to_boundary( numeric::xyzVector< core::Real > &idxX ) const;

	core::Real
	move_to_boundary( numeric::xyzVector< core::Real > &idxX, numeric::xyzVector< core::Real > &dpen );

	/// @brief fast etable energies
	inline void
	fast_eval_etable_split_fasol(
		core::scoring::etable::Etable const &etable,
		core::conformation::Atom const &atom1,
		core::conformation::Atom const &atom2,
		core::Real & lj_atrE,
		core::Real & lj_repE,
		core::Real & fa_solE1,
		core::Real & fa_solE2
	);

	/// @brief fast lk_ball
	inline core::Real
	fast_get_lk_fractional_contribution(
		numeric::xyzVector< core::Real > const &at2,
		core::Size const atom1_n_attached_waters,
		core::scoring::lkball::WaterCoords const & atom1_waters,
		core::Real ramp_width_A2,
		core::Real d2_low,
		core::Real multiwater_fade
	);

	/// @brief fast lk_br
	inline core::Real
	fast_get_lkbr_fractional_contribution(
		core::Size const atom1_n_attached_waters,
		core::Size const atom2_n_attached_waters,
		core::scoring::lkball::WaterCoords const & atom1_waters,
		core::scoring::lkball::WaterCoords const & atom2_waters,
		core::Real overlap_gap_A2,
		core::Real overlap_width_A2,
		core::Real multiwater_fade
	);

private:
	// raw data
	core::scoring::ScoreFunctionOP sfxn_, sfxn_clash_, sfxn_soft_, sfxn_cart_, sfxn_cst_;
	core::scoring::ScoreFunctionOP sfxn_1b_, sfxn_1b_clash_, sfxn_1b_soft_;
	core::scoring::etable::Etable etable_;
	core::pose::PoseOP ref_pose_;
	bool has_cst_energies_;

	// needed for scoring
	core::scoring::lkball::LK_BallEnergyOP LKBe_;
	core::scoring::etable::coulomb::CoulombOP coulomb_;
	core::scoring::hbonds::HBondDatabaseCOP hb_database_;
	core::energy_methods::CartesianBondedEnergyOP cartbonded_;

	// alternate modes
	bool exact_, debug_, force_exact_min_;

	// reporting runtime
	std::chrono::duration<core::Real> min_time_, pack_time_;

	//unique atoms grid calculated for
	std::map< int, core::conformation::Atom > uniq_atoms_; // KEEP AS INT, IT IS NOT A SIZE


	// gridded spline coeffs
	numeric::xyzVector< core::Real > origin_;  // cartesian coordinates of index (1,1,1)
	core::Real bbox_padding_, voxel_spacing_, hash_grid_;
	core::Size hash_subgrid_;
	numeric::xyzVector< core::Size > dims_;

	// etable function values per atom type
	std::map< int, ObjexxFCL::FArray3D< float > > raw_faatr_, raw_farep_, raw_fasol_, raw_lkball_, raw_lkbridge_;
	ObjexxFCL::FArray3D< float > raw_faelec_;

	// etable coeffs per atom type
	std::map< int, ObjexxFCL::FArray3D< float > > coeffs_faatr_, coeffs_farep_, coeffs_fasol_, coeffs_lkball_, coeffs_lkbridge_;
	ObjexxFCL::FArray3D< float > coeffs_faelec_;

	// fast lookup table for polars
	GridHash3D<hbDon> hbdonors_;
	GridHash3D<hbAcc> hbacceptors_;

	core::Real maxdis_;

	// penalty factor when ligands hits boundary
	core::Real out_of_bound_e_;

	// fade function for lk ball
	core::Real LK_fade_;

	// weight on LJ repulsion b/w grid
	core::Real w_rep_;

	// original weight on elec & hbond /
	core::Real w_fa_elec_, w_hbond_sc_, w_hbond_bb_sc_;
	core::Real elec_scale_; //temporary scaling factor

	// global smoothing factor
	core::Real smoothing_;

	// # packer cycles
	core::Size packer_cycles_;

	// for minimize
	core::Size maxiter_minimize_;

	// for adjusting grid box
	numeric::xyzVector< core::Real > lig_com_;
	core::Real maxRad_;

	utility::vector1< core::Size > ligids_;
};

typedef utility::pointer::shared_ptr< GridScorer > GridScorerOP;
typedef utility::pointer::shared_ptr< GridScorer const > GridScorerCOP;


// helper function 1
inline void
GridScorer::fast_eval_etable_split_fasol(
	core::scoring::etable::Etable const &etable,
	core::conformation::Atom const &at1,
	core::conformation::Atom const &at2,
	core::Real & lj_atrE,
	core::Real & lj_repE,
	core::Real & fa_solE1,
	core::Real & fa_solE2
) {
	lj_atrE = lj_repE = fa_solE1 = fa_solE2 = 0;

	core::Real dis2 = at1.xyz().distance_squared( at2.xyz() );
	int atype1 = at1.type() < at2.type() ? at1.type() : at2.type();
	int atype2 = at1.type() < at2.type() ? at2.type() : at1.type();

	core::Real ljE = 0., atrE = 0., repE = 0.;

	if ( dis2 > etable.max_dis2() + 1e-4 ) return;

	core::scoring::etable::EtableParamsOnePair const & p = etable.analytic_params_for_pair( atype1, atype2 );
	bool inorder = at1.type() < at2.type(); // otherwise, atype1 and atype2 are swapped

	if ( dis2 > p.maxd2 ) return;
	if ( dis2 < etable.min_dis2() ) dis2 = etable.min_dis2();

	core::Real const dis = std::sqrt(dis2);
	core::Real const inv_dis2 = 1.0/dis2;

	if ( dis2 < p.ljrep_linear_ramp_d2_cutoff ) {
		ljE = dis*p.lj_switch_slope + p.lj_switch_intercept;
	} else if ( dis < p.ljatr_cubic_poly_xlo ) {
		core::Real const inv_dis6  = inv_dis2 * inv_dis2 * inv_dis2;
		ljE = ( p.lj_r12_coeff * inv_dis6 + p.lj_r6_coeff ) * inv_dis6;
	} else if ( dis < p.ljatr_cubic_poly_xhi ) {
		ljE = numeric::eval_cubic_polynomial( dis, p.ljatr_cubic_poly_parameters );
	}

	if ( dis < p.lj_minimum ) {
		atrE = p.lj_val_at_minimum;
		repE = ljE - p.lj_val_at_minimum;
	} else {
		atrE = ljE;
	}

	lj_atrE = atrE * p.ljatr_final_weight;
	lj_repE = repE;

	core::Real &fa_solE1sort = inorder? fa_solE1 : fa_solE2;
	core::Real &fa_solE2sort = inorder? fa_solE2 : fa_solE1;

	if ( dis < p.fasol_cubic_poly_close_start ) {
		fa_solE1sort = p.fasol_cubic_poly1_close_flat * p.fasol_final_weight;
		fa_solE2sort = p.fasol_cubic_poly2_close_flat * p.fasol_final_weight;
	} else if ( dis < p.fasol_cubic_poly_close_end ) {
		fa_solE1sort = p.fasol_final_weight * numeric::eval_cubic_polynomial( dis, p.fasol_cubic_poly1_close );
		fa_solE2sort = p.fasol_final_weight * numeric::eval_cubic_polynomial( dis, p.fasol_cubic_poly2_close );
	} else if ( dis < etable.fasol_cubic_poly_far_xlo() ) {
		/// exponential evaluation
		core::Real dis_rad1 = dis - etable.lj_radius(atype1);
		core::Real x1 = ( dis_rad1 * dis_rad1 ) * etable.lk_inv_lambda2(atype1);
		core::Real dis_rad2 = dis - etable.lj_radius(atype2);
		core::Real x2 = ( dis_rad2 * dis_rad2 ) * etable.lk_inv_lambda2(atype2);

		fa_solE1sort = p.fasol_final_weight * inv_dis2 * std::exp(-x1) * p.lk_coeff1;
		fa_solE2sort = p.fasol_final_weight * inv_dis2 * std::exp(-x2) * p.lk_coeff2;

	} else if ( dis < etable.fasol_cubic_poly_far_xhi() ) {
		fa_solE1sort = p.fasol_final_weight * numeric::eval_cubic_polynomial( dis, p.fasol_cubic_poly1_far );
		fa_solE2sort = p.fasol_final_weight * numeric::eval_cubic_polynomial( dis, p.fasol_cubic_poly2_far );
	}
}

// helper function 2
// for speed do non-soft variant
//    since this populates grid, derivs don't have to be continuous
inline core::Real
GridScorer::fast_get_lk_fractional_contribution(
	numeric::xyzVector< core::Real > const &atom2_xyz,
	core::Size const atom1_n_attached_waters,
	core::scoring::lkball::WaterCoords const & atom1_waters,
	core::Real ramp_width_A2,
	core::Real d2_low,
	core::Real //temp
) {
	core::Real weighted_d2_water_delta = ramp_width_A2;
	for ( core::Size idx = 1; idx <= atom1_n_attached_waters; ++idx ) {
		core::Real d2_delta = atom2_xyz.distance_squared( atom1_waters[idx] ) - d2_low;
		//weighted_d2_water_delta += exp( -d2_delta/multiwater_fade );
		weighted_d2_water_delta = std::min( weighted_d2_water_delta, d2_delta );
	}
	//weighted_d2_water_delta = -multiwater_fade * log( weighted_d2_water_delta );

	core::Real frac( 0.0 );
	if ( weighted_d2_water_delta < 0.0 ) {
		frac = 1.0;
	} else if ( weighted_d2_water_delta < ramp_width_A2 ) {
		core::Real xprime( weighted_d2_water_delta/ramp_width_A2 );
		frac = ( 1 - xprime*xprime );
		frac *= frac;
	}
	return frac;
}

// helper function 3
// for speed do non-soft variant
//    since this populates grid, derivs don't have to be continuous
inline core::Real
GridScorer::fast_get_lkbr_fractional_contribution(
	core::Size const atom1_n_attached_waters,
	core::Size const atom2_n_attached_waters,
	core::scoring::lkball::WaterCoords const & atom1_waters,
	core::scoring::lkball::WaterCoords const & atom2_waters,
	core::Real overlap_gap_A2,
	core::Real overlap_width_A2,
	core::Real //temp
) {
	core::Real weighted_d2_water_delta = overlap_width_A2;
	for ( core::Size idx1 = 1; idx1 <= atom1_n_attached_waters; ++idx1 ) {
		for ( core::Size idx2 = 1; idx2 <= atom2_n_attached_waters; ++idx2 ) {
			core::Real d2_delta = (atom1_waters[idx1] - atom2_waters[idx2]).length_squared() - overlap_gap_A2;
			weighted_d2_water_delta = std::min( weighted_d2_water_delta, d2_delta );
		}
	}


	core::Real frac( 0.0 );
	if ( weighted_d2_water_delta < 0.0 ) {
		frac = 1.0;
	} else if ( weighted_d2_water_delta < overlap_width_A2 ) {
		core::Real xprime( weighted_d2_water_delta/overlap_width_A2 );
		frac = ( 1 - xprime*xprime );
		frac *= frac;
	}

	return frac;
}

// helper function for getting hbond energies
inline
core::Real
get_hbond_score_weighted (
	core::scoring::ScoreFunctionOP sf,
	core::scoring::hbonds::HBondDatabaseCOP database,
	core::scoring::hbonds::HBondOptions const & hbopt,
	core::scoring::hbonds::HBEvalTuple const &hbt,
	numeric::xyzVector<core::Real> const &D,
	numeric::xyzVector<core::Real> const &H,
	numeric::xyzVector<core::Real> const &A,
	numeric::xyzVector<core::Real> const &B,
	numeric::xyzVector<core::Real> const &B_0
) {
	core::Real energy=0;
	core::scoring::hbonds::hb_energy_deriv(
		*database, hbopt, hbt,
		D,H,A,B,B_0,
		energy, false, core::scoring::hbonds::DUMMY_DERIVS
	);
	if ( energy<hbopt.max_hb_energy() ) {
		switch(core::scoring::hbonds::get_hbond_weight_type(hbt.eval_type())){
		case core::scoring::hbonds::hbw_SR_BB_SC :
		case core::scoring::hbonds::hbw_LR_BB_SC :
			energy *= sf->get_weight( core::scoring::hbond_bb_sc );
			break;
		case core::scoring::hbonds::hbw_SC :
		default :
			energy *= sf->get_weight( core::scoring::hbond_sc );
		}
		return energy;
	}
	return 0.0;
}


} // ga_ligand_dock
} // ligand_docking
} // protocols

#endif
