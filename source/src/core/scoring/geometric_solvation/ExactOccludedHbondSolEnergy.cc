// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk
/// @author Andrea Bazzoli (bazzoli@ku.edu)

// Project Headers
#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh>
#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergyCreator.hh>

#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Atom.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <basic/Tracer.hh>

#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/constants.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

#include <basic/options/keys/score.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>
#include <vector>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <basic/options/option.hh>

// Utility headers
#include <utility/Bound.hh>
#include <utility/string_util.hh>
#include <utility/vector0_bool.hh>
#include <utility/file/FileName.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.geometric_solvation.ExactOccludedHbondSolEnergy" );

namespace core {
namespace scoring {
namespace geometric_solvation {


/// @details This must return a fresh instance of the ExactOccludedHbondSolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ExactOccludedHbondSolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {

	return create_ExactSHOEnergy_from_cmdline(options);
}


ScoreTypes
ExactOccludedHbondSolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( occ_sol_exact );
	return sts;
}


using namespace core;
using namespace core::scoring;
using namespace core::scoring::hbonds;

core::Real const geosol_kT = { 0.593 };

// slope and intercept of the linear function to scale SHO values to the EEF1 range. They were obtained by
// fitting average SHO energies to average LK energies over an ensemble of non-H-bonded heavy atoms from
// 207 high-resolutionstructures.
core::Real const LK_MATCHING_SLOPE = { 0.4775 };
core::Real const LK_MATCHING_INTERCEPT = { 0 }; // to be splitted among polar groups forming the same LK heavy atom

// private constructor
GridInfo::GridInfo() {

	// Setup dimensions for water grids - use the same dimensions for all grids
	core::Real const water_grid_width = 10.;
	core::Real const water_grid_depth = 8.;

	// Use55 thousand points per grid
	xnum_points_ = 41;
	ynum_points_ = 41;
	znum_points_ = 33;

	xstep_ = water_grid_width / ( xnum_points_ - 1 );
	ystep_ = water_grid_width / ( ynum_points_ - 1 );
	zstep_ = water_grid_depth / ( znum_points_ - 1 );

	// Note: the point at the origin will NOT be considered in calculations - the grid starts AFTER the origin!!
	xorigin_ = -water_grid_width/2 - xstep_;
	yorigin_ = -water_grid_width/2 - ystep_;
	zorigin_ = 0.;
}

// private constructor
WaterWeightGridSet::WaterWeightGridSet() :
	hbondoptions_( hbonds::HBondOptionsOP( new HBondOptions ) ),
	hb_database_(HBondDatabase::get_database())
{
	using namespace hbonds;

	// override command line settings
	hbondoptions_->use_sp2_chi_penalty(false);
	hbondoptions_->measure_sp3acc_BAH_from_hvy(false);

	// We need water grids for each donor and acceptor type
	// We could read them in, but let's just compute them from scratch instead
	// We'll store them as a map keyed on hbonds::HBEvalType, with each value a 3D vector of the weights

	TR << "computing and storing water weight grids for acceptor types." << std::endl;
	for ( Size i = 1; i <= hbacc_MAX; i++ ) {
		HBEvalTuple const hbe( hbdon_H2O, HBAccChemType(i), seq_sep_other );
		sum_all_water_weights_[hbe.eval_type()] = fill_water_grid( all_water_weights_[hbe.eval_type()], hbe, *GridInfo::get_instance(), true /*water is donor*/);
	}

	TR << "computing and storing water weight grids for donor types." << std::endl;
	for ( Size i = 1; i <= hbdon_MAX; i++ ) {
		HBEvalTuple const hbe( HBDonChemType(i), hbacc_H2O, seq_sep_other );
		sum_all_water_weights_[hbe.eval_type()] = fill_water_grid( all_water_weights_[hbe.eval_type()], hbe, *GridInfo::get_instance(), false /*water is acceptor*/);
	}
}


// Fill in the water grid
core::Real WaterWeightGridSet::fill_water_grid(
	WaterWeightGridSet::Grid & water_weights,
	hbonds::HBEvalTuple const & hbond_eval_type,
	GridInfo const & grid_info, bool const water_is_donor)
{

	static Vector const base_to_outer(0,0,1);
	core::Real const entropy_scaling = 1.0;
	core::Real const water_O_H_distance = 0.958;

	// Setup grid, initialize to zero
	water_weights.clear();
	water_weights.resize(grid_info.xnum_points());
	for ( core::Size tx=0; tx<grid_info.xnum_points(); tx++ ) {
		water_weights[tx].resize(grid_info.ynum_points());
		for ( core::Size ty=0; ty<grid_info.ynum_points(); ty++ ) {
			water_weights[tx][ty].resize(grid_info.znum_points(), 0.);
		}
	}

	// Fill in the water weight grid
	core::Real sum_grid_water_weight = 0.;
	core::Vector water_position(grid_info.xorigin(),grid_info.yorigin(),grid_info.zorigin());

	for ( core::Size tx=0; tx<grid_info.xnum_points(); tx++ ) {
		water_position.x() += grid_info.xstep();
		water_position.y() = grid_info.yorigin();
		for ( core::Size ty=0; ty<grid_info.ynum_points(); ty++ ) {
			water_position.y() += grid_info.ystep();
			water_position.z() = grid_info.zorigin();
			for ( core::Size tz=0; tz<grid_info.znum_points(); tz++ ) {
				water_position.z() += grid_info.zstep();

				// Compute the current geometry
				core::Real AHdis, xD, xH;
				if ( water_is_donor ) {

					// water is the donor, give it perfect geometry
					xD = 0.9999;

					// compute the distance to the accepting water proton
					// subtract the water's OH distance to get the AHdis,
					// since the distance computed was from the acceptor to the water oxygen
					// note: water proton lies on the line between the acceptor and the water oxygen
					AHdis = water_position.length();
					AHdis -= water_O_H_distance; // water O-H distance

					// find cosine of the base-acceptor-water_proton angle (xH)
					// note: this is the same as the base-acceptor-water_oxygen angle
					// note: careful to normalize-by-copy, since water_position is reused around the loop
					xH = dot( base_to_outer, water_position.normalized() );

				} else {

					// water is the acceptor, give it perfect geometry (lone pair along the donor-acceptor line)
					xH = 1./3.;  // perfect geometry is cos( 180 - 109.5 degrees), which is 1/3

					// compute the distance to the accepting water
					AHdis = water_position.length();

					// find the cosine of the base-proton-water angle (xD)
					// note: careful to normalize-by-copy, since water_position is reused around the loop
					xD = dot( base_to_outer, water_position.normalized() );

				}

				if ( xH < MIN_xH ) continue;
				if ( xH > MAX_xH ) continue;
				if ( xD < MIN_xD ) continue;
				if ( xD > MAX_xD ) continue;
				if ( AHdis < MIN_R ) continue;
				if ( AHdis > MAX_R ) continue;

				// Get the Hbond energy
				core::Real curr_water_hbond;
				core::Real dummy_chi( 0.0 );
				debug_assert( ! hbondoptions_->use_sp2_chi_penalty() ); // APL avoid the new sp2 chi term.
				debug_assert( ! hbondoptions_->measure_sp3acc_BAH_from_hvy() );
				hbond_compute_energy(*hb_database_, *hbondoptions_, hbond_eval_type,
					AHdis, xD, xH, dummy_chi, curr_water_hbond );

				// Save the Hbond energy
				curr_water_hbond *= entropy_scaling;
				if ( curr_water_hbond < 0 ) {
					core::Real curr_water_weight = exp( - curr_water_hbond / geosol_kT );
					water_weights[tx][ty][tz] = curr_water_weight;
					sum_grid_water_weight += curr_water_weight;
				}
			}
		}
	}

	return sum_grid_water_weight;
}


WaterWeightGridSet::Grid const &
WaterWeightGridSet::get_water_weight_grid( hbonds::HBEvalType const & hbond_eval_type ) const {

	// Check that we have weights for this Hbond type
	all_water_weights_iterator curr_water_weights_iter = all_water_weights_.find( hbond_eval_type );
	if ( curr_water_weights_iter == all_water_weights_.end( ) ) {
		TR << "Could not look up map element" << std::endl;
		TR << "get_water_weight_grid hbond_eval_type " << hbond_eval_type << std::endl;
		debug_assert(false);
		exit(1);
	}
	return curr_water_weights_iter->second;
}


core::Real
WaterWeightGridSet::get_sum_water_weight_grid( hbonds::HBEvalType const & hbond_eval_type ) const {

	// Check that we have weights for this Hbond type
	sum_water_weights_iterator curr_sum_water_weights_iter = sum_all_water_weights_.find( hbond_eval_type );
	if ( curr_sum_water_weights_iter == sum_all_water_weights_.end( ) ) {
		TR << "Could not look up map element" << std::endl;
		TR << "get_sum_water_weight_grid hbond_eval_type " << hbond_eval_type << std::endl;
		debug_assert(false);
		exit(1);
	}
	return curr_sum_water_weights_iter->second;
}


///
/// @brief prints a given xz-plane of a water grid
///
/// @param[in] hbond_eval_type HBEvalType of the interaction between the grid's polar group and a water molecule
/// @param[in] y y-coordinate of the xz-plane
///
/// @details grid values are printed, from top to bottom, by decreasing z-coordinate, and, from left to right, by\n
///  increasing x-coordinate. Namely, the ith output line contains the values for z = N-1-i (i=0,...,N-1,\n
///  where N is the number of values assumed by the z-coordinate); within each line, the jth column contains the\n
///  value for x = j (j=0,...,M-1, where M is the number of values assumed by the x-coordinate).
///
void WaterWeightGridSet::print_water_weight_grid_xz_plane(
	hbonds::HBEvalType const & hbond_eval_type,
	int const y) const {

	Grid const & grid = get_water_weight_grid(hbond_eval_type);

	TR << "printing xz-plane at y = " << y << " for water grid of HBEvalType " << hbond_eval_type
		<< ". Top-to-bottom: decreasing z; left-to-right: increasing x" << std::endl;
	for ( int z=GridInfo::get_instance()->znum_points()-1; z >= 0 ; z-- ) {
		for ( int x=0; x<(int)(GridInfo::get_instance()->xnum_points()); x++ ) {
			TR.width(12);
			TR << grid[x][y][z];
		}
		TR << std::endl;
	}
}


ExactOccludedHbondSolEnergy::~ExactOccludedHbondSolEnergy() {}


void ExactOccludedHbondSolEnergy::allocate_grid_of_occluded_sites() {
	occluded_sites_.clear();
	occluded_sites_.resize(GridInfo::get_instance()->xnum_points());
	for ( core::Size tx=0; tx<GridInfo::get_instance()->xnum_points(); tx++ ) {
		occluded_sites_[tx].resize(GridInfo::get_instance()->ynum_points());
		for ( core::Size ty=0; ty<GridInfo::get_instance()->ynum_points(); ty++ ) {
			occluded_sites_[tx][ty].resize(GridInfo::get_instance()->znum_points());
		}
	}
}


ExactOccludedHbondSolEnergy::ExactOccludedHbondSolEnergy(
	etable::Etable const & etable_in,
	bool const analytic_etable_evaluation,
	bool const exact_occ_skip_Hbonders,
	bool const exact_occ_pairwise,
	bool const exact_occ_pairwise_by_res,
	bool const exact_occ_split_between_res,
	bool const exact_occ_self_res_occ,
	core::Real const occ_radius_scaling,
	bool const verbose
) :
	parent( methods::EnergyMethodCreatorOP( new ExactOccludedHbondSolEnergyCreator ) ),
	exact_occ_skip_Hbonders_( exact_occ_skip_Hbonders ),
	exact_occ_pairwise_( exact_occ_pairwise ),
	exact_occ_pairwise_by_res_( exact_occ_pairwise_by_res ),
	exact_occ_split_between_res_( exact_occ_split_between_res ),
	exact_occ_self_res_occ_( exact_occ_self_res_occ ),
	occ_radius_scaling_( occ_radius_scaling ),
	hbondoptions_( hbonds::HBondOptionsOP( new HBondOptions ) ),
	hb_database_( HBondDatabase::get_database() ),
	hbond_set_( hbonds::HBondSetOP( new hbonds::HBondSet( *hbondoptions_ ) ) ),
	lk_safe_max_dis2_( etable_in.get_safe_max_dis2() ),
	verbose_( verbose )
{
	hbondoptions_->use_sp2_chi_penalty( false ); // apl preserve old behavior
	hbondoptions_->measure_sp3acc_BAH_from_hvy( false );

	if ( verbose_ ) TR <<"ExactOccludedHbondSolEnergy constructor" << std::endl;
	if ( exact_occ_split_between_res_ && ! exact_occ_pairwise_ ) {
		TR << "Error - cannot split occ energy between residues unless pairwise calculations are used!" << std::endl;
		exit(1);
	}

	// Keep a copy of the AtomTypeSet (to lookup atomic radii)
	atom_type_set_ptr_ = chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD );

	// force allocation of water grids
	WaterWeightGridSet::get_instance();

	allocate_grid_of_occluded_sites();

	if ( analytic_etable_evaluation ) {
		etable_evaluator_ = etable::EtableEvaluatorOP( new etable::AnalyticEtableEvaluator( etable_in ) );
	} else {
		etable_evaluator_ = etable::EtableEvaluatorOP( new etable::TableLookupEvaluator( etable_in ) );
	}
}


ExactOccludedHbondSolEnergy::ExactOccludedHbondSolEnergy( ExactOccludedHbondSolEnergy const & src ):
	parent( src ),
	exact_occ_skip_Hbonders_( src.exact_occ_skip_Hbonders_ ),
	exact_occ_pairwise_( src.exact_occ_pairwise_ ),
	exact_occ_pairwise_by_res_( src.exact_occ_pairwise_by_res_ ),
	exact_occ_split_between_res_( src.exact_occ_split_between_res_ ),
	exact_occ_self_res_occ_( src.exact_occ_self_res_occ_ ),
	occ_radius_scaling_( src.occ_radius_scaling_ ),
	hbondoptions_(src.hbondoptions_),
	hb_database_(src.hb_database_),
	atom_type_set_ptr_( src.atom_type_set_ptr_ ),
	hbond_set_(hbonds::HBondSetOP(new hbonds::HBondSet(*src.hbond_set_))),
	lk_safe_max_dis2_( src.lk_safe_max_dis2_ ),
	verbose_( src.verbose_ )
{
	allocate_grid_of_occluded_sites();
}


methods::EnergyMethodOP
ExactOccludedHbondSolEnergy::clone() const
{
	return methods::EnergyMethodOP( new ExactOccludedHbondSolEnergy( *this ) );
}


void
ExactOccludedHbondSolEnergy::init_hbond_data( pose::Pose const& pose) const
{
	hbond_set_->resize_bb_donor_acceptor_arrays( pose.size() );
	core::scoring::hbonds::fill_hbond_set(pose, false, *hbond_set_);
}


void
ExactOccludedHbondSolEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const &) const
{
	pose.update_residue_neighbors();
	init_hbond_data(pose);
}


void
ExactOccludedHbondSolEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	pose.update_residue_neighbors();
}


void
ExactOccludedHbondSolEnergy::setup_for_derivatives( pose::Pose & , ScoreFunction const & ) const
{
	TR << "Error - cannot compute derivatives for ExactOccludedHbondSolEnergy (occ_sol_exact)" << std::endl;
	debug_assert(false);
	exit(1);
}


void
ExactOccludedHbondSolEnergy::setup_for_minimizing( pose::Pose & , ScoreFunction const & , kinematics::MinimizerMapBase const & ) const
{
	TR << "Error - cannot compute derivatives for ExactOccludedHbondSolEnergy (occ_sol_exact)" << std::endl;
	debug_assert(false);
	exit(1);
}


void ExactOccludedHbondSolEnergy::residue_energy(
	conformation::Residue const & polar_rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const {

	core::Size polar_resnum = (core::Size) polar_rsd.seqpos();
	core::Real residue_geosol(0.);

	// loop over donors in polar_rsd
	for ( Size const polar_atom : polar_rsd.Hpos_polar() ) {
		residue_geosol += compute_donor_atom_energy(polar_rsd, polar_resnum, polar_atom, pose);
	}

	// loop over acceptors in polar_rsd
	for ( Size const polar_atom : polar_rsd.accpt_pos() ) {
		residue_geosol += compute_acceptor_atom_energy(polar_rsd, polar_resnum, polar_atom, pose);
	}

	if ( exact_occ_split_between_res_ ) {

		// Here we need to add contributions from this residue occluding all neighboring polar groups then divide everything by two
		// note: this will make the code run twice as slow as it otherwise would, but that's okay because it's just for parameterization / analysis
		utility_exit_with_message( "PAIRWISE OUTPUT FORMAT IS NOT YET SUPPORTED" );
	}

	emap[ occ_sol_exact ] += residue_geosol;
}


/// @brief computes energy of fully buried polar group
///
/// @details The energy is +5 kcal/mole in any case
///
core::Real ExactOccludedHbondSolEnergy::compute_fully_buried_ene() const {

	return 5;
}


///
/// @brief computes the desolvation energy (i.e., either SHO or LK) of a donor atom
///
/// @param[in] polar_rsd residue to which the atom belongs
/// @param[in] polar_resnum index of the residue in its pose
/// @param[in] polar_atom index of the atom in the residue
/// @param[in] pose the pose
///
core::Real ExactOccludedHbondSolEnergy::compute_donor_atom_energy(
	conformation::Residue const & polar_rsd,
	core::Size polar_resnum,
	core::Size const polar_atom,
	pose::Pose const & pose) const {

	core::id::AtomID aid(polar_atom, polar_resnum);
	Size const base = polar_rsd.atom_base(polar_atom);

	// split factor is at least the number of hydrogens attached to the base
	Size split_factor = polar_rsd.attached_H_end(base) - polar_rsd.attached_H_begin(base) + 1;
	if ( polar_rsd.heavyatom_is_an_acceptor( base ) ) {
		// base is also an acceptor, so it contributes to the splitting
		split_factor++;
	}

	if ( hbond_set_->atom_hbonds(aid).size() > 0 ) {
		// atom is H-bonded: use LK energy of base atom divided by split factor
		core::Real lk_ene = get_atom_lk_energy(base, polar_rsd, pose);
		return lk_ene / split_factor;
	} else {
		// atom is not H-bonded: use SHO energy
		core::Real sho_ene = compute_sho_donor_atom_energy(polar_rsd, polar_resnum, polar_atom, pose);
		return ( LK_MATCHING_SLOPE * sho_ene ) + ( LK_MATCHING_INTERCEPT / split_factor );
	}
}


/// @brief computes the SHO energy of a donor atom
///
/// @param[in] polar_rsd residue to which the atom belongs
/// @param[in] polar_resnum index of the residue in its pose
/// @param[in] polar_atom index of the atom in the residue
/// @param[in] pose the pose
///
core::Real ExactOccludedHbondSolEnergy::compute_sho_donor_atom_energy(
	conformation::Residue const & polar_rsd,
	core::Size polar_resnum,
	core::Size const polar_atom,
	pose::Pose const & pose) const {

	Size const base_atom( polar_rsd.atom_base( polar_atom ) );
	HBEvalTuple const curr_hbond_eval_type( get_hb_don_chem_type( base_atom, polar_rsd ), hbacc_H2O, seq_sep_other );

	core::Real const fully_buried_ene = compute_fully_buried_ene();
	core::Real const grid_constant = compute_grid_constant(curr_hbond_eval_type, fully_buried_ene);

	core::Real polar_group_energy = 0.;
	if ( ! exact_occ_pairwise_ && ! exact_occ_pairwise_by_res_ ) {
		polar_group_energy = compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
			grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type.eval_type() ) );
	} else {
		// loop over all atoms of neighboring residues, INCLUDING SELF
		core::scoring::TenANeighborGraph const & graph = pose.energies().tenA_neighbor_graph();
		utility::vector1 <core::Size> neighborlist;
		neighborlist.push_back( polar_resnum);
		for ( utility::graph::Graph::EdgeListConstIter
				neighbor_iter = graph.get_node( polar_resnum )->const_edge_list_begin(),
				neighbor_iter_end = graph.get_node( polar_resnum )->const_edge_list_end();
				neighbor_iter != neighbor_iter_end; ++neighbor_iter ) {
			neighborlist.push_back( (*neighbor_iter)->get_other_ind( polar_resnum ) );
		}
		for ( Size const occ_resnum : neighborlist ) {
			conformation::Residue const & occ_rsd = pose.residue(occ_resnum);
			if ( exact_occ_pairwise_by_res_ ) {
				polar_group_energy += compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
					grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type.eval_type() ),
					true, occ_resnum );
			} else {
				for ( Size occ_atomno = 1; occ_atomno <= occ_rsd.natoms(); ++occ_atomno ) {
					polar_group_energy += compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
						grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type.eval_type() ),
						true, occ_resnum, true, occ_atomno );
				}
			}
		}
	}

	return  polar_group_energy;
}


/// @brief computes the desolvation energy (i.e., either SHO or LK) of an acceptor atom
///
/// @param[in] polar_rsd residue to which the atom belongs
/// @param[in] polar_resnum index of the residue in its pose
/// @param[in] polar_atom index of the atom in the residue
/// @param[in] pose the pose
///
core::Real ExactOccludedHbondSolEnergy::compute_acceptor_atom_energy(
	conformation::Residue const & polar_rsd,
	core::Size polar_resnum,
	core::Size const polar_atom,
	pose::Pose const & pose) const {

	core::id::AtomID aid(polar_atom, polar_resnum);
	Size split_factor = 1;
	if ( polar_rsd.heavyatom_has_polar_hydrogens( polar_atom ) ) {
		// atom is also a donor: divide by 1 + number of attached hydrogens
		Size const NH = polar_rsd.attached_H_end(polar_atom) - polar_rsd.attached_H_begin(polar_atom) + 1;
		split_factor += NH;
	}

	if ( hbond_set_->atom_hbonds(aid).size() > 0 ) {
		// atom is H-bonded: use LK energy
		core::Real lk_ene = get_atom_lk_energy(polar_atom, polar_rsd, pose);
		return lk_ene / split_factor;
	} else {
		// atom is not H-bonded: use SHO energy
		core::Real sho_ene = compute_sho_acceptor_atom_energy( polar_rsd, polar_resnum, polar_atom, pose );
		return ( LK_MATCHING_SLOPE * sho_ene ) + ( LK_MATCHING_INTERCEPT / split_factor );
	}
}


///
/// @brief computes the SHO energy of an acceptor atom
///
/// @param[in] polar_rsd residue to which the atom belongs
/// @param[in] polar_resnum index of the residue in its pose
/// @param[in] polar_atom index of the atom in the residue
/// @param[in] pose the pose
///
core::Real ExactOccludedHbondSolEnergy::compute_sho_acceptor_atom_energy(
	conformation::Residue const & polar_rsd,
	core::Size polar_resnum,
	core::Size const polar_atom,
	pose::Pose const & pose) const {

	chemical::AtomTypeSetCOP atom_type_set_ptr( atom_type_set_ptr_ );
	hbonds::HBEvalTuple const curr_hbond_eval_type( hbdon_H2O, get_hb_acc_chem_type( polar_atom, polar_rsd ), seq_sep_other);

	core::Real const fully_buried_ene = compute_fully_buried_ene();
	core::Real const grid_constant = compute_grid_constant(curr_hbond_eval_type, fully_buried_ene);

	core::Real polar_group_energy = 0.;
	if ( ! exact_occ_pairwise_ && ! exact_occ_pairwise_by_res_ ) {
		polar_group_energy = compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
			grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type.eval_type() ) );
	} else {
		// loop over all atoms of neighboring residues, INCLUDING SELF
		core::scoring::TenANeighborGraph const & graph = pose.energies().tenA_neighbor_graph();
		utility::vector1 <core::Size> neighborlist;
		neighborlist.push_back( polar_resnum);
		for ( utility::graph::Graph::EdgeListConstIter
				neighbor_iter = graph.get_node( polar_resnum )->const_edge_list_begin(),
				neighbor_iter_end = graph.get_node( polar_resnum )->const_edge_list_end();
				neighbor_iter != neighbor_iter_end; ++neighbor_iter ) {
			neighborlist.push_back( (*neighbor_iter)->get_other_ind( polar_resnum ) );
		}
		for ( Size const occ_resnum : neighborlist ) {
			conformation::Residue const & occ_rsd = pose.residue(occ_resnum);
			if ( exact_occ_pairwise_by_res_ ) {
				polar_group_energy += compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
					grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type.eval_type() ),
					true, occ_resnum );
			} else {
				for ( Size occ_atomno = 1; occ_atomno <= occ_rsd.natoms(); ++occ_atomno ) {
					polar_group_energy += compute_polar_group_sol_energy(pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
						grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type.eval_type() ),
						true, occ_resnum, true, occ_atomno );
				}
			}
		}
	}

	return polar_group_energy;
}


Real ExactOccludedHbondSolEnergy::compute_polar_group_sol_energy(
	pose::Pose const & pose,
	conformation::Residue const & polar_rsd,
	Size const polar_atom,
	bool const restrict_to_single_occluding_residue, // = false
	Size const single_occluding_resinx, // = 0
	bool const restrict_to_single_occluding_atom, // = false
	Size const single_occluding_atominx // = 0
) const {
	Size const base_atom( polar_rsd.atom_base( polar_atom ) );

	HBEvalTuple curr_hbond_eval_type;
	bool const atom_is_donor = polar_rsd.atom_type( base_atom ).is_donor();
	if ( atom_is_donor ) {
		curr_hbond_eval_type = HBEvalTuple(get_hb_don_chem_type( base_atom, polar_rsd ), hbacc_H2O, seq_sep_other);
	} else if ( polar_rsd.atom_type( polar_atom).is_acceptor() ) {
		curr_hbond_eval_type = HBEvalTuple( hbdon_H2O, get_hb_acc_chem_type( polar_atom, polar_rsd ), seq_sep_other);
	} else {
		debug_assert( false ); // Not a donor or an acceptor, don't know what to do.
	}

	Real const fully_buried_ene = compute_fully_buried_ene();
	Real const grid_constant = compute_grid_constant( curr_hbond_eval_type, fully_buried_ene);

	return compute_polar_group_sol_energy(
		pose, polar_rsd, polar_atom, *GridInfo::get_instance(),
		grid_constant, WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_type.eval_type() ),
		restrict_to_single_occluding_residue,
		single_occluding_resinx,
		restrict_to_single_occluding_atom,
		single_occluding_atominx);
}


core::Real ExactOccludedHbondSolEnergy::compute_polar_group_sol_energy(
	pose::Pose const & pose,
	conformation::Residue const & polar_rsd,
	core::Size const polar_atomno,
	GridInfo const & grid_info,
	core::Real const & grid_constant,
	WaterWeightGridSet::Grid const & water_weights,
	bool const restrict_to_single_occluding_residue, // = false
	core::Size const single_occluding_resinx, // = 0
	bool const restrict_to_single_occluding_atom, // = false
	core::Size const single_occluding_atominx // = 0
) const {

	// note: the "restrict_to_single_occluding_residue / atom" options are so that pairwise additivity can be enforced (for parameterization / analysis)
	core::Size const polar_resnum = (core::Size) polar_rsd.seqpos();

	// jumpout immediately if we're not including self but the single requested residue is self
	if ( restrict_to_single_occluding_residue && ! exact_occ_self_res_occ_ && ( polar_resnum == single_occluding_resinx ) ) {
		return 0;
	}

	Size const base_atom( polar_rsd.atom_base( polar_atomno ) );

	core::Real polar_group_hb_energy(0.);

	bool const hydrogens_can_occlude = false;
	core::Real const water_radius = 1.4;

	// Reset grid of occluded sites, by setting everything to false (for "not occluded")
	for ( core::Size tx=0; tx<grid_info.xnum_points(); tx++ ) {
		for ( core::Size ty=0; ty<grid_info.ynum_points(); ty++ ) {
			for ( core::Size tz=0; tz<grid_info.znum_points(); tz++ ) {
				occluded_sites_[tx][ty][tz] = false;
			}
		}
	}

	// Find the transformation which puts the donor/acceptor of interest onto the existing grid
	// The plan is to apply this transformation to bring each occluding atom onto the existing grid (translation then matrix multiplication)
	core::Vector const & orig_polar_atom_xyz( polar_rsd.atom( polar_atomno ).xyz() );
	core::Vector const & orig_base_atom_xyz( polar_rsd.atom( base_atom ).xyz() );
	core::Vector const translation_vector = -orig_polar_atom_xyz;
	core::Vector const translated_base_atom = orig_base_atom_xyz + translation_vector;

	// We want to transform the positions of occluding atoms from the canonical cartesian reference frame to
	// a cartesian reference frame centered on the polar group
	core::Vector const cartesian_x(1,0,0);
	core::Vector const cartesian_y(0,1,0);
	core::Vector const cartesian_z(0,0,1);

	// There's not a unique solution for this, so we'll arbitrarily pick a second basis vector, requiring that the dot product with desired_z is zero
	core::Vector desired_z = -translated_base_atom.normalized();

	// Treat the case where polar_atom.z == base_atom.z ; this leads to desired_z.z == 0
	core::Real arbitrary_x, arbitrary_y, arbitrary_z;

	// jk note: we get numerical problems if z-coor is about zero (we divide by it below) - if so then use y to divide by
	if ( std::abs( desired_z.z() ) > 0.01 ) {
		arbitrary_x = 1;
		arbitrary_y = 1;
		arbitrary_z = -1. * ((arbitrary_x * desired_z.x()) + (arbitrary_y * desired_z.y())) / desired_z.z();
	} else {
		if ( std::abs( desired_z.y() ) > 0.01 ) {
			arbitrary_x = 1;
			arbitrary_z = 1;
			arbitrary_y = -1. * ((arbitrary_x * desired_z.x()) + (arbitrary_z * desired_z.z())) / desired_z.y();
		} else {
			// jk note: we get numerical problems if z-coor AND y-coor are about zero (it does happen!) - if so then use x to divide by
			arbitrary_y = 1;
			arbitrary_z = 1;
			arbitrary_x = -1. * ((arbitrary_y * desired_z.y()) + (arbitrary_z * desired_z.z())) / desired_z.x();
		}
	}
	core::Vector desired_x(arbitrary_x,arbitrary_y,arbitrary_z);
	desired_x.normalize();
	core::Vector desired_y = cross_product( desired_z, desired_x );

	// The transformation matrix to do this is i.i'  i.j', etc. where i,j,k are the unit vectors of the starting system
	// and i',j',k' are the unit vectors of the target system
	// for reference see http://kwon3d.com/theory/transform/transform.html
	numeric::xyzMatrix< Length > transformation_matrix;
	transformation_matrix.xx( desired_x.dot(cartesian_x) );
	transformation_matrix.xy( desired_x.dot(cartesian_y) );
	transformation_matrix.xz( desired_x.dot(cartesian_z) );
	transformation_matrix.yx( desired_y.dot(cartesian_x) );
	transformation_matrix.yy( desired_y.dot(cartesian_y) );
	transformation_matrix.yz( desired_y.dot(cartesian_z) );
	transformation_matrix.zx( desired_z.dot(cartesian_x) );
	transformation_matrix.zy( desired_z.dot(cartesian_y) );
	transformation_matrix.zz( desired_z.dot(cartesian_z) );

	// Double-check transformation matrix
	core::Vector new_base_atom_location = transformation_matrix * translated_base_atom;
	debug_assert( std::abs(new_base_atom_location.normalized().x()) < 0.001 );
	debug_assert( std::abs(new_base_atom_location.normalized().y()) < 0.001 );
	debug_assert( std::abs(new_base_atom_location.normalized().z() + 1.) < 0.001 );

	utility::vector1 <core::Size> neighborlist;
	if ( restrict_to_single_occluding_residue ) {

		// consider only a single occluding residue (for pairwise additivity)
		debug_assert ( single_occluding_resinx > 0 );
		neighborlist.push_back( single_occluding_resinx );
	} else {

		// loop over all atoms of neighboring residues, INCLUDING SELF
		core::scoring::TenANeighborGraph const & graph = pose.energies().tenA_neighbor_graph();
		neighborlist.push_back( polar_resnum);
		for ( utility::graph::Graph::EdgeListConstIter
				neighbor_iter = graph.get_node( polar_resnum )->const_edge_list_begin(),
				neighbor_iter_end = graph.get_node( polar_resnum )->const_edge_list_end();
				neighbor_iter != neighbor_iter_end; ++neighbor_iter ) {
			neighborlist.push_back( (*neighbor_iter)->get_other_ind( polar_resnum ) );
		}
	}

	chemical::AtomTypeSetCOP atom_type_set_ptr( atom_type_set_ptr_ );

	for ( Size const occ_resnum : neighborlist ) {
		if ( ! exact_occ_self_res_occ_ && ( polar_resnum == occ_resnum ) ) {
			continue;
		}
		conformation::Residue const& occ_rsd = pose.residue(occ_resnum);

		// loop over all occluding atoms in this residue
		Size atom_startinx = 1;
		Size atom_lastinx = occ_rsd.natoms();
		if ( restrict_to_single_occluding_atom ) {

			// consider only a single occluding atom (for pairwise additivity)
			debug_assert ( single_occluding_atominx > 0 );
			atom_startinx = single_occluding_atominx;
			atom_lastinx = single_occluding_atominx;
		}

		//  TR << "computing occlusion of polar residue " << polar_resnum << " by residue " << occ_resnum << std::endl;

		core::id::AtomID const polat_id(polar_atomno, polar_resnum);

		for ( Size occ_atomno = atom_startinx; occ_atomno <= atom_lastinx; ++occ_atomno ) {

			bool const occ_atom_is_hydrogen = occ_rsd.atom_is_hydrogen( occ_atomno );
			if ( occ_atom_is_hydrogen && ! hydrogens_can_occlude ) continue;

			// can be occluded by atoms directly bonded to this group, but not by self
			if ( polar_resnum == occ_resnum ) {
				if ( polar_atomno == occ_atomno ) continue;
				if ( base_atom == occ_atomno ) continue;
			}

			core::Real occ_radius = (*atom_type_set_ptr)[occ_rsd.atom_type_index( occ_atomno )].lj_radius();

			// catch proline NV here (and other virtual atoms, etc.)
			if ( occ_radius < 0.1 ) continue;
			occ_radius *= occ_radius_scaling_;
			core::Real const sq_dist_cut = ( occ_radius + water_radius ) * ( occ_radius + water_radius );

			// Apply the transformation to put this atom onto the current grid
			core::Vector const & orig_occ_atom_xyz( occ_rsd.atom( occ_atomno ).xyz() );
			core::Vector const transformed_occ_atom_xyz = transformation_matrix * ( orig_occ_atom_xyz + translation_vector );

			// Double-check transformations
			debug_assert( std::abs(orig_polar_atom_xyz.distance( orig_occ_atom_xyz ) - transformed_occ_atom_xyz.magnitude()) < 0.001 );
			debug_assert( std::abs(orig_base_atom_xyz.distance( orig_occ_atom_xyz ) - new_base_atom_location.distance( transformed_occ_atom_xyz )) < 0.001 );

			// Loop over all water positions, mark those which are occluded by this atom
			core::Vector water_position(grid_info.xorigin(),grid_info.yorigin(),grid_info.zorigin());
			for ( core::Size wx=0; wx<grid_info.xnum_points(); wx++ ) {
				water_position.x() += grid_info.xstep();
				core::Real sq_xdist = ( water_position.x() - transformed_occ_atom_xyz.x() ) * ( water_position.x() - transformed_occ_atom_xyz.x() );
				if ( sq_xdist > sq_dist_cut ) continue;
				water_position.y() = grid_info.yorigin();
				for ( core::Size wy=0; wy<grid_info.ynum_points(); wy++ ) {
					water_position.y() += grid_info.ystep();
					core::Real sq_ydist = ( water_position.y() - transformed_occ_atom_xyz.y() ) * ( water_position.y() - transformed_occ_atom_xyz.y() );
					if ( sq_ydist > sq_dist_cut ) continue;
					water_position.z() = grid_info.zorigin();
					for ( core::Size wz=0; wz<grid_info.znum_points(); wz++ ) {
						water_position.z() += grid_info.zstep();
						core::Real sq_zdist = ( water_position.z() - transformed_occ_atom_xyz.z() ) * ( water_position.z() - transformed_occ_atom_xyz.z() );
						if ( sq_zdist > sq_dist_cut ) continue;
						core::Real sq_curr_dist = sq_xdist + sq_ydist + sq_zdist;
						if ( sq_curr_dist < sq_dist_cut ) {

							// this atom occludes this water site
							occluded_sites_[wx][wy][wz] = true;
						}
					}
				}
			}

		}
	}

	// Compute and store the solvation energy for the grid occluded by all nearby atoms
	// Compute the numerator (the sum of occluded weights)
	core::Real sum_occluded_weights(0.);
	for ( core::Size tx=0; tx<grid_info.xnum_points(); tx++ ) {
		for ( core::Size ty=0; ty<grid_info.ynum_points(); ty++ ) {
			for ( core::Size tz=0; tz<grid_info.znum_points(); tz++ ) {
				if ( occluded_sites_[tx][ty][tz] ) {
					sum_occluded_weights += water_weights[tx][ty][tz];
				}
			}
		}
	}

	// Compute the energetic cost of occluding this polar group
	core::Real geometric_solvation_energy = - geosol_kT * log( 1 - ( sum_occluded_weights / grid_constant ) );
	core::Real desired_hb_weight = 0.;

	geometric_solvation_energy = geometric_solvation_energy + (desired_hb_weight*polar_group_hb_energy );
	return geometric_solvation_energy;
}


/// @brief computes the grid constant for a given polar group (i.e., the denominator in the solvation energy equation)
///
/// @param[in] hbond_eval_type HBEvalTuple describing the hydrogen bond of the polar group to water
/// @param[in] fully_buried_ene energy of the polar group when fully buried
///
core::Real ExactOccludedHbondSolEnergy::compute_grid_constant(
	core::scoring::hbonds::HBEvalTuple const & hbond_eval_type,
	core::Real fully_buried_ene ) const
{
	// Compute energy of water molecule in bulk water
	core::Real const Emax_weight = exp( -fully_buried_ene / geosol_kT );
	core::Real const sum_water_weights = WaterWeightGridSet::get_instance()->get_sum_water_weight_grid( hbond_eval_type.eval_type() );
	core::Real const Ebulk_weight = ( sum_water_weights * Emax_weight ) / ( 1. - Emax_weight);

	return sum_water_weights + Ebulk_weight;
}


/// @brief returns the LK energy of a given atom due to all its neighboring residues
///
/// @param[in] atom_idx index of the atom in its residue
/// @param[in] res the atom's residue
/// @param[in] ps the atom's pose
///
core::Real ExactOccludedHbondSolEnergy::get_atom_lk_energy(
	core::Size const atom_idx,
	core::conformation::Residue const& res,
	core::pose::Pose const& ps
) const {

	core::Size residx = res.seqpos();

	core::scoring::EnergyGraph const & energy_graph( ps.energies().energy_graph() );

	core::Real lk_tot = 0;
	for ( utility::graph::Graph::EdgeListConstIter
			NITB = energy_graph.get_node(residx)->const_edge_list_begin(),
			NITE = energy_graph.get_node(residx)->const_edge_list_end(),
			nit = NITB; nit != NITE; ++nit ) {

		core::Size nri( ( *nit )->get_other_ind( residx ) );
		core::conformation::Residue const& occ_res = ps.residue(nri);
		lk_tot += get_atom_lk_energy(atom_idx, res, occ_res);
	}

	return lk_tot;
}


/// @brief returns the LK energy of a given atom due to a given residue
///
/// @param[in] atom_idx index of the atom in its residue
/// @param[in] res the atom's residue
/// @param[in] occ_res the occluding residue (i.e, the one to which the LK energy is due)
///
/// @details the code implements the same kind of approximation as
///  core::scoring::methods::LK_PolarNonPolarEnergy::get_residue_pair_energy_one_way()
///
core::Real ExactOccludedHbondSolEnergy::get_atom_lk_energy(
	core::Size atom_idx,
	conformation::Residue const& res,
	conformation::Residue const& occ_res) const {

	etable::count_pair::CountPairFunctionOP cpfxn =
		etable::count_pair::CountPairFactory::create_count_pair_function(res, occ_res, etable::count_pair::CP_CROSSOVER_4);

	core::conformation::Atom const& atom = res.atom(atom_idx);
	core::Vector atom_xyz = res.xyz(atom_idx);
	core::Real tot_lk_ene = 0;

	core::Size const NHVY = occ_res.nheavyatoms();
	for ( core::Size j = 1; j <= NHVY; ++j ) {

		core::Real cp_weight = 1.0;
		Size path_dist( 0 );
		if ( cpfxn->count( atom_idx, j, cp_weight, path_dist ) ) {

			core::Vector occ_atom_xyz = occ_res.xyz(j);
			core::Vector diff = atom_xyz - occ_atom_xyz;
			core::Real const d2 = diff.length_squared();
			if ( ( d2 > 0 ) && ( d2 < lk_safe_max_dis2_) ) {

				core::conformation::Atom const& occ_atom = occ_res.atom( j );
				core::Real lk_ene = 0;
				core::Real dummy_deriv = 0;
				etable_evaluator_->atom_pair_lk_energy_and_deriv_v( atom, occ_atom, lk_ene, dummy_deriv, false );
				tot_lk_ene += ( cp_weight*lk_ene );
			}
		}
	}

	return tot_lk_ene;
}


core::Size
ExactOccludedHbondSolEnergy::version() const
{
	return 1; // Initial versioning
}


/// @brief creates an ExactOccludedHbondSolEnergy object according to command-line options
///
/// @param[in] options options not from the command line
ExactOccludedHbondSolEnergyOP create_ExactSHOEnergy_from_cmdline(methods::EnergyMethodOptions const & options) {

	etable::EtableOptions etable_options = options.etable_options();
	etable_options.no_lk_polar_desolvation = false;
	etable_options.proline_N_is_lk_nonpolar = false;

	ExactOccludedHbondSolEnergyOP sho_op(
		new ExactOccludedHbondSolEnergy(
		*(ScoringManager::get_instance()->etable( etable_options ).lock()),
		options.analytic_etable_evaluation(),
		basic::options::option[ basic::options::OptionKeys::score::exact_occ_skip_Hbonders ],
		basic::options::option[ basic::options::OptionKeys::score::exact_occ_pairwise ],
		basic::options::option[ basic::options::OptionKeys::score::exact_occ_pairwise_by_res ],
		basic::options::option[ basic::options::OptionKeys::score::exact_occ_split_between_res ],
		! basic::options::option[ basic::options::OptionKeys::score::exact_occ_self_res_no_occ ],
		basic::options::option[basic::options::OptionKeys::score::exact_occ_radius_scaling]
		)
	);

	return sho_op;
}


} // geometric_solvation
} // scoring
} // core
