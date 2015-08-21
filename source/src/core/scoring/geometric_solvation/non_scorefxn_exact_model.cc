// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk

// Project Headers
#include <core/scoring/geometric_solvation/ExactOccludedHbondSolEnergy.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Atom.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <basic/Tracer.hh>

#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/types.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>

//#include <core/scoring/ScoreFunction.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>


// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <cmath>
#include <iostream>
#include <iomanip>
#include <map>

#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>
#include <numeric/xyz.functions.hh>


//Vector dummy_res_energy_vector_;

static thread_local basic::Tracer TR( "core.scoring.geometric_solvation.exact_model" );

namespace core {
namespace scoring {
namespace geometric_solvation {

using namespace core;
using namespace core::scoring;
using namespace core::scoring::hbonds;

core::Real const geosol_kT = { 0.593 };

// apply this weight to everything, so that scale will match LK
core::Real const LK_MATCHING_WEIGHT_OLD_EXACT = { 0.23968 };


void add_to_individual_sol_energies(
	pose::Pose & input_pose,
	core::Size const polar_resnum,
	core::Size const polar_atomno,
	core::scoring::etable::EtableOP etable_ptr,
	GridInfo const & grid_info,
	core::Real const & grid_constant,
	std::vector < std::vector < std::vector <core::Real> > > const & water_weights,
	std::vector < std::vector < std::vector <bool> > > & occluded_sites,
	bool const hydrogens_can_occlude,
	bool const pairwise_additive,
	bool const pairwise_additive_output,
	utility::vector1 <core::Real> & residue_energies ) {

	core::Real const water_radius = 1.4;

	if ( pairwise_additive_output && ! pairwise_additive ) {
		TR << "Error - pairwise additive output doesn't make sense when calculations are not pairwise additive!!" << std::endl;
		exit(1);
	}

	// Reset grid of occluded sites, by setting everything to false (for "not occluded")
	for ( core::Size tx=0; tx<grid_info.xnum_points(); tx++ ) {
		for ( core::Size ty=0; ty<grid_info.ynum_points(); ty++ ) {
			for ( core::Size tz=0; tz<grid_info.znum_points(); tz++ ) {
				occluded_sites[tx][ty][tz] = false;
			}
		}
	}

	// Find the transformation which puts the donor/acceptor of interest onto the existing grid
	// The plan is to apply this transformation to bring each occluding atom onto the existing grid (translation then matrix multiplication)
	core::Size const base_atomno( input_pose.residue(polar_resnum).atom_base( polar_atomno ) );
	core::Vector const & orig_polar_atom_xyz( input_pose.residue(polar_resnum).atom( polar_atomno ).xyz() );
	core::Vector const & orig_base_atom_xyz( input_pose.residue(polar_resnum).atom( base_atomno ).xyz() );
	core::Vector translation_vector = -1.0 * orig_polar_atom_xyz;
	core::Vector translated_base_atom = orig_base_atom_xyz + translation_vector;

	// We want to translate positions of occluding atoms _from_ a cartesian basis set into one using the reference frame of the polar group
	core::Vector cartesian_x(1,0,0);
	core::Vector cartesian_y(0,1,0);
	core::Vector cartesian_z(0,0,1);

	// There's not a unique solution for this, so we'll arbitrarily pick a second basis vector, requiring that the dot product with desired_z is zero
	core::Vector desired_z = -1. * translated_base_atom.normalized();

	// Treat the case where polar_atom.z == base_atom.z ; this leads to desired_z.z == 0
	core::Real arbitrary_x, arbitrary_y, arbitrary_z;
	if ( std::abs( desired_z.z() ) > 0.01 ) {
		arbitrary_x = 1;
		arbitrary_y = 1;
		arbitrary_z = -1. * ((arbitrary_x * desired_z.x()) + (arbitrary_y * desired_z.y())) / desired_z.z();
	} else {
		arbitrary_x = 1;
		arbitrary_z = 1;
		arbitrary_y = -1. * ((arbitrary_x * desired_z.x()) + (arbitrary_z * desired_z.z())) / desired_z.y();
	}
	core::Vector desired_x(arbitrary_x,arbitrary_y,arbitrary_z);
	desired_x.normalize();
	core::Vector desired_y = cross_product( desired_x, desired_z );

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

	// Loop over all atoms of neighboring residues, INCLUDING SELF
	core::scoring::TenANeighborGraph const & graph = input_pose.energies().tenA_neighbor_graph();
	utility::vector1 <core::Size> neighborlist;
	neighborlist.push_back( polar_resnum);
	for ( core::graph::Graph::EdgeListConstIter
			neighbor_iter = graph.get_node( polar_resnum )->const_edge_list_begin(),
			neighbor_iter_end = graph.get_node( polar_resnum )->const_edge_list_end();
			neighbor_iter != neighbor_iter_end; ++neighbor_iter ) {
		neighborlist.push_back( (*neighbor_iter)->get_other_ind( polar_resnum ) );
	}

	for ( Size occ_inx = 1; occ_inx <= neighborlist.size(); ++occ_inx ) {
		core::Size const occ_resnum( neighborlist[occ_inx] );
		conformation::Residue const occ_rsd = input_pose.residue(occ_resnum);
		//  TR << "jk computing occlusion of polar residue " << polar_resnum << " by residue " << occ_resnum << std::endl;

		for ( Size occ_atomno = 1; occ_atomno <= occ_rsd.natoms(); ++occ_atomno ) {

			bool const occ_atom_is_hydrogen = occ_rsd.atom_is_hydrogen( occ_atomno );
			if ( occ_atom_is_hydrogen && ! hydrogens_can_occlude ) continue;

			// can be occluded by atoms directly bonded to this group, but not by self
			if ( polar_resnum == occ_resnum ) {
				if ( polar_atomno == occ_atomno ) continue;
				if ( base_atomno == occ_atomno ) continue;
			}

			// If pairwise-additive, reset grid of occluded sites for every atom
			if ( pairwise_additive ) {
				for ( core::Size tx=0; tx<grid_info.xnum_points(); tx++ ) {
					for ( core::Size ty=0; ty<grid_info.ynum_points(); ty++ ) {
						for ( core::Size tz=0; tz<grid_info.znum_points(); tz++ ) {
							occluded_sites[tx][ty][tz] = false;
						}
					}
				}
			}

			// Apply the transformation to put this atom onto the current grid
			core::Vector const & orig_occ_atom_xyz( occ_rsd.atom( occ_atomno ).xyz() );
			core::Vector const translated_occ_atom_xyz = orig_occ_atom_xyz + translation_vector;
			core::Vector const transformed_occ_atom_xyz = transformation_matrix * ( orig_occ_atom_xyz + translation_vector );

			// Double-check transformations
			debug_assert( std::abs(orig_polar_atom_xyz.distance( orig_occ_atom_xyz ) - transformed_occ_atom_xyz.magnitude()) < 0.001 );
			debug_assert( std::abs(orig_base_atom_xyz.distance( orig_occ_atom_xyz ) - new_base_atom_location.distance( transformed_occ_atom_xyz )) < 0.001 );
			// Find sites occluded by this atom, set these to true
			core::Real const occ_radius = etable_ptr->lj_radius( occ_rsd.atom_type_index( occ_atomno ) );
			core::Real const sq_dist_cut = ( occ_radius + water_radius ) * ( occ_radius + water_radius );

			// Loop over all water positions, mark those which are occluded
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
							occluded_sites[wx][wy][wz] = true;
						}
					}
				}
			}

			if ( pairwise_additive ) {
				// For pairwise additive version, compute and save energy for every occluding atom separately
				core::Real sum_occluded_weights(0.);
				for ( core::Size tx=0; tx<grid_info.xnum_points(); tx++ ) {
					for ( core::Size ty=0; ty<grid_info.ynum_points(); ty++ ) {
						for ( core::Size tz=0; tz<grid_info.znum_points(); tz++ ) {
							if ( occluded_sites[tx][ty][tz] ) {
								core::Real const curr_water_weight = water_weights[tx][ty][tz];
								sum_occluded_weights += curr_water_weight;
							}
						}
					}
				}
				core::Real const geometric_solvation_energy = - geosol_kT * log( 1 - ( sum_occluded_weights / grid_constant ) );
				if ( pairwise_additive_output ) {
					// For pairwise additive output, split the energy between polar and occluding residue (eg. to match form used in scorefxn)
					residue_energies[ polar_resnum ] += geometric_solvation_energy / 2.;
					residue_energies[ occ_resnum ] += geometric_solvation_energy / 2.;
				} else {
					// For testing, write output in the non-pairwise additive format (eg. to compare to non-pairwise additive model)
					residue_energies[ polar_resnum ] += geometric_solvation_energy;
				}
			}

		}
	}

	if ( ! pairwise_additive ) {
		// Compute and store the solvation energy for the grid occluded by all nearby atoms
		// Compute the numerator (the sum of occluded weights)
		core::Real sum_occluded_weights(0.);
		for ( core::Size tx=0; tx<grid_info.xnum_points(); tx++ ) {
			for ( core::Size ty=0; ty<grid_info.ynum_points(); ty++ ) {
				for ( core::Size tz=0; tz<grid_info.znum_points(); tz++ ) {
					if ( occluded_sites[tx][ty][tz] ) {
						core::Real const curr_water_weight = water_weights[tx][ty][tz];
						sum_occluded_weights += curr_water_weight;
					}
				}
			}
		}
		// If we have a non-pairwise additive model, we can't split the energy between the polar and occluding residues
		core::Real const geometric_solvation_energy = - geosol_kT * log( 1 - ( sum_occluded_weights / grid_constant ) );
		residue_energies[ polar_resnum ] += geometric_solvation_energy;
	}

	return;
}


core::Real compute_exact_geosol(
	pose::Pose & input_pose,
	bool const hydrogens_can_occlude,
	bool const pairwise_additive,
	bool const pairwise_additive_output,
	utility::vector1<core::Real> & residue_energies ) {

	TR << "jk geometric solvation exact scoring" << std::endl;

	residue_energies.clear();
	residue_energies.resize( input_pose.total_residue(), 0.);

	// Get a copy of the Etable (to lookup atomic radii)
	core::scoring::etable::EtableOP etable_ptr( new core::scoring::etable::Etable( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ), core::scoring::etable::EtableOptions() ) );

	// Allocate memory for grid of occluded sites
	std::vector < std::vector < std::vector <bool> > > occluded_sites;
	occluded_sites.clear();
	occluded_sites.resize(GridInfo::get_instance()->xnum_points());
	for ( core::Size tx=0; tx<GridInfo::get_instance()->xnum_points(); tx++ ) {
		occluded_sites[tx].resize(GridInfo::get_instance()->ynum_points());
		for ( core::Size ty=0; ty<GridInfo::get_instance()->ynum_points(); ty++ ) {
			occluded_sites[tx][ty].resize(GridInfo::get_instance()->znum_points());
		}
	}

	// get exact geometric solvation scores as a fxn of residue number
	TR << "jk computing exact solvation scores" << std::endl;
	for ( Size polar_resnum = 1; polar_resnum <= input_pose.total_residue(); polar_resnum++ ) {

		conformation::Residue const polar_rsd = input_pose.residue(polar_resnum);

		// loop over donors in polar_rsd
		for ( chemical::AtomIndices::const_iterator
				hnum  = polar_rsd.Hpos_polar().begin(),
				hnume = polar_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			Size const polar_atom( *hnum );
			Size const base_atom( polar_rsd.atom_base( polar_atom ) );
			hbonds::HBEvalTuple const curr_hbond_eval_tuple(
				get_hb_don_chem_type( polar_atom, polar_rsd ),
				hbacc_H2O, seq_sep_other);

			// Figure out max LK energy
			std::string const base_atom_name = polar_rsd.atom_name( base_atom );
			core::Real max_possible_LK = etable_ptr->lk_dgfree( polar_rsd.atom_type_index( base_atom ) );
			if ( ( base_atom_name == " N  " ) && polar_rsd.is_lower_terminus() ) max_possible_LK /= 3; // charged N-terminus
			if ( base_atom_name == " NZ " ) max_possible_LK /= 3; // Lys
			if ( base_atom_name == " ND2" ) max_possible_LK /= 2; // Asn
			if ( base_atom_name == " NE2" ) max_possible_LK /= 2; // Gln
			if ( base_atom_name == " NH1" ) max_possible_LK /= 2; // Arg
			if ( base_atom_name == " NH2" ) max_possible_LK /= 2; // Arg
			// Note: inner nitrogen of Arg (NE) is extra strong, since it's the same atom type as the other two but doesn't get
			// cut in half because there's only one proton...
			//   TR << "jk max LK for donor with base " << base_atom_name << " is  " << max_possible_LK << std::endl;

			// Compute Ebulk (using the LK energy)
			core::Real const Emax_weight = exp( max_possible_LK / geosol_kT );
			core::Real const sum_water_weights = WaterWeightGridSet::get_instance()->get_sum_water_weight_grid( curr_hbond_eval_tuple.eval_type() );
			core::Real const Ebulk_weight = ( sum_water_weights * Emax_weight ) / ( 1. - Emax_weight);
			// This grid constant is the denominator in computing solvation energies,
			// it depends on the grid dimensions, and sets the max possible solvation energy (in this case to match LK)
			core::Real const grid_constant = sum_water_weights + Ebulk_weight;
			// Setup then call compute_individual_sol_energies
			add_to_individual_sol_energies(input_pose, polar_resnum, polar_atom, etable_ptr, *GridInfo::get_instance(), grid_constant,
				WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbond_eval_tuple.eval_type() ), occluded_sites, hydrogens_can_occlude,
				pairwise_additive, pairwise_additive_output, residue_energies );
		}

		// loop over acceptors in polar_rsd
		for ( chemical::AtomIndices::const_iterator
				anum  = polar_rsd.accpt_pos().begin(),
				anume = polar_rsd.accpt_pos().end(); anum != anume; ++anum ) {
			Size const polar_atom( *anum );
			//Size const base_atom ( polar_rsd.atom_base( polar_atom ) );
			hbonds::HBEvalType const curr_hbeval_type = hbonds::HBEval_lookup( hbdon_H2O, get_hb_acc_chem_type( polar_atom, polar_rsd ), seq_sep_other);

			// Figure out max LK energy
			//std::string const base_atom_name = polar_rsd.atom_name( base_atom );
			core::Real max_possible_LK = etable_ptr->lk_dgfree( polar_rsd.atom_type_index( polar_atom ) );
			//   TR << "jk max LK for acceptor " << polar_rsd.atom_name(polar_atom) << " is  " << max_possible_LK << std::endl;
			// Compute Ebulk (using the LK energy)
			core::Real const Emax_weight = exp( max_possible_LK / geosol_kT );
			core::Real const sum_water_weights = WaterWeightGridSet::get_instance()->get_sum_water_weight_grid( curr_hbeval_type );
			core::Real const Ebulk_weight = ( sum_water_weights * Emax_weight ) / ( 1. - Emax_weight);
			// This grid constant is the denominator in computing solvation energies,
			// it depends on the grid dimensions, and sets the max possible solvation energy (in this case to match LK)
			core::Real const grid_constant = sum_water_weights + Ebulk_weight;
			// Setup then call compute_individual_sol_energies
			add_to_individual_sol_energies(input_pose, polar_resnum, polar_atom, etable_ptr, *GridInfo::get_instance(), grid_constant,
				WaterWeightGridSet::get_instance()->get_water_weight_grid( curr_hbeval_type ), occluded_sites, hydrogens_can_occlude,
				pairwise_additive, pairwise_additive_output, residue_energies );
		}

	}

	TR << "jk finished computing exact geometric solvation scores" << std::endl;

	core::Real total_solvation_energy(0.);
	for ( Size i = 1; i <= input_pose.total_residue(); i++ ) {
		residue_energies[i] *= LK_MATCHING_WEIGHT_OLD_EXACT;
		total_solvation_energy += residue_energies[i];
	}

	return total_solvation_energy;
}


} // geometric_solvation
} // scoring
} // core

