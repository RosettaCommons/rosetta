// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/DarcParticleSwarmMinimizer.cc
/// @brief
/// @author Karen R. Khar
/// @author Ragul Gowthaman


#include <protocols/pockets/DarcParticleSwarmMinimizer.hh>
#include <core/optimization/ParticleSwarmMinimizer.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <utility/exit.hh>
#include <numeric/random/random.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <ObjexxFCL/format.hh>
#include <algorithm>
#include <utility/vector1.hh>
#include <cmath>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>

namespace protocols {
namespace pockets {

using namespace ObjexxFCL::format;

// note: f_fitness is not used, rather we hard-code the objective function here to enable parellization across particles (for GPU)
//void DarcParticleSwarmMinimizer::score_all_particles(core::optimization::Multifunc & f_fitness, core::optimization::ParticleOPs & particles) {
void DarcParticleSwarmMinimizer::score_all_particles(core::optimization::Multifunc & , core::optimization::ParticleOPs & particles) {

	core::Size const N = particles.size();
	core::Size const nconformers = pfp_.compute_ligand_nconformers();

#ifdef USEOPENCL
	if(nfp_.gpu_.use()){
		std::vector<basic::gpu::float4> particles_rotation_offset(N);
		std::vector<basic::gpu::float4> particles_translation_offset(N);
		//In this GPU block, we compute only the arrays required to be passed for GPU calculations
		//No CPU computations goes in here
		for (Size j = 1; j <= N; ++j) {
			core::optimization::Multivec const vars = particles[j]->p_;//particles[particles_num]->p_;
			core::Size conformer=((core::Size)(floor(vars[7])) % nconformers);
			//fill array with translation offset; Note: particles indexed from 0.
			particles_translation_offset[j-1].x = vars[1];
			particles_translation_offset[j-1].y = vars[2];
			particles_translation_offset[j-1].z = vars[3];
			particles_translation_offset[j-1].w = conformer;//Here we pass the conformer number of ligand; indexed from 0
			//fill array with rotational angles
			particles_rotation_offset[j-1].x = vars[4];
			particles_rotation_offset[j-1].y = vars[5];
			particles_rotation_offset[j-1].z = vars[6];
			particles_rotation_offset[j-1].w = 0.0;
		}
		//call to GPU for DARC calculations
		nfp_.gpu_calculate_particle_scores(particles, particles_rotation_offset, particles_translation_offset);
	}
	else {
#endif

		//All the CPU calculations starts from here
		using namespace core::optimization;
		using namespace basic::options;
		ligand_natoms_shapecalc_ = pfp_.compute_ligand_natoms();
		ligand_natoms_elstscalc_ = pfp_.compute_ligand_natoms_with_hydrogens();
		core::Size num_origins_ = pfp_.num_origins();
		utility::vector1< std::vector<basic::gpu::float4> >atoms(num_origins_, std::vector<basic::gpu::float4>(ligand_natoms_shapecalc_ * N));
		utility::vector1< std::vector<basic::gpu::float4> >atom_maxmin_phipsi(num_origins_, std::vector<basic::gpu::float4>(ligand_natoms_shapecalc_ * N));
		std::vector<basic::gpu::float4> eatoms(ligand_natoms_elstscalc_ * N);
		//Initialize arrays requird for CPU calculations
		for (Size j = 1; j <= N; ++j) {
			core::optimization::Multivec const vars = particles[j]->p_;//particles[particles_num]->p_;
			numeric::xyzVector<core::Real> origin_offset;
			origin_offset.x() = vars[1];
			origin_offset.y() = vars[2];
			origin_offset.z() = vars[3];
			// note: conformer is indexed starting at 0
			core::Size conformer=((core::Size)(floor(vars[7])) % nconformers);
			core::conformation::ResidueCOP ligand_rsd = pfp_.select_conf_and_move_ligand_( nfp_, origin_offset, vars[4], vars[5], vars[6], conformer );
			fill_atom_arrays_(j-1, ligand_rsd, atoms, atom_maxmin_phipsi);
			if (!option[ OptionKeys::fingerprint::darc_shape_only ].user()){
				core::pose::Pose ligand_pose_for_elec_calc = pfp_.get_oriented_pose( nfp_, vars[4], vars[5], vars[6], origin_offset, conformer );
				fill_atom_arrays_for_electrostatics_(j-1, ligand_pose_for_elec_calc, eatoms );
			}
		}
		//CPU DARC calculations for each particle
		for (core::Size j = 1; j <= N; ++j) {
			core::Real particle_score(0.);
			particle_score = DarcPSO_fp_compare_( j-1, missing_pt_, steric_, extra_pt_, atoms, atom_maxmin_phipsi );
			//calculate Electorstatics score and add to DARC-shape-only score
			if (!option[ OptionKeys::fingerprint::darc_shape_only ].user()){
				std::vector < std::vector < std::vector <core::Real> > > espGrid = nfp_.espGrid_;
				std::vector < std::vector < std::vector <ElectrostaticpotentialGrid::PtType> > > typGrid = nfp_.typGrid_;
				core::Real spacing = nfp_.esp_spacing_;
				core::Size dim_x = nfp_.esp_dim_.x();
				core::Size dim_y = nfp_.esp_dim_.y();
				core::Size dim_z = nfp_.esp_dim_.z();
				core::Real mid_x = nfp_.esp_mid_.x();
				core::Real mid_y = nfp_.esp_mid_.y();
				core::Real mid_z = nfp_.esp_mid_.z();
				particle_score += DarcPSO_elsts_score_( j-1, dim_x, dim_y, dim_z, mid_x, mid_y, mid_z, spacing, espGrid, typGrid, eatoms );
			}
			//std::cout<<"Particle score(cpu): "<<particle_score<<std::endl;
			particles[j]->set_score(particle_score);
		}
		//END CPU DARC calculations
#ifdef USEOPENCL
  }
#endif
}

	void DarcParticleSwarmMinimizer::fill_atom_arrays_( core::Size particle_inx, core::conformation::ResidueCOP ligand_rsd, utility::vector1< std::vector<basic::gpu::float4> > & atoms, utility::vector1< std::vector<basic::gpu::float4> > & atom_maxmin_phipsi ) {
		using namespace basic::options;
		core::Real const radius_scale = option[ OptionKeys::fingerprint::atom_radius_scale ];
		core::Real const atom_buffer = option[ OptionKeys::fingerprint::atom_radius_buffer ];
		core::Size origin_num = 1;
		for (utility::vector1<numeric::xyzVector<core::Real> >::const_iterator i_mori = pfp_.multi_origin_list_.begin(); i_mori != pfp_.multi_origin_list_.end(); ++i_mori, ++origin_num) {

			for (Size i = 1; i <= ligand_natoms_shapecalc_; ++i) {
				core::Size const curr_array_inx = ( ligand_natoms_shapecalc_ * particle_inx ) + ( i-1 );
				numeric::xyzVector<core::Real> this_atomcoors = ligand_rsd->atom(i).xyz() - *i_mori;
				core::Real const this_atom_radius = ( ligand_rsd->atom_type(i).lj_radius() - atom_buffer ) * radius_scale;
				atoms[origin_num][curr_array_inx].x = this_atomcoors.x();
				atoms[origin_num][curr_array_inx].y = this_atomcoors.y();
				atoms[origin_num][curr_array_inx].z = this_atomcoors.z();
				atoms[origin_num][curr_array_inx].w = this_atom_radius;

				// find the atom center (in spherical coors)
				spherical_coor_triplet atom_center;
				convert_cartesian_to_spherical_coor_triplet( this_atomcoors, atom_center );
				core::Real const tmp_atomx=this_atomcoors.x();
				core::Real const tmp_atomy=this_atomcoors.y();
				core::Real const tmp_atomz=this_atomcoors.z();

				core::Real curr_max_phi( 999 );
				core::Real curr_min_phi( -999 );
				if ( std::abs(tmp_atomz) > 0.00001 ) {
					core::Real const inside_phi_asin = this_atom_radius / tmp_atomz;
					if ( std::abs(inside_phi_asin) < 1. ) {
						core::Real const max_angular_displacement_phi = std::abs( asin( inside_phi_asin ) );
						curr_max_phi = atom_center.phi + max_angular_displacement_phi;
						curr_min_phi = atom_center.phi - max_angular_displacement_phi;
					}
				}

				core::Real curr_max_psi( 999 );
				core::Real curr_min_psi( -999 );
				if ( ( std::abs(tmp_atomx) > 0.00001 ) || ( std::abs(tmp_atomx) > 0.00001 ) ) {
					core::Real const inside_psi_asin = this_atom_radius / sqrt( (tmp_atomx*tmp_atomx) + (tmp_atomy*tmp_atomy) );
					if ( std::abs(inside_psi_asin) < 1. ) {
						core::Real const max_angular_displacement_psi = std::abs( asin( inside_psi_asin ) );
						curr_max_psi = atom_center.psi + max_angular_displacement_psi;
						curr_min_psi = atom_center.psi - max_angular_displacement_psi;
					}
				}

				//We currently do not use ligand_maxmin_phipsi
				//if ( curr_max_phi > ligand_maxmin_phipsi[particle_inx].x ) ligand_maxmin_phipsi[particle_inx].x = curr_max_phi;
				//if ( curr_max_psi > ligand_maxmin_phipsi[particle_inx].y ) ligand_maxmin_phipsi[particle_inx].y = curr_max_psi;
				//if ( curr_min_phi < ligand_maxmin_phipsi[particle_inx].z ) ligand_maxmin_phipsi[particle_inx].z = curr_min_phi;
				//if ( curr_min_psi < ligand_maxmin_phipsi[particle_inx].w ) ligand_maxmin_phipsi[particle_inx].w = curr_min_psi;
				atom_maxmin_phipsi[origin_num][curr_array_inx].x = curr_max_phi;
				atom_maxmin_phipsi[origin_num][curr_array_inx].y = curr_max_psi;
				atom_maxmin_phipsi[origin_num][curr_array_inx].z = curr_min_phi;
				atom_maxmin_phipsi[origin_num][curr_array_inx].w = curr_min_psi;

			}
		}
	}

	void DarcParticleSwarmMinimizer::fill_atom_arrays_for_electrostatics_( core::Size particle_inx, core::pose::Pose ligand_pose_for_elec_calc, std::vector<basic::gpu::float4> & eatoms ) {

		core::Size lig_res_num = 0;
		for ( int j = 1, resnum = ligand_pose_for_elec_calc.total_residue(); j <= resnum; ++j ) {
			if (!ligand_pose_for_elec_calc.residue(j).is_protein()){
				lig_res_num = j;
				break;
			}
		}
		core::conformation::Residue const & ligand_rsd = ligand_pose_for_elec_calc.residue(lig_res_num);

		for (Size i = 1; i <= ligand_natoms_elstscalc_; ++i) {

		core::Size const curr_array_inx = ( ligand_natoms_elstscalc_ * particle_inx ) + ( i-1 );

    numeric::xyzVector<core::Real> this_atomcoors = ligand_rsd.atom(i).xyz();
    core::Real const this_atom_charge = ligand_rsd.atomic_charge(i);

    eatoms[curr_array_inx].x = this_atomcoors.x();
    eatoms[curr_array_inx].y = this_atomcoors.y();
    eatoms[curr_array_inx].z = this_atomcoors.z();
    eatoms[curr_array_inx].w = this_atom_charge;
		}
	}

	core::Real
	DarcParticleSwarmMinimizer::DarcPSO_fp_compare_(
																									core::Size particle_inx,
																									core::Real const & missing_point_weight,
																									core::Real const & steric_weight,
																									core::Real const & extra_point_weight,
																									utility::vector1< std::vector<basic::gpu::float4> > & atoms,
																									utility::vector1< std::vector<basic::gpu::float4> > & atom_maxmin_phipsi ) {

		core::Real Total_score = 0;
		core::Size num_rays = 0;
		core::Real underpack_dist = 0, steric_dist = 0, rays_missing_ligand = 0, rays_missing_pocket = 0;

		for (std::list<spherical_coor_triplet>::const_iterator ni = nfp_.triplet_fingerprint_data().begin();
				 ni != nfp_.triplet_fingerprint_data().end(); ++ni) {

			core::Real curr_phi = ni->phi;
			core::Real curr_psi = ni->psi;
			core::Size curr_ori = ni->ori;
			core::Real best_rho_sq(9999.);
			//    core::Size best_intersecting_atom(0);
			for (Size i = 1, i_end = ligand_natoms_shapecalc_; i <= i_end; ++i) {
				core::Size const curr_array_inx = ( ligand_natoms_shapecalc_ * particle_inx ) + ( i-1 );

				if ( atoms[curr_ori][curr_array_inx].w < 0.001 ) continue;

				while ( curr_phi < atom_maxmin_phipsi[curr_ori][curr_array_inx].z ) {
					curr_phi += numeric::constants::r::pi_2;
				}

				while ( curr_phi > atom_maxmin_phipsi[curr_ori][curr_array_inx].x ) {
					curr_phi -= numeric::constants::r::pi_2;
				}
				if ( curr_phi < atom_maxmin_phipsi[curr_ori][curr_array_inx].z ) continue;
				if ( curr_phi > atom_maxmin_phipsi[curr_ori][curr_array_inx].x ) continue;

				while ( curr_psi < atom_maxmin_phipsi[curr_ori][curr_array_inx].w ) {
					curr_psi += numeric::constants::r::pi_2;
				}
				while ( curr_psi > atom_maxmin_phipsi[curr_ori][curr_array_inx].y ) {
					curr_psi -= numeric::constants::r::pi_2;
				}
				if ( curr_psi < atom_maxmin_phipsi[curr_ori][curr_array_inx].w ) continue;
				if ( curr_psi > atom_maxmin_phipsi[curr_ori][curr_array_inx].y ) continue;

				core::Real const min_intersect_SQ = Find_Closest_Intersect_SQ(curr_phi, curr_psi, atoms[curr_ori][curr_array_inx].x, atoms[curr_ori][curr_array_inx].y, atoms[curr_ori][curr_array_inx].z, atoms[curr_ori][curr_array_inx].w);

				if ( min_intersect_SQ < best_rho_sq ) {
					best_rho_sq = min_intersect_SQ;
				}
			}
			core::Real plaid_rho = 9999.;
			if ( best_rho_sq < 9998. ) {
				plaid_rho = sqrt(best_rho_sq);
			}

			if ( (plaid_rho > 9998.) && (ni->rho > 0.001) ) {
				Total_score += missing_point_weight;
				rays_missing_ligand++;
				num_rays++;
			}
			if ( (plaid_rho < 9999.) && (ni->rho < 0.001 ) ) {
				Total_score += extra_point_weight;
				rays_missing_pocket++;
				num_rays++;
			}
			if ( (plaid_rho < 9999.) && (ni->rho > 0.001 ) ) {
				core::Real distance_deviation = std::abs( ni->rho - plaid_rho );
				if (plaid_rho > ni->rho) underpack_dist += distance_deviation;
				if (plaid_rho < ni->rho) steric_dist += distance_deviation;
				core::Real dist_deviation = ( plaid_rho - ni->rho );
				if (dist_deviation < 0.0) dist_deviation = ( ni->rho - plaid_rho ) * steric_weight;
				Total_score += dist_deviation;
				num_rays++;
			}
		}
		using namespace basic::options;
		bool get_darc_components = option[ OptionKeys::fingerprint::darc_components ]();
		if (get_darc_components){
			std::cout<<"OPT "<<underpack_dist<<" "<<steric_dist<< " "<<rays_missing_ligand<<" "<<rays_missing_pocket<<" "<<num_rays;
			//std::cout<<"CMPSCORE "<<" "<<Total_score/num_rays<<" "<<(underpack_dist+(steric_dist*steric_weight)+(rays_missing_ligand*missing_point_weight)+(rays_missing_pocket*extra_point_weight))/num_rays<<std::endl;
		}
		return (Total_score/num_rays);
	}

core::Real
DarcParticleSwarmMinimizer::DarcPSO_elsts_score_(
		core::Size particle_inx,
		core::Size dim_x,
		core::Size dim_y,
		core::Size dim_z,
		core::Real mid_x,
		core::Real mid_y,
		core::Real mid_z,
		core::Real spacing,
		std::vector < std::vector < std::vector <core::Real> > > espGrid,
		std::vector < std::vector < std::vector <ElectrostaticpotentialGrid::PtType> > > typGrid,
		std::vector<basic::gpu::float4> & atom_coors_charge ) {

	using namespace basic::options;
  core::Real const esp_wt = option[ OptionKeys::fingerprint::esp_weight ];

	core::Real curr_xcoor, curr_ycoor, curr_zcoor, curr_charge;
	core::Real lig_elsts_energy(0.);

	for (core::Size i = 1, i_end = ligand_natoms_elstscalc_; i <= i_end; ++i) {
		core::Real atm_esp_energy(0.);
		core::Size const curr_array_inx = ( ligand_natoms_elstscalc_ * particle_inx ) + ( i-1 );

		curr_xcoor = atom_coors_charge[curr_array_inx].x;
		curr_ycoor = atom_coors_charge[curr_array_inx].y;
		curr_zcoor = atom_coors_charge[curr_array_inx].z;
		curr_charge = atom_coors_charge[curr_array_inx].w;

		core::Real X = (curr_xcoor - ( mid_x - ((static_cast<core::Real>(dim_x-1)/2) * spacing) ))/ spacing;
		core::Real Y = (curr_ycoor - ( mid_y - ((static_cast<core::Real>(dim_y-1)/2) * spacing) ))/ spacing;
		core::Real Z = (curr_zcoor - ( mid_z - ((static_cast<core::Real>(dim_z-1)/2) * spacing) ))/ spacing;

		core::Real X1 = std::floor(X+1);
		core::Real Y1 = std::floor(Y+1);
		core::Real Z1 = std::floor(Z+1);

		core::Real X0 = std::ceil(X-1);
		core::Real Y0 = std::ceil(Y-1);
		core::Real Z0 = std::ceil(Z-1);

		//if the ligand atom goes out of the grid
		if( ((X0 >= dim_x)||(X0 < 0)) || ((Y0 >= dim_y)||(Y0 < 0)) || ((Z0 >= dim_z)||(Z0 < 0)) || ((X1 >= dim_x)||(X1 < 0)) || ((Y1 >= dim_y)||(Y1 < 0)) || ((Z1 >= dim_z)||(Z1 < 0)) ) {
			atm_esp_energy = 100.0;
		} else {
			core::Real Xd = (X-X0)/(X1-X0);
			core::Real Yd = (Y-Y0)/(Y1-Y0);
			core::Real Zd = (Z-Z0)/(Z1-Z0);


			core::Real C00 = espGrid[Size(X0)][Size(Y0)][Size(Z0)]*(1-Xd) + espGrid[Size(X1)][Size(Y0)][Size(Z0)]*Xd;
			core::Real C10 = espGrid[Size(X0)][Size(Y1)][Size(Z0)]*(1-Xd) + espGrid[Size(X1)][Size(Y1)][Size(Z0)]*Xd;
			core::Real C01 = espGrid[Size(X0)][Size(Y0)][Size(Z1)]*(1-Xd) + espGrid[Size(X1)][Size(Y0)][Size(Z1)]*Xd;
			core::Real C11 = espGrid[Size(X0)][Size(Y1)][Size(Z1)]*(1-Xd) + espGrid[Size(X1)][Size(Y1)][Size(Z1)]*Xd;
			core::Real C0 = C00*(1-Yd) + C10*Yd;
			core::Real C1 = C01*(1-Yd) + C11*Yd;
			core::Real C = C0*(1-Zd) + C1*Zd;

			atm_esp_energy = C * curr_charge;

			if ( ( typGrid[Size(X0)][Size(Y0)][Size(Z0)] == ElectrostaticpotentialGrid::PROTEIN ) ||
					 ( typGrid[Size(X0)][Size(Y1)][Size(Z0)] == ElectrostaticpotentialGrid::PROTEIN ) ||
					 ( typGrid[Size(X0)][Size(Y0)][Size(Z1)] == ElectrostaticpotentialGrid::PROTEIN ) ||
					 ( typGrid[Size(X0)][Size(Y1)][Size(Z1)] == ElectrostaticpotentialGrid::PROTEIN ) ||
					 ( typGrid[Size(X1)][Size(Y0)][Size(Z0)] == ElectrostaticpotentialGrid::PROTEIN ) ||
					 ( typGrid[Size(X1)][Size(Y1)][Size(Z0)] == ElectrostaticpotentialGrid::PROTEIN ) ||
					 ( typGrid[Size(X1)][Size(Y0)][Size(Z1)] == ElectrostaticpotentialGrid::PROTEIN ) ||
					 ( typGrid[Size(X1)][Size(Y1)][Size(Z1)] == ElectrostaticpotentialGrid::PROTEIN ) ) {
				// return zero if the energy is favourable and the atom goes into protein
				if (atm_esp_energy < 0.) atm_esp_energy = 0.;
			}
		}
		lig_elsts_energy += atm_esp_energy;
	}
	using namespace basic::options;
	bool get_darc_components = option[ OptionKeys::fingerprint::darc_components ]();
	if (get_darc_components){
		std::cout<<" "<<lig_elsts_energy<<std::endl;
	}

	return lig_elsts_energy * esp_wt;
}


} // namespace pockets
} // namespace protocols
