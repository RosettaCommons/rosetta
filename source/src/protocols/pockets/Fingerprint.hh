// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pockets/MultiFingerprint.hh
/// @brief  protocols::pockets::MultiFingerprint header
/// @author Ragul Gowthaman

#ifndef INCLUDED_protocols_pockets_Fingerprint_hh
#define INCLUDED_protocols_pockets_Fingerprint_hh

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/pockets/Fingerprint.fwd.hh>
#include <protocols/pockets/PocketGrid.fwd.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <protocols/pockets/DarcParticleSwarmMinimizer.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1_bool.hh>
#include <list>
#include <cmath>
#include <iostream>
#include <utility/vector1.hh>
#include <basic/gpu/GPU.hh>

namespace protocols {
namespace pockets {

struct spherical_coor_triplet {
	// note: these are in radians
	core::Real phi;
	core::Real psi;
	core::Real rho;
	core::Size ori;
};

struct triplet_and_originnum {
	// note: these are in radians
	core::Real phi;
	core::Real psi;
	core::Real rho;
	core::Size originnum;
};

struct ray_distance_derivs {
	core::Real dDist_dv1;
	core::Real dDist_dv2;
	core::Real dDist_dv3;
	core::Real dDist_dv4;
	core::Real dDist_dv5;
	core::Real dDist_dv6;
};

class FingerprintBase : public utility::pointer::ReferenceCount {

	friend class FingerprintMultifunc;
	friend class DarcParticleSwarmMinimizer;

public:

	FingerprintBase();

	void print_to_file(std::string const & output_filename) const;
	//void print_to_file(std::string const & output_filename, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset) const;

	void print_to_pdb(std::string const & output_pdbname) const;

	void print_to_pdb(std::string const & output_pdbname, numeric::xyzVector<core::Real> const & translation ) const;

	// const accessor functions

	core::Size num_origins() const { return num_origins_; };
	numeric::xyzVector<core::Real> origin() const { return origin_; };
	utility::vector1< numeric::xyzVector<core::Real> > multi_origin_list() const { return multi_origin_list_; };

	numeric::xyzVector<core::Real> CoM() const { return CoM_; };
	void CHEAT_CoM( numeric::xyzVector<core::Real> const & inp_CoM ) { CoM_ = inp_CoM; }; // CHEAT!!

	std::list< spherical_coor_triplet > const & triplet_fingerprint_data() const { return triplet_fingerprint_data_; };

	numeric::xyzVector<core::Real> pocketGrid_mid_;
	numeric::xyzVector<core::Real> pocketGrid_dim_;
	core::Real pocketGrid_spacing_;


protected:
	core::Size num_origins_;
	numeric::xyzVector<core::Real> origin_;
	utility::vector1< numeric::xyzVector<core::Real> > multi_origin_list_;
	numeric::xyzVector<core::Real> CoM_;
	std::list< spherical_coor_triplet > triplet_fingerprint_data_;

public:

	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~FingerprintBase();

};

class NonPlaidFingerprint : public FingerprintBase {
public:
	NonPlaidFingerprint();
	~NonPlaidFingerprint();
	numeric::xyzVector<core::Real> pocket_CoM_;
	std::list< numeric::xyzVector<core::Real> > egg_and_ext_list_;
	std::list< numeric::xyzVector<core::Real> > eggshell_list_;
	std::list< numeric::xyzVector<core::Real> > extshell_list_;
	void change_CoM_to_ligandCoM(numeric::xyzVector<core::Real> const & ligandCoM);

	void setup_from_espGrid( std::string const & input_filename, core::pose::Pose const & protein_pose, bool delphi );

	void setup_from_PocketGrid( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid );

	void setup_from_PocketGrid( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, PocketGrid const & grid_for_extshell );

	void setup_from_PocketGrid_and_known_ligand( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, PocketGrid const & grid_for_extshell, core::pose::Pose const & known_ligand_pose, core::Real const & trim_dist );

	void setup_from_PocketGrid_using_bound_ligand( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, PocketGrid const & grid_for_extshell, core::pose::Pose const & known_ligand_pose );

	void setup_from_EggshellGrid();

	void setup_from_ConnollySurface( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, PocketGrid const & grid_for_extshell );

	void write_eggshell_to_pdb_file( std::string const & output_eggshell_name ) const;


	void setup_from_eggshell_pdb_file( std::string const & input_filename);

	void trim_based_on_known_ligand( core::pose::Pose const & known_ligand_pose );

	void include_eggshell_points_based_on_known_ligand( core::pose::Pose const & known_ligand_pose, core::Real const & trim_dist);

	void setup_from_eggshell_triplet_file(std::string const & input_filename);

	void setup_from_PlaidFingerprint( PlaidFingerprint const & pfp );

	void choose_origin_with_lowest_eggshell_ruggedness( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell );

	void set_origin ( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell );

	void set_origin_from_option_( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, Size const & set_origin_option );

	void set_multiple_origin( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell );

	numeric::xyzVector<core::Real> place_origin_point( core::Real const & angle );

	void set_origin_away_from_protein_center ( core::pose::Pose const & protein_pose );

	void set_origin_from_residue ( core::pose::Pose const & protein_pose );

	core::Size get_pose_resnum(int const pdbnum, char const pdbchn, core::pose::Pose const & ps);

	void set_origin_away_from_eggshell( std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, core::pose::Pose const & protein_pose );

	void set_origin_away_from_eggshell_plane( std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, core::pose::Pose const & protein_pose, Size const & set_origin_option );

	core::Real get_Rvalue (core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, Size const & set_origin_option);

	numeric::xyzVector<core::Real> calculate_protein_CoM( core::pose::Pose const & protein_pose);

	std::list< numeric::xyzVector<core::Real> > combine_xyz_lists (std::list< numeric::xyzVector<core::Real> > const & xyz_list_1 , std::list< numeric::xyzVector<core::Real> > const & xyz_list_2);

	std::list<spherical_coor_triplet> remove_duplicate_phi_psi(std::list<spherical_coor_triplet> const & rounded_triplet);

	std::list<spherical_coor_triplet> convert_cart_to_spherical_and_round(std::list< numeric::xyzVector<core::Real> > const & xyz_list);

	std::list<numeric::xyzVector<core::Real> > convert_spherical_list_to_cartesian_list(std::list<spherical_coor_triplet> const & unique_triplet);

	std::list<spherical_coor_triplet> set_rho_to_zero(std::list<spherical_coor_triplet> const & rounded_triplet);

	core::Real get_electrostatics_energy( core::pose::Pose const & ligand_pose);


	core::Real get_nearest_neighbour_esp_energy(numeric::xyzVector<core::Real> const & ligand_atom, core::Real const & atom_charge);

	core::Real get_interpolated_esp_energy(numeric::xyzVector<core::Real> const & ligand_atom, core::Real const & atom_charge);

	core::Real get_interpolated_esp_energy_with_type(numeric::xyzVector<core::Real> const & ligand_atom, core::Real const & atom_charge);

	core::Real get_surface_esp( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list, std::string const & inp_espGrid_fname );

	std::vector < std::vector < std::vector <core::Real> > > espGrid_;
	std::vector < std::vector < std::vector <ElectrostaticpotentialGrid::PtType> > > typGrid_;
	core::Real esp_spacing_;
	numeric::xyzVector<core::Real> esp_mid_;
	numeric::xyzVector<core::Size> esp_dim_;

#ifdef USEOPENCL
	public:
    basic::gpu::GPU gpu_;
	private:
    struct {
			cl_mem rays;
			cl_mem RAYorigins;
			cl_mem atomcoords_shapecalc;
			cl_mem atomcoords_elstscalc;
			cl_mem particles_rotation_offset;
			cl_mem particles_translation_offset;
			cl_mem ray_scores;
			cl_mem raytest_scores;
			cl_mem particle_scores;
			cl_mem particletest_scores;
			cl_mem weights;
			cl_mem ligCoM;
			unsigned int rays_size;
			unsigned int ligCoM_size;
			unsigned int RAYorigins_size;
			unsigned int atomcoords_shapecalc_size;
			unsigned int atomcoords_elstscalc_size;
			unsigned int particles_translation_offset_size;
			unsigned int particles_rotation_offset_size;
			unsigned int ray_scores_size;
			unsigned int raytest_scores_size;
			unsigned int particle_scores_size;
			unsigned int particletest_scores_size;
			int num_rays;

			int num_ligatoms_elstscalc;
			int num_totatoms_elstscalc;
			int num_ligatoms_shapecalc;
			int num_totatoms_shapecalc;
			int num_particles;
			int num_eatoms;
			cl_mem elec_weights;
			cl_mem elec_scores;
			unsigned int elec_scores_size;
			cl_mem gpu_espGrid;
			unsigned int gpu_espGrid_size;
			cl_mem gpu_typGrid;
			unsigned int gpu_typGrid_size;
			cl_mem eatoms;
			unsigned int eatoms_size;
			cl_mem grid_mid;
			cl_mem grid_dim;

    } gpu_memory_;

	public:

		void gpu_setup( core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, int & num_particles, core::Real const & electrostatics_weight );
    int gpu_calculate_particle_scores(core::optimization::ParticleOPs & particles, std::vector<basic::gpu::float4> & particles_rotation_offset, std::vector<basic::gpu::float4> & particles_translation_offset);
    void gpu_setup_rays();
		void gpu_setup_atomcoords(PlaidFingerprint & pf, int & num_particles);
    void gpu_setup_espGrid();

#endif

};

class PlaidFingerprint : public FingerprintBase {

	friend class FingerprintMultifunc;
	friend class DarcParticleSwarmMinimizer;

public:

	PlaidFingerprint( core::pose::Pose const & input_pose, FingerprintBase & fp );

	core::Real find_optimal_rotation( FingerprintBase & fp, core::Real const & angle_increment, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight );

	core::Real find_optimal_rotation( FingerprintBase & fp, core::Real const & angle_increment, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, numeric::xyzVector<core::Real> const & no_CoM_offset );

	core::Real search_random_poses( FingerprintBase & fp, core::Size const & num_pose_search, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight );

	core::Real search_random_poses( FingerprintBase & fp, core::Size const & num_pose_search, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight,  core::Real const & extra_point_weight, numeric::xyzVector<core::Real> const & no_CoM_offset );

	core::Real fp_compare( FingerprintBase & fp, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight ) const;

	void fp_compare_deriv( FingerprintBase & fp, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, core::Real & dE_dx, core::Real & dE_dy, core::Real & dE_dz, core::Real & dE_dv4, core::Real & dE_dv5, core::Real & dE_dv6 ) const;

	void dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset );

	void dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform );

	void dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset,  numeric::xyzVector<core::Real> const & CoM_offset );

	void dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform, numeric::xyzVector<core::Real> const & CoM_offset );

	core::pose::Pose get_oriented_pose( FingerprintBase & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, numeric::xyzVector<core::Real> const & CoM_offset, core::Size const conformer );

	core::pose::Pose get_oriented_pose( FingerprintBase & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform, numeric::xyzVector<core::Real> const & CoM_offset, core::Size const conformer );

	core::Real rmsd( core::pose::Pose const & original_pose,  core::pose::Pose const & oriented_pose );

	// void move_origin(numeric::xyzVector<core::Real> const & new_origin );

	numeric::xyzVector<core::Real> calculate_ligand_CoM( core::pose::Pose const & ligand_pose );

	core::pose::Pose & pose() { return pose_; };
	core::Size compute_ligand_resnum( core::pose::Pose const & pose ) const;
	core::Size compute_ligand_resnum() const { return compute_ligand_resnum(pose_); };
	core::Size compute_ligand_natoms( core::pose::Pose const & pose ) const;
	core::Size compute_ligand_natoms() const { return compute_ligand_natoms(pose_); };
	core::Size compute_ligand_natoms_with_hydrogens( core::pose::Pose const & pose ) const;
	core::Size compute_ligand_natoms_with_hydrogens() const { return compute_ligand_natoms_with_hydrogens(pose_); };
	core::Size compute_ligand_nconformers( core::pose::Pose const & pose ) const;
	core::Size compute_ligand_nconformers() const { return compute_ligand_nconformers(pose_); };

	void move_ligand_and_update_rhos_(FingerprintBase & fp, numeric::xyzVector<core::Real> const & CoM_offset, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, core::Size const & conformer, bool const update_derivatives = false ) {
		core::conformation::ResidueCOP ligand_rsd = select_conf_and_move_ligand_( fp, CoM_offset, angle1_offset, angle2_offset, angle3_offset, conformer );
		update_rhos_( fp, ligand_rsd, update_derivatives );
	}

private:
	PlaidFingerprint(); // no default constructor

	core::pose::Pose pose_;

	// derivatives are optionally filled by build_from_pose_ , used in fp_compare_deriv
	std::list< ray_distance_derivs > derivs_of_ray_distances_; // note: refers to rays in the same order as the triplet data

	void update_rhos_( FingerprintBase & fp, core::conformation::ResidueCOP curr_ligand_rsd, bool const update_derivatives = false );

	core::conformation::ResidueCOP select_conf_and_move_ligand_( FingerprintBase & fp, numeric::xyzVector<core::Real> const & CoM_offset, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, core::Size const & conformer );

	void apply_rotation_offset_to_pose_( core::pose::Pose & pose, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset ) const;

};//class PlaidFingerprint

//void correct_phi_psi( core::Real & phi, core::Real & psi );

inline void convert_cartesian_to_spherical_coor_triplet( numeric::xyzVector<core::Real> const & coord, spherical_coor_triplet & triplet ){
	core::Real const triplet_rho = sqrt((coord.x()*coord.x())+(coord.y()*coord.y())+(coord.z()*coord.z()));
	triplet.rho = triplet_rho;
	triplet.phi = acos(coord.z()/triplet_rho);
	triplet.psi = atan2((coord.y()),(coord.x()));
}

inline void convert_spherical_coor_triplet_to_cartesian( spherical_coor_triplet const & triplet, numeric::xyzVector<core::Real> & coord ) {
	core::Real const triplet_phi = triplet.phi;
	core::Real const triplet_psi = triplet.psi;
	core::Real const triplet_rho = triplet.rho;
	core::Real const rho_times_sin_triplet_phi = triplet_rho*sin(triplet_phi);
	coord.x() = rho_times_sin_triplet_phi*cos(triplet_psi);
	coord.y() = rho_times_sin_triplet_phi*sin(triplet_psi);
	coord.z() = triplet_rho*cos(triplet_phi);
}

inline void convert_cartesian_to_grid( numeric::xyzVector<core::Real> const & cart_coord, numeric::xyzVector<core::Real> const & mid, numeric::xyzVector<core::Real> const & dim, core::Real const & spacing, numeric::xyzVector<core::Real> & grid_coord ) {
	grid_coord.x() = (cart_coord.x() - ( mid.x() - ((static_cast<core::Real>(dim.x()-1)/2) * spacing) ))/spacing;
	grid_coord.y() = (cart_coord.y() - ( mid.y() - ((static_cast<core::Real>(dim.y()-1)/2) * spacing) ))/spacing;
	grid_coord.z() = (cart_coord.z() - ( mid.z() - ((static_cast<core::Real>(dim.z()-1)/2) * spacing) ))/spacing;
}

inline void convert_grid_to_cartesian( numeric::xyzVector<core::Real> const & grid_coord, numeric::xyzVector<core::Real> const & mid, numeric::xyzVector<core::Real> const & dim, core::Real const & spacing, numeric::xyzVector<core::Real> & cart_coord ) {
	cart_coord.x() = grid_coord.x() * spacing + ( mid.x() - ((static_cast<core::Real>(dim.x()-1)/2) * spacing) );
	cart_coord.y() = grid_coord.y() * spacing + ( mid.y() - ((static_cast<core::Real>(dim.y()-1)/2) * spacing) );
	cart_coord.z() = grid_coord.z() * spacing + ( mid.z() - ((static_cast<core::Real>(dim.z()-1)/2) * spacing) );
}

// helper functions to compute derivatives
double dD_dv1(const double,const double,const double,const double,const double,const double,const double,const double,const double,const double) ;
double dD_dv2(const double,const double,const double,const double,const double,const double,const double,const double,const double,const double) ;
double dD_dv3(const double,const double,const double,const double,const double,const double,const double,const double,const double,const double) ;
double dD_dv4(const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double) ;
double dD_dv5(const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double) ;
double dD_dv6(const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double,const double) ;

// another helper function
core::Real Find_Closest_Intersect_SQ(core::Real const & phiAngle, core::Real const & psiAngle, core::Real const & atomX, core::Real const & atomY, core::Real const & atomZ, core::Real const & atom_radius );


}//pockets
}//protocols

#endif
