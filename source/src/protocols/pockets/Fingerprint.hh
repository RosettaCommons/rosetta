// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/Fingerprint.hh
/// @brief  protocols::pockets::Fingerprint header
/// @author Ragul Gowthaman

#ifndef INCLUDED_protocols_pockets_Fingerprint_hh
#define INCLUDED_protocols_pockets_Fingerprint_hh

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/pockets/Fingerprint.fwd.hh>
#include <protocols/pockets/PocketGrid.fwd.hh>
#include <protocols/pockets/FingerprintMultifunc.hh>
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

  typedef struct {
    // note: these are in radians
    core::Real phi;
    core::Real psi;
    core::Real rho;
  } spherical_coor_triplet;

  typedef struct {
    core::Real dDist_dv1;
    core::Real dDist_dv2;
    core::Real dDist_dv3;
    core::Real dDist_dv4;
    core::Real dDist_dv5;
    core::Real dDist_dv6;
  } ray_distance_derivs;

  class FingerprintBase : public utility::pointer::ReferenceCount {

    friend class FingerprintMultifunc;
    friend class DarcParticleSwarmMinimizer;

  public:

    FingerprintBase();

    ///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
    virtual ~FingerprintBase();

    void print_to_file(std::string const & output_filename) const;
    //void print_to_file(std::string const & output_filename, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset) const;

    void print_to_pdb(std::string const & output_pdbname) const;

    void print_to_pdb(std::string const & output_pdbname, numeric::xyzVector<core::Real> const & translation ) const;

    // const accessor functions
    numeric::xyzVector<core::Real> origin() const { return origin_; };

    numeric::xyzVector<core::Real> CoM() const { return CoM_; };

    std::list< spherical_coor_triplet > const & triplet_fingerprint_data() const { return triplet_fingerprint_data_; };

    // CHEAT!!
    void CHEAT_CoM( numeric::xyzVector<core::Real> const & inp_CoM ) { CoM_ = inp_CoM; };

  protected:
    numeric::xyzVector<core::Real> origin_;
    std::list< spherical_coor_triplet > triplet_fingerprint_data_;
    numeric::xyzVector<core::Real> CoM_;

  };

  class NonPlaidFingerprint : public FingerprintBase {
  public:
    NonPlaidFingerprint();
    ~NonPlaidFingerprint();

    std::list< numeric::xyzVector<core::Real> > egg_and_ext_list_;
    std::list< numeric::xyzVector<core::Real> > eggshell_list_;
    std::list< numeric::xyzVector<core::Real> > extshell_list_;

    void setup_from_PocketGrid( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid );

    void setup_from_PocketGrid( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, PocketGrid const & grid_for_extshell );

		void setup_from_PocketGrid_and_known_ligand( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, PocketGrid const & grid_for_extshell, core::pose::Pose const & known_ligand_pose, core::Real const & trim_dist );

		void setup_from_PocketGrid_using_bound_ligand( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, PocketGrid const & grid_for_extshell, core::pose::Pose const & known_ligand_pose );
    void setup_from_EggshellGrid();

    void write_eggshell_to_pdb_file( std::string const & output_eggshell_name ) const;

#ifdef USEOPENCL
    void gpu_setup( core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, int & num_particles, PlaidFingerprint & pf );
    int gpu_calculate_particle_scores(core::optimization::ParticleOPs & particles, std::vector<basic::gpu::float4> &atoms, std::vector<basic::gpu::float4> &atom_maxmin_phipsi);
    void gpu_setup_rays();
#endif

    void setup_from_eggshell_pdb_file( std::string const & input_filename);

    void trim_based_on_known_ligand( core::pose::Pose const & known_ligand_pose );

		void include_eggshell_points_based_on_known_ligand( core::pose::Pose const & known_ligand_pose, core::Real const & trim_dist);

    void setup_from_eggshell_triplet_file(std::string const & input_filename);

    void setup_from_PlaidFingerprint( PlaidFingerprint const & pfp );

    void set_origin ( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell );

    void set_origin_from_option_( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, Size const & set_origin_option );

    void set_origin_away_from_protein_center ( core::pose::Pose const & protein_pose );

		void set_origin_from_residue ( core::pose::Pose const & protein_pose );

    void set_origin_away_from_eggshell( std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, core::pose::Pose const & protein_pose );

    void set_origin_away_from_eggshell_plane( std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, core::pose::Pose const & protein_pose, Size const & set_origin_option );

    core::Real get_Rvalue (core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, Size const & set_origin_option);

    numeric::xyzVector<core::Real> calculate_protein_CoM( core::pose::Pose const & protein_pose);

    std::list< numeric::xyzVector<core::Real> > combine_xyz_lists (std::list< numeric::xyzVector<core::Real> > const & xyz_list_1 , std::list< numeric::xyzVector<core::Real> > const & xyz_list_2);

		std::list<spherical_coor_triplet> remove_duplicate_phi_psi(std::list<spherical_coor_triplet> const & rounded_triplet);

		std::list<spherical_coor_triplet> convert_cart_to_spherical_and_round(std::list< numeric::xyzVector<core::Real> > const & xyz_list);

		std::list<numeric::xyzVector<core::Real> > convert_spherical_list_to_cartesian_list(std::list<spherical_coor_triplet> const & unique_triplet);

		std::list<spherical_coor_triplet> set_rho_to_zero(std::list<spherical_coor_triplet> const & rounded_triplet);

#ifdef USEOPENCL
	private:
    basic::gpu::GPU gpu_;
    struct {
    	cl_mem rays;
    	cl_mem atoms;
    	cl_mem atom_maxmin_phipsi;
    	cl_mem ray_scores;
    	cl_mem particle_scores;
    	cl_mem weights;
    	unsigned int rays_size;
    	unsigned int atoms_size;
    	unsigned int atom_maxmin_phipsi_size;
    	unsigned int ray_scores_size;
    	unsigned int particle_scores_size;
	    int num_rays;
	    int num_atoms;
	    int num_particles;
    } gpu_memory_;
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

    //	void move_origin(numeric::xyzVector<core::Real> const & new_origin );

    numeric::xyzVector<core::Real> calculate_ligand_CoM( core::pose::Pose const & ligand_pose );

    core::pose::Pose & pose() { return pose_; };
    core::Size compute_ligand_resnum( core::pose::Pose const & pose ) const;
    core::Size compute_ligand_resnum() const { return compute_ligand_resnum(pose_); };
    core::Size compute_ligand_natoms( core::pose::Pose const & pose ) const;
    core::Size compute_ligand_natoms() const { return compute_ligand_natoms(pose_); };
    core::Size compute_ligand_nconformers( core::pose::Pose const & pose ) const;
    core::Size compute_ligand_nconformers() const { return compute_ligand_nconformers(pose_); };

  private:
    PlaidFingerprint(); // no default constructor

    core::pose::Pose pose_;

    // derivatives are optionally filled by build_from_pose_ , used in fp_compare_deriv
    std::list< ray_distance_derivs > derivs_of_ray_distances_; // note: refers to rays in the same order as the triplet data

    void move_ligand_and_update_rhos_(FingerprintBase & fp, numeric::xyzVector<core::Real> const & CoM_offset, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, core::Size const & conformer, bool const update_derivatives = false ) {
      core::conformation::ResidueCOP ligand_rsd = select_conf_and_move_ligand_( fp, CoM_offset, angle1_offset, angle2_offset, angle3_offset, conformer );
      update_rhos_( fp, ligand_rsd, update_derivatives );
    }

    void update_rhos_( FingerprintBase & fp, core::conformation::ResidueCOP curr_ligand_rsd, bool const update_derivatives = false );

    core::conformation::ResidueCOP select_conf_and_move_ligand_( FingerprintBase & fp, numeric::xyzVector<core::Real> const & CoM_offset, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, core::Size const & conformer );

    void apply_rotation_offset_to_pose_( core::pose::Pose & pose, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset ) const;

  };

//void correct_phi_psi( core::Real & phi, core::Real & psi );
inline void
convert_cartesian_to_spherical_coor_triplet( numeric::xyzVector<core::Real> const & coord,
		spherical_coor_triplet & triplet )
{
	core::Real const triplet_rho = sqrt((coord.x()*coord.x()) + (coord.y()*coord.y()) + (coord.z()*coord.z()));
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
