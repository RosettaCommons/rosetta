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
// AUTO-REMOVED #include <protocols/pockets/FingerprintMultifunc.fwd.hh>
#include <protocols/pockets/PocketGrid.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1_bool.hh>
#include <list>
#include <cmath>

#include <utility/vector1.hh>




namespace protocols {
namespace pockets {

typedef struct {
	core::Real phi;
	core::Real psi;
	core::Real rho;
} spherical_coor_triplet;

class FingerprintBase : public utility::pointer::ReferenceCount {

	friend class FingerprintMultifunc;

public:

  FingerprintBase();

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
  NonPlaidFingerprint() {};

  void setup_from_PocketGrid( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid );

  void setup_from_EggshellGrid( core::pose::Pose const & protein_pose, EggshellGrid const & pocket_grid );

  void trim_based_on_known_ligand( core::pose::Pose const & known_ligand_pose );

  void setup_from_file(std::string const & input_filename);

  void setup_from_PlaidFingerprint( PlaidFingerprint const & pfp );



};

class PlaidFingerprint : public FingerprintBase {

	friend class FingerprintMultifunc;

public:


  PlaidFingerprint( core::pose::Pose const & input_pose, FingerprintBase const & fp );

	core::Real Find_Intersect(core::Real const & phiAngle, core::Real const & psiAngle, core::Real const & atomX, core::Real const & atomY, core::Real const & atomZ, core::Real const & atom_radius, core::Real const & desired_rho );

	core::Real find_optimal_rotation( FingerprintBase const & fp, core::Real const & angle_increment, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight );

	core::Real find_optimal_rotation( FingerprintBase const & fp, core::Real const & angle_increment, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, numeric::xyzVector<core::Real> const & no_CoM_offset );

	core::Real search_random_poses( FingerprintBase const & fp, core::Size const & num_pose_search, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight );

	core::Real search_random_poses( FingerprintBase const & fp, core::Size const & num_pose_search, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight,  core::Real const & extra_point_weight, numeric::xyzVector<core::Real> const & no_CoM_offset );


	core::Real fp_compare( FingerprintBase const & fp, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight ) const;
	void dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase const & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset );

	void dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase const & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform );

	void dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase const & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform, numeric::xyzVector<core::Real> const & CoM_offset );

	core::pose::Pose get_oriented_pose( FingerprintBase const & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform, numeric::xyzVector<core::Real> const & CoM_offset );

	core::Real rmsd( core::pose::Pose const & original_pose,  core::pose::Pose const & oriented_pose );

	void move_origin(numeric::xyzVector<core::Real> const & new_origin );


private:
	PlaidFingerprint(); // no default constructor

	core::pose::Pose pose_;

	void build_from_pose_( FingerprintBase const & fp);

	void build_from_pose_( FingerprintBase const & fp, numeric::xyzVector<core::Real> const & CoM_offset, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset);

	void apply_rotation_offset_to_pose_( core::pose::Pose & pose, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset ) const;

};

void correct_phi_psi( core::Real & phi, core::Real & psi );

inline void convert_cartesian_to_spherical_coor_triplet( numeric::xyzVector<core::Real> const & coord, spherical_coor_triplet & triplet ){
	triplet.rho = sqrt((coord.x()*coord.x())+(coord.y()*coord.y())+(coord.z()*coord.z()));
	triplet.phi = acos(coord.z()/triplet.rho)*(1/numeric::constants::f::pi_over_180) ;
	triplet.psi = atan2((coord.y()),(coord.x()))*(1/numeric::constants::f::pi_over_180);
}

inline void convert_spherical_coor_triplet_to_cartesian( spherical_coor_triplet const & triplet, numeric::xyzVector<core::Real> & coord ) {
	coord.x() = triplet.rho*sin(triplet.phi*numeric::constants::f::pi_over_180)*cos(triplet.psi*numeric::constants::f::pi_over_180);
	coord.y() = triplet.rho*sin(triplet.phi*numeric::constants::f::pi_over_180)*sin(triplet.psi*numeric::constants::f::pi_over_180);
	coord.z() = triplet.rho*cos(triplet.phi*numeric::constants::f::pi_over_180);
}


}//pockets
}//protocols


#endif
