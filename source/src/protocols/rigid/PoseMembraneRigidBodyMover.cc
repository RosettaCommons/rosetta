// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PoseMembraneRigidBodyMover.cc
/// @brief
/// @author Yifan Song


#include <protocols/rigid/PoseMembraneRigidBodyMover.hh>

// Rosetta Headers
#include <core/pose/Pose.hh>

#include <core/types.hh>

#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/MembranePotential.hh>

#include <core/id/AtomID.hh>
#include <core/kinematics/Jump.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/pose/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>

//Auto Headers
//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


namespace protocols {
namespace rigid {

using namespace core;

static basic::Tracer TR( "protocols.moves.PoseMembraneRigidBodyMover" );

/// move pose into a membrane/*
MovePoseToMembraneCenterMover::MovePoseToMembraneCenterMover() : moves::Mover("MovePoseToMembraneCenterMover")
{
}

void MovePoseToMembraneCenterMover::apply( core::pose::Pose & pose )
{

	Vector const membrane_normal( scoring::MembraneEmbed_from_pose( pose ).normal() );
	Vector const membrane_center( scoring::MembraneEmbed_from_pose( pose ).center() );

	Vector const current_pose_center( estimate_membrane_center(pose) );

	Real const pose_shifted_from_membrane_center( membrane_normal.dot( current_pose_center-membrane_center ) );
	using namespace ObjexxFCL::format;
	//TR << "MovePoseToMembraneCenterMover: estimate membrane center " << F(8,3,current_pose_center.x()) << F(8,3,current_pose_center.y()) << F(8,3,current_pose_center.z()) << " shifted by:" << F(8,3,pose_shifted_from_membrane_center) << std::endl;
	if ( fabs(pose_shifted_from_membrane_center) > 20. ) {
		Vector shift_back = -pose_shifted_from_membrane_center * membrane_normal;

		WholeBodyTranslationMover whole_body_translation_mover(shift_back);
		whole_body_translation_mover.apply(pose);
		//TR << "MovePoseToMembraneCenterMover: shifting                 " << F(8,3,shift_back.x()) << F(8,3,shift_back.y()) << F(8,3,shift_back.z()) << std::endl;
	}

	//Vector const current_pose_center2( estimate_membrane_center(pose) );
	//Real const pose_shifted_from_membrane_center2( membrane_normal.dot( current_pose_center-membrane_center ) );
	//TR << "MovePoseToMembraneCenterMover: after moving center      " << F(8,3,current_pose_center2.x()) << F(8,3,current_pose_center2.y()) << F(8,3,current_pose_center2.z()) << " shifted by:" << F(8,3,pose_shifted_from_membrane_center2) << std::endl;
}

/// use membrane topology information to estimate the membrane center of the current pose position
Vector
MovePoseToMembraneCenterMover::estimate_membrane_center( core::pose::Pose & pose )
{
	core::scoring::MembraneTopology const & topology( core::scoring::MembraneTopology_from_pose(pose) );
	//Define vectors for inside and outside cap residue
	Vector sum_membrane_anchor(0);
	Vector center(0);

	for ( Size i=1; i<=topology.tmhelix(); ++i ) {
		//using namespace ObjexxFCL::format;
		//TR << "MovePoseToMembraneCenterMover: after moving center      " << I(4,i) << " " << topology.allow_tmh_scoring(i) << std::endl;

		if ( !topology.allow_tmh_scoring(i) ) continue;
		Vector const & start( pose.residue( topology.span_begin(i) ).atom( 2 ).xyz());
		Vector const & end( pose.residue( topology.span_end(i) ).atom( 2 ).xyz());

		sum_membrane_anchor+=start;
		sum_membrane_anchor+=end;
	}
	center = 0.5*(sum_membrane_anchor)/(float)topology.tmh_inserted();

	// if symmetrical pose, add subunit information
	if ( core::pose::symmetry::is_symmetric( pose ) ) {

		using namespace core::conformation::symmetry;
		auto const & symm_conf (
			dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

		for ( Size i_clone = 1; i_clone <= symm_info->num_bb_clones(); ++i_clone ) {

			for ( Size i=1; i<=topology.tmhelix(); ++i ) {
				if ( !topology.allow_tmh_scoring(i) ) continue;

				Vector const & start( pose.residue( symm_info->bb_clones( topology.span_begin(i) )[i_clone] ).atom( 2 ).xyz());
				Vector const & end( pose.residue( symm_info->bb_clones( topology.span_end(i) )[i_clone] ).atom( 2 ).xyz());

				sum_membrane_anchor+=start;
				sum_membrane_anchor+=end;
			}
		}
		center = 0.5*(sum_membrane_anchor)/(float)(symm_info->num_bb_clones() + 1)/(float)topology.tmh_inserted();
	}

	//using namespace ObjexxFCL::format;
	//TR << "Yifan debug: current origin          " << F(8,3,pose.conformation().atom_tree().root()->xyz().x()) << F(8,3,pose.conformation().atom_tree().root()->xyz().y()) << F(8,3,pose.conformation().atom_tree().root()->xyz().z()) << std::endl;
	//TR << "Yifan debug: current membrane center " << F(8,3,center.x()) << F(8,3,center.y()) << F(8,3,center.z()) << std::endl;
	//TR << "Yifan debug: current res 29 position " << F(8,3,pose.residue(29).atom("CA").xyz().x()) << F(8,3,pose.residue(29).atom("CA").xyz().y()) << F(8,3,pose.residue(29).atom("CA").xyz().z()) << std::endl;

	return center;
}

std::string
MovePoseToMembraneCenterMover::get_name() const {
	return "MovePoseToMembraneCenterMover";
}

/// perturb the pose along membrane normal
MembraneCenterPerturbationMover::MembraneCenterPerturbationMover() : moves::Mover("MembraneCenterPerturbationMover"),
	trans_mag_(10.)
{
}

MembraneCenterPerturbationMover::MembraneCenterPerturbationMover( core::Real const & trans_mag_in ) : moves::Mover() {
	trans_mag_ = trans_mag_in;
}

void MembraneCenterPerturbationMover::apply( core::pose::Pose & pose )
{
	Vector const membrane_normal( scoring::MembraneEmbed_from_pose( pose ).normal() );
	Real random_translation (trans_mag_ * (2.*numeric::random::rg().uniform()-1.));
	Vector const trans_vect = random_translation * membrane_normal;

	MovePoseToMembraneCenterMover move_pose_to_membrane_mover;
	move_pose_to_membrane_mover.apply(pose);

	WholeBodyTranslationMover whole_body_translation_mover(trans_vect);
	whole_body_translation_mover.apply(pose);
}

std::string MembraneCenterPerturbationMover::get_name() const {
	return "MembraneCenterPerturbationMover";
}

/// randomly rotate the pose along an axis on membrane plane
MembraneNormalPerturbationMover::MembraneNormalPerturbationMover() :
	Mover("MembraneNormalPerturbationMover"),
	rotation_mag_(15.)
{
}

MembraneNormalPerturbationMover::MembraneNormalPerturbationMover( core::Real const & rotation_mag_in ) : moves::Mover() {
	rotation_mag_ = rotation_mag_in;
}

void MembraneNormalPerturbationMover::apply( core::pose::Pose & pose )
{
	Vector const membrane_normal( scoring::MembraneEmbed_from_pose( pose ).normal() );
	Vector const membrane_center( scoring::MembraneEmbed_from_pose( pose ).center() );

	//find a random axis through center, perpendicular to normal
	Vector test_vec ;
	test_vec.x() = 1.; test_vec.y() = 1.; test_vec.z() = 1.;
	Vector x_axis (test_vec.cross(membrane_normal).normalize());
	Vector y_axis (membrane_normal.cross(x_axis).normalize());

	Real random_angle (numeric::conversions::radians(360.) * numeric::random::rg().uniform());
	Vector random_axis (cos(random_angle) * x_axis + sin(random_angle) * y_axis);
	Real random_rotation_angle (rotation_mag_ * numeric::random::rg().uniform());

	MovePoseToMembraneCenterMover move_pose_to_membrane_mover;
	move_pose_to_membrane_mover.apply(pose);

	WholeBodyRotationMover whole_body_rotation_mover(random_axis, membrane_center, random_rotation_angle);
	whole_body_rotation_mover.apply(pose);
}

std::string MembraneNormalPerturbationMover::get_name() const {
	return "MembraneNormalPerturbationMover";
}

/// Whole Body Translation
WholeBodyTranslationMover::WholeBodyTranslationMover( core::Vector const & trans_in ) : moves::Mover() {
	trans_ = trans_in;
}

void WholeBodyTranslationMover::apply( core::pose::Pose & pose )
{
	using namespace ObjexxFCL::format;

	// MovePoseToMembraneCenterMover move_pose_to_membrane_mover;

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		//utility::vector1 <kinematics::Jump> backup_jumps;
		//for (Size i =1; i<= pose.num_jump(); ++i) {
		// backup_jumps.push_back(pose.jump(i));
		//}

		core::kinematics::Jump first_jump = pose.jump( 1 );

		//TR << "Translate: " << "Jump (before): " << first_jump << std::endl;
		core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( 1 );
		first_jump.translation_along_axis( upstream_stub, trans_, trans_.length() );
		//TR << "Translate: " << "Jump (after):  " << first_jump << std::endl;
		pose.set_jump( 1, first_jump );

		/*
		id::AtomID root_id(pose.conformation().atom_tree().root()->atom_id());
		for (core::Size iatom=1; iatom <= pose.residue(root_id.rsd()).natoms(); ++iatom) {
		using namespace ObjexxFCL::format;
		id::AtomID atom_id(iatom, root_id.rsd());

		Vector const curr_position(pose.xyz(atom_id));
		Vector const new_position (curr_position + trans_);
		pose.set_xyz(atom_id, new_position);

		//TR << "Yifan debug: curr_position     " << I(4, iatom) << F(8,3,curr_position.x()) << F(8,3,curr_position.y()) << F(8,3,curr_position.z()) << std::endl;
		//TR << "Yifan debug: trans_            " << I(4, iatom) << F(8,3,trans_.x()) << F(8,3,trans_.y()) << F(8,3,trans_.z()) << std::endl;
		//TR << "Yifan debug: new_position      " << I(4, iatom) << F(8,3,new_position.x()) << F(8,3,new_position.y()) << F(8,3,new_position.z()) << std::endl;
		}

		for (Size i =1; i<= pose.num_jump(); ++i) {
		using namespace core::conformation::symmetry;
		SymmetricConformation const & symm_conf (
		dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
		if (!symm_info->jump_is_independent(i)) continue;

		pose.set_jump(i, backup_jumps[i]);
		}
		*/
	} else {
		//numeric::xyzMatrix< Real > const & R( numeric::xyzMatrix< Real >::identity() );
		//pose.apply_transform_Rx_plus_v(R,trans_);

		//core::kinematics::Stub upstream_stub root_stub = pose.conformation().atom_tree().root().get_stub();

		core::kinematics::Jump first_jump = pose.jump( 1 );

		//TR << "Translate: " << "Jump (before): " << first_jump << std::endl;
		core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( 1 );
		first_jump.translation_along_axis( upstream_stub, trans_, trans_.length() );
		//TR << "Translate: " << "Jump (after):  " << first_jump << std::endl;
		pose.set_jump( 1, first_jump );

	}


	//protocols::moves::AddPyMOLObserver(pose, true); // Lets ask PyMOL to store history...

	//TR << "Yifan debug: root origin    " << F(8,3,pose.conformation().atom_tree().root()->xyz().x()) << F(8,3,pose.conformation().atom_tree().root()->xyz().y()) << F(8,3,pose.conformation().atom_tree().root()->xyz().z()) << std::endl;
	//TR << "Yifan debug: new root position" << pose.conformation().atom_tree().root()->xyz() << std::endl;
	//test_center = move_pose_to_membrane_mover.estimate_membrane_center(pose);
	//TR << "Yifan debug: re-estimate membrane center " << F(8,3,test_center.x()) << F(8,3,test_center.y()) << F(8,3,test_center.z()) << std::endl;
	//pose.conformation().atom_tree().root()->show(1);
}

std::string
WholeBodyTranslationMover::get_name() const {
	return "WholeBodyTranslationMover";
}

/// Whole Body Rotation
WholeBodyRotationMover::WholeBodyRotationMover(
	Vector const & axis,
	Vector const & center,
	Real const & alpha /* degrees */
)
: moves::Mover()
{
	axis_   = axis;
	center_ = center;
	alpha_  = alpha;
}

/*
WholeBodyRotationMover::WholeBodyRotationMover( numeric::xyzMatrix< Real > const & R_in ) : Mover() {
rotation_m_ = R_in;
}
*/


void WholeBodyRotationMover::apply( core::pose::Pose & pose )
{
	debug_assert (! core::pose::symmetry::is_symmetric( pose ) );
	core::kinematics::Jump first_jump = pose.jump( 1 );

	//TR << "Rotate: " << "Jump (before): " << first_jump << std::endl;
	core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( 1 );
	first_jump.rotation_by_axis( upstream_stub, axis_, center_, alpha_ /* degrees */ );
	//TR << "Rotate: " << "Jump (after):  " << first_jump << std::endl;
	pose.set_jump( 1, first_jump );


	//Vector translation;
	//translation.zero();
	//pose.apply_transform_Rx_plus_v(rotation_m_, translation);
}

std::string
WholeBodyRotationMover::get_name() const {
	return "WholeBodyRotationMover";
}

} //ns rigid
} //ns protocols
