// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/docking/EllipsoidalRandomizationMover.cc
/// @brief
/// @author Nick Marze (nickmarze@gmail.com)

// Unit Headers
#include <protocols/docking/EllipsoidalRandomizationMover.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/types.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <numeric/PCA.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

// Utility Headers
#include <utility/vector1.hh>

#include <numeric/random/DistributionSampler.hh>
#include <numeric/trig.functions.hh>

#include <protocols/docking/util.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>

#include <protocols/scoring/Interface.hh>

#include <basic/Tracer.hh>


using namespace protocols::moves;
using namespace core;

static basic::Tracer TR( "protocols.docking.EllipsoidalRandomizationMover" );

namespace protocols {
namespace docking {

EllipsoidalRandomizationMover::EllipsoidalRandomizationMover() : moves::Mover()
{
	set_default();
	init_from_options();
}

EllipsoidalRandomizationMover::EllipsoidalRandomizationMover( EllipsoidalRandomizationMover const & object_to_copy ) : moves::Mover()
{
	copy_data( *this, object_to_copy );
}

EllipsoidalRandomizationMover::EllipsoidalRandomizationMover( core::Size rb_jump, bool ellipsoid_is_first_partner, bool autofoldtree ) : moves::Mover()
{
	set_default();
	init_from_options();

	ellipsoid_is_first_partner_ = ellipsoid_is_first_partner;
	rb_jump_ = rb_jump;
	autofoldtree_ = autofoldtree;
}

EllipsoidalRandomizationMover &
EllipsoidalRandomizationMover::operator=( EllipsoidalRandomizationMover const & object_to_copy )
{
	if ( this == &object_to_copy ) {
		return *this;
	}

	Mover::operator=( object_to_copy );
	copy_data( *this, object_to_copy );
	return *this;
}

EllipsoidalRandomizationMover::~EllipsoidalRandomizationMover()= default;

void
EllipsoidalRandomizationMover::set_default()
{

	a_axis_ = 1.0;
	b_axis_ = 1.0;
	c_axis_ = 1.0;

	c_alpha_centroid_.zero();
	c_alpha_plane_centroid_.zero();
	c_alpha_non_ellipsoid_centroid_.zero();

	rb_jump_ = 1;
	slide_axis_.zero();

	random_point_on_ellipsoid_.zero();
	normal_to_plane_.zero();

	ellipsoid_is_first_partner_ = true;
	autofoldtree_ = true;
	partners_ = "_";
	return;
}

void
EllipsoidalRandomizationMover::init_from_options()
{
	using namespace basic::options;

	// This defaults to "_"
	if ( option[ OptionKeys::docking::partners ].user() ) {
		set_partners(option[ OptionKeys::docking::partners ]());
	}

	if ( option[ OptionKeys::docking::randomize1 ].user() ) {
		ellipsoid_is_first_partner_ = true;
	}

	if ( option[ OptionKeys::docking::randomize2 ].user() ) {
		ellipsoid_is_first_partner_ = false;
	}
}

void
EllipsoidalRandomizationMover::copy_data( EllipsoidalRandomizationMover object_to_copy_to, EllipsoidalRandomizationMover object_to_copy_from )
{
	object_to_copy_to.c_alpha_centroid_ = object_to_copy_from.c_alpha_centroid_;
	object_to_copy_to.c_alpha_plane_centroid_ = object_to_copy_from.c_alpha_plane_centroid_;
	object_to_copy_to.c_alpha_non_ellipsoid_centroid_ = object_to_copy_from.c_alpha_non_ellipsoid_centroid_;
	object_to_copy_to.slide_axis_ = object_to_copy_from.slide_axis_;
	object_to_copy_to.a_axis_ = object_to_copy_from.a_axis_;
	object_to_copy_to.b_axis_ = object_to_copy_from.b_axis_;
	object_to_copy_to.c_axis_ = object_to_copy_from.c_axis_;
	object_to_copy_to.rb_jump_ = object_to_copy_from.rb_jump_;
	object_to_copy_to.ellipsoid_is_first_partner_ = object_to_copy_from.ellipsoid_is_first_partner_;
	object_to_copy_to.partners_ = object_to_copy_from.partners_;
}

void
EllipsoidalRandomizationMover::apply( core::pose::Pose & pose )
{
	TR << "Calculating geometry for rotation and translation..." << std::endl;

	std::pair< numeric::xyzMatrix< core::Real >, Vector > rotation_and_translation = calculate_geometry( pose );

	Vector translation_vector = rotation_and_translation.second;
	numeric::xyzMatrix< core::Real > rotation_matrix = rotation_and_translation.first;

	Real translation_length = translation_vector.length();

	rigid::RigidBodyTransMover superimpose_normal_takeoffs;

	TR << "Performing translation" << std::endl;

	superimpose_normal_takeoffs.trans_axis( translation_vector );
	superimpose_normal_takeoffs.step_size( translation_length );
	superimpose_normal_takeoffs.apply( pose );

	TR << "Performing rotation" << std::endl;

	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
	core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( rb_jump_ );
	flexible_jump.rotation_by_matrix( upstream_stub, c_alpha_plane_centroid_, rotation_matrix );
	pose.set_jump( rb_jump_, flexible_jump );

	TR << "Sliding partners away" << std::endl;

	//move partners apart to remove clashes
	Real const slide_away_factor = 0.5;
	//use normal to point on ellipsoid surface as slide axis, rather than line between centroids,
	//to preserve inter-partner contacts
	rigid::RigidBodyTransMover slide_away( slide_axis_, rb_jump_ );
	slide_away.step_size( slide_away_factor * c_axis_ );
	slide_away.apply( pose );

}

std::pair< numeric::xyzMatrix< core::Real >, Vector >
EllipsoidalRandomizationMover::calculate_geometry( core::pose::Pose & pose )
{
	//How ellipsoidal randomization is achieved:
	//1. Define an ellipsoid based on PCA of one partner's c-alpha coordinates
	//   (3 principal components).
	//2. Define a plane based on PCA of other partner's c-alpha coordinates at
	//   initial protein-protein interface (2 principal components).
	//3. Select a random point on ellipsoid surface.
	//4. Calculate the (empirical) outward normal to the random point.
	//5. Calculate the normal to the non-ellipsoid plane, inward to the non-ellipsoid
	//   partner.
	//6. Calculate the translation and rotation necessary to align corresponding
	//   points on the two normal vectors.

	//ellipsoid partner
	numeric::xyzMatrix< Real > axes = calculate_axes( pose );
	//non-ellipsoid partner
	utility::vector1< Vector > plane_axes = calculate_plane_axes( pose );

	random_point_on_ellipsoid_ = point_on_ellipsoid( axes );

	Vector normal_to_point = normal_to_ellipsoid( random_point_on_ellipsoid_ );
	Vector normal_terminus_from_origin = normal_to_point + random_point_on_ellipsoid_;

	Vector unit_z_axis;
	unit_z_axis.x( 0.0 );
	unit_z_axis.y( 0.0 );
	unit_z_axis.z( 1.0 );

	numeric::xyzMatrix< Real > rotation_matrix = get_rotation_matrix( axes.col(1), unit_z_axis );

	Vector rotated_point_on_ellipsoid = rotation_matrix * random_point_on_ellipsoid_;
	Vector rotated_normal_terminus = rotation_matrix * normal_terminus_from_origin;

	//Reorients point to original ellipsoid
	rotated_point_on_ellipsoid += c_alpha_centroid_;
	rotated_normal_terminus += c_alpha_centroid_;

	normal_to_plane_ = inward_normal_to_plane( plane_axes );
	Vector normal_terminus_from_plane = normal_to_plane_ + c_alpha_plane_centroid_;

	Vector placeholder_vector = rotated_normal_terminus - rotated_point_on_ellipsoid;
	Vector long_normal = rotated_point_on_ellipsoid + 2 * placeholder_vector;

	//Putting normal vectors in pdb format for visualization
	TR << "ATOM      1  C   ALA C   1" << rotated_point_on_ellipsoid << std::endl;
	TR << "ATOM      2  C   ALA C   2" << long_normal << std::endl;

	//4 points to align:
	//rotated_point_on_ellipsoid -> c_alpha_plane_centroid_
	//rotated_normal_terminus -> normal_terminus_from_plane

	Vector translation_vector;

	if ( ellipsoid_is_first_partner_ ) {
		//translate plane partner to ellipsoid partner
		translation_vector = rotated_point_on_ellipsoid - c_alpha_plane_centroid_;
		slide_axis_ = rotated_normal_terminus - rotated_point_on_ellipsoid;
		normal_terminus_from_plane += translation_vector;
		c_alpha_plane_centroid_ += translation_vector;
	} else {
		//translate ellipsoid partner to plane partner
		translation_vector = c_alpha_plane_centroid_ - rotated_point_on_ellipsoid;
		slide_axis_ = c_alpha_plane_centroid_ - normal_terminus_from_plane;
		rotated_normal_terminus += translation_vector;
		rotated_point_on_ellipsoid += translation_vector;
	}

	Vector rotation_source_vector;
	Vector rotation_target_vector;

	if ( ellipsoid_is_first_partner_ ) {
		//rotate plane partner onto ellipsoid normal
		rotation_source_vector = normal_terminus_from_plane - c_alpha_plane_centroid_;
		rotation_target_vector = rotated_normal_terminus - c_alpha_plane_centroid_;
	} else {
		//rotate ellipsoid partner onto plane normal
		rotation_source_vector = rotated_normal_terminus - c_alpha_plane_centroid_;
		rotation_target_vector = normal_terminus_from_plane - c_alpha_plane_centroid_;
	}

	numeric::xyzMatrix< Real > rotation_mover_matrix = get_rotation_matrix( rotation_target_vector, rotation_source_vector );

	std::pair<numeric::xyzMatrix< Real >, Vector > rotation_and_translation;
	rotation_and_translation = ( std::make_pair( rotation_mover_matrix, translation_vector ) );

	return rotation_and_translation;
}

std::string
EllipsoidalRandomizationMover::get_name() const
{
	return "EllipsoidalRandomizationMover";
}

void
EllipsoidalRandomizationMover::set_partners( std::string const & partners )
{
	partners_ = partners;
	return;
}

numeric::xyzMatrix< Real >
EllipsoidalRandomizationMover::calculate_axes( core::pose::Pose & pose_in )
{

	utility::vector1< Size > partner_residue_start_stop = get_partner_residue_start_stop( pose_in, ellipsoid_is_first_partner_ );

	utility::vector1< Vector > c_alpha_coords;
	for ( Size i=partner_residue_start_stop[1]; i<=partner_residue_start_stop[2]; ++i ) {
		c_alpha_coords.push_back( pose_in.residue( i ).xyz( "CA" ) );
	}

	Size n_res = c_alpha_coords.size();

	for ( Size i = 1; i <= n_res; ++i ) {
		c_alpha_centroid_ += c_alpha_coords[i];
	}
	c_alpha_centroid_ /= n_res;

	numeric::xyzMatrix< Real > eigenvectors = numeric::principal_components( c_alpha_coords );

	Vector eigenvalues = numeric::principal_component_eigenvalues( c_alpha_coords );

	//2.0 works well for roughly ellipsoid proteins
	//Should be increased for highly non-ellipsoidal proteins
	Real axis_mult_factor = 2.0;

	numeric::xyzMatrix< Real > ellipsoid_axes;
	ellipsoid_axes.col_x( axis_mult_factor * eigenvectors.col_x() * std::sqrt( std::fabs( eigenvalues.x() ) ) );
	ellipsoid_axes.col_y( axis_mult_factor * eigenvectors.col_y() * std::sqrt( std::fabs( eigenvalues.y() ) ) );
	ellipsoid_axes.col_z( axis_mult_factor * eigenvectors.col_z() * std::sqrt( std::fabs( eigenvalues.z() ) ) );

	return ellipsoid_axes;
}

utility::vector1< core::Size >
EllipsoidalRandomizationMover::get_partner_residue_start_stop( core::pose::Pose & pose_in, bool first_partner )
{
	utility::vector1< core::Size > start_stop;

	pose::PDBInfoCOP pdb_info = pose_in.pdb_info();

	if ( !pdb_info ) {
		utility_exit_with_message("Attempting to identify docking partners, however, the pdb_info object associated "
			"with the pose does not exist.");
	}

	Size first_residue_second_partner = 0;
	Size last_residue_first_partner = 0;

	if ( partners_ == "_" ) {

		utility_exit_with_message("Attempting to identify docking partners, however, no partners are defined. Please use the -partners option flag");

	} else {
		char first_chain_second_partner = char();
		for ( Size i=1; i<=partners_.length()-1; ++i ) {
			if ( partners_[i-1] == '_' ) {
				first_chain_second_partner = partners_[i];
			}
		}
		for ( Size i=2; i<= pose_in.size(); ++i ) {
			if ( pdb_info->chain( i ) == first_chain_second_partner ) {
				first_residue_second_partner = i;
				last_residue_first_partner = i-1;
				break;
			}
		}
	}

	if ( first_partner ) {
		start_stop.push_back( 1 );
		start_stop.push_back( last_residue_first_partner );
	} else {
		start_stop.push_back( first_residue_second_partner );
		start_stop.push_back( pose_in.size() );
	}

	return start_stop;
}

Real
EllipsoidalRandomizationMover::single_beta_sample( double alpha_param, double beta_param )
{
#ifndef  __native_client__
	using Mybeta = boost::math::beta_distribution<double>;

	using core::Real;
	using numeric::random::DistributionSampler;

	Mybeta ellipsoid_dist( alpha_param, beta_param );
	DistributionSampler< Mybeta > sampler( ellipsoid_dist );

	Real beta_sample = sampler.sample();

	return beta_sample;
#endif
}


Vector
EllipsoidalRandomizationMover::point_on_ellipsoid( numeric::xyzMatrix< core::Real > axes )
{
	Vector first_pca_vector = axes.col(1);
	Vector second_pca_vector = axes.col(2);
	Vector third_pca_vector = axes.col(3);

	//Strips vectors of directionality; defines three axes of ellipsoid centered at origin
	c_axis_ = first_pca_vector.length();
	a_axis_ = second_pca_vector.length();
	b_axis_ = third_pca_vector.length();

	//Sample elliptical shell from ellipsoid, weighted by surface area of shell
	Real const elliptical_beta_dist_param = 1.5;
	Real const uniform_beta_dist_param = 1.0;

	Real ellipsoid_sample_z = single_beta_sample( elliptical_beta_dist_param, elliptical_beta_dist_param );
	Real z_coord = ( ( ellipsoid_sample_z - 0.5 ) * 2 * c_axis_ );
	Real ellipse_a = ( a_axis_ * std::sqrt( 1.0 - ( ( z_coord * z_coord ) / ( c_axis_ * c_axis_ ) ) ) );
	Real ellipse_b = ( b_axis_ * std::sqrt( 1.0 - ( ( z_coord * z_coord ) / ( c_axis_ * c_axis_ ) ) ) );

	Real ab_ratio = ( ellipse_a / ellipse_b );

	//Empirical correlation for beta distribution parameters that evenly sample ellipse
	Real const linear_regime_boundary = 16.0;
	Real const correlation_param_one = 0.18;
	Real const correlation_param_two = 0.5;
	Real beta_parameter;

	if ( ab_ratio <= linear_regime_boundary ) {
		beta_parameter = ( ( correlation_param_one * std::log( ab_ratio ) ) + correlation_param_two );
	} else {
		beta_parameter = 1.0;
	}

	Real ellipse_sample_x = single_beta_sample( beta_parameter, beta_parameter );
	Real x_coord = ( ( ellipse_sample_x - 0.5 ) * 2 * ellipse_a );

	Real pos_or_neg = single_beta_sample( uniform_beta_dist_param, uniform_beta_dist_param );

	//Account for both halves of ellipse
	Real y_coord;
	if ( pos_or_neg <= 0.5 ) {
		y_coord = std::sqrt( ( ellipse_b * ellipse_b * ( 1.0 - ( ( x_coord * x_coord ) / ( ellipse_a * ellipse_a ) ) ) ) );
	} else {
		y_coord = ( -1.0 ) * std::sqrt( ( ellipse_b * ellipse_b * ( 1.0 - ( ( x_coord * x_coord ) / ( ellipse_a * ellipse_a ) ) ) ) );
	}


	Vector random_point_on_surface;
	random_point_on_surface.x( x_coord );
	random_point_on_surface.y( y_coord );
	random_point_on_surface.z( z_coord );

	return random_point_on_surface;
}

Vector
EllipsoidalRandomizationMover::normal_to_ellipsoid( Vector ellipsoid_point )
{
	Vector point_plus_dx;
	Vector point_minus_dx;
	Vector point_plus_dy;
	Vector point_minus_dy;


	//Perturb random point in 4 directions: +/- x, +/- y
	Real const delta = 0.00001;

	point_plus_dx.x( ellipsoid_point.x() + delta );
	point_plus_dx.y( ellipsoid_point.y() );
	point_plus_dx.z( ellipsoid_point.z() );

	point_minus_dx.x( ellipsoid_point.x() - delta );
	point_minus_dx.y( ellipsoid_point.y() );
	point_minus_dx.z( ellipsoid_point.z() );

	point_plus_dy.x( ellipsoid_point.x() );
	point_plus_dy.y( ellipsoid_point.y() + delta );
	point_plus_dy.z( ellipsoid_point.z() );

	point_minus_dy.x( ellipsoid_point.x() );
	point_minus_dy.y( ellipsoid_point.y() - delta );
	point_minus_dy.z( ellipsoid_point.z() );


	//Adjust z coordinate to keep perturbed points on ellipsoid surface
	point_plus_dx.z( recalculate_z_coordinate( point_plus_dx ) );
	point_minus_dx.z( recalculate_z_coordinate( point_minus_dx ) );
	point_plus_dy.z( recalculate_z_coordinate( point_plus_dy ) );
	point_minus_dy.z( recalculate_z_coordinate( point_minus_dy ) );

	Vector dx_vector = point_plus_dx - point_minus_dx;
	Vector dy_vector = point_plus_dy - point_minus_dy;

	Vector cross_forward = dx_vector.cross( dy_vector );
	Vector cross_reverse = dy_vector.cross( dx_vector );

	Vector normal_terminus_forward = ellipsoid_point + cross_forward.normalize_any();
	Vector normal_terminus_reverse = ellipsoid_point + cross_reverse.normalize_any();

	if ( normal_terminus_reverse.length() < normal_terminus_forward.length() ) {
		return cross_forward.normalize_any();
	} else {
		return cross_reverse.normalize_any();
	}

}

Vector
EllipsoidalRandomizationMover::inward_normal_to_plane( utility::vector1< Vector > two_plane_vectors )
{
	Vector normal_to_plane;

	Vector cross_forward = two_plane_vectors[1].cross( two_plane_vectors[2] );
	Vector cross_reverse = two_plane_vectors[2].cross( two_plane_vectors[1] );

	Vector distance_to_centroid_forward = c_alpha_plane_centroid_ + cross_forward.normalize_any() - c_alpha_non_ellipsoid_centroid_;
	Vector distance_to_centroid_reverse = c_alpha_plane_centroid_ + cross_reverse.normalize_any() - c_alpha_non_ellipsoid_centroid_;

	if ( distance_to_centroid_reverse.length() < distance_to_centroid_forward.length() ) {
		return cross_reverse.normalize_any();
	} else {
		return cross_forward.normalize_any();
	}
}

Real
EllipsoidalRandomizationMover::recalculate_z_coordinate( Vector perturbed_ellipsoid_point )
{
	//z^2 = c^2(1-x^2/a^2-y^2/b^2)
	Real new_z_coord = 0.0;

	Real inside_root = ( ( c_axis_ * c_axis_ ) * ( 1.0 - ( ( perturbed_ellipsoid_point.x() * perturbed_ellipsoid_point.x() ) / ( a_axis_ * a_axis_ ) ) - ( ( perturbed_ellipsoid_point.y() * perturbed_ellipsoid_point.y() ) / ( b_axis_ * b_axis_ ) ) ) );
	if ( inside_root < 0.0 ) {
		new_z_coord = ( -1.0 ) * perturbed_ellipsoid_point.z() ;
	} else if ( perturbed_ellipsoid_point.z() < 0.0 ) {
		new_z_coord = ( -1.0 ) * std::sqrt( inside_root );
	} else /* */ {
		new_z_coord = std::sqrt( inside_root );
	}

	return new_z_coord;
}

numeric::xyzMatrix< core::Real >
EllipsoidalRandomizationMover::get_rotation_matrix( Vector target_vector, Vector source_vector )
{
	Vector scaled_target_vector;
	Vector scaled_source_vector;

	scaled_target_vector.x( target_vector.x() );
	scaled_target_vector.y( target_vector.y() );
	scaled_target_vector.z( target_vector.z() );
	scaled_target_vector.normalize_any();

	scaled_source_vector.x( source_vector.x() );
	scaled_source_vector.y( source_vector.y() );
	scaled_source_vector.z( source_vector.z() );
	scaled_source_vector.normalize_any();

	Real rotation_angle = numeric::arccos( scaled_source_vector.dot( scaled_target_vector ) );
	Vector rotation_axis = scaled_source_vector.cross( scaled_target_vector );

	numeric::xyzMatrix< Real > rotation_matrix = numeric::rotation_matrix( rotation_axis, rotation_angle );

	return rotation_matrix;
}

void
EllipsoidalRandomizationMover::set_foldtree( core::pose::Pose & pose_in )
{
	//just instantiating to match setup_foldtree arguments;
	//setup_foldtree clears this variable and defaults it to 1,
	//which is the correct number of jumps in this case
	DockJumps movable_jumps;

	docking::setup_foldtree( pose_in, partners_, movable_jumps );
	return;
}

utility::vector1< bool >
EllipsoidalRandomizationMover::get_interface_residues( core::pose::Pose & pose_in, core::Real interface_distance_cutoff, bool autofoldtree )
{
	using namespace core::scoring;
	using ObjexxFCL::FArray1D_bool;

	if ( autofoldtree ) {
		set_foldtree( pose_in );
		FArray1D_bool partner1_( pose_in.size(), false );
		pose_in.fold_tree().partition_by_jump( rb_jump_, partner1_);
	}

	protocols::scoring::Interface interface_obj( rb_jump_ );
	pose_in.update_residue_neighbors();
	interface_obj.distance( interface_distance_cutoff );
	interface_obj.calculate( pose_in );

	utility::vector1< bool > interface_residues;

	for ( core::Size i = 1; i <= pose_in.size(); ++i ) {
		if ( interface_obj.is_interface( i ) ) {
			interface_residues.push_back( true );
		} else {
			interface_residues.push_back( false );
		}
	}

	return interface_residues;
}

utility::vector1< Vector >
EllipsoidalRandomizationMover::calculate_plane_axes( core::pose::Pose & pose_in )
{
	bool non_ellipsoid_is_first_partner;
	if ( ellipsoid_is_first_partner_ ) {
		non_ellipsoid_is_first_partner = false;
	} else {
		non_ellipsoid_is_first_partner = true;
	}

	utility::vector1< Size > partner_residue_start_stop = get_partner_residue_start_stop( pose_in, non_ellipsoid_is_first_partner );

	utility::vector1< Vector > c_alpha_plane_coords;
	utility::vector1< Vector > c_alpha_non_ellipsoid_coords;

	core::Real interface_distance_cutoff = 8.0;
	Size n_res_plane = 0;

	//A several-residue interface is necessary to maintain good contacts
	while ( n_res_plane < 5 ) {
		if ( interface_distance_cutoff > 20.0 ) {
			utility_exit_with_message("Docking partners are not within 20.0 Angstroms. Please move partners closer together in input file");
		}
		utility::vector1< bool > is_interface = get_interface_residues( pose_in, interface_distance_cutoff, autofoldtree_ );
		TR << "Getting interface residues at " << interface_distance_cutoff << " Angstroms" << std::endl;
		for ( Size i=partner_residue_start_stop[1]; i<=partner_residue_start_stop[2]; ++i ) {
			// check for CAs first -- RNA for example would not have this
			if ( pose_in.residue_type(i).has("CA") ) {
				c_alpha_non_ellipsoid_coords.push_back( pose_in.residue( i ).xyz( "CA" ) );
				if ( is_interface[i] ) {
					c_alpha_plane_coords.push_back( pose_in.residue( i ).xyz( "CA" ) );
				}
			}
		}
		n_res_plane = c_alpha_plane_coords.size();
		interface_distance_cutoff += 1.0;
	}

	Size n_res_non_ellipsoid = c_alpha_non_ellipsoid_coords.size();

	for ( Size i = 1; i <= n_res_plane; ++i ) {
		c_alpha_plane_centroid_ += c_alpha_plane_coords[i];
	}
	c_alpha_plane_centroid_ /= n_res_plane;

	for ( Size i = 1; i <= n_res_non_ellipsoid; ++i ) {
		c_alpha_non_ellipsoid_centroid_ += c_alpha_non_ellipsoid_coords[i];
	}
	c_alpha_non_ellipsoid_centroid_ /= n_res_non_ellipsoid;

	numeric::xyzMatrix< Real > eigenvectors = numeric::principal_components( c_alpha_plane_coords );

	Vector eigenvalues = numeric::principal_component_eigenvalues( c_alpha_plane_coords );

	utility::vector1< Vector > plane_axes;
	plane_axes.push_back( eigenvectors.col_x() * std::sqrt( std::fabs( eigenvalues.x() ) ) );
	plane_axes.push_back( eigenvectors.col_y() * std::sqrt( std::fabs( eigenvalues.y() ) ) );

	return plane_axes;
}

Vector
EllipsoidalRandomizationMover::get_slide_axis() const
{
	return slide_axis_;
}

Vector
EllipsoidalRandomizationMover::get_spin_center() const
{
	return c_alpha_plane_centroid_;
}

Vector
EllipsoidalRandomizationMover::get_c_alpha_centroid() const
{
	return c_alpha_centroid_;
}

Vector
EllipsoidalRandomizationMover::get_c_alpha_non_ellipsoid_centroid() const
{
	return c_alpha_non_ellipsoid_centroid_;
}

core::Real
EllipsoidalRandomizationMover::get_a_axis() const
{
	return a_axis_;
}

core::Real
EllipsoidalRandomizationMover::get_b_axis() const
{
	return b_axis_;
}

core::Real
EllipsoidalRandomizationMover::get_c_axis() const
{
	return c_axis_;
}

Vector
EllipsoidalRandomizationMover::get_random_point_on_ellipsoid() const
{
	return random_point_on_ellipsoid_;
}

Vector
EllipsoidalRandomizationMover::get_normal_to_plane() const
{
	return normal_to_plane_;
}

protocols::moves::MoverOP
EllipsoidalRandomizationMover::clone() const
{
	return protocols::moves::MoverOP( new EllipsoidalRandomizationMover( *this ) );
}

protocols::moves::MoverOP
EllipsoidalRandomizationMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new EllipsoidalRandomizationMover() );
}

void
EllipsoidalRandomizationMover::show( std::ostream & output ) const
{
	using namespace std;

	Mover::show( output );

	output << "Current slide axis" << endl << slide_axis_ << endl;
	output << "Current rb jump" << endl << rb_jump_ << endl;

}

std::ostream &
operator<<( std::ostream & output, EllipsoidalRandomizationMover const & object_to_output )
{
	object_to_output.show( output );
	return output;
}

}
}
