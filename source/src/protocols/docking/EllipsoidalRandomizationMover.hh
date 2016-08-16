// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/docking/EllipsoidalRandomizationMover.hh
/// @brief
/// @author Nick Marze (nickmarze@gmail.com)
#ifndef INCLUDED_protocols_docking_EllipsoidalRandomizationMover_hh
#define INCLUDED_protocols_docking_EllipsoidalRandomizationMover_hh

// Unit Headers
#include <protocols/docking/EllipsoidalRandomizationMover.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <protocols/moves/Mover.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <utility/vector1.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;
namespace protocols {
namespace docking {

/// @details This mover approximates one docking partner as an ellipsoid and docks the other partner
/// along a random normal to said ellipsoid; sliding the partners into contact is necessary to complete
/// this docking move.
class EllipsoidalRandomizationMover : public moves::Mover {
public:
	/// @brief Default constructor
	EllipsoidalRandomizationMover();

	/// @brief Copy constructor
	EllipsoidalRandomizationMover( EllipsoidalRandomizationMover const & object_to_copy );

	/// @brief Constructor with arguments
	EllipsoidalRandomizationMover( core::Size, bool );

	/// @brief Assignment operator
	EllipsoidalRandomizationMover & operator=( EllipsoidalRandomizationMover const & object_to_copy );

	/// @brief Destructor
	~EllipsoidalRandomizationMover();

	/// @brief Sets up the default values for the object including the movemap and minimization type.
	void set_default();
	/// @brief Initializes data members from options.
	void init_from_options();

	/// @brief Applies EllipsoidalRandomizationMover mover.
	void apply( core::pose::Pose & );

	/// @brief Returns mover name.
	virtual std::string get_name() const;

	/// @brief Generates string representation of EllipsoidalRandomizationMover.
	virtual void show( std::ostream & output=std::cout ) const;

	virtual protocols::moves::MoverOP clone() const;

	virtual protocols::moves::MoverOP fresh_instance() const;

	/// @brief Defines geometry required for docking move.
	std::pair< numeric::xyzMatrix< core::Real >, Vector > calculate_geometry( core::pose::Pose & );

	/// @brief Calculates ellipsoid axes from atomic coordinates.
	numeric::xyzMatrix< core::Real > calculate_axes( core::pose::Pose & );

	/// @brief Returns a sample from a beta distribution.
	Real single_beta_sample( double, double );

	/// @brief Selects a random point on ellipsoid surface.
	Vector point_on_ellipsoid( numeric::xyzMatrix< core::Real > );

	/// @brief Calculates 3rd coordinate of a point on ellipsoid surface given first 2 coordinates.
	Real recalculate_z_coordinate( Vector );

	/// @brief Calculates normal vector to a point on ellipsoid.
	Vector normal_to_ellipsoid( Vector );

	/// @brief Calculates a rotation matrix to superimpose source vector to target vector.
	numeric::xyzMatrix< core::Real > get_rotation_matrix( Vector, Vector );

	/// @brief Identifies first and last residues of a docking partner.
	utility::vector1< core::Size > get_partner_residue_start_stop( core::pose::Pose &, bool );

	/// @brief Sets pose fold tree based on -partners flag.
	void set_foldtree( core::pose::Pose & );

	/// @brief Identifies residues at interface of docking partners.
	utility::vector1< bool > get_interface_residues( core::pose::Pose &, core::Real );

	/// @brief Calculates two vectors defining a plane at the docking interface.
	utility::vector1< Vector > calculate_plane_axes( core::pose::Pose & );

	/// @brief Calculates normal vector from a plane at the docking interface.
	Vector inward_normal_to_plane( utility::vector1< Vector > );

	/// @brief Copies data members from one object to another.
	void copy_data( EllipsoidalRandomizationMover object_to_copy_to, EllipsoidalRandomizationMover object_to_copy_from );

	/// @brief Set the partners_ string.
	void set_partners( std::string const & );

	/// @brief Get the axis for sliding partners into contact.
	Vector get_slide_axis() const;

	/// @brief Get the point about which to spin the docking partner.
	Vector get_spin_center() const;

	//For unit tests
	Vector get_c_alpha_centroid() const;
	Vector get_c_alpha_non_ellipsoid_centroid() const;
	core::Real get_a_axis() const;
	core::Real get_b_axis() const;
	core::Real get_c_axis() const;
	Vector get_random_point_on_ellipsoid() const;
	Vector get_normal_to_plane() const;

private:

	//Centroid of all ellipsoid c-alpha coordinates
	Vector c_alpha_centroid_;
	//Centroid of non-ellipsoid interface residue c-alpha coordinates
	Vector c_alpha_plane_centroid_;
	//Centroid of all non-ellipsoid c-alpha coordinates
	Vector c_alpha_non_ellipsoid_centroid_;

	Vector slide_axis_;

	//Ellipsoid axes; c > a > b
	core::Real a_axis_;
	core::Real b_axis_;
	core::Real c_axis_;

	Vector random_point_on_ellipsoid_;
	Vector normal_to_plane_;

	core::Size rb_jump_;

	bool ellipsoid_is_first_partner_;

	std::string partners_;

};

std::ostream & operator<<( std::ostream & output, EllipsoidalRandomizationMover const & object_to_output );

}
}
#endif
