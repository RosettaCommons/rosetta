// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/ui_protocols/HelicalBundlePoseDrawOpenGLWidget.h
/// @brief  Headers for a derived class for drawing parametrically-generated helical bundles.  Derived from
/// ui::ui_core::SimplePoseDrawOpenGLWidget class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef HELICALBUNDLEPOSEDRAWOPENGLWIDGET_H
#define HELICALBUNDLEPOSEDRAWOPENGLWIDGET_H

//Qt headers:
#include <QObject>

// Associated headers:
#include <ui/ui_protocols/helical_bundle/HelicalBundleDialogueWidget.fwd.h>

// Rosetta numeric headers:
#include <numeric/xyzVector.hh>
#include <numeric/xyzTransform.hh>

// Rosetta core headers:
#include <core/types.hh>

// UI core headers:
#include <ui/ui_core/pose_draw/SimplePoseDrawOpenGLWidget.h>

namespace ui {
namespace ui_protocols {
namespace helical_bundle {

/// @brief The ways in which a user can manipulate the display: rotate everything, drag non-parametric geometry, rotate non-parametric
/// geometry, etc.
enum HBPDOGLW_extended_drag_mode {
	HBPDOGLW_base_class_mode = 1, //Keep this first
	HBPDOGLW_rotate_nonparametric_geometry,
	HBPDOGLW_drag_nonparametric_geometry, //Keep this second-to-last
	HBPDOGLW_end_of_list = HBPDOGLW_drag_nonparametric_geometry
};

class HelicalBundlePoseDrawOpenGLWidget : public ui::ui_core::pose_draw::SimplePoseDrawOpenGLWidget
{
	Q_OBJECT

public:

	/// @brief Constructor
	HelicalBundlePoseDrawOpenGLWidget( QWidget* parent = nullptr );

	/// @brief Set the extended drag mode.
	/// @details Set this to HBPDOGLW_base_class_mode to use the base class drag mode.
	void set_extended_drag_mode( HBPDOGLW_extended_drag_mode const mode_in );

	/// @brief Given the current orientation of the viewport, determine unit vectors in the viewport X (horizontal), Y (vertical), and Z (out of screen) directions.
	void determine_xyz_unit_vectors();

	/// @brief Get the current nonparametric translation vector.
	inline numeric::xyzVector< core::Real > const & get_nonparametric_translation_vect() const { return nonparametric_translation_vect_; }

	/// @brief Get the current nonparametric rotation transform.
	inline numeric::xyzTransform< core::Real > const & get_nonparametric_rotation_xform() const { return nonparametric_rotation_xform_; }

	/// @brief Reset the current nonparametric translation vector.
	void reset_nonparametric_translation_vect();

	/// @brief Reset the current nonparametric rotation transform.
	void reset_nonparametric_rotation_xform();

	/// @brief Set the current nonparametric centre of mass.
	void set_nonparametric_COM_vect( numeric::xyzVector< core::Real > const &vect_in );

Q_SIGNALS:

	/// @brief Signal that the nonparametric geometry has moved, and that it's time to update it.
	void nonparametric_geometry_has_moved();

protected:

	/// @brief What to do when the user moves the mouse.
	/// @details Overrides the base class mouseMoveEvent() function.  Adds functionality for the
	/// extended drag modes, while calling the base class function for the HBPDOGLW_base_class_mode.
	void mouseMoveEvent( QMouseEvent* event ) override;

private:

	/// @brief The drag mode, with extended options for translating the nonparametric geometry or rotating it.
	HBPDOGLW_extended_drag_mode extended_drag_mode_;

	/// @brief A unit vector in the X-direction, used for translating in the screen plane.
	numeric::xyzVector< GLfloat > x_unit_vect_;

	/// @brief A unit vector in the Y-direction, used for translating in the screen plane.
	numeric::xyzVector< GLfloat > y_unit_vect_;

	/// @brief A unit vector in the Z-direction, used for translating perpendicular to the screen plane.
	numeric::xyzVector< GLfloat > z_unit_vect_;

	/// @brief The translation vector for nonparametric geometry.
	numeric::xyzVector< core::Real > nonparametric_translation_vect_;

	/// @brief The centre-of-mass of nonparametric geometry, for rotation.
	numeric::xyzVector< core::Real > nonparametric_COM_vect_;

	/// @brief The rotation vector for nonparametric geometry.
	numeric::xyzTransform< core::Real > nonparametric_rotation_xform_;

};

} //namespace helical_bundle
} //namespace ui_protocols
} //namespace ui

#endif // HELICALBUNDLEPOSEDRAWOPENGLWIDGET_H
