// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/ui_protocols/HelicalBundlePoseDrawOpenGLWidget.cpp
/// @brief  A derived class for drawing parametrically-generated helical bundles.  Derived from
/// ui::ui_core::SimplePoseDrawOpenGLWidget class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Associated headers:
#include <ui/ui_protocols/helical_bundle/HelicalBundlePoseDrawOpenGLWidget.h>

// Numeric headers:

// Core headers:
#include <core/types.hh>

//Qt headers:
#include <QMouseEvent>
#include <QDebug>

namespace ui {
namespace ui_protocols {
namespace helical_bundle {

HelicalBundlePoseDrawOpenGLWidget::HelicalBundlePoseDrawOpenGLWidget( QWidget* parent ) :
	ui::ui_core::pose_draw::SimplePoseDrawOpenGLWidget( parent ),
	extended_drag_mode_( HBPDOGLW_base_class_mode ),
	x_unit_vect_( numeric::xyzVector< GLfloat >(1,0,0) ),
	y_unit_vect_( numeric::xyzVector< GLfloat >(0,1,0) ),
	z_unit_vect_( numeric::xyzVector< GLfloat >(0,0,1) ),
	nonparametric_translation_vect_( numeric::xyzVector< core::Real >(0,0,0) ),
	nonparametric_COM_vect_( numeric::xyzVector< core::Real >(0,0,0) ),
	nonparametric_rotation_xform_( numeric::xyzTransform< core::Real >::identity() )
	//TODO -- initialize vars here
{}

/// @brief Set the extended drag mode.
/// @details Set this to HBPDOGLW_base_class_mode to use the base class drag mode.
void
HelicalBundlePoseDrawOpenGLWidget::set_extended_drag_mode(
	HBPDOGLW_extended_drag_mode const mode_in
) {
	debug_assert( mode_in > 0  && mode_in <= HBPDOGLW_end_of_list );
	extended_drag_mode_ = mode_in;
}

/// @brief Given the current orientation of the viewport, determine unit vectors in the viewport X (horizontal), Y (vertical), and Z (out of screen) directions.
void
HelicalBundlePoseDrawOpenGLWidget::determine_xyz_unit_vectors() {
	numeric::xyzTransform< GLfloat > rotationx( numeric::xyzTransform< GLfloat >::rot_deg( numeric::xyzVector<GLfloat>(1, 0, 0), -rotx(), numeric::xyzVector<GLfloat>(0,0,0) ) );
	numeric::xyzTransform< GLfloat > rotationz( numeric::xyzTransform< GLfloat >::rot_deg( numeric::xyzVector<GLfloat>(0, 0, 1), -rotz(), numeric::xyzVector<GLfloat>(0,0,0) ) );
	x_unit_vect_ = rotationz*rotationx*numeric::xyzVector< GLfloat >(1,0,0);
	y_unit_vect_ = rotationz*rotationx*numeric::xyzVector< GLfloat >(0,1,0);
	z_unit_vect_ = rotationz*rotationx*numeric::xyzVector< GLfloat >(0,0,1);

	/*qDebug() << "SCREEN ROT: " << rotx() << " " << roty() << " " << rotz();
	qDebug() << "X UNIT: " << x_unit_vect_.x() << " " << x_unit_vect_.y() << " " << x_unit_vect_.z();
	qDebug() << "Y UNIT: " << y_unit_vect_.x() << " " << y_unit_vect_.y() << " " << y_unit_vect_.z();
	qDebug() << "Z UNIT: " << z_unit_vect_.x() << " " << z_unit_vect_.y() << " " << z_unit_vect_.z();
	*/
}

/// @brief Reset the current nonparametric translation vector.
void
HelicalBundlePoseDrawOpenGLWidget::reset_nonparametric_translation_vect() {
	nonparametric_translation_vect_.x(0);
	nonparametric_translation_vect_.y(0);
	nonparametric_translation_vect_.z(0);
}

/// @brief Reset the current nonparametric rotation transform.
void
HelicalBundlePoseDrawOpenGLWidget::reset_nonparametric_rotation_xform() {
	nonparametric_rotation_xform_ = numeric::xyzTransform< core::Real >::identity();
}

/// @brief Set the current nonparametric centre of mass.
void
HelicalBundlePoseDrawOpenGLWidget::set_nonparametric_COM_vect(
	numeric::xyzVector< core::Real > const &vect_in
) {
	nonparametric_COM_vect_ = vect_in;
}



//////////////////////////PROTECTED/////////////////////////////////////////

/// @brief What to do when the user moves the mouse.
/// @details Overrides the base class mouseMoveEvent() function.  Adds functionality for the
/// extended drag modes, while calling the base class function for the HBPDOGLW_base_class_mode.
void
HelicalBundlePoseDrawOpenGLWidget::mouseMoveEvent(
	QMouseEvent* event
) {
	if( extended_drag_mode_ == HBPDOGLW_base_class_mode ) {
		ui::ui_core::pose_draw::SimplePoseDrawOpenGLWidget::mouseMoveEvent( event );
		return;
	} else if ( extended_drag_mode_ == HBPDOGLW_drag_nonparametric_geometry ) {
		if(event->buttons() == Qt::LeftButton ) {
			GLfloat const dx( 75 * static_cast<GLfloat>( event->x() - lastpos().x() ) / static_cast<GLfloat>( width() ) );
			GLfloat const dy( -75 * static_cast<GLfloat>( event->y() - lastpos().y() ) / static_cast<GLfloat>( height() ) );

			numeric::xyzVector< core::Real > const delta_x( static_cast< numeric::xyzVector < core::Real > >( x_unit_vect_ * dx ) );
			numeric::xyzVector< core::Real > const delta_y( static_cast< numeric::xyzVector < core::Real > >( y_unit_vect_ * dy ) );

			nonparametric_translation_vect_ += (delta_x + delta_y);

			Q_EMIT nonparametric_geometry_has_moved();

			set_lastpos( event->pos() );
		} else if(event->buttons() == Qt::RightButton) {
			GLfloat const dz( 75 * static_cast<GLfloat>( event->y() - lastpos().y() ) / static_cast<GLfloat>( height() ) );
			numeric::xyzVector< core::Real > const delta_z( static_cast< numeric::xyzVector< core::Real > >( z_unit_vect_ * dz ) );
			nonparametric_translation_vect_ += delta_z;
			Q_EMIT nonparametric_geometry_has_moved();
			set_lastpos( event->pos() );
		}
		return;
	} else if ( extended_drag_mode_ == HBPDOGLW_rotate_nonparametric_geometry ) { //If we're rotating the nonparametric geometry
		if( event->buttons() == Qt::LeftButton ) {
			GLfloat const dx( static_cast<GLfloat>( event->x() - lastpos().x() ) / static_cast<GLfloat>( width() ) );
			GLfloat const dy( static_cast<GLfloat>( event->y() - lastpos().y() ) / static_cast<GLfloat>( height() ) );

			nonparametric_rotation_xform_ =
					numeric::xyzTransform<core::Real>::rot_deg( y_unit_vect_, 180*dx, nonparametric_COM_vect_-nonparametric_translation_vect_ )*
					numeric::xyzTransform<core::Real>::rot_deg( x_unit_vect_, 180*dy, nonparametric_COM_vect_-nonparametric_translation_vect_ )*
					nonparametric_rotation_xform_;

			set_lastpos( event->pos() );

			Q_EMIT nonparametric_geometry_has_moved();
		} else if( event->buttons() == Qt::RightButton ) {
			GLfloat const dy( static_cast<GLfloat>( event->y() - lastpos().y() ) / static_cast<GLfloat>( height() ) );
			nonparametric_rotation_xform_ = numeric::xyzTransform<core::Real>::rot_deg( z_unit_vect_, -180*dy, nonparametric_COM_vect_-nonparametric_translation_vect_ )*nonparametric_rotation_xform_;

			set_lastpos( event->pos() );

			Q_EMIT nonparametric_geometry_has_moved();
		}
	}
}


} //namespace helical_bundle
} //namespace ui_protocols
} //namespace ui
