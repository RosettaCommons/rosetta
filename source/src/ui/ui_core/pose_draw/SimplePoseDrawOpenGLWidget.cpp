// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/ui_core/SimplePoseDrawOpenGLWidget.cpp
/// @brief  A base class for drawing a simple Van der Waals representation of a pose, with
/// the option to colour atoms based on a residue selection or on per-residue score.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <ui/ui_core/pose_draw/SimplePoseDrawOpenGLWidget.h>

//Rosetta numeric includes
#include <numeric/xyzVector.hh>

// Rosetta core includes
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#if defined(MAC) || defined(__APPLE__) || defined (__OSX__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <QDebug>
#include <QMouseEvent>

namespace ui {
namespace ui_core {
namespace pose_draw {

///////////////////PUBLIC METHODS///////////////////////////////////

/// @brief Constructor.
SimplePoseDrawOpenGLWidget::SimplePoseDrawOpenGLWidget(QWidget* parent) :
	QOpenGLWidget(parent),
	pose_(nullptr),
	residue_selector_(nullptr),
	rotx_(90.0),
	roty_(0.0),
	rotz_(0.0),
	forward_back_(-40.0),
	lastpos_( 0, 0 ),
	colour_mode_( SPDOGLW_default_colours ),
	energy_range_max_(5.0),
	energy_range_min_(-3.0),
	score_type_( core::scoring::fa_rep ),
	drag_mode_( SPDOGLW_rotate_and_zoom_viewport )
{
	QSurfaceFormat fmt;
	fmt.setDepthBufferSize(24);
	fmt.setStencilBufferSize(8);
	fmt.setVersion(2,1);
	fmt.setProfile(QSurfaceFormat::CoreProfile);
	fmt.setSamples(4);
	setFormat( fmt );
}

/// @brief Set the pose that we'll be drawing.
void SimplePoseDrawOpenGLWidget::set_pose( core::pose::PoseCOP pose ) {
	pose_ = pose;
}

/// @brief Set the residue selector that we'll use for colouring the pose.
/// @details Only used if colour_mode_ == SPODGLW_colour_by_selection.
void
SimplePoseDrawOpenGLWidget::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	residue_selector_ = selector_in;
}

/// @brief Set the colour mode that we'll use for colouring the pose.
void
SimplePoseDrawOpenGLWidget::set_colour_mode(
	SPDOGLW_colour_mode const mode_in
) {
	runtime_assert_string_msg( mode_in >= SPDOGLW_default_colours && mode_in <= SPDOGLW_end_of_list, "Invalid colour mode selected!" );
	colour_mode_ = mode_in;
}

/// @brief Set the lower and upper ends of the energy range in the gradient of colour values used.
/// @details There is no check that upper > lower.
void
SimplePoseDrawOpenGLWidget::set_energy_range(
	core::Real const &lower,
	core::Real const &upper
) {
	energy_range_max_ = upper;
	energy_range_min_ = lower;
}

/// @brief Set the score term to use if we're colouring by a single score type.
/// @details Only used if colour_mode_ == SPDOGLW_colour_by_score_term.
void
SimplePoseDrawOpenGLWidget::set_score_type(
	core::scoring::ScoreType const scoretype_in
) {
	score_type_ = scoretype_in;
}

/// @brief The pose has changed; redraw it.
void
SimplePoseDrawOpenGLWidget::update_pose_draw() {
	update();
}

/// @brief Set what happens when the user drags in the viewport.
void
SimplePoseDrawOpenGLWidget::set_drag_mode(
	SPDOGLW_drag_mode const mode_in
) {
	debug_assert( mode_in > 0 && mode_in <= SPDOGLW_end_of_drag_mode_list );
	drag_mode_ = mode_in;
}


////////////////PROTECTED METHODS///////////////////////////////////

/// @brief Initialize the OpenGL rendering environment.
void
SimplePoseDrawOpenGLWidget::initializeGL() {
	glClearColor( 0.1, 0.2, 0.3, 1.0 );
	glShadeModel(GL_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	//glEnable(GL_CULL_FACE);

	glEnable(GL_LIGHT1);
	GLfloat lightAmbient1[] = { 0.1f, 0.15f, 0.2f, 1.0f };
	GLfloat lightDiffuse1[] = { 1.0f, 0.95f, 0.92f, 1.0f };
	GLfloat lightPosition1[] = { -200.0f, 200.0f, -55.0f, 1.0f};
	glLightfv(GL_LIGHT1, GL_AMBIENT, lightAmbient1);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightDiffuse1);
	glLightfv(GL_LIGHT1, GL_SPECULAR, lightDiffuse1);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPosition1);

	glEnable(GL_LIGHT2);
	GLfloat lightAmbient2[] = { 0.0f, 0.0f, 0.0f, 1.0f };
	GLfloat lightDiffuse2[] = { 0.3f, 0.45f, 0.5f, 1.0f };
	GLfloat lightPosition2[] = { 200.0f, -200.0f, -15.0f, 1.0f};
	glLightfv(GL_LIGHT2, GL_AMBIENT, lightAmbient2);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, lightDiffuse2);
	glLightfv(GL_LIGHT2, GL_SPECULAR, lightDiffuse2);
	glLightfv(GL_LIGHT2, GL_POSITION, lightPosition2);
}

/// @brief What to do when the window has been resized.
void
SimplePoseDrawOpenGLWidget::resizeGL(int width, int height) {
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	GLfloat const aspectratio( GLfloat(width) / GLfloat(height) );
	GLfloat const zNear(0.2), zFar(1000.0);
	GLfloat const fH( tan( 70.0 / 360.0 * 3.141592654 ) * zNear );
	GLfloat const fW( fH * aspectratio);

	glFrustum( -fW, fW, -fH, fH, zNear, zFar );

	glMatrixMode(GL_MODELVIEW);
}

/// @brief How to actually draw a pose.
void
SimplePoseDrawOpenGLWidget::paintGL() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(pose_ == nullptr) {
		glFlush();
		return; //Do nothing if no pose
	}

	core::select::residue_selector::ResidueSubset selection( pose_->total_residue(), false);
	if( colour_mode_ == SPDOGLW_colour_by_selection  && residue_selector_ != nullptr ) {
		selection = residue_selector_->apply(*pose_);
	}

	glEnable(GL_LIGHTING);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, forward_back_);
	glRotatef(rotx_, 1.0f, 0.0f, 0.0f);
	glRotatef(roty_, 0.0f, 1.0f, 0.0f);
	glRotatef(rotz_, 0.0f, 0.0f, 1.0f);

	GLfloat colorvect[] = { 0.25f, 0.25f, 0.3f, 1.0f};
	GLfloat colorvectselected[] = { 0.98f, 0.95f, 0.04f, 1.0f };
	GLfloat colorvectN[] = { 0.15f, 0.15f, 0.75f, 1.0f};
	GLfloat colorvectO[] = { 0.75f, 0.15f, 0.15f, 1.0f};
	GLfloat colorvectS[] = { 0.85f, 0.82f, 0.07f, 1.0f};
	GLfloat colorvectH[] = { 0.75f, 0.75f, 0.75f, 1.0f};
	// Colors for less-common elements (in proteins) that are 
	// nonetheless quite common in nucleic acids. Editorial
	// decision to echo PyMOL colors, which may or may not
	// work out in the lighting environment available.
	// Phosphorus
	GLfloat colorvectP[] = { 1.00f, 0.502f, 0.00f, 1.0f};
	// Halogens
	GLfloat colorvectF[] = { 0.702f, 1.00f, 1.00f, 1.0f};
	GLfloat colorvectCl[] = { 0.122f, 0.941f, 0.122f, 1.0f};
	GLfloat colorvectBr[] = { 0.651f, 0.161f, 0.161f, 1.0f};
	GLfloat colorvectI[] = { 0.580f, 0.00f, 0.580f, 1.0f};
	
	// Ca, Mg, Na
	GLfloat colorvectCa[] = { 0.239f, 1.00f, 0.00f, 1.0f};
	GLfloat colorvectMg[] = { 0.541f, 1.00f, 0.00f, 1.0f};
	GLfloat colorvectNa[] = { 0.671f, 0.361f, 0.949f, 1.0f};
	
	GLfloat colorvectspec[] = { 0.2f, 0.2f, 0.2f, 1.0f};
	GLfloat colorvectenergy_max[] = { 0.9f, 0.6f, 0.1f, 1.0f };
	GLfloat colorvectenergy_min[] = { 0.01f, 0.4f, 0.8f, 1.0f };
	glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvect );
	glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, colorvectspec );
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 15.0f );

	GLfloat const energyrange( static_cast<GLfloat>( energy_range_max_ - energy_range_min_ ) );

	for(core::Size ir(1), irmax(pose_->total_residue()); ir<=irmax; ++ir) {
			for(core::Size ia(1), iamax(pose_->residue(ir).natoms()); ia<=iamax; ++ia) {
					numeric::xyzVector< core::Real > xyz( pose_->xyz( core::id::AtomID(ia, ir) ) );

					glTranslatef( xyz.x(), xyz.y(), xyz.z() );
					std::string const elem( pose_->residue(ir).atom_type(ia).element() );
					if( !elem.compare("N") ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectN );
					} else if( !elem.compare("S") ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectS );
					} else if( !elem.compare("O") ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectO );
					} else if( !elem.compare("P") ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectO );
					} else if( !elem.compare("F") ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectF );
					} else if( !elem.compare("Cl") ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectCl );
					} else if( !elem.compare("Br") ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectBr );
					} else if( !elem.compare("I") ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectI );
					} else if( !elem.compare("Ca") ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectCa );
					} else if( !elem.compare("Mg") ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectMg );
					} else if( !elem.compare("Na") ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectNa );
					} else if( pose_->residue(ir).atom_type(ia).is_polar_hydrogen() ) {
							glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectH );
					} else {
							if( colour_mode_ == SPDOGLW_colour_by_total_energy || colour_mode_ == SPDOGLW_colour_by_score_term ) {
								GLfloat const totalenergy(
									colour_mode_ == SPDOGLW_colour_by_total_energy ?
									static_cast<GLfloat>( pose_->energies().residue_total_energies(ir)[ core::scoring::total_score ] ) :
									static_cast<GLfloat>( pose_->energies().residue_total_energies(ir)[ score_type_ ] )
								);
								GLfloat energyfract( (totalenergy - static_cast<GLfloat>( energy_range_min_ )) / energyrange );
								if(energyfract > 1) energyfract = 1.0f;
								if(energyfract < 0) energyfract = 0.0f;
								GLfloat const invenergyfract( 1.0f - energyfract );
								GLfloat colorvectenergy[] = { energyfract * colorvectenergy_max[0] + invenergyfract * colorvectenergy_min[0], energyfract * colorvectenergy_max[1] + invenergyfract * colorvectenergy_min[1], energyfract * colorvectenergy_max[2] + invenergyfract * colorvectenergy_min[2], 1.0f  };
								glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, colorvectenergy );
							} else {
								glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE,
									colour_mode_ == SPDOGLW_colour_by_selection && selection[ir] ? colorvectselected : colorvect
								);
							}
					}
					//qDebug() << "About to draw sphere at " << xyz.x() << ", " << xyz.y() << ", " << xyz.z();
					glutSolidSphere( 0.65*pose_->residue(ir).atom_type(ia).lj_radius(), 32, 16 );
					glTranslatef( -xyz.x(), -xyz.y(), -xyz.z() );
			}
	}

	glDisable(GL_LIGHTING);
	glFinish();
}

/// @brief What to do when the user clicks.
void
SimplePoseDrawOpenGLWidget::mousePressEvent( QMouseEvent* event ) {
	lastpos_ = event->pos();
}

/// @brief What do do when the user drags.
void
SimplePoseDrawOpenGLWidget::mouseMoveEvent(QMouseEvent* event) {
	if(event->buttons() == Qt::LeftButton ) {
		GLfloat const dx( static_cast<GLfloat>( event->x() - lastpos_.x() ) / static_cast<GLfloat>( width() ) );
		GLfloat const dy( static_cast<GLfloat>( event->y() - lastpos_.y() ) / static_cast<GLfloat>( height() ) );
		rotx_ += 180*dy;
		rotz_ -= 180*dx;
		if(rotx_ > 360.0) rotx_ -= 360.0;
		else if(rotx_ < 0.0) rotx_ += 360.0;
		if(rotz_ > 360.0) rotz_ -= 360.0;
		else if(rotz_ < 0.0) rotz_ += 360.0;
		lastpos_ = event->pos();
		update();
	} else if(event->buttons() == Qt::RightButton) {
		GLfloat const dz( static_cast<GLfloat>( event->y() - lastpos_.y() ) / static_cast<GLfloat>( height() ) );
		forward_back_ -= 50*dz;
		if(forward_back_ > -1.0) forward_back_ = -1.0;
		lastpos_ = event->pos();
		update();
	}
}



} //namespace pose_draw
} //namespace ui_core
} //namespace ui
