// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/ui_core/SimplePoseDrawOpenGLWidget.h
/// @brief  Headers for a base class for drawing a simple Van der Waals representation of a pose, with
/// the option to colour atoms based on a residue selection or on per-residue score.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef SIMPLEPOSEDRAWOPENGLWIDGET_H
#define SIMPLEPOSEDRAWOPENGLWIDGET_H

#include <ui/ui_core/pose_draw/SimplePoseDrawOpenGLWidget.fwd.h>

#include <QObject>
#include <QOpenGLWidget>
#include <QtOpenGL/QGL>

//Rosetta core includes
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreType.hh>

//OpenGL includes
#if defined(MAC) || defined(__APPLE__) || defined (__OSX__)
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

namespace ui {
namespace ui_core {
namespace pose_draw {

class SimplePoseDrawOpenGLWidget : public QOpenGLWidget
{
	Q_OBJECT

public:

	/// @brief Enum for the different ways that the pose could be coloured.
	enum class ColorMode {
						   none, //Keep first
						   selection,
						   input_domain,
						   total_energy,
						   score_term, //Keep second-to-last
						   list = score_term //Keep last
	};

	enum class DragMode {
						 rotate_and_zoom_viewport, // = 1, //Keep first
						 list = rotate_and_zoom_viewport //Keep last.  TODO: if another is added, keep it second-to-last and set this to equal its value.
	};


	/// @brief Constructor
	SimplePoseDrawOpenGLWidget( QWidget* parent=nullptr );

	/// @brief Set the pose that we'll be drawing.
	void set_pose( core::pose::PoseCOP pose );

	/// @brief Get currently assigned pose object
	core::pose::PoseCOP pose() const { return pose_; }

	/// @brief Set the residue selector that we'll use for colouring the pose.
	/// @details Only used if colour_mode_ == SPODGLW_colour_by_selection.
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Set the colour mode that we'll use for colouring the pose.
	void set_color_mode( ColorMode mode);

	/// @brief Set the lower and upper ends of the energy range in the gradient of colour values used.
	/// @details There is no check that upper > lower.
	void set_energy_range( core::Real const &lower, core::Real const &upper );

	/// @brief Set the score term to use if we're colouring by a single score type.
	/// @details Only used if colour_mode_ == SPDOGLW_colour_by_score_term.
	void set_score_type( core::scoring::ScoreType scoretype_in );

	/// @brief The pose has changed; redraw it.
	void update_pose_draw();

	/// @brief Set what happens when the user drags in the viewport.
	void set_drag_mode( DragMode mode);

protected:

	/// @brief Initialize the OpenGL rendering environment.
	void initializeGL() override;

	/// @brief What to do when the window has been resized.
	void resizeGL(int width, int height) override;

	/// @brief How to actually draw a pose.
	void paintGL() override;

	/// @brief What to do when the user clicks.
	void mousePressEvent( QMouseEvent* event ) override;

	/// @brief What do do when the user drags.
	void mouseMoveEvent(QMouseEvent* event) override;

	/// @brief Read the x-rotation value.
	inline GLfloat const & rotx() const { return rotx_; }

	/// @brief Read the y-rotation value.
	inline GLfloat const & roty() const { return roty_; }

	/// @brief Read the z-rotation value.
	inline GLfloat const & rotz() const { return rotz_; }

	/// @brief The start of a mouse click-and-drag event.
	inline QPoint const & lastpos() const { return lastpos_; }

	/// @brief Set the value of lastpos.
	inline void set_lastpos( QPoint const &point_in ) { lastpos_ = point_in; }

private: //Private variables

	/// @brief Const owning pointer to the pose that we'll be drawing.
	///
	core::pose::PoseCOP pose_;

	/// @brief Const owning pointer to the residue selector that we'll use for
	/// colouring the pose.
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;

	/// @brief Rotation of the scene about the x-axis.
	GLfloat rotx_;

	/// @brief Rotation of the scene about the y-axis.
	GLfloat roty_;

	/// @brief Rotation of the scene about the z-axis.
	GLfloat rotz_;

	/// @brief How far have I moved backwards or forwards?
	GLfloat forward_back_;

	/// @brief The start of a mouse click-and-drag event.
	QPoint lastpos_;

	/// @brief The way in which we're currently colouring the pose.
	///
	ColorMode color_mode_ = ColorMode::none;

	/// @brief The energy value at the far upper end of the energy colour gradient.
	///
	core::Real energy_range_max_;

	/// @brief The energy value at the far lower end of the energy colour gradient.
	///
	core::Real energy_range_min_;

	/// @brief The score term to use if we're colouring by a single score type.
	/// @details Only used if colour_mode_ == SPDOGLW_colour_by_score_term.
	core::scoring::ScoreType score_type_;

	/// @brief The behaviour of the viewport when the user clicks and drags.
	DragMode drag_mode_ = DragMode::rotate_and_zoom_viewport;

};

} //namespace pose_draw
} //namespace ui_core
} //namespace ui

#endif // SIMPLEPOSEDRAWOPENGLWIDGET_H
