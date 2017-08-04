// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/ui_protocols/HelicalBundlePoseDrawOpenGLWidget.h
/// @brief  Headers for a simple widget that shows a slider and other options that can be set
/// for a single parameter during helical bundle design.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef HELIXOPTIONWIDGET_H
#define HELIXOPTIONWIDGET_H

#include <ui/ui_protocols/helical_bundle/HelixOptionWidget.fwd.h>
#include <QWidget>

//Rosetta core includes
#include <core/types.hh>

namespace Ui {
class HelixOptionWidget;
}

namespace ui {
namespace ui_protocols {
namespace helical_bundle {

/// @brief Is the widget being used to set the value for a parameter for a helix, or to
/// copy it from another control?
enum HOW_widget_mode {
	HOW_set_value=1, //Keep first
	HOW_copy_value,
	HOW_copy_pitch, //Keep second-to-last
	HOW_end_of_list = HOW_copy_pitch //Keep last
};

/// @brief Is an update being made from the slider, the spinner, or neither?
enum HOW_update_source {
	HOW_update_from_neither=1, //Keep first
	HOW_update_from_slider,
	HOW_update_from_spinner, //Keep second-to-last
	HOW_end_of_update_source_list = HOW_update_from_spinner //Keep last
};

class HelixOptionWidget : public QWidget
{
	Q_OBJECT

public:
	explicit HelixOptionWidget(QWidget *parent = 0);
	~HelixOptionWidget();

	/// @brief Set the options for this dialogue.
	/// @param[in] name The name of the option.
	/// @param[in] range_min The min value of the range of allowed values.
	/// @param[in] range_max The max value of the range of allowed values.
	/// @param[in] show_pitch_option Iff true, the option to copy pitch from helix is shown.
	void configure( std::string const &name, core::Real const &range_min, core::Real const &range_max, bool const show_pitch_option );

	/// @brief Set the value of the slider and the spinner to some value.
	/// @details If reset_mode is true, then the HelixOptionWidget is also reset so that
	/// it is in direct-setting mode.
	void set_value( core::Real const &value, bool const reset_mode=false );

	/// @brief Get the value of the slider and spinner.
	core::Real value() const;

	/// @brief Get the current widget mode.
	HOW_widget_mode widget_mode() const;

	/// @brief Get the value of the spinner that indicates which helix should be copied.
	core::Size helix_to_copy() const;

	/// @brief Set whether the "copy from helix" option is allowed or not.
	void set_copying_enabled( bool const setting, core::Size const max_index );

	/// @brief If the "copy_from_helix" option is enabled, reduce the index of the helix from which to
	/// copy by one.
	void decrement_copy_from_helix();

	/// @brief Set this widget to copy helix 1.
	void set_copies_helix1();

	/// @brief Set this widget to copy pitch from helix 1.
	/// @details Assumes that this is an omega0 widget.
	void set_copies_helix1_pitch();

Q_SIGNALS:

	/// @brief Emitted when a control is changed.
	void control_has_changed();

private: //Functions

	/// @brief Update the slider from the spinner or the spinner from the slider.
	void update_slider_spinner_values();

private Q_SLOTS:

	/// @brief If the user has selected "set to value".
	void on_radioButton_1_clicked();

	/// @brief If the user has selected "copy from helix".
	void on_radioButton_2_clicked();

	/// @brief If the user has selected "copy pitch from helix".
	void on_radioButton_3_clicked();

	/// @brief The user has changed the value in the spin box.  Update the corresponding slider.
	void on_doubleSpinBox_valueChanged(double arg1);

	/// @brief The user has moved the slider.  Update the corresponding spin box.
	void on_horizontalSlider_valueChanged(int value);

	/// @brief The user has changed which helix we're copying from.  Emit a control_has_changed() signal.
	void on_spinBox_valueChanged(int arg1);

private:
	Ui::HelixOptionWidget *ui;

	/// @brief What is the source of a changed value -- the slider or the spinner?
	HOW_update_source update_source_;
};

} //namespace helical_bundle
} //namespace ui_protocols
} //namespace ui

#endif // HELIXOPTIONWIDGET_H
