// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/ui_protocols/HelicalBundlePoseDrawOpenGLWidget.cpp
/// @brief  A simple widget that shows a slider and other options that can be set
/// for a single parameter during helical bundle design.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include "HelixOptionWidget.h"
#include "ui_HelixOptionWidget.h"

//C++ includes
#include <math.h>

namespace ui {
namespace ui_protocols {
namespace helical_bundle {

HelixOptionWidget::HelixOptionWidget(QWidget *parent) :
	QWidget(parent),
	ui(new Ui::HelixOptionWidget),
	update_source_( HOW_update_from_neither )
{
	ui->setupUi(this);
}

HelixOptionWidget::~HelixOptionWidget()
{
	delete ui;
}

/// @brief Set the options for this dialogue.
/// @param[in] name The name of the option.
/// @param[in] range_min The min value of the range of allowed values.
/// @param[in] range_max The max value of the range of allowed values.
/// @param[in] show_pitch_option Iff true, the option to copy pitch from helix is shown.
void
HelixOptionWidget::configure(
	std::string const &name,
	core::Real const &range_min,
	core::Real const &range_max,
	bool const show_pitch_option
) {
	ui->groupBox->setTitle(QString( name.c_str() ));
	ui->doubleSpinBox->setMinimum(range_min);
	ui->doubleSpinBox->setMaximum(range_max);
	ui->doubleSpinBox->setDecimals(2);
	if(!show_pitch_option) {
		ui->radioButton_3->hide();
	}
}

/// @brief Set the value of the slider and the spinner to some value.
/// @details If reset_mode is true, then the HelixOptionWidget is also reset so that
/// it is in direct-setting mode.
void
HelixOptionWidget::set_value(
	core::Real const &value,
	bool const reset_mode/*=false*/
) {
	ui->doubleSpinBox->setValue(value);
	update_source_ = HOW_update_from_spinner;
	update_slider_spinner_values();

	if( reset_mode ) { //If we're also resetting:
		ui->spinBox->setValue( ui->radioButton_2->isEnabled() ? 1 : 0 );
		ui->radioButton_1->setChecked(true);
		on_radioButton_1_clicked();
	}
}

/// @brief Get the value of the slider and spinner.
core::Real
HelixOptionWidget::value() const {
	return ui->doubleSpinBox->value();
}

/// @brief Get the current widget mode.
HOW_widget_mode
HelixOptionWidget::widget_mode() const {
	if(ui->radioButton_1->isChecked()) {
		return HOW_set_value;
	} else if ( ui->radioButton_2->isChecked() ) {
		return HOW_copy_value;
	} else if ( ui->radioButton_3->isChecked() ) {
		return HOW_copy_pitch;
	}
	return HOW_set_value;
}

/// @brief Get the value of the spinner that indicates which helix should be copied.
core::Size
HelixOptionWidget::helix_to_copy() const {
	return static_cast<core::Size>( ui->spinBox->value() );
}

/// @brief Set whether the "copy from helix" option is allowed or not.
void
HelixOptionWidget::set_copying_enabled(
	bool const setting,
	core::Size const max_index
) {
	if(setting==false) {
		ui->radioButton_1->setChecked(true);
	}
	ui->radioButton_2->setEnabled( setting );
	ui->radioButton_3->setEnabled( setting );

	ui->spinBox->setMinimum(1);
	ui->spinBox->setMaximum(max_index);
	if(ui->spinBox->value() < 1 || static_cast<core::Size>( ui->spinBox->value() ) > max_index) { ui->spinBox->setValue(1);}

	ui->radioButton_1->update();
	ui->radioButton_2->update();
	if(ui->radioButton_3->isVisible()) ui->radioButton_3->update();
	ui->spinBox->update();
}

/// @brief If the "copy_from_helix" option is enabled, reduce the index of the helix from which to
/// copy by one.
void
HelixOptionWidget::decrement_copy_from_helix() {
	if( !ui->radioButton_1->isChecked() ) {
		int const newval( ui->spinBox->value() - 1 );
		if( newval > 0 ) ui->spinBox->setValue( newval );
	}
	ui->spinBox->update();
}

/// @brief Set this widget to copy helix 1.
void
HelixOptionWidget::set_copies_helix1() {
	ui->radioButton_2->setChecked(true);
	ui->spinBox->setValue(1);
	ui->radioButton_2->update();
	ui->spinBox->update();
	on_radioButton_2_clicked();
}


/// @brief Set this widget to copy pitch from helix 1.
/// @details Assumes that this is an omega0 widget.
void
HelixOptionWidget::set_copies_helix1_pitch() {
	ui->radioButton_3->setChecked(true);
	ui->spinBox->setValue(1);
	ui->radioButton_3->update();
	ui->spinBox->update();
	on_radioButton_3_clicked();
}

/// @brief Update the slider from the spinner or the spinner from the slider.
void
HelixOptionWidget::update_slider_spinner_values() {
	if(update_source_ == HOW_update_from_slider) {
		double const absval( static_cast<double>( ui->horizontalSlider->value() ) / 1000.0 * (ui->doubleSpinBox->maximum() - ui->doubleSpinBox->minimum() ) + ui->doubleSpinBox->minimum() );
		ui->doubleSpinBox->setValue(absval);
	} else if (update_source_ == HOW_update_from_spinner) {
		int const fractval( round( 1000.0*(ui->doubleSpinBox->value() - ui->doubleSpinBox->minimum()) / (ui->doubleSpinBox->maximum() - ui->doubleSpinBox->minimum() )) );
		ui->horizontalSlider->setValue( fractval );
	}
	update_source_ = HOW_update_from_neither;
}



/////////////////////PRIVATE SLOTS/////////////////

/// @brief If the user has selected "set to value".
void HelixOptionWidget::on_radioButton_1_clicked()
{
	ui->spinBox->setEnabled(false);
	ui->doubleSpinBox->setEnabled(true);
	ui->horizontalSlider->setEnabled(true);
	Q_EMIT control_has_changed();
}

} //namespace helical_bundle
} //namespace ui_protocols
} //namespace ui

/// @brief If the user has selected "copy from helix".
void ui::ui_protocols::helical_bundle::HelixOptionWidget::on_radioButton_2_clicked()
{
	ui->spinBox->setEnabled(true);
	ui->doubleSpinBox->setEnabled(false);
	ui->horizontalSlider->setEnabled(false);
	Q_EMIT control_has_changed();
}

/// @brief If the user has selected "copy pitch from helix".
void ui::ui_protocols::helical_bundle::HelixOptionWidget::on_radioButton_3_clicked()
{
	ui->spinBox->setEnabled(true);
	ui->doubleSpinBox->setEnabled(false);
	ui->horizontalSlider->setEnabled(false);
	Q_EMIT control_has_changed();
}

/// @brief The user has changed the value in the spin box.  Update the corresponding slider.
void ui::ui_protocols::helical_bundle::HelixOptionWidget::on_doubleSpinBox_valueChanged(double /*arg1*/)
{
	if(update_source_ == HOW_update_from_neither) {
		update_source_ = HOW_update_from_spinner;
		update_slider_spinner_values();
	}

	/*int const fractval( round( 1000.0*(ui->doubleSpinBox->value() - ui->doubleSpinBox->minimum()) / (ui->doubleSpinBox->maximum() - ui->doubleSpinBox->minimum() )) );
	if(fractval != ui->horizontalSlider->value()) {
		ui->horizontalSlider->setValue( fractval );
	}*/
	Q_EMIT control_has_changed();
}

/// @brief The user has moved the slider.  Update the corresponding spin box.
void ui::ui_protocols::helical_bundle::HelixOptionWidget::on_horizontalSlider_valueChanged(int /*value*/)
{
	if(update_source_ == HOW_update_from_neither) {
		update_source_ = HOW_update_from_slider;
		update_slider_spinner_values();
	}

	/*double const absval( static_cast<double>( ui->horizontalSlider->value() ) / 1000.0 * (ui->doubleSpinBox->maximum() - ui->doubleSpinBox->minimum() ) + ui->doubleSpinBox->minimum() );
	if( abs(ui->doubleSpinBox->value() - absval) > 0.001 ) {
		ui->doubleSpinBox->setValue(absval);
	}*/
	Q_EMIT control_has_changed();
}

/// @brief The user has changed which helix we're copying from.  Emit a control_has_changed() signal.
void ui::ui_protocols::helical_bundle::HelixOptionWidget::on_spinBox_valueChanged(int /*arg1*/)
{
	Q_EMIT control_has_changed();
}
