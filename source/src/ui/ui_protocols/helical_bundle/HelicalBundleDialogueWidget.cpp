// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/ui_protocols/HelicalBundleDialogueWidget.cpp
/// @brief  A dialogue widget that displays controls for manipulating
/// Crick parameters during helical bundle parametric design.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include "HelicalBundleDialogueWidget.h"
#include "ui_HelicalBundleDialogueWidget.h"

#include <ui/ui_protocols/helical_bundle/HelixOptionWidget.h>

//Rosetta utility headers
#include <utility/vector1.hh>

//Rosetta basic headers
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//C++ headers
#include <sstream>

//Qt headers
#include <QFileDialog>
#include <QFileInfo>

namespace ui {
namespace ui_protocols {
namespace helical_bundle {

HelicalBundleDialogueWidget::HelicalBundleDialogueWidget(QWidget *parent) :
	QWidget(parent),
	ui(new Ui::HelicalBundleDialogueWidget),
	helix_index_(0),
	r0_widget_( new HelixOptionWidget(this) ),
	omega0_widget_( new HelixOptionWidget(this) ),
	delta_omega0_widget_( new HelixOptionWidget(this) ),
	delta_omega1_widget_( new HelixOptionWidget(this) ),
	delta_t_widget_( new HelixOptionWidget(this) ),
	z0_offset_widget_( new HelixOptionWidget(this) ),
	z1_offset_widget_( new HelixOptionWidget(this) ),
	epsilon_widget_( new HelixOptionWidget(this) ),
	custom_params_file_( "14_helix" )
{
	ui->setupUi(this);

	ui->verticalLayout_3->addWidget(r0_widget_);
	ui->verticalLayout_3->addWidget(omega0_widget_);
	ui->verticalLayout_3->addWidget(delta_omega0_widget_);
	ui->verticalLayout_3->addWidget(delta_omega1_widget_);
	ui->verticalLayout_3->addWidget(delta_t_widget_);
	ui->verticalLayout_3->addWidget(z0_offset_widget_);
	ui->verticalLayout_3->addWidget(z1_offset_widget_);
	ui->verticalLayout_3->addWidget(epsilon_widget_);

	r0_widget_->configure( "Radius from z-axis (r0)", 0, 30, false );
	r0_widget_->show();
	r0_widget_->set_value(5.0);
	omega0_widget_->configure( "Twist per residue about z-axis (omega0)", -25.0, 25.0, true );
	omega0_widget_->show();
	omega0_widget_->set_value(0.0);
	delta_omega0_widget_->configure( "Rotation about z-axis (delta_omega0)", 0, 360, false );
	delta_omega0_widget_->show();
	delta_omega0_widget_->set_value(0.0);
	delta_omega1_widget_->configure( "Rotation about helix axis (delta_omega1)", 0, 360, false );
	delta_omega1_widget_->show();
	delta_omega1_widget_->set_value(0.0);
	delta_t_widget_->configure( "Offset along peptide backbone (delta_t)", -10, 10, false );
	delta_t_widget_->show();
	delta_t_widget_->set_value(0.0);
	z0_offset_widget_->configure( "Offset along z-axis (z0_offset)", -10, 10, false );
	z0_offset_widget_->show();
	z0_offset_widget_->set_value(0.0);
	z1_offset_widget_->configure( "Offset along helix axis (z1_offset)", -10, 10, false );
	z1_offset_widget_->show();
	z1_offset_widget_->set_value(0.0);
	epsilon_widget_->configure( "Bundle eccentricity (epsilon)", 0, 2, false );
	epsilon_widget_->show();
	epsilon_widget_->set_value(1.0);

	connect( r0_widget_, SIGNAL(control_has_changed()), this, SLOT(on_control_changed()) );
	connect( omega0_widget_, SIGNAL(control_has_changed()), this, SLOT(on_control_changed()) );
	connect( delta_omega0_widget_, SIGNAL(control_has_changed()), this, SLOT(on_control_changed()) );
	connect( delta_omega1_widget_, SIGNAL(control_has_changed()), this, SLOT(on_control_changed()) );
	connect( delta_t_widget_, SIGNAL(control_has_changed()), this, SLOT(on_control_changed()) );
	connect( z0_offset_widget_, SIGNAL(control_has_changed()), this, SLOT(on_control_changed()) );
	connect( z1_offset_widget_, SIGNAL(control_has_changed()), this, SLOT(on_control_changed()) );
	connect( epsilon_widget_, SIGNAL(control_has_changed()), this, SLOT(on_control_changed()) );

	update();
	Q_EMIT control_has_changed();
}

HelicalBundleDialogueWidget::~HelicalBundleDialogueWidget()
{
	delete ui;
}

/// @brief Get r0
core::Real HelicalBundleDialogueWidget::r0( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const {
	if ( r0_widget_->widget_mode() == HOW_copy_value ) {
		core::Size const helix_to_copy( r0_widget_->helix_to_copy() );
		debug_assert( helix_to_copy > 0 && helix_to_copy < helix_index_ );
		return prev_widgets[ helix_to_copy ]->r0( prev_widgets );
	}
	return r0_widget_->value();
}

/// @brief Get omega0
core::Real HelicalBundleDialogueWidget::omega0( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const {
	if ( omega0_widget_->widget_mode() == HOW_copy_value ) {
		core::Size const helix_to_copy( omega0_widget_->helix_to_copy() );
		debug_assert( helix_to_copy > 0 && helix_to_copy < helix_index_ );
		return prev_widgets[ helix_to_copy ]->omega0( prev_widgets );
	}
	return omega0_widget_->value();
}

/// @brief Get delta_omega0
core::Real HelicalBundleDialogueWidget::delta_omega0( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const {
	if ( delta_omega0_widget_->widget_mode() == HOW_copy_value ) {
		core::Size const helix_to_copy( delta_omega0_widget_->helix_to_copy() );
		debug_assert( helix_to_copy > 0 && helix_to_copy < helix_index_ );
		return prev_widgets[ helix_to_copy ]->delta_omega0( prev_widgets );
	}
	return delta_omega0_widget_->value();
}

/// @brief Get delta_omega1
core::Real HelicalBundleDialogueWidget::delta_omega1( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const {
	if ( delta_omega1_widget_->widget_mode() == HOW_copy_value ) {
		core::Size const helix_to_copy( delta_omega1_widget_->helix_to_copy() );
		debug_assert( helix_to_copy > 0 && helix_to_copy < helix_index_ );
		return prev_widgets[ helix_to_copy ]->delta_omega1( prev_widgets );
	}
	return delta_omega1_widget_->value();
}

/// @brief Get delta_t
core::Real HelicalBundleDialogueWidget::delta_t( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const {
	if ( delta_t_widget_->widget_mode() == HOW_copy_value ) {
		core::Size const helix_to_copy( delta_t_widget_->helix_to_copy() );
		debug_assert( helix_to_copy > 0 && helix_to_copy < helix_index_ );
		return prev_widgets[ helix_to_copy ]->delta_t( prev_widgets );
	}
	return delta_t_widget_->value();
}

/// @brief Get z0_offset
core::Real HelicalBundleDialogueWidget::z0_offset( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const {
	if ( z0_offset_widget_->widget_mode() == HOW_copy_value ) {
		core::Size const helix_to_copy( z0_offset_widget_->helix_to_copy() );
		debug_assert( helix_to_copy > 0 && helix_to_copy < helix_index_ );
		return prev_widgets[ helix_to_copy ]->z0_offset( prev_widgets );
	}
	return z0_offset_widget_->value();
}

/// @brief Get z1_offset
core::Real HelicalBundleDialogueWidget::z1_offset( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const {
	if ( z1_offset_widget_->widget_mode() == HOW_copy_value ) {
		core::Size const helix_to_copy( z1_offset_widget_->helix_to_copy() );
		debug_assert( helix_to_copy > 0 && helix_to_copy < helix_index_ );
		return prev_widgets[ helix_to_copy ]->z1_offset( prev_widgets );
	}
	return z1_offset_widget_->value();
}

/// @brief Get epsilon
core::Real HelicalBundleDialogueWidget::epsilon( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const {
	if ( epsilon_widget_->widget_mode() == HOW_copy_value ) {
		core::Size const helix_to_copy( epsilon_widget_->helix_to_copy() );
		debug_assert( helix_to_copy > 0 && helix_to_copy < helix_index_ );
		return prev_widgets[ helix_to_copy ]->epsilon( prev_widgets );
	}
	return epsilon_widget_->value();
}

/// @brief Get helix length (residues).
core::Size HelicalBundleDialogueWidget::helix_length() const { return ui->helix_length_spin_box->value(); }

/// @brief Get whether to flip the helix direction.
bool HelicalBundleDialogueWidget::invert_helix() const { return ui->invert_checkbox->isChecked(); }

/// @brief Get the Crick params file to use for the current helix.
/// @details Currently hard-coded.  Could change later.
std::string
HelicalBundleDialogueWidget::params_file() const {
	if( ui->radioButton_l_alpha_helix->isChecked() ) { return "alpha_helix.crick_params"; }
	else if( ui->radioButton_beta_strand->isChecked() ) { return "beta_strand.crick_params"; }
	else if( ui->radioButton_custom_helix->isChecked() ) { return custom_params_file_; }

	return("alpha_helix.crick_params"); //Default
}

/// @brief Set the index of the helix that this control controls.
void
HelicalBundleDialogueWidget::set_helix_index(
	core::Size const new_index
) {
	helix_index_ = new_index;
	std::stringstream newtext;
	newtext << "Helix " << new_index;
	ui->label->setText(QString( newtext.str().c_str() ));
	ui->label->update();

	bool const copying_allowed( new_index > 1 );
	r0_widget_->set_copying_enabled(copying_allowed, new_index - 1);
	omega0_widget_->set_copying_enabled(copying_allowed, new_index - 1);
	delta_omega0_widget_->set_copying_enabled(copying_allowed, new_index - 1);
	delta_omega1_widget_->set_copying_enabled(copying_allowed, new_index - 1);
	delta_t_widget_->set_copying_enabled(copying_allowed, new_index - 1);
	z1_offset_widget_->set_copying_enabled(copying_allowed, new_index - 1);
	z0_offset_widget_->set_copying_enabled(copying_allowed, new_index - 1);
	epsilon_widget_->set_copying_enabled(copying_allowed, new_index - 1);
}

/// @brief Get the index of the helix that this control controls.
core::Size HelicalBundleDialogueWidget::helix_index() const { return helix_index_; }

/// @brief If this helix copies its pitch from another, get the index of the helix
/// from which we're copying.
/// @details Returns 0 if we're not copying pitch (even if we're copying omega0
/// directly from another helix).
core::Size
HelicalBundleDialogueWidget::omega0_pitchcopy_helix() const {
	if(omega0_widget_->widget_mode() != HOW_copy_pitch) return 0;
	return omega0_widget_->helix_to_copy();
}

/// @brief If any controls are set to track another helix with helix index other_index,
/// reset these controls to set a value directly.
void
HelicalBundleDialogueWidget::reset_controls_dependent_on_other_helix(
	core::Size const other_index
) {
	if(r0_widget_->widget_mode() == HOW_copy_value && r0_widget_->helix_to_copy() == other_index) {
		r0_widget_->set_value(5.0, true);
	}
	if( (omega0_widget_->widget_mode() == HOW_copy_value || omega0_widget_->widget_mode() == HOW_copy_pitch) && omega0_widget_->helix_to_copy() == other_index) {
		omega0_widget_->set_value(0.0, true);
	}
	if( delta_omega0_widget_->widget_mode() == HOW_copy_value && delta_omega0_widget_->helix_to_copy() == other_index) {
		delta_omega0_widget_->set_value(0.0, true);
	}
	if( delta_omega1_widget_->widget_mode() == HOW_copy_value && delta_omega1_widget_->helix_to_copy() == other_index) {
		delta_omega1_widget_->set_value(0.0, true);
	}
	if( delta_t_widget_->widget_mode() == HOW_copy_value && delta_t_widget_->helix_to_copy() == other_index) {
		delta_t_widget_->set_value(0.0, true);
	}
	if( z1_offset_widget_->widget_mode() == HOW_copy_value && z1_offset_widget_->helix_to_copy() == other_index) {
		z1_offset_widget_->set_value(0.0, true);
	}
	if( z0_offset_widget_->widget_mode() == HOW_copy_value && z0_offset_widget_->helix_to_copy() == other_index) {
		z0_offset_widget_->set_value(0.0, true);
	}
	if( epsilon_widget_->widget_mode() == HOW_copy_value && epsilon_widget_->helix_to_copy() == other_index) {
		epsilon_widget_->set_value(1.0, true);
	}
}

/// @brief Is a custom params file currently selected?
bool
HelicalBundleDialogueWidget::custom_params_file_selected() const {
	return ui->radioButton_custom_helix->isChecked();
}

/// @brief Get the custom params file.
std::string const &
HelicalBundleDialogueWidget::custom_params_file() const {
	return custom_params_file_;
}

/// @brief Set r0.
void
HelicalBundleDialogueWidget::set_r0(
	core::Real const &setting
) {
	r0_widget_->set_value( setting );
	Q_EMIT control_has_changed();
}

/// @brief Set omega0.
void
HelicalBundleDialogueWidget::set_omega0(
	core::Real const &setting
) {
	omega0_widget_->set_value( setting );
	Q_EMIT control_has_changed();
}

/// @brief Set delta_omega0.
void
HelicalBundleDialogueWidget::set_delta_omega0(
	core::Real const &setting
) {
	delta_omega0_widget_->set_value( setting );
	Q_EMIT control_has_changed();
}

/// @brief Set delta_omega1.
void
HelicalBundleDialogueWidget::set_delta_omega1(
	core::Real const &setting
) {
	delta_omega1_widget_->set_value( setting );
	Q_EMIT control_has_changed();
}

/// @brief Set delta_t.
void
HelicalBundleDialogueWidget::set_delta_t(
	core::Real const &setting
) {
	delta_t_widget_->set_value( setting );
	Q_EMIT control_has_changed();
}

/// @brief Set z1_offset.
void
HelicalBundleDialogueWidget::set_z1_offset(
	core::Real const &setting
) {
	z1_offset_widget_->set_value( setting );
	Q_EMIT control_has_changed();
}

/// @brief Set z0_offset.
void
HelicalBundleDialogueWidget::set_z0_offset(
	core::Real const &setting
) {
	z0_offset_widget_->set_value( setting );
	Q_EMIT control_has_changed();
}

/// @brief Set epsilon.
void
HelicalBundleDialogueWidget::set_epsilon(
	core::Real const &setting
) {
	epsilon_widget_->set_value( setting );
	Q_EMIT control_has_changed();
}

/// @brief Set everything except delta_omega0 to copy helix 1.
/// @details Does nothing if this is helix 1.
void
HelicalBundleDialogueWidget::set_everything_except_delta_omega0_copies_helix1() {
	if( helix_index_ == 1 ) return;
	r0_widget_->set_copies_helix1();
	omega0_widget_->set_copies_helix1_pitch();
	delta_omega1_widget_->set_copies_helix1();
	delta_t_widget_->set_copies_helix1();
	z1_offset_widget_->set_copies_helix1();
	z0_offset_widget_->set_copies_helix1();
	epsilon_widget_->set_copies_helix1();
	Q_EMIT control_has_changed();
}

/// @brief Set whether this helix is inverted.
void
HelicalBundleDialogueWidget::set_invert(
	bool const inverted
) {
	ui->invert_checkbox->setChecked(inverted);
	ui->invert_checkbox->update();
	Q_EMIT control_has_changed();
}

/// @brief Set the length of this helix.
void
HelicalBundleDialogueWidget::set_helix_length(
	core::Size const length
) {
	ui->helix_length_spin_box->setValue(length);
	ui->helix_length_spin_box->update();
	Q_EMIT control_has_changed_requiring_pose_rebuild();
}

/// @brief Set that this "helix" is actually a beta-strand.
void
HelicalBundleDialogueWidget::set_to_beta_strand() {
	ui->radioButton_beta_strand->setChecked(true);
	ui->radioButton_beta_strand->update();
	Q_EMIT control_has_changed_requiring_pose_rebuild();
}

/// @brief Given the index of a helix that is to be removed, decrement the indices of helices PAST that helix
/// on which THIS helix depends.
void
HelicalBundleDialogueWidget::update_control_dependencies_given_helix_to_be_removed(
	core::Size const helix_to_remove
) {
	if(r0_widget_->widget_mode() == HOW_copy_value && r0_widget_->helix_to_copy() > helix_to_remove) {
		r0_widget_->decrement_copy_from_helix();
	}
	if( (omega0_widget_->widget_mode() == HOW_copy_value || omega0_widget_->widget_mode() == HOW_copy_pitch) && omega0_widget_->helix_to_copy() > helix_to_remove) {
		omega0_widget_->decrement_copy_from_helix();
	}
	if(delta_omega0_widget_->widget_mode() == HOW_copy_value && delta_omega0_widget_->helix_to_copy() > helix_to_remove) {
		delta_omega0_widget_->decrement_copy_from_helix();
	}
	if(delta_omega1_widget_->widget_mode() == HOW_copy_value && delta_omega1_widget_->helix_to_copy() > helix_to_remove) {
		delta_omega1_widget_->decrement_copy_from_helix();
	}
	if(delta_t_widget_->widget_mode() == HOW_copy_value && delta_t_widget_->helix_to_copy() > helix_to_remove) {
		delta_t_widget_->decrement_copy_from_helix();
	}
	if(z1_offset_widget_->widget_mode() == HOW_copy_value && z1_offset_widget_->helix_to_copy() > helix_to_remove) {
		z1_offset_widget_->decrement_copy_from_helix();
	}
	if(z0_offset_widget_->widget_mode() == HOW_copy_value && z0_offset_widget_->helix_to_copy() > helix_to_remove) {
		z0_offset_widget_->decrement_copy_from_helix();
	}
	if(epsilon_widget_->widget_mode() == HOW_copy_value && epsilon_widget_->helix_to_copy() > helix_to_remove) {
		epsilon_widget_->decrement_copy_from_helix();
	}
}


//////////////////////////////PRIVATE SLOTS/////////////////////////////

/// @brief Sub-control has signaled that it has changed.  What should I do?  (Fire a control_has_changed() signal).
void HelicalBundleDialogueWidget::on_control_changed() {
	Q_EMIT control_has_changed();
}


} //helical_bundle
} //protocols
} //ui

/// @brief The helix length spin box value has changed.
void ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget::on_helix_length_spin_box_valueChanged(int /*arg1*/)
{
	Q_EMIT control_has_changed_requiring_pose_rebuild();
}

/// @brief The first radio button has been toggled.
void ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget::on_radioButton_l_alpha_helix_toggled(bool /*checked*/)
{
	Q_EMIT control_has_changed_requiring_pose_rebuild();
}

/// @brief The second radio button has been toggled.
void ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget::on_radioButton_beta_strand_toggled(bool /*checked*/)
{
	Q_EMIT control_has_changed_requiring_pose_rebuild();
}

/// @brief The invert button has been toggled.
void ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget::on_invert_checkbox_clicked()
{
	Q_EMIT control_has_changed_requiring_pose_rebuild();
}

/// @brief What to do when the "remove helix" button is clicked.
/// @details This is public to allow other things to remove helices.
void ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget::on_remove_button_clicked()
{
	Q_EMIT helix_has_been_removed( helix_index() );
}

void ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget::on_radioButton_custom_helix_clicked()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	QString const filename( QFileDialog::getOpenFileName(this, "Load crick params file", QString( (option[ in::path::database ](1).name() + "/protocol_data/crick_parameters/").c_str() ), "Crick params file (*.crick_params)") );
	if( filename.isEmpty() ) return;
	QFileInfo fileinfo(filename);

	custom_params_file_ = fileinfo.baseName().toStdString();
	ui->radioButton_custom_helix->setText( QString( std::string( "Custom: " + custom_params_file_ ).c_str() ) );

	Q_EMIT control_has_changed_requiring_pose_rebuild();
}
