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
	subwidgets_(), //Initialized below
	custom_params_file_( "14_helix" )
{
	using namespace protocols::helical_bundle;

	for( core::Size i(1); i<=static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++i ) {
		subwidgets_[i] = new HelixOptionWidget(this);
	}

	ui->setupUi(this);

	for( core::Size i(1); i<=static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++i ) {
		ui->verticalLayout_3->addWidget(subwidgets_[i]);
	}

	subwidgets_[BPC_r0]->configure( "Radius from z-axis (r0)", 0, 30, false );
	subwidgets_[BPC_omega0]->configure( "Twist per residue about z-axis (omega0)", -25.0, 25.0, true );
	subwidgets_[BPC_delta_omega0]->configure( "Rotation about z-axis (delta_omega0)", 0, 360, false );
	subwidgets_[BPC_delta_omega1]->configure( "Rotation about helix axis (delta_omega1)", 0, 360, false );
	subwidgets_[BPC_delta_t]->configure( "Offset along peptide backbone (delta_t)", -10, 10, false );
	subwidgets_[BPC_z0_offset]->configure( "Offset along z-axis (z0_offset)", -10, 10, false );
	subwidgets_[BPC_z1_offset]->configure( "Offset along helix axis (z1_offset)", -10, 10, false );
	subwidgets_[BPC_epsilon]->configure( "Bundle eccentricity (epsilon)", 0, 2, false );

	for( core::Size i(1); i<=static_cast<core::Size>(BPC_last_parameter_to_be_sampled); ++i ) {
		subwidgets_[i]->show();
		if( i == static_cast< core::Size >( BPC_r0 ) ) {
			subwidgets_[i]->set_value(5.0);
		} else if( i == static_cast< core::Size >( BPC_epsilon ) ) {
			subwidgets_[i]->set_value(1.0);
		} else {
			subwidgets_[i]->set_value(0.0);
		}
		connect( subwidgets_[i], SIGNAL(control_has_changed()), this, SLOT(on_control_changed()) );
	}

	update();
	Q_EMIT control_has_changed();
}

HelicalBundleDialogueWidget::~HelicalBundleDialogueWidget()
{
	delete ui;
}

/// @brief Get a subwidget corresponding to a parameter.  (Const access).
HelixOptionWidget *
HelicalBundleDialogueWidget::get_subwidget_nonconst(
	core::Size const subwidget_index
) {
	debug_assert( subwidget_index > 0 );
	debug_assert( subwidget_index <= subwidgets_.size() );
	return subwidgets_[subwidget_index];
}

/// @brief Get a subwidget corresponding to a parameter.  (Const access).
HelixOptionWidget const *
HelicalBundleDialogueWidget::get_subwidget(
	core::Size const subwidget_index
) const {
	debug_assert( subwidget_index > 0 );
	debug_assert( subwidget_index <= subwidgets_.size() );
	return subwidgets_[subwidget_index];
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
	for( core::Size i(1); i<=static_cast<core::Size>(protocols::helical_bundle::BPC_last_parameter_to_be_sampled); ++i ) {
		subwidgets_[i]->set_copying_enabled( copying_allowed, new_index - 1 );
	}
}

/// @brief Get the index of the helix that this control controls.
core::Size HelicalBundleDialogueWidget::helix_index() const { return helix_index_; }

/// @brief If this helix copies its pitch from another, get the index of the helix
/// from which we're copying.
/// @details Returns 0 if we're not copying pitch (even if we're copying omega0
/// directly from another helix).
core::Size
HelicalBundleDialogueWidget::omega0_pitchcopy_helix() const {
	if(subwidgets_[protocols::helical_bundle::BPC_omega0]->widget_mode() != HOW_copy_pitch) return 0;
	return subwidgets_[protocols::helical_bundle::BPC_omega0]->helix_to_copy();
}

/// @brief If any controls are set to track another helix with helix index other_index,
/// reset these controls to set a value directly.
void
HelicalBundleDialogueWidget::reset_controls_dependent_on_other_helix(
	core::Size const other_index
) {
	for( core::Size i(1); i<=static_cast<core::Size>(protocols::helical_bundle::BPC_last_parameter_to_be_sampled); ++i ) {
		if(subwidgets_[i]->widget_mode() == HOW_copy_value && subwidgets_[i]->helix_to_copy() == other_index) {
			if( i == static_cast<core::Size>(protocols::helical_bundle::BPC_r0 ) ) {
				subwidgets_[i]->set_value(5.0, true);
			} else if ( i == static_cast<core::Size>(protocols::helical_bundle::BPC_epsilon ) ) {
				subwidgets_[i]->set_value(1.0, true);
			} else {
				subwidgets_[i]->set_value(0.0, true);
			}
		}
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

/// @brief Set the value of a subwidget.
void
HelicalBundleDialogueWidget::set_real_subwidget(
	core::Size const subwidget_index,
	core::Real const & value
) {
	debug_assert( subwidget_index > 0 );
	debug_assert( subwidget_index <= static_cast< core::Size >( protocols::helical_bundle::BPC_last_parameter_to_be_sampled ) );
	subwidgets_[subwidget_index]->set_value( value, true );
}

/// @brief Set everything except delta_omega0 to copy helix 1.
/// @details Does nothing if this is helix 1.
void
HelicalBundleDialogueWidget::set_everything_except_delta_omega0_copies_helix1() {
	if( helix_index_ == 1 ) return;

	for( core::Size i(1); i<=static_cast<core::Size>(protocols::helical_bundle::BPC_last_parameter_to_be_sampled); ++i ) {
		if( i == static_cast<core::Size>(protocols::helical_bundle::BPC_delta_omega0) ) {
			continue;
		} else if( i == static_cast<core::Size>(protocols::helical_bundle::BPC_omega0) ) {
			subwidgets_[i]->set_copies_helix1_pitch();
		} else {
			subwidgets_[i]->set_copies_helix1();
		}
	}
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
	for( core::Size i(1); i<=static_cast<core::Size>(protocols::helical_bundle::BPC_last_parameter_to_be_sampled); ++i ) {
		if(subwidgets_[i]->widget_mode() == HOW_copy_value && subwidgets_[i]->helix_to_copy() > helix_to_remove) {
			subwidgets_[i]->decrement_copy_from_helix();
		}
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
