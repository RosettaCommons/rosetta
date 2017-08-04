// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/ui_protocols/HelicalBundleDialogueWidget.h
/// @brief  Headers for a dialogue widget that displays controls for manipulating
/// Crick parameters during helical bundle parametric design.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef HELICALBUNDLEDIALOGUEWIDGET_H
#define HELICALBUNDLEDIALOGUEWIDGET_H

#include <QWidget>
#include <ui/ui_protocols/helical_bundle/HelicalBundleDialogueWidget.fwd.h>
#include <ui/ui_protocols/helical_bundle/HelixOptionWidget.fwd.h>

//Rosetta utility includes
#include <utility/vector1.fwd.hh>

//Rosetta core includes
#include <core/types.hh>

namespace Ui {
class HelicalBundleDialogueWidget; //Forward declaration of class used for UI.
}

namespace ui {
namespace ui_protocols {
namespace helical_bundle {

class HelicalBundleDialogueWidget : public QWidget
{
	Q_OBJECT

public:
	explicit HelicalBundleDialogueWidget(QWidget *parent = 0);
	~HelicalBundleDialogueWidget();

	/// @brief Get r0
	core::Real r0( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const;

	/// @brief Get omega0
	core::Real omega0( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const;

	/// @brief Get delta_omega0
	core::Real delta_omega0( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const;

	/// @brief Get delta_omega1
	core::Real delta_omega1( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const;

	/// @brief Get delta_t
	core::Real delta_t( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const;

	/// @brief Get z0_offset
	core::Real z0_offset( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const;

	/// @brief Get z1_offset
	core::Real z1_offset( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const;

	/// @brief Get epsilon
	core::Real epsilon( utility::vector1< HelicalBundleDialogueWidget* > prev_widgets ) const;

	/// @brief Get helix length (residues).
	core::Size helix_length() const;

	/// @brief Get whether to flip the helix direction.
	bool invert_helix() const;

	/// @brief Get the Crick params file to use for the current helix.
	/// @details Currently hard-coded.  Could change later.
	std::string params_file() const;

	/// @brief Set the index of the helix that this control controls.
	void set_helix_index( core::Size const new_index );

	/// @brief Get the index of the helix that this control controls.
	core::Size helix_index() const;

	/// @brief If this helix copies its pitch from another, get the index of the helix
	/// from which we're copying.
	/// @details Returns 0 if we're not copying pitch (even if we're copying omega0
	/// directly from another helix).
	core::Size omega0_pitchcopy_helix() const;

	/// @brief If any controls are set to track another helix with helix index other_index,
	/// reset these controls to set a value directly.
	void reset_controls_dependent_on_other_helix( core::Size const other_index );

	/// @brief Given the index of a helix that is to be removed, decrement the indices of helices PAST that helix
	/// on which THIS helix depends.
	void update_control_dependencies_given_helix_to_be_removed( core::Size const helix_to_remove );

	/// @brief Is a custom params file currently selected?
	bool custom_params_file_selected() const;

	/// @brief Get the custom params file.
	std::string const & custom_params_file() const;

	/// @brief Set r0.
	void set_r0( core::Real const &setting );

	/// @brief Set omega0.
	void set_omega0( core::Real const &setting );

	/// @brief Set delta_omega0.
	void set_delta_omega0( core::Real const &setting );

	/// @brief Set delta_omega1.
	void set_delta_omega1( core::Real const &setting );

	/// @brief Set delta_t.
	void set_delta_t( core::Real const &setting );

	/// @brief Set z1_offset.
	void set_z1_offset( core::Real const &setting );

	/// @brief Set z0_offset.
	void set_z0_offset( core::Real const &setting );

	/// @brief Set epsilon.
	void set_epsilon( core::Real const &setting );

	/// @brief Set everything except delta_omega0 to copy helix 1.
	/// @details Does nothing if this is helix 1.
	void set_everything_except_delta_omega0_copies_helix1();

	/// @brief Set whether this helix is inverted.
	void set_invert( bool const inverted );

	/// @brief Set the length of this helix.
	void set_helix_length( core::Size const length );

	/// @brief Set that this "helix" is actually a beta-strand.
	void set_to_beta_strand();

Q_SIGNALS:

	/// @brief One of the controls has changed; signal that this object has changed.
	void control_has_changed();

	/// @brief One of the controls has changed, requiring a complete rebuild of the pose.
	void control_has_changed_requiring_pose_rebuild();

	/// @brief This helix has been removed.  Signal to remove it!
	void helix_has_been_removed( core::Size const helix_index );

public Q_SLOTS:

	/// @brief What to do when the "remove helix" button is clicked.
	/// @details This is public to allow other things to remove helices.
	void on_remove_button_clicked();

private Q_SLOTS:

	/// @brief Sub-control has signaled that it has changed.  What should I do?  (Fire a control_has_changed() signal).
	void on_control_changed();

	/// @brief The first radio button has been toggled.
	void on_radioButton_l_alpha_helix_toggled(bool checked);

	/// @brief The second radio button has been toggled.
	void on_radioButton_beta_strand_toggled(bool checked);

	/// @brief The invert button has been toggled.
	void on_invert_checkbox_clicked();

	/// @brief The helix length spin box value has changed.
	void on_helix_length_spin_box_valueChanged(int arg1);

	void on_radioButton_custom_helix_clicked();

private:
	Ui::HelicalBundleDialogueWidget *ui;

	/// @brief What's the index of this helix?
	core::Size helix_index_;

	/// @brief Controls for setting r0.
	HelixOptionWidget* r0_widget_;

	/// @brief Controls for setting omega0.
	HelixOptionWidget* omega0_widget_;

	/// @brief Controls for setting delta_omega0.
	HelixOptionWidget* delta_omega0_widget_;

	/// @brief Controls for setting delta_omega1.
	HelixOptionWidget* delta_omega1_widget_;

	/// @brief Controls for setting delta_t.
	HelixOptionWidget* delta_t_widget_;

	/// @brief Controls for setting z0_offset.
	HelixOptionWidget* z0_offset_widget_;

	/// @brief Controls for setting z1_offset.
	HelixOptionWidget* z1_offset_widget_;

	/// @brief Controls for setting epsilon.
	HelixOptionWidget* epsilon_widget_;

	/// @brief Name of a custom params file to use.
	std::string custom_params_file_;

};

} //helical_bundle
} //protocols
} //ui

#endif // HELICALBUNDLEDIALOGUEWIDGET_H
