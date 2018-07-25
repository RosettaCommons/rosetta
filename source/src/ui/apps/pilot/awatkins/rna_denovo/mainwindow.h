// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/apps/pilot/awatkins/rna_denovo/mainwindow.h
/// @brief  Headers for main window class for the stepwise GUI project.
/// @author Andrew M. Watkins (amw579@stanford.edu)

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QTimer>

//Rosetta numeric headres
#include <numeric/xyzVector.hh>
#include <numeric/xyzTransform.hh>

//Rosetta core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

//Rosetta ui headres
#include <ui/ui_protocols/helical_bundle/HelicalBundlePoseDrawOpenGLWidget.fwd.h>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private: //Private functions

	/// @brief Use the PerturbBundle mover to update the pose.
	void update_pose();

	/// @brief Use the MakeBundle mover to build the pose.
	void rebuild_pose_from_scratch();

private Q_SLOTS:

    /// @brief What do to when a control reports that the user has changed a value.  (Update the pose).
	void on_control_changed();

	/// @brief What do do when a control reports that the user has changed a value, AND a full rebuild of the pose is required.
	void on_control_changed_requiring_pose_rebuild();

	/// @brief Periodically, check whether we need to update the pose and the display, and do so if we need to.
	void on_timer_refresh();

	/// @brief The user has selected that the pose be coloured by total score (per residue).  Update accordingly.
    void on_radiobutton_colourby_totalscore_clicked();

private: //DATA:

    Ui::MainWindow *ui;

	/// @brief A pose on which to operate.
    core::pose::PoseOP pose_;

	/// @brief A scorefunction, used for colouring the pose.
	core::scoring::ScoreFunctionOP sfxn_;

	/// @brief Do we need to update the pose?
	bool need_to_update_pose_;

	/// @brief Do we need to rebuild the pose from scratch?
	bool need_to_rebuild_pose_;

	/// @brief A QTimer used to refresh the display
	QTimer* timer_;

	/// @brief A widget to draw the pose.
	ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget* pose_draw_widget_;

};

#endif // MAINWINDOW_H
