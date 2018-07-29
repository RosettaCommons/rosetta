// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/apps/pilot/awatkins/rna_denovo/mainwindow.cpp
/// @brief  Main window class for the stepwise interface project.
/// @author Andrew Watkins (amw579@stanford.edu)

#include "mainwindow.h"
#include "ui_mainwindow.h"

//Qt includes:
#include <QScrollArea>
#include <QDebug>
#include <QVBoxLayout>
#include <QFileDialog>
#include <QDebug>

//Rosetta numeric headers
#include <numeric/conversions.hh>

//Rosetta core headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/AA.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/options/RNA_DeNovoProtocolOptions.hh>
#include <core/types.hh>
#include <utility/options/OptionCollection.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//Rosetta protocols headers
#include <protocols/rna/denovo/RNA_DeNovoProtocol.hh>
#include <core/import_pose/RNA_DeNovoSetup.hh>
#include <core/import_pose/options/RNA_DeNovoProtocolOptions.hh>

//Rosetta ui headers
#include <ui/ui_protocols/helical_bundle/HelicalBundlePoseDrawOpenGLWidget.h>

//C++ headers
#include <sstream>


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
	ui(new Ui::MainWindow),
	pose_(new core::pose::Pose),
	sfxn_(nullptr),
	need_to_update_pose_(false),
	need_to_rebuild_pose_(false),
	timer_(new QTimer),
	pose_draw_widget_( new ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget(this) )
{

    ui->setupUi(this);
	connect( timer_, SIGNAL(timeout()), this, SLOT( on_timer_refresh() ) );

    //ui->tab_1->deleteLater();
    //doubleSpinBoxPhi_->setMinimum(range_min);
    //doubleSpinBoxPsi_->setMaximum(range_max);
	ui->verticalLayout_opengl->addWidget(pose_draw_widget_);
	pose_draw_widget_->set_pose(pose_);
	pose_draw_widget_->set_color_mode( ui::ui_core::pose_draw::SimplePoseDrawOpenGLWidget::ColorMode::selection );
	pose_draw_widget_->show();

	timer_->start(50);

}

MainWindow::~MainWindow()
{
    delete ui;
}

/////////////PRIVATE FUNCTIONS/////////////////////////////


void
MainWindow::update_pose() {

	if( sfxn_ != nullptr ) (*sfxn_)(*pose_);

	pose_draw_widget_->update_pose_draw();
}

void
MainWindow::rebuild_pose_from_scratch() {
	pose_->clear();


    pose_draw_widget_->update_pose_draw();
}


/////////////PRIVATE Q_SLOTS///////////////////////////////

/// @brief What do to when a control reports that the user has changed a value.  (Update the pose).
void
MainWindow::on_control_changed() {
	need_to_update_pose_ = true;
}

/// @brief What do do when a control reports that the user has changed a value, AND a full rebuild of the pose is required.
void
MainWindow::on_control_changed_requiring_pose_rebuild() {
	need_to_rebuild_pose_ = true;
}

/// @brief Periodically, check whether we need to update the pose and the display, and do so if we need to.
void MainWindow::on_timer_refresh() {
	if( need_to_rebuild_pose_ ) {
		need_to_rebuild_pose_ = false;
		need_to_update_pose_ = false;
		rebuild_pose_from_scratch();
	} else if(need_to_update_pose_) {
		need_to_update_pose_ = false;
		need_to_rebuild_pose_ = false;
		update_pose();
		//pose_->dump_pdb("temp.pdb"); //DELETE ME -- for testing only.
	}
}

/// @brief The user has selected that the pose be coloured by total score (per residue).  Update accordingly.
void MainWindow::on_radiobutton_colourby_totalscore_clicked()
{
	//ui->radiobutton_colourby_selection->setChecked(false);
	sfxn_ = core::scoring::get_score_function();
	(*sfxn_)(*pose_); //Score the pose.
	pose_draw_widget_->set_color_mode( ui::ui_core::pose_draw::SimplePoseDrawOpenGLWidget::ColorMode::total_energy );
	pose_draw_widget_->update_pose_draw();
}
