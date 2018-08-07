// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/apps/pilot/awatkins/stepwise_interface/mainwindow.cpp
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
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/LayerSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/chemical/AA.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/import_pose/import_pose.hh>

//Rosetta protocols headers
#include <protocols/helical_bundle/MakeBundle.hh>
#include <protocols/helical_bundle/PerturbBundle.hh>
#include <protocols/helical_bundle/MakeBundleHelix.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>

//Rosetta ui headers
#include <ui/ui_protocols/helical_bundle/HelicalBundleDialogueWidget.h>
#include <ui/ui_protocols/helical_bundle/HelicalBundlePoseDrawOpenGLWidget.h>

//C++ headers
#include <sstream>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
	ui(new Ui::MainWindow),
	pose_(new core::pose::Pose),
	nonparametric_pose_(nullptr), //Null unless a pose is loaded.
	nonparametric_pose_offset_(numeric::xyzVector<core::Real>(0,0,0)),
	nonparametric_pose_rotation_( numeric::xyzTransform< core::Real >::identity() ),
	sfxn_(nullptr),
	layer_selector_( new core::select::residue_selector::LayerSelector ),
	residue_selector_( layer_selector_ ),
	need_to_update_pose_(false),
	need_to_rebuild_pose_(false),
	timer_(new QTimer),
	pose_draw_widget_( new ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget(this) )
{

    ui->setupUi(this);
	connect( timer_, SIGNAL(timeout()), this, SLOT( on_timer_refresh() ) );

	ui->tab_1->deleteLater();

	update_layer_selector_settings();

	ui->radiobutton_colourby_selection->setAutoExclusive(true);
	ui->radiobutton_colourby_totalscore->setAutoExclusive(true);

	ui->verticalLayout_opengl->addWidget(pose_draw_widget_);
	pose_draw_widget_->set_pose(pose_);
	pose_draw_widget_->set_colour_mode( ui::ui_core::pose_draw::SPDOGLW_colour_by_selection );
	pose_draw_widget_->set_residue_selector(residue_selector_);
	pose_draw_widget_->show();

	connect( pose_draw_widget_, SIGNAL(nonparametric_geometry_has_moved()), this, SLOT(on_nonparametric_geometry_moved()) );

	timer_->start(50);

}

MainWindow::~MainWindow()
{
    delete ui;
}

/////////////PRIVATE FUNCTIONS/////////////////////////////


/// @brief Use the BundleGridSampler to update the pose.
void
MainWindow::update_pose() {
	protocols::helical_bundle::PerturbBundle pertbundle;
	pertbundle.set_use_degrees(true);

	core::Size nhelices(0);
	utility::vector1< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* > previous_widgets;

	for(core::Size i(0), imax(ui->tabWidget->count()); i<imax; ++i) {
		QScrollArea* curscrollarea(  dynamic_cast< QScrollArea* >( ui->tabWidget->widget(i) ) );
		if(curscrollarea == nullptr) continue;
		ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* curwidget_ptr(  dynamic_cast< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* >( curscrollarea->widget() ) );
		if(curwidget_ptr == nullptr) continue;
		ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget& curwidget( *curwidget_ptr );

		++nhelices;
		debug_assert(curwidget.helix_index() == nhelices); //Should be true

		pertbundle.add_helix(nhelices);

		pertbundle.r0(nhelices)->set_being_set(true);
		pertbundle.r0(nhelices)->set_default_value( curwidget.r0(previous_widgets) );

		// Special-case logic for omega0 copying pitch instead of value:
		core::Size const omega0_pitchcopy_helix( curwidget.omega0_pitchcopy_helix() );
		if(omega0_pitchcopy_helix > 0) {
			pertbundle.omega0(nhelices)->set_use_defaults(false);
			pertbundle.omega0(nhelices)->set_helix_to_copy(omega0_pitchcopy_helix);
			pertbundle.omega0(nhelices)->set_omega0_copies_pitch_instead(true);
		} else {
			pertbundle.omega0(nhelices)->set_being_set(true);
			pertbundle.omega0(nhelices)->set_default_value( numeric::conversions::radians( curwidget.omega0(previous_widgets) ) );
		}

		pertbundle.delta_omega0(nhelices)->set_being_set(true);
		pertbundle.delta_omega0(nhelices)->set_default_value( numeric::conversions::radians( curwidget.delta_omega0(previous_widgets) ) );
		pertbundle.delta_omega1(nhelices)->set_being_set(true);
		pertbundle.delta_omega1(nhelices)->set_default_value( numeric::conversions::radians( curwidget.delta_omega1(previous_widgets) ) );
		pertbundle.delta_t(nhelices)->set_being_set(true);
		pertbundle.delta_t(nhelices)->set_default_value( curwidget.delta_t(previous_widgets) );
		pertbundle.z0_offset(nhelices)->set_being_set(true);
		pertbundle.z0_offset(nhelices)->set_default_value( curwidget.z0_offset(previous_widgets) );
		pertbundle.z1_offset(nhelices)->set_being_set(true);
		pertbundle.z1_offset(nhelices)->set_default_value( curwidget.z1_offset(previous_widgets) );
		pertbundle.epsilon(nhelices)->set_being_set(true);
		pertbundle.epsilon(nhelices)->set_default_value( curwidget.epsilon(previous_widgets) );

		previous_widgets.push_back(curwidget_ptr);
	}
	if(nhelices > 0) {
		pertbundle.apply(*pose_);
		//qDebug() << "Updating a helix"; //DELETE ME
	}

	if( sfxn_ != nullptr ) (*sfxn_)(*pose_);

	pose_draw_widget_->update_pose_draw();
}

void
MainWindow::rebuild_pose_from_scratch() {
	pose_->clear();

	if( nonparametric_pose_ != nullptr ) {
		(*pose_) = (*nonparametric_pose_);
		on_nonparametric_geometry_moved(); //Translate the geometry appropriately.
	}

	protocols::helical_bundle::MakeBundle mkbundle;
	mkbundle.set_use_degrees(true);
	mkbundle.set_reset_pose(pose_->total_residue() == 0);

	core::Size nhelices(0);
	utility::vector1< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* > previous_widgets;

	for(core::Size i(0), imax(ui->tabWidget->count()); i<=imax; ++i) {
		QScrollArea* curscrollarea(  dynamic_cast< QScrollArea* >( ui->tabWidget->widget(i) ) );
		if(curscrollarea == nullptr) continue;
		ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* curwidget_ptr(  dynamic_cast< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* >( curscrollarea->widget() ) );
		if(curwidget_ptr == nullptr) continue;
		ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget &curwidget( *curwidget_ptr );

		++nhelices;
		debug_assert(curwidget.helix_index() == nhelices); //Should be true.
		mkbundle.add_helix();

		protocols::helical_bundle::MakeBundleHelix &helix( *(mkbundle.helix(nhelices)) );

		helix.set_helix_length(curwidget.helix_length());
		helix.set_minor_helix_params_from_file(curwidget.params_file()); //NEED TO SET UP CACHING FOR THIS!

		helix.set_r0( curwidget.r0( previous_widgets ) );
		helix.set_omega0( numeric::conversions::radians( curwidget.omega0( previous_widgets ) ) );
		helix.set_delta_omega0( numeric::conversions::radians( curwidget.delta_omega0( previous_widgets ) ) );
		helix.set_delta_omega1_all( numeric::conversions::radians( curwidget.delta_omega1( previous_widgets ) ) );
		helix.set_delta_t( curwidget.delta_t( previous_widgets ) );
		helix.set_z0_offset( curwidget.z0_offset( previous_widgets ) );
		helix.set_z1_offset( curwidget.z1_offset( previous_widgets ) );
		helix.set_epsilon( curwidget.epsilon( previous_widgets ) );
		helix.set_invert_helix( curwidget.invert_helix() );

		previous_widgets.push_back( curwidget_ptr );
	}
	if(nhelices > 0) {
			mkbundle.apply(*pose_);
			update_pose(); //Needed for omega0 pitch copying.  Calls pose_draw_widget_->update_pose_draw().
	} else {
		pose_draw_widget_->update_pose_draw();
	}
}

/// @brief Update the indices of helices in UI elements.
void
MainWindow::update_helix_indices() {
	core::Size nhelices(0);
	for(core::Size i(0), imax(ui->tabWidget->count()); i<imax; ++i) {
		QScrollArea* curscrollarea(  dynamic_cast< QScrollArea* >( ui->tabWidget->widget(i) ) );
		if(curscrollarea == nullptr) continue;
		ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* curwidget(  dynamic_cast< ui::ui_protocols::helical_bundle::HelicalBundleDialogueWidget* >( curscrollarea->widget() ) );
		if(curwidget == nullptr) continue;
		++nhelices;
		curwidget->set_helix_index(nhelices);

		std::stringstream curname;
		curname << "Helix " << nhelices;
		ui->tabWidget->setTabText( i, QString( curname.str().c_str() ) );
	}
}

/// @brief Update the layer selector's settings from the UI.
void
MainWindow::update_layer_selector_settings() {
	layer_selector_->set_layers( ui->checkbox_select_core->isChecked() , ui->checkbox_select_boundary->isChecked(), ui->checkbox_select_surface->isChecked() );
}

/// @brief Switch the residue selector to a composite selector that does not select the imported geometry,
/// or just use the layer selector, depending on whether we have imported geometry.
void
MainWindow::update_residue_selector() {
	if( nonparametric_pose_ == nullptr ) {
		residue_selector_ = layer_selector_;
	} else {
		core::select::residue_selector::AndResidueSelectorOP andsel( new core::select::residue_selector::AndResidueSelector );
		core::select::residue_selector::NotResidueSelectorOP notsel( new core::select::residue_selector::NotResidueSelector );
		core::select::residue_selector::ResidueIndexSelectorOP indexsel( new core::select::residue_selector::ResidueIndexSelector );
		for(core::Size i(1), imax( nonparametric_pose_->total_residue() ); i<=imax; ++i) indexsel->append_index(i);
		notsel->set_residue_selector(indexsel);
		andsel->add_residue_selector(notsel);
		andsel->add_residue_selector(layer_selector_);
		residue_selector_ = andsel;
	}

	pose_draw_widget_->set_residue_selector(residue_selector_);
}

/// @brief Given a centre-of-mass vector, calculate the centre of mass of the nonparametric
/// geometry in the pose.
/// @details Does nothing (returns 0, 0, 0) if no nonparametric geometry is loaded.
void
MainWindow::compute_com_vect(
	numeric::xyzVector< core::Real > & comvect
) const {
	comvect = numeric::xyzVector< core::Real >(0, 0, 0);
	if( nonparametric_pose_ == nullptr ) return;

	core::pose::Pose const & pose( *pose_ ); //Get a reference to the pose for speed

	core::Size counter(0);
	for( core::Size ir(1), irmax(nonparametric_pose_->total_residue()); ir<=irmax; ++ir ) {
		for( core::Size ia(1), iamax( pose.residue_type(ir).natoms() ); ia<=iamax; ++ia ) {
			if( pose.residue_type(ir).atom(ia).is_hydrogen() ) continue;
			comvect += pose.xyz( core::id::AtomID(ia,ir) );
			++counter;
		}
	}
	if(counter > 0) {
		comvect /= static_cast<core::Real>(counter);
	}
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

/// @brief The "export PDB" button is clicked.  Dump a PDB file.
void MainWindow::on_pushButton_clicked()
{
	QString const filename( QFileDialog::getSaveFileName(this, "Save PDB file", "", "PDB file (*.pdb)") );
	if( filename.isEmpty() ) return;

	std::string pdb_filename( filename.toStdString() );

	/*if( !pdb_filename.substr(pdb_filename.length()-5, 4).compare(".pdb") ) {
		pdb_filename += ".pdb";
	}*/
	pose_->dump_pdb( pdb_filename );
}

/// @brief The user has checked or unchecked the allow design checkbox.
void MainWindow::on_checkbox_allow_design_clicked()
{
	bool const checked( ui->checkbox_allow_design->isChecked() );
	ui->checkbox_ala->setEnabled( checked );
	ui->checkbox_met->setEnabled( checked );
	ui->checkbox_aliphatic->setEnabled( checked );
	ui->checkbox_gly->setEnabled( checked );
	ui->checkbox_pro->setEnabled( checked );
	ui->checkbox_aromatic->setEnabled( checked );
	ui->checkbox_charged->setEnabled( checked );
	ui->checkbox_polar->setEnabled( checked );
	ui->checkbox_his->setEnabled( checked );
	ui->checkbox_cys->setEnabled( checked );
}

/// @brief The user has clicked on the "run minimizer" button.
void MainWindow::on_button_run_minimizer_clicked()
{
	core::select::residue_selector::ResidueSubset selection( residue_selector_->apply(*pose_) );
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb(false);
	movemap->set_chi(false);
	movemap->set_jump(false);
	for(core::Size i(1), imax(selection.size()); i<=imax; ++i) {
		movemap->set_chi( i, selection[i] );
	}

	protocols::minimization_packing::MinMover minmove( movemap, core::scoring::get_score_function(), "lbfgs_armijo_nonmonotone", 0.0000001, false );
	minmove.apply(*pose_);
	pose_draw_widget_->update_pose_draw();
}

/// @brief The user has selected that the pose be coloured by what's selected.  Update accordingly.
void MainWindow::on_radiobutton_colourby_selection_clicked()
{
	//ui->radiobutton_colourby_totalscore->setChecked(false);
	sfxn_ = nullptr;
	pose_draw_widget_->set_colour_mode( ui::ui_core::pose_draw::SPDOGLW_colour_by_selection );
	pose_draw_widget_->update_pose_draw();
}

/// @brief The user has selected that the pose be coloured by total score (per residue).  Update accordingly.
void MainWindow::on_radiobutton_colourby_totalscore_clicked()
{
	//ui->radiobutton_colourby_selection->setChecked(false);
	sfxn_ = core::scoring::get_score_function();
	(*sfxn_)(*pose_); //Score the pose.
	pose_draw_widget_->set_colour_mode( ui::ui_core::pose_draw::SPDOGLW_colour_by_total_energy );
	pose_draw_widget_->update_pose_draw();
}

/// @brief When the user clicks the "clear imported pose" button, clear the imported pose.
void MainWindow::on_clearimported_button_clicked()
{
    pose_.reset();
    pose_ = nullptr;
}

/// @brief When the user clicks the "import pose" button, import a pose.
void MainWindow::on_importPDB_button_clicked()
{
	QString const filename( QFileDialog::getOpenFileName(this, "Load PDB file", "", "PDB file (*.pdb)") );
	if( filename.isEmpty() ) return;

	std::string pdb_filename( filename.toStdString() );
    core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
    pose_ = core::import_pose::get_pdb_with_full_model_info( pdb_filename, rsd_set );
}

/// @brief When the user specifies that the mouse should be used for dragging the
/// nonparametric geometry around.
void MainWindow::on_tool_drag_nonparametric_clicked()
{
	pose_draw_widget_->set_extended_drag_mode( ui::ui_protocols::helical_bundle::HBPDOGLW_drag_nonparametric_geometry );
	pose_draw_widget_->determine_xyz_unit_vectors();
}

/// @brief When the user specifies that the mouse should be used for rotating the
/// viewing direction around.
void MainWindow::on_tool_rotate_zoom_view_clicked()
{
	pose_draw_widget_->set_drag_mode( ui::ui_core::pose_draw::SPDOGLW_rotate_and_zoom_viewport );
	pose_draw_widget_->set_extended_drag_mode( ui::ui_protocols::helical_bundle::HBPDOGLW_base_class_mode );
}

/// @brief When the user specifies that the mouse should be used for rotating the
/// nonparametric geometry around.
void MainWindow::on_tool_rotate_nonparametric_clicked()
{
	pose_draw_widget_->set_extended_drag_mode( ui::ui_protocols::helical_bundle::HBPDOGLW_rotate_nonparametric_geometry );
	pose_draw_widget_->determine_xyz_unit_vectors();
	numeric::xyzVector< core::Real > comvect(0,0,0);
	compute_com_vect( comvect );
	pose_draw_widget_->set_nonparametric_COM_vect( comvect );
}


/// @brief When the user indicates that the nonparametric geometry has moved.
void
MainWindow::on_nonparametric_geometry_moved() {
	if(nonparametric_pose_ == nullptr) return; //Do nothing if there's no nonparametric geometry.

	//Get Pose references for speed:
	core::pose::Pose const & nppose( *nonparametric_pose_ );
	core::pose::Pose & pose( *pose_ );

	nonparametric_pose_offset_ = pose_draw_widget_->get_nonparametric_translation_vect(); //Get the translation vector.
	nonparametric_pose_rotation_ = pose_draw_widget_->get_nonparametric_rotation_xform(); //Get the rotation transform.

	for(core::Size ir(1), irmax(nppose.total_residue()); ir<=irmax; ++ir) { //Loop through all residues
		for(core::Size ia(1), iamax( nppose.residue(ir).natoms() ); ia<=iamax; ++ia ) { //Loop through all atoms
			core::id::AtomID const curat( ia, ir );
			pose.set_xyz( curat, nonparametric_pose_rotation_*nppose.xyz( curat ) + nonparametric_pose_offset_ );
		}
	}
	pose.update_residue_neighbors();

	if( sfxn_ != nullptr ) {
		(*sfxn_)(pose);
	}

	pose_draw_widget_->update_pose_draw();
}
