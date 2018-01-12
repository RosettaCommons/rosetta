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
#include <ui/ui_protocols/helical_bundle/HelicalBundleDialogueWidget.h>
#include <ui/ui_protocols/helical_bundle/HelicalBundlePoseDrawOpenGLWidget.h>

//C++ headers
#include <sstream>







// Unit headers
#include <protocols/viewer/ConformationViewer.hh>

// Package headers
#include <protocols/viewer/viewers.hh>

#include <core/chemical/rna/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>

#include <core/scoring/electron_density/ElectronDensity.hh>

#include <core/types.hh>

// Project headers
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
// using -edensity::mapfile option to determine whether to try to display density contours or not

// C++ Headers

// GLUT

#include <pthread.h>










/// @brief observer that attaches to a Conformation and displays graphics
class ConformationViewer : public utility::pointer::ReferenceCount {

private: // typedefs

    typedef utility::pointer::ReferenceCount Super;

public: // typedefs

    typedef utility::vector1< core::conformation::ResidueCOP > ResidueCOPs;

public: // construct/destruct

    /// @brief default constructor
    ConformationViewer( ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget * pdw );

    /// @brief constructor
    ConformationViewer( std::string const & name, ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget * pdw );

    /// @brief constructor
    ConformationViewer( std::string const & name, int length, int width, bool debug_pause, ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget * pdw );

    // @brief default destructor
    ~ConformationViewer();


private: // disallow copy

    /// @brief disallow copy constructor
    // NOTE: if implementing copy constructor, remember to set 'conf_' to NULL
    //       as there is no transferal of subject Conformation on copy construct
    ConformationViewer( ConformationViewer const & rval );

    /// @brief disallow copy assignment
    // NOTE: if ConformationViewer copy assignment, remember to leave 'conf_' untouched
    //       as any current subject Conformation is kept on copy assign
    ConformationViewer & operator =( ConformationViewer const & rval );


public: // window management

    /// called by glutDisplayFunc
    void
    display_func();

    /// called by glutIdleFunc
    void
    display_if_necessary();


    int
    window() const
    {
        return my_window_;
    }


    void
    window( int const setting )
    {
        my_window_ = setting;
    }

    int
    get_width() const
    {
        return width_;
    }

    int
    get_length() const
    {
        return length_;
    }

    protocols::viewer::GraphicsState &
    get_gs() {
        return current_gs_;
    }


    std::string const &
    name() const
    {
        return name_;
    }

public: // observer interface

    /// @brief is currently observing a Conformation?
    /// @return the Conformation, otherwise NULL
    core::conformation::Conformation const *
    is_observing() const;

    /// @brief attach to Conformation
    void
    attach_to( core::conformation::Conformation const & conf );

    /// @brief detach from Conformation
    void
    detach_from();

    /// @brief upon receiving a ConnectionEvent do...
    void
    on_connection_change( core::conformation::signals::ConnectionEvent const & event );

    /// @brief upon receiving a GeneralEvent update the residues and atom tree root
    void
    on_xyz_change( core::conformation::signals::XYZEvent const & event );

    void
    set_center_vector( core::Vector const & setting ){ center_vector_defined_ = true; center_vector_ = setting; }

private:
    ResidueCOPs residues_;
    utility::vector1< char > secstruct_;
    protocols::viewer::GraphicsState current_gs_;

    // density isosurface
    utility::vector1< protocols::viewer::triangle > triangles_;

    core::id::AtomID anchor_id_;

    std::string const name_;

    bool new_conformation_;

    int my_window_;
    int length_;
    int width_;

    bool use_debug_pause_;

    bool center_vector_defined_;
    core::Vector center_vector_;

    utility::signals::Link connection_event_link_;
    utility::signals::Link xyz_event_link_;

    /// @brief the Conformation being observed, we need this to ensure
    ///  debug_pause is reset upon detachment
    core::conformation::Conformation const * conf_;

    ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget * pose_draw_widget_ = nullptr;

    pthread_mutex_t residues_mut_;
};

typedef std::shared_ptr< ConformationViewer > ConformationViewerOP;
typedef std::shared_ptr< ConformationViewer const > ConformationViewerCOP;


using namespace core;

// signal that our window hasn't been created yet
int const BAD_WINDOW( -999 );

/// @brief default constructor
ConformationViewer::ConformationViewer( ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget * pdw ) :
    Super(),
    new_conformation_( true ),
    use_debug_pause_( false ),
    center_vector_defined_( false ),
    conf_( NULL ),
    pose_draw_widget_( pdw )
{}

// CALLED BY WORKER THREAD
ConformationViewer::ConformationViewer(std::string const & name_in, ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget * pdw ):
    Super(),
    name_( name_in ),
    new_conformation_( true ),
    my_window_( BAD_WINDOW ),
    length_( 900 ),
    width_( 900 ),
    use_debug_pause_( false ),
    center_vector_defined_( false ),
    conf_( NULL ),
    pose_draw_widget_( pdw )
{
    pthread_mutex_init( &residues_mut_, NULL );
}

ConformationViewer::ConformationViewer(std::string const & name_in, int length, int width, bool debug_pause, ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget * pdw ):
    Super(),
    name_( name_in ),
    new_conformation_( true ),
    my_window_( BAD_WINDOW ),
    length_( length ),
    width_( width ),
    use_debug_pause_( debug_pause ),
    center_vector_defined_( false ),
    conf_( NULL ),
    pose_draw_widget_( pdw )
{
    pthread_mutex_init( &residues_mut_, NULL );
}

/// @brief default destructor
ConformationViewer::~ConformationViewer()
{
    detach_from();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// called by our window's glutDisplayFunc, conformation_viewer_display, when glut registers need for
// redrawing our window
//
// CALLED BY GRAPHICS THREAD

void
ConformationViewer::display_func()
{

    if ( residues_.empty() ) return;

    // lock the residues vector
    pthread_mutex_lock( &residues_mut_ );

    pose_draw_widget_->update_pose_draw();


    pthread_mutex_unlock( &residues_mut_ );
}

/////////////////////////////////////////////////////////////////

// called by glutIdleFunc -- GRAPHICS THREAD
void
ConformationViewer::display_if_necessary()
{

    if ( new_conformation_ && my_window_ != BAD_WINDOW ) {
        //std::cout << "display_if_necessary new_conf=true" << std::endl;
        new_conformation_ = false;
        display_func();

    }

}

/////////////////////////////////////////////////////////////////
/// @brief is currently observing a Conformation?
/// @return the Conformation, otherwise NULL
core::conformation::Conformation const *
ConformationViewer::is_observing() const
{
    return conf_;
}

/////////////////////////////////////////////////////////////////
/// @brief attach to Conformation
void
ConformationViewer::attach_to(
    core::conformation::Conformation const & conf
)
{
    detach_from();

    xyz_event_link_ = conf.attach_xyz_obs( &ConformationViewer::on_xyz_change, this );
    connection_event_link_ = conf.attach_connection_obs( &ConformationViewer::on_connection_change, this );

    conf.debug_pause( use_debug_pause_ );

    conf_ = &conf; // keeping pointer only so debug_pause() can be cleared later
}

/////////////////////////////////////////////////////////////////
/// @brief detach from Conformation
void
ConformationViewer::detach_from()
{
    xyz_event_link_.invalidate();
    connection_event_link_.invalidate();

    if ( conf_ ) {
        conf_->debug_pause( false );
    }

    conf_ = 0;

// next line is commented out under the assumption that we may want
// to draw the last state of the Conformation even after detachment
//	residues_.clear();
}

/////////////////////////////////////////////////////////////////
/// @brief upon receiving a ConnectionEvent do...
void
ConformationViewer::on_connection_change(
    core::conformation::signals::ConnectionEvent const & event
)
{
    using core::conformation::signals::ConnectionEvent;

    switch ( event.tag ) {
        case ConnectionEvent::DISCONNECT:
            detach_from();
            break;
        case ConnectionEvent::TRANSFER:
            // connection is being transferred, so swap the Conformation pointer
            attach_to( *event.conformation );
            break;
        default: // do nothing
            break;
    }
}

/////////////////////////////////////////////////////////////////
/// @brief upon receiving a GeneralEvent update the residues and atom tree root
void
ConformationViewer::on_xyz_change(
    core::conformation::signals::XYZEvent const & event
)
{
    pthread_mutex_lock( &residues_mut_ );

    // Grab the residues.  This is a roundabout way of doing it, but
    // at the moment it's the only const access to the residues and
    // it doesn't trigger any residue updates inside the Conformation.
    // The alternative is to change all viewer function calls to using
    // CAPs instead, but for now we don't do that in case display of
    // the residues persist after a Conformation is destroyed.
    core::conformation::ResidueCAPs res_caps = event.conformation->const_residues();
    residues_.resize( res_caps.size() );
    secstruct_.resize( res_caps.size() );
    for ( core::Size i = 1, ie = res_caps.size(); i <= ie; ++i ) {
        residues_[ i ] = core::conformation::ResidueCOP( res_caps[ i ] );
        secstruct_[ i ] = event.conformation->secstruct( i );
    }

    core::kinematics::tree::AtomCOP root_atom = event.conformation->atom_tree().root();
    if ( root_atom ) {
        anchor_id_ = root_atom->id();
        core::conformation::ResidueCOP rsd_root = residues_[ anchor_id_.rsd() ];
        if ( rsd_root->is_RNA() ) {
            anchor_id_ = id::AtomID( rsd_root->atom_index( core::chemical::rna::default_jump_atom( rsd_root->type() ) ), anchor_id_.rsd() );
        }
    }

    new_conformation_ = true;

    // always set the debug pause in case it gets reset elsewhere
    event.conformation->debug_pause( use_debug_pause_ );

    pose_draw_widget_->update_pose_draw();

    pthread_mutex_unlock( &residues_mut_ );
}

void
add_conformation_viewer(
    conformation::Conformation & conformation,
    ui::ui_protocols::helical_bundle::HelicalBundlePoseDrawOpenGLWidget * pdw,
    std::string const & name_in, // = ""
    int const length,
    int const width ,
    bool const debug_pause,
    bool const set_center_vector ,
    core::Vector const center_vector
)
{

    //pthread_mutex_lock( &new_conformation_viewers_mut );

    // create a new viewer
    std::string const window_name
        ( name_in.empty() ? "conformation" : name_in );

    ConformationViewerOP viewer( new ConformationViewer( window_name, length, width, debug_pause, pdw ) );
    if ( set_center_vector ) viewer->set_center_vector( center_vector );

    viewer->attach_to( conformation );

    //new_conformation_viewers.push_back( viewer );

    //pthread_mutex_unlock( &new_conformation_viewers_mut );

    // allow main to start if this is the 1st window
    //pthread_cond_broadcast( &start_cond );

}











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
    QString const sequence( ui->lineEdit->text() );
    if( sequence.isEmpty() ) return;
    QString const secstruct( ui->lineEdit_2->text() );
    if( secstruct.isEmpty() ) return;

    // Set up a RNA_DeNovoProtocol with these parameters.
    utility::options::OptionKeyList opts;
    // OK, gotta use list_options_read functions provided by RNA denovo stuff, maybe???
    using namespace protocols::rna::denovo;
	core::import_pose::options::RNA_DeNovoProtocolOptions::list_options_read( opts );

    utility::options::OptionCollection faux_cl_opts = basic::options::option;
    using namespace basic::options::OptionKeys;
    using namespace core::import_pose;
    // TODO: enable provision of sequence and secstruct file paths.
    faux_cl_opts[ rna::denovo::minimize_rna ].set_cl_value( "true" );
    faux_cl_opts[ rna::denovo::sequence ].set_cl_value( sequence.toStdString() );
    faux_cl_opts[ rna::denovo::secstruct ].set_cl_value( secstruct.toStdString() );
    RNA_DeNovoSetup rna_de_novo_setup;
    rna_de_novo_setup.initialize_from_options( faux_cl_opts );
    pose_ = rna_de_novo_setup.pose();
    pose_draw_widget_->set_pose( pose_ );

    RNA_DeNovoProtocol protocol( rna_de_novo_setup.options(), rna_de_novo_setup.rna_params() );
    pose_draw_widget_->update_pose_draw();

    add_conformation_viewer( pose_->conformation(), pose_draw_widget_, "foo", 900, 900, false, false, core::Vector(0,0,0) );

    protocol.apply( *pose_ );
    pose_draw_widget_->update_pose_draw();
    //pose_->dump_pdb( pdb_filename );
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
