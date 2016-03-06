// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

// Unit headers
#include <protocols/viewer/viewers.hh>

#include <protocols/viewer/SilentObserver.hh>
#include <protocols/viewer/SilentObserver.fwd.hh>
#include <core/chemical/AtomType.hh>

// Package headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarlo.tmpl.hh>
#include <core/id/AtomID.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

#if defined(WIN32) || defined(BOINC)
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <protocols/viewer/triangleIterator.hh>
#endif

#include <core/conformation/Residue.hh>


// Project headers

#ifdef GL_GRAPHICS
#include <protocols/viewer/triangleIterator.hh>
#include <protocols/viewer/ConformationViewer.hh>
#include <protocols/viewer/ConformationViewer.fwd.hh>
#include <protocols/viewer/triangleIterator.hh> //should not be auto-removed!! needed for graphics!
#include <protocols/viewer/ConformationViewer.hh>  //should not be auto-removed!! needed for graphics!
#include <protocols/viewer/ConformationViewer.fwd.hh>  //should not be auto-removed!! needed for graphics!
#include <ObjexxFCL/string.functions.hh>  //should not be auto-removed!! needed for graphics!
#include <core/chemical/ResidueTypeSet.hh>   //should not be auto-removed!! needed for graphics!
#include <core/chemical/ChemicalManager.hh>  //should not be auto-removed!! needed for graphics!
#include <core/chemical/AtomTypeSet.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#endif

// Utility headers
#include <utility/vector1.hh>

#include <numeric/NumericTraits.hh>

// GLUT
#if defined GL_GRAPHICS || defined BOINC_GRAPHICS

#ifdef MAC
#include <GLUT/glut.h>
#elif _WIN32
#include <glut/glut.h>
#else
#include "GL/glut.h"
#endif

#endif

#if defined( MAC ) || defined( __APPLE__ ) || defined( __OSX__ )
#define THREAD_STACK_SIZE 16 * 1024 * 1024 // 16 MB
#endif


using namespace core; /////////////////////////////////////////// DANGER

namespace protocols {
namespace viewer {

#ifdef GL_GRAPHICS
// prototypes
void check_for_new_conformation_viewers();

// global data
typedef std::map< int, ConformationViewerOP > ConformationViewers;
ConformationViewers conformation_viewers;
utility::vector1< ConformationViewerOP > new_conformation_viewers;

pthread_mutex_t new_conformation_viewers_mut = PTHREAD_MUTEX_INITIALIZER;

pthread_mutex_t start_mut = PTHREAD_MUTEX_INITIALIZER;

pthread_cond_t start_cond = PTHREAD_COND_INITIALIZER;

#endif

#if defined GL_GRAPHICS || defined BOINC_GRAPHICS

namespace graphics {
int window_size = 30;
int specialKey = 0;
bool click;
bool clicked_button;
int click_x;
int click_y;

// pointer to a graphics_state object for each ConformationViewer
std::map< int , GraphicsState* > gs_map_;

// Vector bg_color( 1.0f, 1.0f, 1.0f ); // white
Vector bg_color( 0.01f, 0.03f, 0.4f ); // dark blue
Vector bg_color2( 0.0f, 0.0f, 0.2f ); // darker blue
Vector border_color( 0.02f, 0.08f, 0.75f ); // lighter blue
Vector ghost_color_vect( 0.95f, 0.85f, 0.7f); //Orange-grey
Vector atom_specular_color(0.4f, 0.37f, 0.35f); //Grey-orange

// rhiju parameters for cartoons
const int NUM_SEGMENTS = 5;
const int NUM_SEGMENTS_COIL = 10;
const core::Real HELIX_HALF_WIDTH = 1.5;
const core::Real STRAND_HALF_WIDTH = 1.0;
const core::Real COILRADIUS = 0.2;
const core::Real HELIX_HERMITE_FACTOR = 4.7;
const core::Real STRAND_HERMITE_FACTOR = 4.7;
const core::Real COIL_HERMITE_FACTOR = 5.0;
const core::Real NA_HERMITE_FACTOR = 7.7;
const core::Real CHAINBREAK_CUTOFF2 = 4.5*4.5;
const core::Real CHAINBREAK_CUTOFF2_NA = 7.5*7.5;
const core::Real SHOWBONDCUTOFF2 = (5.0*5.0)/(NUM_SEGMENTS*NUM_SEGMENTS);
const core::Real SHOWBONDCUTOFF2_COIL = (5.0*5.0)/(NUM_SEGMENTS_COIL*NUM_SEGMENTS_COIL);
const core::Real SHOWBONDCUTOFF2_NA = (8.0*8.0)/(NUM_SEGMENTS*NUM_SEGMENTS);

//lin parameters for spacefill
core::Real const ligand_sphere_opacity( 1.0 );
core::Real const protein_sphere_opacity( 1.0 );
core::Real const ghost_sphere_opacity( 0.33 );
core::Real const ligand_sphere_shininess( 25 );
core::Real const protein_sphere_shininess( 25 );
//lin parameters for ball and stick
core::Real const protein_wireframeScale( 0.2 );
core::Real const protein_stickScale( 0.2 );
core::Real const protein_sphereScale( 0.2 );
int sphereDisplayList = 0;
core::Real const BOND_LENGTH_CUTOFF2 = 6.0*6.0;
}

#endif

#if defined GL_GRAPHICS


void processMouse(int button, int state, int x, int y) {
	using namespace graphics;

	//std::cout << "processMouse: " << glutGetWindow() << ' ' << button << ' ' << state << ' ' << x << ' ' << y << ' ' <<
	// " graphicsvars: " << click << ' ' << clicked_button << ' ' << click_x << ' ' << click_y << std::endl;

	specialKey = glutGetModifiers();

	click_x = x;
	click_y = y;
	if ( state == GLUT_DOWN ) {
		clicked_button = button;
		if ( !click ) {
			click = true;
			//      std::cout << "click button: " << button << std::endl;
		}
	}
	if ( state == GLUT_UP ) {
		clicked_button = -1;
		click = false;
	}
}


// rhiju
// Jack/Phil's rotate with mouse routine. Now with more
// intuitive rotating.
void processMouseActiveMotion(int x, int y) {
	using namespace graphics;

	//std::cout << "processMouseActiveMotion: " << glutGetWindow() << ' ' << x << ' ' << y << ' ' <<
	// " graphicsvars: " << click << ' ' << clicked_button << ' ' << click_x << ' ' << click_y << std::endl;

	static int old_x;
	static int old_y;
	if ( click ) {
		old_x = click_x;
		old_y = click_y;
		click = false;
	}

	float delta_x = old_x - x;
	float delta_y = old_y - y;

	// std::cout << "clicked_button " << clicked_button << "    specialKey " << specialKey << std::endl;

	if ( ((specialKey == GLUT_ACTIVE_SHIFT) & (clicked_button == 0)) || clicked_button == 1 ) { // Zoom in/out
		double s = exp( -1.0* (double) delta_y*0.01);
		glScalef(s,s,s);
	} else if ( (specialKey == GLUT_ACTIVE_CTRL) & (clicked_button == 0) ) { // Recontour
		GraphicsState* current_gs = gs_map_[ glutGetWindow() ];
		if ( !current_gs ) {
			std::cerr << "ignoring processKeyPress for window id " << glutGetWindow() << std::endl;
			return;
		}
		current_gs->density_sigma += delta_y*0.02;
		current_gs->density_redraw = true;
	} else if ( (specialKey == GLUT_ACTIVE_SHIFT) & (clicked_button > 0) ) { //Rotate around z-axis
		// See below for explanation of premultiplication.
		GLfloat currentrotation[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, currentrotation);
		glLoadIdentity();
		glRotatef(delta_x,0,0,1);
		glMultMatrixf(currentrotation);
	} else if ( specialKey == GLUT_ACTIVE_ALT && clicked_button == 0 ) { //Pan
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT,viewport);
		//glTranslatef( -1.0*delta_x * (_right-_left)/(viewport[2]),
		//  -1.0*delta_y * (_bottom-_top)/(viewport[3]), 0);
		//Scale factors map from screen coordinates to molecule coordinates.
		glMatrixMode(GL_MODELVIEW);
		double xscale = (graphics::window_size * 2.0)/(viewport[2]);
		double yscale = (graphics::window_size * 2.0)/(viewport[3]);

		GLfloat currentrotation[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, currentrotation);
		glLoadIdentity();
		glTranslatef( -delta_x * xscale, delta_y * yscale, 0);
		glMultMatrixf(currentrotation);
	} else { //Rotate the sucker.
		//double axis_z = 0;
		double axis_x = -delta_y;
		double axis_y = -delta_x;
		double userangle = sqrt(delta_x*delta_x + delta_y*delta_y);

		glMatrixMode(GL_MODELVIEW);
		// Standard GLUT rotation is a postmultiplication - rotation around
		// molecule's inertial frame --  and leads to
		// non-intuitive behavior.
		//glRotatef(userangle,axis_x,axis_y,0.0);

		//A premultiplication -- rotates around the axis the user actually sees.
		// A little more complicated to code; unfortunately GLUT doesn't have a one-line
		// function for it.
		GLfloat currentrotation[16];
		glGetFloatv(GL_MODELVIEW_MATRIX, currentrotation);
		glLoadIdentity();
		glRotatef(userangle,axis_x,axis_y,0.0);
		glMultMatrixf(currentrotation);
	}

	glutPostRedisplay();

	old_x = x;
	old_y = y;
}


void processKeyPress(unsigned char key, int /* x */, int /* y */) {
	using namespace graphics;
	using namespace protocols::viewer;

	GraphicsState* current_gs = gs_map_[ glutGetWindow() ];
	if ( !current_gs ) {
		std::cerr << "ignoring processKeyPress for window id " << glutGetWindow() << std::endl;
		return;
	}

	if ( key == 67 || key == 99 ) { //'c' control color
		current_gs->Color_mode = ColorMode ( current_gs->Color_mode + 1 );
		if ( current_gs->Color_mode > RESIDUE_CPK_COLOR ) current_gs->Color_mode = RAINBOW_COLOR;
	}
	if ( key == 66 || key == 98 ) { //'b' control backbone display
		current_gs->BBdisplay_state = BBdisplayState ( current_gs->BBdisplay_state + 1 );
		if ( current_gs->BBdisplay_state > SHOW_BACKBONE ) current_gs->BBdisplay_state = SHOW_NOBB;
	}
	if ( key == 72 || key == 104 ) { //'H' or 'h': toggle hydrogens
		current_gs->show_H_state = ShowHState ( current_gs->show_H_state + 1 );
		if ( current_gs->show_H_state > SHOW_H ) current_gs->show_H_state = SHOW_NO_H;
	}
	if ( key == 83 || key == 115 ) { //'s' control sidechain display
		current_gs->SCdisplay_state = SCdisplayState ( current_gs->SCdisplay_state + 1 );
		if ( current_gs->SCdisplay_state > SHOW_WIREFRAME ) current_gs->SCdisplay_state = SHOW_NOSC;
	}

	glutPostRedisplay();
}
/**/


/////////////////////////////////////////////////////////////////////////
// GRAPHICS THREAD
//
// this is the displayFunc for the windows of ConformationViewer objects
void
conformation_viewer_display( void )
{
	// which viewer?
	int const window( glutGetWindow() );

	//std::cout << "conformation_viewer_display: " << window << std::endl;

	if ( conformation_viewers.count( window ) ) {
		conformation_viewers.find( window )->second->display_func();
		glutSwapBuffers();
		//glFlush();
	}

}

/////////////////////////////////////////////////////////////////////////
//
// GRAPHICS THREAD
void
idle_func( void )
{

	check_for_new_conformation_viewers();

	// pthread_mutex_lock( &conformation_viewers_mut );

	for ( ConformationViewers::const_iterator iter = conformation_viewers.begin(), iter_end = conformation_viewers.end();
			iter != iter_end; ++iter ) {
		iter->second->display_if_necessary();
	}

	// pthread_mutex_unlock( &conformation_viewers_mut );

}


/////////////////////////////////////////////////////////////////////////
// GRAPHICS THREAD


int
conformation_viewer_window_init( GraphicsState& gs, std::string const & window_name, int length, int width )
{
	using namespace graphics;

	glutInitWindowSize (length, width);
	glutInitWindowPosition (370, 10);

	int const window = glutCreateWindow ( window_name.c_str() );

	// register gs object in gs map
	gs.BBdisplay_state = SHOW_BACKBONE;
	gs.SCdisplay_state = SHOW_WIREFRAME;
	gs.Color_mode = RAINBOW_COLOR;
	gs.Trajectory_state = SHOW_ALL_TRIALS;
	gs.show_H_state = SHOW_NO_H;

	graphics::gs_map_[window] = &gs;

	glClearColor( bg_color.x(), bg_color.y(), bg_color.z(), 1.0 );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ) ;

	glutDisplayFunc( conformation_viewer_display );

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	using graphics::window_size;
	glOrtho( -window_size, window_size, -window_size, window_size, -10*window_size, 10*window_size );
	glEnable( GL_DEPTH_TEST );
	glMatrixMode(GL_MODELVIEW);

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POINT_SMOOTH);
	glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	glutMouseFunc( processMouse );
	//glutMouseWheelFunc( processMouseWheel );
	glutMotionFunc( processMouseActiveMotion );
	glutKeyboardFunc( processKeyPress );
	//glutSpecialFunc(specialKeyPressFunc);

	runtime_assert( glutGetWindow() == window );
	return window;

} // conformation_viewer_window_init( ... )


int
conformation_viewer_window_init( GraphicsState& gs, std::string const & window_name )
{

	using namespace graphics;
	glutInitWindowSize (900, 900);
	int window = conformation_viewer_window_init( gs, window_name, 900, 900 );
	return window;
}

/////////////////////////////////////////////////////////////////////////
// GRAPHICS THREAD
//
void
check_for_new_conformation_viewers()
{
	pthread_mutex_lock( &new_conformation_viewers_mut );
	if ( !new_conformation_viewers.empty() ) {
		for ( Size i=1; i<= new_conformation_viewers.size(); ++i ) {
			ConformationViewerOP viewer( new_conformation_viewers[i] );
			// create the new window

			int width  = viewer->get_width();
			int length = viewer->get_length();

			int const new_window( conformation_viewer_window_init( viewer->get_gs(), viewer->name(), length, width ) );
			viewer->window( new_window );
			conformation_viewers[ new_window ] = viewer;
		}
		new_conformation_viewers.clear();
	}
	pthread_mutex_unlock( &new_conformation_viewers_mut );
}


/////////////////////////////////////////////////////////////////////////
//
// WORKER THREAD (conformation owner)
void
add_conformation_viewer(
	conformation::Conformation & conformation,
	std::string const & name_in, // = ""
	int const length /* = 900 */,
	int const width /* = 900 */,
	bool const debug_pause /* = false */,
	bool const set_center_vector /* = false */,
	core::Vector const center_vector /* = empty_vector*/
)
{

	pthread_mutex_lock( &new_conformation_viewers_mut );

	// create a new viewer
	std::string const window_name
		( name_in.empty() ? "conformation"+ObjexxFCL::string_of( conformation_viewers.size() + new_conformation_viewers.size() ) :
		name_in );

	ConformationViewerOP viewer( new ConformationViewer( window_name, length, width, debug_pause ) );
	if ( set_center_vector ) viewer->set_center_vector( center_vector );

	viewer->attach_to( conformation );

	new_conformation_viewers.push_back( viewer );

	pthread_mutex_unlock( &new_conformation_viewers_mut );

	// allow main to start if this is the 1st window
	pthread_cond_broadcast( &start_cond );

}


/////////////////////////////////////////////////////////////////////////
//
// WORKER THREAD (conformation owner)
//
void
add_monte_carlo_viewer(
	moves::MonteCarlo & mc,
	std::string const & name_in, //= ""
	int const length,
	int const width,
	bool debug_pause
)
{

	pthread_mutex_lock( &new_conformation_viewers_mut );

	// create a new viewer
	std::string const tag
		( name_in.empty() ?
		"MC"+ObjexxFCL::string_of( conformation_viewers.size() + new_conformation_viewers.size() ) :
		name_in );

	ConformationViewerOP
		viewer1( new ConformationViewer(tag+"_last_accepted", length, width, debug_pause) ),
		viewer2( new ConformationViewer(tag+"_best_accepted", length, width, debug_pause) );

	std::cerr << "attaching viewers!!!!" << std::endl;
	mc.attach_observer_to_last_accepted_conformation( *viewer1 );
	mc.attach_observer_to_lowest_score_conformation ( *viewer2 );

	new_conformation_viewers.push_back( viewer1 );
	new_conformation_viewers.push_back( viewer2 );

	pthread_mutex_unlock( &new_conformation_viewers_mut );

	// allow main to start if this is the 1st window
	pthread_cond_broadcast( &start_cond );
}


void
silly_window_display()
{

	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ) ;
	glutWireTeapot(0.5);

	glFlush();
	glutSwapBuffers();
	glFlush();

}


/// testing/hacking
void
silly_window_init() {


	glutInitWindowSize (500, 500);
	glutInitWindowPosition (100, 100);

	glutCreateWindow ( "silly_window" );

	glutDisplayFunc( silly_window_display );

	//glutKeyboardFunc(keyPressFunc);

	//glutSpecialFunc(specialKeyPressFunc);

	glClearColor (0.0, 0.0, 0.0, 0.0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
viewer_main( VoidFunc worker_main )
{
	// initialize attributes with default values
	pthread_attr_t stack_size_attr;
	pthread_attr_init( &stack_size_attr );

#ifdef THREAD_STACK_SIZE
	// set stack size allocated to thread
	size_t thread_stack_size;
	int error = pthread_attr_getstacksize( &stack_size_attr, &thread_stack_size );
	if ( !error && thread_stack_size < THREAD_STACK_SIZE ) {
		pthread_attr_setstacksize( &stack_size_attr, THREAD_STACK_SIZE );
	}
#endif

	// launch rosetta thread (worker)
	pthread_t p;
	pthread_create ( &p, &stack_size_attr, worker_main, NULL );

	// start glut
	int argc(1);
	char **argv = new char*[1];
	argv[0] = new char[1];
	argv[0][0] = ' ';
	glutInit( &argc, argv );

	glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH );

	//silly_window_init();

	// wait for a window to get created::
	pthread_mutex_lock( &start_mut );
	if ( new_conformation_viewers.empty() )  {
		pthread_cond_wait( &start_cond, &start_mut );
	}
	pthread_mutex_unlock( &start_mut );

	check_for_new_conformation_viewers();

	glutIdleFunc( idle_func );

	glutMainLoop();
	return 0;

}

#endif

#if defined GL_GRAPHICS || defined BOINC_GRAPHICS


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// helper function
void
glVertex3fxyz( Vector const & coord )
{
	glVertex3f((float)coord.x(), (float)coord.y(), (float)coord.z() );
}

void
set_bg_color( Vector new_bg_color ) {
	runtime_assert( new_bg_color.x() >= 0 );
	runtime_assert( new_bg_color.y() >= 0 );
	runtime_assert( new_bg_color.z() >= 0 );
	graphics::bg_color = new_bg_color;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// helper function
void
glColor3fxyz( Vector const & coord )
{
	glColor3f((float)coord.x(), (float)coord.y(), (float)coord.z() );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
Vector
atom_color_by_element( std::string const & element )
{

	if ( element == "N" ) {
		return Vector(0.0, 0.0, 1.0 ); // blue
	} else if ( element == "O" ) {
		return Vector(1.0, 0.0, 0.0 ); // red
	} else if ( element == "S" ) {
		return Vector(1.0, 1.0, 0.0 ); // yellow
	} else if ( element == "P" ) {
		return Vector(1.0, 0.5, 0.0 ); // orange
	} else if ( element == "H" ) {
		return Vector(1.0, 1.0, 1.0 ); // white
	}

	return Vector(0.5, 0.5, 0.5 ); // gray (default)

}

inline float dgaussian( float x, float mean, float sd )
{
	float sqrt_2pi = 2.50662721600161f;
	float sd_sq = sd * sd;
	float mean_diff_sq = (x - mean) * (x - mean);
	return 1 / (sd * sqrt_2pi) * std::exp( -1 * mean_diff_sq / ( 2 * sd_sq ) );
}

Vector
residue_color_by_group( core::conformation::Residue const & res, int total_residue )
{
	float sd     = 0.10;
	float max    = 0.8;
	float min    = 0.2;
	float factor = (float) res.seqpos() / (float) total_residue;
	// rescale factor according to min and max
	factor = (max - min) * factor + min;

	float red    = dgaussian( factor, 0.75, sd );
	float green  = dgaussian( factor, 0.50, sd );
	float blue   = dgaussian( factor, 0.25, sd );

	Vector color( red, green, blue );
	return color;
}

// from rosetta++ protein_graphics.cc
void rainbow_color( float frac , float & red, float & green, float & blue , bool mute_color ) {
	//float muting = .7;
	float my_color = frac;
	red = my_color;
	blue = 1.0 - my_color ;
	green = (my_color < .5) ? 2.0*my_color : 2.0-2.0*my_color;
	if ( mute_color ) {
		float saturation = sqrt(red*red + green*green + blue*blue);
		float muting = .7;
		red = muting*red/saturation;
		green = muting*green/saturation;
		blue = muting*blue/saturation;
	}
}
// chu enable color by chain
void chain_color( int const chain, float & red, float & green, float & blue ) {
	static int const num_color = 5;
	int chain_local = chain%num_color;
	switch ( chain_local ) {
	case 1 : // blue
		red = 0.0;
		green = 0.0;
		blue = 1.0;
		return;
	case 2 : // green
		red = 0.0;
		green = 1.0;
		blue = 0.0;
		return;
	case 3 : // yellow
		red = 0.7;
		green = 0.7;
		blue = 0.0;
		return;
	case 4 : // orange
		red = 0.7;
		green = 0.5;
		blue = 0.0;
		return;
	case 0 : // red
		red = 1.0;
		green = 0.0;
		blue = 0.0;
		return;
	default :
		red = 0.0;
		green = 0.0;
		blue = 1.0;
		return;
	}
}

// from rosetta++ protein_graphics.cc
void get_residue_color( float i, float & red, float & green, float & blue, bool mute_color, int total_residue ) {
	float i_local = i;
	if ( i > total_residue ) i_local = total_residue;
	rainbow_color ( float(i_local) / float(total_residue), red, green, blue, mute_color);
}

// from rosetta++ protein_graphics.cc
std::map<std::string, Vector>  get_sidechain_color_rhiju() {

	std::map<std::string, Vector> sidechain_color_rhiju;

	sidechain_color_rhiju[ "ALA" ] = Vector( 0.3, 0.3, 0.3); //gray
	sidechain_color_rhiju[ "CYS" ] = Vector( 0.7, 0.7, 0.0); //yellow:
	sidechain_color_rhiju[ "ASP" ] = Vector( 0.7, 0.0, 0.0); //red
	sidechain_color_rhiju[ "GLU" ] = Vector( 0.7, 0.0, 0.0); //red
	sidechain_color_rhiju[ "PHE" ] = Vector( 0.3, 0.3, 0.3);
	sidechain_color_rhiju[ "GLY" ] = Vector( 0.7, 0.5, 0.0); //orange; this shouldn't happen from sidechain.
	sidechain_color_rhiju[ "HIS" ] = Vector( 0.0, 0.0, 0.7); //blue
	sidechain_color_rhiju[ "ILE" ] = Vector( 0.3, 0.3, 0.3);
	sidechain_color_rhiju[ "LYS" ] = Vector( 0.0, 0.0, 0.7); //blue
	sidechain_color_rhiju[ "LEU" ] = Vector( 0.3, 0.3, 0.3);
	sidechain_color_rhiju[ "MET" ] = Vector( 0.3, 0.3, 0.3);
	sidechain_color_rhiju[ "ASN" ] = Vector( 0.0, 0.5, 0.0); //green
	sidechain_color_rhiju[ "PRO" ] = Vector( 0.3, 0.3, 0.3);
	sidechain_color_rhiju[ "GLN" ] = Vector( 0.0, 0.5, 0.0); //green
	sidechain_color_rhiju[ "ARG" ] = Vector( 0.0, 0.0, 0.7); //blue
	sidechain_color_rhiju[ "SER" ] = Vector( 0.0, 0.5, 0.0); //green
	sidechain_color_rhiju[ "THR" ] = Vector( 0.0, 0.5, 0.0); //green
	sidechain_color_rhiju[ "VAL" ] = Vector( 0.3, 0.3, 0.3);
	sidechain_color_rhiju[ "TRP" ] = Vector( 0.3, 0.3, 0.3);
	sidechain_color_rhiju[ "TYR" ] = Vector( 0.0, 0.5, 0.0); //green
	sidechain_color_rhiju[ "SEP" ] = Vector( 0.5, 0.5, 0.0); //orange
	sidechain_color_rhiju[ "GUA" ] = Vector( 0.5, 0.0, 0.0); //red [now matching EteRNA]
	sidechain_color_rhiju[ "ADE" ] = Vector( 0.5, 0.5, 0.0); //yellow
	sidechain_color_rhiju[ "CYT" ] = Vector( 0.0, 0.5, 0.0); //green
	sidechain_color_rhiju[ "THY" ] = Vector( 0.0, 0.0, 0.5); //blue [now matching EteRNA]
	sidechain_color_rhiju[ "RGU" ] = Vector( 0.5, 0.0, 0.0); //red [used to be blue, now matching EteRNA]
	sidechain_color_rhiju[ "RAD" ] = Vector( 0.5, 0.5, 0.0); //yellow
	sidechain_color_rhiju[ "RCY" ] = Vector( 0.0, 0.5, 0.0); //green
	sidechain_color_rhiju[ "URA" ] = Vector( 0.0, 0.0, 0.5); //blue  [now matching EteRNA]
	sidechain_color_rhiju[ "  G" ] = Vector( 0.5, 0.0, 0.0); //red [now matching EteRNA]
	sidechain_color_rhiju[ "  A" ] = Vector( 0.5, 0.5, 0.0); //yellow
	sidechain_color_rhiju[ "  C" ] = Vector( 0.0, 0.5, 0.0); //green
	sidechain_color_rhiju[ "  U" ] = Vector( 0.0, 0.0, 0.5); //blue [now matching EteRNA]
	sidechain_color_rhiju[ " MG" ] = Vector( 0.0, 1.0, 0.0); //bright green

	return sidechain_color_rhiju;
}

////////////////////////////////////////////////////////////////////
Vector get_atom_color(
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	int const & r,
	int const & i ) {

	float red,green,blue;
	static std::map<std::string, Vector> sidechain_color_rhiju = get_sidechain_color_rhiju();

	switch ( gs.Color_mode ) {

	case CPK_COLOR :
		return atom_color_by_element( residues[r]->atom_type(i).element());

	case RAINBOW_COLOR :
		rainbow_color( float(r)/ float(gs.nres_for_graphics), red, green, blue, true /*mute_color*/);
		return Vector(red, green, blue);

	case RESIDUE_COLOR :
		if ( sidechain_color_rhiju.find(residues[r]->name3()) != sidechain_color_rhiju.end() ) {
			return sidechain_color_rhiju[residues[r]->name3()];
		}
		return Vector( 1.0, 0.5, 0.0); //orange

	case CHAIN_COLOR :
		chain_color( residues[r]->chain(), red, green, blue );
		return Vector(red, green, blue);

	case RAINBOW_CPK_COLOR :
		if ( !residues[r]->atom_is_backbone(i)  ) { //non carbon atoms
			return atom_color_by_element( residues[r]->atom_type(i).element());
		}
		rainbow_color( float(r)/ float(gs.nres_for_graphics), red, green, blue, true /*mute_color*/);
		return Vector(red, green, blue);

	case RESIDUE_CPK_COLOR :
		if ( !residues[r]->atom_is_backbone(i)  ) { //non carbon atoms
			return atom_color_by_element( residues[r]->atom_type(i).element());
		}
		if ( sidechain_color_rhiju.find(residues[r]->name3()) != sidechain_color_rhiju.end() ) {
			return sidechain_color_rhiju[residues[r]->name3()];
		}
		return Vector( 1.0, 1.0, 1.0);

	case RHIJU_COLOR :
		if ( residues[r]->is_virtual(i) || residues[r]->is_repulsive(i) ) {
			return Vector( 1.0, 1.0, 1.0 );
		} else if ( residues[r]->atom_is_backbone(i)  ) {
			rainbow_color( float(r)/ float(gs.nres_for_graphics), red, green, blue, false /*mute_color*/);
			return Vector(red, green, blue);
		} else if ( sidechain_color_rhiju.find(residues[r]->name3()) != sidechain_color_rhiju.end() ) {
			return sidechain_color_rhiju[ residues[r]->name3() ];
		}
		break;
	case GHOST_COLOR :
		if ( residues[r]->atom_is_hydrogen(i) ) {
			return Vector( graphics::ghost_color_vect.x()+0.1 , graphics::ghost_color_vect.y()+0.1, graphics::ghost_color_vect.z()+0.1 );
		}
		return Vector(graphics::ghost_color_vect.x(), graphics::ghost_color_vect.y(), graphics::ghost_color_vect.z());
	}

	return Vector( 1.0, 1.0, 1.0);
}


////////////////////////////////////////////////////////////////////////////////////////////////
void
display_residues_wireframe(
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	Vector const & center )
{
	using namespace conformation;
	using namespace chemical;

	// could get these from some runtime-configurable options set
	//Real const bond_width( 0.1 );

	// In case the view has been rotated... set z to be the axis pointing out of the screen
	glMatrixMode(GL_MODELVIEW);
	GLfloat currentrotation[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, currentrotation);
	Vector const z( currentrotation[2], currentrotation[6], currentrotation[10] );

	Size const nres( residues.size() );

	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( *(residues[i] ) );
		Size const natoms( rsd.natoms() );

		// draw connection to previous residue
		if ( i>1 && !rsd.is_lower_terminus() && rsd.is_polymer() ) {

			Residue const & prev_rsd( *(residues[i-1]));
			if ( prev_rsd.is_polymer() && prev_rsd.n_mainchain_atoms() > 0 ) {
				int const atom1( prev_rsd.mainchain_atoms()[ prev_rsd.n_mainchain_atoms() ] );
				int const atom2( rsd.mainchain_atoms()[ 1 ] );

				Vector const color1( get_atom_color( gs, residues, i-1,   atom1 ) );
				Vector const color2( get_atom_color( gs, residues, i  , atom2 ) );

				Vector const xyz1( prev_rsd.xyz( atom1 ) - center );
				Vector const xyz2(      rsd.xyz( atom2 ) - center );

				Vector const bond( xyz2 - xyz1 );
				if ( !prev_rsd.is_virtual( atom1 ) &&
						!rsd.is_virtual( atom2 ) &&
						bond.length_squared() <= graphics::BOND_LENGTH_CUTOFF2 )  {

					Vector width( cross( bond, z ) );
					if ( width.length_squared() ) width.normalize();
					width *= graphics::protein_wireframeScale;

					// also need to draw the elbow?
					//   if (bond.length_squared() < 9.0 ){
					glColor3fxyz( color1 );
					glBegin(GL_POLYGON);
					glVertex3fxyz ( xyz1 + width );
					glVertex3fxyz ( xyz1 - width );
					glColor3fxyz( color2 );
					glVertex3fxyz ( xyz2 - width );
					glVertex3fxyz ( xyz2 + width );
					glEnd();
				}
			}
		}

		// draw the atom bonds
		utility::vector1< Vector > prev1( natoms ), prev2( natoms );
		utility::vector1< bool > prev_set( natoms, false );

		for ( Size m=1; m<= natoms; ++m ) {

			AtomIndices const & nbrs( rsd.bonded_neighbor(m) );

			for ( Size jj=1; jj<= nbrs.size(); ++jj ) {
				Size const n( nbrs[jj] );
				if ( n < m ) continue;
				//if ( rsd.atom_type(m).is_hydrogen() && graphics::exclude_hydrogens ) continue;


				Vector const color1( get_atom_color( gs, residues, i, m ) );
				Vector const color2( get_atom_color( gs, residues, i, n ) );

				Vector const xyz1( rsd.xyz(m) - center );
				Vector const xyz2( rsd.xyz(n) - center );

				Vector const bond( xyz2 - xyz1 );

				// check for chainbreaks
				//if (bond.length_squared() > graphics::BOND_LENGTH_CUTOFF2 ) break;

				Vector width( cross( bond, z ) );
				if ( width.length_squared() ) width.normalize();
				width *= graphics::protein_wireframeScale;

				if ( rsd.atom_type(m).is_hydrogen() || rsd.atom_type(n).is_hydrogen() || rsd.is_virtual(m) || rsd.is_virtual(n) ) width *= 0.5;

				glColor3fxyz( color1 );


				if ( prev_set[ m ] ) {
					// draw the elbow
					glBegin(GL_POLYGON);
					glVertex3fxyz ( prev1[m] );
					glVertex3fxyz ( xyz1 + width );
					glVertex3fxyz ( xyz1 - width );
					glVertex3fxyz ( prev2[m] );
					glEnd();
				} else {
					prev_set[m] = true;
					prev1[m] = xyz1 - width;
					prev2[m] = xyz1 + width;
				}


				glBegin(GL_POLYGON);
				glVertex3fxyz ( xyz1 + width );
				glVertex3fxyz ( xyz1 - width );
				glColor3fxyz( color2 ); // change color
				glVertex3fxyz ( xyz2 - width );
				glVertex3fxyz ( xyz2 + width );
				glEnd();

				if ( !prev_set[n] ) {
					// for drawing the elbow at atomn
					prev_set[n] = true;
					prev1[n] = xyz2 + width;
					prev2[n] = xyz2 - width;
				} else {
					// draw the elbow
					glBegin(GL_POLYGON);
					glVertex3fxyz ( prev1[n] );
					glVertex3fxyz ( xyz2 - width );
					glVertex3fxyz ( xyz2 + width );
					glVertex3fxyz ( prev2[n] );
					glEnd();
				}
			} // jj
		} // i
	} // nres
} // void display_residues_wireframe


////////////////////////////////////////////////////////////////////////////////////////////////
// Secondary structure display methods from Rhiju's code in rosetta++
// protein_graphics.cc


// placeholder function
//bool check_occupancy( int, const core::pose::Pose & ) { return true; }

void get_direction( Vector & direction, const int & next_res, const int & prior_res,
	utility::vector1< core::conformation::ResidueCOP > const & residues ) {
	direction = residues[ next_res ]->xyz( "CA" ) - residues[ prior_res ]->xyz( "CA" );
	if ( direction.length_squared() > 0.00001 ) direction.normalize();
}


void get_normal( Vector & normal, const int n, utility::vector1< core::conformation::ResidueCOP > const & residues ) {
	normal = cross( (residues[ n ]->xyz( "CA" )-residues[ n-1 ]->xyz( "CA" )),
		(residues[ n+1 ]->xyz( "CA" )-residues[ n ]->xyz( "CA" )) );
#ifdef BOINC
	if (normal.length_squared() > 0.00001) normal.normalize();
#else
	if ( normal.length_squared() > 0.00001 ) normal.normalize(15);
#endif
}


void get_axis_and_tangent( Vector & axis, Vector & tangent,
	const Vector & direction, const Vector & normal) {
	//Magic linear combinations from Per Kraulis' molscript.
	const core::Real HELIX_ALPHA = 0.5585;
	axis = std::cos(HELIX_ALPHA) * normal + std::sin(HELIX_ALPHA) * direction;
	const core::Real HELIX_BETA = -0.1920;
	tangent = std::cos(HELIX_BETA) * direction + std::sin(HELIX_BETA) * normal;
}


Vector get_CA_segment( const Vector & prev_CA, const Vector & current_CA,
	const Vector & prev_tangent, const Vector & tangent,
	const float & p, const float & hermite_factor) {
	//Hermitean interpolation.
	const core::Real p2 = p*p;
	const core::Real p3 = p*p*p;
	const core::Real h1 = 2*p3 - 3*p2 + 1;
	const core::Real h2 =-2*p3 + 3*p2;
	const core::Real h3 =   p3 - 2*p2 + p;
	const core::Real h4 =   p3 -   p2;
	return (h1*prev_CA + h2*current_CA + h3*hermite_factor*prev_tangent + h4*hermite_factor*tangent);
}


float get_half_width( const std::string & taper, const float & secstruct_half_width, const float & p ) {
	const Real PI = numeric::NumericTraits<Real>::pi();

	core::Real half_width = secstruct_half_width;
	if ( taper == "start" ) {
		half_width = graphics::COILRADIUS + (secstruct_half_width - graphics::COILRADIUS)*
			0.5 * ( - cos( PI * p )  + 1.0 );
	}
	if ( taper == "end" ) {
		half_width = graphics::COILRADIUS + (secstruct_half_width - graphics::COILRADIUS)*
			0.5 * ( cos( PI * p )  + 1.0 );
	}
	if ( taper == "strand_ultimate" ) {
		half_width = graphics::COILRADIUS + (2*secstruct_half_width - graphics::COILRADIUS) * (1 - p);
	}
	return half_width;
}


void set_initial_polygon_vertices(const Vector & vec1, const Vector & vec2, GraphicsState & gs) {
	gs.previous_vertex1 = vec1;
	gs.previous_vertex2 = vec2;
	gs.previous_width_vector = 0.0;
}


void draw_next_polygon( const Vector & vec1, const Vector & vec2,
	const float & red, const float & green, const float & blue, const int & aa,
	GraphicsState & gs,
	bool is_coil = false,
	bool darken_inside = false) {
	if ( aa<1 ) return;

	//Give the cartoon some thickness?
	const bool show_thickness = true;
	Vector width_vector;

	core::Real showbondcutoff2 = (is_coil) ? graphics::SHOWBONDCUTOFF2_COIL : graphics::SHOWBONDCUTOFF2;

	const core::Real bond_length2 = ((vec1 + vec2)/2.0 - (gs.previous_vertex1 + gs.previous_vertex2)/2.0).length_squared();
	if ( !show_thickness ) {
		glColor3f(red,green,blue);
		glBegin(GL_POLYGON);
		glVertex3fxyz ( gs.previous_vertex1 );
		glVertex3fxyz ( gs.previous_vertex2 );
		glVertex3fxyz ( vec2 );
		glVertex3fxyz ( vec1 );
		glVertex3fxyz ( gs.previous_vertex1 );
		glEnd();
	} else {
		width_vector = cross( (vec1 + vec2) - (gs.previous_vertex1 + gs.previous_vertex2),
			(vec1 + gs.previous_vertex1) - (vec2 + gs.previous_vertex2) );
		if ( width_vector.length_squared() > 0.000001 ) width_vector.normalize();
		const float cartoon_width = 0.3;
		width_vector *= -1.0f * cartoon_width;

		if ( gs.previous_width_vector == 0.0 ) gs.previous_width_vector = width_vector;

		if ( bond_length2 < showbondcutoff2 ) {
			//outside
			glColor3f(red,green,blue);
			glBegin(GL_POLYGON);
			glVertex3fxyz ( gs.previous_vertex1 + gs.previous_width_vector );
			glVertex3fxyz ( gs.previous_vertex2 + gs.previous_width_vector );
			glVertex3fxyz ( vec2 + width_vector );
			glVertex3fxyz ( vec1 + width_vector );
			glEnd();

			//inside
			if ( darken_inside ) glColor3f(0.5*red,0.5*green,0.5*blue);
			glBegin(GL_POLYGON);
			glVertex3fxyz ( gs.previous_vertex1 - gs.previous_width_vector );
			glVertex3fxyz ( gs.previous_vertex2 - gs.previous_width_vector );
			glVertex3fxyz ( vec2 - width_vector );
			glVertex3fxyz ( vec1 - width_vector );
			glEnd();

			//edges
			glColor3f(0.7*red,0.7*green,0.7*blue);
			glBegin(GL_POLYGON);
			glVertex3fxyz ( gs.previous_vertex1 + gs.previous_width_vector );
			glVertex3fxyz ( vec1 + width_vector );
			glVertex3fxyz ( vec1 - width_vector );
			glVertex3fxyz ( gs.previous_vertex1 - gs.previous_width_vector );
			glEnd();
			glBegin(GL_POLYGON);
			glVertex3fxyz ( gs.previous_vertex2 + gs.previous_width_vector );
			glVertex3fxyz ( vec2 + width_vector );
			glVertex3fxyz ( vec2 - width_vector );
			glVertex3fxyz ( gs.previous_vertex2 - gs.previous_width_vector );
			glEnd();
		}

	}

	gs.previous_vertex1 = vec1;
	gs.previous_vertex2 = vec2;
	gs.previous_width_vector = width_vector;
}


void draw_secstruct_chunk(
	const Vector & prev_CA,
	const Vector & current_CA,
	const Vector & prev_tangent,
	const Vector & tangent,
	const Vector & prev_axis,
	const Vector & axis,
	const int &  n,
	const char & secstruct_res,
	const std::string & taper,
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues )
{
	float hermite_factor;
	float red( 0.0 ), green( 0.0 ), blue( 0.0 );
	float secstruct_half_width;
	bool darken_inside;
	if ( secstruct_res == 'H' ) {
		secstruct_half_width = graphics::HELIX_HALF_WIDTH;
		darken_inside = true;
		hermite_factor = graphics::HELIX_HERMITE_FACTOR;
	} else if ( secstruct_res == 'E' ) {
		secstruct_half_width = graphics::STRAND_HALF_WIDTH;
		darken_inside = false;
		hermite_factor = graphics::STRAND_HERMITE_FACTOR;
	} else {
		runtime_assert( secstruct_res == 'N' );
		secstruct_half_width = graphics::STRAND_HALF_WIDTH;
		darken_inside = true;
		hermite_factor = graphics::NA_HERMITE_FACTOR;
	}
	Vector axis_segment, CA_segment;
	for ( int s = 1; s <= graphics::NUM_SEGMENTS; s++ ) {
		const core::Real p = s / static_cast<core::Real>(graphics::NUM_SEGMENTS);
		axis_segment = p*axis + (1-p)*prev_axis;
		CA_segment = get_CA_segment( prev_CA, current_CA, prev_tangent, tangent, p, hermite_factor);
		if ( gs.Color_mode == CHAIN_COLOR ) {
			chain_color( residues[n]->chain(), red, green, blue );
		} else if ( gs.Color_mode == GHOST_COLOR ) {
			red=graphics::ghost_color_vect.x();
			green=graphics::ghost_color_vect.y();
			blue=graphics::ghost_color_vect.z();
		} else  {
			get_residue_color( static_cast<float>(n - 1) + p, red, green, blue, false, gs.nres_for_graphics );
		}
		const core::Real half_width = get_half_width( taper, secstruct_half_width, p );
		draw_next_polygon(CA_segment - half_width*axis_segment, CA_segment + half_width*axis_segment,
			red,green,blue, n, gs, false, darken_inside);
	}
}


void draw_helix(
	const int & start,
	const int & end,
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues )
{

	const int total_residue = residues.size();

	//Starting point.
	int prior_res( start-1 );
	if ( prior_res < 1 ) prior_res = 1;
	Vector direction, normal, tangent, axis, current_CA;
	Vector prev_CA, prev_tangent, prev_axis;

	//For the starting point, need to figure out axis from next residue...
	int next_res = start + 1;
	if ( next_res > total_residue-1 ) next_res = total_residue-1;
	get_direction( direction, next_res+1, next_res-1, residues);
	get_normal( normal, next_res, residues );
	get_axis_and_tangent( axis, tangent, direction, normal );

	get_direction(direction, start+1, prior_res, residues);
	tangent = direction;
	current_CA = residues[ start ]->xyz( "CA" );
	if ( start==1 ) {
		set_initial_polygon_vertices( current_CA - graphics::HELIX_HALF_WIDTH*axis,
			current_CA + graphics::HELIX_HALF_WIDTH*axis, gs);
	}

	//previous residue's helix geometry
	prev_CA = current_CA;
	prev_tangent = tangent;
	prev_axis = axis;
	std::string taper = "start";

	// Draw the body of the helix.
	for ( int n = start+1; n<=end-1; n++ ) {
		//new residue's helix geometry
		get_direction( direction, n+1, n-1, residues );
		get_normal( normal, n, residues );
		get_axis_and_tangent( axis, tangent, direction, normal);
		current_CA = residues[ n ]->xyz( "CA" );
		draw_secstruct_chunk( prev_CA, current_CA, prev_tangent, tangent, prev_axis, axis, n,
			'H', taper, gs, residues);

		//previous residue's helix geometry
		prev_CA = current_CA;
		prev_tangent = tangent;
		prev_axis = axis;
		taper = "none";
	}

	//last piece.
	//  next_res = end + 1;
	//  if (next_res > nres)  next_res = nres;
	//  get_direction( direction, next_res, end - 1, xyz_full);
	get_direction( direction, end, end - 1, residues);
	tangent = direction;
	current_CA = residues[ end ]->xyz( "CA" );
	taper = "end";
	draw_secstruct_chunk( prev_CA, current_CA, prev_tangent, tangent, prev_axis, axis, end,
		'H', taper, gs, residues);
}


void draw_strand(
	const int & start,
	const int & end,
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues ) {

	const int total_residue = residues.size();

	//Pre-smooth?  priestle_smooth
	//core::pose::Pose smooth_pose( pose );

	//Starting point.
	int prior_res( start-1 );
	if ( prior_res < 1 ) prior_res = 1;
	Vector direction, normal, tangent, axis, current_CA, prev_CA, prev_direction, prev_normal;

	//For the starting point, need to figure out axis from next residue...
	int next_res = start + 1;
	if ( next_res > total_residue-1 ) next_res = total_residue-1;
	get_direction( direction, next_res+1, next_res-1, residues);
	get_normal( normal, next_res, residues );

	current_CA = residues[ start ]->xyz( "CA" );
	if ( start==1 ) {
		set_initial_polygon_vertices( current_CA - graphics::STRAND_HALF_WIDTH*normal,
			current_CA + graphics::STRAND_HALF_WIDTH*normal, gs);
	}

	//previous residue's strand geometry
	prev_CA = current_CA;
	prev_normal = normal;
	prev_direction = direction;
	std::string taper = "start";

	// Draw the body of the strand.
	for ( int n = start+1; n<=end-1; n++ ) {
		//new residue's strand geometry
		get_direction( direction, n+1, n-1, residues);
		get_normal( normal, n, residues );
		if ( dot(normal, prev_normal) < 0.0 ) normal *= -1.0;
		current_CA = residues[ n ]->xyz( "CA" );
		draw_secstruct_chunk( prev_CA, current_CA, prev_direction, direction, prev_normal, normal, n,
			'E', taper, gs, residues);

		//previous residue's strand geometry
		prev_CA = current_CA;
		prev_normal = normal;
		prev_direction = direction;
		taper = "none";
	}

	//last piece.
	next_res = end + 1;
	if ( next_res > total_residue )  next_res = total_residue;
	get_direction( direction, next_res, end - 1, residues);
	tangent = direction;
	current_CA = residues[ end ]->xyz( "CA" );
	taper = "strand_ultimate";
	draw_secstruct_chunk( prev_CA, current_CA, prev_direction, direction, prev_normal, normal, end,
		'E', taper, gs, residues);
}


void draw_coil_chunk(
	const Vector & prev_CA,
	const Vector & current_CA,
	const Vector & prev_tangent,
	const Vector & tangent,
	const int & n,
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues ) {
	float red( 0.0 ), green( 0.0 ), blue( 0.0 );

	Vector  axis_segment, CA_segment, prev_CA_segment, bond, prev_bond;

	GLfloat currentrotation[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, currentrotation);
	Vector z(currentrotation[2],currentrotation[6],currentrotation[10]); // z pointing out of window in current view.

	prev_CA_segment = prev_CA;
	prev_bond = current_CA - prev_CA;
	for ( int s = 1; s <= graphics::NUM_SEGMENTS_COIL; s++ ) {
		const float p = s / static_cast<float>(graphics::NUM_SEGMENTS_COIL);
		CA_segment = get_CA_segment( prev_CA, current_CA, prev_tangent, tangent, p, graphics::COIL_HERMITE_FACTOR);
		if ( gs.Color_mode == CHAIN_COLOR ) {
			chain_color( residues[n]->chain(), red, green, blue );
		} else if ( gs.Color_mode == GHOST_COLOR ) {
			red=graphics::ghost_color_vect.x();
			green=graphics::ghost_color_vect.y();
			blue=graphics::ghost_color_vect.z();
		} else  {
			get_residue_color ( static_cast<float>(n) + p, red, green, blue, false, gs.nres_for_graphics );
		}
		//Need to replace the following with a real cylinder.
		bond = CA_segment - prev_CA_segment;
		Vector width =  cross( bond, z );
		if ( width.length_squared() > 0.0001 )  width.normalize();
		width = (width - 0.5*z).normalized(); //nice shadow effect
		width *= graphics::COILRADIUS;

		draw_next_polygon( CA_segment + width, CA_segment - width, red, green, blue, n, gs, false );

		prev_CA_segment = CA_segment;
		prev_bond = bond;
	}
}


void draw_coil(
	const int & start,
	const int & end,
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues ) {
	const int total_residue = residues.size();
	// Draw the body of the coil.

	Vector prev_CA, prev_direction, direction, current_CA;
	prev_CA = residues[ start ]->xyz( "CA" );
	get_direction( prev_direction, start+1, start, residues );
	if ( start==1 ) {
		set_initial_polygon_vertices( prev_CA, prev_CA, gs );
	}
	for ( int n = start; n <= end-1; n++ ) {
		if ( n < end-1 && end < total_residue ) {
			get_direction( direction, n+2, n, residues);
		} else {
			get_direction( direction, n+1, n, residues);
		}
		current_CA = residues[ n+1 ]->xyz( "CA" );
		draw_coil_chunk( prev_CA, current_CA, prev_direction, direction, n, gs, residues);
		//previous residue's coil geometry
		prev_CA = current_CA;
		prev_direction = direction;
	}
}


void draw_segment(
	const int & start_segment,
	const int & end_segment,
	const char & prev_secstruct,
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues ) {
	const int size_segment = end_segment - start_segment + 1;
	if ( start_segment >= end_segment ) return; //not drawable
	if ( prev_secstruct=='H' && size_segment >= 4 ) {
		draw_helix( start_segment, end_segment, gs, residues );
	} else if ( prev_secstruct=='E' && size_segment >= 2 ) {
		draw_strand( start_segment, end_segment, gs, residues );
	} else {
		draw_coil( start_segment, end_segment, gs, residues );
	}
}


bool check_chainbreak(const int & i, utility::vector1< core::conformation::ResidueCOP > const & residues) {
	if ( i==1 ) return false;

	float chainbreak_cutoff2 = graphics::CHAINBREAK_CUTOFF2;
	Vector vec = residues[i]->xyz("CA");
	Vector vec_prev = residues[i-1]->xyz("CA");
	const float dist2 = (vec-vec_prev).length_squared();
	if ( dist2 > chainbreak_cutoff2 ) {
		return true;
	}
	return false;
}


void draw_secstruct(
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	utility::vector1< char > const & ss,
	int const begin, int const end)
{

	char prev_secstruct = ss[ begin ];

	int start_segment( begin );
	//bool is_chainbreak = false;

	int protein_end = 0;
	for ( int i = begin+1; i<= end; i++ ) {
		if ( ! (residues[i-1]->is_protein() && residues[i]->is_protein()) ) continue;
		char current_secstruct = ss[ i ];
		protein_end = i;
		//if (current_secstruct != prev_secstruct || !check_occupancy(i, pose) || check_chainbreak(i, pose)) {
		if ( current_secstruct != prev_secstruct || check_chainbreak(i, residues) ) {
			int const end_segment = i - 1;
			draw_segment( start_segment, end_segment, prev_secstruct, gs, residues );
			start_segment = i;

			// if (!check_occupancy(i, pose)) {
			//  start_segment++;
			// } else
			if ( !check_chainbreak(i, residues) ) {
				draw_segment( end_segment, end_segment+1, 'L', gs,  residues ); // connector region
			}
		}
		prev_secstruct = current_secstruct;
	}
	//if (check_occupancy(end, pose))
	draw_segment( start_segment, protein_end, prev_secstruct, gs, residues );
}


void draw_Calpha_trace(
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	const int & start, const int & end, float xwidth = 0.5)
{

	GLfloat currentrotation[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, currentrotation);
	Vector z(currentrotation[2],currentrotation[6],currentrotation[10]); // z pointing out of window in current view.
	const core::Real z_offset = .1;
	const Vector z_halo = z * z_offset;

	Vector prev1( 0.0 ), prev2( 0.0 );
	bool last_bonded = false;
	for ( int i = start; i<end ; i++ ) {
		if ( ! (residues[i]->is_protein() && residues[i+1]->is_protein()) ) continue;
		Vector ca_pos1, ca_pos2;
		Vector ca_pos1tmp = residues[i]->xyz("CA");
		Vector ca_pos2tmp = residues[i+1]->xyz("CA");
		ca_pos1 = ca_pos1tmp;
		ca_pos2 = ca_pos2tmp;

		Vector bond;
		bond = ca_pos2 - ca_pos1;
		Vector width( cross( bond, z ));
		if ( width.length_squared() > 0.0001 )  width.normalize();
		width = width * (core::Real)xwidth;

		float red, green, blue;
		if ( gs.Color_mode == CHAIN_COLOR ) {
			chain_color( residues[i]->chain(), red, green, blue );
		} else if ( gs.Color_mode == GHOST_COLOR ) {
			red=graphics::ghost_color_vect.x();
			green=graphics::ghost_color_vect.y();
			blue=graphics::ghost_color_vect.z();
		} else {
			get_residue_color( i, red, green, blue, false,  gs.nres_for_graphics );
		}
		if ( i > 1 && last_bonded ) {
			glColor3f(red, green, blue);
			glBegin(GL_POLYGON);
			glVertex3fxyz ( prev1 );
			glVertex3fxyz ( ca_pos1 + width );
			glVertex3fxyz ( ca_pos1 - width );
			glVertex3fxyz ( prev2 );
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINES);
			glVertex3fxyz ( prev1 );
			glVertex3fxyz ( ca_pos1 + width );
			glVertex3fxyz ( ca_pos1 - width );
			glVertex3fxyz ( prev2 );
			glEnd();
		}
		last_bonded = false;
		if ( bond.length_squared() < 16 ) {
			last_bonded = true;
			glColor3f(red, green, blue);
			glBegin(GL_POLYGON);
			glVertex3fxyz ( ca_pos1 + width );
			glVertex3fxyz ( ca_pos2 + width );
			glVertex3fxyz ( ca_pos2 - width );
			glVertex3fxyz ( ca_pos1 - width );
			glEnd();

			glColor3f(0,0,0);
			glBegin(GL_LINES);
			glVertex3fxyz ( z_halo + ca_pos1 + width );
			glVertex3fxyz ( z_halo + ca_pos2 + width );
			glVertex3fxyz ( z_halo + ca_pos2 - width );
			glVertex3fxyz ( z_halo + ca_pos1 - width );
			glEnd();
		}
		prev1 = ca_pos2 + width;
		prev2 = ca_pos2 - width;
	}

}


void
draw_sidechains( GraphicsState & gs, utility::vector1< core::conformation::ResidueCOP > const & residues, const int & start, const int & end ) {

	if ( gs.SCdisplay_state == SHOW_NOSC ) return;
	else if ( gs.SCdisplay_state == SHOW_SCSPHERES ) {
		draw_sphere( gs, residues, SPHERE_MODE_SC );
		return;
	}

	GLfloat currentrotation[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, currentrotation);
	Vector z(currentrotation[2],currentrotation[6],currentrotation[10]); // z pointing out of window in current view.

	float xwidth = graphics::protein_stickScale;
	for ( int r = start; r<=end ; ++r ) {

		if ( residues[r]->is_ligand() && gs.SCdisplay_state != SHOW_STICK
				&& gs.SCdisplay_state != SHOW_WIREFRAME ) continue;

		if ( !residues[r]->is_ligand() && gs.SCdisplay_state == SHOW_WIREFRAME ) {
			xwidth = graphics::protein_wireframeScale;
		}

		if ( residues[r]->residue_type_set()->name() == chemical::CENTROID
				&& residues[r]->aa() != core::chemical::aa_unk ) continue;

		Size const natoms( residues[r]->natoms() );

		// draw the atom bonds
		utility::vector1< Vector > prev1( natoms ), prev2( natoms );
		utility::vector1< bool > prev_set( natoms, false );

		for ( Size i=1; i<= natoms; ++i ) {
			core::chemical::AtomIndices const & nbrs( residues[r]->bonded_neighbor(i) );

			for ( Size jj=1; jj<= nbrs.size(); ++jj ) {
				Size const j( nbrs[jj] );
				if ( j < i ) continue;
				if ( residues[r]->is_virtual(j) ) continue; //no virtual atoms
				if ( residues[r]->atom_type(j).is_hydrogen() && gs.show_H_state == SHOW_NO_H ) continue;
				//if ( pose.residue(r).atom_is_backbone(j) && pose.residue(r).atom_name(j) != "CA" ) continue;

				Vector const color1 = get_atom_color( gs, residues, r, i );
				Vector const color2 = get_atom_color( gs, residues, r, j );

				Vector const xyz1( residues[r]->xyz(i));
				Vector const xyz2( residues[r]->xyz(j));

				Vector const bond( xyz2 - xyz1 );

				if ( bond.length_squared() > graphics::BOND_LENGTH_CUTOFF2 ) break;

				Vector width( cross( bond, z ) );
				if ( width.length_squared() ) width.normalize();
				width *= xwidth;

				if ( residues[r]->atom_type(j).is_hydrogen() ) width *= 0.5;

				glColor3fxyz( color1 );

				if ( prev_set[ i ] ) {
					// draw the elbow
					glBegin(GL_POLYGON);
					glVertex3fxyz ( prev1[i] );
					glVertex3fxyz ( xyz1 + width );
					glVertex3fxyz ( xyz1 - width );
					glVertex3fxyz ( prev2[i] );
					glEnd();
				} else {
					prev_set[i] = true;
					prev1[i] = xyz1 - width;
					prev2[i] = xyz1 + width;
				}

				glBegin(GL_POLYGON);
				glVertex3fxyz ( xyz1 + width );
				glVertex3fxyz ( xyz1 - width );
				glColor3fxyz( color2 ); // change color
				glVertex3fxyz ( xyz2 - width );
				glVertex3fxyz ( xyz2 + width );
				glEnd();

				if ( !prev_set[j] ) {
					// for drawing the elbow at atomj
					prev_set[j] = true;
					prev1[j] = xyz2 + width;
					prev2[j] = xyz2 - width;
				} else {
					// draw the elbow
					glBegin(GL_POLYGON);
					glVertex3fxyz ( prev1[j] );
					glVertex3fxyz ( xyz2 - width );
					glVertex3fxyz ( xyz2 + width );
					glVertex3fxyz ( prev2[j] );
					glEnd();
				}
			} // jj
		} // i
	} // nres

	// detect disulfides...
	Real const DISULFIDE_LENGTH_CUTOFF2( 3.0 * 3.0 );
	for ( int m = start; m <= end; m++ ) {
		if ( residues[ m ]->aa() != core::chemical::aa_cys || !residues[ m ]->type().has( "SG" ) || residues[ m ]->type().atom_type_set().name() != "fa_standard" ) continue;
		Size const i = residues[ m ]->atom_index( " SG " );
		Vector const & xyz1 = residues[ m ]->xyz( i );

		for ( int n = start; n <= end; n++ ) {
			if ( residues[ n ]->aa() != core::chemical::aa_cys || !residues[ n ]->type().has( "SG" ) || residues[ n ]->type().atom_type_set().name() != "fa_standard" ) continue;
			Size const j = residues[ n ]->atom_index( " SG " );
			Vector const & xyz2 = residues[ n ]->xyz( j );

			Vector const bond( xyz2 - xyz1 );
			if ( bond.length_squared() > DISULFIDE_LENGTH_CUTOFF2 ) continue;

			Vector width( cross( bond, z ) );
			if ( width.length_squared() ) width.normalize();
			width *= xwidth;

			Vector const color1 = get_atom_color( gs, residues, m, i );
			Vector const color2 = get_atom_color( gs, residues, n, j );
			glColor3fxyz( color1 );
			glBegin(GL_POLYGON);
			glVertex3fxyz ( xyz1 + width );
			glVertex3fxyz ( xyz1 - width );
			glColor3fxyz( color2 ); // change color
			glVertex3fxyz ( xyz2 - width );
			glVertex3fxyz ( xyz2 + width );
			glEnd();
		} // n
	} // m

} // void draw_sidechains


void draw_backbone(
	GraphicsState & gs,
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	utility::vector1< char > const & ss ) {

	switch( gs.BBdisplay_state ) {
	case SHOW_BACKBONE :
		draw_Calpha_trace( gs, residues, 1, residues.size() );
		return;
	case SHOW_CARTOON :
		draw_secstruct( gs, residues, ss, 1, residues.size() );
		return;
	case SHOW_NOBB :
		return;
	case SHOW_BBSPHERES :
		draw_sphere( gs, residues, SPHERE_MODE_BB );
		return;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//chu copied from rosetta++ protein_graphics.cc
void
draw_sphere( GraphicsState & gs, utility::vector1< core::conformation::ResidueCOP > const & residues, spheremode const &sphere_mode )
{
	using namespace graphics;
	//lin currently only ligand
	//lin should be easy to extend to non ligand
	const int nres = residues.size();
	glPushAttrib(GL_ENABLE_BIT);
	glEnable(GL_LIGHTING);

	GLUquadricObj *sphereObj = gluNewQuadric();
	// Create a glu quadric object
	gluQuadricDrawStyle(sphereObj, GLU_FILL);
	gluQuadricNormals(sphereObj, GLU_SMOOTH);

	for ( int i = 1; i <= nres; i++ ) {
		conformation::Residue const & rsd = *(residues[i]);

		if ( sphere_mode == SPHERE_MODE_LIGAND ) {
			if ( rsd.aa() != chemical::aa_unk || rsd.is_protein() ) continue;  // ligand residues only
		} else {
			if ( rsd.aa() == chemical::aa_unk && !rsd.is_protein() ) continue; // NOT ligand residues
		}

		float const sphere_opacity ( gs.Color_mode == GHOST_COLOR? ghost_sphere_opacity : (rsd.is_protein() ? protein_sphere_opacity : ligand_sphere_opacity) );
		float const sphere_shininess ( rsd.is_protein() ? static_cast<float>(protein_sphere_shininess) : static_cast<float>(ligand_sphere_shininess) );

		// loop through each heavy atom
		int atom_begin = 1 ;
		int atom_end = rsd.natoms();

		for ( int j = atom_begin; j <= atom_end; ++j ) {
			conformation::Atom const & atom( rsd.atom(j) );

			if ( sphere_mode == SPHERE_MODE_BB && !rsd.atom_is_backbone(j) ) continue; //Skip if atom is not backbone if we're doing backbone.
			if ( sphere_mode == SPHERE_MODE_SC && rsd.atom_is_backbone(j) ) continue; //Skip if atom is backbone if we're NOT doing backbone.

			if ( rsd.is_virtual(j) ) continue; //no virtual atoms

			Vector const xyz( rsd.xyz(j) );

			// Figure out the material based on the atom type.
			Vector const atom_color ( get_atom_color( gs, residues, i, j ) );

			//GLfloat mat_shininess[] = { sphere_shininess };
			GLfloat atom_material[4] = {
				static_cast<GLfloat>(atom_color[0])*sphere_opacity,
				static_cast<GLfloat>(atom_color[1])*sphere_opacity,
				static_cast<GLfloat>(atom_color[2])*sphere_opacity,
				sphere_opacity,
				};
			GLfloat specular_material[4] = {
				static_cast<GLfloat>(atom_specular_color.x()),
				static_cast<GLfloat>(atom_specular_color.y()),
				static_cast<GLfloat>(atom_specular_color.z()),
				1.0f,
				};

			// Highlight the nonprotein
			if ( ! rsd.is_protein() ) {
				for ( int color = 0 ; color < 3 ; color++ ) {
					atom_material[color] *= 1.5;
					if ( atom_material[color] > 1.0 ) {
						atom_material[color] = 1.0;
					}
				}
				if ( atom_material[3] < 0.1 ) {
					atom_material[3] = 0.1;
				}
			}

			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,  atom_material);
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  atom_material);
			glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specular_material);
			glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &sphere_shininess);

			glPushMatrix();
			glTranslatef(xyz(1), xyz(2), xyz(3));

			// draw a sphere using the sphere display list
			if ( graphics::sphereDisplayList == 0 ) {
				const int SPHERE_SLICES(16), SPHERE_STACKS(16);
				float sphereRadius(1.0);
				graphics::sphereDisplayList = glGenLists(1);
				runtime_assert(graphics::sphereDisplayList != 0);

				glNewList(graphics::sphereDisplayList, GL_COMPILE);
				gluSphere(sphereObj, sphereRadius, SPHERE_SLICES, SPHERE_STACKS);
				glEndList();
			}
			if ( graphics::sphereDisplayList != 0 ) {
				//float sphereRadius, sphereScale(1.0);
				const float scale_for_display_list(1.0);
				float sphereScale = 0.8 * scale_for_display_list ;
				float sphereRadius = sphereScale *  rsd.atom_type_set()[ atom.type() ].lj_radius();
				glScalef(sphereRadius, sphereRadius, sphereRadius);
			}
			glCallList(graphics::sphereDisplayList);

			glPopMatrix();
		}
	}
	// clean up
	glPopAttrib();
	gluDeleteQuadric(sphereObj);
	glDisable(GL_LIGHTING);// Turn lighting off
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////
////// electron density display functions
//////
void render_density(
	GraphicsState &gs,
	utility::vector1< triangle > &triangles) {
	// "on demand"
	if ( !gs.density_redraw ) return;

	triangles.clear();

	const core::scoring::electron_density::ElectronDensity& edm = core::scoring::electron_density::getDensityMap();
	const float thresh = edm.getMean() + gs.density_sigma*edm.getStdev();
	triangleIterator tri_it( edm.get_data(), thresh );

	numeric::xyzVector_float vertex[3], normal[3];
	numeric::xyzVector< core::Real> cart_vi[3], cart_ni[3];
	while ( tri_it.hasNext() ) {
		tri_it.next( vertex, normal );
		for ( int j=0; j<3; ++j ) {
			edm.idx2cart( vertex[j] , cart_vi[j] );
			edm.idxoffset2cart( normal[j] , cart_ni[j] );
			cart_ni[j].normalize();

			//vertices.push_back( cart_vi );
			//normals.push_back( cart_ni );
		}
		triangles.push_back( triangle( cart_vi, cart_ni ) );
	}
}

// for sorting
bool operator < (const triangle&  left, const triangle&  right) {
	return (left.depth_ < right.depth_);
}


void
display_density(
	GraphicsState & /*gs*/,
	utility::vector1< triangle > &triangles ) {
	using namespace conformation;
	using namespace chemical;

	if ( triangles.size() == 0 ) return;
	//glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );  // wireframe

	glMatrixMode(GL_MODELVIEW);
	GLfloat currentrotation[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, currentrotation);
	numeric::xyzVector_float const z( currentrotation[2], currentrotation[6], currentrotation[10] );

	glPushAttrib(GL_ENABLE_BIT);
	glEnable(GL_BLEND);                    //activate blending mode
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);    //define blending factors
	glEnable(GL_LIGHTING);
	glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE ) ;
	glEnable ( GL_COLOR_MATERIAL ) ;
	glDisable( GL_DEPTH_TEST );

	// sort triangles back to front
	core::Size const ntris = triangles.size();
	for ( core::Size i=1; i<=ntris; ++i ) {
		triangles[i].update( z );
	}

	std::sort( triangles.begin(), triangles.end() );

	glBegin(GL_TRIANGLES);
	Vector const color( 0.8,0.8,0.8 );
	glColor4f( color[0], color[1], color[2], 0.3 );

	for ( Size i=1; i<=ntris; ++i ) {
		for ( Size j=0; j<3; ++j ) {
			glNormal3f ( triangles[i].normals_[j][0] , triangles[i].normals_[j][1] , triangles[i].normals_[j][2] );
			glVertex3f ( triangles[i].vertices_[j][0], triangles[i].vertices_[j][1], triangles[i].vertices_[j][2] );
		}
	}
	glEnd();

	glPopAttrib();
	glDisable(GL_LIGHTING);
}

void
draw_conformation_and_density(
	utility::vector1< core::conformation::ResidueCOP > const & residues,
	utility::vector1< char > const & ss,
	utility::vector1< triangle > &triangles,
	GraphicsState & gs,
	Vector const & center )
{

	using namespace graphics;
	const int total_residue = residues.size();

	//Set background color
	glClearColor( bg_color.x(), bg_color.y(), bg_color.z(), 1.0 );

	// clear
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	// In case the view has been rotated... set z to be the axis pointing out of the screen
	glMatrixMode(GL_MODELVIEW);
	GLfloat currentrotation[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, currentrotation);

	Vector const x( currentrotation[0], currentrotation[4], currentrotation[ 8] );
	Vector const y( currentrotation[1], currentrotation[5], currentrotation[ 9] );
	Vector const z( currentrotation[2], currentrotation[6], currentrotation[10] );
	Vector light = z + y ;
	Vector light2 = z - y - x;

	GLfloat light_position[] = { static_cast<GLfloat>(light[0]), static_cast<GLfloat>(light[1]), static_cast<GLfloat>(light[2]), 0.0 };
	GLfloat light2_position[] = { static_cast<GLfloat>(light2[0]), static_cast<GLfloat>(light2[1]), static_cast<GLfloat>(light2[2]), 0.0 };
	//GLfloat light_ambient[]  = { 0.0, 0.0, 0.0, 1.0 };
	//GLfloat light_diffuse[]  = { 1.0, 1.0, 1.0, 1.0 };
	//GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
	//GLfloat mat_emission[]   = { 0.0, 0.0, 0.0, 1.0 };
	//GLfloat mat_specular[]   = { 1.0, 1.0, 1.0, 1.0 };

	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);

	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT1, GL_POSITION, light2_position);
	//glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	//glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	//glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	//glMaterialfv( GL_FRONT_AND_BACK, GL_EMISSION, mat_emission );
	//glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular );


	if ( total_residue > 0 ) {
		glPushMatrix();
		glTranslatef(-center.x(), -center.y(), -center.z());

		draw_backbone( gs, residues, ss );
		draw_sidechains( gs, residues, 1, total_residue );
		draw_sphere( gs, residues, SPHERE_MODE_LIGAND );

		render_density( gs, triangles );  // renders "on-demand" depending on gs object
		display_density( gs, triangles );
		glPopMatrix();
	}

	//glDisable(GL_LIGHT0);// Turn lighting off
	//glDisable(GL_LIGHT1);// Turn lighting off
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
fill_nres_for_graphics( GraphicsState & gs, utility::vector1< conformation::ResidueCOP > const & residues ) {
	Size nres( 0 );
	for ( Size n = 1; n <= residues.size(); n++ ) {
		if ( residues[ n ]->is_protein() || residues[ n ]->is_NA() ) {
			nres++;
		}
	}
	if ( nres == 0 ) nres = 1;
	gs.nres_for_graphics = nres;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void draw_conformation(
	utility::vector1< conformation::ResidueCOP > const & residues,
	utility::vector1< char > const & ss,
	GraphicsState & gs,
	Vector const & center
) {

	using namespace graphics;
	//const int total_residue = residues.size();

#ifndef BOINC_GRAPHICS

	//Set background color
	glClearColor( bg_color.x(), bg_color.y(), bg_color.z(), 1.0 );

	// clear
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

#endif

	GLfloat light_position[] = { 100.0, 100.0, 100.0, 0.0 };
	GLfloat light_color[] = {1.0f, 1.0f, 1.0f, 1.0f };
	GLfloat light_position2[] = { -100.0, -100.0, -100.0, 0.0 };
	GLfloat light_color2[] = {0.3f, 0.35f, 0.5f, 1.0f };
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_color);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_color);
	glLightfv(GL_LIGHT1, GL_POSITION, light_position2);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_color2);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_color2);

	glPushMatrix();
	glTranslatef(-center.x(), -center.y(), -center.z());

	utility::vector1< conformation::ResidueCOP > residues_protein, other_residues, residues_sphere;
	for ( Size n = 1; n <= residues.size(); n++ ) {
		conformation::ResidueCOP rsd = residues[ n ];
		if ( rsd->is_protein() ) {
			residues_protein.push_back( rsd );
		} else if ( rsd->is_metal() ) {
			residues_sphere.push_back( rsd );
		} else {
			other_residues.push_back( rsd );
		}
	}

	fill_nres_for_graphics( gs, residues );
	if  ( residues_protein.size() > 0 ) {
		draw_backbone( gs, residues_protein, ss );
		draw_sidechains( gs, residues_protein, 1, residues_protein.size() );
		draw_sphere( gs, residues_protein, SPHERE_MODE_LIGAND );
	}

	if ( other_residues.size() > 0 ) {
		ColorMode colormode_save = gs.Color_mode;
		gs.Color_mode = RHIJU_COLOR;
		display_residues_wireframe( gs, other_residues, Vector( 0.0, 0.0, 0.0) );
		gs.Color_mode = colormode_save;
	}

	if ( residues_sphere.size() > 0 ) {
		ColorMode colormode_save = gs.Color_mode;
		gs.Color_mode = RESIDUE_COLOR;
		draw_sphere( gs, residues_sphere, SPHERE_MODE_LIGAND );
		gs.Color_mode = colormode_save;
	}

	glPopMatrix();
	glDisable(GL_LIGHT0);// Turn lighting off
	glDisable(GL_LIGHT1);// Turn lighting off
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Vector
get_center( utility::vector1< conformation::ResidueCOP > const & residues ){

	core::Size nres = residues.size();
	Vector center_of_mass( 0, 0, 0 );

	for ( int i = 1; i <= (int) nres; ++i ) {
		conformation::Residue const & rsd = *residues[i];
		if ( !rsd.is_virtual( rsd.nbr_atom() ) ) center_of_mass += rsd.nbr_atom_xyz();
	}
	if ( nres > 0 ) center_of_mass /= nres;

	return center_of_mass;
}

/// @brief Clear the background and fill it with the background colour.
///
void clear_bg()
{
	using namespace graphics;
	glClearColor( bg_color.x(), bg_color.y(), bg_color.z(), 1.0 );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT ) ;
	return;
}

/// @brief Draw a frame for a window.
///
void draw_frame()
{
	using namespace graphics;

	glLineWidth(2.5);
	glColor3f(border_color.x(), border_color.y(), border_color.z());
	glBegin(GL_LINES);
	glVertex2f(-1.0, -1.0);
	glVertex2f(1.0, -1.0);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(-1.0, 1.0);
	glVertex2f(1.0, 1.0);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(1.0, 1.0);
	glVertex2f(1.0, -1.0);
	glEnd();
	glBegin(GL_LINES);
	glVertex2f(-1.0, 1.0);
	glVertex2f(-1.0, -1.0);
	glEnd();

	return;
}

/// @brief Fill window background with black.
///
void draw_black_bg()
{
	using namespace graphics;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glBegin(GL_QUADS);
	glColor3f(0.0f, 0.0f, 0.0f);
	glVertex2f(1.0,1.0);
	glVertex2f(-1.0,1.0);
	glVertex2f(-1.0,-1.0);
	glVertex2f(1.0,-1.0);
	glEnd();

	draw_frame();
}

/// @brief Draw a gradient for the window background.
///
void draw_gradient_bg()
{
	using namespace graphics;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glBegin(GL_QUADS);
	glColor3f(bg_color.x(), bg_color.y(), bg_color.z());
	glVertex2f(1.0,1.0);
	glVertex2f(-1.0,1.0);
	glColor3f(bg_color2.x(), bg_color2.y(), bg_color2.z());
	glVertex2f(-1.0,-1.0);
	glVertex2f(1.0,-1.0);
	glEnd();

	draw_frame();

	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void draw_pose(
	const core::pose::Pose & pose,
	GraphicsState & gs
) {
	core::Size nres = pose.total_residue();
	utility::vector1< char > ss(nres);
	utility::vector1< conformation::ResidueCOP > residues(nres);

	for ( int i=1; i<=(int)pose.total_residue(); ++i ) {
		ss[i] = pose.secstruct(i);
		residues[i] = pose.residue(i).get_self_ptr();
	}

	Vector center_of_mass = get_center( residues );

	draw_conformation( residues, ss, gs, center_of_mass );
}

#endif

/// @brief
utility::vector1< SilentObserverOP > silent_observers; // quick hack to eliminate OPs being destroyed!
void
add_monte_carlo_silent_viewer(
	moves::MonteCarlo & mc,
	std::string const & name_in,
	bool fullatom = false
)
{
	SilentObserverOP viewer1( new SilentObserver( name_in + "_last_accepted", fullatom ) );
	SilentObserverOP viewer2( new SilentObserver( name_in + "_best_accepted", fullatom ) );

	std::cout << "attaching viewers!" << std::endl;
	mc.attach_observer_to_last_accepted_pose( *viewer1 );
	mc.attach_observer_to_lowest_score_pose ( *viewer2 );

	silent_observers.push_back( viewer1 );
	silent_observers.push_back( viewer2 );
}

///////////////////////////////////////
void
clear_conformation_viewers()
{

#ifdef GL_GRAPHICS
	conformation_viewers.clear();
#endif
}


} // viewer
} // protocols
