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

#ifdef BOINC_GRAPHICS
#include <utility/io/izstream.hh>
#ifdef _WIN32
#include "boinc_win.h"
#else
#include <math.h>
#endif

#if _MSC_VER >= 1400
#endif


#include <protocols/boinc/boinc.hh>
#include <core/scoring/rms_util.hh>


#include <core/types.hh>
#include <core/init/init.hh>

// avoid having to create the static protocol movers
//#include <protocols/init/init.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/boinc.OptionKeys.gen.hh>

#include <utility/exit.hh>

#include "parse.h"
#include "util.h"
#include "gutil.h"
#include "boinc_gl.h"
#include "app_ipc.h"
#include "boinc_api.h"
#include "graphics2.h"
#include "graphics_data.h"

#include <core/io/serialization/serialize_pose.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.fwd.hh>
#include <protocols/boinc/boinc_shmem.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/viewer/GraphicsState.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/boinc.OptionKeys.gen.hh>

#undef read

// Project headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <vector>

// New boinc specific options

#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/Fstring.hh>
#include <ObjexxFCL/string.functions.hh>

#ifdef MAC
#include <GLUT/glut.h>
#elif _WIN32
#include <glut/glut.h>
#else
#include "GL/glut.h"
#endif

#include "txf_util.h"

#ifdef __APPLE__
#include "boinc_mac_app_icon.hh"
#endif

namespace graphics {
	float window_size = { 28 };
	float native_window_size = { 28 };
	int default_aspect_width = 6;
	int default_aspect_height = 4;
	float aspect_width, aspect_height;
	float aspect;
	int small_box;
	int mouse_x, mouse_y;
	bool mouse_down;
	int window_height;
	std::string where_in_window;
	GLfloat nativerotation[16], lowrotation[16], bestrotation[16], currentrotation[16];
	bool const allow_rotation = true;

	float default_text_color[4] = {1., 1., 1., 1.};
	float default_structure_text_color[4] = {.5f, .5f, .5f, 1.0f};

	float wu_text_box_height = 0.0;
	std::vector<std::string> wu_desc_rows;
	// tinker with these to modify how the text is displayed
	float wu_desc_rows_per_small_box = 6.0;

	// boinc data
	APP_INIT_DATA app_init_data;
	// shared memory
	protocols::boinc::BoincSharedMemory* shmem = NULL;

	double cpu_time=0;

	std::vector<float> low_rmsd_vector;
	std::vector<float> low_energy_vector;
	std::vector<float> last_accepted_rmsd_vector;
	std::vector<float> last_accepted_energy_vector;
	std::vector<float> model_rmsd_vector;
	std::vector<float> model_energy_vector;

	float last_accepted_rmsd;
	float low_energy_rmsd;
	unsigned int last_low_energy_update_cnt;
	unsigned int last_model_cnt;

	enum GraphicsPoseType {
	  CURRENT,
	  ACCEPTED,
	  LOW,
	  NATIVE
	};

	protocols::viewer::GraphicsState current_gs;

	static core::Size max_pose_nres = 0;
	static core::Size current_pose_nres = 0;
	static core::Size native_pose_nres = 0;

	bool native_exists = false;
	static core::pose::PoseOP nativeposeOP;
	static core::pose::PoseOP currentposeOP;
	static core::pose::PoseOP lastacceptedposeOP;
	static core::pose::PoseOP	lowenergyposeOP;

	// should probably make these maps, indexes by "low", etc.
	int low_viewport_x, low_viewport_y, low_viewport_width, low_viewport_height;
	int native_viewport_x,  native_viewport_y,  native_viewport_width,  native_viewport_height;
	int best_viewport_x,  best_viewport_y,  best_viewport_width,  best_viewport_height;
	int current_viewport_x, current_viewport_y, current_viewport_width, current_viewport_height;

}

void scale(const int & iw, const int & ih) {
	using namespace graphics;
	float aspect_ratio = (float)aspect_width/(float)aspect_height;
	float w=(float)iw, h=(float)ih;
	float xs, ys;
	if (h*aspect_ratio > w) {
		xs = 1.0f;
		ys = (w/aspect_ratio)/h;
	} else {
		xs = (h*aspect_ratio)/w;
		ys = 1.0f;
	}
	glScalef(xs, ys*aspect_ratio, 1.0f);
}

// call this to render 2D stuff, with 0..1x0..1 getting mapped
// to the full window.  You must call ortho_done() when done.
//
void mode_ortho_start() {
	glDisable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0, 1, 0, 1, 0, 1);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	gluLookAt(
		0.0, 0.0, 1.0,        // eye position
		0, 0, 0,              // where we're looking
		0.0, 1.0, 0.          // up is in positive Y direction
	);
	int viewport[4];
	glMatrixMode(GL_MODELVIEW);
	glGetIntegerv(GL_VIEWPORT, (GLint*)viewport);
	scale(viewport[2], viewport[3]);
}

void mode_ortho_done() {
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
	glLoadIdentity();
}

void get_bounds( std::vector< float > const & t, float & mn, float & mx ) {
	if ( !t.size() ) return;

	//OK, now look back -- if there was a huge jump at some point,
	// ignore points in the trajectory before then.
	int count(0);
	int LOOKBACK = 500;
	mn=t[t.size()-1];
	mx=mn;
	for ( int i = t.size()-1; i >= std::max(int(0),int( t.size()-LOOKBACK)); --i ) {
		mn = std::min( mn, t[i]);
		mx = std::max( mx, t[i]);
		count++;
	}
	mn-=0.00;
	mx+=0.00;
	return;
/*
	//How about using some estimate of the variance. This
	// will help us figure out if a big jump occurred recently.
	int const NUM_LOOK_BACK = 10000;
	int const i_start = std::max( static_cast<int>(t.size()) - NUM_LOOK_BACK, 0);
	int const i_end = t.size() - 1;
	int const numpoints = i_end - i_start + 1;

	float x (0.0);
	for ( int i=i_start+1; i<= i_end; ++i ) {
		float const score_jump = std::abs(t[i+1] - t[i]);
		if (score_jump < 10.0 ) x += score_jump;
	}
	float const avg_jump  = x / numpoints; //This should be of order 1 (for energy).

	//OK, now look back -- if there was a huge jump at some point,
	// ignore points in the trajectory before then.
	mn = t[ i_end ];
	mx = mn;
	int count(0);
	for ( int i = i_end-1; i >= i_start; --i ) {
		if (std::abs(t[i+1] - t[i]) > 50.0*avg_jump) break;
		mn = std::min( mn, t[i]);
		mx = std::max( mx, t[i]);
		count++;
	}

	*/
}


void
plot_2D(
	std::vector<float> const & xdata,
	std::vector<float> const & ydata,
	float const xmin,
	float const xmax,
	float const ymin,
	float const ymax
) {
	using namespace graphics;

	unsigned int total_steps = xdata.size();
	if (ydata.size() < total_steps) total_steps = ydata.size();

	if ( total_steps > 0 ) {

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		float maxdatax =log(1.0+xmax-xmin);
		float mindatax =log(1.0);
		float maxdatay = log(1.0+ymax-ymin);
		float mindatay = log(1.0);

		gluOrtho2D( mindatax, maxdatax, mindatay, maxdatay );
		glMatrixMode(GL_MODELVIEW);
		glDisable( GL_DEPTH_TEST );
		glLoadIdentity();

		glPointSize( 2.0 );
		glBegin( GL_POINTS );
		glColor3f(.5,.5,.5);
		for (unsigned int step = 0; step < total_steps; step++) {
			float x = xdata[step];
			float y = ydata[step];
			x -= xmin;
			if( x > 0 ) x = log(1.0+x);
			else continue;
			y -= ymin;
			if( y > 0 ) y = log(1.0+y);
			else continue;
			glVertex2f( x, y ) ;
		}
		glEnd();

		int decoys = (native_exists) ? model_rmsd_vector.size() :  model_energy_vector.size();
		if (decoys > 0) {
			glPointSize( 4.0 );
			glBegin( GL_POINTS );
			glColor3f(0.8f,0.13f,0.0f); //rust
			for (int step = 0; step < decoys; step++) {
				float x = (native_exists) ? model_rmsd_vector[step] : model_energy_vector[step];
				float y = model_energy_vector[step];
				if (x && y){
					x -= xmin;
					if( x > 0 ) x = log(1.0+x);
					else continue;
					y -= ymin;
					if( y > 0 ) y = log(1.0+y);
					else continue;

					glVertex2f( x, y );
				}
			}
			glEnd();
		}
/*
		// cross hairs
		glColor3f(.5,.5,.5);
		glBegin( GL_LINES );
		glVertex2f( xdata[total_steps-1], ymax ) ;
		glVertex2f( xdata[total_steps-1], ymin );
		glVertex2f( xmin, ydata[total_steps-1] );
		glVertex2f( xmax, ydata[total_steps-1] );
		glEnd();
*/
		glEnable( GL_DEPTH_TEST );
		glLoadIdentity();
	}
}

void start_rotate( GLfloat decoyrotation[16]) {
	glPushMatrix();
	glMultMatrixf(decoyrotation);
}

void end_rotate() {
	glPopMatrix();
}

void do_the_rotation(	const int & x, const int & y, const int & left, const int & middle,
											const int & right, GLfloat decoyrotation[16], const float & this_window_size ) {
	using namespace graphics;

	double delta_y = y-mouse_y;
	double delta_x = x-mouse_x;

	if (left) {
		double axis_x = delta_y;
		double axis_y = delta_x;
		double userangle = sqrt(delta_x*delta_x + delta_y*delta_y);

		glMatrixMode(GL_MODELVIEW);
		// Standard GLUT rotation is a postmultiplication - rotation around
		// molecule's inertial frame --  and leads to
		// non-intuitive behavior.
		//glRotatef(userangle,axis_x,axis_y,0.0);

		//A premultiplication -- rotates around the axis the user actually sees.
		// A little more complicated to code; unfortunately GLUT doesn't have a one-line
		// function for it.
		glLoadIdentity();
		glRotatef(userangle,axis_x,axis_y,0.0);
		glMultMatrixf(decoyrotation);
		glGetFloatv(GL_MODELVIEW_MATRIX, decoyrotation);  // Store current model view in decoyrotation
		mouse_y = y;
		mouse_x = x;
	}  else if (right){ //Zoom
		double s = exp( 1.0* (double) delta_y*0.01);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glMultMatrixf(decoyrotation);
		glScalef(s,s,s);
		glGetFloatv(GL_MODELVIEW_MATRIX, decoyrotation); // Store current model view in decoyrotation
		mouse_y = y;
		mouse_x = x;
	}
	else if (middle) { //Translate
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT,viewport);
		glMatrixMode(GL_MODELVIEW);
		double xscale = (this_window_size * 2.0)/(viewport[2]);
		double yscale = (this_window_size * 2.0)/(viewport[3]);
		glLoadIdentity();
		glTranslatef( delta_x * xscale, -delta_y * yscale, 0);
		glMultMatrixf(decoyrotation);
		glGetFloatv(GL_MODELVIEW_MATRIX, decoyrotation);
		mouse_y = y;
		mouse_x = x;
	}
	else {
		mouse_down = false;
	}
}

// BOINC GRAPHICS API
void app_graphics_resize(int w, int h){
	using namespace graphics;
	glViewport(0, 0, w, h);
}

// BOINC GRAPHICS API
void boinc_app_mouse_move(int x, int y, int left, int middle, int right) {
	using namespace graphics;

	if (!allow_rotation) return;

	//if (!native_centered) return; // just rotate native for now
	// Where are we? In low, native, accepted, or current?
	if ( where_in_window == "native" )
	  do_the_rotation( x, y, left, middle, right, nativerotation, native_window_size);
	if ( where_in_window == "low" )
	  do_the_rotation( x, y, left, middle, right, lowrotation, window_size);
	if ( where_in_window == "current" )
	  do_the_rotation( x, y, left, middle, right, currentrotation, window_size);
	if ( where_in_window == "best" )
	  do_the_rotation( x, y, left, middle, right, bestrotation, window_size);
}

// BOINC GRAPHICS API
void boinc_app_mouse_button(int x, int y, int, int is_down) {
	using namespace graphics;

	if (!allow_rotation) return;

	if (is_down) {
	      mouse_down = true;
	      mouse_x = x;
	      mouse_y = y;
	  } else {
	      mouse_down = false;
	  }

	if ( mouse_x > native_viewport_x && mouse_x < native_viewport_x + native_viewport_width &&
	     mouse_y > window_height -  native_viewport_y - native_viewport_height&& mouse_y < window_height - native_viewport_y ) {
	  where_in_window = "native"; //native
	} else if
	    ( mouse_x > low_viewport_x && mouse_x < low_viewport_x + low_viewport_width &&
	      mouse_y > window_height - low_viewport_y -low_viewport_height && mouse_y < window_height - low_viewport_y ) {
	  where_in_window = "low";
	} else if
	    ( mouse_x > current_viewport_x && mouse_x < current_viewport_x + current_viewport_width &&
	      mouse_y > window_height - current_viewport_y - current_viewport_height && mouse_y < window_height - current_viewport_y ) {
	  where_in_window = "current";
	} else if
	    ( mouse_x > best_viewport_x && mouse_x < best_viewport_x + best_viewport_width &&
	      mouse_y > window_height - best_viewport_y - best_viewport_height && mouse_y < window_height - best_viewport_y ) {
	  where_in_window = "best";
	} else {
	  where_in_window = ""; //somewhere else
	}
}

// BOINC GRAPHICS API
void boinc_app_key_press(int key, int //lParam           // system-specific key encodings
){
	using namespace graphics;
	using namespace protocols::viewer;
	if (key == 67 || key == 99) { //'c' control color
		current_gs.Color_mode = ColorMode ( current_gs.Color_mode + 1 );
		if ( current_gs.Color_mode > RESIDUE_CPK_COLOR ) current_gs.Color_mode = RAINBOW_COLOR;
	}
	if (key == 66 || key == 98) { //'b' control backbone display
		current_gs.BBdisplay_state = BBdisplayState ( current_gs.BBdisplay_state + 1 );
		if ( current_gs.BBdisplay_state > SHOW_BACKBONE ) current_gs.BBdisplay_state = SHOW_NOBB;
	}
	if (key == 83 || key == 115) { //'s' control sidechain display
		current_gs.SCdisplay_state = SCdisplayState ( current_gs.SCdisplay_state + 1 );
		if ( current_gs.SCdisplay_state > SHOW_WIREFRAME ) current_gs.SCdisplay_state = SHOW_NOSC;
	}
}

// BOINC GRAPHICS API
void boinc_app_key_release(int, int){}

// BOINC GRAPHICS API
void app_graphics_init() {
	using namespace graphics;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	/* Open the Semaphore */
	// for shared memory data sychronization
	protocols::boinc::Boinc::get_semaphore();

	// boinc_graphics_get_shmem() must be called after
	// boinc_parse_init_data_file()
	// Put this in the main loop to allow retries if the
	// worker application has not yet created shared memory
	//
	if (!option[ OptionKeys::boinc::noshmem ]()) {
		protocols::boinc::Boinc::attach_shared_memory();
		shmem =  protocols::boinc::Boinc::get_shmem();
	}

	// get workunit description

	// display work init name first
	std::string wu_name = app_init_data.wu_name;
	if (wu_name.length() > 0) wu_desc_rows.push_back( wu_name );

	if (shmem && shmem->wu_desc_exists) {
		if(!protocols::boinc::Boinc::wait_semaphore()){
			std::vector<std::string> tmp_wu_desc_rows;
			core::io::serialization::BUFFER b((char*)(&shmem->wu_desc_buf ),protocols::boinc::WU_DESC_TEXT_BUFSIZE);
			core::io::serialization::read_binary(tmp_wu_desc_rows, b);
			wu_desc_rows.insert(wu_desc_rows.end(), tmp_wu_desc_rows.begin(), tmp_wu_desc_rows.end());
			protocols::boinc::Boinc::unlock_semaphore();
		}
	}

	aspect_width = float(default_aspect_width);
	aspect_height = float(default_aspect_height);
	if (wu_desc_rows.size() > 0) {
		//small box is 1 aspect_height unit
		wu_text_box_height = float(wu_desc_rows.size())/wu_desc_rows_per_small_box;
		aspect_height += wu_text_box_height;
	}

#ifdef _WIN32
	boinc_set_windows_icon("boinc_rosetta32","boinc_rosetta32");
#endif

	// read project specific prefs
	protocols::boinc::Boinc::read_and_set_project_prefs();
	boinc_max_fps = protocols::boinc::Boinc::get_project_pref_max_gfx_fps();
	boinc_max_gfx_cpu_frac = protocols::boinc::Boinc::get_project_pref_max_gfx_cpu()/100.0;

	// for testing
	if (option[ OptionKeys::boinc::max_fps ].user())
		boinc_max_fps = option[ OptionKeys::boinc::max_fps ]();
	if (option[ OptionKeys::boinc::max_cpu ].user())
	    boinc_max_gfx_cpu_frac = (double)option[ OptionKeys::boinc::max_cpu ]()/100.0;

	glClearColor (0.0, 0.0, 0.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	int size = 30;
	glOrtho(-size, size, -size, size, -size, size);
	glMatrixMode(GL_MODELVIEW);
	glGetFloatv(GL_MODELVIEW_MATRIX, nativerotation);
	glGetFloatv(GL_MODELVIEW_MATRIX, lowrotation);
	glGetFloatv(GL_MODELVIEW_MATRIX, currentrotation);
	glGetFloatv(GL_MODELVIEW_MATRIX, bestrotation);

	glEnable( GL_DEPTH_TEST );
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	// Expects a .txf file to exist in run directory
	// TEXT WILL NOT DISPLAY IF .txf FILE IS MISSING!!
	txf_load_fonts(".");
}

void
writeStrokeString( const std::string & text_string, GLfloat *col, const float & xpos, const float & ypos, int scalefactor=350 ) {
	// USE NICE FONTS
	// Requires .txf file in run directory
	char buf[256];
	sprintf(buf, "%s", text_string.c_str());
	// FROM BOINC API fscale bigger is smaller
	txf_render_string(.1, xpos, ypos, 0, scalefactor, col, 0, buf);
}

void display_text() {
	using namespace graphics;
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::format;
	///////////////////////////////////////////////////////////////
	// START LEFT COLUMN

	// STAGE TEXT
	std::string mode_string = "Stage: unknown";
	if (shmem) {
		if (shmem->total_mc_trial_count == 0 && last_model_cnt == 0 ) {
			// initializing?
			// show it's doing something
			static std::string initializing;
			mode_string = " Initializing";
			initializing = (initializing.length() >= 10) ? " ." : initializing + " .";
			mode_string += initializing;
		} else if (shmem->total_mc_trial_count > 0) {
			// Stage: mover type info
			std::string mover_type;
			if(!protocols::boinc::Boinc::wait_semaphore()){
				core::io::serialization::BUFFER b((char*)(&shmem->mover_type_text ),protocols::boinc::TEXT_BUFSIZE);
				core::io::serialization::read_binary(mover_type,b);
				protocols::boinc::Boinc::unlock_semaphore();
			}
			mode_string = "Stage: " + mover_type;
		}
		// update cpu run time
		cpu_time = shmem->cpu_time;
		float dt = dtime() - shmem->update_time;
		if (dt > 10) {
			boinc_close_window_and_quit("shmem not updated");
		} else if (dt > 5) {
			mode_string = "App not running - exiting in 5 seconds";
		} else if (shmem->status.suspended) {
			mode_string = "App suspended";
    }
	}
	mode_ortho_start();
	glTranslatef( .01, 0.55, 0.0f );
	writeStrokeString( mode_string, default_text_color, 0.0f, 0.0f );
	mode_ortho_done();

	// CPU TIME TEXT
	int hours = int (cpu_time / 3600.0);
	int minutes = int ((cpu_time-3600*hours) / 60.0);
	int seconds = int ((cpu_time-3600*hours-60*minutes));
	std::string cpu_time_str = "CPU time: " + ObjexxFCL::string_of(hours) + " hr "
		+ ObjexxFCL::string_of(minutes) + " min "
		+ ObjexxFCL::string_of(seconds) + " sec ";

	mode_ortho_start();
	glTranslatef( .01, 0.45, 0.0f );
	writeStrokeString( cpu_time_str, default_text_color, 0.0f, 0.0f);
	mode_ortho_done();

	// USER TEXT
	std::string user_string = app_init_data.user_name;
	std::string user_total_credit_string = ObjexxFCL::string_of(app_init_data.user_total_credit,  6);
	std::string user_expavg_credit_string = ObjexxFCL::string_of(app_init_data.user_expavg_credit , 6);
	user_string = user_string + " - Total credit: " + user_total_credit_string
		+ " - RAC: " + user_expavg_credit_string;

	mode_ortho_start();
	glTranslatef( .01, 0.35, 0 );
	writeStrokeString( user_string, default_text_color, 0.0f, 0.0f );
	mode_ortho_done();

	// TEAM TEXT
	std::string team_string = app_init_data.team_name;
	if (team_string != "") {
		mode_ortho_start();
		glTranslatef( .01, .25, 0.0f );
		writeStrokeString( team_string, default_text_color, 0.0f, 0.0f);
		mode_ortho_done();
	}

	// LOGO TEXT
	std::string appver = ObjexxFCL::string_of(app_init_data.app_version/100.0, 3);
	std::string logo_string = "Rosetta@home v" + appver + " http://boinc.bakerlab.org/rosetta/";

	mode_ortho_start();
	glTranslatef( .01, .01, 0.0 );
	writeStrokeString( logo_string, default_text_color, 0.0f, 0.0f, 275 );
	mode_ortho_done();

	///////////////////////////////////////////////////////////////
	// START RIGHT COLUMN
	float xshift = 2.7;

	// RUN STATUS -  PERCENT COMPLETE TEXT
	std::string run_status_str;
	if (shmem) {
		float boinc_pct_complete = 100.0 * shmem->fraction_done;
		run_status_str = F( 7, 2, boinc_pct_complete ) + "% Complete";
	} else {
		run_status_str = "No shared mem";
	}

	mode_ortho_start();
	glTranslatef( xshift, .55, 0.0f );
	writeStrokeString( run_status_str, default_text_color, 0.0f, 0.0f );
	mode_ortho_done();

	// MODEL AND STEP TEXT
	if (shmem) {
		std::string ntrials_string = "Model: " + string_of(shmem->model_count, 4) +
			"  Step: " + string_of(shmem->total_mc_trial_count);
		mode_ortho_start();
		glTranslatef( xshift, .45, 0.0f );
		writeStrokeString( ntrials_string, default_text_color, 0.0f, 0.0f );
		mode_ortho_done();

	  // LAST ACCEPTED ENERGY TEXT
	  std::string score_string2 = "Accepted Energy: " + string_of(shmem->last_accepted_energy, 7);
	  mode_ortho_start();
	  glTranslatef( xshift, .35, 0.0f );
	  writeStrokeString( score_string2, default_text_color, 0.0f, 0.0f );
	  mode_ortho_done();

		double lowenergytextypos = 0.25;

		// LAST ACCEPTED RMSD TEXT
		if ( native_exists ) {
			lowenergytextypos = 0.15;
			std::string acc_rmsd_string = "Accepted RMSD:  " + string_of(last_accepted_rmsd, 4);
			mode_ortho_start();
			glTranslatef( xshift, .25, 0.0f );
			writeStrokeString( acc_rmsd_string, default_text_color, 0.0f, 0.0f );
			mode_ortho_done();
		}

		// LOW ENERGY TEXT
		std::string score_string3 = "Low Energy: " + string_of(shmem->low_energy, 7);
		mode_ortho_start();
		glTranslatef( xshift, lowenergytextypos, 0.0f );
		writeStrokeString( score_string3, default_text_color, 0.0f, 0.0f );
		mode_ortho_done();

		// LOW ENERGY RMSD TEXT
		if ( native_exists ) {
			std::string low_rmsd_string = "Low RMSD:  " + string_of(low_energy_rmsd, 4);
			mode_ortho_start();
			glTranslatef( xshift, .05, 0.0f );
			writeStrokeString( low_rmsd_string, default_text_color, 0.0f, 0.0f );
			mode_ortho_done();
		}
	}
}

void
Structure_display ( const graphics::GraphicsPoseType & graphics_pose_type, const float & this_window_size )
{
	using namespace graphics;
	if (!shmem || !current_pose_nres) return;

	float y_min_screen = - this_window_size;
	float y_max_screen = + this_window_size;
	float x_min_screen = - this_window_size;
	float x_max_screen = + this_window_size;
	float zmax = 300.0;

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(x_min_screen, x_max_screen, y_min_screen, y_max_screen, -zmax, zmax);
  glMatrixMode(GL_MODELVIEW);

	switch( graphics_pose_type ){
		case CURRENT:
			start_rotate(currentrotation);
			break;
		case ACCEPTED:
			start_rotate(bestrotation);
			break;
		case LOW:
			start_rotate(lowrotation);
			break;
		case NATIVE:
			start_rotate(nativerotation);
			break;
	}

	switch( graphics_pose_type ){
		case CURRENT:
			protocols::viewer::draw_pose( *currentposeOP, current_gs );
			break;
		case ACCEPTED:
			protocols::viewer::draw_pose( *lastacceptedposeOP, current_gs );
			break;
		case LOW:
			protocols::viewer::draw_pose( *lowenergyposeOP, current_gs );
			break;
		case NATIVE:
			protocols::viewer::draw_pose( *nativeposeOP, current_gs );
			break;
	}
	glLoadIdentity();
	end_rotate();
}

void
display_wu_desc( const float & height ) {
	using namespace graphics;
	float rowheight = float(small_box)/wu_desc_rows_per_small_box;
	for(unsigned int i=1;i <= wu_desc_rows.size(); i++) {
		glViewport( 0, height - 3.0*small_box-rowheight*(float)i, small_box*aspect_width, rowheight);
		mode_ortho_start();
		writeStrokeString( wu_desc_rows.at(i-1), default_text_color, 0.01, 0.0, 55);
		mode_ortho_done();
	}
}


void plot_timeseries(
	std::vector<float> const & data,
	std::vector<float> const & data2,
	bool const & vertical,
	float const & min,
	float const & max
) {
	using namespace graphics;

	//int logsteps = int( ((float)log(float(data.size() )*5.0f) ) );
	//unsigned int total_steps = pow( log(1),

	unsigned int total_steps;
	total_steps = 10*(1+int(data.size()/10));
	total_steps = std::max( (unsigned int) (50), (unsigned int)( pow(4, 1+floor(log( float(data.size()) )/log(4.0f)) ) ) );
	int start_step = 0;
	unsigned int LOOK_BACK =500;
	if( total_steps > LOOK_BACK ){
		total_steps = LOOK_BACK;
		start_step = data.size() - total_steps;
	}
	start_step = data.size() - total_steps;
	if( start_step < 0 ) start_step = 0;

	if ( data.size() > 0) {
		float mindata,maxdata;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		if (vertical) {
			//gluOrtho2D( min, max, 0, total_steps-1  );
			maxdata =log(1.0+max-min);
			mindata =log(1.0);
			gluOrtho2D( mindata, maxdata, 0, total_steps-1  );

		} else {
			//gluOrtho2D( 0, total_steps-1, min, max );
			maxdata = log(1.0+max-min);
			mindata = log(1.0);
			//maxdata = max; //log(1+max-min);
			//mindata = min; //log(1.0);
			gluOrtho2D( 0, total_steps-1, mindata, maxdata );
		}
		glMatrixMode(GL_MODELVIEW);
	  //  glColor3f(0.1f,0.7f,0.7f); //light blue
	  //  glColor3f(1.0f,0.0f,0.0f); // red
		glDisable( GL_DEPTH_TEST );
		glLoadIdentity();
		glColor3f(.5,.5,.5);

		glColor3f(0.0f,1.0f,0.0f); // yellow
		float min_point = 1e12;
		glLineWidth( 2.0 );
		glBegin( GL_LINE_STRIP );
		for (unsigned int step = start_step; step < data.size(); step++) {
			float x, y;
			float caged_data,point = data[step];
			if ( point < min_point ) min_point = point;
			if (vertical) {
				x = point;
				x -= min;
				if( x > 0 ) x = log(1.0+x);
				y = (total_steps-1) - step+start_step;
				caged_data = std::min(maxdata, std::max( (float)mindata, (float)x ) ) / (maxdata+0.001);
				glColor3f(1.0-caged_data,1.0-caged_data,caged_data); // yellow
		} else {
				x = step-start_step;
				y = point;
				y -= min;
				if( y > 0 ) y = log(1.0+y);

				caged_data = std::min(maxdata, std::max( (float)mindata, (float)y ) ) / (maxdata+0.001);
				glColor3f(1.0-caged_data,1.0-caged_data,caged_data); // yellow

			}
			glVertex2f( x, y ) ;
		}
		glEnd();
		glLineWidth( 1.0 );
		glLoadIdentity();
	}
}

void clear_trajectory() {
	using namespace graphics;
	static bool init = false;
	if ( ! init ) {
		low_rmsd_vector.reserve(300000);
		low_energy_vector.reserve(300000);
		last_accepted_rmsd_vector.reserve(300000);
		last_accepted_energy_vector.reserve(300000);
		init = true;
	}
	low_rmsd_vector.clear();
	low_energy_vector.clear();
	last_accepted_energy_vector.clear();
	last_accepted_rmsd_vector.clear();
}

void get_shmem_structures() {
	using namespace graphics;
	if (!shmem) return;

	// need to make sure we are in sync with the worker app
	if(!protocols::boinc::Boinc::wait_semaphore()){

		// get native pose
		if (shmem->native_pose_exists) {
			// get native
			core::io::serialization::BUFFER bn((char*)(&shmem->native_pose_buf ),protocols::boinc::POSE_BUFSIZE);
			core::io::serialization::read_binary(*nativeposeOP,bn);
			native_pose_nres = (*nativeposeOP).total_residue();
			native_exists = true;
		}

		// get current pose
		if (shmem->current_pose_exists) {
			core::io::serialization::BUFFER bc((char*)(&shmem->current_pose_buf ),protocols::boinc::POSE_BUFSIZE);
			core::io::serialization::read_binary(*currentposeOP,bc);
			current_pose_nres =  (*currentposeOP).total_residue();
		}

		// get last accepted pose and calculate rmsd
		if (shmem->last_accepted_pose_exists) {
			core::io::serialization::BUFFER bla((char*)(&shmem->last_accepted_pose_buf ),protocols::boinc::POSE_BUFSIZE);
			core::io::serialization::read_binary(*lastacceptedposeOP,bla);
			if (native_exists && (*lastacceptedposeOP).total_residue() > 0)
				last_accepted_rmsd = core::scoring::native_CA_rmsd( *nativeposeOP, *lastacceptedposeOP);
		}

		// get low energy pose and calculate rmsd
		if (shmem->low_energy_pose_exists) {
			core::io::serialization::BUFFER ble((char*)(&shmem->low_energy_pose_buf ),protocols::boinc::POSE_BUFSIZE);
			core::io::serialization::read_binary(*lowenergyposeOP,ble);
			last_low_energy_update_cnt = shmem->low_energy_update_cnt;
			if (native_exists && (*lowenergyposeOP).total_residue() > 0)
				low_energy_rmsd = core::scoring::native_CA_rmsd( *nativeposeOP, *lowenergyposeOP);
		}

		// next model? clear trajectory
		if (shmem->model_count != last_model_cnt) {
			// get low energy model rmsd and energy
			model_energy_vector.push_back( shmem->model_low_energy );
			if (native_exists)
				model_rmsd_vector.push_back( shmem->model_low_energy_rmsd );
			clear_trajectory();
			last_model_cnt = shmem->model_count;
		}

		// get monte carlo trial info (resolved to graphics frame rate)
		if (shmem->total_mc_trial_count > 0) {
			low_energy_vector.push_back( shmem->low_energy );
			low_rmsd_vector.push_back( low_energy_rmsd );
			last_accepted_energy_vector.push_back( shmem->last_accepted_energy );
			last_accepted_rmsd_vector.push_back( last_accepted_rmsd );
		}

		max_pose_nres = (native_pose_nres > current_pose_nres) ? native_pose_nres : current_pose_nres;

		protocols::boinc::Boinc::unlock_semaphore();
	}

}


void draw_rosetta_screensaver( int & width, int & height )
{
	using namespace graphics;
	static int last_time_graphic_switch = time(NULL);

	glViewport(0 , 0, width, height );
	window_height = height; // needed to track mouse position.
	aspect = float(width)/float(height);

	get_shmem_structures(); // !!

	// users complaining that they can't see the whole protein...
	if (max_pose_nres > 0) {
		window_size = (current_pose_nres > 100) ? 0.7*(28 + ( current_pose_nres - 100)*0.15) : 28;
		native_window_size = (native_pose_nres > 100) ? 0.7*(28 + ( native_pose_nres - 100)*0.15) : 28;
	}

	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glViewport(0 , 0, width, height );
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0 , 1, 0, 1);
	glMatrixMode(GL_MODELVIEW);
	glEnable( GL_DEPTH_TEST );

	if (aspect >= aspect_width / aspect_height) {
	  // add space to right side (wider than tall)
	  small_box = int(height/aspect_height);
	} else {
	  // add space to bottom (taller than wide)
	  small_box = int(width/aspect_width);
	}
	int dim_main = 2*small_box;


	// generate box above text box (need this to prevent lines from
	// disappearing in graph and plot boxes for some odd reason)
	glViewport( 0, height-3*small_box, small_box*6, small_box*3 );
	gluOrtho2D(0,1,0,aspect);
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	glColor3f( 0.5f, 0.5f, 0.5f );
	glTranslatef( 0.0, 0.0, 0.0 );
	glBegin( GL_LINE_STRIP ) ;
	glVertex2f(0,0); glVertex2f(0,1); glVertex2f(1,1); glVertex2f(1,0); glVertex2f(0,0);
	glEnd();
	glLoadIdentity();

	// generate box above text box (need this to prevent lines from
	// disappearing in graph and plot boxes for some odd reason)
	glViewport( 0, height-3*small_box, small_box*6, small_box );
	gluOrtho2D(0,1,0,aspect);
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	glColor3f( 0.5f, 0.5f, 0.5f );
	glTranslatef( 0.0, 0.0, 0.0 );
	glBegin( GL_LINE_STRIP ) ;
	glVertex2f(0,0); glVertex2f(0,1); glVertex2f(1,1); glVertex2f(1,0); glVertex2f(0,0);
	glEnd();
	glLoadIdentity();

	// generate boxes for plots
	glViewport( small_box*5, height-3*small_box, small_box, small_box );
	gluOrtho2D(0,1,0,aspect);
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	glColor3f( 0.5f, 0.5f, 0.5f );
	glTranslatef( 0.0, 0.0, 0.0 );
	glBegin( GL_LINE_STRIP ) ;
	glVertex2f(0,0); glVertex2f(0,1); glVertex2f(1,1); glVertex2f(1,0); glVertex2f(0,0);
	glEnd();
	glLoadIdentity();

	glViewport( 0, height-3*small_box, 5*small_box, small_box );
	gluOrtho2D(0,1,0,aspect);
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();
	glColor3f( 0.5f, 0.5f, 0.5f );
	glTranslatef( 0.0, 0.0, 0.0 );
	glBegin( GL_LINE_STRIP ) ;
	glVertex2f(0,0); glVertex2f(0,1); glVertex2f(1,1); glVertex2f(1,0); glVertex2f(0,0);
	glEnd();
	glLoadIdentity();

	// Work Unit Description?
	if (wu_desc_rows.size() > 0) {

	  // generate box for description
		glViewport( 0, height - 3*small_box - wu_text_box_height*small_box, small_box*aspect_width, wu_text_box_height*small_box );
	  gluOrtho2D(0,1,0,aspect);
	  glMatrixMode( GL_MODELVIEW );
	  glLoadIdentity();
	  glColor3f( 0.5f, 0.5f, 0.5f );
	  glTranslatef( 0.0, 0.0, 0.0 );
	  glBegin( GL_LINE_STRIP ) ;
	  glVertex2f(0,1); glVertex2f(0,0); glVertex2f(1,0); glVertex2f(1,1); //glVertex2f(0,0);
	  glEnd();
	  glLoadIdentity();

		display_wu_desc( height );

	  // Main text box
	  glViewport( 0, height - 4*small_box-int(wu_text_box_height*small_box),
				int(small_box*aspect_width), small_box );
	  display_text();

	} else {
		glViewport( 0, height-4*small_box, small_box*6, small_box );
		display_text();
	}

	// SEARCHING BOX
	current_viewport_x      = 0;
	current_viewport_y      = height-dim_main;
	current_viewport_width  = current_viewport_height = dim_main;

	glViewport( current_viewport_x, current_viewport_y, current_viewport_width, current_viewport_height );
	mode_ortho_start();
  glColor3f( 0.5f, 0.5f, 0.5f );
  glTranslatef( 0.0, 0.0, 0.0 );
  glBegin( GL_LINE_STRIP ) ;
  glVertex2f(0,1); glVertex2f(1,1); glVertex2f(1,0);
  glEnd();


	using namespace graphics;
	using namespace protocols::viewer;


	// change view_type randomly
	if( (time(NULL) - last_time_graphic_switch) > 100 ){
		last_time_graphic_switch = time(NULL);

		if( rand() % 2 == 0 ) current_gs.BBdisplay_state = SHOW_CARTOON;
		else                  current_gs.BBdisplay_state = SHOW_BACKBONE;

		if( rand() % 3 == 0 && max_pose_nres < 500 ) current_gs.SCdisplay_state = SHOW_STICK;
		else                  current_gs.SCdisplay_state = SHOW_NOSC;

		if( current_gs.SCdisplay_state != SHOW_NOSC ){
			int randnum = rand() % 4;
			if( randnum == 0 ) current_gs.Color_mode = RAINBOW_COLOR;
			if( randnum == 1 ) current_gs.Color_mode = RAINBOW_CPK_COLOR;
			if( randnum == 2 ) current_gs.Color_mode = RESIDUE_CPK_COLOR;
			if( randnum == 3 ) current_gs.Color_mode = CPK_COLOR;
		}else{
			current_gs.Color_mode = RAINBOW_COLOR;
		}
	}

	writeStrokeString( "Searching...", default_structure_text_color, 0.0f, 0.0f, 280);
	mode_ortho_done();
	// do slow changes
	float slowrotate = 0.8;
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(slowrotate,0,1,0.0);
	glMultMatrixf(graphics::nativerotation);
	glGetFloatv(GL_MODELVIEW_MATRIX, graphics::nativerotation);  // Store current model view in decoyrotation
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(slowrotate,0,1,0.0);
	glMultMatrixf(graphics::bestrotation);
	glGetFloatv(GL_MODELVIEW_MATRIX, graphics::bestrotation);  // Store current model view in decoyrotation
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(slowrotate,0,1,0.0);
	glMultMatrixf(graphics::currentrotation);
	glGetFloatv(GL_MODELVIEW_MATRIX, graphics::currentrotation);  // Store current model view in decoyrotation
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glRotatef(slowrotate,0,1,0.0);
	glMultMatrixf(graphics::lowrotation);
	glGetFloatv(GL_MODELVIEW_MATRIX, graphics::lowrotation);  // Store current model view in decoyrotation

	Structure_display(CURRENT, window_size);

	// ACCEPTED BOX
	best_viewport_x      = dim_main;
	best_viewport_y      = height-dim_main;
	best_viewport_width  = best_viewport_height = dim_main;
	glViewport( best_viewport_x, best_viewport_y, best_viewport_width, best_viewport_height );
	mode_ortho_start();
	glColor3f( 0.5f, 0.5f, 0.5f );
	glTranslatef( 0.0, 0.0, 0.0 );
	glBegin( GL_LINE_STRIP ) ;
	glVertex2f(0,1); glVertex2f(1,1); glVertex2f(1,0);
	glEnd();
	writeStrokeString( "Accepted", default_structure_text_color, 0.0f, 0.0f, 280);
	mode_ortho_done();
	Structure_display(ACCEPTED, window_size);

	float rms_min = 0.0;
	float rms_max = 0.0;
	get_bounds( last_accepted_rmsd_vector, rms_min, rms_max );
	if (rms_max > 20) rms_max = 20.00;
	rms_max += 1.00;
	rms_min -= 0.20;
	if (rms_min < 1) rms_min = 0.01;

	float energy_min = 0;
	if( last_accepted_energy_vector.size() >= 1 )
		energy_min = last_accepted_energy_vector[last_accepted_energy_vector.size()-1];
	float energy_max = energy_min;
	get_bounds( last_accepted_energy_vector, energy_min, energy_max );
	energy_max += 10.00;
	energy_min -= 1.00;

	if ( native_exists ) {

		// LOW ENERGY BOX (SMALL)
	  low_viewport_x      = 2 * dim_main;
	  low_viewport_y      = height - small_box;
	  low_viewport_width  = low_viewport_height = small_box;
	  glViewport( low_viewport_x, low_viewport_y, low_viewport_width, low_viewport_height );
		mode_ortho_start();
		glColor3f( 0.5f, 0.5f, 0.5f );
		glTranslatef( 0.0, 0.0, 0.0 );
		glBegin( GL_LINE_STRIP ) ;
		glVertex2f(0,1); glVertex2f(1,1); glVertex2f(1,0);
		glEnd();
		writeStrokeString( "Low Energy", default_structure_text_color, 0.0f, 0.0f, 200);
		mode_ortho_done();
	  Structure_display(LOW, window_size);

	  // NATIVE BOX (SMALL)
	  native_viewport_x      = 2 * dim_main;
	  native_viewport_y      = height - 2*small_box;
	  native_viewport_width  = native_viewport_height = small_box;
	  glViewport( native_viewport_x, native_viewport_y, native_viewport_width, native_viewport_height );
		mode_ortho_start();
    glColor3f( 0.5f, 0.5f, 0.5f );
    glTranslatef( 0.0, 0.0, 0.0 );
    glBegin( GL_LINE_STRIP ) ;
    glVertex2f(0,1); glVertex2f(1,1); glVertex2f(1,0);
    glEnd();
		writeStrokeString( "Native", default_structure_text_color, 0.0f, 0.0f, 200);
		mode_ortho_done();
	  Structure_display(NATIVE, native_window_size);

	  // RMSD TIME SERIES GRAPH BOX
	  glViewport( 5*small_box, height-2*small_box, small_box, 2*small_box );
		mode_ortho_start();
		writeStrokeString( "RMSD", default_structure_text_color, 0.0f, 0.0f, 180 );
		mode_ortho_done();
		plot_timeseries( last_accepted_rmsd_vector,low_rmsd_vector ,true, rms_min, rms_max );

	} else {

		// LOW ENERGY BOX (LARGE WITHOUT NATIVE)
		low_viewport_x      = 2 * dim_main;
		low_viewport_y      = height - dim_main;
		low_viewport_width  = low_viewport_height = dim_main;
	  glViewport( 2*dim_main, height-dim_main, dim_main, dim_main );
		mode_ortho_start();
		writeStrokeString( "Low Energy", default_structure_text_color, 0.0f, 0.0f, 280);
		mode_ortho_done();
	  Structure_display(LOW, window_size);

	}

	// ENERGY TIME SERIES GRAPH BOX
	glViewport( 0, height-3*small_box, 5*small_box, small_box );
	mode_ortho_start();
	writeStrokeString( "Accepted Energy", default_structure_text_color, 0.0f, 0.0f, 250 );
	mode_ortho_done();
	plot_timeseries( last_accepted_energy_vector, low_energy_vector, false, energy_min, energy_max );

	// ENERGY VS RMSD 2D PLOT
	glViewport( 5*small_box, height-3*small_box, small_box, small_box );
	if( native_exists)
		plot_2D( last_accepted_rmsd_vector, last_accepted_energy_vector,	rms_min, rms_max, energy_min, energy_max );


	glFlush();
}

// BOINC GRAPHICS API
// This will be called periodically in the graphics thread. It
// should generate the current graphic. xs and ys are the X and Y
// sizes of the window, and time_of_day is the relative time in
// seconds. The function should return true if it actually drew
// anything. It can refer to the user name, CPU time etc. obtained
// from boinc_get_init_data(). Applications that don't do graphics
// must also supply a dummy app_graphics_render() to link with the
// API.
void app_graphics_render(int xs, int ys, double) { draw_rosetta_screensaver( xs, ys ); }

#if _MSC_VER >= 1400
void
RosettaInvalidParameterHandler(
	const wchar_t* expression,
	const wchar_t* function,
	const wchar_t* file,
	unsigned int line,
	uintptr_t pReserved
)
{
	fprintf(
		stderr,
		"Invalid parameter detected in function %s. File: %s Line: %d\n",
		function,
		file,
		line
	);
	fprintf(
		stderr,
		"Expression: %s\n",
		expression
	);
	DebugBreak();
}
#endif


int main(int argc, char** argv) {
	try {

	using namespace graphics;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::cerr << "Starting graphics.." << std::endl;
	// new boinc client specific options
	option.add_relevant( OptionKeys::boinc::graphics );
	option.add_relevant( OptionKeys::boinc::fullscreen );
	option.add_relevant( OptionKeys::boinc::max_fps );
	option.add_relevant( OptionKeys::boinc::max_cpu );
	option.add_relevant( OptionKeys::boinc::noshmem );

#ifdef _DEBUG
	// DIAGNOSTICS
	// http://boinc.berkeley.edu/trac/wiki/DiagnosticsApi
	boinc_init_graphics_diagnostics(
		BOINC_DIAG_DUMPCALLSTACKENABLED |
		BOINC_DIAG_HEAPCHECKENABLED |
		BOINC_DIAG_MEMORYLEAKCHECKENABLED |
		BOINC_DIAG_REDIRECTSTDERR |
		BOINC_DIAG_REDIRECTSTDOUT |
		BOINC_DIAG_TRACETOSTDERR
	);
#else
 boinc_init_graphics_diagnostics(BOINC_DIAG_DEFAULTS);
#endif


// FROM rosetta++ main.cc
#if _MSC_VER >= 1400
  // Every once and awhile something looks at a std::vector or some other
  // CRT/STL construct that throws an exception.  In this case we should
  // dump whatever information we can and then bail.  When we bail we
  // should dump as much data as possible.
  _set_invalid_parameter_handler(RosettaInvalidParameterHandler);
#endif


#ifdef __APPLE__
	// mac icon
	setMacIcon(argv[0], MacAppIconData, sizeof(MacAppIconData));
#endif

	// create boinc object
	protocols::boinc::Boinc boinc_wu = protocols::boinc::Boinc::instance();

  core::init::init( argc, argv );

	// avoid having to create the static protocol movers
	//protocols::init::init( argc, argv );

	// override database option and set to current directory
	option[in::path::database].value("minirosetta_database");
#ifndef _DEBUG
	option[out::mute].value("all");
#endif
	// BOINC GRAPHICS!
	// get BOINC app init data - user, team, etc. info
	boinc_parse_init_data_file();
	boinc_get_init_data(app_init_data);

	nativeposeOP =  core::pose::PoseOP( new core::pose::Pose() );
	currentposeOP = core::pose::PoseOP( new core::pose::Pose() );
	lowenergyposeOP = core::pose::PoseOP( new core::pose::Pose() );
	lastacceptedposeOP = core::pose::PoseOP( new core::pose::Pose() );

	boinc_graphics_loop(argc, argv);

#ifdef _DEBUG
	  boinc_finish_diag();
#endif
	 } catch ( utility::excn::EXCN_Base const & e ) {
		 std::cout << "caught exception " << e.msg() << std::endl;
		 return -1;
	}
}

#else // BOINC_GRAPHICS!

#include <iostream>

int main(int /*argc*/, char** /*argv*/) {
	std::cout << "Build with boinc api (scons extras=boinc,static)!" << std::endl;
	return -1;
}

#endif
