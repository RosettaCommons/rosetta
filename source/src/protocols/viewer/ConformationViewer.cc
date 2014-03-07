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

#ifdef GL_GRAPHICS //////////////////////////////////////////////////////////////////---------------------

// Unit headers
#include <protocols/viewer/ConformationViewer.hh>

// Package headers
#include <protocols/viewer/viewers.hh>

#include <core/chemical/rna/RNA_Util.hh> // for silly centering of RNA.
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
#ifdef MAC
#include <GLUT/glut.h>
#else
#include "GL/glut.h"
#endif

#include <pthread.h>

using namespace core; ////////////////////////////////// !!!!!!!!!!!!!!!!!!!!!!!!!! DANGER



namespace protocols {
namespace viewer {


// signal that our window hasn't been created yet
int const BAD_WINDOW( -999 );

/// @brief default constructor
ConformationViewer::ConformationViewer() :
	Super(),
	new_conformation_( true ),
	use_debug_pause_( false ),
	conf_( NULL )
{}

// CALLED BY WORKER THREAD
ConformationViewer::ConformationViewer(std::string const & name_in ):
	Super(),
	name_( name_in ),
	new_conformation_( true ),
	my_window_( BAD_WINDOW ),
	length_( 900 ),
	width_( 900 ),
	use_debug_pause_( false ),
	conf_( NULL )
{
	pthread_mutex_init( &residues_mut_, NULL );
}

ConformationViewer::ConformationViewer(std::string const & name_in, int length, int width, bool debug_pause ):
	Super(),
	name_( name_in ),
	new_conformation_( true ),
	my_window_( BAD_WINDOW ),
	length_( length ),
	width_( width ),
	use_debug_pause_( debug_pause ),
	conf_( NULL )
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

	Vector center_vector = residues_[ anchor_id_.rsd() ]->xyz( anchor_id_.atomno() );
	//	if ( true ) center_vector = get_center( residues_ );

	//display_residues( residues_, anchor_id_ ); // in viewer.cc
	if ( basic::options::option[ basic::options::OptionKeys::edensity::mapfile ].user() && false) {

		const core::scoring::electron_density::ElectronDensity& edm = core::scoring::electron_density::getDensityMap();
		runtime_assert ( edm.isMapLoaded() );
		draw_conformation_and_density( residues_, secstruct_, triangles_ , current_gs_, center_vector );

	} else {
		draw_conformation( residues_, secstruct_, current_gs_, center_vector ); // in viewer.cc
	}

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
		glutSetWindow( my_window_ );
		glutPostRedisplay();

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
		residues_[ i ] = res_caps[ i ].get();
		secstruct_[ i ] = event.conformation->secstruct( i );
	}

	core::kinematics::tree::AtomCOP root_atom = event.conformation->atom_tree().root();
	if ( root_atom ) {
		anchor_id_ = root_atom->id();
		core::conformation::ResidueCOP rsd_root = residues_[ anchor_id_.rsd() ];
		if ( rsd_root->is_RNA() ) {
			anchor_id_ = id::AtomID( rsd_root->atom_index( core::chemical::rna::default_jump_atom( *rsd_root ) ), anchor_id_.rsd() );
		}
	}

	new_conformation_ = true;

	// always set the debug pause in case it gets reset elsewhere
	event.conformation->debug_pause( use_debug_pause_ );

	pthread_mutex_unlock( &residues_mut_ );
}


} // namespace conformation
} // namespace core


#endif // ifdef GL_GRAPHICS
