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

#ifdef GL_GRAPHICS

#ifndef INCLUDED_protocols_viewer_ConformationViewer_hh
#define INCLUDED_protocols_viewer_ConformationViewer_hh


// Unit headers
#include <protocols/viewer/ConformationViewer.fwd.hh>

// Package headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.hh>

// Project headers

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/signals/Link.hh>

#include <protocols/viewer/GraphicsState.hh>
#include <protocols/viewer/triangle.hh>

#include <utility/vector1_bool.hh>

//Auto Headers


// C++ Headers

#include <pthread.h>

namespace protocols {
namespace viewer {

/// @brief observer that attaches to a Conformation and displays graphics
class ConformationViewer : public utility::pointer::ReferenceCount {

private: // typedefs

	typedef utility::pointer::ReferenceCount Super;

public: // typedefs

	typedef utility::vector1< core::conformation::ResidueCOP > ResidueCOPs;

public: // construct/destruct

	/// @brief default constructor
	ConformationViewer();

	/// @brief constructor
	ConformationViewer( std::string const & name );

	/// @brief constructor
	ConformationViewer( std::string const & name, int length, int width, bool debug_pause );

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
	utility::vector1< triangle > triangles_;

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

	pthread_mutex_t residues_mut_;
};

} // viewer
} // protocols


#endif

#endif 
