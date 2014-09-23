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
#include <protocols/viewer/SilentObserver.hh>

// Package headers
// AUTO-REMOVED #include <protocols/viewer/viewers.hh>

#include <core/types.hh>

#include <core/io/silent/SilentFileData.hh>
//#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>

#include <ObjexxFCL/string.functions.hh>

#include <core/id/types.hh>
#include <utility/vector1.hh>

using ObjexxFCL::string_of;


// Project headers

// C++ Headers


namespace protocols {
namespace viewer {


/// @brief default constructor
SilentObserver::SilentObserver() {}

/// @brief constructor
SilentObserver::SilentObserver( std::string const & name_in, bool fullatom = false ) :
	frame_count_( 0 ),
	fullatom_( fullatom ),
	silent_file_name_( name_in + ".out" )
{
	sfd_ = core::io::silent::SilentFileDataOP( new core::io::silent::SilentFileData() );
}


/// @brief default destructor
SilentObserver::~SilentObserver() {
	detach_from();
}


/// @brief attach to a Pose
void
SilentObserver::attach_to( Pose & pose ) {
	detach_from();
	energy_event_link_ = pose.attach_energy_obs( &SilentObserver::on_energy_change, this );
}


/// @grief detach from Pose
void
SilentObserver::detach_from() {
	energy_event_link_.invalidate();
}


/// @brief upon receiving an EnergyEvent write to silent file
void
SilentObserver::on_energy_change(
	core::pose::signals::EnergyEvent const & event
) {
	frame_count_++;
	// std::cout << "frame(" << frame_count_ << ")" << std::endl;
	//core::io::silent::ProteinSilentStruct pss( *event.pose, "frame" + string_of(frame_count_)  );
	using namespace core::io::silent;
	SilentStructOP ss = SilentStructFactory::get_instance()->get_silent_struct_out();
	ss->fill_struct( *event.pose, "frame" + string_of(frame_count_) );

	sfd_->write_silent_struct( *ss, silent_file_name_ );
}


} // namespace conformation
} // namespace core


