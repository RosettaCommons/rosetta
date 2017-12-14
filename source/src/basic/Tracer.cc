// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/basic/unti/tracer.cc
/// @brief  Tracer IO system
/// @author Sergey Lyskov
/// @author Rocco Moretti (rmorettiase@gmail.com)


#include <basic/Tracer.hh>

#include <utility>
#include <utility/thread/backwards_thread_local.hh> //For THREAD_LOCAL

namespace basic {

static basic::Tracer TR_Error( "Error", basic::t_error );

Tracer & Error() {
	return TR_Error;
}

static basic::Tracer TR_Warning( "Warning", basic::t_warning );

Tracer & Warning() {
	return TR_Warning;
}

///////////////////////  Tracer //////////////////////////////////////

utility::CSI::CSI_Enum Tracer::Reset(utility::CSI::Reset), Tracer::Bold(utility::CSI::Bold), Tracer::Underline(utility::CSI::Underline),
Tracer::Black(utility::CSI::Black),   Tracer::Red(utility::CSI::Red),   Tracer::Green(utility::CSI::Green),   Tracer::Yellow(utility::CSI::Yellow),   Tracer::Blue(utility::CSI::Blue),   Tracer::Magenta(utility::CSI::Magenta),   Tracer::Cyan(utility::CSI::Cyan),   Tracer::White(utility::CSI::White),
Tracer::bgBlack(utility::CSI::bgBlack), Tracer::bgRed(utility::CSI::bgRed), Tracer::bgGreen(utility::CSI::bgGreen), Tracer::bgYellow(utility::CSI::bgYellow), Tracer::bgBlue(utility::CSI::bgBlue), Tracer::bgMagenta(utility::CSI::bgMagenta), Tracer::bgCyan(utility::CSI::bgCyan), Tracer::bgWhite(utility::CSI::bgWhite);

Tracer::Tracer(
	std::string const & channel,
	TracerPriority priority,
	bool muted_by_default
):
	Fatal(   *this, t_fatal   ),
	Error(   *this, t_error   ),
	Warning( *this, t_warning ),
	Info(    *this, t_info    ),
	Debug(   *this, t_debug   ),
	Trace(   *this, t_trace   ),
	channel_( channel ),
	channel_color_( utility::CSI::Nothing ),
	channel_name_color_( utility::CSI::Nothing ),
	priority_( priority ),
	muted_by_default_( muted_by_default )
{}

Tracer::Tracer(
	std::string const & channel,
	utility::CSI::CSI_Enum const & channel_color,
	utility::CSI::CSI_Enum const & channel_name_color,
	TracerPriority priority,
	bool muted_by_default
):
	Fatal(   *this, t_fatal   ),
	Error(   *this, t_error   ),
	Warning( *this, t_warning ),
	Info(    *this, t_info    ),
	Debug(   *this, t_debug   ),
	Trace(   *this, t_trace   ),
	channel_( channel ),
	channel_color_( channel_color ),
	channel_name_color_( channel_name_color ),
	priority_( priority ),
	muted_by_default_( muted_by_default )
{}

/// @details Ideally, we'd clean ourselves up from the tracer_impl_map,
/// but that's complicated by it being a) thread-dependent, and
/// b) subject to static initialization order fiasco (static destruction order fiasco) issues.
/// Regardless, the associated TracerImpl will be flushed when the static tracer_impl_map gets destroyed.
/// The one tricky bit is that the TracerImpl for non-static tracers stick around, and may accidentally
/// be re-used if we re-use the same pointer ... so avoid non-static `Tracer` objects: use a TracerImpl directly instead.
Tracer::~Tracer() = default;

TracerImpl &
Tracer::tracer_impl() {
	// As the map is thread_local, we don't have to worry about locks here, as each thread gets it's own.
	static THREAD_LOCAL std::unordered_map< Tracer*, std::unique_ptr< TracerImpl > > impl_map;

	std::unique_ptr< TracerImpl > & impl( impl_map[ this ] ); // Will initialize as nullptr if it doesn't exist
	if ( impl == nullptr ) {
		impl = create_impl();
	}
	return *impl;
}

std::unique_ptr< TracerImpl >
Tracer::create_impl() {
	return std::unique_ptr< TracerImpl >( new TracerImpl(channel_, utility::CSI_Sequence(channel_color_), utility::CSI_Sequence(channel_name_color_), priority_, muted_by_default_) );
}

///////////////////////  PyTracer //////////////////////////////////////

void PyTracer::t_flush(std::string const &str)
{
	buf_ += str;
	output_callback(str);
}


} // namespace basic
