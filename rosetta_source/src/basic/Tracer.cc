// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/basic/unti/tracer.cc
/// @brief  Tracer IO system
/// @author Sergey Lyskov


#ifdef USEMPI
#include <mpi.h>
#endif


#include <basic/Tracer.hh>

//#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>

#include <utility/string_util.hh>
#include <utility/tools/make_vector.hh>

#include <iostream>
#include <algorithm>

#include <platform/types.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <cassert>
#include <cstddef>
#include <iosfwd>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>
#include <basic/Tracer.fwd.hh>
#include <boost/algorithm/string/erase.hpp>


namespace basic {

std::ostream *final_channel = &std::cout;

otstreamOP Tracer::ios_hook_;
bool Tracer::ios_hook_raw_;
utility::vector1<std::string> Tracer::monitoring_list_;

std::string const Tracer::AllChannels("_Really_Unique_String_Object_To_Identify_All_Tracer_Channels__qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM");

TracerOptions Tracer::tracer_options_;
bool Tracer::super_mute_(false);
int Tracer::mpi_rank_( 0 );


/// return visibility of TracerProxy object.
bool Tracer::TracerProxy::visible() const
{
	if( !visibility_calculated_ ) {

		#ifndef EXPERIMENTAL_TRACER_FEATURES
			bool visible, muted;
			calculate_visibility(channel_, priority_, visible, muted, tracer_.muted_by_default_);
		#else
			int mute_level;
			experimental_calculate_visibility(channel_, priority_, visible_, mute_level, tracer_.muted_by_default_);
		#endif // EXPERIMENTAL_TRACER_FEATURES

	}
	return visible_;
}


/// Flush inner buffer: send it to bound Tracer object, and clean it.
void Tracer::TracerProxy::t_flush(std::string const & s)
{
	int pr = tracer_.priority();
	tracer_.priority(priority_);
	tracer_ << s;
	tracer_.flush();
	tracer_.priority(pr);
}


Tracer::TracerProxy::~TracerProxy()
{
	/// Do nothing here - contents get flushed in Tracer destructor.
}


/// @details static collection of all Tracer objects
std::vector< Tracer * > & Tracer::all_tracers()
{
	static std::vector< Tracer * > * allTr = new std::vector< Tracer * >;
	return *allTr;
}

void Tracer::flush_all_tracers()
{
	for(std::vector<Tracer *>::iterator it=all_tracers().begin(); it < all_tracers().end(); ++it) {
		(*it)->flush_all_channels();
	}
}


/// @details Constructor of Tracer object.
/// Since most of the Tracer object will be created as static - they Constuctor will be called before
/// Option system is initialized. So we can't really calculate any vizibility or priority here.
/// Such calculation should be done later, whe first IO operation happend.
/// @todo Default Tracer level should probably be modified to t_info here and in options defn.
Tracer::Tracer(std::string const & channel, TracerPriority priority, bool muted_by_default) :
	Fatal(*this,   t_fatal, channel),
	Error(*this,   t_error, channel),
	Warning(*this, t_warning, channel),
	Info(*this,    t_info, channel),
	Debug(*this,   t_debug, channel),
	Trace(*this,   t_trace, channel),
	visible_(true), muted_( false ), muted_by_default_(muted_by_default), begining_of_the_line_(true),
	visibility_calculated_(false)
{
	channel_ = channel;
	priority_ = priority;

	all_tracers().push_back(this);
}

Tracer::~Tracer()
{
	/// We could not gurantee that option system was not deleted already, so we must disable
	/// options access. However we don't want channel to withhold any output since it can be important.
    /// We check if anything still present in channel buffer, and if it is - print it contents with warning.

	std::vector< otstream* > v = utility::tools::make_vector< otstream* >(this, &Fatal, &Error, &Warning,
															 &Info, &Debug, &Trace);

	//bool need_flush = false;
	for(size_t i=0; i<v.size(); i++) {
		if( !v[i]->is_flushed() ) {
			//v[i]->flush();
			(*v[i]) << std::endl;
			(*v[i]) << "WARNING: Message(s) above was printed in the end instead of proper place because this Tracer object has some contents left in inner buffer when destructor was called. Explicit call Tracer::flush() or end your IO with std::endl to disable this warning.\n" << std::endl;
		}
	}

	for(int i = all_tracers().size()-1; i >= 0; --i) {
		if( this == all_tracers()[i] ) all_tracers().erase( all_tracers().begin()+i);
	}

	/* for(std::vector< Tracer * >::iterator it=all_tracers().begin(); it < all_tracers().end(); ++it) {
		if( this == *it ) { all_tracers().erase( it ); break; }
	} */
}


///  @details re-init using data from another tracer object.
void Tracer::init(Tracer const & tr)
{
	channel_ = tr.channel_;
	priority_ = tr.priority_;

	visible_ = true;
	muted_ = false;
	begining_of_the_line_ = true;
	visibility_calculated_ = false;
}

void Tracer::flush_all_channels()
{
	std::vector< otstream* > v = utility::tools::make_vector< otstream* >(this, &Fatal, &Error, &Warning,
															 &Info, &Debug, &Trace);

	for(size_t i=0; i<v.size(); i++) {
		v[i]->flush();
	}
}

bool Tracer::visible( int priority ) const {
	if (!visibility_calculated_) calculate_visibility();

	#ifdef EXPERIMENTAL_TRACER_FEATURES
		muted_ = priority >= mute_level_;
	#endif // EXPERIMENTAL_TRACER_FEATURES

	if ( muted_ ) return false;
	if ( priority > tracer_options_.level ) return false;
	return true;
}

bool Tracer::visible() const {
	if (!visibility_calculated_) calculate_visibility();
	return visible_;
}

Tracer & Tracer::operator ()(int priority)
{
	this->priority(priority);
	return *this;
}

void Tracer::priority(int priority)
{
	priority_ = priority;
	if (visibility_calculated_) {
		#ifdef EXPERIMENTAL_TRACER_FEATURES
			muted_ = priority >= mute_level_;
		#endif // EXPERIMENTAL_TRACER_FEATURES

		visible_ = !muted_ && ( priority <= tracer_options_.level );
	}
}


/// @details Calculate visibility of current Tracer object.
void Tracer::calculate_visibility() const
{
	#ifndef EXPERIMENTAL_TRACER_FEATURES
		calculate_visibility(channel_, priority_, visible_, muted_, muted_by_default_);
	#else
		experimental_calculate_visibility(channel_, priority_, visible_, mute_level_, muted_by_default_);
	#endif // EXPERIMENTAL_TRACER_FEATURES

	visibility_calculated_ = true;

}


/// @details Calculate visibility (static version) of current Tracer object using channel name and priority.
/// result stored in 'muted' and 'visible'.
void Tracer::calculate_visibility(std::string const &channel, int priority, bool &visible, bool &muted, bool muted_by_default)
{
	visible = false;
	if( in(tracer_options_.muted, "all", true) ) {
		if( in(tracer_options_.unmuted, channel, false) ) visible = true;
		else visible = false;
	}
	else {
		if( in(tracer_options_.unmuted, "all", true) ) {
			if( in(tracer_options_.muted, channel, false) ) visible = false;
			else visible = true;
		}
		else {  /// default bechavior: unmute unless muted_by_default is true
			if( muted_by_default ) {
				if( in(tracer_options_.unmuted, channel, false) ) visible = true;
				else visible = false;
			}
			else {
				if( in(tracer_options_.muted, channel, false) ) visible = false;
				else visible = true;
			}
		}
	}

	//if we are in MPI mode --- most of the time one doesn't want to see output from all nodes just the master and 1st client is plenty ..
#ifdef USEMPI
	int already_initialized = 0;
	int already_finalized = 0;
	MPI_Initialized( &already_initialized );
	MPI_Finalized( &already_finalized );
	if ( already_initialized != 0 && already_finalized == 0 ) {
		int mpi_rank, mpi_nprocs;
		MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);/* get current process id */
		MPI_Comm_size (MPI_COMM_WORLD, &mpi_nprocs);/* get number of processes */
		mpi_rank_=mpi_rank;
	}

	if ( in(tracer_options_.muted, "all_high_mpi_rank", true ) ) {
		if ( mpi_rank_>=2 ) visible = false; //* visible: master and 1st client: rank 0 and rank1
	}

	if ( in(tracer_options_.muted, "all_high_mpi_rank_filebuf", true ) ) {
		if ( mpi_rank_>=4 ) visible = false; //* visible: master, filebuf and 1st client: rank 0, 1, 2
	}

#endif
	muted = !visible;
	if( priority > tracer_options_.level ) visible = false;
}


/// @details Check if string representing channel 'ch' is in vector<string> v. Return true if channel
/// is in vector, false other wise.
/// Two mode of operation:
/// Strict:  strict==true - channels compared verbatim.
/// Regular: strict==false - comparing with hierarchy in mind,
///                          ie: ch='basic.pose' v[0]='basic' --> will yield true.
bool Tracer::in(utility::vector1<std::string> const & v, std::string const ch, bool strict)
{
	for(size_t i=1; i<=v.size(); i++) {
		if( v[i] == ch ) return true;

		if( !strict ) {
			if( ch.size() > v[i].size() ) {
				std::string s(ch);  s.resize(v[i].size());
				if( s == v[i] ) return true;
			}
		}
	}
	return false;
}

/// @dtails Write the contents of str to sout prepending the channel
/// name on each line if the print_channel_name flag is set.
template <class out_stream>
void Tracer::prepend_channel_name( out_stream & sout, std::string const &str){
	if ( tracer_options_.print_channel_name ){
		std::string s = str;
		begining_of_the_line_ = true;
		for (size_t i=0; i<s.size(); i++) {
			if ( begining_of_the_line_ ) {
				sout << channel_ << ": ";
#ifdef USEMPI
				sout << "(" << mpi_rank_ << ") ";
#endif
				begining_of_the_line_ = false;
			}
			sout << s[i];
		}
	} else {
		sout << str;
	}
}




/// @details Inform Tracer that is contents was modified, and IO is in order.
void Tracer::t_flush(std::string const &str)
{
	if( !visibility_calculated_ ) calculate_visibility();

	if( ios_hook_ && ios_hook_.get()!=this &&
			( in(monitoring_list_, channel_, false) || in(monitoring_list_, AllChannels, true ) ) ) {
		if (ios_hook_raw_ || visible_){
			prepend_channel_name<otstream>( *ios_hook_, str );
			ios_hook_->flush();
		}
	}

	if ( !super_mute_ && visible_ ){
		prepend_channel_name<std::ostream>( *final_channel, str );
	}
}


/// @details Return reference to static Tracer object (after setting it channel and priority).
Tracer & T(std::string const & channel, TracerPriority priority)
{
	static Tracer * t = new Tracer();
	t->channel_ = channel;
	t->priority_ = priority;
	t->calculate_visibility();
	t->begining_of_the_line_ = true;
	return *t;
}




/// @details Set OStringStream object to which all Tracers output
/// listed in the monitoring_channels_list should be copied.  Note
/// this copies the output of channels even if they are invisible or
/// muted.
void Tracer::set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list)
{
	ios_hook_ = tr;
	monitoring_list_ = utility::split(monitoring_channels_list);
	ios_hook_raw_ = true;
}

/// @details Same as above above but gives the option get only the
/// visible and unmuted tracers.  It can be useful to get the raw
/// output for applications like the comparing tracers, where the
/// output should not change with command line parameters.  It can be
/// useful to get the non-raw output in applications like using the
/// jd2 with MPI, where the output each job should match the output if
/// it was run on a single processor.
	void Tracer::set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list, bool raw)
{
	ios_hook_ = tr;
	monitoring_list_ = utility::split(monitoring_channels_list);
	ios_hook_raw_ = raw;
}


#ifdef EXPERIMENTAL_TRACER_FEATURES
	/// @details Calculate visibility (static version) of current Tracer object using channel name and priority.
	/// result stored in 'muted' and 'visible'.
	void Tracer::experimental_calculate_visibility(std::string const &channel, int priority, bool &visible, int &mute_level, bool muted_by_default)
	{
			visible = false;
			if( in(tracer_options_.muted, "all", true) ) {
					if( in(tracer_options_.unmuted, channel, false) ) visible = true;
					else visible = false;
			}
			else {
					if( in(tracer_options_.unmuted, "all", true) ) {
							if( in(tracer_options_.muted, channel, false) ) visible = false;
							else visible = true;
					}
					else {  /// default bechavior: unmute unless muted_by_default is true
							if( muted_by_default ) {
									if( in(tracer_options_.unmuted, channel, false) ) visible = true;
									else visible = false;
							}
							else {
									if( in(tracer_options_.muted, channel, false) ) visible = false;
									else visible = true;
							}
					}
			}


			//      muted = !visible;
			//comes at end now...   if( priority > tracer_options_.level ) visible = false;
	  if ( !visible ) mute_level = t_fatal;
			else {
					mute_level = 1000;
					if ( in(tracer_options_.muted_trace, "all", true) ) {
							mute_level = t_trace;
					}
					if ( in(tracer_options_.muted_trace, channel, false ) ) {
							mute_level = t_trace;
					}
					if ( in(tracer_options_.muted_debug, "all", true) ) {
							mute_level = t_debug;
					}
					if ( in(tracer_options_.muted_debug, channel, false ) ) {
							mute_level = t_debug;
					}

					if ( in(tracer_options_.muted_info, "all", true) ) {
							mute_level = t_info;
					}
					if ( in(tracer_options_.muted_info, channel, false ) ) {
							mute_level = t_info;
					}

					if ( in(tracer_options_.muted_warning, "all", true) ) {
							mute_level = t_warning;
					}
					if ( in(tracer_options_.muted_warning, channel, false ) ) {
							mute_level = t_warning;
					}
			}
			if ( mute_level < t_warning ) {
					if ( in(tracer_options_.unmuted_error, "all", true) ) {
							mute_level = t_warning;
					}
					if ( in(tracer_options_.unmuted_error, channel, false ) ) {
							mute_level = t_warning;
					}
			}

			if ( mute_level < t_info ) {
					if ( in(tracer_options_.unmuted_warning, "all", true) ) {
							mute_level = t_info;
					}
					if ( in(tracer_options_.unmuted_warning, channel, false ) ) {
							mute_level = t_info;
					}
			}

			if ( mute_level < t_debug ) {
					if ( in(tracer_options_.unmuted_info, "all", true) ) {
							mute_level = t_debug;
					}
					if ( in(tracer_options_.unmuted_info, channel, false ) ) {
							mute_level = t_debug;
					}

			}

			if ( mute_level < t_trace ) {
					if ( in(tracer_options_.unmuted_debug, "all", true) ) {
							mute_level = t_trace;
					}
					if ( in(tracer_options_.unmuted_debug, channel, false ) ) {
							mute_level = t_trace;
					}
			}

	//if we are in MPI mode --- most of the time one doesn't want to see output from all nodes just the master and 1st client is plenty ..
	#ifdef USEMPI
			int already_initialized = 0;
			int already_finalized = 0;
			MPI_Initialized( &already_initialized );
			MPI_Finalized( &already_finalized );
			if ( already_initialized != 0 && already_finalized == 0 ) {
					int mpi_rank, mpi_nprocs;
					MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);/* get current process id */
					MPI_Comm_size (MPI_COMM_WORLD, &mpi_nprocs);/* get number of processes */
					mpi_rank_=mpi_rank;
			}

			if ( in(tracer_options_.muted, "all_high_mpi_rank", true ) ) {
					if ( mpi_rank_>=2 ) {
							mute_level = std::min( (int)t_warning, mute_level );
					}
			}

			if ( in(tracer_options_.muted, "all_high_mpi_rank_filebuf", true ) ) {
					if ( mpi_rank_>=4 ) {
							visible = false; //* visible: master, filebuf and 1st client: rank 0, 1, 2
							mute_level = std::min( (int)t_warning, mute_level );
					}
			}
	#endif

			visible = (priority < mute_level) && ( priority <= tracer_options_.level );
	}

#endif // EXPERIMENTAL_TRACER_FEATURES



} // namespace basic

