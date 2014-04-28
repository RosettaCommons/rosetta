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
#include <utility/basic_sys_util.hh>
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

#include <ObjexxFCL/string.functions.hh>

#ifndef WIN32
#include <unistd.h>
#else
#include <io.h>
#endif

#include <cstdio>

namespace basic {

CSI_Sequence::CSI_Sequence(std::string sequence_)
{
#ifdef WIN32
	if( _isatty(fileno(stdout)) ) sequence = sequence_;
#else
	if( isatty(fileno(stdout)) ) sequence = sequence_;
#endif
}

Tracer::OstreamPointer &Tracer::final_stream()
{
	static std::ostream * final_stream_ = &std::cout;
	return final_stream_;
}

void Tracer::set_new_final_stream(std::ostream *new_final_stream)
{
	final_stream() = new_final_stream;
}

void Tracer::set_default_final_stream()
{
	final_stream() = &std::cout;
}


otstreamOP &Tracer::ios_hook()
{
	static otstreamOP hook;
	return hook;
}

bool &Tracer::ios_hook_raw_()	// uninitilized, we will set correct output during set_ios_hook(...)
{
	static bool raw = true;
	return raw;
}

utility::vector1<std::string> Tracer::monitoring_list_;


//std::string const Tracer::AllChannels("_Really_Unique_String_Object_To_Identify_All_Tracer_Channels__qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM");
std::string const & Tracer::get_all_channels_string() {
	static std::string const all_channels("_Really_Unique_String_Object_To_Identify_All_Tracer_Channels__qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM");
	return all_channels;
}


TracerOptions Tracer::tracer_options_;

bool & Tracer::super_mute_()
{
	static bool mute = false;
	return mute;
}

int Tracer::mpi_rank_( 0 );


/// @details Static objects to hold various ASCII CSI codes, treat them as almost-const's
///          Codes below is all Hogwarts-approved magic numbers, so do not modify them. For reference see: http://en.wikipedia.org/wiki/ANSI_escape_code#CSI_codes
CSI_Sequence Tracer::Reset("\x1b[0m"), Tracer::Bold("\x1b[1m"), Tracer::Underline("\x1b[4m"),
	  Tracer::Black("\x1b[30m"),   Tracer::Red("\x1b[31m"),   Tracer::Green("\x1b[32m"),   Tracer::Yellow("\x1b[33m"),   Tracer::Blue("\x1b[34m"),   Tracer::Magenta("\x1b[35m"),   Tracer::Cyan("\x1b[36m"),   Tracer::White("\x1b[37m"),
	Tracer::bgBlack("\x1b[40m"), Tracer::bgRed("\x1b[41m"), Tracer::bgGreen("\x1b[42m"), Tracer::bgYellow("\x1b[43m"), Tracer::bgBlue("\x1b[44m"), Tracer::bgMagenta("\x1b[45m"), Tracer::bgCyan("\x1b[46m"), Tracer::bgWhite("\x1b[47m");

/// return visibility of TracerProxy object.
bool Tracer::TracerProxy::visible() const
{
	if( !visibility_calculated_ ) {
		bool muted;  int mute_level;
		calculate_visibility(channel_, priority_, visible_, muted, mute_level, tracer_.muted_by_default_);
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
	/// Do nothing here - contents will get flushed in Tracer destructor.
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
	mute_level_(-1), visible_(true), muted_( false ), muted_by_default_(muted_by_default),
	begining_of_the_line_(true), visibility_calculated_(false)
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
	mute_level_ = -1;
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

	if ( muted_ ) return false;
	if ( priority > mute_level_ ) return false;
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
		/*
		#ifdef EXPERIMENTAL_TRACER_FEATURES
			muted_ = priority >= mute_level_;
		#endif // EXPERIMENTAL_TRACER_FEATURES
		*/
		//visible_ = !muted_ && ( priority <= tracer_options_.level );
		visible_ = !muted_ && ( priority <= mute_level_ );
	}
}


/// @details Calculate visibility of current Tracer object.
void Tracer::calculate_visibility() const
{
	calculate_visibility(channel_, priority_, visible_, muted_, mute_level_, muted_by_default_);
	visibility_calculated_ = true;

}


/// @details Calculate visibility (static version) of current Tracer object using channel name and priority.
/// result stored in 'muted' and 'visible'.
void Tracer::calculate_visibility(std::string const &channel, int priority, bool &visible, bool &muted, int &mute_level_, bool muted_by_default)
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

	mute_level_ = tracer_options_.level;
	calculate_tracer_level(tracer_options_.levels, channel, false, mute_level_);
	//std::cout << "levels:" << tracer_options_.levels <<" ch:" << channel << " mute_level:" << mute_level_ << " priority:" << priority << std::endl;

	if( priority > mute_level_ ) visible = false;
}


/// @details Check if string representing channel 'ch' is in vector<string> v. Return true if channel
/// is in vector, false otherwise.
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

/// Same as before but return integer value for matched channel or closest match (we asume that 'v' in levels format, ie like: <channel name>:level )
bool Tracer::calculate_tracer_level(utility::vector1<std::string> const & v, std::string const ch, bool strict, int &res)
{
	unsigned int len = 0;
	bool math = false;
	//std::cout << "Entring:calculate_tracer_level: v=" << v << " ch:" << ch << std::endl;
	for(size_t i=1; i<=v.size(); i++) {
		bool flag = false;
		utility::vector1< std::string > spl = utility::string_split(v[i], ':');
		//std::cout << "Split:" << spl << " size:" << spl[1].size() << std::endl;

		if( spl[1] == "all" && len == 0 ) flag = true;  // we can asume that 'all' is shorter then any valid core/protocol path... but we don't!

		if( spl[1] == ch ) flag=true;

		if( !strict ) {
			if( ch.size() > spl[1].size() ) {
				std::string s(ch);  s.resize(spl[1].size());
				if( s == spl[1] ) flag=true;
			}
		}
		if(flag  && ( len < spl[1].size() ) ) {
			math = true;
			len = spl[1].size();
			res = utility::string2int(spl[2]);
			std::string const spl2_lower = ObjexxFCL::lowercased( spl[2] );
			if( spl2_lower == "fatal" )   res = t_fatal;
			if( spl2_lower == "error"   || spl2_lower == "errors" )   res = t_error;
			if( spl2_lower == "warning" || spl2_lower == "warnings" ) res = t_warning;
			if( spl2_lower == "info" )    res = t_info;
			if( spl2_lower == "debug" )   res = t_debug;
			if( spl2_lower == "trace" )   res = t_trace;
			if( res == -1 ) {
				std::string message( "WARNING:: The setting '" + spl[2] + "' is not recognized as a valid tracer level." );
				if( ios_hook() ) {
					*ios_hook() << message << std::endl;
				} else if ( final_stream() ) {
					*final_stream() << message << std::endl;
				}
				res = t_info; // Set such that you get standard amount of output instead of none.
			}
			//std::cout << "Match:" << spl << " ch:" << ch << " res:"<< res << std::endl;
		}
		else {
			//std::cout << "Fail:" << spl << " ch:" << ch << " res:"<< res << std::endl;
		}

	}
	//std::cout << "Leaving:calculate_tracer_level: v=" << v << " ch:" << ch << " match:" << math <<" res:"<< res << std::endl;
	return math;
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
				if ( tracer_options_.timestamp ) {
					sout << utility::timestamp() << " ";
				}
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

	if( ios_hook() && ios_hook().get()!=this &&
		( in(monitoring_list_, channel_, false) || in(monitoring_list_, get_all_channels_string(), true ) ) ) {
		if (ios_hook_raw_() || visible() ){
			prepend_channel_name<otstream>( *ios_hook(), str );
			ios_hook()->flush();
		}
	}

	if ( !super_mute_() && visible() ){
		prepend_channel_name<std::ostream>( *final_stream(), str );
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
// void Tracer::set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list)
// {
// 	ios_hook_ = tr;
// 	monitoring_list_ = utility::split(monitoring_channels_list);
// 	ios_hook_raw_ = true;
// }
/// When raw==false same as above above but gives the option get only the
/// visible and unmuted tracers.  It can be useful to get the raw
/// output for applications like the comparing tracers, where the
/// output should not change with command line parameters.  It can be
/// useful to get the non-raw output in applications like using the
/// jd2 with MPI, where the output each job should match the output if
/// it was run on a single processor.
void Tracer::set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list, bool raw)
{
	ios_hook() = tr;
	monitoring_list_ = utility::split(monitoring_channels_list);
	ios_hook_raw_() = raw;
}

void PyTracer::t_flush(std::string const &str)
{
	buf_ += str;
	output_callback(str);
}


} // namespace basic
