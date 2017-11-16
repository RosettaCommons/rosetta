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


#ifdef USEMPI
#include <mpi.h>
#endif

#include <basic/TracerImpl.hh>
#include <ObjexxFCL/string.functions.hh>  // for lowercased
#include <algorithm>                      // for find
#include <cassert>                        // for assert
#include <cstddef>                        // for size_t
#include <iosfwd>                         // for string, ostream
#include <iostream>                       // for cout, cerr
#include <ostream>                        // for operator<<, basic_ostream
#include <platform/types.hh>              // for Size
#include <string>                         // for allocator, operator==, basi...
#include <utility/CSI_Sequence.hh>        // for CSI_Sequence, CSI_Black
#include <utility/sys_util.hh>      // for timestamp
#include <utility/string_util.hh>         // for split, string2int, string_s...
#include <utility/tools/make_vector.hh>   // for make_vector
#include <utility/vector1.hh>             // for vector1
#include <utility/vectorL.hh>             // for vectorL
#include <vector>                         // for vector, vector<>::iterator

#ifndef WIN32
#include <unistd.h>
#else
#include <io.h>
#endif

#ifdef MULTI_THREADED

#include <mutex>

#endif


namespace basic {

bool
TracerOptions::operator==( TracerOptions const & other ) const {
	return level == other.level &&
		print_channel_name == other.print_channel_name &&
		timestamp == other.timestamp &&
		muted == other.muted &&
		unmuted == other.unmuted &&
		levels == other.levels;
}

bool
TracerOptions::operator!=( TracerOptions const & other ) const {
	return !(*this == other);
}

///////////////////////  TracerImpl //////////////////////////////////////

std::string TracerImpl::output_prefix_;


#ifdef MULTI_THREADED

/// @brief This mutex ensures that any time static data is read from or written to by
/// a tracer object, that the thread gets exclusive access to it.
/// See the TracerImpl datamembers for items which need to be mutex protected
std::mutex & tracer_static_data_mutex() {
	// Construct On First Use Idiom, to avoid static initialization order fiasco with static tracers.
	// In C++11, static locals are constructed only once, in a threadsafe manner.
	static std::mutex * mutex = new std::mutex;
	return *mutex;
}

#endif


TracerImpl::OstreamPointer &TracerImpl::final_stream()
{
	static std::ostream * final_stream_ = &std::cout;
	return final_stream_;
}

void TracerImpl::set_new_final_stream(std::ostream *new_final_stream)
{
	if ( dynamic_cast< TracerImpl* >(new_final_stream) != nullptr || dynamic_cast< TracerImpl::TracerProxyImpl* >(new_final_stream) != nullptr ) {
		utility_exit_with_message("Error: Setting the final_stream on a Tracer to be another Tracer or a TracerProxy is only going to end in grief! (Try a PyTracer instead.)");
	}
	final_stream() = new_final_stream;
}

void TracerImpl::set_default_final_stream()
{
	final_stream() = &std::cout;
}


otstreamOP &TracerImpl::ios_hook()
{
	static otstreamOP hook;
	return hook;
}

bool &TracerImpl::ios_hook_raw_()
{
	static bool raw = true;
	return raw;
}

utility::vector1< std::string > &
TracerImpl::monitoring_list_()
{
	static utility::vector1< std::string > monitoring_list;
	return monitoring_list;
}


//std::string const TracerImpl::AllChannels("_Really_Unique_String_Object_To_Identify_All_Tracer_Channels__qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM");
std::string const & TracerImpl::get_all_channels_string() {
	static std::string const all_channels("_Really_Unique_String_Object_To_Identify_All_Tracer_Channels__qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM");
	return all_channels;
}


TracerOptionsOP TracerImpl::tracer_options_;

void TracerImpl::set_tracer_options( TracerOptions const & to) {
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > lock( tracer_static_data_mutex() );
#endif

	if ( tracer_options_ && to != *tracer_options_ ) {
		// Unit tests reset the options a lot, but it's always with an identical value
		std::cout << "[ WARNING ] Resetting the global tracer options, which have already been set." << std::endl;
		std::cout << "[ WARNING ] This will not affect Tracers which are already initialized." << std::endl;
	}
	tracer_options_ = TracerOptionsOP( new TracerOptions( to ) );
}

bool & TracerImpl::super_mute_()
{
	static bool mute = false;
	return mute;
}

int TracerImpl::mpi_rank_( 0 );


/// @details Static objects holding various ASCII CSI codes (see utility/CSI_Sequence.hh)
utility::CSI_Sequence TracerImpl::Reset(utility::CSI_Reset()), TracerImpl::Bold(utility::CSI_Bold()), TracerImpl::Underline(utility::CSI_Underline()),
TracerImpl::Black(utility::CSI_Black()),   TracerImpl::Red(utility::CSI_Red()),   TracerImpl::Green(utility::CSI_Green()),   TracerImpl::Yellow(utility::CSI_Yellow()),   TracerImpl::Blue(utility::CSI_Blue()),   TracerImpl::Magenta(utility::CSI_Magenta()),   TracerImpl::Cyan(utility::CSI_Cyan()),   TracerImpl::White(utility::CSI_White()),
TracerImpl::bgBlack(utility::CSI_bgBlack()), TracerImpl::bgRed(utility::CSI_bgRed()), TracerImpl::bgGreen(utility::CSI_bgGreen()), TracerImpl::bgYellow(utility::CSI_bgYellow()), TracerImpl::bgBlue(utility::CSI_bgBlue()), TracerImpl::bgMagenta(utility::CSI_bgMagenta()), TracerImpl::bgCyan(utility::CSI_bgCyan()), TracerImpl::bgWhite(utility::CSI_bgWhite());

TracerImpl::TracerProxyImpl::TracerProxyImpl(
	TracerImpl & tracer,
	int priority,
	std::string const & channel
) :
	tracer_(tracer),
	priority_(priority),
	channel_(channel),
	visible_(true)
{}

void
TracerImpl::TracerProxyImpl::calculate_visibility()
{
	bool muted;  int mute_level;
	TracerImpl::calculate_visibility( channel_, priority_, visible_, muted, mute_level, tracer_.muted_by_default_ );
}

TracerImpl::TracerProxyImpl &
TracerImpl::get_proxy_by_priority( TracerPriority priority) {
	if ( priority <= t_fatal ) {
		return Fatal;
	} else if ( priority <= t_error ) {
		return Error;
	} else if ( priority <= t_warning ) {
		return Warning;
	} else if ( priority <= t_info ) {
		return Info;
	} else if ( priority <= t_debug ) {
		return Debug;
	} else {
		return Trace;
	}
}

/// Flush inner buffer: send it to bound Tracer object, and clean it.
void TracerImpl::TracerProxyImpl::t_flush(std::string const & s)
{
	int pr = tracer_.priority();
	tracer_.priority(priority_);
	tracer_ << s;
	tracer_.flush();
	tracer_.priority(pr);
}


TracerImpl::TracerProxyImpl::~TracerProxyImpl()
{
	/// Do nothing here - contents will get flushed in Tracer destructor.
}

/// @details Constructor of Tracer object.
/// Since most of the Tracer object will be created as static - they Constuctor will be called before
/// Option system is initialized. So we can't really calculate any vizibility or priority here.
/// Such calculation should be done later, whe first IO operation happend.
/// @todo Default Tracer level should probably be modified to t_info here and in options defn.
TracerImpl::TracerImpl(std::string const & channel, TracerPriority priority, bool muted_by_default) :
	Fatal(   *this, t_fatal,   channel ),
	Error(   *this, t_error,   channel ),
	Warning( *this, t_warning, channel ),
	Info(    *this, t_info,    channel ),
	Debug(   *this, t_debug,   channel ),
	Trace(   *this, t_trace,   channel )
{
	init(channel, utility::CSI_Nothing(), utility::CSI_Nothing(), priority, muted_by_default);
}


TracerImpl::TracerImpl(std::string const & channel, utility::CSI_Sequence const & channel_color, utility::CSI_Sequence const & channel_name_color, TracerPriority priority, bool muted_by_default) :
	Fatal(   *this, t_fatal,   channel ),
	Error(   *this, t_error,   channel ),
	Warning( *this, t_warning, channel ),
	Info(    *this, t_info,    channel ),
	Debug(   *this, t_debug,   channel ),
	Trace(   *this, t_trace,   channel )
{
	init(channel, channel_color, channel_name_color, priority, muted_by_default);
}

void TracerImpl::init(std::string const & channel, utility::CSI_Sequence const & channel_color, utility::CSI_Sequence const & channel_name_color, TracerPriority priority, bool muted_by_default)
{
	mute_level_ = -1;
	visible_ = true;
	muted_ = false;
	muted_by_default_ = muted_by_default;

	channel_ = channel;
	priority_ = priority;

	channel_color_ = channel_color;
	channel_name_color_ = channel_name_color;

	// TracerImpl should only be initialized after the TracerOptions object is initialized.
	calculate_visibility();
}

TracerImpl::~TracerImpl()
{
	/// We could not gurantee that option system was not deleted already, so we must disable
	/// options access. However we don't want channel to withhold any output since it can be important.
	/// We check if anything still present in channel buffer, and if it is - print it contents with warning.

	std::vector< otstream* > v = utility::tools::make_vector< otstream* >(
		this, &Fatal, &Error, &Warning,
		&Info, &Debug, &Trace);

	/// PyRosetta final stream could be redirected to some Python object which might be already got destroyed during Python exit
#ifndef PYROSETTA
	//set_new_final_stream( &std::cerr );
	//set_ios_hook(otstreamOP(), "");

	//bool need_flush = false;
	for ( size_t i=0; i<v.size(); i++ ) {
		if ( !v[i]->is_flushed() ) {
			//v[i]->flush();
			(*v[i]) << std::endl;
			(*v[i]) << "[ WARNING ] Message(s) above was printed in the end instead of proper place because this Tracer object has some contents left in inner buffer when destructor was called. Explicitly call Tracer::flush() or end your IO with std::endl to disable this warning.\n" << std::endl;
		}
	}
#endif
}


///  @details re-init using data from another tracer object.
void TracerImpl::init( TracerImpl const & tr )
{
	channel_ = tr.channel_;
	priority_ = tr.priority_;

	visible_ = true;
	muted_ = false;
	mute_level_ = -1;

	// TracerImpl is only likely to be called after the options are set
	calculate_visibility();
}

/// @brief set/get globale string-prefix for all Tracer output strings
void TracerImpl::output_prefix(std::string const &new_output_prefix)
{
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > lock( tracer_static_data_mutex() );
#endif

	output_prefix_ = new_output_prefix;
}

std::string TracerImpl::output_prefix()
{
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > lock( tracer_static_data_mutex() );
#endif
	return output_prefix_;
}


void TracerImpl::flush_all_channels()
{
	std::vector< otstream* > v = utility::tools::make_vector< otstream* >(
		this, &Fatal, &Error, &Warning,
		&Info, &Debug, &Trace);

	for ( size_t i=0; i<v.size(); i++ ) {
		v[i]->flush();
	}
}

bool TracerImpl::visible( int priority ) {
	if ( muted_ ) return false;
	if ( priority > mute_level_ ) return false;
	return true;
}

TracerImpl &
TracerImpl::operator ()(int priority)
{
	this->priority(priority);
	return *this;
}

void TracerImpl::priority(int priority)
{
	priority_ = priority;
	visible_ = !muted_ && ( priority <= mute_level_ );
}


/// @details Calculate visibility of current Tracer object and all of its proxies
void TracerImpl::calculate_visibility()
{
	bool options_set = false; // Indirection is so we don't have to lock around the sub-calculate_visibility calls
	{
#ifdef MULTI_THREADED
		std::lock_guard< std::mutex > lock( tracer_static_data_mutex() );
#endif
		options_set = tracer_options_ != nullptr;
	}
	if ( options_set ) { // Small chance of race conditions with resetting to nullptr, but that's handled with asserts below.
		calculate_visibility(channel_, priority_, visible_, muted_, mute_level_, muted_by_default_);
		Fatal.calculate_visibility();
		Error.calculate_visibility();
		Warning.calculate_visibility();
		Info.calculate_visibility();
		Debug.calculate_visibility();
		Trace.calculate_visibility();
	} else {
		utility_exit_with_message("Tried to calculate the visibility of Tracer "+channel_+" before the init() was called!.");
	}
}


/// @details Calculate visibility (static version) of current Tracer object using channel name and priority.
/// result stored in 'muted' and 'visible'.
void TracerImpl::calculate_visibility(
	std::string const & channel,
	int    priority,
	bool & visible,
	bool & muted,
	int  & mute_level,
	bool   muted_by_default
)
{
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > lock( tracer_static_data_mutex() );
#endif

	assert( tracer_options_ );
	visible = false;
	if ( in(tracer_options_->muted, "all", true) ) {
		if ( in(tracer_options_->unmuted, channel, false) ) visible = true;
		else visible = false;
	} else {
		if ( in(tracer_options_->unmuted, "all", true) ) {
			if ( in(tracer_options_->muted, channel, false) ) visible = false;
			else visible = true;
		} else {  /// default bechavior: unmute unless muted_by_default is true
			if ( muted_by_default ) {
				if ( in(tracer_options_->unmuted, channel, false) ) visible = true;
				else visible = false;
			} else {
				if ( in(tracer_options_->muted, channel, false) ) visible = false;
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

	if ( in(tracer_options_->muted, "all_high_mpi_rank", true ) ) {
		if ( mpi_rank_>=2 ) visible = false; //* visible: master and 1st client: rank 0 and rank1
	}

	if ( in(tracer_options_->muted, "all_high_mpi_rank_filebuf", true ) ) {
		if ( mpi_rank_>=4 ) visible = false; //* visible: master, filebuf and 1st client: rank 0, 1, 2
	}

#endif
	muted = !visible;

	mute_level = tracer_options_->level;
	calculate_tracer_level(tracer_options_->levels, channel, false, mute_level);
	//std::cout << "levels:" << tracer_options_->levels <<" ch:" << channel << " mute_level:" << mute_level << " priority:" << priority << std::endl;

	if ( priority > mute_level ) visible = false;
}


/// @details Check if string representing channel 'ch' is in vector<string> v. Return true if channel
/// is in vector, false otherwise.
/// Two mode of operation:
/// Strict:  strict==true - channels compared verbatim.
/// Regular: strict==false - comparing with hierarchy in mind,
///                          ie: ch='basic.pose' v[0]='basic' --> will yield true.
bool TracerImpl::in(utility::vector1<std::string> const & v, std::string const & ch, bool strict)
{
	for ( size_t i=1; i<=v.size(); i++ ) {
		if ( v[i] == ch ) return true;

		if ( !strict ) {
			if ( ch.size() > v[i].size() ) {
				std::string s(ch);  s.resize(v[i].size());
				if ( s == v[i] ) return true;
			}
		}
	}
	return false;
}

/// @brief Output a message in a manner that is safe if the Tracers/output are poorly initialized.
void TracerImpl::safe_output(std::string const & message ) {
	if ( ios_hook() ) {
		*ios_hook() << message << std::endl;
	} else if ( final_stream() ) {
		*final_stream() << message << std::endl;
	} else {
		std::cerr << "Tracer Error: " << message << std::endl;
	}
}

/// Same as before but return integer value for matched channel or closest match (we asume that 'v' in levels format, ie like: <channel name>:level )
bool TracerImpl::calculate_tracer_level(utility::vector1<std::string> const & v, std::string const & ch, bool strict, int &res)
{
	unsigned int len = 0;
	bool math = false;
	//std::cout << "Entring:calculate_tracer_level: v=" << v << " ch:" << ch << std::endl;
	for ( size_t i=1; i<=v.size(); i++ ) {
		bool flag = false;
		utility::vector1< std::string > spl = utility::string_split(v[i], ':');
		if ( spl.size() != 2 ) {
			safe_output("WARNING: Cannot parse -out:levels setting '"+v[i]+"'. Does not follow the format of 'tracer:level'. Ignoring.");
			continue;
		}
		//std::cout << "Split:" << spl << " size:" << spl[1].size() << std::endl;

		if ( spl[1] == "all" && len == 0 ) flag = true;  // we can asume that 'all' is shorter then any valid core/protocol path... but we don't!

		if ( spl[1] == ch ) flag=true;

		if ( !strict ) {
			if ( ch.size() > spl[1].size() ) {
				std::string s(ch);  s.resize(spl[1].size());
				if ( s == spl[1] ) flag=true;
			}
		}
		if ( flag  && ( len <= spl[1].size() ) ) { // If specified twice, use the later one.
			math = true;
			len = spl[1].size();
			res = utility::string2int(spl[2]);
			std::string const spl2_lower = ObjexxFCL::lowercased( spl[2] );
			if ( spl2_lower == "fatal" )   res = t_fatal;
			if ( spl2_lower == "error"   || spl2_lower == "errors" )   res = t_error;
			if ( spl2_lower == "warning" || spl2_lower == "warnings" ) res = t_warning;
			if ( spl2_lower == "info" )    res = t_info;
			if ( spl2_lower == "debug" )   res = t_debug;
			if ( spl2_lower == "trace" )   res = t_trace;
			if ( res == -1 ) {
				safe_output( "WARNING: The setting '" + spl[2] + "' is not recognized as a valid tracer level." );
				res = t_info; // Set such that you get standard amount of output instead of none.
			}
			//std::cout << "Match:" << spl << " ch:" << ch << " res:"<< res << std::endl;
		} else {
			//std::cout << "Fail:" << spl << " ch:" << ch << " res:"<< res << std::endl;
		}

	}
	//std::cout << "Leaving:calculate_tracer_level: v=" << v << " ch:" << ch << " match:" << math <<" res:"<< res << std::endl;
	return math;
}


/// @details Write the contents of str to sout prepending the channel
/// name on each line if the print_channel_name flag is set.
template <class out_stream>
void TracerImpl::prepend_channel_name( out_stream & sout, std::string const &str )
{
	// This is only called from TracerImpl::t_flush() -- That function holds the tracer_static_data_mutex

	if ( str.empty() ) { return; } // Don't bother printing nothing (an empty line will have at least a carriage return.)

	// If we're calling this function, we're at the begining of the line.

	sout << this->Reset;

	sout << output_prefix_;

	if ( tracer_options_ && tracer_options_->print_channel_name ) {
		sout << channel_name_color_;
		sout << channel_ << ": ";
	}

#ifdef USEMPI
	sout << "(" << mpi_rank_ << ") ";
#endif

	if ( tracer_options_ && tracer_options_->timestamp ) {
		sout << utility::timestamp() << " ";
	}

	sout << this->Reset;

	// If the priority levels warrant it, add additional labeling.
	// Note that the stylings here are somewhat arbitrary.
	// (More important priorities get lower numbers.)
	if ( priority_ <=  t_warning ) { // Quick short-circuit for most output.
		if ( priority_ <= t_fatal ) {
			sout << this->Red << this->Bold << "[ FATAL ]" << this->Reset << ' ';
		} else if ( priority_ <= t_error ) {
			sout << this->Red << this->Bold << "[ ERROR ]" << this->Reset << ' ';
		} else if ( priority_ <= t_warning ) {
			sout << this->Bold << "[ WARNING ]" << this->Reset << ' ';
		}
	}

	sout << channel_color_;

	sout << str;
}


/// @details Inform Tracer that is contents was modified, and IO is in order.
void TracerImpl::t_flush(std::string const &str)
{
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > lock( tracer_static_data_mutex() );
#endif
	bool use_ios_hook = ios_hook() && ios_hook().get()!=this;
	use_ios_hook = use_ios_hook && ( in(monitoring_list_(), channel_, false) || in(monitoring_list_(), get_all_channels_string(), true ) );
	use_ios_hook = use_ios_hook && ( ios_hook_raw_() || visible() );

	if ( use_ios_hook ) {
		prepend_channel_name<otstream>( *ios_hook(), str );
		ios_hook()->flush();
	}

	if ( !super_mute_() && visible() ) {
		prepend_channel_name<std::ostream>( *final_stream(), str );
	}
}

/// @details Set OStringStream object to which all Tracers output
/// listed in the monitoring_channels_list should be copied.  Note
/// this copies the output of channels even if they are invisible or
/// muted.
///
/// When raw==false same as above above but gives the option get only the
/// visible and unmuted tracers.  It can be useful to get the raw
/// output for applications like the comparing tracers, where the
/// output should not change with command line parameters.  It can be
/// useful to get the non-raw output in applications like using the
/// jd2 with MPI, where the output each job should match the output if
/// it was run on a single processor.
void TracerImpl::set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list, bool raw)
{
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > lock( tracer_static_data_mutex() );
#endif
	if ( utility::pointer::dynamic_pointer_cast< TracerImpl >(tr) != nullptr || utility::pointer::dynamic_pointer_cast< TracerImpl::TracerProxyImpl >(tr) != nullptr ) {
		utility_exit_with_message("Error: Setting the ios_hook() to be a true Tracer or a TracerProxy is only going to end in grief! (Use a PyTracer instead.)");
	}
	ios_hook() = tr;
	monitoring_list_() = utility::split(monitoring_channels_list);
	ios_hook_raw_() = raw;
}

} // namespace basic
