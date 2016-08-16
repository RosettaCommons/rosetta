// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/basic/Tracer.hh
/// @brief  Tracer IO system
/// @author Sergey Lyskov


#ifndef INCLUDED_basic_Tracer_hh
#define INCLUDED_basic_Tracer_hh

#include <cassert>                                // for assert
#include <memory>                                 // for allocator, shared_ptr
#include <sstream>                                // for string, basic_strin...
#include <vector>                                 // for vector
#include <utility/pointer/ReferenceCount.hh>      // for ReferenceCount
#include <utility/vector1.hh>                     // for vector1

#include <utility/thread/backwards_thread_local.hh> // for THREAD_LOCAL

#ifdef WIN32
#include <utility/CSI_Sequence.hh>
#else
#include <utility/CSI_Sequence.fwd.hh>
#endif

namespace basic {

/// @brief
/// Priority levels for T() and Tracer object, modeled on the log4j project and its offspring.
/// Priorities in Tracer are still ints so users can pass other arbitrary integer values (for now).
enum TracerPriority {
	t_fatal   = 0,   //< The FATAL level designates very severe error events that will presumably lead the application to abort.
	t_error   = 100, //< The ERROR level designates error events that might still allow the application to continue running.
	t_warning = 200, //< The WARN level designates potentially harmful situations.
	t_info    = 300, //< The INFO level designates informational messages that highlight the progress of the application at coarse-grained level.
	t_debug   = 400, //< The DEBUG level designates fine-grained informational events that are most useful to debug an application.
	t_trace   = 500  //< The TRACE level designates finer-grained informational events than the DEBUG level.
};


/// @brief Base class for Tracer, TracerProxy and UTracer objects.
template <class CharT, class Traits = std::char_traits<CharT> >
class basic_otstream : public std::basic_ostream<CharT, Traits>, public utility::pointer::ReferenceCount
{
protected: /// Inner class declaration

	/// @brief Wrapper class for std::stringbuf
	template <class _CharT, class _Traits = std::char_traits<_CharT> >
	class basic_tstringbuf : public std::basic_stringbuf<_CharT, _Traits> {
	public:
		basic_tstringbuf(basic_otstream *ot) : otsream_(ot) {}
		virtual ~basic_tstringbuf() {}

	protected:
		virtual int sync() {
			otsream_->t_flush( this->str() ); //std::basic_stringbuf<CharT, Traits>::str() );
			//std::basic_stringbuf<CharT, Traits>::str("");
			this->str("");
			return 0;
		}
	private:
		basic_otstream *otsream_;
	};


public:
	basic_otstream() : std::basic_ostream<CharT, Traits> ( new basic_tstringbuf<CharT, Traits> (this) ) {}
	virtual ~basic_otstream() { delete this->rdbuf(); }


	/// @brief Return true if inner string buffer is empty.
	bool is_flushed() const {
		basic_tstringbuf<char> * buf = dynamic_cast< basic_tstringbuf<char> * >( this->rdbuf() );
		return buf->str().size() == 0;
	}

protected:

	/// @brief notification that flush function was called and inner buffer should be outputed.
	/// This is the mechanims by which the std::basic_stringbuf base class communicates with the
	/// Tracer and TracerProxy objects.
	virtual void t_flush(std::string const &) { assert("basic_otstream::t_flush"); };

private:
	basic_otstream(basic_otstream const & );


	/// Data members
	/// @brief inner string buffer
	//std::basic_stringbuf<CharT, Traits> * tstringbuf_;
};


typedef basic_otstream<char> otstream;

typedef utility::pointer::shared_ptr< otstream > otstreamOP;


/// @brief data structure to store all system level options for Tracer system.
struct TracerOptions
{
	/// @brief system priority level
	int level;

	/// @brief should channel name be printed during the IO?
	bool print_channel_name;

	/// @brief should a timestamp be added to the channel name?
	bool timestamp;

	/// @brief list of muted channels
	utility::vector1<std::string> muted;

	/// @brief list of unmuted channels
	utility::vector1<std::string> unmuted;

	/// @brief list of muted channels
	utility::vector1<std::string> levels;
};


/// @brief Class for handling user debug/warnings/errors.
///  Use instance of this class instead of 'std::cout' for all your regular io.
///  Channel argument must be related to the location of the source file. For example if you create
///  Tracer object in src/basic/scoring/myfile.cc,
///    then channel must be something like 'src.basic.scoring.myfile'
class Tracer :  public otstream
{
	/// @brief init Tracer object with given parameters. This is a helper function to be called from various constructors
	void init(
		std::string const & channel,
		std::string const & channel_color,
		std::string const & channel_name_color,
		TracerPriority priority,
		bool muted_by_default
	);

public:

	/// @brief Create Tracer object with given channel and priority
	Tracer(
		std::string const & channel = "",
		TracerPriority priority = t_info,
		bool muted_by_default = false
	);

	/// @brief Create Tracer object with channel color, channel name color and given channel and priority
	/// @details
	///  Ex:
	///  static THREAD_LOCAL basic::Tracer     Blue("blue",       CSI_Blue);
	///
	Tracer(
		std::string const & channel,
		std::string const & channel_color,
		std::string const & channel_name_color = "",
		TracerPriority priority = t_info,
		bool muted_by_default = false
	);



	virtual ~Tracer();

	/// @brief re-init using data from another tracer object.
	void init( Tracer const & tr );

	/// @brief flush tracer buffer and flush buffers of all
	///        sub-channels ie: Fatal, Error, Warning, Info, Debug, Trace
	void flush_all_channels();

	typedef std::ostream * OstreamPointer;

	/// @brief set ios hook for final tracer stream (deafult is std::cout).
	static OstreamPointer &final_stream();
	static void set_new_final_stream( std::ostream *new_final_stream );
	static void set_default_final_stream();


	/// @brief set ios hook for all tracer io operation.
	/// @param monitoring_channels_list is space separated list of channels.
	//static void set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list);
	static void set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list, bool raw=true);

	static std::string const & get_all_channels_string();  // PyRosetta helper function

	/// @brief Is this tracer currently visible?.
	bool visible() const { return visible_; }

	/// @brief is this tracer visible, if it used the given priority value?
	bool visible( int priority ) const;

	/// @brief get/set tracer priority level.
	int priority() const { return priority_; }
	Tracer & operator () (int priority);
	void priority(int priority);

	std::string const & channel() const { return channel_; }

	///@brief Get the channel color.
	std::string const &channel_color() { return channel_color_; }

	///@brief Set the channel color.
	///
	///@details
	/// This can be done in a stream like this:
	///  TR << TR.bgWhite << TR.Black << "Example" << TR.Reset << std::endl;
	void channel_color(std::string const &color) { channel_color_ = color; }

	std::string const &channel_name_color() { return channel_name_color_; }
	void channel_name_color(std::string const &color) { channel_name_color_ = color; }

	/// @brief get/set tracer options - global options for Tracer IO.
	static TracerOptions & tracer_options() { return tracer_options_; }

	/// @brief global super mute flag that allow to mute all io no matter what.
	static bool super_mute() { return super_mute_(); }
	static void super_mute(bool f) { super_mute_() = f; }

	static void flush_all_tracers();

	/// @brief This function should be invoked after the options system has been
	/// initialized, so that the visibility for all tracers that have so far been
	/// constructed and have been waiting for the options system to be initialized
	/// can now have their visibility calculated.  After this function completes,
	/// all newly-constructed Tracers will calculate their visibility in their
	/// constructors.  Visibility is no longer be calculated on a just-in-time
	/// basis and stored in mutable data members.
	static void calculate_tracer_visibilities();

public: /// Inner Classes
	/// @brief Small inner class acting as a proxy to an object that hold it.
	class TracerProxy : public otstream // std::ostringstream //
	{
	public:
		TracerProxy( Tracer & tracer, int priority, std::string const & channel );

		virtual ~TracerProxy();

		/// @brief determine the visibility of the proxy
		void calculate_visibility();
		/// @brief Adding this function to get around unused class data member warnings; you should
		/// never have to worry about whether the visibility for a TracerProxy has been calculated.
		bool visibility_calculated() const { return visibility_calculated_; }
		bool visible() const { return visible_; }

	protected:

		virtual void t_flush( std::string const & );

	private:
		Tracer & tracer_;
		int priority_;

		/// @brief We need to copy channel name here so we can generate appropriate 'warning' message
		/// in destructor, where tracer_ object is no longer valid.
		std::string channel_;

		/// @brief is channel visible?
		bool visible_;

		/// @brief is channel visibility already calculated?
		bool visibility_calculated_;
	};

	/// @brief channels with predefined priority levels.
	TracerProxy Fatal, Error, Warning, Info, Debug, Trace;

	/// @details Static objects holding various ASCII CSI codes (see utility/CSI_Sequence.hh)
	static utility::CSI_Sequence Reset, Bold, Underline,
		Black,   Red,   Green,   Yellow,   Blue,   Magenta,   Cyan,   White,
		bgBlack, bgRed, bgGreen, bgYellow, bgBlue, bgMagenta, bgCyan, bgWhite;

protected:
	/// @brief overload member function.
	virtual void t_flush(std::string const &);

private: /// Functions
	/// @brief copy constructor.
	Tracer( Tracer const & tr );


	/// @brief return true if channel is inside vector, some logic apply.
	static
	bool
	in( utility::vector1<std::string> const &, std::string const & channel, bool strict );

	/// @brief calculate channel priority with hierarchy in mind.
	static
	bool
	calculate_tracer_level(
		utility::vector1<std::string> const & v,
		std::string const & ch,
		bool strict,
		int &res
	);

	/// @brief Tracers must register themselves with the static array of all tracers so
	/// that they can be flushed en masse if need be.  This function is thread safe.
	static
	void
	register_tracer( Tracer * tracer );

	template <class out_stream>
	void prepend_channel_name( out_stream & sout, std::string const &str );

	/// @brief calcualte visibility of the current object depending of the channel name and priority.
	void calculate_visibility();

	/// @brief Adding this function to get around unused class data member warnings; you should
	/// never have to worry about whether the visibility for a Tracer has been calculated.
	bool visibility_calculated() const { return visibility_calculated_; }

	static void calculate_visibility(
		std::string const & channel,
		int    priority,
		bool & visible,
		bool & muted,
		int  & mute_level_,
		bool   muted_by_default
	);

	/// @brief Output a message in a manner that is safe if the Tracers/output are poorly initialized.
	static void safe_output(std::string const &);

private: /// Data members

	/// @brief channel name
	std::string channel_;

	/// @brief default colors for tracer output and tracer channel-name string (ie color of string such as: 'core.pose:')
	std::string channel_color_, channel_name_color_;

	/// @brief channel output priority level
	int priority_;

	/// @brief channel muted priority level (above which level is channel muted), calculated using user suppied -level and -levels options
	int mute_level_;

	/// @brief is channel visible?
	bool visible_;

	/// @brief is channel muted ?
	bool muted_;

	/// @brief is channel muted by default?
	bool muted_by_default_;

	/// @brief is current printing position a begining of the line?
	bool begining_of_the_line_;

	/// @brief is channel visibility already calculated?
	bool visibility_calculated_;

	/// static data members
	/// @brief link to Tracer like object where all output for selecting channels should go.
	static otstreamOP & ios_hook();

	/// @brief should the ios_hook_ the raw output?
	static bool & ios_hook_raw_();

	/// @brief list of channels for which outout should be redirected.
	static utility::vector1< std::string > monitoring_list_;

	/// @brief global option collection for Tracer IO.
	static TracerOptions tracer_options_;

	static bool initial_tracers_visibility_calculated_;

	/// @brief global super mute flag that allow to mute all io no matter what.
	static bool & super_mute_();

	/// @which Mpi rank is this process
	static int mpi_rank_;

	/// @brief T is special function for assign tracer property on the static object.
	friend Tracer & T(std::string const &, TracerPriority);
};

/// @brief Simple singleton class to hold the all_tracers_ array, which
/// otherwise suffers from funky double-construction problems when declared
/// as a static data member of Tracer.
class TracerManager {
public:
	static TracerManager * get_instance();
	std::vector< Tracer * > & all_tracers();

private:
	TracerManager();

private:
	static TracerManager * instance_;
	std::vector< Tracer * > all_tracers_;
};

/// @brief T is special function for assign tracer property on the static object.
Tracer & T(std::string const & channel, TracerPriority priority=t_info);

/// @brief Predefined Error tracer.
inline Tracer & Error(TracerPriority priority=t_error) { return T("Error", priority); }

/// @brief Predefined Warning tracer.
inline Tracer & Warning(TracerPriority priority=t_warning) { return T("Warning", priority); }


/// Special PyRosetta friendly Tracer like buffer. Use it to capture Tracer output with set_ios_hook
class PyTracer :  public otstream
{
public:
	//PyTracer(void) {}
	//virtual ~PyTracer() {}

	std::string buf() { return buf_; }
	void buf(std::string b) { buf_ = b; }

	virtual void output_callback(std::string) {};

protected:
	/// @brief overload member function.
	virtual void t_flush(std::string const &);

private:
	std::string buf_;
};

} // namespace basic



#ifdef NDEBUG  // faster version of Tracer IO for Release version

#ifdef CXX11

#include <utility/stream_util.hh>

namespace basic {

template <class T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type * = nullptr>
Tracer & operator <<( Tracer & TR, T const & entry ) {
	std::ostream &t(TR);
    if( TR.visible() ) { t << entry; }
    return TR;
}

template <class T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type * = nullptr>
Tracer::TracerProxy & operator <<( Tracer::TracerProxy & TR, T const & entry ) {
	std::ostream &t(TR);
    if( TR.visible() ) { t << entry; }
    return TR;
}

} // namespace basic

#else // CXX11

namespace basic {

template <class T>
Tracer & operator <<( Tracer & TR, T const & entry ) {
	std::ostream &t(TR);
    if( TR.visible() ) { t << entry; }
    return TR;
}

template <class T>
Tracer::TracerProxy & operator <<( Tracer::TracerProxy & TR, T const & entry ) {
	std::ostream &t(TR);
    if( TR.visible() ) { t << entry; }
    return TR;
}

} // namespace basic

#endif // CXX11

#endif // NDEBUG


#endif // INCLUDED_basic_tracer_hh
