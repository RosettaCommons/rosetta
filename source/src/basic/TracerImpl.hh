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
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_basic_TracerImpl_hh
#define INCLUDED_basic_TracerImpl_hh

#include <basic/Tracer.fwd.hh>

#include <cassert>                                // for assert
#include <sstream>                                // for basic_stringbuf
#include <utility/vector1.hh>                     // for vector1

#include <utility/CSI_Sequence.hh>


namespace basic {

/// @brief Base class for TracerImpl, TracerProxyImpl and UTracer objects.
template <class CharT, class Traits = std::char_traits<CharT> >
class basic_otstream : public std::basic_ostream<CharT, Traits>
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
	/// TracerImpl and TracerProxyImpl objects.
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
	int level = 300;

	/// @brief should channel name be printed during the IO?
	bool print_channel_name = true;

	/// @brief should a timestamp be added to the channel name?
	bool timestamp = false;

	/// @brief list of muted channels
	utility::vector1<std::string> muted;

	/// @brief list of unmuted channels
	utility::vector1<std::string> unmuted;

	/// @brief list of muted channels
	utility::vector1<std::string> levels;

	bool operator==( TracerOptions const & other ) const;
	bool operator!=( TracerOptions const & other ) const;
};

typedef utility::pointer::shared_ptr< TracerOptions > TracerOptionsOP;

/// @brief Class for handling user debug/warnings/errors.
///  Use instance of this class instead of 'std::cout' for all your regular io.
///  Channel argument must be related to the location of the source file. For example if you create
///  Tracer object in src/basic/scoring/myfile.cc,
///    then channel must be something like 'src.basic.scoring.myfile'
class TracerImpl :  public otstream
{
	/// @brief init Tracer object with given parameters. This is a helper function to be called from various constructors
	void init(
		std::string const & channel,
		utility::CSI_Sequence const & channel_color,
		utility::CSI_Sequence const & channel_name_color,
		TracerPriority priority,
		bool muted_by_default
	);

public:

	/// @brief Create Tracer object with given channel and priority
	TracerImpl(
		std::string const & channel = "",
		TracerPriority priority = t_info,
		bool muted_by_default = false
	);

	/// @brief Create Tracer object with channel color, channel name color and given channel and priority
	/// @details
	///  Ex:
	///  static basic::Tracer     Blue("blue",       CSI_Blue());
	///
	TracerImpl(
		std::string const & channel,
		utility::CSI_Sequence const & channel_color,
		utility::CSI_Sequence const & channel_name_color = utility::CSI_Nothing(),
		TracerPriority priority = t_info,
		bool muted_by_default = false
	);

	virtual ~TracerImpl();

	/// @brief re-init using data from another tracer object.
	void init( TracerImpl const & tr );

	/// @brief flush tracer buffer and flush buffers of all
	///        sub-channels ie: Fatal, Error, Warning, Info, Debug, Trace
	void flush_all_channels();

	typedef std::ostream * OstreamPointer;

	/// @brief Is this tracer currently visible?.
	//Inlined for speed.
	bool visible() {
		return visible_;
	}

	/// @brief is this tracer visible, if it used the given priority value?
	bool visible( int priority );

	/// @brief get/set tracer priority level.
	int priority() const { return priority_; }
	TracerImpl & operator () (int priority);
	void priority(int priority);

	std::string const & channel() const { return channel_; }

	///@brief Get the channel color.
	utility::CSI_Sequence const &channel_color() { return channel_color_; }

	///@brief Set the channel color.
	///
	///@details
	/// This can be done in a stream like this:
	///  TR << TR.bgWhite << TR.Black << "Example" << TR.Reset << std::endl;
	void channel_color(utility::CSI_Sequence const &color) { channel_color_ = color; }

	utility::CSI_Sequence const &channel_name_color() { return channel_name_color_; }
	void channel_name_color(utility::CSI_Sequence const &color) { channel_name_color_ = color; }

public: // Public static member functions

	/// @brief set ios hook for final tracer stream (deafult is std::cout).
	static void set_new_final_stream( std::ostream *new_final_stream );
	static void set_default_final_stream();

	/// @brief set ios hook for all tracer io operation.
	/// @param monitoring_channels_list is space separated list of channels.
	//static void set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list);
	static void set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list, bool raw=true);

	static std::string const & get_all_channels_string();  // PyRosetta helper function

	/// @brief set tracer options - global options for Tracer IO.
	static void set_tracer_options( TracerOptions const & to );

	/// @brief global super mute flag that allow to mute all io no matter what.
	static bool super_mute() { return super_mute_(); }
	static void super_mute(bool f) { super_mute_() = f; }

	/// @brief set/get globale string-prefix for all Tracer output strings
	static void output_prefix(std::string const &);
	static std::string output_prefix();

public: /// Inner Classes
	/// @brief Small inner class acting as a proxy to an object that hold it.
	class TracerProxyImpl : public otstream // std::ostringstream //
	{
	public:
		TracerProxyImpl( TracerImpl & tracer, int priority, std::string const & channel );

		virtual ~TracerProxyImpl();

		/// @brief Is the output from this TracerProxy visible?
		bool visible() { return visible_; }

		/// @brief determine the visibility of the proxy.
		void calculate_visibility();

	protected:

		virtual void t_flush( std::string const & );

	private:
		TracerImpl & tracer_;
		int priority_;

		/// @brief We need to copy channel name here so we can generate appropriate 'warning' message
		/// in destructor, where tracer_ object is no longer valid.
		std::string channel_;

		/// @brief is channel visible?
		bool visible_;
	};

	/// @brief channels with predefined priority levels.
	TracerProxyImpl Fatal, Error, Warning, Info, Debug, Trace;

	TracerProxyImpl & get_proxy_by_priority( TracerPriority priority );

	/// @details Static objects holding various ASCII CSI codes (see utility/CSI_Sequence.hh)
	static utility::CSI_Sequence Reset, Bold, Underline,
		Black,   Red,   Green,   Yellow,   Blue,   Magenta,   Cyan,   White,
		bgBlack, bgRed, bgGreen, bgYellow, bgBlue, bgMagenta, bgCyan, bgWhite;

protected:
	/// @brief overload member function.
	virtual void t_flush(std::string const &);

private: /// Functions
	/// @brief copy constructor.
	TracerImpl( TracerImpl const & tr ) = delete;

	static OstreamPointer &final_stream();

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

	/// @brief helper function for t_flush() -- do not call from elsewhere (without caution)
	template <class out_stream>
	void prepend_channel_name( out_stream & sout, std::string const &str );

	/// @brief calculate visibility of the current object depending of the channel name and priority.
	void calculate_visibility();

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
	utility::CSI_Sequence channel_color_, channel_name_color_;

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

	/// static data members
	/// @brief link to Tracer like object where all output for selecting channels should go.
	static otstreamOP & ios_hook(); // Not necessarily mutex protected - I/O mangling is

	/// @brief should the ios_hook_ the raw output?
	static bool & ios_hook_raw_(); // Use should be `tracer_static_data_mutex()` protected

	/// @brief list of channels for which outout should be redirected.
	static utility::vector1< std::string > & monitoring_list_(); // Use should be `tracer_static_data_mutex()` protected

	/// @brief global option collection for Tracer IO.
	static TracerOptionsOP tracer_options_; // Use should be `tracer_static_data_mutex()` protected

	/// @brief global super mute flag that allow to mute all io no matter what.
	static bool & super_mute_(); // Use should be `tracer_static_data_mutex()` protected

	/// @which Mpi rank is this process
	static int mpi_rank_; // Use should be `tracer_static_data_mutex()` protected

	/// @brief global prefix for all tracer output
	static std::string output_prefix_; // Use should be `tracer_static_data_mutex()` protected

};

} // namespace basic


#ifdef NDEBUG  // faster version of Tracer IO for Release version


#include <utility/stream_util.hh>

namespace basic {

template <class T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type * = nullptr>
TracerImpl & operator <<( TracerImpl & TR, T const & entry ) {
	std::ostream &t(TR);
	if( TR.visible() ) { t << entry; }
	return TR;
}

template <class T, typename std::enable_if< utility::has_insertion_operator_s<T>::value >::type * = nullptr>
TracerImpl::TracerProxyImpl & operator <<( TracerImpl::TracerProxyImpl & TR, T const & entry ) {
	std::ostream &t(TR);
	if( TR.visible() ) { t << entry; }
	return TR;
}

} // namespace basic


#endif // NDEBUG


#endif // INCLUDED_basic_tracer_hh
