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


#ifndef INCLUDED_basic_Tracer_hh
#define INCLUDED_basic_Tracer_hh

#include <basic/Tracer.fwd.hh>
#include <basic/TracerImpl.hh> // Non-fwd is intentional, for inlining efficiency.

#include <unordered_map>

#include <utility/CSI_Sequence.hh>

namespace basic {

/// @brief Predefined Error-level tracer, for use in headers
Tracer & Error();

/// @brief Predefined Warning tracer, for use in headers.
Tracer & Warning();


/// @brief Class for handling user debug/warnings/errors.
///  Use instance of this class instead of 'std::cout' for all your regular io.
///  Channel argument must be related to the location of the source file. For example if you create
///  Tracer object in src/basic/scoring/myfile.cc,
///  then channel must be something like 'basic.scoring.myfile'
///
/// Intended usage is as a global-scope static object.
/// If you need a heap/stack allocated object, use the TracerImpl class directly.
//
/// @details The funky indirection to TracerImpl here is so that the heavy-weight TracerImpl object
/// is only constructed on first use, and isn't static or thread_local
class Tracer
{
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
	///  static basic::Tracer     Blue("blue",       utility::CSI::Blue);
	///
	Tracer(
		std::string const & channel,
		utility::CSI::CSI_Enum const & channel_color,
		utility::CSI::CSI_Enum const & channel_name_color = utility::CSI::Nothing,
		TracerPriority priority = t_info,
		bool muted_by_default = false
	);

	virtual ~Tracer();

	Tracer( Tracer const & ) = delete;
	Tracer & operator=( Tracer const & ) = delete;

	/// @brief output operator
	/// We return a TracerImpl instead of ourself to make subsequent access faster.
	template< typename T >
	TracerImpl & operator <<( T const & entry ) {
		TracerImpl & TR( tracer_impl() );
		TR << entry;
		return TR;
	}

	// Due to funky definition, have to special case std::endl and other stream manipulators
	TracerImpl & operator <<( std::ostream& (*entry)(std::ostream&) ) {
		TracerImpl & TR( tracer_impl() );
		TR << entry;
		return TR;
	}

	/// @brief Allow Tracer object to be passed to something expecting a std::ostream
	operator std::ostream&() { return dynamic_cast<std::ostream&>( tracer_impl() ); }

	///// @brief Allow Tracer object to be passed to something expecting a TracerImpl object
	//operator TracerImpl&() { return tracer_impl(); }

	void flush() {
		tracer_impl().flush();
	}

	/// @brief flush tracer buffer and flush buffers of all
	///        sub-channels ie: Fatal, Error, Warning, Info, Debug, Trace
	void flush_all_channels() {
		tracer_impl().flush_all_channels();
	}

	// want this to be inline-able
	bool visible() {
		return tracer_impl().visible();
	}

	bool visible( int priority ) {
		return tracer_impl().visible(priority);
	}

	std::string const & channel() {
		return tracer_impl().channel();
	}

	TracerImpl & operator () (int priority) {
		return (tracer_impl())(priority); // Call operator() on the TracerImpl class
	}

	//// The IOS manipulators

	std::streamsize width() { return tracer_impl().width(); }
	std::streamsize width (std::streamsize wide) { return tracer_impl().width(wide); }

	std::streamsize precision() { return tracer_impl().precision(); }
	std::streamsize precision (std::streamsize prec) { return tracer_impl().precision(prec); }

	std::ios_base::fmtflags flags() { return tracer_impl().flags(); }
	std::ios_base::fmtflags flags (std::ios_base::fmtflags fmtfl) { return tracer_impl().flags(fmtfl); }

	//// The priority proxy objects
	class TracerProxy
	{
	public:
		TracerProxy( Tracer & tracer, TracerPriority priority ):
			tracer_( tracer ),
			priority_( priority )
		{}

		TracerProxy( TracerProxy const & ) = delete;
		TracerProxy & operator=( TracerProxy const & ) = delete;

		/// @brief output operator
		/// We return a TracerProxyImpl instead of ourself to make subsequent access faster.
		template< typename T >
		TracerImpl::TracerProxyImpl & operator <<( T const & entry ) {
			TracerImpl::TracerProxyImpl & TR( tracer_proxy_impl() );
			TR << entry;
			return TR;
		}

		// Due to funky definition, have to special case std::endl and other stream manipulators
		// Have to special case std::endl, due to it's funky definition.
		TracerImpl::TracerProxyImpl & operator <<( std::ostream& (*entry)(std::ostream&) ) {
			TracerImpl::TracerProxyImpl & TR( tracer_proxy_impl() );
			TR << entry;
			return TR;
		}

		/// @brief Allow Tracer object to be passed to something expecting a std::ostream
		operator std::ostream&() { return dynamic_cast<std::ostream&>( tracer_proxy_impl() ); }

		///// @brief Allow Tracer object to be passed to something expecting a TracerImpl object
		//operator TracerImpl::TracerProxyImpl&() { return tracer_proxy_impl(); }

		void flush() {
			tracer_proxy_impl().flush();
		}

		bool visible() {
			return tracer_proxy_impl().visible();
		}

	private:

		TracerImpl::TracerProxyImpl & tracer_proxy_impl() {
			return tracer_.tracer_impl().get_proxy_by_priority( priority_ );
		}

	private:
		Tracer & tracer_;
		TracerPriority priority_;
	};

public:

	/// @brief channels with predefined priority levels.
	TracerProxy Fatal, Error, Warning, Info, Debug, Trace;

	/// @details These are just convenience references to the enum entries
	static utility::CSI::CSI_Enum Reset, Bold, Underline,
		Black,   Red,   Green,   Yellow,   Blue,   Magenta,   Cyan,   White,
		bgBlack, bgRed, bgGreen, bgYellow, bgBlue, bgMagenta, bgCyan, bgWhite;

	/// TODO: See if we can kill any usage of this, except for the esoteric
	static void set_ios_hook(otstreamOP tr, std::string const & monitoring_channels_list, bool raw=true) {
		TracerImpl::set_ios_hook(tr, monitoring_channels_list, raw);
	}

	static std::string const & get_all_channels_string() {
		return TracerImpl::get_all_channels_string();
	}

	/// @brief global super mute flag that allow to mute all io no matter what.
	static bool super_mute() { return TracerImpl::super_mute(); }
	static void super_mute(bool f) { TracerImpl::super_mute(f); }


protected:

	/// @brief The function which handles the construct-on-first use
	// Not virtual for speed in the usual (non-creation) case.
	TracerImpl & tracer_impl();

	virtual // virtual for MemTracer
	std::unique_ptr< TracerImpl > create_impl();

private:

	// Data needed to initialize the TracerImpl
	// Keep a copy around for construct-on-first-use.
	std::string channel_;
	utility::CSI::CSI_Enum channel_color_;
	utility::CSI::CSI_Enum channel_name_color_;
	TracerPriority priority_;
	bool muted_by_default_;

};


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



#endif // INCLUDED_basic_tracer_hh
