// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/outputter/FormatStringOutputter.hh
/// @brief An abstract base class that implements format string parsing
/// @author Ken Jung

#ifndef INCLUDED_protocols_outputter_FormatStringOutputter_hh
#define INCLUDED_protocols_outputter_FormatStringOutputter_hh

// Unit Headers
#include <protocols/outputter/FormatStringOutputter.fwd.hh>
#include <protocols/outputter/Outputter.hh>

// Project headers
#include <core/io/serialization/PipeMap.fwd.hh>
#include <core/io/serialization/Pipe.fwd.hh>
#include <boost/unordered_map.hpp>

namespace protocols {
namespace outputter {

		using namespace core::io::serialization;
		using core::pose::PoseSP;

#ifdef USELUA
		void lregister_FormatStringOutputter( lua_State * lstate );
#endif

class FormatStringOutputter : public Outputter {

	public:
		FormatStringOutputter();
		virtual ~FormatStringOutputter();

		std::string format_string() { return format_string_; }
		void format_string( std::string s ) { format_string_ = s; }

		virtual void write( PipeMap & p );
		virtual void write( Pipe & p );
		virtual void write( Pose & p )=0;

		void parse_format_string( boost::unordered_map< std::string, std::string> & filenameparts, std::string const & format_string, std::string & filename );
#ifdef USELUA
		virtual void parse_def( utility::lua::LuaObject const & def,
						utility::lua::LuaObject const & tasks );
		virtual void lregister( lua_State * lstate );
#endif

		// factory functions but this won't actually be created
		OutputterSP create();
		static std::string name() {
			return "FormatStringOutputter";
		}

	protected:
			boost::unordered_map<std::string, std::string> filenameparts_;
			std::string format_string_;

}; // end 

} // outputter
} // protocols


#endif //INCLUDED_protocols_outputter_FormatStringOutputter_hh
