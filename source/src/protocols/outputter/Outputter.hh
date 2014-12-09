// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/outputter/Outputter.hh
/// @brief An outputter can take a PipeMapSP, a PipeSP, or a PoseSP and write it to a file
/// @author Ken Jung

#ifndef INCLUDED_protocols_outputter_Outputter_hh
#define INCLUDED_protocols_outputter_Outputter_hh

// Unit Headers
#include <protocols/outputter/Outputter.fwd.hh>

// Project headers
#include <core/io/serialization/PipeMap.fwd.hh>
#include <core/io/serialization/Pipe.fwd.hh>
#include <core/pose/Pose.hh>

#include <utility/lua/LuaObject.hh>

namespace protocols {
namespace outputter {

#ifdef USELUA
void lregister_Outputter( lua_State * lstate );
#endif

		using namespace core::io::serialization;
		using core::pose::Pose;

class Outputter {

	public:
		Outputter();
		virtual ~Outputter();

		virtual void write( PipeMap & p )=0;
		virtual void write( Pipe & p )=0;
		virtual void write( Pose & p )=0;

#ifdef USELUA
		virtual void parse_def( utility::lua::LuaObject const & def,
						utility::lua::LuaObject const & tasks ) = 0;
		virtual void lregister( lua_State * lstate )=0;
#endif

		// factory functions
		virtual OutputterSP create() = 0;
		static std::string name() {
			return "UNDEFINED NAME";
		}


}; // end Outputter base class

} // outputter
} // protocols


#endif //INCLUDED_protocols_outputter_Outputter_hh
