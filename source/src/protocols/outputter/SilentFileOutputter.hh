// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/outputter/SilentFileOutputter.hh
/// @brief An outputter can take a PipeMapSP, a PipeSP, or a PoseSP and write it to a silent file( or many! )
/// @author Ken Jung

#ifndef INCLUDED_protocols_outputter_SilentFileOutputter_hh
#define INCLUDED_protocols_outputter_SilentFileOutputter_hh

// Unit Headers
#include <protocols/outputter/SilentFileOutputter.fwd.hh>
#include <protocols/outputter/FormatStringOutputter.hh>

#include <core/io/silent/SilentFileData.hh>

namespace protocols {
namespace outputter {

#ifdef USELUA
		void lregister_SilentFileOutputter( lua_State * lstate );
#endif

using namespace core::io::serialization;
using core::pose::PoseSP;

class SilentFileOutputter : public FormatStringOutputter {

public:
	SilentFileOutputter();
	virtual ~SilentFileOutputter();

	virtual void write( Pose & p );

	bool score_only() { return score_only_; }
	void score_only(bool b) { score_only_ = b; }

	std::string tag_format_string() { return tag_format_string_; }
	void tag_format_string( std::string s ) { tag_format_string_ = s; }

#ifdef USELUA
		void parse_def( utility::lua::LuaObject const & def,
						utility::lua::LuaObject const & tasks );
		virtual void lregister( lua_State * lstate );
#endif

	// factory functions
	OutputterSP create();
	static std::string name() {
		return "SilentFileOutputter";
	}

protected:
	core::io::silent::SilentFileData sfd_;
	bool score_only_;
	bool binary_;
	std::string tag_format_string_;

}; // end

} // outputter
} // protocols


#endif //INCLUDED_protocols_outputter_SilentFileOutputter_hh
