// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/outputter/SilentFileOutputter.cc
/// @brief An outputter can take a PipeMapSP, a PipeSP, or a PoseSP and write it to a pdb
/// @author Ken Jung

// Unit Headers
#include <protocols/outputter/SilentFileOutputter.hh>
#include <core/io/silent/SilentStructFactory.hh>

// tracer
#include <basic/Tracer.hh>

namespace protocols {
namespace outputter {

static basic::Tracer TR("protocols.outputter.SilentFileOutputter");

#ifdef USELUA
void lregister_SilentFileOutputter( lua_State * lstate ) {
	lregister_FormatStringOutputter( lstate );

	luabind::module(lstate, "protocols")
	[
		luabind::namespace_( "outputter")
		[
			luabind::class_<SilentFileOutputter, FormatStringOutputter>("SilentFileOutputter")
				.def("score_only", (void (SilentFileOutputter::*) (bool) ) &SilentFileOutputter::score_only)
				.def("tag_format_string", (void (SilentFileOutputter::*) (std::string) ) &SilentFileOutputter::tag_format_string)
		]
	];
}
#endif

SilentFileOutputter::SilentFileOutputter(){}
SilentFileOutputter::~SilentFileOutputter(){}
OutputterSP SilentFileOutputter::create() {
	return OutputterSP( new SilentFileOutputter () );
}

void SilentFileOutputter::write( Pose & p ) {
	std::string outfilename;
	std::string tag;
	parse_format_string( filenameparts_, format_string_, outfilename );
	parse_format_string( filenameparts_, tag_format_string_, tag);
	core::io::silent::SilentStructOP tmp = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
	tmp->fill_struct( p, tag );
	sfd_.write_silent_struct( *tmp, outfilename, score_only_ );
}

#ifdef USELUA
void SilentFileOutputter::parse_def( utility::lua::LuaObject const & def,
		utility::lua::LuaObject const & tasks ) {
	format_string_ = def["format_string"] ? def["format_string"].to<std::string>() : "%filebasename_%filemultiplier.silent";
	score_only_= def["score_only"] ? def["score_only"].to<bool>() : false;
	binary_ = def["binary"] ? def["binary"].to<bool>() : true;
	tag_format_string_ = def["tag_format_string"] ? def["tag_format_string"].to<std::string>() : "%filebasename_%filemultiplier";
}

void SilentFileOutputter::lregister( lua_State * lstate ) {
	lregister_SilentFileOutputter(lstate);
}
#endif

} // outputter
} // protocols
