// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/outputter/PDBOutputter.cc
/// @brief An outputter can take a PipeMapSP, a PipeSP, or a PoseSP and write it to a pdb
/// @author Ken Jung

// Unit Headers
#include <protocols/outputter/PDBOutputter.hh>

// tracer
#include <basic/Tracer.hh>

namespace protocols {
namespace outputter {

static thread_local basic::Tracer TR( "protocols.outputter.PDBOutputter" );

#ifdef USELUA
void lregister_PDBOutputter( lua_State * lstate ) {
	lregister_FormatStringOutputter( lstate );
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("outputter")
		[
			luabind::class_<PDBOutputter, FormatStringOutputter>("PDBOutputter")
		]
	];
}
#endif

PDBOutputter::PDBOutputter(){}
PDBOutputter::~PDBOutputter(){}
OutputterSP PDBOutputter::create() {
	return OutputterSP( new PDBOutputter () );
}

void PDBOutputter::write( Pose & p ) {
	std::string outfilename;
	parse_format_string( filenameparts_, format_string_, outfilename );
	p.dump_pdb( outfilename );
}

#ifdef USELUA
void PDBOutputter::lregister( lua_State * lstate ) {
	lregister_PDBOutputter( lstate );
}
#endif

} // outputter
} // protocols
