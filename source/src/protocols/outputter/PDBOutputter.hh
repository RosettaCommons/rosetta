// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/outputter/PDBOutputter.hh
/// @brief An outputter can take a PipeMapSP, a PipeSP, or a PoseSP and write it to a pdb
/// @author Ken Jung

#ifndef INCLUDED_protocols_outputter_PDBOutputter_hh
#define INCLUDED_protocols_outputter_PDBOutputter_hh

// Unit Headers
#include <protocols/outputter/PDBOutputter.fwd.hh>
#include <protocols/outputter/FormatStringOutputter.hh>

namespace protocols {
namespace outputter {

#ifdef USELUA
void lregister_PDBOutputter( lua_State * lstate );
#endif

using namespace core::io::serialization;
using core::pose::PoseSP;

class PDBOutputter : public FormatStringOutputter {

	public:
		PDBOutputter();
		virtual ~PDBOutputter();

		virtual void write( Pose & p );

		// factory functions
		OutputterSP create();
		static std::string name() {
			return "PDBOutputter";
		}

#ifdef USELUA
		virtual void lregister( lua_State * lstate );
#endif

}; // end 

} // outputter
} // protocols


#endif //INCLUDED_protocols_outputter_PDBOutputter_hh
