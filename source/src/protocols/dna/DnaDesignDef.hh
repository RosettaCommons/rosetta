// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DnaDesignDef.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_dna_DnaDesignDef_hh
#define INCLUDED_protocols_dna_DnaDesignDef_hh

#include <protocols/dna/DnaDesignDef.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.fwd.hh>

#include <iosfwd>
#include <string>

namespace protocols {
namespace dna {

/// @brief command-line dna_defs are of the format "C.501.ADE"
/// they are parsed here into this little class for convenience
class DnaDesignDef : public utility::pointer::ReferenceCount {
public:
	char chain;
	int pdbpos; // store pdb position (can be negative), convert to rosetta index later
	std::string name3; // store as string, convert to AA or Residue type later

	// string constructor
	DnaDesignDef( std::string const & );
	~DnaDesignDef() override;
};

std::ostream & operator << ( std::ostream &, DnaDesignDef const & );
std::ostream & operator << ( std::ostream &, DnaDesignDefs const & );
std::ostream & operator << ( std::ostream &, DnaDesignDefOPs const & );

} // namespace dna
} // namespace protocols

#endif
