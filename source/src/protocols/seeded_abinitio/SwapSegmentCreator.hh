// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
///
/// @author Eva-Maria Strauch (evas01@u.washington.edu), May 2011


#ifndef INCLUDED_protocols_seeded_abinitio_SwapSegmentCreator_hh
#define INCLUDED_protocols_seeded_abinitio_SwapSegmentCreator_hh

// Project headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace seeded_abinitio {

class SwapSegmentCreator : public moves::MoverCreator
{
public:
	moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static  std::string mover_name();

};

}
}
#endif
