// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/phosphate/PhosphateMove.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_rna_phosphate_PhosphateMove_HH
#define INCLUDED_protocols_stepwise_modeler_rna_phosphate_PhosphateMove_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/modeler/rna/phosphate/PhosphateMove.fwd.hh>
#include <core/types.hh>
#include <string>
#include <ostream>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace phosphate {

enum PhosphateTerminus { NONE, FIVE_PRIME_PHOSPHATE, THREE_PRIME_PHOSPHATE };

std::string
to_string( PhosphateTerminus const & phosphate_terminus );

class PhosphateMove: public utility::pointer::ReferenceCount {

public:

	//constructor
	PhosphateMove( Size const rsd,
		PhosphateTerminus const terminus ):
		rsd_( rsd ),
		terminus_( terminus )
	{
	}

	//destructor
	~PhosphateMove();

	bool operator== ( PhosphateMove const & other) const;

public:

	core::Size rsd() const { return rsd_; }
	PhosphateTerminus terminus() const { return terminus_; }

private:

	core::Size rsd_;
	PhosphateTerminus terminus_;

};

std::ostream &
operator <<( std::ostream & os, PhosphateMove const & phosphate_move );

} //phosphate
} //rna
} //modeler
} //stepwise
} //protocols

#endif
