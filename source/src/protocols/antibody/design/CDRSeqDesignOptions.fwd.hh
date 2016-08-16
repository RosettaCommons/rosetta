// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_CDRSeqDesignOptions_fwd_hh
#define INCLUDED_protocols_antibody_design_CDRSeqDesignOptions_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace antibody {
namespace design {


// Forward
class CDRSeqDesignOptions;

typedef utility::pointer::shared_ptr< CDRSeqDesignOptions > CDRSeqDesignOptionsOP;
typedef utility::pointer::shared_ptr< CDRSeqDesignOptions const > CDRSeqDesignOptionsCOP;


class CDRSeqDesignOptionsParser;

typedef utility::pointer::shared_ptr< CDRSeqDesignOptionsParser > CDRSeqDesignOptionsParserOP;
typedef utility::pointer::shared_ptr< CDRSeqDesignOptionsParser const > CDRSeqDesignOptionsParserCOP;

typedef utility::vector1<CDRSeqDesignOptionsOP> AntibodyCDRSeqDesignOptions;

}
}
}

#endif //INCLUDED_ TestClass2.fwd.hh





