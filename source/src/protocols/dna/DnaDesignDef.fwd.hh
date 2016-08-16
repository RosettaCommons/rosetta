// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DnaDesignDef.fwd.hh
/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_dna_DnaDesignDef_fwd_hh
#define INCLUDED_protocols_dna_DnaDesignDef_fwd_hh

#include <utility/vector1.fwd.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace dna {

class DnaDesignDef;
typedef utility::pointer::shared_ptr< DnaDesignDef > DnaDesignDefOP;
typedef utility::pointer::shared_ptr< DnaDesignDef const > DnaDesignDefCOP;
typedef utility::vector1< DnaDesignDef > DnaDesignDefs;
typedef utility::vector1< DnaDesignDefOP > DnaDesignDefOPs;

} // namespace dna
} // namespace protocols

#endif
