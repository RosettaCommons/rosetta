// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @author Rocco Moretti (rmorettiase@gmail.com)
///

#ifndef INCLUDED_protocols_rotamer_gen_RDKitRotamers_fwd_hh
#define INCLUDED_protocols_rotamer_gen_RDKitRotamers_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace rotamer_gen {


class RDKitRotamers;
typedef utility::pointer::shared_ptr< RDKitRotamers > RDKitRotamersOP;

} // namespace rotamer_gen
} // namespace protocols

#endif
