// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/DenovoProteinDesign/CreateStartingStructureMover.fwd.hh
/// @brief  CreateStartingStructureMover forward declarations header
/// @author

#ifndef INCLUDED_devel_denovo_protein_design_CreateStartingStructureMover_fwd_hh
#define INCLUDED_devel_denovo_protein_design_CreateStartingStructureMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace devel{
namespace denovo_protein_design{

//Forwards and OP typedefs
class CreateStartingStructureMover;
typedef utility::pointer::shared_ptr< CreateStartingStructureMover > CreateStartingStructureMoverOP;
typedef utility::pointer::shared_ptr< CreateStartingStructureMover const > CreateStartingStructureMoverCOP;

}//DenovoProteinDesign
}//devel

#endif //INCLUDED_devel_DenovoProteinDesign_CreateStartingStructureMover_FWD_HH
