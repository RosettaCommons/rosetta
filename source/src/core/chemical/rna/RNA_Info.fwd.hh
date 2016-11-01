// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/chemical/rna/RNA_Info.fwd.hh
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)

#ifndef INCLUDED_core_chemical_rna_RNA_Info_fwd_hh
#define INCLUDED_core_chemical_rna_RNA_Info_fwd_hh

#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace chemical {
namespace rna {

class RNA_Info;

typedef  utility::pointer::weak_ptr< RNA_Info >  RNA_InfoAP;
typedef  utility::pointer::weak_ptr< RNA_Info const >  RNA_InfoCAP;
typedef  utility::pointer::shared_ptr< RNA_Info >  RNA_InfoOP;
typedef  utility::pointer::shared_ptr< RNA_Info const >  RNA_InfoCOP;

}
}
}

#endif
