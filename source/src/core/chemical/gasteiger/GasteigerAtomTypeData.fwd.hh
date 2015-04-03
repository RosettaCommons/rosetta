// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/gasteiger/GasteigerAtomTypeData.fwd.hh
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_gasteiger_GasteigerAtomTypeData_fwd_hh
#define INCLUDED_core_chemical_gasteiger_GasteigerAtomTypeData_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace chemical {
namespace gasteiger {

class GasteigerAtomTypeData;

typedef  utility::pointer::shared_ptr< GasteigerAtomTypeData >  GasteigerAtomTypeDataOP;
typedef  utility::pointer::shared_ptr< GasteigerAtomTypeData const >  GasteigerAtomTypeDataCOP;

} // gasteiger
} // chemical
} // core


#endif
