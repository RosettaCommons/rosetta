// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Yifan Song
/// @author Frank DiMaio


#ifndef INCLUDED_protocols_hybridization_CartesianSampler_fwd_hh
#define INCLUDED_protocols_hybridization_CartesianSampler_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

class CartesianSampler;
typedef utility::pointer::shared_ptr< CartesianSampler > CartesianSamplerOP;
typedef utility::pointer::shared_ptr< CartesianSampler const > CartesianSamplerCOP;


} // hybridize
//} // comparative_modeling
} // protocols

#endif
