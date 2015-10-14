// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/denovo_design/connection/ExtendChain.fwd.hh
/// @brief  ExtendChain forward header
/// @author Tom Linsky


#ifndef INCLUDED_protocols_denovo_design_connection_ExtendChain_fwd_hh
#define INCLUDED_protocols_denovo_design_connection_ExtendChain_fwd_hh


// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>


namespace protocols {
namespace denovo_design {
namespace connection {

// Forward
class ExtendChain;

// Types
typedef  utility::pointer::shared_ptr< ExtendChain >  ExtendChainOP;
typedef  utility::pointer::shared_ptr< ExtendChain const >  ExtendChainCOP;

typedef  utility::pointer::weak_ptr< ExtendChain >  ExtendChainAP;
typedef  utility::pointer::weak_ptr< ExtendChain const >  ExtendChainCAP;


} // namespace connection
} // namespace denovo_design
} // namespace protocols

#endif
