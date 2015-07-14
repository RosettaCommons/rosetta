// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/symmetry/util.cc
///
/// @brief      Quick calculations for symmetrizing membrane proteins
/// @details    Calculate membrane position based on a fully symmetrized spanning
///             topology. Part of an experiment to stabilize the symmetric scoring
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <protocols/membrane/symmetry/util.hh>

// Package Headers
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh> 
#include <core/conformation/Conformation.hh>

#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <protocols/membrane/util.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh> 
#include <basic/Tracer.hh> 

static basic::Tracer TR( "protocols.membrane.symmetry.util" );

namespace protocols {
namespace membrane {
namespace symmetry {

using namespace core;
using namespace core::conformation::membrane;
using namespace core::pose;

/// @brief Symmetrize Spans
/// @details Create a spanning topology to reflect the full symmetric
/// complex instead of just the asymmetric unit
SpanningTopologyOP
symmetrize_spans( Pose & pose, SpanningTopology & topology ) {

    using namespace core::pose::symmetry;
    using namespace core::conformation::symmetry;
    
    // Check before we cast
    if ( !is_symmetric( pose ) ) {
        utility_exit_with_message( "Cannot create a symmetric spanning topology from an asymmetric pose!" );
    }
    SymmetricConformation & symm_conf ( dynamic_cast< SymmetricConformation & > ( pose.conformation()) );
    
    // Create a new symmetric spanning topology object
    SpanningTopologyOP symmetrized_topology = SpanningTopologyOP( new SpanningTopology() );
    
    // Iterate through the subunits and add spans to the spanning topology
    core::Size nsubunits( symm_conf.Symmetry_Info()->subunits() );
    core::Size nres_monomer( symm_conf.Symmetry_Info()->get_nres_subunit() );
    for ( core::Size i = 1; i <= nsubunits; ++i ) {
        for ( core::Size j = 1; j <= topology.nspans(); ++j ) {
            core::Size new_start( topology.span( j )->start() + nres_monomer*(i-1)  );
            core::Size new_end( topology.span( j )->end() + nres_monomer*(i-1) );
            Span new_span = Span( new_start, new_end );
            symmetrized_topology->add_span( new_span );
        }
    }
    
    // Sanity check on number of transmembrane spans
    TR << "Orig spans: " << topology.nspans() << " Symmetrized spans: " << symmetrized_topology->nspans() << std::endl;
    
    return symmetrized_topology;
}

} // symmetry
} // membrane
} // protocols
