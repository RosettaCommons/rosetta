// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @detailed responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/CutBiasClaimer.hh>

// Package Headers

// Project Headers
#include <core/fragment/SecondaryStructure.hh>

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


//#include <basic/options/option.hh>

//// C++ headers
// AUTO-REMOVED #include <fstream>

// option key includes


static basic::Tracer tr("protocols.topo_broker",basic::t_info);
//static numeric::random::RandomGenerator RG(181134);

namespace protocols {
namespace topology_broker {

using namespace core;

CutBiasClaimer::CutBiasClaimer() {}

CutBiasClaimer::CutBiasClaimer( core::fragment::SecondaryStructure const& ss )
{
	cut_bias_.reserve( ss.total_residue() );
	ObjexxFCL::FArray1D_float const& lf = ss.loop_fraction();
	for ( Size i = 1; i <= ss.total_residue(); i ++ ) {
		cut_bias_.push_back( lf( i ) );
	}
}

CutBiasClaimer::CutBiasClaimer( utility::vector1< core::Real > const& set ) {
	cut_bias_ = set;
}

void
CutBiasClaimer::manipulate_cut_bias( utility::vector1< core::Real >& tot_cut_bias ) {
	for ( Size i = 1; i<=cut_bias_.size() && i<=tot_cut_bias.size(); i++ ) {
		tot_cut_bias[ i ] *= cut_bias_[ i ];
	}
}

} //topology_broker
} //protocols
