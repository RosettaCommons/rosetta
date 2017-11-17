// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/CutBiasClaimer.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/SequenceNumberResolver.hh>


// Package Headers
#include <protocols/topology_broker/Exceptions.hh>

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

// option key includes


static basic::Tracer tr( "protocols.topo_broker", basic::t_info );

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

CutBiasClaimer::CutBiasClaimer( core::fragment::SecondaryStructure const& ss, std::string label )
{
	cut_bias_.reserve( ss.total_residue() );
	ObjexxFCL::FArray1D_float const& lf = ss.loop_fraction();
	for ( Size i = 1; i <= ss.total_residue(); i ++ ) {
		cut_bias_.push_back( lf( i ) );
	}
	set_label( label );
}

CutBiasClaimer::CutBiasClaimer( utility::vector1< core::Real > const& set ) {
	cut_bias_ = set;
}

void
CutBiasClaimer::manipulate_cut_bias( utility::vector1< core::Real >& tot_cut_bias ) {


	core::Size offset = broker().sequence_number_resolver().offset( label() );

	if ( tot_cut_bias.size() < offset + cut_bias_.size() ) {
		std::ostringstream msg;
		msg << " CutBiasClaimer with label '" << label() << "' tried to change cut_bias at position " << ( offset + cut_bias_.size()) <<
			" while sequence is only " << tot_cut_bias.size() << " residues long. " << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RangeError,  msg.str() );
	} else {
		for ( core::Size i = 1; i <= cut_bias_.size(); i++ ) {
			tot_cut_bias[ i + offset ] *= cut_bias_[ i ];
		}
		tr.Debug << "Set cut_bias values in range [" << (1+offset) << "," <<  ( cut_bias_.size() + offset ) << "]" <<std::endl;
	}

	/*for ( Size i = 1; i<=cut_bias_.size() && i<=tot_cut_bias.size(); i++ ) {
	tot_cut_bias[ i ] *= cut_bias_[ i ];
	}*/
}

} //topology_broker
} //protocols
