// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/FARNA_Optimizer.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/FARNA_Optimizer.hh>
#include <protocols/farna/RNA_FragmentMonteCarlo.hh>
#include <protocols/farna/options/RNA_FragmentMonteCarloOptions.hh>
#include <protocols/farna/util.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.FARNA_Optimizer" );

using namespace core;

namespace protocols {
namespace farna {

//Constructor
FARNA_Optimizer::FARNA_Optimizer( utility::vector1< pose::PoseOP > const & pose_list,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::Size cycles /* = 0 */):
	pose_list_( pose_list ),
	scorefxn_( scorefxn ),
	cycles_( cycles )
{}

//Destructor
FARNA_Optimizer::~FARNA_Optimizer()
{}

void
FARNA_Optimizer::apply( core::pose::Pose & pose )
{
	for ( Size n = 1; n <= pose_list_.size(); n++ ) {
		pose = *pose_list_[ n ];

		RNA_FragmentMonteCarloOptionsOP options( new RNA_FragmentMonteCarloOptions );
		//   options->initialize_from_command_line();
		if ( cycles_ > 0 ) options->set_monte_carlo_cycles( cycles_ );
		options->set_bps_moves( true );
		options->set_filter_lores_base_pairs( false );
		options->set_filter_chain_closure( false );
		options->set_verbose( false );

		RNA_FragmentMonteCarloOP rna_fragment_monte_carlo( new RNA_FragmentMonteCarlo( options ) );
		rna_fragment_monte_carlo->set_denovo_scorefxn( scorefxn_ );
		if ( options->minimize_structure() ) rna_fragment_monte_carlo->set_hires_scorefxn( get_rna_hires_scorefxn() );
		rna_fragment_monte_carlo->set_refine_pose( true ); // no heating, etc.
		rna_fragment_monte_carlo->apply( pose );
	}
}


} //farna
} //protocols
