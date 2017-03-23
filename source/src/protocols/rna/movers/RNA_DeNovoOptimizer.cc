// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/movers/RNA_DeNovoOptimizer.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/rna/movers/RNA_DeNovoOptimizer.hh>
#include <protocols/rna/denovo/RNA_FragmentMonteCarlo.hh>
#include <protocols/rna/denovo/options/RNA_FragmentMonteCarloOptions.hh>
#include <protocols/rna/denovo/util.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/Tracer.hh>
#include <utility>

static basic::Tracer TR( "protocols.rna.movers.RNA_DeNovoOptimizer" );

using namespace core;
using namespace protocols::rna::denovo;
using namespace protocols::rna::denovo::options;

namespace protocols {
namespace rna {
namespace movers {

//////////////////////////////////////////////////////////////////////////////////////////
//Constructor
RNA_DeNovoOptimizer::RNA_DeNovoOptimizer( utility::vector1< pose::PoseOP > const & pose_list,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::Size cycles /* = 0 */):
	pose_list_( pose_list ),
	scorefxn_(std::move( scorefxn )),
	cycles_( cycles )
{}

//////////////////////////////////////////////////////////////////////////////////////////
//Destructor
RNA_DeNovoOptimizer::~RNA_DeNovoOptimizer() = default;

//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_DeNovoOptimizer::apply( core::pose::Pose & pose )
{
	using namespace core::scoring;
	for ( Size n = 1; n <= pose_list_.size(); n++ ) {
		pose = *pose_list_[ n ];

		RNA_FragmentMonteCarloOptionsOP options( new RNA_FragmentMonteCarloOptions );
		options->initialize_for_farna_optimizer( cycles_ );
		rna_fragment_monte_carlo_ = RNA_FragmentMonteCarloOP( new RNA_FragmentMonteCarlo( options ) );

		if ( options->minimize_structure() ) {
			static ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_lores_for_stepwise.wts" );
			rna_fragment_monte_carlo_->set_denovo_scorefxn( lores_scorefxn );
			rna_fragment_monte_carlo_->set_hires_scorefxn( scorefxn_ );
		} else {
			rna_fragment_monte_carlo_->set_denovo_scorefxn( scorefxn_ );
		}

		rna_fragment_monte_carlo_->set_refine_pose( true ); // no heating, etc.
		rna_fragment_monte_carlo_->apply( pose );
	}
}



} //movers
} //rna
} //protocols
