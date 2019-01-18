// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

// Unit headers
#include <--path--/--class--.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>



static basic::Tracer TR( "--namespace_dot--.--class--" );

--namespace--
using namespace protocols::jd3;

//Constructor
--class--::--class--()
{}

//Destructor
--class--::~--class--()
{}

jd3::CompletedJobOutput
--class--::run() {

        jd3::CompletedJobOutput output;

        //Ex:
        /*
        mover_->apply( *pose_ );
        core::Real const score = sfxn_->score( *pose_ );

        jd3::JobSummaryOP summary( utility::pointer::make_shared< jd3::standard::EnergyJobSummary >( score ) );
        jd3::JobResultOP result( utility::pointer::make_shared< jd3::standard::PoseJobResult >( pose_ ) );
        output.job_results.push_back( std::make_pair( summary, result ) );

        output.status = jd3::jd3_job_status_success;
        */
        return output;
}

--end_namespace--