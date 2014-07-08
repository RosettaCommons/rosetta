// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   sandbox
/// @brief  apps/pilot/yab/buns_test.cc
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// project headers
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/MetricValue.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <sstream>
#include <string>


static std::string BUNS_CALC_NAME( "buns" );

class BUNS_Output : public protocols::moves::Mover {
	typedef basic::MetricValue< utility::vector1< core::Size > > MV_SizeVec;
	typedef core::pose::Pose Pose;
	typedef std::string String;

	virtual void apply( Pose & pose ) {
		using core::Size;
		using protocols::jd2::JobDistributor;

		MV_SizeVec sv;
		pose.metric( BUNS_CALC_NAME, "residue_bur_unsat_polars", sv );

		std::ostringstream ss;
		for ( Size i = 1, ie = sv.value().size(); i <= ie; ++i ) {
			ss << "BUNS RES   " << i << "   " <<  sv.value()[ i ] << '\n';
		}

		JobDistributor::get_instance()->current_job()->add_string( ss.str() );
	}
};


int main( int argc, char * argv [] ) {
	try {
	using core::pose::metrics::CalculatorFactory;
	using core::scoring::ScoreFunctionOP;
	using protocols::simple_moves::ScoreMover;
	using protocols::simple_moves::ScoreMoverOP;
	using protocols::jd2::JobDistributor;
	using protocols::moves::SequenceMover;
	using protocols::moves::SequenceMoverOP;
	using protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator;

	// initialize rosetta
	devel::init( argc, argv );

	CalculatorFactory::Instance().register_calculator(
		BUNS_CALC_NAME,
		new BuriedUnsatisfiedPolarsCalculator(
			"default",
			"default"
		)
	);

	ScoreFunctionOP sfx = core::scoring::get_score_function();
	ScoreMoverOP scoremover = new ScoreMover( sfx );

	BUNS_Output * buns_output = new BUNS_Output();

	SequenceMoverOP seqmover = new SequenceMover;
	seqmover->add_mover( scoremover );
	seqmover->add_mover( buns_output );

	// run job
	JobDistributor::get_instance()->go( seqmover );
	 } catch ( utility::excn::EXCN_Base const & e ) {
		 std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
