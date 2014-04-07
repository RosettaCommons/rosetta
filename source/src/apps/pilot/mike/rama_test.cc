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


#include <protocols/jd2/JobDistributor.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/Energies.hh>
//#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DiagnosticData.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/ScoreMap.hh>
//#include <utility/exit.hh>
//#include <utility/file/FileName.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <string>
//#include <iostream>

// ObjexxFCL headers
//#include <ObjexxFCL/FArray1D.hh>


// libRosetta headers
// option key includes
// project headers
//silly using/typedef
// type headers
// unit headers
// utility headers


using basic::T;
using basic::Error;
using basic::Warning;


using namespace core;
using namespace protocols;
using namespace protocols::moves;
using namespace core::scoring;

using utility::vector1;

static basic::Tracer TR("rama_test");

class RamaTestMover : public moves::Mover {

public:
	RamaTestMover();

	~RamaTestMover();

	virtual MoverOP clone() const;
	virtual MoverOP fresh_instance() const;

	virtual void apply( Pose & pose );
	virtual std::string get_name() const;
	virtual void test_move( Pose & pose )
	{
		apply(pose);
	}

private:
	ScoreFunctionOP score_function_;
	std::map< std::string, core::Real > score_map_;
	bool verbose_;
	std::string scorefile_;
};

RamaTestMover::RamaTestMover() :
	Mover( "rama_test_mover" ),
	score_function_( getScoreFunction() )
{}

RamaTestMover::~RamaTestMover() {}

MoverOP RamaTestMover::clone() const {
	return new RamaTestMover( *this );
}
MoverOP RamaTestMover::fresh_instance() const {
	return new RamaTestMover;
}

void
RamaTestMover::apply( Pose & pose ) {
	using namespace pose;
	//using datacache::CacheableDataType::SCORE_MAP;

	(*score_function_)(pose);
	const core::Real limit=4.0;

	TR << protocols::jd2::JobDistributor::get_instance()->current_output_name();
	for( core::Size ir=1; ir <= pose.total_residue(); ++ ir ){
		EnergyMap emap = pose.energies().onebody_energies(ir);

		if( emap[ rama ] > limit )	TR << " " << ir << "(" << emap[ rama ] << ")";
	}
	TR <<  std::endl;

}

std::string
RamaTestMover::get_name() const {
	return "rama_test";
}

//ScoreFunctionOP RamaTestMover::score_function() const {
//	return score_function_;
//}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
    try {
    	using namespace protocols::jobdist;
    	using namespace protocols::moves;
    	using namespace scoring;

    	devel::init(argc, argv);

    	MoverOP protocol = new RamaTestMover();
    	protocols::jd2::JobDistributor::get_instance()->go( protocol );
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}
