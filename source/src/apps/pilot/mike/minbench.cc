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
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/Energies.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DiagnosticData.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/jd2/ScoreMap.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <string>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <ObjexxFCL/string.functions.hh>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


using namespace ObjexxFCL;

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace protocols;
using namespace protocols::moves;
using namespace core::scoring;

using utility::vector1;

static THREAD_LOCAL basic::Tracer TR( "sctrials" );

class Benchmark : public moves::Mover {

public:
	Benchmark();

	~Benchmark();

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
};

Benchmark::Benchmark() :
	Mover( "benchmark" ),
	score_function_( get_score_function() )
{}

Benchmark::~Benchmark() {}

MoverOP Benchmark::clone() const {
	return new Benchmark( *this );
}
MoverOP Benchmark::fresh_instance() const {
	return new Benchmark;
}

void
Benchmark::apply( Pose & pose ) {
	using namespace pose;
	//using datacache::CacheableDataType::SCORE_MAP;

	core::optimization::MinimizerOptions options( "dfpmin", 0.000, true, false );
	core::kinematics::MoveMap final_mm;
	final_mm.set_chi( true );
	final_mm.set_bb( true );

	int start_time = time(NULL);
	int end_time = time(NULL);
	int i=0;
	do {
		i++;
		Pose run_pose = pose;
		core::optimization::AtomTreeMinimizer().run( run_pose, final_mm, *score_function_, options );
		core::Real final_score = (*score_function_)(run_pose);
		std::cout << "FinalScore: " << final_score << std::endl;
		end_time = time(NULL);
	} while( (end_time-start_time) < 60 );
	end_time = time(NULL);

	std::cout << "TIME/min200: " << float(end_time - start_time)/float(i) << std::endl;

}

std::string
Benchmark::get_name() const {
	return "benchmark";
}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
    try {
    	using namespace protocols::jobdist;
    	using namespace protocols::moves;
    	using namespace scoring;

    	devel::init(argc, argv);

    	MoverOP protocol = new Benchmark();
    	protocols::jd2::JobDistributor::get_instance()->go( protocol );
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}
