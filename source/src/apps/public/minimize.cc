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
// AUTO-REMOVED #include <core/pose/datacache/CacheableDataType.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/DunbrackRotamer.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/RotamerLibrary.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraint.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh>
// AUTO-REMOVED #include <basic/datacache/DiagnosticData.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <protocols/jd2/ScoreMap.hh>
// AUTO-REMOVED #include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <string>
// AUTO-REMOVED #include <utility/string_util.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

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

static basic::Tracer TR("minimize");
static numeric::random::RandomGenerator RG(91240391);

class Minimize : public moves::Mover {

public:
	Minimize();

	~Minimize();

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

Minimize::Minimize() :
	Mover( "benchmark" ),
	score_function_( get_score_function() )
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ in::file::fullatom ]() ){
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *score_function_ );
	}else{
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *score_function_ );
	}


}

Minimize::~Minimize() {}

MoverOP Minimize::clone() const {
	return new Minimize( *this );
}
MoverOP Minimize::fresh_instance() const {
	return new Minimize;
}

void
Minimize::apply( Pose & pose ) {
	using namespace pose;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;


	if ( option[ in::file::fullatom ]() ){
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose );
	}else{
		core::scoring::constraints::add_constraints_from_cmdline_to_pose( pose );
	}


	std::string min_type =  option[ OptionKeys::run::min_type ]();
	core::Real min_tol =  option[ OptionKeys::run::min_tolerance ]();
	core::optimization::MinimizerOptions options( min_type, min_tol, true, false );
	core::kinematics::MoveMap final_mm;
	final_mm.set_chi( true );
	final_mm.set_bb( true );

	/*core::Real start_score =*/ (*score_function_)(pose);
	core::Size repeats = 1;
	for(core::Size i = 0; i < repeats; i++ ){
		core::optimization::AtomTreeMinimizer().run( pose, final_mm, *score_function_, options );
		TR << "Score: " << i << "  " <<  (*score_function_)(pose) << std::endl;
	}

	core::Real final_score = (*score_function_)(pose);
	TR << "FinalScore: " << final_score << std::endl;

}

std::string
Minimize::get_name() const {
	return "Minimize";
}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
	using namespace protocols::moves;
	using namespace scoring;

	devel::init(argc, argv);

	MoverOP protocol = new Minimize();
	protocols::jd2::JobDistributor::get_instance()->go( protocol );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
