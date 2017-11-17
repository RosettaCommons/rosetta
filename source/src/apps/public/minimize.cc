// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/public/minimize.cc
/// @brief Calls the minimizer to energy-relax a pose.
/// @author No idea who wrote this originally.  Modified by Vikram K. Mulligan (vmullig@uw.edu) on 20 June 2016 as part of the 2016 Documentation XRW.

#include <protocols/jd2/JobDistributor.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <numeric/random/random.hh>
#include <string>
#include <ObjexxFCL/string.functions.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


using namespace ObjexxFCL;

using basic::Error;
using basic::Warning;

using namespace core;
using namespace protocols;
using namespace protocols::moves;
using namespace core::scoring;

using utility::vector1;

static basic::Tracer TR( "minimize" );

class Minimize : public moves::Mover {

public:
	Minimize();

	~Minimize() override;

	MoverOP clone() const override;
	MoverOP fresh_instance() const override;

	void apply( Pose & pose ) override;
	std::string get_name() const override;
	void test_move( Pose & pose ) override
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

	if ( option[ in::file::fullatom ]() ) {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *score_function_ );
	} else {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *score_function_ );
	}


}

Minimize::~Minimize() = default;

MoverOP Minimize::clone() const {
	return MoverOP( new Minimize( *this ) );
}
MoverOP Minimize::fresh_instance() const {
	return MoverOP( new Minimize );
}

void
Minimize::apply( Pose & pose ) {
	using namespace pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	if ( option[ in::file::fullatom ]() ) {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( pose );
	} else {
		core::scoring::constraints::add_constraints_from_cmdline_to_pose( pose );
	}


	std::string min_type =  option[ OptionKeys::run::min_type ]();
	core::Real min_tol =  option[ OptionKeys::run::min_tolerance ]();
	core::optimization::MinimizerOptions options( min_type, min_tol, true, false );
	core::kinematics::MoveMap final_mm;

	if ( !option[ OptionKeys::in::file::movemap ].user() ) {
		final_mm.set_chi( option[ OptionKeys::relax::chi_move ].user() ? option[ OptionKeys::relax::chi_move ]() : true  );
		final_mm.set_bb( option[ OptionKeys::relax::bb_move ].user() ? option[ OptionKeys::relax::bb_move ]() : true  );
		final_mm.set_jump( option[ OptionKeys::relax::jump_move ].user() ? option[ OptionKeys::relax::jump_move ]() : true );
	} else {
		TR << "Initializing movemap from file " << option[ OptionKeys::in::file::movemap ]() << std::endl;
		final_mm.init_from_file(option[ OptionKeys::in::file::movemap ]() );
	}

	/*core::Real start_score =*/ (*score_function_)(pose);
	core::Size repeats = 1;
	for ( core::Size i = 0; i < repeats; i++ ) {
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

		MoverOP protocol( new Minimize() );
		protocols::jd2::JobDistributor::get_instance()->go( protocol );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
