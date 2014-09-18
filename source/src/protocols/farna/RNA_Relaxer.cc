// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to RNA_Relaxer
/// @detailed
/// @author Rhiju Das


#include <protocols/farna/RNA_Relaxer.hh>
#include <protocols/farna/RNA_FragmentMover.hh>
#include <protocols/farna/RNA_Minimizer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/Pose.hh>

#include <core/scoring/constraints/ConstraintSet.hh>

#include <basic/options/option.hh>


//Relaxer stuff

// ObjexxFCL Headers

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>

// External library headers


//C++ headers
#include <vector>
#include <string>
#include <sstream>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/kinematics/Jump.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end



using namespace core;
using basic::T;


static thread_local basic::Tracer TR( "protocols.rna.RNA_Relaxer" );

namespace protocols {
namespace farna {

 /////////////////////////////////////////////////////
 // BASIC IDEA -- relax consists of fragment moves
 //  that aren't too perturbing, followed by minimization
 // So we need a fragment mover and a minimizer.

RNA_Relaxer::RNA_Relaxer( RNA_FragmentMoverOP & rna_fragment_mover, RNA_MinimizerOP & rna_minimizer ):
  Mover(),
	rna_fragment_mover_( rna_fragment_mover ),
	rna_minimizer_( rna_minimizer ),
	relax_cycles_( 50 ),
	max_frag_size_( 3 ),
	num_find_fragment_tries_( 100 ),
	rmsd_min_cutoff_( 1.0 ),
	rmsd_max_cutoff_( 3.5 ),
	simple_rmsd_cutoff_relax_( false )
{
	Mover::type("RNA_Relaxer");
}

RNA_Relaxer::~RNA_Relaxer()
{
}

/// @details  Apply the RNA full atom relaxer.
///
void RNA_Relaxer::apply( core::pose::Pose & pose	)
{

	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	//COMMENT THIS OUT!!!
	//	rna_minimizer_->use_coordinate_constraints( false );
	//rna_minimizer_->skip_o2prime_trials( true );

	time_t pdb_start_time = time(NULL);

	ScoreFunctionOP scorefxn;

	if ( option[ score::weights ].user() ) {
		scorefxn = get_score_function();
	} else {
		scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	}

	if (pose.constraint_set()->has_constraints() ){
		scorefxn->set_weight( atom_pair_constraint, 1.0 );
		scorefxn->set_weight( angle_constraint, 1.0 );
	}

	protocols::moves::MonteCarloOP monte_carlo_( new protocols::moves::MonteCarlo( pose, *scorefxn, 2.0 ) );

	for (Size r = 1; r <= relax_cycles_; r++ ) {

		TR << "RNA relaxer ROUND " << r << " of " << relax_cycles_ << std::endl;

		make_fragment_moves( pose );

		//Minimize it.
		rna_minimizer_->apply( pose );
		monte_carlo_->boltzmann( pose, "relax" );

	}

	monte_carlo_->show_counters();

	pose = monte_carlo_->lowest_score_pose();

	time_t pdb_end_time = time(NULL);

	scorefxn->show( std::cout, pose );

	TR << "RNA relaxer finished in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;

}

std::string
RNA_Relaxer::get_name() const {
	return "RNA_Relaxer";
}

///////////////////////////////////////////////////////////////////////////////////////////
void
RNA_Relaxer::make_fragment_moves( pose::Pose & pose )
{
	if ( simple_rmsd_cutoff_relax_ ) {
		find_fragment_by_simple_rmsd_cutoff( pose );
	} else { //default
		lores_monte_carlo( pose );
	}
}

///////////////////////////////////////////////////////////////////////////////////////////
void
RNA_Relaxer::find_fragment_by_simple_rmsd_cutoff( pose::Pose & pose )
{
	//Find a fragment to try... don't want to deviate too much from the input pose.
	pose::Pose start_pose = pose;

	Size const frag_size = static_cast <int> ( numeric::random::rg().uniform() * max_frag_size_ ) + 1;

	for (Size n = 1; n <= num_find_fragment_tries_ ; n++ ) {

		pose = start_pose;

		rna_fragment_mover_->random_fragment_insertion( pose, frag_size );

		Real const rmsd_to_start = scoring::all_atom_rmsd( pose, start_pose );
		TR << "Testing fragment insertion --> Length: " << frag_size << "  rmsd: " << rmsd_to_start << std::endl;
		if ( rmsd_to_start > rmsd_min_cutoff_ && rmsd_to_start < rmsd_max_cutoff_ ) break;

	}

}

///////////////////////////////////////////////////////////////////////////////////////////
void
RNA_Relaxer::lores_monte_carlo( pose::Pose & pose )
{

	using namespace core::scoring;

	pose::Pose start_pose = pose;

	static ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );
	if (pose.constraint_set()->has_constraints() ){
		lores_scorefxn->set_weight( atom_pair_constraint, 1.0 );
		lores_scorefxn->set_weight( angle_constraint, 1.0 );
	}

	Size const lores_monte_carlo_rounds( 1000 );
	protocols::moves::MonteCarloOP lores_monte_carlo_( new protocols::moves::MonteCarlo( pose, *lores_scorefxn, 2.0 ) );

	for( Size i=1; i <= lores_monte_carlo_rounds ; ++i ) {

		Size const frag_size = static_cast <int> ( numeric::random::rg().uniform() * max_frag_size_ ) + 1;
		rna_fragment_mover_->random_fragment_insertion( pose, frag_size );
		lores_monte_carlo_->boltzmann( pose, "frag" + SS(frag_size) );

	}

	lores_monte_carlo_->show_counters();

	// Want some variety? Don't use the lowest score pose...
	//	pose = lores_monte_carlo_->lowest_score_pose();

	Real const rmsd_to_start = scoring::all_atom_rmsd( pose, start_pose );
	TR << " LoRes MonteCarlo -->  rmsd: " << rmsd_to_start << std::endl;

}


} //farna
} //protocols
