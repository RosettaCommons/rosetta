// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Oliver Lange, Mike Tyka
/// @brief


// libRosetta headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/option.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/electron_density/util.hh>

// AUTO-REMOVED #include <protocols/electron_density/util.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/util.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

// C++ headers
//#include <cstdlib>
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/broker.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/rescore.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <basic/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>



static basic::Tracer TR("main");

using namespace protocols::moves;
using namespace core::scoring;
// using namespace basic::options;
// using namespace basic::options::OptionKeys;



class MyScoreMover : public Mover {
public:
	MyScoreMover();

	virtual void apply( core::pose::Pose& pose );
	std::string get_name() const { return "MyScoreMover"; }

	virtual MoverOP clone() const {
		return new MyScoreMover( *this );
	}

	virtual	MoverOP	fresh_instance() const {
		return new MyScoreMover;
	}

	void set_keep_input_scores(){ keep_scores_flag_ = true; }
	void set_skip_scoring(){ skip_scoring_ = true; }
private:
	core::scoring::ScoreFunctionOP sfxn_;

	bool keep_scores_flag_;   // retains the previous scores from the input silent file or whatever
	bool skip_scoring_;       // skips the actual scoring call, calling evaluators only
	bool assign_ss_;          // Call the dssp object to assign a secondary structure.
};

MyScoreMover::MyScoreMover():
	keep_scores_flag_(false),
	skip_scoring_(false)
 {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	// get scorefxn and add constraints if defined
	sfxn_ = core::scoring::getScoreFunction();
	if ( option[ in::file::keep_input_scores ]() ){
		set_keep_input_scores();
	}
	if ( option[ rescore::skip ]() ){
		set_skip_scoring();
	}

	assign_ss_ = false;
	if ( option[ rescore::assign_ss ]() ){
		assign_ss_ = true;
	}

	// add cst scores from cmd line
	if( option[ in::file::fullatom ]() ) {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *sfxn_ );
	} else {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *sfxn_ );
	}

	// now add density scores from cmd line
	if ( option[ edensity::mapfile ].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *sfxn_ );
	}
}

void MyScoreMover::apply( core::pose::Pose& pose ) {
	if ( !keep_scores_flag_ ) {
		pose.energies().clear();
		pose.data().clear();
	}

	if ( assign_ss_ ) {
		core::scoring::dssp::Dssp my_dssp( pose );
		//core::scoring::dssp::DsspOP my_dssp_OP = new core::scoring::dssp::Dssp( pose );
		my_dssp.insert_ss_into_pose( pose );
	}

	sfxn_->set_weight( core::scoring::linear_chainbreak, 4.0/3.0 );
	sfxn_->set_weight( core::scoring::overlap_chainbreak, 1.0 );

	if ( ! skip_scoring_ ) {
		(*sfxn_)( pose );
	}
}

int
main( int argc, char * argv [] )
{
	try {

	using namespace protocols;
	using namespace protocols::jd2;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	jd2::register_options();

	OPT(rescore::assign_ss);
	OPT(rescore::skip);

	// initialize core
	devel::init(argc, argv);


	//The following lines are to ensure one can rescore the pcs energy term (that uses TopologyClaimer)
	if( option[ broker::setup ].user() ){
		protocols::topology_broker::TopologyBrokerOP top_bro_OP = new  topology_broker::TopologyBroker();
		try{
			add_cmdline_claims(*top_bro_OP, false /*do_I_need_fragments */);
		}
		catch ( utility::excn::EXCN_Exception &excn )  {
			excn.show( TR.Error );
			utility_exit();
		}
	}

	//MyScoreMover* scoremover = new MyScoreMover;
	MoverOP scoremover = new MyScoreMover;

	// add constraints from cmd line
	if ( option[ OptionKeys::constraints::cst_fa_file ].user() || option[ OptionKeys::constraints::cst_file ].user() ) {
			protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
			protocols::simple_moves::ConstraintSetMoverOP loadCsts( new protocols::simple_moves::ConstraintSetMover );
			if( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
				loadCsts->constraint_file( core::scoring::constraints::get_cst_fa_file_option() );
			} else {
				loadCsts->constraint_file( core::scoring::constraints::get_cst_file_option() );
			}
			seqmov->add_mover( loadCsts );
			seqmov->add_mover( scoremover );
			scoremover = seqmov;
	}

	// set pose for density scoring if a map was input
	//   + (potentially) dock map into density
	if ( option[ edensity::mapfile ].user() ) {
		protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
		seqmov->add_mover( new protocols::electron_density::SetupForDensityScoringMover );
		seqmov->add_mover( scoremover );
		scoremover = seqmov;
	}

	// set pose for symmetry
	if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
		protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
		seqmov->add_mover( new protocols::simple_moves::symmetry::SetupForSymmetryMover );
		seqmov->add_mover( scoremover );
		scoremover = seqmov;
	}

	using namespace protocols::jd2;

	// Make sure the default JobOutputter is SilentJobOutputter to ensure that when score_jd2
	// is called with default arguments is prints a proper scorefile and not the hacky thing that
	// the  JobOutputter scorefile() function produces (which for example skips Evaluators!!)

	// Set up a job outputter that writes a scorefile and no PDBs and no Silent Files.
	SilentFileJobOutputterOP jobout = new SilentFileJobOutputter;
	jobout->set_write_no_structures();
	jobout->set_write_separate_scorefile(true);

	// If the user chooses something else, then so be it, but by default score(_jd2) should only create a score
	// file and nothing else.
	protocols::jd2::JobDistributor::get_instance()->set_job_outputter( JobDistributorFactory::create_job_outputter( jobout ));

	try{
		JobDistributor::get_instance()->go( scoremover );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
		std::cout << "Exception: " << std::endl;
		excn.show( std::cout ); //so its also seen in a >LOG file
	}
	 } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

