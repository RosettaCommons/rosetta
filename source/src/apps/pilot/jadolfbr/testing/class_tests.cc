// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/jadolfbr/testing/class_tests.cc
/// @brief Just a file for testing various classes
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <devel/init.hh>

#include <protocols/toolbox/task_operations/RestrictToMoveMapChiOperation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>


	
class TestMMOP : public protocols::moves::Mover{

public:	
	
	TestMMOP(){};
	
	virtual ~TestMMOP(){};
	
	virtual
	std::string
	get_name() const {
		return "TestMMOP";
	}
	
	virtual
	void
	apply(core::pose::Pose & pose){
		
		using namespace core::pack::task::operation;
		using namespace protocols::toolbox::task_operations;
		using namespace protocols::simple_moves;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		
		core::pack::task::TaskFactoryOP tf = new core::pack::task::TaskFactory();
		tf->push_back(new InitializeFromCommandline);
		
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
		
		if (option[ OptionKeys::in::file::movemap ].user()) {
			movemap->init_from_file(option[ OptionKeys::in::file::movemap ]() );
		}
		else{
			movemap->set_jump( option[ OptionKeys::relax::jump_move ]() );
			movemap->set_bb( option[ OptionKeys::relax::bb_move ]() );
			movemap->set_chi( option[ OptionKeys::relax::chi_move ]() );	
		}
		RestrictToMoveMapChiOperationOP mmop = new RestrictToMoveMapChiOperation(movemap);
		mmop->set_design(true);
		mmop->set_include_neighbors(true);
		mmop->set_cutoff_distance(10.0);
		
		tf->push_back(mmop);
		
		PackRotamersMoverOP packer = new PackRotamersMover(core::scoring::getScoreFunction(true));
		packer->task_factory(tf);
		packer->apply(pose);
	}
	
};

int main(int argc, char* argv[]){
	
	devel::init(argc, argv);

	
	try{
		protocols::jd2::JobDistributor::get_instance()->go(new TestMMOP);
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cout << "Exception: " << std::endl;
		excn.show( std::cerr );
	}
	
	return(0);
}
