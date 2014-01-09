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


#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/InterfaceFeatures.hh>
#include <protocols/features/StructureFeatures.hh>
#include <protocols/features/ResidueFeatures.hh>
#include <protocols/features/PdbDataFeatures.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/sql_database/types.hh>
#include <utility/file/file_sys_util.hh>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <limits>

// External Headers
#include <cppdb/frontend.h>

//Auto Headers
#include <utility/vector1.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>


basic::options::StringOptionKey const interface("interface");
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


static basic::Tracer TR("protocols.features.InterfaceFeaturesTests.cxxtest");

using namespace protocols::features;

class InterfaceFeaturesTests : public protocols::moves::Mover{

public:
	//Unit test times out - limit not long enough to test.
	InterfaceFeaturesTests(){};
	
	virtual ~InterfaceFeaturesTests(){};
	
	virtual
	std::string
	get_name() const {
		return "InterfaceFeaturesTests";
	}
	
	virtual
	void
	apply(core::pose::Pose & pose){
		
		//setUp();
		//test_interfaces();
		//tearDown();
		
		setUp();
		test_reporter(pose);
		tearDown();
		

		

	}
	//"/home/jadolfbr/Documents/modeling/databases/antibody_databases/PyIgClassify/DBOUT/renumbered_pdbs/1pg7L-1pg7H.pdb"
	void setUp() {
		//core::import_pose::pose_from_pdb(multimer_, );
		db_name_ =  "InterfaceFeaturesTest2.db3";
		reporter_ = new protocols::features::InterfaceFeatures();
		utility::file::file_delete(db_name_);
		db_session_ = basic::database::get_db_session(db_name_);
		TR <<"Setup"<<std::endl;
	}
	void tearDown() {
		multimer_.clear();
	}
	void test_interfaces() {
		TR << "Testing Interface combos" << std::endl;
		utility::vector1<std::string> interfaces;
		reporter_->make_interface_combos(multimer_, interfaces);
		TR << "Interfaces: "<< interfaces.size() << std::endl;
		for (core::Size i = 1; i<=interfaces.size(); ++i){
			TR<<"Interface: " << interfaces[i]<<std::endl;
		}
	}
	void test_reporter(core::pose::Pose & pose) {
		TR << "Testing features reporter" << std::endl;
		utility::vector1<bool> relavant_residues(pose.total_residue(), true);
		
		reporter_->set_dSASA_cutoff(150); //Not real value
		reporter_->set_pack_separated(false); //speed
		reporter_->set_pack_together(false);
		
		std::string i = basic::options::option[ interface ].value();
		vector1<std::string> interfaces;
		//interfaces.push_back("L_A");
		interfaces.push_back(i);
		//interfaces.push_back("H_A");
		//interfaces.push_back("LH_A");
		
		reporter_->set_interface_chains(interfaces);
		reporter_->set_dSASA_cutoff(150);
		reporter_->set_pack_separated(false);
		reporter_->set_pack_together(true);
		
		StructureFeaturesOP structure_reporter = new StructureFeatures();
		structure_reporter->write_schema_to_db(db_session_);
		StructureID parent_id = structure_reporter->report_features(0, db_session_);
		
		TR << "Reported StructureFeatures" << std::endl;
		ResidueFeaturesOP residue_reporter = new ResidueFeatures();
		residue_reporter->write_schema_to_db(db_session_);
		residue_reporter->report_features(pose, relavant_residues, parent_id, db_session_);
		
		TR << "Reported Residue Features" << std::endl;
		PdbDataFeaturesOP pdb_info_features = new PdbDataFeatures();
		pdb_info_features->write_schema_to_db(db_session_);
		pdb_info_features->report_features(pose, relavant_residues, parent_id, db_session_);
		
		TR << "Reported PDBData Features" << std::endl;
		
		reporter_->write_schema_to_db(db_session_);
		reporter_->report_features(pose, relavant_residues, parent_id, db_session_);
		TR<<"Done reporting features" << std::endl;
		
	}
private:
	core::pose::Pose multimer_;
	InterfaceFeaturesOP reporter_;
	std::string db_name_;
	utility::sql_database::sessionOP db_session_;
	
};



int main(int argc, char* argv[]){
	
	

	
	try{
		using basic::options::option;
		option.add( interface, "dock_chains interface definition, optional, ex LH_A.  Can handle any number of chains. ");
		devel::init(argc, argv);
		protocols::jd2::JobDistributor::get_instance()->go(new InterfaceFeaturesTests);
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cout << "Exception: " << std::endl;
		excn.show( std::cerr );
	}
	
	return(0);
}
