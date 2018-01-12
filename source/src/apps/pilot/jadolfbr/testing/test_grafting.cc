// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/test_cart_graft.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/constraints/util.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/grafting/util.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/loops_main.hh>

#include <protocols/minimization_packing/MinMover.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/string_util.hh>
#include <numeric/random/random.hh>

static numeric::random::RandomGenerator RG(1365422);
static basic::Tracer TR("GraftingTest");


namespace myspace {


using namespace protocols::antibody;
using namespace protocols::antibody::design;

class GraftTester : public protocols::moves::Mover{

public:
	//Unit test times out - limit not long enough to test.
	GraftTester(){
		mode_ = "rb_min";
		rounds_ = 1;
	};

	virtual ~GraftTester(){};

	virtual
	std::string
	get_name() const {
		return "GraftTester";
	}

	virtual
	void
	apply(core::pose::Pose & pose){

		//setUp();
		//test_interfaces();
		//tearDown();

		init(pose);
		if ( mode_ == "rb_min" ) {
			test_rb_min_graft(pose);
		} else if ( mode_ == "cart" ) {
			test_cart_graft(pose);
		}
	}

	void
	set_cdrs(utility::vector1<CDRNameEnum> cdrs){
		cdrs_ = cdrs;
	}

	void
	test_cart_graft(core::pose::Pose & start_pose) {
		for ( core::Size i = 1; i <= rounds_; ++i ) {

			core::pose::Pose pose = start_pose;
			TR << "Graft round: " << i <<std::endl;

			CDRNameEnum cdr = cdrs_[RG.random_range(1, cdrs_.size())];
			if ( cdr_set_[cdr].size() ==0 ) {
				continue;
			}
			core::Size cdr_index = RG.random_range(1, cdr_set_[cdr].size());

			CDRDBPose c_p = cdr_set_[cdr][cdr_index];

			core::pose::Pose piece = *(c_p.pose);
			TR <<"Grafting: " << ab_info_->get_cluster_name(c_p.cluster)<<" from "<<c_p.pdb << std::endl;

			protocols::antibody::design::insert_cdr_into_antibody(ab_info_, cdr, pose, piece);

			modeler_->set_cdr_only(cdr, true);
			modeler_->cdr_overhang(cdr, 1);

			modeler_->minimize_cdrs(pose, true, true, false, true, false);

			std::pair<bool, core::Size> cb = protocols::loops::has_severe_pep_bond_geom_issues(pose, modeler_->get_cdr_loop_with_overhang(pose, cdr), true, true, 1.5, 15, 15);

			TR <<"CB "<<cb.first;

			pose.dump_pdb("cart_cb_"+utility::to_string(cb.first)+"_"+utility::to_string(i)+".pdb");
		}
	}
	void
	test_rb_min_graft(core::pose::Pose & start_pose){
		using namespace protocols::simple_moves;
		using namespace core::kinematics;

		core::pose::Pose pose = start_pose;

		CDRNameEnum cdr = cdrs_[RG.random_range(1, cdrs_.size())];
		if ( cdr_set_[cdr].size() ==0 ) {
			utility_exit_with_message("NO cdrs in CDR set!");
		}
		core::Size cdr_index = RG.random_range(1, cdr_set_[cdr].size());

		CDRDBPose c_p = cdr_set_[cdr][cdr_index];

		core::pose::Pose piece = *(c_p.pose);
		TR <<"Grafting: " << ab_info_->get_cluster_name(c_p.cluster)<<" from "<<c_p.pdb << std::endl;

		protocols::antibody::design::insert_cdr_into_antibody(ab_info_, cdr, pose, piece);

		pose.dump_pdb("inserted_pose.pdb");

		MinMoverOP min_mover = new MinMover();
		MoveMapOP movemap = new MoveMap();

		core::Size start = ab_info_->get_CDR_start(cdr, pose);
		core::Size end = ab_info_->get_CDR_end(cdr, pose);

		for ( core::Size i = start; i <= end; ++i ) {
			movemap->set_bb(i, true);
			//movemap->set_chi(i, true);
		}

		movemap->set_bb( start - 1, true );
		movemap->set_bb( start - 2, true );
		movemap->set_bb( start - 3, true);
		movemap->set_bb( end + 1, true );
		movemap->set_bb( end + 2, true);
		movemap->set_bb( end + 3, true);
		movemap->set_jump(true);

		protocols::loops::LoopsOP connections_ = new protocols::loops::Loops();
		protocols::loops::Loop term1 = protocols::loops::Loop(start - 2, start + 2, start-1);
		protocols::loops::Loop term2 = protocols::loops::Loop(end - 2, end + 2, end);
		connections_->add_loop(term1);
		connections_->add_loop(term2);

		core::kinematics::FoldTree original_ft = pose.fold_tree();
		core::kinematics::FoldTree new_ft = pose.fold_tree();
		protocols::loops::fold_tree_from_loops(pose, *connections_, new_ft);
		pose.fold_tree(new_ft);
		protocols::loops::add_cutpoint_variants(pose);
		min_mover->movemap(movemap);
		min_mover->score_function(scorefxn_);
		min_mover->apply(pose);
		pose.fold_tree(original_ft);
		start_pose = pose;

	}

private:


	void
	init(core::pose::Pose & start_pose) {
		ab_info_ = new AntibodyInfo(start_pose, AHO_Scheme, North);
		modeler_ = new AntibodyDesignModeler(ab_info_);
		manager_ = new AntibodyDatabaseManager(ab_info_);

		scorefxn_ = core::scoring::get_score_function();
		scorefxn_->set_weight(chainbreak, 1.0);
		scorefxn_->set_weight(dihedral_constraint, 1.0);
		scorefxn_->set_weight(coordinate_constraint, 1.0);
		modeler_->set_scorefunction(scorefxn_);

		TR << "Classes set" << std::endl;
		cdr_options_ = protocols::antibody::design::get_cdr_set_options();
		TR << "Options set " << std::endl;
		cdr_set_ = manager_->load_cdr_poses(cdr_options_, start_pose, true /*use light chain type*/);
		TR<< "CDRs loaded" << std::endl;
		utility::vector1<CDRNameEnum> cdrs;
		for ( core::Size i = 1; i <= 6; ++i ) {
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			if ( cdr_set_[cdr].size() >=1 ) {
				cdrs_.push_back(cdr);
			}
		}


	}

private:
	AntibodyInfoOP ab_info_;
	AntibodyDesignModelerOP modeler_;
	AntibodyDatabaseManagerOP manager_;
	utility::vector1<CDRSetOptionsOP> cdr_options_;
	utility::vector1<CDRNameEnum> cdrs_;
	core::scoring::ScoreFunctionOP scorefxn_;
	std::string mode_;
	core::Size rounds_;
	std::map<CDRNameEnum, utility::vector1<protocols::antibody::CDRPose> > cdr_set_;
};

} //namespace myspace


int main(int argc, char* argv[]){

	try{
		using basic::options::option;
		//option.add( mode, "rb_min or cart ");
		devel::init(argc, argv);
		protocols::jd2::JobDistributor::get_instance()->go(new myspace::GraftTester());
	} catch (utility::excn::Exception& excn ) {
		std::cout << "Exception: " << std::endl;
		excn.show( std::cerr );
		return -1;
	}

	return(0);
}
