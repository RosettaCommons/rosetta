// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/docking/DockingTest.cxxtest.hh
/// @brief  tests for container Docking Movers classes.
/// @author Sergey Lyskov

// Test headers
#include <test/UTracer.hh>

// Unit headers
#include <protocols/docking/DockFilters.hh>
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/DockingLowRes.fwd.hh>
#include <protocols/docking/DockingHighRes.fwd.hh>
#include <protocols/docking/EllipsoidalRandomizationMover.hh>

// Package headers
#include <protocols/docking/metrics.hh>
#include <protocols/docking/util.hh>

// Project headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_task_operations/RestrictToInterface.hh>

#include <basic/options/option.hh> //needed to set option
// option key includes
#include <basic/options/keys/docking.OptionKeys.gen.hh>

// Utility headers

#include <core/import_pose/import_pose.hh> // AUTO IWYU For pose_from_file, PDB_file
#include <core/init_util.hh> // AUTO IWYU For core_init_with_additional_options
#include <core/kinematics/MoveMap.hh> // AUTO IWYU For MoveMap
#include <core/pose/Pose.hh> // AUTO IWYU For Pose


///////////////////////////////////////////////////////////////////////////
/// @name DockingTest
/// @brief: tests functions specific to the Docking Protocol
/// @author Sid Chaudhury
///////////////////////////////////////////////////////////////////////////
class DockingTest : public CxxTest::TestSuite {

	core::pose::Pose fullatom_pose;
	core::pose::Pose centroid_pose;
	core::Size rb_jump;

public:
	void setUp() {
		core_init_with_additional_options("-out:file:no_chainend_ter -INTEGRATION_TEST");
		core::import_pose::pose_from_file( fullatom_pose, "protocols/docking/DockingTest.pdb" , core::import_pose::PDB_file);
		rb_jump = 1;

		//setting up the fold tree as is used in docking
		core::kinematics::FoldTree fold_tree;

		core::Size jump_pos1 = 197;
		core::Size jump_pos2 = 282;
		core::Size cutpoint = 245;

		fold_tree.clear();
		fold_tree.add_edge( jump_pos1, jump_pos2, rb_jump );
		fold_tree.add_edge( 1, jump_pos1, core::kinematics::Edge::PEPTIDE );
		fold_tree.add_edge( jump_pos1, cutpoint, core::kinematics::Edge::PEPTIDE );
		fold_tree.add_edge( jump_pos2, cutpoint+1, core::kinematics::Edge::PEPTIDE );
		fold_tree.add_edge( jump_pos2, fullatom_pose.size(), core::kinematics::Edge::PEPTIDE );
		fold_tree.reorder( 1 );

		fullatom_pose.fold_tree(fold_tree);

		centroid_pose = fullatom_pose;

		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
		to_centroid.apply(centroid_pose);
	}

	void tearDown() {
		fullatom_pose.clear();
		centroid_pose.clear();
	}

	/// @brief test the docking protocol functions
	void test_DockingProtocolFunctions() {

		using protocols::docking::DockingProtocolOP;
		using protocols::docking::DockingProtocol;
		DockingProtocolOP docking_protocol( new DockingProtocol() ); //Defaults to rb_jump = 1, consistent with the value used

		test::UTracer UT("protocols/docking/DockingProtocolFunctions.u");

		core::scoring::ScoreFunctionOP scorefxn_low, scorefxn_high ;
		core::pose::Pose decoy_pose, multichain_pose;

		UT << "Testing docking low-res scoring function..."<<std::endl;
		scorefxn_low = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" ) ;
		scorefxn_low->show(UT, centroid_pose);

		UT << "Testing docking high-res scoring function..."<<std::endl;
		scorefxn_high = core::scoring::ScoreFunctionFactory::create_score_function( "docking" ) ;
		scorefxn_high->show(UT, fullatom_pose);

		UT << "Testing DockingProtocol.setup_foldtree()..."<< std::endl;
		protocols::docking::setup_foldtree( fullatom_pose, docking_protocol->partners(), docking_protocol->movable_jumps() );
		UT << fullatom_pose.fold_tree() << std::endl;

		UT << "Testing DockingProtocol.setup_foldtree()for multichain..."<< std::endl;
		core::import_pose::pose_from_file( multichain_pose, "protocols/docking/DockingMultiChain.pdb" , core::import_pose::PDB_file);
		DockingProtocolOP docking_protocol2( new DockingProtocol() );
		core::pose::DockingPartners interface = core::pose::DockingPartners::docking_partners_from_string("AB_E");
		protocols::docking::setup_foldtree( multichain_pose, interface, docking_protocol2->movable_jumps() );
		UT << multichain_pose.fold_tree() << std::endl;

		UT << "Testing interface-dependant scoring for multichain docking..."<<std::endl;
		protocols::simple_moves::SwitchResidueTypeSetMover to_centroid( core::chemical::CENTROID );
		to_centroid.apply(multichain_pose);
		scorefxn_low->show(UT, multichain_pose);

		UT << "Testing DockingProtocol.calc_interaction_energy()..."<<std::endl;
		core::Real int_energy = protocols::docking::calc_interaction_energy(fullatom_pose, scorefxn_high, docking_protocol->movable_jumps() );

		UT << int_energy << std::endl;

		UT << "Testing DockingProtocol.recover_sidechains()..."<<std::endl;
		core::import_pose::pose_from_file( decoy_pose, "protocols/docking/DockingDecoy.pdb" , core::import_pose::PDB_file);
		protocols::docking::setup_foldtree( decoy_pose, docking_protocol->partners(), docking_protocol->movable_jumps() );

		UT << "Testing DockingProtocol.calc_Lrmsd()..."<<std::endl;

		UT << protocols::docking::calc_Lrmsd(fullatom_pose, decoy_pose, docking_protocol->movable_jumps() ) << std::endl;

		UT << "Testing DockingProtocol.calc_Irmsd()..."<<std::endl;

		UT << protocols::docking::calc_Irmsd(fullatom_pose, decoy_pose, scorefxn_high, docking_protocol->movable_jumps() ) << std::endl;

		UT << "Testing DockingProtocol.calc_Fnat()..."<<std::endl;

		UT << protocols::docking::calc_Fnat(fullatom_pose, decoy_pose, scorefxn_high, docking_protocol->movable_jumps() ) << std::endl;

		UT << "Testing DockingProtocol.docking_lowres_filter()..."<<std::endl;

		protocols::docking::DockingLowResFilter lowres_filter;
		UT << lowres_filter.apply( centroid_pose ) << std::endl;

		protocols::docking::DockingHighResFilter highres_filter;
		highres_filter.set_score_cutoff( 1000000.0 );
		highres_filter.set_scorefunction( scorefxn_high );
		highres_filter.apply( fullatom_pose );

		UT << "Testing DockingProtocol.docking_highres_filter()..."<<std::endl;
		UT << highres_filter.apply( fullatom_pose ) << std::endl;
		// Note: we really should add logical unit tests for these latter two filters to test the logical cases
	}

	void test_DockingPacking() {
		using protocols::minimization_packing::PackRotamersMover;
		using protocols::minimization_packing::PackRotamersMoverOP;

		core::scoring::ScoreFunctionOP scorefxn_pack = core::scoring::get_score_function();
		core::scoring::ScoreFunctionOP scorefxn_dockmin = core::scoring::ScoreFunctionFactory::create_score_function("docking", "docking_min");
		(*scorefxn_pack)(fullatom_pose);

		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		using namespace protocols::simple_task_operations;

		TaskFactoryOP tf( new TaskFactory );
		tf->push_back( utility::pointer::make_shared< InitializeFromCommandline >() );
		tf->push_back( utility::pointer::make_shared< IncludeCurrent >() );
		tf->push_back( utility::pointer::make_shared< RestrictToRepacking >() );
		tf->push_back( utility::pointer::make_shared< RestrictToInterface >( rb_jump ) );

		core::pack::task::PackerTaskOP task = tf->create_task_and_apply_taskoperations( fullatom_pose );

		core::Real temperature = 0.8;
		protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo(fullatom_pose, *scorefxn_dockmin, temperature) );

		protocols::minimization_packing::PackRotamersMoverOP pack_interface_repack( new protocols::minimization_packing::PackRotamersMover( scorefxn_pack, task ) );

		pack_interface_repack->apply(fullatom_pose);

		mc->reset(fullatom_pose);

		core::Real energy_cut = 0.01;
		protocols::minimization_packing::RotamerTrialsMoverOP rotamer_trials( new protocols::minimization_packing::EnergyCutRotamerTrialsMover( scorefxn_dockmin, *task, mc, energy_cut ) );
		rotamer_trials->apply(fullatom_pose);

		test::UTracer UT("protocols/docking/DockingPacking.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		fullatom_pose.dump_pdb(UT);

	}

	void test_DockingSlideIntoContact() {
		using protocols::docking::DockingSlideIntoContact;

		DockingSlideIntoContact slide( rb_jump );

		protocols::rigid::RigidBodyTransMoverOP trans_mover( new protocols::rigid::RigidBodyTransMover(centroid_pose, rb_jump, false ) );
		trans_mover->step_size(10.26);
		trans_mover->apply(centroid_pose);
		slide.apply( centroid_pose );
		test::UTracer UT("protocols/docking/DockingSlideIntoContact.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		centroid_pose.dump_pdb(UT);
	}

	void test_DockingRigidBodyMinimize() {
		using protocols::minimization_packing::MinMover;
		using protocols::minimization_packing::MinMoverOP;

		core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" ) ;

		core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
		movemap->set_chi(false);
		movemap->set_bb(false);
		movemap->set_jump(rb_jump, true);

		core::Real tolerance = 0.01;
		std::string min_type = "dfpmin_armijo_nonmonotone";
		bool nb_list = true;

		protocols::minimization_packing::MinMoverOP minmover( new protocols::minimization_packing::MinMover(movemap, scorefxn, min_type, tolerance, nb_list) );
		minmover->apply(fullatom_pose);
		test::UTracer UT("protocols/docking/DockingRigidBodyMinimize.pdb");
		UT.abs_tolerance(0.003);
		UT << std::endl;
		fullatom_pose.dump_pdb(UT);
	}

	void test_DockingProtocol_clone(){
		protocols::docking::DockingProtocolOP dockerprot1( new protocols::docking::DockingProtocol() ) ;
		protocols::docking::DockingProtocolOP dockerprot2 = utility::pointer::static_pointer_cast< protocols::docking::DockingProtocol > ( dockerprot1->clone() );
		TS_ASSERT( dockerprot1.get() != dockerprot2.get() );
		TS_ASSERT( dockerprot1->to_centroid() != dockerprot2->to_centroid() );
		TS_ASSERT( dockerprot1->docking_lowres_mover() != dockerprot2->docking_lowres_mover() );
		TS_ASSERT( dockerprot1->to_all_atom() != dockerprot2->to_all_atom() );
		TS_ASSERT( dockerprot1->docking_highres_mover() != dockerprot2->docking_highres_mover() );
		TS_ASSERT( dockerprot1->perturber() != dockerprot2->perturber() );
		TS_ASSERT( dynamic_cast< protocols::simple_moves::SwitchResidueTypeSetMover const * > ( dockerprot2->to_centroid().get() ) );
	}

	void test_EllipsoidalRandomizationMover() {
		using protocols::docking::EllipsoidalRandomizationMover;

		basic::options::option[ basic::options::OptionKeys::docking::partners ].value( "E_I" );

		protocols::docking::EllipsoidalRandomizationMoverOP mover( new EllipsoidalRandomizationMover( rb_jump, false, true ) );
		mover->apply( fullatom_pose );

		core::Real a_axis = mover->get_a_axis();
		core::Real b_axis = mover->get_b_axis();
		core::Real c_axis = mover->get_c_axis();

		core::Vector c_alpha_centroid = mover->get_nbr_atom_centroid();
		core::Vector c_alpha_plane_centroid = mover->get_spin_center();
		core::Vector c_alpha_non_ellipsoid_centroid = mover->get_nbr_atom_non_ellipsoid_centroid();

		core::Vector random_point_on_ellipsoid = mover->get_random_point_on_ellipsoid();
		core::Real should_be_one = ( random_point_on_ellipsoid.x() * random_point_on_ellipsoid.x() ) / ( a_axis * a_axis ) +
			( random_point_on_ellipsoid.y() * random_point_on_ellipsoid.y() ) / ( b_axis * b_axis ) +
			( random_point_on_ellipsoid.z() * random_point_on_ellipsoid.z() ) / ( c_axis * c_axis );

		core::Vector ellipsoid_normal = mover->normal_to_ellipsoid( random_point_on_ellipsoid );
		core::Real should_be_one_too = ( ellipsoid_normal.x() * ellipsoid_normal.x()) +
			( ellipsoid_normal.y() * ellipsoid_normal.y()) +
			( ellipsoid_normal.z() * ellipsoid_normal.z());

		core::Vector normal_to_plane = mover->get_normal_to_plane();
		core::Vector slide_axis = mover->get_slide_axis();

		core::Real beta_sample = mover->single_beta_sample( 1.0, 1.0 );

		core::Real tolerance = 0.001;

		TS_ASSERT_DELTA( a_axis, 11.586, tolerance );
		TS_ASSERT_DELTA( b_axis, 8.331, tolerance );
		TS_ASSERT_DELTA( c_axis, 15.997, tolerance );
		TS_ASSERT_DELTA( c_alpha_centroid.x(), 44.094, tolerance );
		TS_ASSERT_DELTA( c_alpha_centroid.y(), 33.088, tolerance );
		TS_ASSERT_DELTA( c_alpha_centroid.z(), 33.046, tolerance );
		TS_ASSERT_DELTA( c_alpha_plane_centroid.x(), 31.114, tolerance );
		TS_ASSERT_DELTA( c_alpha_plane_centroid.y(), 29.683, tolerance );
		TS_ASSERT_DELTA( c_alpha_plane_centroid.z(), 29.135, tolerance );
		TS_ASSERT_DELTA( c_alpha_non_ellipsoid_centroid.x(), 23.415, tolerance );
		TS_ASSERT_DELTA( c_alpha_non_ellipsoid_centroid.y(), 26.422, tolerance );
		TS_ASSERT_DELTA( c_alpha_non_ellipsoid_centroid.z(), 22.750, tolerance );
		TS_ASSERT_DELTA( should_be_one, 1.000, tolerance );
		TS_ASSERT_DELTA( should_be_one_too, 1.000, tolerance );
		TS_ASSERT_DELTA( normal_to_plane.x(), -1.0 * slide_axis.x(), tolerance );
		TS_ASSERT_DELTA( normal_to_plane.y(), -1.0 * slide_axis.y(), tolerance );
		TS_ASSERT_DELTA( normal_to_plane.z(), -1.0 * slide_axis.z(), tolerance );
		TS_ASSERT( beta_sample >= 0.0 );
		TS_ASSERT( beta_sample <= 1.0 );

	}
};
