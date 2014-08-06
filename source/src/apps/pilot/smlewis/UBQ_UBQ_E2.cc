// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/apps/pilot/smlewis/UBQ_UBQ_E2.cc
/// @brief  this application is a one-shot for modeling the thioester bond between UBQ and an E2
/// @author Steven Lewis

// Unit Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>

//movers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh> //Sequence Mover
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/moves/OutputMovers.hh> //pdbdumpmover
#include <protocols/simple_moves/TorsionDOFMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/rigid/RotateJumpAxisMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>

#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>

// Numeric Headers
#include <numeric/conversions.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <basic/prof.hh>

#include <numeric/xyz.io.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

//local options
basic::options::FileOptionKey const UBQpdb("UBQpdb");
basic::options::FileOptionKey const E2pdb("E2pdb");
basic::options::IntegerOptionKey const E2_residue("E2_residue");
basic::options::RealOptionKey const SASAfilter("SASAfilter");
basic::options::RealOptionKey const scorefilter("scorefilter");

//tracers
using basic::Error;
using basic::Warning;
static basic::Tracer TR("apps.pilot.smlewis.UBQ_UBQ_E2");

class UBQ_E2Mover : public protocols::moves::Mover {
private: //enum for atomID vector
	enum atomID {
		CYS_C  = 1,
		CYS_CA = 2,
		CYS_CB = 3,
		CYS_SG = 4,
		GLY_C  = 5,
		GLY_CA = 6,
		GLY_N  = 7,
		GLY2_C = 8,
		LYS_2HZ= 9,
		LYS_NZ = 10,
		LYS_CE = 11,
		LYS_CD = 12,
		atomID_tot = LYS_CD
	};



public:
	UBQ_E2Mover()
	: init_for_input_yet_(false),
		fullatom_scorefunction_(NULL),
		task_factory_(NULL),
		thioester_mm_(NULL),
		loop_(), //we want default ctor
		atomIDs(atomID_tot, core::id::BOGUS_ATOM_ID ),
		InterfaceSasaDefinition_("InterfaceSasaDefinition_" + 1),
		IAM_(new protocols::analysis::InterfaceAnalyzerMover)
	{
		//set up fullatom scorefunction
		using namespace core::scoring;
		fullatom_scorefunction_ = get_score_function();
		TR << "Using fullatom scorefunction (TALARIS_2013), pair may be modified\n"
			 << *fullatom_scorefunction_;

		using namespace core::pose::metrics;
		using namespace protocols::toolbox::pose_metric_calculators;
		//magic number: chains 1 and 2; set up interface SASA calculator
		if( !CalculatorFactory::Instance().check_calculator_exists( InterfaceSasaDefinition_ ) ){
			CalculatorFactory::Instance().register_calculator( InterfaceSasaDefinition_, new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator(core::Size(1), core::Size(2)));
		}

		IAM_->set_use_centroid_dG(false);

		//set up loop object
		if ( basic::options::option[ basic::options::OptionKeys::loops::loop_file ].user() ) {
			protocols::loops::Loops loops( basic::options::option[ basic::options::OptionKeys::loops::loop_file ].value()[1] );
			loop_ = *(loops.begin());\
		}
	}

	///@brief init_on_new_input system allows for initializing these details the first time apply() is called.  the job distributor will reinitialize the whole mover when the input changes (a freshly constructed mover, which will re-run this on first apply().
	virtual
	void
	init_on_new_input() {
		init_for_input_yet_ = true;


		///////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////Create starting complex pose/////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		TR << "Creating starting pose..." << std::endl;

		//read poses
		core::pose::Pose E2;
		core::import_pose::pose_from_pdb( E2, basic::options::option[E2pdb].value() );
		core::Size const E2length = E2.total_residue();

		core::pose::Pose UBQ;
		core::import_pose::pose_from_pdb( UBQ, basic::options::option[UBQpdb].value() );
		core::Size const UBQlength = UBQ.total_residue();
		core::pose::Pose UBQ_second(UBQ);

		//determine cysteine target
		runtime_assert(E2.conformation().num_chains() == 1);
		char const E2chain(E2.pdb_info()->chain(1));
		core::Size const E2_cys(E2.pdb_info()->pdb2pose(E2chain, basic::options::option[E2_residue].value()));
		runtime_assert(E2.residue_type(E2_cys).aa() == core::chemical::aa_cys);

		//determine c_term target on UBQ
		core::Size const UBQ_term = UBQlength;

		//strip C-term from UBQ - best to do this with a full replace to re-draw the carboxyl oxygen
		//UBQ.dump_pdb("pre-removeUQB.pdb");
		core::chemical::ResidueTypeSetCAP fa_standard(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
		UBQ.conformation().delete_residue_slow( UBQ_term );
		UBQ.append_residue_by_bond( *(core::conformation::ResidueFactory::create_residue(fa_standard->name_map("GLY")) ) );
		//UBQ.dump_pdb("post-removeUQB.pdb");

		//replace cysteine
		core::chemical::ResidueType const & cyx_rsd_type( fa_standard->name_map("CYX") );
		//E2.dump_pdb("prereplace_E2.pdb");
		core::pose::replace_pose_residue_copying_existing_coordinates( E2, E2_cys, cyx_rsd_type );
		//E2.dump_pdb("postreplace_E2.pdb");


		// check safety of connections (from phil)
		core::chemical::ResidueType const & ubq_rsd_type( UBQ.residue_type( UBQ_term ) );
		core::Size const cyx_connid( 3 );
		core::Size const ubq_connid( 2 );

		runtime_assert( cyx_rsd_type.n_residue_connections() == cyx_connid &&
			cyx_rsd_type.lower_connect_id() != cyx_connid &&
			cyx_rsd_type.upper_connect_id() != cyx_connid );

		runtime_assert( ubq_rsd_type.n_residue_connections() == ubq_connid &&
			ubq_rsd_type.lower_connect_id() != ubq_connid);

		//E2.conformation().show_residue_connections();
		//UBQ.conformation().show_residue_connections();

		//cross fingers - making the actual connection!
		/*void
			append_residue_by_bond(
			conformation::Residue const & new_rsd,
			bool const build_ideal_geometry = false,
			int const connection = 0,
			Size const anchor_residue = 0,
			int const anchor_connection = 0,
			bool const start_new_chain = false
			)*/
		core::pose::Pose complex(E2);
		complex.append_residue_by_bond( UBQ.residue( UBQ_term ), true, ubq_connid, E2_cys, cyx_connid );
		//complex.dump_pdb("just1_complex.pdb");

		//not that this does anything
		complex.conformation().insert_ideal_geometry_at_residue_connection( E2_cys, cyx_connid );

		core::Size const ubq_pos( complex.total_residue() );
		core::id::AtomID const atom0( cyx_rsd_type.atom_index( "C" ), E2_cys );
		core::id::AtomID const atom1( cyx_rsd_type.atom_index( "CA" ), E2_cys );
		core::id::AtomID const atom2( cyx_rsd_type.atom_index( "CB" ), E2_cys );
		core::id::AtomID const atom3( cyx_rsd_type.atom_index( "SG" ), E2_cys );
		core::id::AtomID const atom4( ubq_rsd_type.atom_index( "C"  ), ubq_pos );
		core::id::AtomID const atom5( ubq_rsd_type.atom_index( "CA" ), ubq_pos );
		core::id::AtomID const atom6( ubq_rsd_type.atom_index( "N"  ), ubq_pos );


		//starting values derived from 1FXT.pdb
		complex.conformation().set_torsion_angle( atom0, atom1, atom2, atom3, numeric::conversions::radians(106.5) );
		complex.conformation().set_torsion_angle( atom1, atom2, atom3, atom4, numeric::conversions::radians(-60.0) );
		complex.conformation().set_torsion_angle( atom2, atom3, atom4, atom5, numeric::conversions::radians(180.0) );
		complex.conformation().set_torsion_angle( atom3, atom4, atom5, atom6, numeric::conversions::radians(100.5) );
		//complex.dump_pdb("just1_complex2.pdb");

		//now add the rest of ubiquitin
		for( core::Size i=UBQ_term-1; i>= 1; --i ) {
			complex.prepend_polymer_residue_before_seqpos( UBQ.residue(i), E2length+1, false );
		}

		core::Size const complexlength( complex.total_residue());
		complex.conformation().insert_ideal_geometry_at_polymer_bond( complexlength-1 );
		complex.conformation().insert_chain_ending(E2length);
		//complex.dump_pdb("initcomplex.pdb");


		//now add in the second ubiquitin - link to thioester C=O from lys48
		complex.conformation().append_residue_by_jump(UBQ_second.residue(48), complexlength, "C", "NZ", true);
		for( core::Size i(49); i <= UBQlength; ++i) complex.conformation().append_polymer_residue_after_seqpos(UBQ_second.residue(i), complex.total_residue(), false);
		for( core::Size i(47); i >= 1; --i) complex.conformation().prepend_polymer_residue_before_seqpos(UBQ_second.residue(i), complexlength+1, false);

		//check it!
		TR << complex.fold_tree();

		//the combined pose is complete - time to pack the atomID vector
		atomIDs[1] = atom0;
		atomIDs[2] = atom1;
		atomIDs[3] = atom2;
		atomIDs[4] = atom3;
		atomIDs[5] = core::id::AtomID( ubq_rsd_type.atom_index("C" ), complexlength );
		atomIDs[6] = core::id::AtomID( ubq_rsd_type.atom_index("CA" ), complexlength );
		atomIDs[7] = core::id::AtomID( ubq_rsd_type.atom_index("N" ), complexlength );
		atomIDs[8] = core::id::AtomID( ubq_rsd_type.atom_index("C" ), complexlength-1 );

		core::chemical::ResidueType const & lys_rsd_type( fa_standard->name_map("LYS") );

		atomIDs[LYS_2HZ] = core::id::AtomID( lys_rsd_type.atom_index("2HZ"), complexlength + 48);
		atomIDs[LYS_NZ]  = core::id::AtomID( lys_rsd_type.atom_index("NZ" ), complexlength + 48);
		atomIDs[LYS_CE]  = core::id::AtomID( lys_rsd_type.atom_index("CE" ), complexlength + 48);
		atomIDs[LYS_CD]  = core::id::AtomID( lys_rsd_type.atom_index("CD" ), complexlength + 48);

		//ok, jump is in place (numbered 1) - now we have to make a statement about where to put it
		//relevant atoms for placing NZ
		core::Vector const & C_xyz(complex.residue(complexlength).atom("C").xyz());
		core::Vector const & O_xyz(complex.residue(complexlength).atom("O").xyz());
		//TR << "C " << C_xyz << " O " << O_xyz << std::endl;
		core::Vector newpos(C_xyz+(3*(C_xyz-O_xyz)/O_xyz.distance(C_xyz)));
		core::Vector oldpos(complex.residue(complexlength+48).atom("NZ").xyz());

		//core::Vector oldposC(complex.residue(complexlength+48).atom("CE").xyz());
		//core::Vector newposC(newpos+(oldpos-oldposC));

		complex.set_xyz(core::id::AtomID(lys_rsd_type.atom_index("NZ"), complexlength+48), newpos);
		//	complex.set_xyz(core::id::AtomID(lys_rsd_type.atom_index("CE"), complexlength+48), newposC);

		core::kinematics::Jump newjump(complex.atom_tree().jump(core::id::AtomID(lys_rsd_type.atom_index("NZ"), complexlength+48)));

		complex.set_xyz(core::id::AtomID(lys_rsd_type.atom_index("NZ"), complexlength+48), oldpos);
		//	complex.set_xyz(core::id::AtomID(lys_rsd_type.atom_index("CE"), complexlength+48), oldposC);

		complex.conformation().set_jump(1, newjump);

		starting_pose_ = complex;
		starting_pose_.dump_pdb("starting_complex.pdb");
		TR << "starting pose finished." << std::endl;

		///////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////Finish starting complex pose/////////////////////////////////////////////
		//////////////////////Start creating move accessory data///////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		//setup MoveMaps
		//small/shear behave improperly @ the last residue - psi is considered nonexistent and the wrong phis apply.
		thioester_mm_ = new core::kinematics::MoveMap;
		//thioester_mm_->set_bb(complexlength, true);
		//thioester_mm_->set(core::id::TorsionID(complexlength, core::id::BB, core::id::phi_torsion), true);
		//thioester_mm_->set(core::id::TorsionID(complexlength, core::id::BB, core::id::psi_torsion), true);
		thioester_mm_->set_bb(complexlength-1, true);
		thioester_mm_->set_bb(complexlength-2, true);
		//thioester_mm_->set(complex.atom_tree().torsion_angle_dof_id(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5]), false);


		std::set< core::Size > loop_posns;
		//setup loop
		if ( basic::options::option[ basic::options::OptionKeys::loops::loop_file ].user() ) {
			TR << "loop " <<  loop_ << std::endl;
			//set up for interface-plus-neighbors-positions operation
			for (core::Size j(loop_.start()), end(loop_.stop()); j <= end; ++j){
				loop_posns.insert(j);
			}//for each residue in loop
		} //no else needed - default loop is safe enough

		//setup of TaskFactory
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		task_factory_ = new TaskFactory;
		task_factory_->push_back( new InitializeFromCommandline );
		if ( basic::options::option[ basic::options::OptionKeys::packing::resfile ].user() ) {
			task_factory_->push_back( new ReadResfile );
		}
		//task_factory_->push_back( new protocols::toolbox::task_operations::RestrictToInterfaceOperation );
		task_factory_->push_back( new IncludeCurrent );
		//prevent repacking at linkage cysteine!
		PreventRepackingOP prevent(new PreventRepacking);
		prevent->include_residue(E2_cys);
		prevent->include_residue(complexlength+48);
		task_factory_->push_back(prevent);

		std::string const interface_calc12("UBQE2_InterfaceNeighborDefinitionCalculator12");
		std::string const interface_calc23("UBQE2_InterfaceNeighborDefinitionCalculator23");
		std::string const interface_calc13("UBQE2_InterfaceNeighborDefinitionCalculator13");
		std::string const neighborhood_calc("UBQE2_NeighborhoodByDistanceCalculator");
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( interface_calc12, new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( core::Size(1), core::Size(2)) );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( interface_calc23, new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( core::Size(2), core::Size(3)) );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( interface_calc13, new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( core::Size(1), core::Size(3)) );

		core::pose::metrics::CalculatorFactory::Instance().register_calculator( neighborhood_calc, new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( loop_posns ) );

		//this is the constructor parameter for the calculator - pairs of calculators and calculations to perform
		utility::vector1< std::pair< std::string, std::string> > calcs_and_calcns;
		calcs_and_calcns.push_back(std::make_pair(interface_calc12, "interface_residues"));
		calcs_and_calcns.push_back(std::make_pair(interface_calc23, "interface_residues"));
		calcs_and_calcns.push_back(std::make_pair(interface_calc13, "interface_residues"));
		calcs_and_calcns.push_back(std::make_pair(neighborhood_calc, "neighbors"));

		using protocols::toolbox::task_operations::RestrictByCalculatorsOperation;
		task_factory_->push_back(new RestrictByCalculatorsOperation( calcs_and_calcns ));

	}

	virtual ~UBQ_E2Mover(){};

	virtual
	void
	apply( core::pose::Pose & pose ){
		if( !init_for_input_yet_ ) init_on_new_input();

		pose = starting_pose_;

		//TR << "foldtree, movemap: " << std::endl;
		//core::kinematics::simple_visualize_fold_tree_and_movemap( pose.fold_tree(), *movemap_, TR);

		/////////////////fullatom Monte Carlo//////////////////////////////////////////////////////////
		//make the monte carlo object
		using protocols::moves::MonteCarlo;
		using protocols::moves::MonteCarloOP;
		using basic::options::option;
		using namespace basic::options::OptionKeys::AnchoredDesign;
		MonteCarloOP mc( new MonteCarlo( pose, *fullatom_scorefunction_, option[ refine_temp ].value() ) );

		//////////////////////////Small/ShearMovers////////////////////////////////////////////////////////
		protocols::simple_moves::BackboneMoverOP small_mover = new protocols::simple_moves::SmallMover(thioester_mm_, 0.8, 1);
		small_mover->angle_max( 'H', 4.0 );
		small_mover->angle_max( 'E', 4.0 );
		small_mover->angle_max( 'L', 4.0 );

		protocols::simple_moves::BackboneMoverOP shear_mover = new protocols::simple_moves::ShearMover(thioester_mm_, 0.8, 1);
		shear_mover->angle_max( 'H', 4.0 );
		shear_mover->angle_max( 'E', 4.0 );
		shear_mover->angle_max( 'L', 4.0 );

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_chi1(new protocols::simple_moves::TorsionDOFMover);
		DOF_mover_chi1->set_DOF(atomIDs[1], atomIDs[2], atomIDs[3], atomIDs[4]);
		DOF_mover_chi1->check_mmt(true);
		DOF_mover_chi1->temp(0.4);
		DOF_mover_chi1->set_angle_range(-180, 180);
		DOF_mover_chi1->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_chi2(new protocols::simple_moves::TorsionDOFMover);
		DOF_mover_chi2->set_DOF(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5]);
		DOF_mover_chi2->check_mmt(true);
		DOF_mover_chi2->temp(0.4);
		DOF_mover_chi2->set_angle_range(-180, 180);
		DOF_mover_chi2->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_thioester(new protocols::simple_moves::TorsionDOFMover);
		DOF_mover_thioester->set_DOF(atomIDs[3], atomIDs[4], atomIDs[5], atomIDs[6]);
		DOF_mover_thioester->check_mmt(true);
		DOF_mover_thioester->temp(0.4);
		DOF_mover_thioester->set_angle_range(-180, 180);
		DOF_mover_thioester->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_psi(new protocols::simple_moves::TorsionDOFMover);
		DOF_mover_psi->set_DOF(atomIDs[4], atomIDs[5], atomIDs[6], atomIDs[7]);
		DOF_mover_psi->check_mmt(true);
		DOF_mover_psi->temp(0.4);
		DOF_mover_psi->set_angle_range(-180, 180);
		DOF_mover_psi->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_phi(new protocols::simple_moves::TorsionDOFMover);
		DOF_mover_phi->set_DOF(atomIDs[5], atomIDs[6], atomIDs[7], atomIDs[8]);
		DOF_mover_phi->check_mmt(true);
		DOF_mover_phi->temp(0.4);
		DOF_mover_phi->set_angle_range(-180, 180);
		DOF_mover_phi->tries(1000);

		protocols::moves::RandomMoverOP backbone_mover( new protocols::moves::RandomMover() );
		backbone_mover->add_mover(small_mover, 2.0);
		backbone_mover->add_mover(shear_mover, 1.0);
		backbone_mover->add_mover(DOF_mover_chi1, 0.75);
		backbone_mover->add_mover(DOF_mover_chi2, 0.75);
		backbone_mover->add_mover(DOF_mover_thioester, 0.75);
		backbone_mover->add_mover(DOF_mover_psi, 0.75);
		backbone_mover->add_mover(DOF_mover_phi, 0.75);

		///////////////////////////loop movement/////////////////////////////////////////////////////
		if( loop_.stop() - loop_.start() >= 3 ) { //empty loop; skip it!
			//make kinematic mover
			using protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP;
			using protocols::loops::loop_closure::kinematic_closure::KinematicMover;
			KinematicMoverOP kin_mover( new KinematicMover() );
			kin_mover->set_temperature( 0.8 );
			kin_mover->set_vary_bondangles( true );
			kin_mover->set_sample_nonpivot_torsions( true );
			kin_mover->set_rama_check( true );

			using protocols::loops::loop_closure::kinematic_closure::KinematicWrapperOP;
			using protocols::loops::loop_closure::kinematic_closure::KinematicWrapper;
			KinematicWrapperOP kin_wrapper( new KinematicWrapper(kin_mover, loop_));

			backbone_mover->add_mover(kin_wrapper, 5);

		}

		//////////////////////////////RotateJumpAxisMover for second ubiquitin//////////////////////
		protocols::rigid::RotateJumpAxisMoverOP RJAmover(new protocols::rigid::RotateJumpAxisMover(1));
		backbone_mover->add_mover(RJAmover, 1);


		/////////////////////////////sidechainmover for second ubiquitin/////////////////////////////
		core::Size Kres(pose.fold_tree().jump_edge(1).stop());
		core::pack::task::PackerTaskOP task(core::pack::task::TaskFactory::create_packer_task((pose)));
		utility::vector1_bool packable(pose.total_residue(), false); //false = nobody is packable
		packable[Kres] = true;
		task->restrict_to_residues(packable); //now only one position is mobile
		task->restrict_to_repacking();
		task->nonconst_residue_task(Kres).or_ex1_sample_level(core::pack::task::EX_SIX_QUARTER_STEP_STDDEVS);
		task->nonconst_residue_task(Kres).or_ex2_sample_level(core::pack::task::EX_SIX_QUARTER_STEP_STDDEVS);
		task->nonconst_residue_task(Kres).or_ex3_sample_level(core::pack::task::EX_SIX_QUARTER_STEP_STDDEVS);
		task->nonconst_residue_task(Kres).or_ex4_sample_level(core::pack::task::EX_SIX_QUARTER_STEP_STDDEVS);
		//		TR << *task << std::endl;

		protocols::simple_moves::sidechain_moves::SidechainMoverOP SCmover( new protocols::simple_moves::sidechain_moves::SidechainMover() );
		SCmover->set_task(task);
		SCmover->set_prob_uniform(0); //we want only Dunbrack rotamers, 0 percent chance of uniform sampling

		backbone_mover->add_mover(SCmover, 1);

		//////////////////////////////added "chi" for lysine//////////////////////////////////////
		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_lys(new protocols::simple_moves::TorsionDOFMover);
		DOF_mover_lys->set_DOF(atomIDs[LYS_2HZ], atomIDs[LYS_NZ], atomIDs[LYS_CE], atomIDs[LYS_CD]);
		DOF_mover_lys->check_mmt(false);
		DOF_mover_lys->temp(10);
		DOF_mover_lys->set_angle_range(-180, 180);
		DOF_mover_lys->tries(1000);

		backbone_mover->add_mover(DOF_mover_lys, 1);

		/////////////////////////minimize backbone DOFs//////////////////////////////////////////////
		using protocols::simple_moves::MinMoverOP;
		using protocols::simple_moves::MinMover;
		protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover(
																				thioester_mm_,
																				fullatom_scorefunction_,
																				basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
																				0.01,
																				true /*use_nblist*/ );

		/////////////////////////////////rotamer trials mover///////////////////////////////////////////
		using protocols::simple_moves::RotamerTrialsMoverOP;
		using protocols::simple_moves::EnergyCutRotamerTrialsMover;
		protocols::simple_moves::RotamerTrialsMoverOP rt_mover(new protocols::simple_moves::EnergyCutRotamerTrialsMover(
																																	fullatom_scorefunction_,
																																	task_factory_,
																																	mc,
																																	0.01 /*energycut*/ ) );

		///////////////////////package RT/min for JumpOutMover////////////////////////////////////////
		protocols::moves::SequenceMoverOP RT_min_seq( new protocols::moves::SequenceMover );
		RT_min_seq->add_mover(rt_mover);
		RT_min_seq->add_mover(min_mover);

		protocols::moves::JumpOutMoverOP bb_if_RT_min( new protocols::moves::JumpOutMover(
																																											backbone_mover,
																																											RT_min_seq,
																																											fullatom_scorefunction_,
																																											20.0));

		///////////////////////////////repack///////////////////////////////////////////////
		protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;
		pack_mover->task_factory( task_factory_ );
		pack_mover->score_function( fullatom_scorefunction_ );

		protocols::simple_moves::MinMoverOP min_mover_pack = new protocols::simple_moves::MinMover(
																						 thioester_mm_,
																						 fullatom_scorefunction_,
																						 basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
																						 0.01,
																						 true /*use_nblist*/ );

		using protocols::simple_moves::TaskAwareMinMoverOP;
		using protocols::simple_moves::TaskAwareMinMover;
		protocols::simple_moves::TaskAwareMinMoverOP TAmin_mover = new protocols::simple_moves::TaskAwareMinMover(min_mover_pack, task_factory_);

		/////////////////////////////////////////refine loop///////////////////////////////////////////

		core::Size const refine_applies = option[ refine_cycles ].value(); //default 5
		core::Size const repack_cycles = option[ refine_repack_cycles ].value();
		//core::Size const min_cycles = repack_cycles/2;
		TR << "   Current     Low    total cycles =" << refine_applies << std::endl;
		for ( core::Size i(1); i <= refine_applies; ++i ) {

			if( (i % repack_cycles == 0) || (i == refine_applies) ) { //full repack
				pack_mover->apply(pose);
				TAmin_mover->apply(pose);
				//} else if ( i % min_cycles == 0 ) { //minimize
				//min_mover->apply(pose);
			} else {
				bb_if_RT_min->apply(pose);
			}


			mc->boltzmann(pose);
			TR << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
		}//end the exciting for loop
		mc->recover_low( pose );

		//filter on SASAs of 1000?
		(*fullatom_scorefunction_)(pose);
		analyze_and_filter(pose);
		return;
	}

	void
	analyze_and_filter(core::pose::Pose & pose){

		//Filter on total score
		core::Real const score((*fullatom_scorefunction_)(pose));
		if( score > basic::options::option[scorefilter].value() ){
			set_last_move_status(protocols::moves::FAIL_RETRY);
			TR << "total score filter failed; score " << score << std::endl;
			return;
		}

		//filter on interface SASA - requires some hacking to break up thioester
		core::pose::Pose copy(pose);
		//hack the pose up for analysis purposes
		core::Size const cbreak(copy.conformation().chain_end(1));
		using core::kinematics::Edge;
		core::kinematics::FoldTree main_tree(copy.total_residue());
		main_tree.clear();
		main_tree.add_edge(Edge(1, cbreak, Edge::PEPTIDE));
		main_tree.add_edge(Edge(cbreak+1, copy.total_residue(), Edge::PEPTIDE));
		main_tree.add_edge(Edge(cbreak, cbreak+1, 1));
		main_tree.reorder(1);
		//TR << main_tree << std::endl;
		copy.fold_tree(main_tree);

 		//Filter on SASA
		basic::MetricValue< core::Real > mv_delta_sasa;
		copy.metric(InterfaceSasaDefinition_, "delta_sasa", mv_delta_sasa);
		if(mv_delta_sasa.value() < basic::options::option[SASAfilter].value()){
			set_last_move_status(protocols::moves::FAIL_RETRY);
			TR << "interface SASA filter failed; SASA " << mv_delta_sasa.value() << std::endl;
			return;
		}

		//passed filters; run IAM
		IAM_->apply(copy);

		//print mobile region fine-grained data
		protocols::jd2::JobOP job_me(protocols::jd2::JobDistributor::get_instance()->current_job());
		using numeric::conversions::degrees;
		job_me->add_string_real_pair("cysteine_chi1_C-CA-CB-SG", degrees(pose.atom_tree().torsion_angle(atomIDs[1], atomIDs[2], atomIDs[3], atomIDs[4])));
		job_me->add_string_real_pair("cysteine_chi2_CA-CB-SG-C", degrees(pose.atom_tree().torsion_angle(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5])));
		job_me->add_string_real_pair("thioester_CB-SG-C-CA", degrees(pose.atom_tree().torsion_angle(atomIDs[3], atomIDs[4], atomIDs[5], atomIDs[6])));
		job_me->add_string_real_pair("glycine_psi_SG-C-CA-N", degrees(pose.atom_tree().torsion_angle(atomIDs[4], atomIDs[5], atomIDs[6], atomIDs[7])));
		job_me->add_string_real_pair("glycine_phi_C-CA-N-C", degrees(pose.atom_tree().torsion_angle(atomIDs[5], atomIDs[6], atomIDs[7], atomIDs[8])));

		set_last_move_status(protocols::moves::MS_SUCCESS);
		return;
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new UBQ_E2Mover;
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }

	virtual
	std::string
	get_name() const { return "UBQ_E2Mover"; }

private:
	bool init_for_input_yet_;

	core::scoring::ScoreFunctionOP fullatom_scorefunction_;
	core::pack::task::TaskFactoryOP task_factory_;
	core::kinematics::MoveMapOP thioester_mm_;
// 	core::kinematics::MoveMapOP loop_mm_;
// 	core::kinematics::MoveMapOP all_mm_;

	protocols::loops::Loop loop_;

	///@brief vector contains atomIDs for thioester bond and atoms before/after bond to determine various torsions
	utility::vector1< core::id::AtomID > atomIDs;

	core::pose::Pose starting_pose_; //maintained from run to run

	std::string const InterfaceSasaDefinition_; //calculator name

	protocols::analysis::InterfaceAnalyzerMoverOP IAM_;

};

typedef utility::pointer::owning_ptr< UBQ_E2Mover > UBQ_E2MoverOP;

int main( int argc, char* argv[] )
{

	try {


	using basic::options::option;
	using namespace basic::options::OptionKeys;
 	option.add( UBQpdb, "ubiquitin structure" ).def("1UBQ.pdb");
 	option.add( E2pdb, "E2 structure" ).def("2OB4.pdb");
 	option.add( E2_residue, "E2 catalytic cysteine (PDB numbering)").def(85);
	option.add( SASAfilter, "filter out interface dSASA less than this").def(1000);
	option.add( scorefilter, "filter out total score greater than this").def(10);

	//initialize options
	devel::init(argc, argv);
	basic::prof_reset();

	if(basic::options::option[ basic::options::OptionKeys::in::file::s ].active()
		|| basic::options::option[ basic::options::OptionKeys::in::file::l ].active()
		|| basic::options::option[ basic::options::OptionKeys::in::file::silent ].active())
		utility_exit_with_message("do not use an input PDB with UBQ_E2 (program uses internally)");

	protocols::jd2::JobDistributor::get_instance()->go(new UBQ_E2Mover);

	basic::prof_show();
	TR << "************************d**o**n**e**************************************" << std::endl;

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
