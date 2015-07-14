// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/apps/public/scenarios/chemically_conjugated_docking/UBQ_E2_thioester.cc
/// @brief  this application is a one-shot for modeling the thioester bond between UBQ and an E2
/// @author Steven Lewis

// Unit Headers
#include <apps/public/scenarios/chemically_conjugated_docking/Gp_extra_bodies.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>

#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/id/TorsionID.hh> //used for movemap to set omega angle false
#include <core/id/types.hh> //used for movemap to set omega angle false

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

//movers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh> //Sequence Mover
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
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
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/InterGroupNeighborsCalculator.hh>

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
#include <utility/vector0.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/chemically_conjugated_docking.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>
//tracers
using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "apps.public.scenarios.chemically_conjugated_docking.UBQ_E2_thioester" );

class UBQ_E2Mover : public protocols::moves::Mover {

private: //enum for atomID vector for tracking bonds near the thioester
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
		fullatom_scorefunction_(/* NULL */),
		task_factory_(/* NULL */),
		thioester_mm_(/* NULL */),
		loop_(), //we want default ctor
		atomIDs(8, core::id::BOGUS_ATOM_ID ),
		InterfaceSasaDefinition_("InterfaceSasaDefinition_" + 1),
		IAM_(protocols::analysis::InterfaceAnalyzerMoverOP( new protocols::analysis::InterfaceAnalyzerMover )),
		two_ubiquitins_(false),
		extra_bodies_chains_(), //uninitializable
		ubq2_lys_pos_in_complex_(0)
	{
		//set up fullatom scorefunction
		using namespace core::scoring;
		fullatom_scorefunction_ = get_score_function();

		TR << "Using fullatom scorefunction from commandline:\n" << *fullatom_scorefunction_;

		using namespace core::pose::metrics;
		using namespace protocols::toolbox::pose_metric_calculators;
		//magic number: chains 1 and 2; set up interface SASA calculator
		if( !CalculatorFactory::Instance().check_calculator_exists( InterfaceSasaDefinition_ ) ){
			CalculatorFactory::Instance().register_calculator( InterfaceSasaDefinition_, PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator(core::Size(1), core::Size(2)) ));
		}

		IAM_->set_use_centroid_dG(false);

		two_ubiquitins_ = basic::options::option[ basic::options::OptionKeys::chemically_conjugated_docking::two_ubiquitins ].value();
	}

	/// @brief init_on_new_input system allows for initializing these details the first time apply() is called.  the job distributor will reinitialize the whole mover when the input changes (a freshly constructed mover, which will re-run this on first apply().
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
		core::import_pose::pose_from_pdb( E2, basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::E2pdb].value() );
		core::Size const E2length = E2.total_residue();

		core::pose::Pose UBQ;
		core::import_pose::pose_from_pdb( UBQ, basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::UBQpdb].value() );
		core::Size const UBQlength = UBQ.total_residue();
		core::pose::PoseOP UBQ_second;
		if(two_ubiquitins_) UBQ_second = core::pose::PoseOP( new core::pose::Pose(UBQ) );
		if(basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::UBQ2_pdb].user()){
			core::import_pose::pose_from_pdb( *UBQ_second, basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::UBQ2_pdb].value() );
		}

		//determine cysteine target
		runtime_assert(E2.conformation().num_chains() == 1);
		char const E2chain(E2.pdb_info()->chain(1));
		core::Size const E2_cys(E2.pdb_info()->pdb2pose(E2chain, basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::E2_residue].value()));
		runtime_assert(E2.residue_type(E2_cys).aa() == core::chemical::aa_cys);

		//determine c_term target on UBQ
		core::Size const UBQ_term = UBQlength;

		//strip C-term from UBQ - best to do this with a full replace to re-draw the carboxyl oxygen
		//UBQ.dump_pdb("pre-removeUQB.pdb");
		core::chemical::ResidueTypeSetCOP fa_standard(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
		UBQ.conformation().delete_residue_slow( UBQ_term );
		UBQ.append_residue_by_bond( *(core::conformation::ResidueFactory::create_residue(fa_standard->name_map("GLY")) ) );
		UBQ.conformation().insert_ideal_geometry_at_polymer_bond( UBQ_term-1 );
		UBQ.set_omega(UBQ_term-1, 180);
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

		//complex.dump_pdb("isthisthefix.pdb");
		core::Size const complexlength( complex.total_residue());
		complex.conformation().insert_ideal_geometry_at_polymer_bond( complexlength-1 );
		complex.set_omega(complexlength-1, 180);
		complex.conformation().insert_chain_ending(E2length);
		//complex.dump_pdb("initcomplex.pdb");

		//pack atom ID vector
		atomIDs[1] = atom0;
		atomIDs[2] = atom1;
		atomIDs[3] = atom2;
		atomIDs[4] = atom3;
		atomIDs[5] = core::id::AtomID( ubq_rsd_type.atom_index("C" ), complexlength );
		atomIDs[6] = core::id::AtomID( ubq_rsd_type.atom_index("CA" ), complexlength );
		atomIDs[7] = core::id::AtomID( ubq_rsd_type.atom_index("N" ), complexlength );
		atomIDs[8] = core::id::AtomID( ubq_rsd_type.atom_index("C" ), complexlength-1 );

		if(two_ubiquitins_) {
			//now add in the second ubiquitin - link to thioester C=O from UB_lys
			core::Size const ub_lys_pos(basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::UBQ2_lys].value());
			complex.conformation().append_residue_by_jump(UBQ_second->residue(ub_lys_pos), complexlength, "C", "NZ", true);
			//there is a bug lurking here somewhere for C-terminal lysines
			for( core::Size i(ub_lys_pos+1); i <= UBQ_second->total_residue(); ++i) complex.conformation().append_polymer_residue_after_seqpos(UBQ_second->residue(i), complex.total_residue(), false);
			for( core::Size i(ub_lys_pos-1); i >= 1; --i) complex.conformation().prepend_polymer_residue_before_seqpos(UBQ_second->residue(i), complexlength+1, false);

			//check it!
			TR << complex.fold_tree() << std::endl;

			ubq2_lys_pos_in_complex_ = ub_lys_pos+complexlength;

			//continue paking atomIDs vector
			core::chemical::ResidueType const & lys_rsd_type( fa_standard->name_map("LYS") );
			//atomIDs[LYS_2HZ] = core::id::AtomID( lys_rsd_type.atom_index("2HZ"), ubq2_lys_pos_in_complex_);
			//atomIDs[LYS_NZ]  = core::id::AtomID( lys_rsd_type.atom_index("NZ" ), ubq2_lys_pos_in_complex_);
			//atomIDs[LYS_CE]  = core::id::AtomID( lys_rsd_type.atom_index("CE" ), ubq2_lys_pos_in_complex_);
			//atomIDs[LYS_CD]  = core::id::AtomID( lys_rsd_type.atom_index("CD" ), ubq2_lys_pos_in_complex_);

			atomIDs.push_back(core::id::AtomID( lys_rsd_type.atom_index("2HZ"), ubq2_lys_pos_in_complex_));
			atomIDs.push_back(core::id::AtomID( lys_rsd_type.atom_index("NZ" ), ubq2_lys_pos_in_complex_));
			atomIDs.push_back(core::id::AtomID( lys_rsd_type.atom_index("CE" ), ubq2_lys_pos_in_complex_));
			atomIDs.push_back(core::id::AtomID( lys_rsd_type.atom_index("CD" ), ubq2_lys_pos_in_complex_));

			//ok, jump is in place (numbered 1) - now we have to make a statement about where to put it
			//relevant atoms for placing NZ
			core::Vector const & C_xyz(complex.residue(complexlength).atom("C").xyz());
			core::Vector const & O_xyz(complex.residue(complexlength).atom("O").xyz());
			//TR << "C " << C_xyz << " O " << O_xyz << std::endl;
			core::Vector newpos(C_xyz+(3*(C_xyz-O_xyz)/O_xyz.distance(C_xyz)));
			core::Vector oldpos(complex.residue(ubq2_lys_pos_in_complex_).atom("NZ").xyz());

			//core::Vector oldposC(complex.residue(ubq2_lys_pos_in_complex_).atom("CE").xyz());
			//core::Vector newposC(newpos+(oldpos-oldposC));

			//force NZ (alone) to proper position - this breaks the lysine
			complex.set_xyz(core::id::AtomID(lys_rsd_type.atom_index("NZ"), ubq2_lys_pos_in_complex_), newpos);
			//	complex.set_xyz(core::id::AtomID(lys_rsd_type.atom_index("CE"), ubq2_lys_pos_in_complex_), newposC);

			//get a copy of that jump
			core::kinematics::Jump newjump(complex.atom_tree().jump(core::id::AtomID(lys_rsd_type.atom_index("NZ"), ubq2_lys_pos_in_complex_)));

			//reset NZ position to fix the lysine
			complex.set_xyz(core::id::AtomID(lys_rsd_type.atom_index("NZ"), ubq2_lys_pos_in_complex_), oldpos);
			//	complex.set_xyz(core::id::AtomID(lys_rsd_type.atom_index("CE"), ubq2_lys_pos_in_complex_), oldposC);

			//reapply the properly-calculated jump
			complex.conformation().set_jump(1, newjump);
		}

		////////////////////////////extra bodies/////////////////////////////////////////////////
		//The purpose of this code is to allow for static extra things in the system; for its original
		//incarnations, to allow for a RING domain to occlude its site on E2 and a GAP (+MG +GDP) to
		//occlude ras

		//check if extra bodies exist
		if (basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user() == true) {
			extra_bodies_chains_ = apps::public1::scenarios::chemically_conjugated_docking::add_extra_bodies(complex, TR);
		}

		starting_pose_ = complex;
		starting_pose_.dump_pdb("starting_complex.pdb");
		TR << "starting pose finished." << std::endl;

		///////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////Finish starting complex pose/////////////////////////////////////////////
		//////////////////////Start creating move accessory data///////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////

		//setup MoveMaps
		//small/shear behave improperly @ the last residue - psi is considered nonexistent and the wrong phis apply.
		bool const dont_minimize_omega(basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::dont_minimize_omega].value());
		thioester_mm_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
		for( core::Size i(1), ntailres(basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::n_tail_res]); i<ntailres; ++i){ //slightly irregular < comparison because C-terminus is functionally zero-indexed
			thioester_mm_->set_bb((complexlength-i), true);
			if(dont_minimize_omega){
				thioester_mm_->set( core::id::TorsionID(complexlength-i, core::id::BB, core::id::omega_torsion), false);

			}
		}
		//thioester_mm_->set_bb(complexlength, true);
		//thioester_mm_->set(core::id::TorsionID(complexlength, core::id::BB, core::id::phi_torsion), true);
		//thioester_mm_->set(core::id::TorsionID(complexlength, core::id::BB, core::id::psi_torsion), true);
		//thioester_mm_->set_bb(complexlength-1, true);
		//thioester_mm_->set_bb(complexlength-2, true);
		//thioester_mm_->set(complex.atom_tree().torsion_angle_dof_id(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5]), false);

		//setup loop
		std::set< core::Size > loop_posns;
		if ( basic::options::option[ basic::options::OptionKeys::loops::loop_file ].user() ) {
			loop_ = *( protocols::loops::Loops( true ).begin() );
			TR << "loop " <<  loop_ << std::endl;
			//set up interface-plus-neighbors-positions operation
			for (core::Size j(loop_.start()), end(loop_.stop()); j <= end; ++j){
				loop_posns.insert(j);
			}//for each residue in loop
		} //no else needed - default loop is safe enough

		//setup of TaskFactory
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		task_factory_ = core::pack::task::TaskFactoryOP( new TaskFactory );
		task_factory_->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
		if ( basic::options::option[ basic::options::OptionKeys::packing::resfile ].user() ) {
			task_factory_->push_back( TaskOperationCOP( new ReadResfile ) );
		}
		//task_factory_->push_back( new protocols::toolbox::task_operations::RestrictToInterfaceOperation );
		task_factory_->push_back( TaskOperationCOP( new IncludeCurrent ) );
		//prevent repacking at linkage cysteine!
		PreventRepackingOP prevent( new PreventRepacking );
		prevent->include_residue(E2_cys);
		if(two_ubiquitins_) {
			//prevent repacking at the lysine!
			prevent->include_residue(ubq2_lys_pos_in_complex_);
 		}
		task_factory_->push_back(prevent);

		//old way - this is not wrong, but it is inferior to the method below
		using core::pose::metrics::PoseMetricCalculatorOP;
		if(basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::publication].value()){
			using core::pose::metrics::PoseMetricCalculatorOP;
			std::string const interface_calc("UBQE2_InterfaceNeighborDefinitionCalculator");
			std::string const neighborhood_calc("UBQE2_NeighborhoodByDistanceCalculator");
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( interface_calc, PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( core::Size(1), core::Size(2)) ) );
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( neighborhood_calc, PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( loop_posns ) ) );

			//this is the constructor parameter for the calculator - pairs of calculators and calculations to perform
			utility::vector1< std::pair< std::string, std::string> > calcs_and_calcns;
			calcs_and_calcns.push_back(std::make_pair(interface_calc, "interface_residues"));
			calcs_and_calcns.push_back(std::make_pair(neighborhood_calc, "neighbors"));

			using protocols::toolbox::task_operations::RestrictByCalculatorsOperation;
			task_factory_->push_back(TaskOperationCOP( new RestrictByCalculatorsOperation( calcs_and_calcns ) ));

		} else {

			TR << "using new way" << std::endl;
			//new way
			//partially stolen from FloppyTail - I need to go back and extract this code unit
			utility::vector1< std::set < core::Size > > regions; //a set of regions to turn into groups for comparison
			std::set < core::Size > const empty; //easier to add empty sets to the vector than construct-then-add

			//E2
			core::Size const E2_end(complex.conformation().chain_end(1));
			regions.push_back(empty); //insert a new set to work with
			core::Size const E2_index(1);
			for(core::Size i(1); i<=E2_end; ++i) {
				regions[E2_index].insert(i);
			}

			//ubiquitin tail
			//complexlength = end of ubiquitin in ubq+e2 pose
			core::Size const tail_begin(complexlength-basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::n_tail_res]+1); //odd construction accounts for functional zero-indexing of the tail
			regions.push_back(empty); //insert a new set to work with
			core::Size const tail_index(2);
			for(core::Size i(complexlength); i>=tail_begin; --i) {
				regions[tail_index].insert(i);
			}

			//ubiquitin core+tail
			//including the tail in both groups ensures it will always repack
			regions.push_back(empty); //insert a new set to work with
			core::Size const ubq_index(3);
			for(core::Size i(E2_end+1); i<=complexlength; ++i) {
				regions[ubq_index].insert(i);
			}

			if(two_ubiquitins_) {
				//second ubiquitin
				TR << "Assuming second ubiquitin is the third chain" << std::endl;
				runtime_assert(complex.conformation().num_chains() >= 3);
				regions.push_back(empty); //insert a new set to work with
				core::Size const ubq2_index(4);
				core::Size const ubq2_end(complex.conformation().chain_end(3));
				for(core::Size i(complexlength+1); i<=ubq2_end; ++i) {
					regions[ubq2_index].insert(i);
				}
			}

			//this will double-count loop residues with E2 - but that's fine, it just ensures they pack no matter what
			if( !loop_posns.empty() ){
				regions.push_back(loop_posns);
			}

			//if extra bodies exist, we are adding them to the first chain set
			if (basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user() == true) {
				apps::public1::scenarios::chemically_conjugated_docking::pack_extra_bodies(extra_bodies_chains_, complex, regions[E2_index], TR);
			}

			//make all pairs of groups (without replacement
			//if you have 1, 2, 3, 4; make 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
			core::Size const num_regions(regions.size());
			utility::vector1< std::pair< std::set<core::Size>, std::set<core::Size> > > vector_of_pairs;
			for(core::Size first_group(1); first_group < num_regions; ++first_group) {
				for(core::Size second_group(first_group+1); second_group <= num_regions; ++second_group){
					vector_of_pairs.push_back(std::make_pair(regions[first_group], regions[second_group]));
				}
			}

			//check contents of vector_of_pairs
			core::Size const num_pairs(vector_of_pairs.size());
			for(core::Size i(1); i<=num_pairs; ++i){
				core::Size const
					onestart(*(vector_of_pairs[i].first.begin())),
					onestop(*(vector_of_pairs[i].first.rbegin())),
					twostart(*(vector_of_pairs[i].second.begin())),
					twostop(*(vector_of_pairs[i].second.rbegin()));

				TR << "IGNC will compare group " << onestart << "-" << onestop << " with " << twostart << "-" << twostop << std::endl;

				// if (basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user() == true) {
				// 	TR.Error << "upcoming debug me non-contiguous set errors are not really errors if using extra bodies" << std::endl;
				// }

				// core::Size guess(onestart);
				// for(std::set<core::Size>::const_iterator iter(vector_of_pairs[i].first.begin()), end(vector_of_pairs[i].first.end()); iter != end; ++iter) {
				// 	if(guess++ != *iter) TR.Error << "non-contiguous set, debug me!" << std::endl;
				// 	TR << *iter << std::endl;
				// }
				// guess = twostart;
				// for(std::set<core::Size>::const_iterator iter(vector_of_pairs[i].second.begin()), end(vector_of_pairs[i].second.end()); iter != end; ++iter) {
				// 	if(guess++ != *iter) TR.Error << "non-contiguous set, debug me!" << std::endl;
				// 	TR << *iter << std::endl;
				// }

			}

			if (basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user() == true) {
				TR << "Those group labels do not take the piling of extra bodies into the first nonmoving group into account" << std::endl;
			}

			//check if calculator exists; create if not
			std::string const calc("IGNC_UBQ_E2_thioester");
			if(core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists(calc)){
				core::pose::metrics::CalculatorFactory::Instance().remove_calculator(calc);
				TR.Error << "removed a PoseMetricCalculator " << calc << ", track down why" << std::endl;
			}
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( calc, PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::InterGroupNeighborsCalculator(vector_of_pairs) ) );

			//now that calculator exists, add the sucker to the TaskFactory via RestrictByCalculatorsOperation
			utility::vector1< std::pair< std::string, std::string> > calculators_used;
			std::pair< std::string, std::string> IGNC_cmd( calc, "neighbors" );
			calculators_used.push_back( IGNC_cmd );
			task_factory_->push_back( TaskOperationCOP( new protocols::toolbox::task_operations::RestrictByCalculatorsOperation( calculators_used ) ) );

		}

		//calculator for number of neighbors for I44
		if ( basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::publication].value() ) {
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( "I44neighbors", PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::NeighborsByDistanceCalculator( 198 ) ) );
		}

		//add constraints; protected internally if no constraints
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( starting_pose_ );
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *fullatom_scorefunction_ );
		core::scoring::constraints::ConstraintSetCOP cst_set(starting_pose_.constraint_set());
		if(cst_set && false){
			cst_set->show_definition(TR, starting_pose_);
		}
	}

	virtual ~UBQ_E2Mover(){};

	virtual
	void
	apply( core::pose::Pose & pose ){
		if( !init_for_input_yet_ ) init_on_new_input();

		pose = starting_pose_;

		//TR << "foldtree, movemap: " << std::endl;
		//core::kinematics::simple_visualize_fold_tree_and_movemap( pose.fold_tree(), *movemap_, TR);

		TR << "foldtree, " << pose.fold_tree() << std::endl;

		/////////////////fullatom Monte Carlo//////////////////////////////////////////////////////////
		//make the monte carlo object
		using protocols::moves::MonteCarlo;
		using protocols::moves::MonteCarloOP;
		using basic::options::option;
		using namespace basic::options::OptionKeys::AnchoredDesign;
		MonteCarloOP mc( new MonteCarlo( pose, *fullatom_scorefunction_, option[ refine_temp ].value() ) );

		//////////////////////////Small/ShearMovers////////////////////////////////////////////////////////
		protocols::simple_moves::BackboneMoverOP small_mover( new protocols::simple_moves::SmallMover(thioester_mm_, 0.8, 1) );
		small_mover->angle_max( 'H', 4.0 );
		small_mover->angle_max( 'E', 4.0 );
		small_mover->angle_max( 'L', 4.0 );

		protocols::simple_moves::BackboneMoverOP shear_mover( new protocols::simple_moves::ShearMover(thioester_mm_, 0.8, 1) );
		shear_mover->angle_max( 'H', 4.0 );
		shear_mover->angle_max( 'E', 4.0 );
		shear_mover->angle_max( 'L', 4.0 );

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_chi1( new protocols::simple_moves::TorsionDOFMover );
		DOF_mover_chi1->set_DOF(atomIDs[1], atomIDs[2], atomIDs[3], atomIDs[4]);
		DOF_mover_chi1->check_mmt(true);
		DOF_mover_chi1->temp(0.4);
		DOF_mover_chi1->set_angle_range(-180, 180);
		DOF_mover_chi1->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_chi2( new protocols::simple_moves::TorsionDOFMover );
		DOF_mover_chi2->set_DOF(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5]);
		DOF_mover_chi2->check_mmt(true);
		DOF_mover_chi2->temp(0.4);
		DOF_mover_chi2->set_angle_range(-180, 180);
		DOF_mover_chi2->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_thioester( new protocols::simple_moves::TorsionDOFMover );
		DOF_mover_thioester->set_DOF(atomIDs[3], atomIDs[4], atomIDs[5], atomIDs[6]);
		DOF_mover_thioester->check_mmt(true);
		DOF_mover_thioester->temp(0.4);
		DOF_mover_thioester->set_angle_range(-180, 180);
		DOF_mover_thioester->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_psi( new protocols::simple_moves::TorsionDOFMover );
		DOF_mover_psi->set_DOF(atomIDs[4], atomIDs[5], atomIDs[6], atomIDs[7]);
		DOF_mover_psi->check_mmt(true);
		DOF_mover_psi->temp(0.4);
		DOF_mover_psi->set_angle_range(-180, 180);
		DOF_mover_psi->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_phi( new protocols::simple_moves::TorsionDOFMover );
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
			KinematicWrapperOP kin_wrapper( new KinematicWrapper(kin_mover, loop_) );

			backbone_mover->add_mover(kin_wrapper, 5);

		}

		if(two_ubiquitins_){
			//////////////////////////////RotateJumpAxisMover for second ubiquitin//////////////////////
			protocols::rigid::RotateJumpAxisMoverOP RJAmover( new protocols::rigid::RotateJumpAxisMover(1) );
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
			SCmover->set_change_chi_without_replacing_residue(true);

			backbone_mover->add_mover(SCmover, 1);

			//////////////////////////////added "chi" for lysine//////////////////////////////////////
			protocols::simple_moves::TorsionDOFMoverOP DOF_mover_lys( new protocols::simple_moves::TorsionDOFMover );
			DOF_mover_lys->set_DOF(atomIDs[LYS_2HZ], atomIDs[LYS_NZ], atomIDs[LYS_CE], atomIDs[LYS_CD]);
			DOF_mover_lys->check_mmt(false);
			DOF_mover_lys->temp(10);
			DOF_mover_lys->set_angle_range(-180, 180);
			DOF_mover_lys->tries(1000);

			backbone_mover->add_mover(DOF_mover_lys, 1);
		}

		/////////////////////////minimize backbone DOFs//////////////////////////////////////////////
		using protocols::simple_moves::MinMoverOP;
		using protocols::simple_moves::MinMover;
		protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover(
																				thioester_mm_,
																				fullatom_scorefunction_,
																				basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
																				0.01,
																				true /*use_nblist*/ ) );

		/////////////////////////////////rotamer trials mover///////////////////////////////////////////
		using protocols::simple_moves::RotamerTrialsMoverOP;
		using protocols::simple_moves::EnergyCutRotamerTrialsMover;
		protocols::simple_moves::RotamerTrialsMoverOP rt_mover( new protocols::simple_moves::EnergyCutRotamerTrialsMover(
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
																																											20.0) );

		///////////////////////////////repack///////////////////////////////////////////////
		protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover );
		pack_mover->task_factory( task_factory_ );
		pack_mover->score_function( fullatom_scorefunction_ );

		protocols::simple_moves::MinMoverOP min_mover_pack( new protocols::simple_moves::MinMover(
																						 thioester_mm_,
																						 fullatom_scorefunction_,
																						 basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
																						 0.01,
																						 true /*use_nblist*/ ) );

		using protocols::simple_moves::TaskAwareMinMoverOP;
		using protocols::simple_moves::TaskAwareMinMover;
		protocols::simple_moves::TaskAwareMinMoverOP TAmin_mover( new protocols::simple_moves::TaskAwareMinMover(min_mover_pack, task_factory_) );

		/////////////////////////////////////////refine loop///////////////////////////////////////////

		core::Size const refine_applies = option[ refine_cycles ].value(); //default 5
		core::Size const repack_cycles = option[ refine_repack_cycles ].value();
		//core::Size const min_cycles = repack_cycles/2;
		TR << "   Current     Low    total cycles =" << refine_applies << std::endl;
		for ( core::Size i(1); i <= refine_applies; ++i ) {
			//pdb_out1.apply(pose);
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
		if( score > basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::scorefilter].value() ){
			set_last_move_status(protocols::moves::FAIL_RETRY);
			TR << "total score filter failed; score " << score << std::endl;
			return;
		}

		//these interface analyses are less interpretable in the three-body case
		if(!two_ubiquitins_ && !basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user()) {

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
			if(mv_delta_sasa.value() < basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::SASAfilter].value()){
				set_last_move_status(protocols::moves::FAIL_RETRY);
				TR << "interface SASA filter failed; SASA " << mv_delta_sasa.value() << std::endl;
				return;
			}

			//passed filters; run IAM
			IAM_->apply(copy);

		}

		//print mobile region fine-grained data
		protocols::jd2::JobOP job_me(protocols::jd2::JobDistributor::get_instance()->current_job());
		using numeric::conversions::degrees;
		job_me->add_string_real_pair("cysteine_chi1_C-CA-CB-SG", degrees(pose.atom_tree().torsion_angle(atomIDs[1], atomIDs[2], atomIDs[3], atomIDs[4])));
		job_me->add_string_real_pair("cysteine_chi2_CA-CB-SG-C", degrees(pose.atom_tree().torsion_angle(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5])));
		job_me->add_string_real_pair("thioester_CB-SG-C-CA", degrees(pose.atom_tree().torsion_angle(atomIDs[3], atomIDs[4], atomIDs[5], atomIDs[6])));
		job_me->add_string_real_pair("Cterm_psi_SG-C-CA-N", degrees(pose.atom_tree().torsion_angle(atomIDs[4], atomIDs[5], atomIDs[6], atomIDs[7])));
		job_me->add_string_real_pair("Cterm_phi_C-CA-N-C", degrees(pose.atom_tree().torsion_angle(atomIDs[5], atomIDs[6], atomIDs[7], atomIDs[8])));

		if ( basic::options::option[ basic::options::OptionKeys::chemically_conjugated_docking::publication ].value()) {
			//I44 neighbors
			basic::MetricValue< core::Size > I44numn;
			pose.metric("I44neighbors", "num_neighbors", I44numn);
			job_me->add_string_real_pair("I44neighbors", I44numn.value());
		}

		set_last_move_status(protocols::moves::MS_SUCCESS);
		return;
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return protocols::moves::MoverOP( new UBQ_E2Mover );
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

	/// @brief vector contains atomIDs for thioester bond and atoms before/after bond to determine various torsions
	utility::vector1< core::id::AtomID > atomIDs;

	core::pose::Pose starting_pose_; //maintained from run to run

	std::string const InterfaceSasaDefinition_; //calculator name

	protocols::analysis::InterfaceAnalyzerMoverOP IAM_;

	/// @brief used for two-ubiquitins mode
	bool two_ubiquitins_;

	/// @brief used to track which chains are "extra" nonmoving bodies in extra bodies mode
	utility::vector1< core::Size > extra_bodies_chains_;

	/// @brief the lysine attachment position for the second ubiquitin (the second moving chain), in the whole complex numbering
	core::Size ubq2_lys_pos_in_complex_;

};

typedef utility::pointer::shared_ptr< UBQ_E2Mover > UBQ_E2MoverOP;

int main( int argc, char* argv[] )
{
try {
	//initialize options
	devel::init(argc, argv);
	basic::prof_reset();

	if(basic::options::option[ basic::options::OptionKeys::in::file::s ].user()
		|| basic::options::option[ basic::options::OptionKeys::in::file::l ].user()
		|| basic::options::option[ basic::options::OptionKeys::in::file::silent ].user())
		utility_exit_with_message("do not use an input PDB with this protocol (program uses internally); use -UBQpdb and -E2pdb instead");

	protocols::jd2::JobDistributor::get_instance()->go(protocols::moves::MoverOP( new UBQ_E2Mover ));

	basic::prof_show();
	TR << "************************d**o**n**e**************************************" << std::endl;
} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
	return 0;
}
