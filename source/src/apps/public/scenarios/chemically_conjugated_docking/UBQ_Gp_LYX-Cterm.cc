// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/apps/public/scenarios/chemically_conjugated_docking/UBQ_Gp_LYX-Cterm.cc
/// @brief  this application is a one-shot for modeling a ubiquitinated G-protein; this version uses a native lysine to ubiquitin linkage via LYX
/// @author Steven Lewis

// Unit Headers
#include <apps/public/scenarios/chemically_conjugated_docking/Gp_quantification_metrics.hh>
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
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictByCalculatorsOperation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/kinematics/util.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

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
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>

#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/InterGroupNeighborsCalculator.hh>
// AUTO-REMOVED #include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>

// Numeric Headers
#include <numeric/conversions.hh>
// AUTO-REMOVED #include <numeric/xyz.io.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <basic/prof.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/kinematics/AtomTree.hh>

//tracers
using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR("apps.public.scenarios.chemically_conjugated_docking.UBQ_Gp_LYX-Cterm");

class UBQ_GTPaseMover : public protocols::moves::Mover {
public:
	UBQ_GTPaseMover()
	: init_for_input_yet_(false),
		fullatom_scorefunction_(NULL),
		task_factory_(NULL),
		amide_mm_(NULL),
		loop_(), //we want default ctor
		atomIDs(8, core::id::BOGUS_ATOM_ID ),
		InterfaceSasaDefinition_("InterfaceSasaDefinition_" + 1),
		IAM_(new protocols::analysis::InterfaceAnalyzerMover),
		GTPase_lys_(0)
	{
		//set up fullatom scorefunction
		using namespace core::scoring;
		fullatom_scorefunction_ = get_score_function();
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *fullatom_scorefunction_ ); //protected if(option) internally

		TR << "Using fullatom scorefunction from commandline:\n" << *fullatom_scorefunction_;

		using namespace core::pose::metrics;
		using namespace protocols::toolbox::pose_metric_calculators;
		//magic number: chains 1 and 2; set up interface SASA calculator
		if( !CalculatorFactory::Instance().check_calculator_exists( InterfaceSasaDefinition_ ) ){
			CalculatorFactory::Instance().register_calculator( InterfaceSasaDefinition_, new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator(core::Size(1), core::Size(2)));
		}

		IAM_->set_use_centroid_dG(false);
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
		core::pose::Pose GTPase;
		core::import_pose::pose_from_pdb( GTPase, basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::GTPasepdb].value() );
		core::Size const GTPaselength = GTPase.total_residue();

		core::pose::Pose UBQ;
		core::import_pose::pose_from_pdb( UBQ, basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::UBQpdb].value() );
		core::Size const UBQlength = UBQ.total_residue();

		//determine cysteine target
		runtime_assert(GTPase.conformation().num_chains() == 1);
		char const GTPasechain(GTPase.pdb_info()->chain(1));
		GTPase_lys_ = GTPase.pdb_info()->pdb2pose(GTPasechain, basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::GTPase_residue].value());
		//runtime_assert(GTPase.residue_type(GTPase_lys_).aa() == core::chemical::aa_lys);

		//determine c_term target on UBQ
		core::Size const UBQ_term = UBQlength;

		//strip C-term from UBQ - best to do this with a full replace to re-draw the carboxyl oxygen
		//UBQ.dump_pdb("pre-removeUQB.pdb");
		core::chemical::ResidueTypeSetCAP fa_standard(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
		UBQ.conformation().delete_residue_slow( UBQ_term );
		UBQ.append_residue_by_bond( *(core::conformation::ResidueFactory::create_residue(fa_standard->name_map("GLY")) ) );
		//UBQ.dump_pdb("post-removeUQB.pdb");

		//replace lysine
		core::chemical::ResidueType const & lyx_rsd_type( fa_standard->name_map("LYX") );
		//GTPase.dump_pdb("prereplace_GTPase.pdb");
		GTPase.replace_residue( GTPase_lys_, core::conformation::Residue(lyx_rsd_type, true), true);
		//GTPase.dump_pdb("postreplace_GTPase.pdb");

		// check safety of connections (from phil)
		core::chemical::ResidueType const & ubq_rsd_type( UBQ.residue_type( UBQ_term ) );
		core::Size const lyx_connid( 3 );
		core::Size const ubq_connid( 2 );

		runtime_assert( lyx_rsd_type.n_residue_connections() == lyx_connid &&
			lyx_rsd_type.lower_connect_id() != lyx_connid &&
			lyx_rsd_type.upper_connect_id() != lyx_connid );

		runtime_assert( ubq_rsd_type.n_residue_connections() == ubq_connid &&
			ubq_rsd_type.lower_connect_id() != ubq_connid);

		//GTPase.conformation().show_residue_connections();
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
		core::pose::Pose complex(GTPase);
		complex.append_residue_by_bond( UBQ.residue( UBQ_term ), true, ubq_connid, GTPase_lys_, lyx_connid );
		//complex.dump_pdb("just1_complex.pdb");

		//not that this does anything
		complex.conformation().insert_ideal_geometry_at_residue_connection( GTPase_lys_, lyx_connid );

		core::Size const ubq_pos( complex.total_residue() );
		core::id::AtomID const atom0( lyx_rsd_type.atom_index( "CG" ), GTPase_lys_ );
		core::id::AtomID const atom1( lyx_rsd_type.atom_index( "CD" ), GTPase_lys_ );
		core::id::AtomID const atom2( lyx_rsd_type.atom_index( "CE" ), GTPase_lys_ );
		core::id::AtomID const atom3( lyx_rsd_type.atom_index( "NZ" ), GTPase_lys_ );
		core::id::AtomID const atom4( ubq_rsd_type.atom_index( "C"  ), ubq_pos );
		core::id::AtomID const atom5( ubq_rsd_type.atom_index( "CA" ), ubq_pos );
		core::id::AtomID const atom6( ubq_rsd_type.atom_index( "N"  ), ubq_pos );

		//starting values derived from the peptide bond and a straight-out lysine
		for( core::Size chi(1); chi<=4; ++chi) complex.set_chi(chi, GTPase_lys_, 180);
		//complex.conformation().set_torsion_angle( atom0, atom1, atom2, atom3, numeric::conversions::radians(106.5) );
		//complex.conformation().set_torsion_angle( atom1, atom2, atom3, atom4, numeric::conversions::radians(-60.0) );
		complex.conformation().set_torsion_angle( atom2, atom3, atom4, atom5, numeric::conversions::radians(180.0) );
		complex.conformation().set_torsion_angle( atom3, atom4, atom5, atom6, numeric::conversions::radians(135.0) );
		//complex.dump_pdb("just1_complex2.pdb");

		//now add the rest of ubiquitin
		for( core::Size i=UBQ_term-1; i>= 1; --i ) {
			complex.prepend_polymer_residue_before_seqpos( UBQ.residue(i), GTPaselength+1, false );
		}

		core::Size const complexlength( complex.total_residue());
		complex.conformation().insert_ideal_geometry_at_polymer_bond( complexlength-1 );
		complex.conformation().insert_chain_ending(GTPaselength);
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
		amide_mm_ = new core::kinematics::MoveMap;
		for( core::Size i(1), ntailres(basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::n_tail_res]); i<ntailres; ++i){ //slightly irregular < comparison because C-terminus is functionally zero-indexed
			amide_mm_->set_bb((complexlength-i), true);
		}
		//amide_mm_->set_bb(complexlength, true);
		//amide_mm_->set(core::id::TorsionID(complexlength, core::id::BB, core::id::phi_torsion), true);
		//amide_mm_->set(core::id::TorsionID(complexlength, core::id::BB, core::id::psi_torsion), true);
		//amide_mm_->set_bb(complexlength-1, true);
		//amide_mm_->set_bb(complexlength-2, true);
		//amide_mm_->set(complex.atom_tree().torsion_angle_dof_id(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5]), false);

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
		task_factory_ = new TaskFactory;
		task_factory_->push_back( new InitializeFromCommandline );
		if ( basic::options::option[ basic::options::OptionKeys::packing::resfile ].user() ) {
			task_factory_->push_back( new ReadResfile );
		}
		//task_factory_->push_back( new protocols::toolbox::task_operations::RestrictToInterfaceOperation );
		task_factory_->push_back( new IncludeCurrent );
		//prevent repacking at linkage lysine!
		PreventRepackingOP prevent(new PreventRepacking);
		prevent->include_residue(GTPase_lys_);
		task_factory_->push_back(prevent);

		if(false){
			std::string const interface_calc("UBQGTPase_InterfaceNeighborDefinitionCalculator");
			std::string const neighborhood_calc("UBQGTPase_NeighborhoodByDistanceCalculator");
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( interface_calc, new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( core::Size(1), core::Size(2)) );
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( neighborhood_calc, new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( loop_posns ) );

			//this is the constructor parameter for the calculator - pairs of calculators and calculations to perform
			utility::vector1< std::pair< std::string, std::string> > calcs_and_calcns;
			calcs_and_calcns.push_back(std::make_pair(interface_calc, "interface_residues"));
			calcs_and_calcns.push_back(std::make_pair(neighborhood_calc, "neighbors"));

			using protocols::toolbox::task_operations::RestrictByCalculatorsOperation;
			task_factory_->push_back(new RestrictByCalculatorsOperation( calcs_and_calcns ));

		} else {
			//functions, etc here use UBQ/E2 nomenclature until I can extract it out
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
			std::string const calc("IGNC_UBQ_Gp_LYX-Cterm");
			if(core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists(calc)){
				core::pose::metrics::CalculatorFactory::Instance().remove_calculator(calc);
				TR.Error << "removed a PoseMetricCalculator " << calc << ", track down why" << std::endl;
			}
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( calc, new protocols::toolbox::pose_metric_calculators::InterGroupNeighborsCalculator(vector_of_pairs) );

			//now that calculator exists, add the sucker to the TaskFactory via RestrictByCalculatorsOperation
			utility::vector1< std::pair< std::string, std::string> > calculators_used;
			std::pair< std::string, std::string> IGNC_cmd( calc, "neighbors" );
			calculators_used.push_back( IGNC_cmd );
			task_factory_->push_back( new protocols::toolbox::task_operations::RestrictByCalculatorsOperation( calculators_used ) );

		}

		//create constraints
		core::scoring::constraints::add_constraints_from_cmdline_to_pose( starting_pose_ ); //protected internally if no constraints

	}
	virtual ~UBQ_GTPaseMover(){};

	virtual
	void
	apply( core::pose::Pose & pose ){
		if( !init_for_input_yet_ ) init_on_new_input();

		pose = starting_pose_;

		//TR << "foldtree, movemap: " << std::endl;
		//core::kinematics::simple_visualize_fold_tree_and_movemap( pose.fold_tree(), *movemap_, TR);

		TR << "foldtree, " << pose.fold_tree() << std::flush;

		/////////////////fullatom Monte Carlo//////////////////////////////////////////////////////////
		//make the monte carlo object
		using protocols::moves::MonteCarlo;
		using protocols::moves::MonteCarloOP;
		using basic::options::option;
		using namespace basic::options::OptionKeys::AnchoredDesign;
		MonteCarloOP mc( new MonteCarlo( pose, *fullatom_scorefunction_, option[ refine_temp ].value() ) );

		//////////////////////////Small/ShearMovers////////////////////////////////////////////////////////
		protocols::simple_moves::BackboneMoverOP small_mover = new protocols::simple_moves::SmallMover(amide_mm_, 0.8, 1);
		small_mover->angle_max( 'H', 4.0 );
		small_mover->angle_max( 'E', 4.0 );
		small_mover->angle_max( 'L', 4.0 );

		protocols::simple_moves::BackboneMoverOP shear_mover = new protocols::simple_moves::ShearMover(amide_mm_, 0.8, 1);
		shear_mover->angle_max( 'H', 4.0 );
		shear_mover->angle_max( 'E', 4.0 );
		shear_mover->angle_max( 'L', 4.0 );

		// These two DOFs should map to the lysine sidechain; SidechainMover should take care of them, I think.
		// protocols::simple_moves::TorsionDOFMoverOP DOF_mover_chi1(new protocols::simple_moves::TorsionDOFMover);
		// DOF_mover_chi1->set_DOF(atomIDs[1], atomIDs[2], atomIDs[3], atomIDs[4]);
		// DOF_mover_chi1->check_mmt(true);
		// DOF_mover_chi1->temp(0.4);
		// DOF_mover_chi1->set_angle_range(-180, 180);
		// DOF_mover_chi1->tries(1000);

		// protocols::simple_moves::TorsionDOFMoverOP DOF_mover_chi2(new protocols::simple_moves::TorsionDOFMover);
		// DOF_mover_chi2->set_DOF(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5]);
		// DOF_mover_chi2->check_mmt(true);
		// DOF_mover_chi2->temp(0.4);
		// DOF_mover_chi2->set_angle_range(-180, 180);
		// DOF_mover_chi2->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_isopeptide(new protocols::simple_moves::TorsionDOFMover);
		DOF_mover_isopeptide->set_DOF(atomIDs[3], atomIDs[4], atomIDs[5], atomIDs[6]);
		DOF_mover_isopeptide->check_mmt(true);
		DOF_mover_isopeptide->temp(0.4);
		DOF_mover_isopeptide->set_angle_range(-180, 180);
		DOF_mover_isopeptide->tries(1000);

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

		//Also add a SidechainMover for LYX (I hope...)
		//set up "pack only the moving conjugate" packer task
		utility::vector1< bool > repack_residues(pose.total_residue(), false); //this could be member data
		repack_residues[GTPase_lys_] = true;
		core::pack::task::PackerTaskOP SC_task(core::pack::task::TaskFactory::create_packer_task(pose) );
		SC_task->restrict_to_residues(repack_residues);
		SC_task->restrict_to_repacking(); //SCmover will design, oops
		//and the mover
		protocols::simple_moves::sidechain_moves::SidechainMoverOP SC_mover(new protocols::simple_moves::sidechain_moves::SidechainMover() );
		SC_mover->set_change_chi_without_replacing_residue(true);
		SC_mover->set_task(SC_task);

		protocols::moves::RandomMoverOP backbone_mover( new protocols::moves::RandomMover() );
		backbone_mover->add_mover(small_mover, 2.0);
		backbone_mover->add_mover(shear_mover, 1.0);
		// 		backbone_mover->add_mover(DOF_mover_chi1, 0.75); //SC mover will handle this DOF
		// 		backbone_mover->add_mover(DOF_mover_chi2, 0.75); //SC mover will handle this DOF
 		backbone_mover->add_mover(DOF_mover_isopeptide, 0.75);
 		backbone_mover->add_mover(DOF_mover_psi, 0.75);
		backbone_mover->add_mover(DOF_mover_phi, 0.75);
		backbone_mover->add_mover(SC_mover, 1.0);

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

		/////////////////////////minimize backbone DOFs//////////////////////////////////////////////
		using protocols::simple_moves::MinMoverOP;
		using protocols::simple_moves::MinMover;
		protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover(
			amide_mm_,
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
			amide_mm_,
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
		if(!basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user()) {

			//filter on interface SASA - requires some hacking to break up isopeptide
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
		job_me->add_string_real_pair("lysine_chi3_CG-CD-CE-NZ", degrees(pose.atom_tree().torsion_angle(atomIDs[1], atomIDs[2], atomIDs[3], atomIDs[4])));
		job_me->add_string_real_pair("lysine_chi4_CD-CE-NZ-C", degrees(pose.atom_tree().torsion_angle(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5])));
		job_me->add_string_real_pair("amide_CE-NZ-C-CA", degrees(pose.atom_tree().torsion_angle(atomIDs[3], atomIDs[4], atomIDs[5], atomIDs[6])));
		job_me->add_string_real_pair("glycine_psi_NZ-C-CA-N", degrees(pose.atom_tree().torsion_angle(atomIDs[4], atomIDs[5], atomIDs[6], atomIDs[7])));
		job_me->add_string_real_pair("glycine_phi_C-CA-N-C", degrees(pose.atom_tree().torsion_angle(atomIDs[5], atomIDs[6], atomIDs[7], atomIDs[8])));

		set_last_move_status(protocols::moves::MS_SUCCESS);
		if(basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::publication].value()){
			apps::public1::scenarios::chemically_conjugated_docking::create_extra_output(pose, TR, !basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::pdz].value(), GTPase_lys_);
		}
		return;
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new UBQ_GTPaseMover;
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }

	virtual
	std::string
	get_name() const { return "UBQ_GTPaseMover"; }

private:
	bool init_for_input_yet_;

	core::scoring::ScoreFunctionOP fullatom_scorefunction_;
	core::pack::task::TaskFactoryOP task_factory_;
	core::kinematics::MoveMapOP amide_mm_;
// 	core::kinematics::MoveMapOP loop_mm_;
// 	core::kinematics::MoveMapOP all_mm_;

	protocols::loops::Loop loop_;

	///@brief vector contains atomIDs for isopeptide bond and atoms before/after bond to determine various torsions
	utility::vector1< core::id::AtomID > atomIDs;

	core::pose::Pose starting_pose_; //maintained from run to run

	std::string const InterfaceSasaDefinition_; //calculator name

	protocols::analysis::InterfaceAnalyzerMoverOP IAM_;

	core::Size GTPase_lys_; //converted to member data for sharing between setup and apply

	///@brief used to track which chains are "extra" nonmoving bodies in extra bodies mode
	utility::vector1< core::Size > extra_bodies_chains_;
};

typedef utility::pointer::owning_ptr< UBQ_GTPaseMover > UBQ_GTPaseMoverOP;

int main( int argc, char* argv[] )
{
try {
	//initialize options
	devel::init(argc, argv);
	basic::prof_reset();

	if(basic::options::option[ basic::options::OptionKeys::in::file::s ].user()
		|| basic::options::option[ basic::options::OptionKeys::in::file::l ].user()
		|| basic::options::option[ basic::options::OptionKeys::in::file::silent ].user())
		utility_exit_with_message("do not use an input PDB with this protocol (program uses internally); use -UBQpdb and -GTPasepdb instead");

	protocols::jd2::JobDistributor::get_instance()->go(new UBQ_GTPaseMover);

	basic::prof_show();
	TR << "************************d**o**n**e**************************************" << std::endl;
} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
	return 0;
}
