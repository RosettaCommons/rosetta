// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/public/scenarios/chemically_conjugated_docking/UBQ_Gp_CYD-CYD.cc
/// @brief  this application is a one-shot for modeling a ubiquitinated G-protein; this version uses a disulfide linkage to match an experimentally-induced disulfide
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
#include <utility/excn/Exceptions.hh>
// C++ headers
#include <string>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/AnchoredDesign.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <utility/vector0.hh>

//Auto Headers
#include <core/kinematics/AtomTree.hh>

//tracers
using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR("apps.public.scenarios.chemically_conjugated_docking.UBQ_Gp_CYD-CYD");

class UBQ_GTPase_disulfide_Mover : public protocols::moves::Mover {
public:
	UBQ_GTPase_disulfide_Mover()
	: init_for_input_yet_(false),
		fullatom_scorefunction_(/* NULL */),
		task_factory_(/* NULL */),
		disulfide_mm_(/* NULL */),
		loop_(), //we want default ctor
		atomIDs(8, core::id::BOGUS_ATOM_ID ),
		InterfaceSasaDefinition_("InterfaceSasaDefinition_" + 1),
		IAM_(protocols::analysis::InterfaceAnalyzerMoverOP( new protocols::analysis::InterfaceAnalyzerMover ))
	{
		//set up fullatom scorefunction
		using namespace core::scoring;
		fullatom_scorefunction_ = get_score_function();
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *fullatom_scorefunction_ ); //protected if(option) internally

		TR << "Using fullatom scorefunction from commandline:\n" << *fullatom_scorefunction_;

		using namespace core::pose::metrics;
		using namespace protocols::toolbox::pose_metric_calculators;
		//magic number: chains 1 and 2; set up interface SASA calculator
		if ( !CalculatorFactory::Instance().check_calculator_exists( InterfaceSasaDefinition_ ) ) {
			CalculatorFactory::Instance().register_calculator( InterfaceSasaDefinition_, PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator(core::Size(1), core::Size(2)) ));
		}

		IAM_->set_use_centroid_dG(false);
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
		core::pose::Pose GTPase;
		core::import_pose::pose_from_file( GTPase, basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::GTPasepdb].value() , core::import_pose::PDB_file);
		core::Size const GTPaselength = GTPase.total_residue();

		core::pose::Pose UBQ;
		core::import_pose::pose_from_file( UBQ, basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::UBQpdb].value() , core::import_pose::PDB_file);
		core::Size const UBQlength = UBQ.total_residue();

		//determine cysteine target
		runtime_assert(GTPase.conformation().num_chains() == 1);
		char const GTPasechain(GTPase.pdb_info()->chain(1));
		GTPase_cyd_ = GTPase.pdb_info()->pdb2pose(GTPasechain, basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::GTPase_residue].value());
		//runtime_assert(GTPase.residue_type(GTPase_cyd_).aa() == core::chemical::aa_cys);

		//determine c_term target on UBQ
		core::Size const UBQ_term = UBQlength;

		//replace cysteine
		core::chemical::ResidueTypeSetCOP fa_standard(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
		core::chemical::ResidueType const & cyd_rsd_type( fa_standard->name_map("CYS:disulfide") );
		//GTPase.dump_pdb("prereplace_GTPase.pdb");
		GTPase.replace_residue( GTPase_cyd_, core::conformation::Residue(cyd_rsd_type, true), true);
		//GTPase.dump_pdb("postreplace_GTPase.pdb");

		//replace with CYD on ubiquitin too
		core::chemical::ResidueType const & cyd_rsd_term_type(fa_standard->name_map("CYS:CtermProteinFull:disulfide"));
		UBQ.replace_residue( UBQ_term, core::conformation::Residue(cyd_rsd_term_type, true), true);

		// check safety of connections (from phil)
		core::chemical::ResidueType const & ubq_rsd_type( UBQ.residue_type( UBQ_term ) );
		runtime_assert(ubq_rsd_type.name() == cyd_rsd_term_type.name());
		core::Size const cyd_connid( 3 );
		core::Size const ubq_connid( 2 );

		runtime_assert( cyd_rsd_type.n_possible_residue_connections() == cyd_connid &&
			cyd_rsd_type.lower_connect_id() != cyd_connid &&
			cyd_rsd_type.upper_connect_id() != cyd_connid );

		runtime_assert( ubq_rsd_type.n_possible_residue_connections() == ubq_connid &&
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
		complex.append_residue_by_bond( UBQ.residue( UBQ_term ), true, ubq_connid, GTPase_cyd_, cyd_connid );
		//complex.dump_pdb("just1_complex.pdb");

		//not that this does anything
		complex.conformation().insert_ideal_geometry_at_residue_connection( GTPase_cyd_, cyd_connid );
		complex.conformation().insert_ideal_geometry_at_residue_connection( complex.total_residue(), ubq_connid );

		core::Size const ubq_pos( complex.total_residue() );
		core::id::AtomID const atom0( cyd_rsd_type.atom_index( "C" ), GTPase_cyd_ );
		core::id::AtomID const atom1( cyd_rsd_type.atom_index( "CA" ), GTPase_cyd_ );
		core::id::AtomID const atom2( cyd_rsd_type.atom_index( "CB" ), GTPase_cyd_ );
		core::id::AtomID const atom3( cyd_rsd_type.atom_index( "SG" ), GTPase_cyd_ );
		core::id::AtomID const atom4( ubq_rsd_type.atom_index( "SG"  ), ubq_pos );
		core::id::AtomID const atom5( ubq_rsd_type.atom_index( "CB" ), ubq_pos );
		core::id::AtomID const atom6( ubq_rsd_type.atom_index( "CA"  ), ubq_pos );

		//starting values derived from disulfide code (FullatomDisulfidePotential.cc and also $database/scoring/score_functions/disulfides/fa_SS_distance_score
		complex.conformation().set_torsion_angle( atom0, atom1, atom2, atom3, numeric::conversions::radians(52.0) ); //chi1
		complex.conformation().set_torsion_angle( atom1, atom2, atom3, atom4, 1.709486 ); //chi2
		complex.conformation().set_torsion_angle( atom2, atom3, atom4, atom5, 1.641427 ); //S-S bond
		complex.conformation().set_torsion_angle( atom3, atom4, atom5, atom6, 1.709486 ); //chi2
		complex.conformation().set_bond_length( atom3, atom4, 2.020); //S-S bond
		complex.conformation().set_bond_angle( atom2, atom3, atom4, 1.819120); //CB-SG-SG
		complex.conformation().set_bond_angle( atom3, atom4, atom5, 1.819120); //SG-SG-CB

		//   TR << "real disulf dist" << complex.xyz(atom3).distance(complex.xyz(atom4)) << std::endl;
		//   TR << GTPase_cyd_ << complex.xyz(atom3) << std::endl;
		//   TR << ubq_pos << complex.xyz(atom4) << std::endl;
		//   TR << GTPase_cyd_ << " " << ubq_pos << std::endl;
		//   complex.dump_scored_pdb("just1_complex2.pdb", *fullatom_scorefunction_);
		//   TR << "real disulf dist" << complex.xyz(atom3).distance(complex.xyz(atom4)) << std::endl;

		//now add the rest of ubiquitin
		for ( core::Size i=UBQ_term-1; i>= 1; --i ) {
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
		atomIDs[5] = core::id::AtomID( ubq_rsd_type.atom_index("SG" ), complexlength );
		atomIDs[6] = core::id::AtomID( ubq_rsd_type.atom_index("CB" ), complexlength );
		atomIDs[7] = core::id::AtomID( ubq_rsd_type.atom_index("CA" ), complexlength );
		atomIDs[8] = core::id::AtomID( ubq_rsd_type.atom_index("C" ), complexlength );


		////////////////////////////extra bodies/////////////////////////////////////////////////
		//The purpose of this code is to allow for static extra things in the system; for its original
		//incarnations, to allow for a RING domain to occlude its site on E2 and a GAP (+MG +GDP) to
		//occlude ras

		//check if extra bodies exist
		if ( basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user() == true ) {
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
		//small/shear behave fine @ the last residue
		disulfide_mm_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
		for ( core::Size i(0), ntailres(basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::n_tail_res]); i<ntailres; ++i ) { //slightly irregular < comparison because C-terminus is functionally zero-indexed
			disulfide_mm_->set_bb((complexlength-i), true);
		}
		//disulfide_mm_->set(core::id::TorsionID(complexlength, core::id::BB, core::id::phi_torsion), true);
		//disulfide_mm_->set(core::id::TorsionID(complexlength, core::id::BB, core::id::psi_torsion), true);
		//disulfide_mm_->set_bb(complexlength-1, true);
		//disulfide_mm_->set_bb(complexlength-2, true);
		//disulfide_mm_->set(complex.atom_tree().torsion_angle_dof_id(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5]), false);

		//setup loop
		std::set< core::Size > loop_posns;
		if ( basic::options::option[ basic::options::OptionKeys::loops::loop_file ].user() ) {
			loop_ = *( protocols::loops::Loops( true ).begin());
			TR << "loop " <<  loop_ << std::endl;
			//set up interface-plus-neighbors-positions operation
			for ( core::Size j(loop_.start()), end(loop_.stop()); j <= end; ++j ) {
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
		prevent->include_residue(GTPase_cyd_);
		prevent->include_residue(complexlength);
		task_factory_->push_back(prevent);

		if ( false ) {
			using core::pose::metrics::PoseMetricCalculatorOP;
			std::string const interface_calc("UBQGTPase_InterfaceNeighborDefinitionCalculator");
			std::string const neighborhood_calc("UBQGTPase_NeighborhoodByDistanceCalculator");
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( interface_calc, PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator( core::Size(1), core::Size(2)) ) );
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( neighborhood_calc, PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( loop_posns ) ) );

			//this is the constructor parameter for the calculator - pairs of calculators and calculations to perform
			utility::vector1< std::pair< std::string, std::string> > calcs_and_calcns;
			calcs_and_calcns.push_back(std::make_pair(interface_calc, "interface_residues"));
			calcs_and_calcns.push_back(std::make_pair(neighborhood_calc, "neighbors"));

			using protocols::toolbox::task_operations::RestrictByCalculatorsOperation;
			task_factory_->push_back(TaskOperationCOP( new RestrictByCalculatorsOperation( calcs_and_calcns ) ));

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
			for ( core::Size i(1); i<=E2_end; ++i ) {
				regions[E2_index].insert(i);
			}

			//ubiquitin tail
			//complexlength = end of ubiquitin in ubq+e2 pose
			core::Size const tail_begin(complexlength-basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::n_tail_res]+1); //odd construction accounts for functional zero-indexing of the tail
			regions.push_back(empty); //insert a new set to work with
			core::Size const tail_index(2);
			for ( core::Size i(complexlength); i>=tail_begin; --i ) {
				regions[tail_index].insert(i);
			}

			//ubiquitin core+tail
			//including the tail in both groups ensures it will always repack
			regions.push_back(empty); //insert a new set to work with
			core::Size const ubq_index(3);
			for ( core::Size i(E2_end+1); i<=complexlength; ++i ) {
				regions[ubq_index].insert(i);
			}

			//this will double-count loop residues with E2 - but that's fine, it just ensures they pack no matter what
			if ( !loop_posns.empty() ) {
				regions.push_back(loop_posns);
			}

			//if extra bodies exist, we are adding them to the first chain set
			if ( basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user() == true ) {
				apps::public1::scenarios::chemically_conjugated_docking::pack_extra_bodies(extra_bodies_chains_, complex, regions[E2_index], TR);
			}

			//make all pairs of groups (without replacement
			//if you have 1, 2, 3, 4; make 1-2, 1-3, 1-4, 2-3, 2-4, 3-4
			core::Size const num_regions(regions.size());
			utility::vector1< std::pair< std::set<core::Size>, std::set<core::Size> > > vector_of_pairs;
			for ( core::Size first_group(1); first_group < num_regions; ++first_group ) {
				for ( core::Size second_group(first_group+1); second_group <= num_regions; ++second_group ) {
					vector_of_pairs.push_back(std::make_pair(regions[first_group], regions[second_group]));
				}
			}

			//check contents of vector_of_pairs
			core::Size const num_pairs(vector_of_pairs.size());
			for ( core::Size i(1); i<=num_pairs; ++i ) {
				core::Size const
					onestart(*(vector_of_pairs[i].first.begin())),
					onestop(*(vector_of_pairs[i].first.rbegin())),
					twostart(*(vector_of_pairs[i].second.begin())),
					twostop(*(vector_of_pairs[i].second.rbegin()));

				TR << "IGNC will compare group " << onestart << "-" << onestop << " with " << twostart << "-" << twostop << std::endl;

				// if (basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user() == true) {
				//  TR.Error << "upcoming debug me non-contiguous set errors are not really errors if using extra bodies" << std::endl;
				// }

				// core::Size guess(onestart);
				// for(std::set<core::Size>::const_iterator iter(vector_of_pairs[i].first.begin()), end(vector_of_pairs[i].first.end()); iter != end; ++iter) {
				//  if(guess++ != *iter) TR.Error << "non-contiguous set, debug me!" << std::endl;
				//  TR << *iter << std::endl;
				// }
				// guess = twostart;
				// for(std::set<core::Size>::const_iterator iter(vector_of_pairs[i].second.begin()), end(vector_of_pairs[i].second.end()); iter != end; ++iter) {
				//  if(guess++ != *iter) TR.Error << "non-contiguous set, debug me!" << std::endl;
				//  TR << *iter << std::endl;
				// }

			}

			if ( basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user() == true ) {
				TR << "Those group labels do not take the piling of extra bodies into the first nonmoving group into account" << std::endl;
			}

			//check if calculator exists; create if not
			std::string const calc("IGNC_UBQ_Gp_CYD-CYD");
			if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists(calc) ) {
				core::pose::metrics::CalculatorFactory::Instance().remove_calculator(calc);
				TR.Error << "removed a PoseMetricCalculator " << calc << ", track down why" << std::endl;
			}
			using core::pose::metrics::PoseMetricCalculatorOP;
			core::pose::metrics::CalculatorFactory::Instance().register_calculator( calc, PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::InterGroupNeighborsCalculator(vector_of_pairs) ) );

			//now that calculator exists, add the sucker to the TaskFactory via RestrictByCalculatorsOperation
			utility::vector1< std::pair< std::string, std::string> > calculators_used;
			std::pair< std::string, std::string> IGNC_cmd( calc, "neighbors" );
			calculators_used.push_back( IGNC_cmd );
			task_factory_->push_back( TaskOperationCOP( new protocols::toolbox::task_operations::RestrictByCalculatorsOperation( calculators_used ) ) );

		}

		//create constraints
		core::scoring::constraints::add_constraints_from_cmdline_to_pose( starting_pose_ ); //protected internally if no constraints

	}
	virtual ~UBQ_GTPase_disulfide_Mover(){};

	virtual
	void
	apply( core::pose::Pose & pose ){
		if ( !init_for_input_yet_ ) init_on_new_input();

		pose = starting_pose_;

		//TR << "foldtree, movemap: " << std::endl;
		//core::kinematics::simple_visualize_fold_tree_and_movemap( pose.fold_tree(), *movemap_, TR);

		TR << "foldtree, " << '\n' << pose.fold_tree() << std::flush;

		/////////////////fullatom Monte Carlo//////////////////////////////////////////////////////////
		//make the monte carlo object
		using protocols::moves::MonteCarlo;
		using protocols::moves::MonteCarloOP;
		using basic::options::option;
		using namespace basic::options::OptionKeys::AnchoredDesign;
		MonteCarloOP mc( new MonteCarlo( pose, *fullatom_scorefunction_, option[ refine_temp ].value() ) );

		//////////////////////////Small/ShearMovers////////////////////////////////////////////////////////
		protocols::simple_moves::BackboneMoverOP small_mover( new protocols::simple_moves::SmallMover(disulfide_mm_, 0.8, 1) );
		small_mover->angle_max( 'H', 4.0 );
		small_mover->angle_max( 'E', 4.0 );
		small_mover->angle_max( 'L', 4.0 );

		protocols::simple_moves::BackboneMoverOP shear_mover( new protocols::simple_moves::ShearMover(disulfide_mm_, 0.8, 1) );
		shear_mover->angle_max( 'H', 4.0 );
		shear_mover->angle_max( 'E', 4.0 );
		shear_mover->angle_max( 'L', 4.0 );

		//handled by SidechainMover (?)
		// protocols::simple_moves::TorsionDOFMoverOP DOF_mover_chi1(new protocols::simple_moves::TorsionDOFMover);
		// DOF_mover_chi1->set_DOF(atomIDs[1], atomIDs[2], atomIDs[3], atomIDs[4]);
		// DOF_mover_chi1->check_mmt(true);
		// DOF_mover_chi1->temp(0.4);
		// DOF_mover_chi1->set_angle_range(-180, 180);
		// DOF_mover_chi1->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_chi2( new protocols::simple_moves::TorsionDOFMover );
		DOF_mover_chi2->set_DOF(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5]);
		DOF_mover_chi2->check_mmt(true);
		DOF_mover_chi2->temp(0.4);
		DOF_mover_chi2->set_angle_range(-180, 180);
		DOF_mover_chi2->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_disulfide( new protocols::simple_moves::TorsionDOFMover );
		DOF_mover_disulfide->set_DOF(atomIDs[3], atomIDs[4], atomIDs[5], atomIDs[6]);
		DOF_mover_disulfide->check_mmt(true);
		DOF_mover_disulfide->temp(0.4);
		DOF_mover_disulfide->set_angle_range(-180, 180);
		DOF_mover_disulfide->tries(1000);

		protocols::simple_moves::TorsionDOFMoverOP DOF_mover_psi( new protocols::simple_moves::TorsionDOFMover );
		DOF_mover_psi->set_DOF(atomIDs[4], atomIDs[5], atomIDs[6], atomIDs[7]);
		DOF_mover_psi->check_mmt(true);
		DOF_mover_psi->temp(0.4);
		DOF_mover_psi->set_angle_range(-180, 180);
		DOF_mover_psi->tries(1000);

		//Well, I forget why this is disabled, but it was for publication.  You can re-enable it.
		// protocols::simple_moves::TorsionDOFMoverOP DOF_mover_phi(new protocols::simple_moves::TorsionDOFMover);
		// DOF_mover_phi->set_DOF(atomIDs[5], atomIDs[6], atomIDs[7], atomIDs[8]);
		// DOF_mover_phi->check_mmt(true);
		// DOF_mover_phi->temp(0.4);
		// DOF_mover_phi->set_angle_range(-180, 180);
		// DOF_mover_phi->tries(1000);

		//Also add a SidechainMover for CYD (I hope...)
		//set up "pack only the moving conjugates" packer task
		utility::vector1< bool > repack_residues(pose.total_residue(), false); //this could be member data
		repack_residues[GTPase_cyd_] = true;
		repack_residues[pose.total_residue()] = true;
		core::pack::task::PackerTaskOP SC_task(core::pack::task::TaskFactory::create_packer_task(pose) );
		SC_task->restrict_to_residues(repack_residues);
		SC_task->restrict_to_repacking(); //SCmover will design, oops
		//and the mover
		protocols::simple_moves::sidechain_moves::SidechainMoverOP SC_mover( new protocols::simple_moves::sidechain_moves::SidechainMover() );
		SC_mover->set_change_chi_without_replacing_residue(true);
		SC_mover->set_task(SC_task);

		protocols::moves::RandomMoverOP backbone_mover( new protocols::moves::RandomMover() );
		backbone_mover->add_mover(small_mover, 2.0);
		backbone_mover->add_mover(shear_mover, 1.0);
		//   backbone_mover->add_mover(DOF_mover_chi1, 0.75);
		backbone_mover->add_mover(DOF_mover_chi2, 0.75);
		backbone_mover->add_mover(DOF_mover_disulfide, 0.75);
		backbone_mover->add_mover(DOF_mover_psi, 0.75);
		//   backbone_mover->add_mover(DOF_mover_phi, 0.75);
		backbone_mover->add_mover(SC_mover, 3.0);

		///////////////////////////loop movement/////////////////////////////////////////////////////
		if ( loop_.stop() - loop_.start() >= 3 ) { //empty loop; skip it!
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

		/////////////////////////minimize backbone DOFs//////////////////////////////////////////////
		using protocols::simple_moves::MinMoverOP;
		using protocols::simple_moves::MinMover;
		protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover(
			disulfide_mm_,
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
			disulfide_mm_,
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
			if ( (i % repack_cycles == 0) || (i == refine_applies) ) { //full repack
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
		if ( score > basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::scorefilter].value() ) {
			set_last_move_status(protocols::moves::FAIL_RETRY);
			TR << "total score filter failed; score " << score << std::endl;
			return;
		}

		//these interface analyses are less interpretable in the three-body case
		if ( !basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::extra_bodies].user() ) {

			//filter on interface SASA - requires some hacking to break up disulfide
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
			if ( mv_delta_sasa.value() < basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::SASAfilter].value() ) {
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
		job_me->add_string_real_pair("Ncysteine_chi1_C-CA-CB-SG", degrees(pose.atom_tree().torsion_angle(atomIDs[1], atomIDs[2], atomIDs[3], atomIDs[4])));
		job_me->add_string_real_pair("Ncysteine_chi2_CA-CB-SG-SG", degrees(pose.atom_tree().torsion_angle(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5])));
		job_me->add_string_real_pair("disulfide_CB-SG-SG-CB", degrees(pose.atom_tree().torsion_angle(atomIDs[3], atomIDs[4], atomIDs[5], atomIDs[6])));
		job_me->add_string_real_pair("Ccysteine_psi_SG-SG-CB-CA", degrees(pose.atom_tree().torsion_angle(atomIDs[4], atomIDs[5], atomIDs[6], atomIDs[7])));
		job_me->add_string_real_pair("Ccysteine_phi_SG-CB-CA-C", degrees(pose.atom_tree().torsion_angle(atomIDs[5], atomIDs[6], atomIDs[7], atomIDs[8])));
		//job_me->add_string_real_pair("disulf_dist", pose.xyz(atomIDs[4]).distance(pose.xyz(atomIDs[5])));

		set_last_move_status(protocols::moves::MS_SUCCESS);
		if ( basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::publication].value() ) {
			apps::public1::scenarios::chemically_conjugated_docking::create_extra_output(pose, TR, !basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::pdz].value(), GTPase_cyd_);
		}
		return;
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return protocols::moves::MoverOP( new UBQ_GTPase_disulfide_Mover );
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }

	virtual
	std::string
	get_name() const { return "UBQ_GTPase_disulfide_Mover"; }

private:
	bool init_for_input_yet_;

	core::scoring::ScoreFunctionOP fullatom_scorefunction_;
	core::pack::task::TaskFactoryOP task_factory_;
	core::kinematics::MoveMapOP disulfide_mm_;
	//  core::kinematics::MoveMapOP loop_mm_;
	//  core::kinematics::MoveMapOP all_mm_;

	protocols::loops::Loop loop_;

	/// @brief vector contains atomIDs for disulfide bond and atoms before/after bond to determine various torsions
	utility::vector1< core::id::AtomID > atomIDs;

	core::pose::Pose starting_pose_; //maintained from run to run

	std::string const InterfaceSasaDefinition_; //calculator name

	protocols::analysis::InterfaceAnalyzerMoverOP IAM_;

	core::Size GTPase_cyd_; //converted to member data for sharing between setup and apply

	/// @brief used to track which chains are "extra" nonmoving bodies in extra bodies mode
	utility::vector1< core::Size > extra_bodies_chains_;
};

typedef utility::pointer::shared_ptr< UBQ_GTPase_disulfide_Mover > UBQ_GTPase_disulfide_MoverOP;

int main( int argc, char* argv[] )
{
	try {
		//initialize options
		devel::init(argc, argv);
		basic::prof_reset();

		if ( basic::options::option[ basic::options::OptionKeys::in::file::s ].user()
				|| basic::options::option[ basic::options::OptionKeys::in::file::l ].user()
				|| basic::options::option[ basic::options::OptionKeys::in::file::silent ].user() ) {
			utility_exit_with_message("do not use an input PDB with this protocol (program uses internally); use -UBQpdb and -GTPasepdb instead");
		}

		protocols::jd2::JobDistributor::get_instance()->go(protocols::moves::MoverOP( new UBQ_GTPase_disulfide_Mover ));

		basic::prof_show();
		TR << "NOTE on interpreting results: the interface energies are somewhat broken due to there being a disulfide across the interface; the bond is still scored in the separated state, which leads to enormously bad bond-length energies.  This biases by about 6000 energy units.  Relative ranking (which is all you should do anyway) is still correct." << std::endl;
		TR << "************************d**o**n**e**************************************" << std::endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
