// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/chemically_conjugated_docking/UBQ_GTPaseMover.cc
/// @brief A mover typically used for binding ubiquitin to a substrate protein.
/// @author Steven Lewis and Hope Anderson

// Unit headers
#include <protocols/chemically_conjugated_docking/UBQ_GTPaseMover.hh>
#include <protocols/chemically_conjugated_docking/UBQ_GTPaseMoverCreator.hh>
#include <protocols/chemically_conjugated_docking/Gp_quantification_metrics.hh>
#include <protocols/chemically_conjugated_docking/Gp_extra_bodies.hh>

// Core headers
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
#include <protocols/task_operations/RestrictByCalculatorsOperation.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/FileJobOutputter.hh>
#include <protocols/jd2/wwPDBJobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>

//movers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/moves/MoverContainer.hh> //Sequence Mover
#include <protocols/minimization_packing/RotamerTrialsMover.hh>
#include <protocols/minimization_packing/TaskAwareMinMover.hh>
#include <protocols/moves/OutputMovers.hh> //pdbdumpmover
#include <protocols/simple_moves/TorsionDOFMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>

#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <protocols/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/pose_metric_calculators/InterGroupNeighborsCalculator.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>

// Basic/Utility headers
#include <numeric/conversions.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <core/select/util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/prof.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/keys/chemically_conjugated_docking.OptionKeys.gen.hh>
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
static basic::Tracer TR( "protocols.chemically_conjugated_docking.UBQ_GTPaseMover" );

namespace protocols {
namespace chemically_conjugated_docking {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
UBQ_GTPaseMover::UBQ_GTPaseMover():
	protocols::moves::Mover( UBQ_GTPaseMover::mover_name() ),
	fullatom_scorefunction_(/* NULL */),
	task_factory_(/* NULL */),
	amide_mm_(/* NULL */),
	loop_(), //we want default ctor
	atomIDs(8, core::id::AtomID::BOGUS_ATOM_ID() ),
	InterfaceSasaDefinition_("InterfaceSasaDefinition_" + 1),
	IAM_(protocols::analysis::InterfaceAnalyzerMoverOP( new protocols::analysis::InterfaceAnalyzerMover )),
	extra_bodies_(false),
	n_tail_res_(3),
	scorefilter_(),
	SASAfilter_(),
	UBQpdb_(),
	selector_()
{
	//set up fullatom scorefunction
	using namespace core::scoring;
	fullatom_scorefunction_ = get_score_function();
	core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *fullatom_scorefunction_ ); //protected if(option) internally

	TR << "Using fullatom scorefunction from commandline:\n" << *fullatom_scorefunction_;

	using namespace core::pose::metrics;
	using namespace protocols::pose_metric_calculators;
	//magic number: chains 1 and 2; set up interface SASA calculator
	if ( !CalculatorFactory::Instance().check_calculator_exists( InterfaceSasaDefinition_ ) ) {
		CalculatorFactory::Instance().register_calculator( InterfaceSasaDefinition_, PoseMetricCalculatorOP( new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator(core::Size(1), core::Size(2)) ));
	}

	IAM_->set_use_centroid_dG(false);
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
UBQ_GTPaseMover::~UBQ_GTPaseMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief init_on_new_input system allows for initializing these details the first time apply() is called.  the job distributor will reinitialize the whole mover when the input changes (a freshly constructed mover, which will re-run this on first apply().
// modified by Hope Anderson so that init_on_new_input is called every time; this is so the residueSelector chooses a new lysine on GTPase with every apply.
void
UBQ_GTPaseMover::initialize(core::pose::Pose & GTPase) {

	TR << "Creating new starting pose." << std::endl;

	TR << "GTPase pose has size " << GTPase.size() << std::endl;

	// selects lysine
	utility::vector1< bool > subset = selector_ -> apply(GTPase);
	utility::vector1< core::Size > res = core::select::get_residues_from_subset(subset);

	core::Size GTPase_residue = res[1];

	TR << "Lysine number " << GTPase_residue << " will be bound to UBQ." << std::endl;

	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////Create starting complex pose/////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	TR << "Creating starting pose..." << std::endl;

	//read poses
	core::Size const GTPaselength = GTPase.size();

	core::pose::Pose UBQ;
	try{
		core::import_pose::pose_from_file( UBQ, UBQpdb_ , core::import_pose::PDB_file);
	}
catch ( utility::excn::Exception & e ) {
	std::stringstream error_msg;
	error_msg << "Could not create pose for UBQ file. Did you remember to provide a file path for UBQ?\n";
	error_msg << e.msg();
	throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
}

	core::Size const UBQlength = UBQ.size();

	//determine cysteine target
	runtime_assert(GTPase.conformation().num_chains() == 1);
	char const GTPasechain(GTPase.pdb_info()->chain(1));
	GTPase_lys_ = GTPase.pdb_info()->pdb2pose(GTPasechain, GTPase_residue);
	//runtime_assert(GTPase.residue_type(GTPase_lys_).aa() == core::chemical::aa_lys);

	//determine c_term target on UBQ
	core::Size const UBQ_term = UBQlength;

	//strip C-term from UBQ - best to do this with a full replace to re-draw the carboxyl oxygen
	//UBQ.dump_pdb("pre-removeUQB.pdb");
	core::chemical::ResidueTypeSetCOP fa_standard(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
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

	runtime_assert( lyx_rsd_type.n_possible_residue_connections() == lyx_connid &&
		lyx_rsd_type.lower_connect_id() != lyx_connid &&
		lyx_rsd_type.upper_connect_id() != lyx_connid );

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
	complex.append_residue_by_bond( UBQ.residue( UBQ_term ), true, ubq_connid, GTPase_lys_, lyx_connid );
	//complex.dump_pdb("just1_complex.pdb");

	//not that this does anything
	complex.conformation().insert_ideal_geometry_at_residue_connection( GTPase_lys_, lyx_connid );

	core::Size const ubq_pos( complex.size() );
	core::id::AtomID const atom0( lyx_rsd_type.atom_index( "CG" ), GTPase_lys_ );
	core::id::AtomID const atom1( lyx_rsd_type.atom_index( "CD" ), GTPase_lys_ );
	core::id::AtomID const atom2( lyx_rsd_type.atom_index( "CE" ), GTPase_lys_ );
	core::id::AtomID const atom3( lyx_rsd_type.atom_index( "NZ" ), GTPase_lys_ );
	core::id::AtomID const atom4( ubq_rsd_type.atom_index( "C"  ), ubq_pos );
	core::id::AtomID const atom5( ubq_rsd_type.atom_index( "CA" ), ubq_pos );
	core::id::AtomID const atom6( ubq_rsd_type.atom_index( "N"  ), ubq_pos );

	//starting values derived from the peptide bond and a straight-out lysine
	for ( core::Size chi(1); chi<=4; ++chi ) complex.set_chi(chi, GTPase_lys_, 180);
	//complex.conformation().set_torsion_angle( atom0, atom1, atom2, atom3, numeric::conversions::radians(106.5) );
	//complex.conformation().set_torsion_angle( atom1, atom2, atom3, atom4, numeric::conversions::radians(-60.0) );
	complex.conformation().set_torsion_angle( atom2, atom3, atom4, atom5, numeric::conversions::radians(180.0) );
	complex.conformation().set_torsion_angle( atom3, atom4, atom5, atom6, numeric::conversions::radians(135.0) );
	//complex.dump_pdb("just1_complex2.pdb");

	//now add the rest of ubiquitin
	for ( core::Size i=UBQ_term-1; i>= 1; --i ) {
		complex.prepend_polymer_residue_before_seqpos( UBQ.residue(i), GTPaselength+1, false );
	}

	core::Size const complexlength( complex.size());
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
	if ( extra_bodies_ == true ) {
		extra_bodies_chains_ = protocols::chemically_conjugated_docking::add_extra_bodies(complex);
	}

	starting_pose_ = complex;


	protocols::jd2::JobOP job_me = protocols::jd2::JobDistributor::get_instance()->current_job();
	protocols::jd2::JobOutputterOP job_out(protocols::jd2::JobDistributor::get_instance()->job_outputter());
	//protocols::jd2::FileJobOutputterOP job_files = utility::pointer::static_pointer_cast< protocols::jd2::FileJobOutputter >(job_out);
	//protocols::jd2::wwPDBJobOutputterOP job_pdbs = utility::pointer::static_pointer_cast< protocols::jd2::wwPDBJobOutputter >(job_files);
	//protocols::jd2::wwPDBJobOutputterOP job_pdbs(protocols::jd2::JobDistributor::get_instance()->job_outputter());

	job_out -> other_pose(job_me, starting_pose_, "starting_pose");

	//starting_pose_.dump_pdb("starting_complex.pdb");
	//TR << "starting pose finished." << std::endl;

	///////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////Finish starting complex pose/////////////////////////////////////////////
	//////////////////////Start creating move accessory data///////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////////////////

	//setup MoveMaps
	//small/shear behave improperly @ the last residue - psi is considered nonexistent and the wrong phis apply.
	amide_mm_ = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
	for ( core::Size i(1), ntailres(n_tail_res_); i<ntailres; ++i ) { //slightly irregular < comparison because C-terminus is functionally zero-indexed
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
	//task_factory_->push_back( new protocols::task_operations::RestrictToInterfaceOperation );
	task_factory_->push_back( TaskOperationCOP( new IncludeCurrent ) );
	//prevent repacking at linkage lysine!
	PreventRepackingOP prevent( new PreventRepacking );
	prevent->include_residue(GTPase_lys_);
	task_factory_->push_back(prevent);

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
	core::Size const tail_begin(complexlength-n_tail_res_+1); //odd construction accounts for functional zero-indexing of the tail
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
	if ( extra_bodies_ == true ) {
		protocols::chemically_conjugated_docking::pack_extra_bodies(extra_bodies_chains_, complex, regions[E2_index]);
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
	}

	if ( extra_bodies_ == true ) {
		TR << "Those group labels do not take the piling of extra bodies into the first nonmoving group into account" << std::endl;
	}

	//check if calculator exists; create if not
	std::string const calc("IGNC_chemically_conjugated_docking");
	if ( core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists(calc) ) {
		core::pose::metrics::CalculatorFactory::Instance().remove_calculator(calc);
		//TR.Error << "removed a PoseMetricCalculator " << calc << ", track down why" << std::endl;
	}
	using core::pose::metrics::PoseMetricCalculatorOP;
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( calc, PoseMetricCalculatorOP( new protocols::pose_metric_calculators::InterGroupNeighborsCalculator(vector_of_pairs) ) );

	//now that calculator exists, add the sucker to the TaskFactory via RestrictByCalculatorsOperation
	utility::vector1< std::pair< std::string, std::string> > calculators_used;
	std::pair< std::string, std::string> IGNC_cmd( calc, "neighbors" );
	calculators_used.push_back( IGNC_cmd );
	task_factory_->push_back( TaskOperationCOP( new protocols::task_operations::RestrictByCalculatorsOperation( calculators_used ) ) );

	//create constraints
	core::scoring::constraints::add_constraints_from_cmdline_to_pose( starting_pose_ ); //protected internally if no constraints

}

/// @brief Apply the mover
void
UBQ_GTPaseMover::apply( core::pose::Pose& pose){
	initialize(pose);

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
	protocols::simple_moves::BackboneMoverOP small_mover( new protocols::simple_moves::SmallMover(amide_mm_, 0.8, 1) );
	small_mover->angle_max( 'H', 4.0 );
	small_mover->angle_max( 'E', 4.0 );
	small_mover->angle_max( 'L', 4.0 );

	protocols::simple_moves::BackboneMoverOP shear_mover( new protocols::simple_moves::ShearMover(amide_mm_, 0.8, 1) );
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

	protocols::simple_moves::TorsionDOFMoverOP DOF_mover_isopeptide( new protocols::simple_moves::TorsionDOFMover );
	DOF_mover_isopeptide->set_DOF(atomIDs[3], atomIDs[4], atomIDs[5], atomIDs[6]);
	DOF_mover_isopeptide->check_mmt(true);
	DOF_mover_isopeptide->temp(0.4);
	DOF_mover_isopeptide->set_angle_range(-180, 180);
	DOF_mover_isopeptide->tries(1000);

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

	//Also add a SidechainMover for LYX (I hope...)
	//set up "pack only the moving conjugate" packer task
	utility::vector1< bool > repack_residues(pose.size(), false); //this could be member data
	repack_residues[GTPase_lys_] = true;
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
	//   backbone_mover->add_mover(DOF_mover_chi1, 0.75); //SC mover will handle this DOF
	//   backbone_mover->add_mover(DOF_mover_chi2, 0.75); //SC mover will handle this DOF
	backbone_mover->add_mover(DOF_mover_isopeptide, 0.75);
	backbone_mover->add_mover(DOF_mover_psi, 0.75);
	backbone_mover->add_mover(DOF_mover_phi, 0.75);
	backbone_mover->add_mover(SC_mover, 1.0);

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
	using protocols::minimization_packing::MinMoverOP;
	using protocols::minimization_packing::MinMover;
	protocols::minimization_packing::MinMoverOP min_mover( new protocols::minimization_packing::MinMover(
		amide_mm_,
		fullatom_scorefunction_,
		basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
		0.01,
		true /*use_nblist*/ ) );

	/////////////////////////////////rotamer trials mover///////////////////////////////////////////
	using protocols::minimization_packing::RotamerTrialsMoverOP;
	using protocols::minimization_packing::EnergyCutRotamerTrialsMover;
	protocols::minimization_packing::RotamerTrialsMoverOP rt_mover( new protocols::minimization_packing::EnergyCutRotamerTrialsMover(
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
	protocols::minimization_packing::PackRotamersMoverOP pack_mover( new protocols::minimization_packing::PackRotamersMover );
	pack_mover->task_factory( task_factory_ );
	pack_mover->score_function( fullatom_scorefunction_ );

	protocols::minimization_packing::MinMoverOP min_mover_pack( new protocols::minimization_packing::MinMover(
		amide_mm_,
		fullatom_scorefunction_,
		basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
		0.01,
		true /*use_nblist*/ ) );

	using protocols::minimization_packing::TaskAwareMinMoverOP;
	using protocols::minimization_packing::TaskAwareMinMover;
	protocols::minimization_packing::TaskAwareMinMoverOP TAmin_mover( new protocols::minimization_packing::TaskAwareMinMover(min_mover_pack, task_factory_) );

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
UBQ_GTPaseMover::analyze_and_filter(core::pose::Pose & pose){

	//Filter on total score
	core::Real const score((*fullatom_scorefunction_)(pose));
	if ( score > scorefilter_ ) {
		set_last_move_status(protocols::moves::FAIL_RETRY);
		TR << "total score filter failed; score " << score << std::endl;
		return;
	}

	//these interface analyses are less interpretable in the three-body case
	if ( !extra_bodies_ ) {

		//filter on interface SASA - requires some hacking to break up isopeptide
		core::pose::Pose copy(pose);
		//hack the pose up for analysis purposes
		core::Size const cbreak(copy.conformation().chain_end(1));
		using core::kinematics::Edge;
		core::kinematics::FoldTree main_tree(copy.size());
		main_tree.clear();
		main_tree.add_edge(Edge(1, cbreak, Edge::PEPTIDE));
		main_tree.add_edge(Edge(cbreak+1, copy.size(), Edge::PEPTIDE));
		main_tree.add_edge(Edge(cbreak, cbreak+1, 1));
		main_tree.reorder(1);
		//TR << main_tree << std::endl;
		copy.fold_tree(main_tree);

		//Filter on SASA
		basic::MetricValue< core::Real > mv_delta_sasa;
		copy.metric(InterfaceSasaDefinition_, "delta_sasa", mv_delta_sasa);
		if ( mv_delta_sasa.value() < SASAfilter_ ) {
			set_last_move_status(protocols::moves::FAIL_RETRY);
			TR << "interface SASA filter failed; SASA " << mv_delta_sasa.value() << std::endl;
			return;
		}

		//passed filters; run IAM
		IAM_->apply(copy);
	}

	//print mobile region fine-grained data
	using numeric::conversions::degrees;
	protocols::jd2::JobOP job_me = protocols::jd2::JobDistributor::get_instance()->current_job();
	job_me->add_string_real_pair("lysine_chi3_CG-CD-CE-NZ", degrees(pose.atom_tree().torsion_angle(atomIDs[1], atomIDs[2], atomIDs[3], atomIDs[4])));
	job_me->add_string_real_pair("lysine_chi4_CD-CE-NZ-C", degrees(pose.atom_tree().torsion_angle(atomIDs[2], atomIDs[3], atomIDs[4], atomIDs[5])));
	job_me->add_string_real_pair("amide_CE-NZ-C-CA", degrees(pose.atom_tree().torsion_angle(atomIDs[3], atomIDs[4], atomIDs[5], atomIDs[6])));
	job_me->add_string_real_pair("glycine_psi_NZ-C-CA-N", degrees(pose.atom_tree().torsion_angle(atomIDs[4], atomIDs[5], atomIDs[6], atomIDs[7])));
	job_me->add_string_real_pair("glycine_phi_C-CA-N-C", degrees(pose.atom_tree().torsion_angle(atomIDs[5], atomIDs[6], atomIDs[7], atomIDs[8])));

	set_last_move_status(protocols::moves::MS_SUCCESS);
	if ( basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::publication].value() ) {
		protocols::chemically_conjugated_docking::create_extra_output(pose, !basic::options::option[basic::options::OptionKeys::chemically_conjugated_docking::pdz].value(), GTPase_lys_);
	}
	return;
}

void
UBQ_GTPaseMover::make_index_selector(core::Size const index){
	core::select::residue_selector::ResidueSelectorCOP res_op(new core::select::residue_selector::ResidueIndexSelector(index));
	selector_ = res_op;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
UBQ_GTPaseMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
UBQ_GTPaseMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const &)
{

	// simple tags
	set_UBQpdb(tag -> getOption < std::string > ("UBQpdb",""));
	set_extra_bodies(tag -> getOption < bool > ("extraBodies", false));
	set_n_tail_res(tag -> getOption < core::Size > ("numTailRes", 3));
	set_scorefilter(tag -> getOption < core::Real > ("scorefilter", 0));
	set_SASAfilter(tag -> getOption < core::Real > ("SASAfilter", 1000));

	// setting the selector
	std::string const selectorname = tag -> getOption < std::string > ("binding_lysine","");
	core::Size const lys_num = tag -> getOption <core::Size> ("lys_num", 0);

	if ( selectorname=="" && lys_num==0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception, "Please enter either a residue selector or the location for binding.");
	} else if ( lys_num==0 ) { // setting selector from selector created in script
		try {
			set_selector(data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selectorname ));
			TR << "Using residue selector " << selectorname << std::endl;
		} catch ( utility::excn::Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selectorname << "' in the DataMap.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}
	} else { // setting selector from a given residue index
		//set_selector((new core::select::residue_selector::ResidueIndexSelector(lys_num));
		make_index_selector(lys_num);

		TR << "Created selector for index " << lys_num << std::endl;
	}
	debug_assert(get_selector());
}

void UBQ_GTPaseMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.
	attlist
		+XMLSchemaAttribute("UBQpdb", xs_string, "File path to UBQ pdb file")
		+XMLSchemaAttribute("binding_lysine", xs_string, "A residue selector giving the lysine residue on the substrate to which UBQ will bind.")
		+XMLSchemaAttribute("extraBodies", xs_boolean, "Whether there are extra bodies in the simulation.")
		+XMLSchemaAttribute("numTailRes", xsct_non_negative_integer, "Number of tail residues on the ubiquitin.")
		+XMLSchemaAttribute("scorefilter", xs_decimal, "Filter for total score")
		+XMLSchemaAttribute("SASAfilter", xs_decimal, "Filter for SASA")
		+XMLSchemaAttribute("lys_num", xs_integer, "An integer giving the lysine residue to which UBQ will bind. Use this or binding_lysine.");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "A mover for modeling docking via a LYX bond.", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
UBQ_GTPaseMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new UBQ_GTPaseMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
UBQ_GTPaseMover::clone() const
{
	return protocols::moves::MoverOP( new UBQ_GTPaseMover( *this ) );
}

std::string UBQ_GTPaseMover::get_name() const {
	return mover_name();
}

std::string UBQ_GTPaseMover::mover_name() {
	return "UBQ_GTPaseMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
UBQ_GTPaseMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new UBQ_GTPaseMover );
}

std::string
UBQ_GTPaseMoverCreator::keyname() const
{
	return UBQ_GTPaseMover::mover_name();
}

void UBQ_GTPaseMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	UBQ_GTPaseMover::provide_xml_schema( xsd );
}

// getters
//

core::scoring::ScoreFunctionOP UBQ_GTPaseMover::get_fullatom_scorefunction() const{
	return fullatom_scorefunction_;
}

void UBQ_GTPaseMover::set_fullatom_scorefunction(core::scoring::ScoreFunctionOP const & sfxn){
	fullatom_scorefunction_ = sfxn;
}

core::pack::task::TaskFactoryOP UBQ_GTPaseMover::get_task_factory() const{
	return task_factory_;
}

void UBQ_GTPaseMover::set_task_factory(core::pack::task::TaskFactoryOP const & tasks){
	task_factory_ = tasks;
}

core::kinematics::MoveMapOP UBQ_GTPaseMover::get_amide_mm() const{
	return amide_mm_;
}

void UBQ_GTPaseMover::set_amide_mm(core::kinematics::MoveMapOP const & mm){
	amide_mm_ = mm;
}

protocols::loops::Loop const & UBQ_GTPaseMover::get_loop() const{
	return loop_;
}

void UBQ_GTPaseMover::set_loop(protocols::loops::Loop const & loop){
	loop_ = loop;
}

utility::vector1<core::id::AtomID> const & UBQ_GTPaseMover::get_atom_IDs() const{
	return atomIDs;
}

void UBQ_GTPaseMover::set_atom_IDs(utility::vector1<core::id::AtomID> const & ids){
	atomIDs = ids;
}

core::pose::Pose const & UBQ_GTPaseMover::get_starting_pose() const{
	return starting_pose_;
}

std::string const & UBQ_GTPaseMover::get_InterfaceSasaDefinition() const{
	return InterfaceSasaDefinition_;
}

protocols::analysis::InterfaceAnalyzerMoverOP UBQ_GTPaseMover::get_IAM() const{
	return IAM_;
}

void UBQ_GTPaseMover::set_IAM(protocols::analysis::InterfaceAnalyzerMoverOP const & iam){
	IAM_ = iam;
}

core::Size UBQ_GTPaseMover::get_GTPase_lys() const{
	return GTPase_lys_;
}

void UBQ_GTPaseMover:: set_GTPase_lys(core::Size const lys){
	GTPase_lys_ = lys;
}

utility::vector1<core::Size> const & UBQ_GTPaseMover::get_extra_bodies_chains() const{
	return extra_bodies_chains_;
}

void UBQ_GTPaseMover::set_extra_bodies_chains(utility::vector1<core::Size> const & chains){
	extra_bodies_chains_ = chains;
}

bool UBQ_GTPaseMover::get_extra_bodies() const{
	return extra_bodies_;
}

void UBQ_GTPaseMover::set_extra_bodies(bool const hasExtra){
	extra_bodies_ = hasExtra;
}

core::Size UBQ_GTPaseMover::get_n_tail_res() const{
	return n_tail_res_;
}

void UBQ_GTPaseMover::set_n_tail_res(core::Size const numTails){
	n_tail_res_ = numTails;
}

core::Real UBQ_GTPaseMover::get_scorefilter() const{
	return scorefilter_;
}

void UBQ_GTPaseMover::set_scorefilter(core::Real const maxScore){
	scorefilter_ = maxScore;
}

core::Real UBQ_GTPaseMover::get_SASAfilter() const{
	return SASAfilter_;
}

void UBQ_GTPaseMover::set_SASAfilter(core::Real const minSASA){
	SASAfilter_ = minSASA;
}

std::string const & UBQ_GTPaseMover::get_UBQpdb() const{
	return UBQpdb_;
}

void UBQ_GTPaseMover::set_UBQpdb(std::string const & ubqFile){
	UBQpdb_ = ubqFile;
}

core::select::residue_selector::ResidueSelectorCOP UBQ_GTPaseMover::get_selector() const{
	return selector_;
}

void UBQ_GTPaseMover::set_selector(core::select::residue_selector::ResidueSelectorCOP const & select){
	selector_ = select;
}

//////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, UBQ_GTPaseMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //chemically_conjugated_docking
