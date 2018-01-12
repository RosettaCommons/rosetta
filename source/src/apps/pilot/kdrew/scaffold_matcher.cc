// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


//kdrew: this app attempts to find the orientation of a scaffold that best satisifies the positions of hotspot residues (provided as inverse rotamer libraries).
//scaffold in this case is generally meant as a small noncanonical backbone with a handful of designable positions.


//Headers are generally organized by either what they do or where they come from.  This organization is first core library headers, then protocols library, then utility stuff.


// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_trials.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/BackboneStubConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/PyMOLMover.hh>

#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/ncbb/oop/OopPatcher.hh>
#include <protocols/ncbb/hbs/HbsPatcher.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/ncbb/oop/OopDockDesignProtocol.hh>
#include <protocols/ncbb/NcbbDockDesignProtocol.hh>
#include <protocols/ncbb/util.hh>

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>


using namespace basic::options;
using namespace protocols;


// tracer - used to replace cout
static basic::Tracer TR( "Hotspot Placement" );

namespace scaffold_matcher {
BooleanOptionKey const pymol( "scaffold_matcher::pymol" );
BooleanOptionKey const keep_history( "scaffold_matcher::keep_history" );

StringOptionKey const primary_hs_stub_lib( "scaffold_matcher::primary_hs_stub_lib" );
FileVectorOptionKey const ancillary_hs_stub_libs( "scaffold_matcher::ancillary_hs_stub_libs" );
StringOptionKey const hstarget( "scaffold_matcher::hstarget" );

BooleanOptionKey const hs_repack_only( "scaffold_matcher::hs_repack_only" );
RealOptionKey const stub_constraint_strength("scaffold_matcher::stub_constraint_strength");
}

class HotspotPlacementMover : public moves::Mover {

public:

	//default ctor
	HotspotPlacementMover(): Mover("HotspotPlacementMover"){}

	//default dtor
	~HotspotPlacementMover() override= default;

	//methods
	void setup_pert_foldtree( core::pose::Pose & pose);
	void setup_pert_foldtree_byres( core::pose::Pose & pose, Size dock_jump_pos_pro, Size dock_jump_pos_pep );

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override { return "OopHotspotPlacementMover"; }

};

using HotspotPlacementMoverOP = utility::pointer::shared_ptr<HotspotPlacementMover>;
using HotspotPlacementMoverCOP = utility::pointer::shared_ptr<const HotspotPlacementMover>;


int
main( int argc, char* argv[] )
{
	try {
		option.add( scaffold_matcher::pymol, "Set up pymol mover. Default false" ).def(false);
		option.add( scaffold_matcher::keep_history, "Keep history in pymol. Requires ohp::pymol set to true. Default false" ).def(false);
		option.add( scaffold_matcher::primary_hs_stub_lib, "The path to primary hotspot stub pdb" ).def( "" );
		option.add( scaffold_matcher::ancillary_hs_stub_libs, "The path to additional hotspot stub pdbs" ).def( "" );
		option.add( scaffold_matcher::hstarget, "Target pdb" ).def( "" );
		option.add( scaffold_matcher::hs_repack_only, "Keep scaffold fixed, only repack sidechains" ).def( true );
		option.add( scaffold_matcher::stub_constraint_strength, "strength of stub constraints (reasonable: 0.5 - 5.0)" ).def( 5.0 );

		devel::init(argc, argv);

		//create mover instance
		HotspotPlacementMoverOP HP_mover( new HotspotPlacementMover() );

		//call job distributor
		protocols::jd2::JobDistributor::get_instance()->go( HP_mover );

	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main


void
HotspotPlacementMover::apply(
	core::pose::Pose & scaffold_pose
)
{

	//core::pose::Pose start_pose = core::pose::Pose(scaffold_pose);
	core::pose::Pose target_pose;
	core::import_pose::pose_from_file(target_pose, option[ scaffold_matcher::hstarget].value() , core::import_pose::PDB_file);

	//kdrew: create starting pose by combining the target pose and the scaffold pose
	core::pose::Pose combined_pose = core::pose::Pose(target_pose);
	for ( Size i = 1; i <= scaffold_pose.size(); ++i ) {
		if ( 1 == i ) {
			combined_pose.append_residue_by_jump(scaffold_pose.residue(i), combined_pose.size(), "", "", true );
		} else {
			combined_pose.append_residue_by_bond( scaffold_pose.residue(i) );
		}
	}

	std::stringstream startpdbname;
	startpdbname << "start_pose.pdb";
	combined_pose.dump_pdb( startpdbname.str() );

	core::pose::Pose start_pose;
	core::import_pose::pose_from_file(start_pose, startpdbname.str() , core::import_pose::PDB_file);

	setup_pert_foldtree(start_pose);

	if ( option[ scaffold_matcher::pymol ].value() ) {
		protocols::moves::PyMOLObserverOP pymover = protocols::moves::AddPyMOLObserver(start_pose, option[ scaffold_matcher::keep_history ].value() );
	}

	core::scoring::ScoreFunctionOP place_hs_score_fxn( new core::scoring::ScoreFunction() );

	core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();
	core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

	/*
	core::scoring::ScoreFunctionOP bump_scorefxn = new core::scoring::ScoreFunction;
	bump_scorefxn->reset();
	bump_scorefxn->set_weight( core::scoring::fa_rep, 1.0 );
	*/

	// get peptide start and end positions
	Size pep_start( start_pose.conformation().chain_begin( 2 ) ); Size pep_end( start_pose.size() );
	TR << "pep_start: " << pep_start << " pep_end: " << pep_end << std::endl;

	////kdrew: automatically find ncbb positions
	utility::vector1< core::Size > ncbb_seq_positions = core::pose::ncbb::initialize_ncbbs(start_pose);

	//kdrew: turn on atom pair constraints if not already on
	if ( ncbb_seq_positions.size() > 0 && score_fxn->has_zero_weight( core::scoring::atom_pair_constraint ) ) {
		score_fxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	}


	core::kinematics::MoveMapOP place_hs_mm( new core::kinematics::MoveMap() );
	//kdrew: set backbone of target false and backbone of scaffold true
	place_hs_mm->set_bb( false );
	place_hs_mm->set_chi( false );
	place_hs_mm->set_jump( 1, true );
	protocols::minimization_packing::MinMoverOP place_hs_min( new protocols::minimization_packing::MinMover( place_hs_mm, place_hs_score_fxn, "lbfgs_armijo_nonmonotone", 0.001, true ) );

	//Primary Hotspot Setup
	protocols::hotspot_hashing::HotspotStubSetOP hotspot_stub_setOP( new protocols::hotspot_hashing::HotspotStubSet );
	hotspot_stub_setOP->read_data( option[ scaffold_matcher::primary_hs_stub_lib ].value() );

	//Secondary Hotspot Setup
	utility::vector1<protocols::hotspot_hashing::HotspotStubSetOP> ancillary_stubs;
	utility::vector1<std::string> ancillary_stub_locations =  option[ scaffold_matcher::ancillary_hs_stub_libs ]();

	for ( core::Size i = 1; i <= ancillary_stub_locations.size(); ++i ) {
		protocols::hotspot_hashing::HotspotStubSetOP hotspot_stub_set2OP( new protocols::hotspot_hashing::HotspotStubSet );
		hotspot_stub_set2OP->read_data(ancillary_stub_locations[i]);
		ancillary_stubs.push_back (hotspot_stub_set2OP);
	}

	//kdrew: packer used before final scoring of scaffold placement
	using core::pack::task::operation::TaskOperationCOP;
	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
	tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
	//kdrew: do not do design, makes NATAA if res file is not specified
	core::pack::task::operation::RestrictToRepackingOP rtrp( new core::pack::task::operation::RestrictToRepacking() );
	tf->push_back( rtrp );
	protocols::minimization_packing::PackRotamersMoverOP packer( new protocols::minimization_packing::PackRotamersMover() );
	packer->task_factory( tf );
	packer->score_function( score_fxn );
	//packer->score_function( bump_scorefxn );

	//kdrew: mc_temp = 1.0, turn off rotation 0.0 and translation 0.0, outer loop number = 10, no_design = true, final_design = false, pymol=false, keep_history=false
	protocols::ncbb::NcbbDockDesignProtocolOP NDDP_mover( new protocols::ncbb::NcbbDockDesignProtocol( score_fxn, 0.0, 0.0, 0.0, 10, true, false, false, false ) );
	protocols::ncbb::setup_filter_stats();

	// create move map for minimization
	core::kinematics::MoveMapOP min_mm( new core::kinematics::MoveMap() );
	min_mm->set_bb( true );
	min_mm->set_chi( true );
	min_mm->set_jump( 1, true );

	// create minimization mover
	minimization_packing::MinMoverOP min( new minimization_packing::MinMover( min_mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01, true ) );

	//kdrew: setup default packer task for both L and D residues
	utility::vector1< bool > aas(20,false);
	core::chemical::ResidueTypeSetCOP rs( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
	core::pack::task::PackerTaskOP packer_task = core::pack::task::TaskFactory::create_packer_task( start_pose );
	for ( core::Size resnum=pep_start; resnum <= pep_end; ++resnum ) {
		core::chemical::ResidueType rtype = start_pose.residue_type( resnum );
		std::list< core::chemical::ResidueTypeCOP > allowed_aas = packer_task->residue_task(resnum).allowed_residue_types();
		if ( rtype.is_d_aa() ) {
			packer_task->nonconst_residue_task(resnum).restrict_absent_canonical_aas(aas);
			for ( std::list< core::chemical::ResidueTypeCOP >::const_iterator restype = allowed_aas.begin();
					restype != allowed_aas.end(); ++restype ) {

				core::chemical::ResidueType allowed_rtype =  **restype;
				core::chemical::ResidueType const & chiral_aa = protocols::simple_moves::chiral::get_chiral_residue_type(allowed_rtype, protocols::simple_moves::chiral::D_CHIRALITY, *rs);
				TR << "chiral_aa.name(): "<< chiral_aa.name() << " chiral_aa.name3(): " << chiral_aa.name3()<< std::endl;
				packer_task->nonconst_residue_task(resnum).allow_noncanonical_aa(chiral_aa.name3(),*rs);
			}
		}
	}

	for ( core::Size resnum=pep_start; resnum <= pep_end; ++resnum ) {
		core::chemical::ResidueType rtype = start_pose.residue_type( resnum );

		// Check that this position is allowed to be used for stub constraints
		if ( ! packer_task->pack_residue(resnum) ) continue;


		std::list< core::chemical::ResidueTypeCOP > allowed_aas = packer_task->residue_task( resnum ).allowed_residue_types();
		for ( std::list< core::chemical::ResidueTypeCOP >::const_iterator restype = allowed_aas.begin();
				restype != allowed_aas.end(); ++restype ) {

			TR << "residues allowed: "<< (*restype)->name3() << " is allowed AA" << std::endl;
			// Loop over all stubs with this restype
			Size hs_index = 0;

			protocols::hotspot_hashing::HotspotStubSet::Hotspots res_stub_set( hotspot_stub_setOP->retrieve( (*restype )->name3() ) );

			for ( auto & hs_stub : res_stub_set ) {
				Pose pose( start_pose );

				//kdrew: setup fold tree centered around resnum
				Size pro_start( pose.conformation().chain_begin( 1 ) );
				Size pro_end( pose.conformation().chain_end( 1 ) );
				Size dock_jump_pos_pro( core::pose::residue_center_of_mass( pose, pro_start, pro_end ) );
				setup_pert_foldtree_byres(pose, dock_jump_pos_pro, resnum);

				hs_index++;

				TR << "resnum: " << resnum << " stub: "<< hs_stub.second->residue()->name() << " is residue stub #" << hs_index << std::endl;

				core::Size fixed_res(1);
				core::id::AtomID fixed_atom = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res );

				//kdrew: setup atom pair constraint between stub and resnum, CA, CB, N, C
				core::Real distance = 0.0;
				core::Real stdev = 0.05;
				core::scoring::func::HarmonicFuncOP harm_func( new core::scoring::func::HarmonicFunc( distance, stdev ) );

				//kdrew: probably cannot do atom pair constraint because stub is not really part of pose,
				//kdrew: try xyz constraint?
				core::id::AtomID aid_scaffold_CA( pose.residue( resnum ).atom_index("CA"), resnum );
				core::id::AtomID aid_scaffold_CB( pose.residue( resnum ).atom_index("CB"), resnum );
				core::id::AtomID aid_scaffold_C( pose.residue( resnum ).atom_index("C"), resnum );
				core::id::AtomID aid_scaffold_N( pose.residue( resnum ).atom_index("N"), resnum );

				core::scoring::constraints::ConstraintCOP coord_CA_cst( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::CoordinateConstraint( aid_scaffold_CA, fixed_atom, hs_stub.second->residue()->xyz( hs_stub.second->residue()->atom_index("CA") ), harm_func ) ) );
				core::scoring::constraints::ConstraintCOP coord_CB_cst( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::CoordinateConstraint( aid_scaffold_CB, fixed_atom, hs_stub.second->residue()->xyz( hs_stub.second->residue()->atom_index("CB") ), harm_func ) ) );
				core::scoring::constraints::ConstraintCOP coord_C_cst( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::CoordinateConstraint( aid_scaffold_C, fixed_atom, hs_stub.second->residue()->xyz( hs_stub.second->residue()->atom_index("C") ), harm_func ) ) );
				core::scoring::constraints::ConstraintCOP coord_N_cst( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::CoordinateConstraint( aid_scaffold_N, fixed_atom, hs_stub.second->residue()->xyz( hs_stub.second->residue()->atom_index("N") ), harm_func ) ) );

				pose.add_constraint( coord_CA_cst );
				pose.add_constraint( coord_CB_cst );
				pose.add_constraint( coord_C_cst );
				pose.add_constraint( coord_N_cst );

				//kdrew: only use coord constraint score term to place scaffold
				TR << "Place_hs full score without coord constraints: " << (*place_hs_score_fxn)(pose) << std::endl;

				place_hs_score_fxn->set_weight( core::scoring::coordinate_constraint, 1.0 );

				core::Real place_hs_full_coord_score = (*place_hs_score_fxn)(pose);

				TR << "Place_hs full score with coord cst: " << place_hs_full_coord_score << std::endl;
				TR << "Place_hs total score with coord constraints is: " << pose.energies().total_energies()[ core::scoring::total_score ] << std::endl;
				TR << "Place_hs coord cst score : " << pose.energies().total_energies()[ core::scoring::coordinate_constraint] << std::endl;

				place_hs_min->apply(pose);

				place_hs_full_coord_score = (*place_hs_score_fxn)(pose);
				TR << "after min: Place_hs full score with coord cst: " << place_hs_full_coord_score << std::endl;
				TR << "after min: Place_hs total score with coord constraints is: " << pose.energies().total_energies()[ core::scoring::total_score ] << std::endl;
				TR << "after min: Place_hs coord cst score : " << pose.energies().total_energies()[ core::scoring::coordinate_constraint] << std::endl;

				/*
				std::stringstream pdbname;
				pdbname << "placed_min_pose_" << resnum << "_" << hs_index << ".pdb";
				pose.dump_scored_pdb( pdbname.str(), *place_hs_score_fxn );
				*/

				place_hs_score_fxn->set_weight( core::scoring::coordinate_constraint, 0.0 );
				TR << "coord constraint removal successful: " << pose.remove_constraint( coord_CA_cst, true ) << std::endl;
				TR << "coord constraint removal successful: " << pose.remove_constraint( coord_CB_cst, true ) << std::endl;
				TR << "coord constraint removal successful: " << pose.remove_constraint( coord_C_cst, true ) << std::endl;
				TR << "coord constraint removal successful: " << pose.remove_constraint( coord_N_cst, true ) << std::endl;

				//kdrew: now use full score term to evaluate "goodness" of scaffold placement
				TR << "Full score without stub: " << (*score_fxn)(pose) << std::endl;

				if ( score_fxn->has_zero_weight( core::scoring::backbone_stub_constraint) ) {
					score_fxn->set_weight( core::scoring::backbone_stub_constraint, 1.0 );
				}

				core::Real stub_bonus_value = hs_stub.second->bonus_value();
				core::Real const force_constant = option[ scaffold_matcher::stub_constraint_strength ].value();

				//core::scoring::constraints::ConstraintCOP bbstub_cst = new core::scoring::constraints::BackboneStubConstraint( pose, resnum, fixed_atom, *(hs_stub->second->residue()), stub_bonus_value, force_constant, "CB","CA","N","C" );

				//kdrew: TODO: make target_residue a pose
				core::pose::Pose stub_pose = core::pose::Pose();
				stub_pose.append_residue_by_bond( *(hs_stub.second->residue()) );
				core::scoring::constraints::ConstraintCOP bbstub_cst( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::BackboneStubConstraint( pose, resnum, fixed_atom, stub_pose, 1, stub_bonus_value, force_constant) ) );

				pose.add_constraint( bbstub_cst );

				core::scoring::constraints::ConstraintCOPs additional_hs_constraints;
				//kdrew: put additional constraints to satisify other hotspot residues
				//kdrew: for each residue in scaffold (peptide)
				for ( core::Size resnum2=pep_start; resnum2 <= pep_end; ++resnum2 ) {

					core::chemical::ResidueType rtype2 = pose.residue_type( resnum2 );

					if ( resnum != resnum2 ) {
						utility::vector1< core::scoring::constraints::ConstraintCOP > ambig_csts;
						//kdrew: Check packer (resfile) that this position is allowed to be used for stub constraints
						if ( ! packer_task->pack_residue(resnum2) ) continue;

						//kdrew: for each allowed type at this position
						std::list< core::chemical::ResidueTypeCOP > allowed_aas2 = packer_task->residue_task( resnum2 ).allowed_residue_types();
						for ( std::list< core::chemical::ResidueTypeCOP >::const_iterator restype2 = allowed_aas2.begin();
								restype2 != allowed_aas2.end(); ++restype2 ) {
							//kdrew: for each inverse rotamer in stub set
							for ( core::Size jjj = 1; jjj <= ancillary_stubs.size(); ++jjj ) {
								protocols::hotspot_hashing::HotspotStubSet::Hotspots res_stub_set2( ancillary_stubs[jjj]->retrieve( (*restype2 )->name3() ) );
								for ( auto & hs_stub2 : res_stub_set2 ) {
									core::pose::Pose stub2_pose = core::pose::Pose();
									stub2_pose.append_residue_by_bond( *(hs_stub2.second->residue()) );
									core::Real stub_bonus_value2 = hs_stub2.second->bonus_value();
									ambig_csts.push_back(
										core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::BackboneStubConstraint( pose,
										resnum2,
										fixed_atom,
										stub2_pose,
										1, //kdrew: assuming stub is in pose with only one residue, change for peptoids
										stub_bonus_value2,
										force_constant ) ) )
									);
								}
							}
						}
						if ( ambig_csts.size() > 0 ) {
							TR << "adding ambiguous constraints to resnum: "<<resnum2<<std::endl;
							additional_hs_constraints.push_back( core::scoring::constraints::ConstraintCOP( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AmbiguousConstraint(ambig_csts) ) ) );
						}
					}
				}

				additional_hs_constraints = pose.add_constraints( additional_hs_constraints );
				core::Real full_stub_score = (*score_fxn)(pose);

				TR << "Full score with stub: " << full_stub_score << std::endl;
				TR << "Total score with long-range constraints is: " << pose.energies().total_energies()[ core::scoring::total_score ] << std::endl;
				TR << "Stub score with long-range constraints is: " << pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;

				if ( option[ scaffold_matcher::hs_repack_only].value() ) {
					packer->apply( pose );
				} else { // sample scaffold conformations
					NDDP_mover->apply( pose );
				}

				//min->apply( pose );

				//kdrew: need to use restrict to repack to properly do rotamer trials
				//core::pack::task::PackerTaskOP bump_packer_task = core::pack::task::TaskFactory::create_packer_task( pose );
				//core::pack::rotamer_trials( pose , *bump_scorefxn, bump_packer_task );

				full_stub_score = (*score_fxn)(pose);

				std::stringstream repack_pdbname;
				utility::file::FileName input_pdbname( scaffold_pose.pdb_info()->name() );
				repack_pdbname << input_pdbname.base() << "_repack_pose_" << resnum << "_" << hs_index << ".pdb";
				pose.dump_scored_pdb( repack_pdbname.str(), *score_fxn );

				TR << "after repack: Full score with stub: " << full_stub_score << std::endl;
				TR << "after repack: Total score with long-range constraints is: " << pose.energies().total_energies()[ core::scoring::total_score ] << std::endl;
				TR << "after repack: Stub score with long-range constraints is: " << pose.energies().total_energies()[ core::scoring::backbone_stub_constraint ] << std::endl;

				score_fxn->set_weight( core::scoring::backbone_stub_constraint, 0.0 );
				TR << "constraint removal successful: " << pose.remove_constraint( bbstub_cst, true ) << std::endl;
				std::cout << std::endl;
			}
		}
	}

}

void
HotspotPlacementMover::setup_pert_foldtree(
	core::pose::Pose & pose
)
{
	// get the start and end for both chains
	Size pro_start( pose.conformation().chain_begin( 1 ) );
	Size pro_end( pose.conformation().chain_end( 1 ) );
	Size pep_start( pose.conformation().chain_begin( 2 ) );
	Size pep_end( pose.conformation().chain_end( 2 ) );

	// get jump positions based on the center of mass of the chains
	Size dock_jump_pos_pro( core::pose::residue_center_of_mass( pose, pro_start, pro_end ) );
	Size dock_jump_pos_pep( core::pose::residue_center_of_mass( pose, pep_start, pep_end ) );

	setup_pert_foldtree_byres(pose, dock_jump_pos_pro, dock_jump_pos_pep);

}
void
HotspotPlacementMover::setup_pert_foldtree_byres(
	core::pose::Pose & pose,
	Size dock_jump_pos_pro,
	Size dock_jump_pos_pep
)
{
	using namespace core;
	using namespace core::kinematics;

	// get current fold tree
	FoldTree f( pose.fold_tree() );
	f.clear();

	// get the start and end for both chains
	Size pro_start( pose.conformation().chain_begin( 1 ) );
	Size pro_end( pose.conformation().chain_end( 1 ) );
	Size pep_start( pose.conformation().chain_begin( 2 ) );
	Size pep_end( pose.conformation().chain_end( 2 ) );

	// build fold tree
	Size jump_index( f.num_jump() + 1 );
	f.add_edge( pro_start, dock_jump_pos_pro, Edge::PEPTIDE );
	f.add_edge( dock_jump_pos_pro, pro_end, Edge::PEPTIDE );
	f.add_edge( pep_start, dock_jump_pos_pep, Edge::PEPTIDE );
	f.add_edge( dock_jump_pos_pep, pep_end, Edge::PEPTIDE );
	f.add_edge( dock_jump_pos_pro, dock_jump_pos_pep, jump_index );

	//kdrew: make the jump go through the CB atom
	core::chemical::AtomType pro_atom_type = pose.residue( dock_jump_pos_pro ).atom_type( pose.residue( dock_jump_pos_pro ).nbr_atom() );
	std::string pro_atom_name = pro_atom_type.element();
	//kdrew: need to change this if using peptoids or anything that does not have a traditional CB
	//awatkins: this works for peptoids and scaffolds containnig beta aas because they at least still name it CB
	f.set_jump_atoms( jump_index, pro_atom_name, "CB" );
	//f.set_jump_atoms( jump_index, pose.residue( dock_jump_pos_pro ).atom_type( pose.residue( dock_jump_pos_pro ).nbr_atom() ).element(), "CB" );

	// set pose foldtree to foldtree we just created
	f.reorder(1);
	f.check_fold_tree();
	assert( f.check_fold_tree() );

	std::cout << "AFTER: " << f << std::endl;

	pose.fold_tree( f );
}
