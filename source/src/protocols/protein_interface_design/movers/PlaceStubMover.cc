// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/PlaceStubMover.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu), Eva-Maria Strauch (evas01@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/PlaceStubMover.hh>
#include <protocols/protein_interface_design/movers/PlaceStubMoverCreator.hh>
#include <boost/foreach.hpp>

// Project Headers
#include <utility/tag/Tag.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>
#include <numeric/xyzVector.hh>
#include <protocols/protein_interface_design/movers/PlaceUtils.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/datacache/DataMapObj.hh>
#include <protocols/moves/ResId.hh>

#include <core/scoring/constraints/BackboneStubConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/pack_rotamers.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>

#include <protocols/moves/Mover.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/moves/MoverStatus.hh>

#include <protocols/protein_interface_design/movers/SetAtomTree.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

// Utility Headers
#include <utility/exit.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Unit Headers
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/protein_interface_design/design_utils.hh>
#include <protocols/protein_interface_design/movers/BuildAlaPose.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <core/pack/task/operation/TaskOperations.hh>

// C++ headers
#include <map>
#include <algorithm>

#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/NumericTraits.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>
#include <protocols/simple_filters/EnergyPerResidueFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_moves/DesignRepackMover.hh>
using namespace core::scoring;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.PlaceStubMover" );
static thread_local basic::Tracer stats_TR( "STATS.PlaceStubMover" );
static thread_local basic::Tracer TR_debug( "DEBUG.PlaceStubMover" );


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;
using namespace core;

std::string
PlaceStubMoverCreator::keyname() const
{
	return PlaceStubMoverCreator::mover_name();
}

protocols::moves::MoverOP
PlaceStubMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PlaceStubMover );
}

std::string
PlaceStubMoverCreator::mover_name()
{
	return "PlaceStub";
}

PlaceStubMover::PlaceStubMover() :
	simple_moves::DesignRepackMover( PlaceStubMoverCreator::mover_name() ),
	score_threshold_( 0.0 ),
	host_chain_( 2 ),
	stub_set_( /* NULL */ ),
	add_constraints_( false ),
	after_placement_filter_( /* NULL */ ),
	final_filter_( /* NULL */ ),
	coord_cst_func_( /* NULL */ ),
	default_fold_tree_( /* NULL */ ),
	leave_coord_csts_after_placement_( false ),
	place_scaffold_( false ),
	max_cb_cb_dist_( 4.0 ),
	hurry_( true ),
	triage_positions_( true ),
	stub_energy_threshold_( 1.0 ),
	residue_level_tasks_for_placed_hotspots_( /* NULL */ ),
	residue_numbers_( /* NULL */ )
{
	coord_cst_std_.clear();
	disallowed_host_pos_.clear();
	stub_minimize_movers_.clear();
	design_movers_.clear();
	placed_stubs_.clear();
	curr_coordinate_constraints_.clear();
	previous_coordinate_constraints_.clear();
	saved_bb_constraints_.clear();
	saved_placed_stubs_.clear();
}

PlaceStubMover::PlaceStubMover(
	protocols::hotspot_hashing::HotspotStubSetOP stub_set,
	core::Real score_threshold,
	core::Size const host_chain,
	protocols::filters::FilterOP final_filter,
	bool const hurry/*=false*/,
	bool const triage_positions/*=true*/,
	core::Real stub_energy_threshold /*= 1.0*/
) :
	simple_moves::DesignRepackMover( PlaceStubMoverCreator::mover_name() ),
	score_threshold_( score_threshold ),
	host_chain_( host_chain ),
	hurry_(hurry),
	triage_positions_(triage_positions),
	stub_energy_threshold_(stub_energy_threshold)
{
	if( stub_set ) stub_set_ = protocols::hotspot_hashing::HotspotStubSetOP( new protocols::hotspot_hashing::HotspotStubSet( *stub_set ) );
	final_filter_ = final_filter->clone();
}

PlaceStubMover::~PlaceStubMover() {}

protocols::moves::MoverOP
PlaceStubMover::clone() const {
	return( protocols::moves::MoverOP( new PlaceStubMover( *this ) ) );
}


/// @details utility function that places a stub at an acceptor position and cleans up the pose
void
PlaceStubMover::place_stub( core::pose::Pose & pose, core::conformation::Residue const res_stub, core::Size const res_num )
{
	runtime_assert( res_num <= pose.total_residue() );
	core::Size const chain_begin( pose.conformation().chain_begin( host_chain_ ) );
	core::Size const chain_end  ( pose.conformation().chain_end  ( host_chain_ ) );
/// Using a default fold tree, rather than the existing atom tree, b/c the atom
/// tree might assume the existence of sidechain atoms that will no longer be available
/// after the residue replacement.
	pose.fold_tree( *default_fold_tree_ );
	std::pair< core::Size, bool > const res_cst( std::make_pair( res_num, add_constraints_ ) );
	placed_stubs_.push_back( res_cst );
	pose.replace_residue( res_num, res_stub, true );
	using namespace core::chemical;
//removing variant types can only be done if the residue would still be connected to the chain
	if( res_num < chain_end )
		core::pose::remove_upper_terminus_type_from_pose_residue( pose, res_num );
	if( res_num > chain_begin )
		core::pose::remove_lower_terminus_type_from_pose_residue( pose, res_num );
	pose.conformation().update_polymeric_connection( res_num , true); // o/w residues connections mess up
	if( res_num > chain_begin )
		pose.conformation().update_polymeric_connection( res_num - 1 , true);
/// return to stub_based atom tree
	if( place_scaffold_ ){
	// currently only supports host being on chain2
	// we append the stub residue to the target by jump and then impose that jump
	// on the target-scaffold complex
		runtime_assert( host_chain_==2 );// other options not supported yet

		core::pose::Pose partner_chain( *pose.split_by_chain( 1 ) );
		TR_debug<<"foldtree before append: "<<partner_chain.fold_tree()<<std::endl;
		partner_chain.append_residue_by_jump( res_stub, 1 );
		TR_debug<<"foldtree after append: "<<partner_chain.fold_tree()<<std::endl;
		core::kinematics::Jump const saved_jump( partner_chain.jump( partner_chain.num_jump() ) );

		core::kinematics::FoldTree pose_stub_tree;
		core::Size const rb_jump( 1 );
		{//setup pose-stub ft
			pose_stub_tree.clear();
			using namespace core::kinematics;
			pose_stub_tree.add_edge( 1, partner_chain.total_residue()-1, Edge::PEPTIDE );
			pose_stub_tree.add_edge( 1, partner_chain.total_residue(), rb_jump );
			pose_stub_tree.set_jump_atoms( rb_jump, partner_chain.residue( 1 ).atom_type( partner_chain.residue( 1 ).nbr_atom()).element(), "CB" );
			partner_chain.fold_tree( pose_stub_tree );
			TR_debug<<"pose stub tree: "<<partner_chain.fold_tree()<<std::endl;
		}//pose stub ft
		core::kinematics::Jump const pose_stub_jump( partner_chain.jump( rb_jump ) );
		core::kinematics::FoldTree new_ft;
		{//setup new ft for target-scaffold
			core::Size const host_chain_begin( pose.conformation().chain_begin( host_chain_ ) );
			core::Size const host_chain_end  ( pose.conformation().chain_end  ( host_chain_ ) );
			new_ft.clear();
			core::Size const jump_pos1( 1 );
			core::Size const jump_pos2( res_num );
			new_ft.clear();
			new_ft.add_edge( jump_pos1, jump_pos2, rb_jump );
			new_ft.add_edge( jump_pos1, pose.conformation().chain_end( 1 ), kinematics::Edge::PEPTIDE );
			if( jump_pos2 > host_chain_begin )
				new_ft.add_edge( pose.conformation().chain_begin( 2 ), jump_pos2, kinematics::Edge::PEPTIDE );
			if( jump_pos2 < host_chain_end )
				new_ft.add_edge( jump_pos2, pose.total_residue(), kinematics::Edge::PEPTIDE );
			new_ft.set_jump_atoms( rb_jump, pose.residue( 1 ).atom_type( pose.residue( 1 ).nbr_atom()).element(), "CB" );
			new_ft.reorder( 1 );
			new_ft.delete_self_edges();
			TR_debug<<"pose fold tree before imposing new jump "<<pose.fold_tree()<<std::endl;
			pose.fold_tree( new_ft );
			TR_debug<<"pose fold tree after imposing new jump "<<pose.fold_tree()<<std::endl;
		}//new ft
		pose.set_jump( 1, pose_stub_jump );
	}
	stub_based_atom_tree( pose, res_stub, 0.5/*cst_sdev*/ );
	pose.update_residue_neighbors();
}

/// @details minimize the rb orientation in the presence of a strong bb_stub_constraint potential, while reducing
/// all other attractive scores. fa_sol is also reduced to avoid blowing up the structure. Returns false if current
/// pose has no backbone_stub_constraint score. Removes stub constraints after minimization.
/// If more than one stub has already been placed then the jumps are held rigid during minimization
/// If a stub is specified, then a custom-made hotspot constraint is used
bool
PlaceStubMover::StubMinimize( core::pose::Pose & pose, protocols::hotspot_hashing::HotspotStubCOP stub/* = NULL*/, core::Size const host_residue/*= 0*/, bool const hurry /*=false*/ ){
	core::Size fixed_res(1);
	if( host_chain_ == 1 ) fixed_res = pose.total_residue();
	core::id::AtomID const fixed_atom_id = core::id::AtomID( pose.residue(fixed_res).atom_index("CA"), fixed_res );
	core::pack::task::PackerTaskOP stub_task = stub_set_->prepare_hashing_packer_task_( pose, host_chain_ );
	core::Size const host_chain_begin( pose.conformation().chain_begin( host_chain_ ) );
	core::Size const host_chain_end  ( pose.conformation().chain_end  ( host_chain_ ) );

	using namespace core::scoring::constraints;
	ConstraintCOPs const csts_before_min = pose.constraint_set()->get_all_constraints();
	core::Size const num_csts_before_min( csts_before_min.size() );

	using namespace core::scoring;
	ScoreFunctionOP scorefxn = get_score_function();
	ScoreFunctionOP stub_scorefxn( scorefxn->clone() );
	if( hurry ) {
		TR << "Speeding up StubMinimize..." << std::endl;
		stub_scorefxn->reset();
		stub_scorefxn->set_weight( backbone_stub_constraint, 10.0 );//is getting reset, if stub minimize mover is specified
		stub_scorefxn->set_weight( fa_rep, 0.44 );
		stub_scorefxn->set_weight( fa_dun, 0.56 );
		stub_scorefxn->set_weight( coordinate_constraint, 1.0 );
		//adding additional score types to ensure proper backbone parameters
	}
	else {stub_scorefxn->set_weight( fa_atr, 0.01 );
		stub_scorefxn->set_weight( backbone_stub_constraint, 10.0 );//is getting reset, if stub minimize mover is specified
		stub_scorefxn->set_weight( fa_sol, 0.01 );
		stub_scorefxn->set_weight( fa_pair, 0.01 );
		stub_scorefxn->set_weight( hbond_lr_bb, 0.01 );
		stub_scorefxn->set_weight( hbond_sr_bb, 0.01 );
		stub_scorefxn->set_weight( hbond_bb_sc, 0.01 );
		stub_scorefxn->set_weight( hbond_sc, 0.01 );
		stub_scorefxn->set_weight( coordinate_constraint, 1.0 );
	}
	using namespace protocols::hotspot_hashing;
	core::scoring::constraints::ConstraintCOPs stub_constraints;
	if( stub != NULL ){ //one stub-based constraint
		runtime_assert( host_residue );
		core::conformation::Residue const host_res( pose.conformation().residue( host_residue ) );
		core::Real dummy1, dummy2;
		if( !test_res_res_aln( host_res, *stub->residue(), dummy1, dummy2 ) ){
			TR_debug<<"Not minimizing towards this stub, b/c of large orientational discrepancy\n";
			return( false );
		}
	// setup the constraint based on a single stub <-> single host residue.
	// The cb force is computed analytically so as to (almost) maximise the steepness
	// of the potential while ascertaining that the potential is felt by the host residue
	// at the distance where it is positioned to begin with.
		TR<<"making single-stub based constraint\n";
		core::Real const distance( pose.residue( host_residue ).xyz( "CB" ).distance( stub->residue()->xyz( "CB" ) ) );
		core::Real const bonus( stub->bonus_value() );
	//I'm capping the value of cb_force to [0.05 .. 10.0]
		core::Real const cb_force( std::max( std::min( 10.0, -bonus / (distance * distance) - 0.05 ), 0.05 ) );
		TR << "Cb force: "<<cb_force<<std::endl;
		// I'm circumventing add_hotspot_constraints to pose and adding the constraint directly
		// since there's no ambiguity here and no need to switch to ala pose etc. And I don't
		// want all the quality control machinery to be applied to this stub; I know it's good.
		stub_constraints.push_back( core::scoring::constraints::ConstraintOP( new BackboneStubConstraint( pose, host_residue, fixed_atom_id, host_res, bonus, cb_force ) ) );
		stub_constraints = pose.add_constraints( stub_constraints );
	}// stub() != NULL
	else{ //multiple stubs
	/// Now circumventing add_hotspot_constraints b/c of the huge computational load of
	/// computing sasa etc.
		using namespace core::scoring::constraints;
		if( stub_set_->constraints().empty() ){
			for( core::Size resi=host_chain_begin; resi<=host_chain_end; ++resi )
				if( std::find( prevent_repacking_.begin(), prevent_repacking_.end(),resi ) != prevent_repacking_.end() )
					stub_task->nonconst_residue_task( resi ).prevent_repacking();
			TR<<"adding multiple stub constraints to pose\n";
			stub_set_->add_hotspot_constraints_to_pose( pose, fixed_atom_id, stub_task, stub_set_, 0.7/*cb force*/, 0/*worst allowed stub bonus*/, false/*apply self energies*/, 10.0/*bump cutoff*/, true/*apply ambiguous constraints*/ );
			(*stub_scorefxn)(pose);//for ->active_constraint to be set below.
		}
		else{//constraints not empty
			ConstraintCOPs constraints;
			constraints = pose.add_constraints( stub_set_->constraints() );
			core::Size const constraint_num( constraints.size() );
			TR<<"adding "<<constraint_num<<" stub constraints to pose"<<std::endl;
			(*stub_scorefxn)(pose);//for ->active_constraint to be set below.
			if( !prevent_repacking_.empty() ){
				ConstraintCOPs to_be_removed;
				for( ConstraintCOPs::const_iterator it = constraints.begin(); it != constraints.end(); ++it ){
					AmbiguousConstraintCOP cst = utility::pointer::dynamic_pointer_cast< core::scoring::constraints::AmbiguousConstraint const > ( (*it) );
					runtime_assert( cst != 0 );

					using namespace core::scoring::constraints;
					ConstraintCOP active_constraint = cst->active_constraint();

					if( active_constraint->type() == "BackboneStub" ){
						BackboneStubConstraintCOP bb_cst = utility::pointer::dynamic_pointer_cast< core::scoring::constraints::BackboneStubConstraint const > ( active_constraint );
						runtime_assert( bb_cst != 0 );
						if( std::find( prevent_repacking_.begin(), prevent_repacking_.end(), bb_cst->seqpos() ) != prevent_repacking_.end() )
							to_be_removed.push_back( *it ); //remove the entire ambiguous constraint, if the active constraint points to a non-repackable residue
					}
				}
				if( !to_be_removed.empty() ){
					pose.remove_constraints( to_be_removed );
					TR<<"removed "<<to_be_removed.size()<<" constraints\n";
				}
			}//!prevent_repacking_.empty()
		}//constraints not empty
	}//multiple stubs

	using namespace core::scoring;
	core::Size num_rigid_jumps( 0 );
	for( utility::vector1< std::pair< core::Size, bool > >::const_iterator it=placed_stubs_.begin(), end=placed_stubs_.end(); it!=end; ++it )
		if( !it->second ) num_rigid_jumps++;
	if( num_rigid_jumps > 1 ){
		TR<<"****WARNING WARNING**** Activating chainbreak weight, implying two rigid jumps.\n";
		TR<<"Not fully supported by PlaceStub\n";
		stub_scorefxn->set_weight( chainbreak, 1.0 ); //***** Is this the correct chainbreak weight to activate here?
	}
	//currently we only support one rigid jump. The remainder have to be set
	//by constraints. There are no explicit checks in the code for this and
	//it is up to the user to make sure that all's well...

	simple_filters::ScoreTypeFilter const stf( stub_scorefxn, backbone_stub_constraint, 1.0 );
	core::Real const before_min( stf.compute( pose ) );
	if( before_min >= -0.0001 ){
		if( stub != NULL ){
			TR<<"no bb stub constraint score even though I computed it analytically! Ask Sarel what's wrong here."<<std::endl;;
			runtime_assert( stub != 0 );
		}
		protocols::hotspot_hashing::remove_hotspot_constraints_from_pose( pose );
		TR<<"removing stub constraints from pose\n";
		return false;
	}
	//for minimization (rb and sc of previous placed stubs)
	utility::vector1< bool > sc_min( pose.total_residue(), false );
	utility::vector1< bool > const no_min( pose.total_residue(), false );
	utility::vector1< core::Size > const no_targets;
	// minimize the sc of all placed stubs
	for( utility::vector1< std::pair< core::Size, bool > >::const_iterator iter=placed_stubs_.begin(), end=placed_stubs_.end(); iter!=end; ++iter ) sc_min[ iter->first ] = true;
	//getting angle alignments, CB distances and CB force
	using namespace protocols::hotspot_hashing;
	core::Real C_N_angle_before_minmize_mover(0);
	core::Real CB_CA_angle_before_minmize_mover(0);
	core::Real before_min_distance(0);
	core::Real curr_bonus(0);
	core::Size before_apply_stub_minimize_mover(0);
	core::Real before_min_cb_force(0);

	if( stub != 0 ){
		before_min_distance = pose.residue( host_residue ).xyz( "CB" ).distance( stub->residue()->xyz( "CB" ) );
		curr_bonus = stub->bonus_value();
		before_apply_stub_minimize_mover = pose.constraint_set()->get_all_constraints().size();//coordinate constraint shouldnt change
		before_min_cb_force = -curr_bonus / (before_min_distance *before_min_distance );

		//scoring the pose before testing the alignment to get the angles....
		test_res_res_aln(  pose.residue( host_residue ), *stub->residue(), C_N_angle_before_minmize_mover, CB_CA_angle_before_minmize_mover );
		(*stub_scorefxn )( pose ); //to update values

		//minimizing stub using user-defined movers or a default minimization (rb and sc of placed stubs)
		typedef std::pair< simple_moves::DesignRepackMoverOP, core::Real > DesignMoverRealPair;
		if ( stub_minimize_movers_.size() ){
			TR<<"entering movers for stub minimization....\n";
			//minimize rb and sc of previous place
			TR<<"performing an initial rb minimization as well as sc minimization of previously placed stubs (if present) \n";
			MinimizeInterface( pose, stub_scorefxn, no_min/*bb*/, sc_min, min_rb()/*rb*/, false /*optimize foldtree*/, no_targets, true/*simultaneous optimization*/);
			//starting mover list
			for( utility::vector1< DesignMoverRealPair >::const_iterator it=stub_minimize_movers_.begin(); it!=stub_minimize_movers_.end(); ++it ){
				simple_moves::DesignRepackMoverOP const curr_mover( it->first );
				core::Real const bb_cst_weight( it->second );
				TR<<"applying mover: "<<curr_mover->get_name()<<'\n';
				//restricting movers for stub minimization
				curr_mover->prevent_repacking( prevent_repacking() );
				curr_mover->optimize_foldtree( false );
				curr_mover->design( false ); //we dont want any design to take place within any mover for stub minimization
				TR<<" design and repacking during stub minimization is prevented" << std::endl;
				TR<<"using weight: "<<bb_cst_weight<<" for the stub bb constraints "<<std::endl;
				ScoreFunctionOP minimize_mover_scorefxn_repack( curr_mover->scorefxn_repack() );
				if( minimize_mover_scorefxn_repack )
					minimize_mover_scorefxn_repack->set_weight( backbone_stub_constraint, bb_cst_weight );
				ScoreFunctionOP minimize_mover_scorefxn_minimize( curr_mover->scorefxn_minimize() );
				if( minimize_mover_scorefxn_minimize )
					minimize_mover_scorefxn_minimize->set_weight( backbone_stub_constraint, bb_cst_weight );
				curr_mover->apply( pose );
			} //end of user movers
			MinimizeInterface( pose, stub_scorefxn, no_min/*bb*/, sc_min, min_rb()/*rb*/, false /*optimize foldtree*/, no_targets, true/*simultaneous optimization*/ );
			utility::vector1< bool > min_host( pose.total_residue(), false );
			for( core::Size i=host_chain_begin; i<=host_chain_end; ++i ) min_host[ i ] = true;
			utility::vector1< bool > const no_min_rb( pose.num_jump(), false );
			MinimizeInterface( pose, scorefxn, min_host/*bb*/, sc_min, no_min_rb/*rb*/, false /*optimize foldtree*/, no_targets, true/*simultaneous optimization*/ );
		}//size()
		else{
			TR<<"Doing rb minimization towards the stub and sc minimization of placed stubs" <<std::endl;
			MinimizeInterface( pose, stub_scorefxn, no_min/*bb*/, sc_min, min_rb()/*rb*/, false /*optimize foldtree*/, no_targets, true/*simultaneous optimization*/ );
		}
		//Reporting statistics for single-stub minimization
		//getting angle alignments, CB distances and CB force
		core::Real C_N_angle_after_minmize_mover;
		core::Real CB_CA_angle_after_minmize_mover;
		test_res_res_aln(  pose.residue( host_residue ), *(stub->residue()), C_N_angle_after_minmize_mover, CB_CA_angle_after_minmize_mover );
		//get distances and cb-force -- this is not the applied one, this is just for stats
		core::Real const after_min_distance( pose.residue( host_residue ).xyz( "CB" ).distance( stub->residue()->xyz( "CB" ) ) );
		Real const pi(numeric::NumericTraits<Real>::pi());
		core::Real const curr_bonus( stub->bonus_value() );
		core::Real const after_min_cb_force( -curr_bonus / (after_min_distance * after_min_distance) );
		stats_TR<<"CB distances, before minimize movers: "<<before_min_distance<<" and after: "<<after_min_distance<<'\n';
		stats_TR<<"C_N angle alignment, before minimization: "<<( C_N_angle_before_minmize_mover *180 )/pi<<" and after: "<<( C_N_angle_after_minmize_mover *180 )/pi<<'\n';
		stats_TR<<"CB_CA angle alignment, before minimization: "<<( CB_CA_angle_before_minmize_mover *180 )/pi<<" and after: "<<( CB_CA_angle_after_minmize_mover *180 )/pi<<'\n';
		stats_TR<<"cb force, before minimization: "<<"before entering minimize movers: "<<before_min_cb_force <<" and after minimization: "<<after_min_cb_force<<'\n';
		//quality checks
		core::Size const after_apply_stub_minimize_mover( pose.constraint_set()->get_all_constraints().size() );
		TR_debug<<"before applying stub minimization mover "<<before_apply_stub_minimize_mover<<" constraints. After: "<<after_apply_stub_minimize_mover<<std::endl;
		if( before_apply_stub_minimize_mover != after_apply_stub_minimize_mover ){
			TR<<" ***ERROR: stub minimization has changed the number of constraints on the pose. Before: "<< before_apply_stub_minimize_mover<<" after: "<<after_apply_stub_minimize_mover<<". This behaviour is unsupported."<<std::endl;
			runtime_assert( before_apply_stub_minimize_mover != after_apply_stub_minimize_mover );
		}
	}//stub!=NULL
	else{
		TR<<"Doing rb minimization towards the stub and sc minimization of placed stubs" <<std::endl;
		MinimizeInterface( pose, stub_scorefxn, no_min/*bb*/, sc_min, min_rb()/*rb*/, false /*optimize foldtree*/, no_targets, true/*simultaneous optimization*/ );
	}

	////now completly done with minimization, checking that nothign unwanted has changed and clean up
	core::Real const after_min( stf.compute( pose ) );
	TR<<"backbone_stub_constraint score changed from "<<before_min<<" to "<<after_min<<" during minimization\n";
	TR<<"removing hotspot constraints from pose" <<std::endl;
	protocols::hotspot_hashing::remove_hotspot_constraints_from_pose( pose );

	ConstraintCOPs const csts_after_min = pose.constraint_set()->get_all_constraints();
	core::Size const num_csts_after_min( csts_after_min.size() );
	TR_debug<<"Csts before min "<<num_csts_before_min<<" after min "<<num_csts_after_min<<std::endl;
	runtime_assert( num_csts_after_min == num_csts_before_min );//constraints shouldnt have changed

	TR<<"removed all hotspot constraints\n";

	if( after_min >= -0.0001 ){
		TR<<"bb stub constraint score too high. Skipping this stub\n";

		return false;
	}
	return true;
}

/// @details selects stubs by iterating over the stub_set_. Returns the status of stub selection
bool
PlaceStubMover::SelectStubIteratively( protocols::hotspot_hashing::HotspotStubSet::Hs_vec::const_iterator stub_it ) {
	using namespace protocols::hotspot_hashing;

	bool accepted;
	do {
		if( stub_it == stub_set_->end() )
			return false;
		HotspotStubOP stub = stub_it->second.second;
		core::Real const bonus( stub->bonus_value() );
		TR<<"trying stub of bonus "<<bonus<<'\n';
		accepted = bonus <= score_threshold_;
		if( !accepted ){
			TR<<"failed this stub\n";
			++stub_it;
		}
	} while( !accepted );
	return true;
}

/// @details remove all coordinate constraints from pose and then reapply them, changing
/// the HarmonicFunc's sdev to a new value
/// Nov09 Changing of previous logic. Now, only changing the sdev associated with the
/// coordinate constraints.
void
PlaceStubMover::refresh_coordinate_constraints( core::pose::Pose &, core::Real const sdev )
{
//	remove_coordinate_constraints_from_pose( pose );
	coord_cst_func_->sd( sdev );
//	previous_coordinate_constraints_ = pose.add_constraints( previous_coordinate_constraints_ );
//	curr_coordinate_constraints_ = pose.add_constraints( curr_coordinate_constraints_ );
}

/// @details Placing a stub in the context of a complex by putting a stub on top of
///  the scaffold. The following steps are taken:
///0. minimization of rb dofs in a bb_stub_constraint-dominated force field. Only the
///   constraints implied by the stubset associated with this placestub mover will be
///   applied.
///1. stub selection, by iterating over a randomly shuffled stubset.
///2. finding a residue on the host chain where the stub may be placed (distance cutoff).
///3. testing the stub's repulsive energy on that residue in the context of a poly-alanine host chain.
///4. repack/minimize the rb orientation of the pose (plus a little bb minimization around the stub, but no repacking of host_chain_).
///5. testing the stub's total_score in the context of the full pose.
///6. a vector of user-defined movers derived from DesignRepackMover. These may
///   include further placestub movers. If this is the case, these movers can be constrained
///   with coordinate constraints that are placed on the functional groups of the stub
///   sidechain. These movers can also define the force constant of these harmonic constraints.
///   Notice that this design should in
///   principle be fully extensible to as many stub placement and design movers as the
///   user chooses. If any of the user defined design movers (including subsequent stub
///   placement) fail, another stub placement is attempted and the user defined movers are
///   called again.
///7. user_defined final_filter_ is called. If the pose does not pass this, place stub return to iterate over another stub. If all stubs fail, the mover signals failure and exits.
void
PlaceStubMover::apply( core::pose::Pose & pose )
{
	using namespace protocols::hotspot_hashing;

	residue_level_tasks_for_placed_hotspots_->clear();
	curr_coordinate_constraints_.clear();
	// rescore the pose to ensure that backbone stub cst's are populated
	core::scoring::ScoreFunctionOP bbcst_scorefxn( new core::scoring::ScoreFunction );
	bbcst_scorefxn->reset();
	bbcst_scorefxn->set_weight( core::scoring::backbone_stub_constraint, 1.0 );
	(*bbcst_scorefxn)( pose );

	// Remove old hotspot constraints from pose
	core::Size count_two_sided( 0 );
	core::Size const before_removal( pose.constraint_set()->get_all_constraints().size() );
	saved_bb_constraints_ = protocols::hotspot_hashing::remove_hotspot_constraints_from_pose( pose );

	core::Size const after_removal( pose.constraint_set()->get_all_constraints().size() );
	core::Size const bb_constraints_size( saved_bb_constraints_.size() );
	if( bb_constraints_size > 0 ){
		TR_debug<<"removed "<<bb_constraints_size<<" bb constraints from pose\n";
		TR_debug<<"before removal "<<before_removal<<" after removal "<<after_removal<<'\n';
	}

	core::Size const host_chain_begin( pose.conformation().chain_begin( host_chain_ ));
	core::Size const host_chain_end( pose.conformation().chain_end( host_chain_ ));
	if( !default_fold_tree_ )
		default_fold_tree_ = core::kinematics::FoldTreeOP( new kinematics::FoldTree( pose.fold_tree() ) );

	// change all interface residues on host_chain_ to ala, except for those that
	// might already contain stubs or that are prevented from repacking for some other reason.
	{// change to ala pose
		// Switch to Ala unless we are doing place scaffold as a replacement for docking
		if( !(place_scaffold_ && !triage_positions_) ){
			BuildAlaPose toAla( host_chain_ == 1/*partner1*/, host_chain_ == 2 /*partner2*/ );
			utility::vector1< core::Size > no_repack;
			if( !prevent_repacking().empty() ) no_repack = prevent_repacking();
			if( !placed_stubs_.empty() ){
				for( utility::vector1< std::pair< core::Size, bool > >::const_iterator it=placed_stubs_.begin(), end=placed_stubs_.end(); it!=end; ++it )
				no_repack.push_back( it->first );
			}
			if( !no_repack.empty() ){
				std::sort( no_repack.begin(), no_repack.end() );
				std::unique( no_repack.begin(), no_repack.end() );
				toAla.prevent_repacking( no_repack );
			}

			TR<<"switching interface to alanine\n";
			toAla.apply( pose );
			bool const minimization_status( StubMinimize( pose, NULL/*HotspotStubCOP*/, 0/*host_res*/, hurry_  ) );

			// Optional Triaging
			if( triage_positions_ ) {
				if( !minimization_status ){ // requires partners to already be placed correctly
					TR<<"constraint score evaluates to 0 in this pose. Failing\n";
					set_last_move_status( protocols::moves::FAIL_RETRY );
					final_cleanup( pose );
					return;
				}
			}
			// End Triaging
		}
	}//change to ala pose

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	core::scoring::ScoreFunctionCOP soft_rep = core::scoring::ScoreFunctionFactory::create_score_function( SOFT_REP_DESIGN_WTS );

	core::pose::Pose const saved_pose( pose ); //saved after minimization
	saved_prevent_repacking_ = prevent_repacking();
	saved_placed_stubs_ = placed_stubs_;
	numeric::random::random_permutation( stub_set_->begin(), stub_set_->end(), numeric::random::rg() );// randomly shuffling stubs so that the selection doesn't repeat the same stub order each time
	HotspotStubSet::Hs_vec::iterator stub_it( stub_set_->begin() );
	std::vector< core::Size > host_positions;
	if( task_factory() ){
		utility::vector1< core::Size > const designable( protocols::rosetta_scripts::residue_packer_states( pose, task_factory(), true/*designable*/, false/*packable*/ ) );
		BOOST_FOREACH( core::Size const d, designable )
			host_positions.push_back( d );
	}
	else{
		for( core::Size i=host_chain_begin+1; i<=host_chain_end-1; ++i ) // avoid placements on chain termini
			host_positions.push_back( i );
	}

	do { // for each two sided trial
		//select a stub
		if( !SelectStubIteratively( stub_it ) ){
			TR<<"Stub selection failed" <<std::endl;
			set_last_move_status( protocols::moves::FAIL_RETRY );
			final_cleanup( pose );
			return;
		}
		HotspotStubOP stub( stub_it->second.second );
		++stub_it;

		// next, place the stub on a position in the host chain that is close enough
		// randomize the order of placement on host residues to encourage diversity
		numeric::random::random_permutation( host_positions.begin(), host_positions.end(), numeric::random::rg() );
		for( std::vector< core::Size >::const_iterator host_pos_it = host_positions.begin(); host_pos_it!=host_positions.end(); ++host_pos_it) {
			using namespace core::conformation;
			core::Size const res( *host_pos_it );
			Residue const res_host( pose.residue( res ) );
			Residue const res_stub( *stub->residue() );

			// Obligatory Triaging
			if( res_host.aa() != chemical::aa_ala ) continue; // position already designed
			if( res_host.aa() == chemical::aa_gly
					|| ( res_host.aa() == chemical::aa_pro && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline] )
					|| res_host.type().is_disulfide_bonded() )
				continue; // disallowed host amino acids
			if( !prevent_repacking().empty()
					&& std::find( prevent_repacking().begin(), prevent_repacking().end(), res ) != prevent_repacking().end() )
				continue; //not allowed to repack
			if( !disallowed_host_pos_.empty()
					&& std::find( disallowed_host_pos_.begin(), disallowed_host_pos_.end(), res ) != disallowed_host_pos_.end() )
				continue; //not allowed to be a stub

			core::Real const distance( res_host.xyz( "CB" ).distance( res_stub.xyz( "CB" ) ) );

			// Optional Triaging of residues too far away
			if( triage_positions_ ) {
				// triage based on distance
				if( distance >= max_cb_cb_dist_ ) continue; // too far to be productive
				core::Real dummy1, dummy2;
				if( !test_res_res_aln( res_host, *stub->residue(), dummy1, dummy2 ) ){
					TR_debug<<"stub failed placement due to alignment"<<std::endl;
					continue; // stub not aligned properly to position
				}
			}
			// End Triaging

			TR << "Trying host position " << res << std::endl;

			protocols::filters::FilterCOP rep_filter( protocols::filters::FilterOP( new simple_filters::EnergyPerResidueFilter(
				res, soft_rep, core::scoring::fa_rep, 20.0 ) ) );
			stub_set_->filter( rep_filter );
			stub->set_filter( rep_filter );
			bool const self_energy_pass( stub->get_scaffold_status( res ) );
			if( self_energy_pass ){

				// the following two lines make sure that t//he pose is at the
				// starting line before each stub placement. These lines belong
				// wherever the pose first changes after minimization.
				pose=saved_pose; // to make sure that we've come back to the beginning
				cst_cleanup( pose );

				++count_two_sided;
				place_stub( pose, res_stub, res );
				{ // build ala pose scope
					BuildAlaPose toAla( host_chain_ == 1/*partner1*/, host_chain_ == 2 /*partner2*/ );
					utility::vector1< core::Size > no_repack;

					if( !prevent_repacking().empty() )
						no_repack = prevent_repacking();

					for( utility::vector1< std::pair< core::Size, bool > >::const_iterator
							it=placed_stubs_.begin(), end=placed_stubs_.end();
							it!=end; ++it ) {
						no_repack.push_back( it->first );
					}
					//place_stub should have inserted the stub into placed_stubs
					assert( std::find(no_repack.begin(), no_repack.end(), res) != no_repack.end() );

					std::sort( no_repack.begin(), no_repack.end() );
					std::unique( no_repack.begin(), no_repack.end() );
					toAla.prevent_repacking( no_repack );

					TR<<"switching interface to alanine\n";
					toAla.apply( pose );
				}

				stats_TR << "Stub of energy "<<stub->bonus_value()
					<<" and type "<<pose.residue( res ).name3()
					<<" placed on residue "<<res
					<<". The distance between stub's Cb and the residue's is "
					<<distance<<"A. Passing rep filter"<<std::endl;

				{ // minimize around hs
					// this will minimize the rb orientation (+ the stub's sc) with a bb constraint based solely on the relevant stub
					core::Size const before_min( pose.constraint_set()->get_all_constraints().size() );
					StubMinimize( pose, stub, res, hurry_ );
					core::Size const after_min( pose.constraint_set()->get_all_constraints().size() );
					TR_debug <<"DEBUG: before minimize around hs we had "<<before_min<<" constraints. after: "<<after_min<<std::endl;
					Residue const res_host_after_min( pose.residue( res ) );
					core::Real const dist_after_min_threshold( 2.0 );
					core::Real const distance_after_min( res_host_after_min.xyz( "CB" ).distance( res_stub.xyz( "CB" ) ) );
					stats_TR<<"stub - host residue Cb distance after stub-based minimization is "<<distance_after_min<<std::endl;
					if( distance_after_min >= dist_after_min_threshold ){
						stats_TR<<"Distance too large. failing stub\n";
						pose = saved_pose;
						cst_cleanup( pose );
						continue;
					}
					/// modify after-placement filter to contain the stub's residue information
					protocols::filters::FilterOP modified_after_placement_filter( after_placement_filter_->clone() );
					protocols::moves::modify_ResId_based_object( modified_after_placement_filter, res );
					bool const pass_after_placement( modified_after_placement_filter->apply( pose ) );
					if( pass_after_placement )
						stats_TR<<"Passed after placement filter"<<std::endl;
					else{
						stats_TR<<"failed after placement filter"<<std::endl;
						pose = saved_pose;
						cst_cleanup( pose );
						continue;
					}
				}//end minimization scope
				kinematics::FoldTree const stub_ft( pose.fold_tree() );
				core::Size const before_remove( pose.constraint_set()->get_all_constraints().size() );
//				remove_coordinate_constraints_from_pose( pose );
				core::Size const after_remove( pose.constraint_set()->get_all_constraints().size() );
				TR_debug <<"DEBUG: before remove coord constraints we had "<<before_remove<<" constraints. after: "<<after_remove<<std::endl;

				// this is to make sure that nothing awful has happened during stub-based minimization.
				protocols::filters::FilterCOP total_energy_filter( protocols::filters::FilterOP( new simple_filters::EnergyPerResidueFilter( res, scorefxn, core::scoring::total_score, stub_energy_threshold_ ) ) );
				bool const interface_energy_pass( total_energy_filter->apply( pose ) );
				if( interface_energy_pass ){
					//If a subsequent placement fails, try placing the stub on another position or
					//iterate over more stubs
					bool subsequent_stub_placement_failure( false );
					TR<<"redesigning remainder of interface with user defined design movers\n";
					prevent_repacking_.push_back( res );
					if( residue_numbers_ ){
						TR<<"Pushing residue number "<<res<<" to residue_numbers_ object"<<std::endl;
						residue_numbers_->obj.clear();
						residue_numbers_->obj.push_back( res ); // in principle should push the number to the stack, but that keeps memory of past placements. Need to find a better way to deal with that.
					}
					utility::vector1< core::Size > empty;
					runtime_assert( coord_cst_std_.size() == design_movers_.size() );
					utility::vector1< core::Real >::const_iterator it_sdev( coord_cst_std_.begin() );
					for( utility::vector1< DesignMoverFoldTreePair >::iterator it=design_movers_.begin(); it!=design_movers_.end() && !subsequent_stub_placement_failure; ++it ){
						TR<<"applying design mover "<<it->first->get_name()<<'\n';
						bool const use_constraints( it->second );
						if( use_constraints ){
							TR<<"using coordinate constraints with sdev "<<*it_sdev<<'\n';
							core::Size const before_refresh( pose.constraint_set()->get_all_constraints().size() );
							refresh_coordinate_constraints( pose, *it_sdev );
							core::Size const after_refresh( pose.constraint_set()->get_all_constraints().size() );
							TR_debug<<"before refreshing coord cst "<<before_refresh<<" constraints. after: "<<after_refresh<<std::endl;
						}//use constraints
						else
							TR<<"no constraints applied\n";
						core::Size const before_apply_design_mover( pose.constraint_set()->get_all_constraints().size() );
						it->first->prevent_repacking( prevent_repacking() );
						it->first->optimize_foldtree( false );
						if( it->first->get_name() == "PlaceStub" ){ //special treatment to place stub movers
							// We want the child place stub movers to contain all the information
							// obtained up to now, including, which residues harbour stubs, where
							// repacking is disallowed, the ft and what were the native sidechains. In some
							// cases, we might want the child placestub movers to remain pristine, so
							// we modify and apply clones of these movers.
							movers::PlaceStubMoverOP modified_place_stub( utility::pointer::dynamic_pointer_cast< protocols::protein_interface_design::movers::PlaceStubMover > ( it->first->clone() ) );
							runtime_assert( modified_place_stub != 0 );
							modified_place_stub->placed_stubs_ = placed_stubs_;
							utility::vector1< core::Size > new_prev_repack( prevent_repacking() );
							for( utility::vector1< std::pair< core::Size, bool > >::const_iterator it=placed_stubs_.begin(), end=placed_stubs_.end(); it!=end; ++it )
								new_prev_repack.push_back( it->first );

							modified_place_stub->prevent_repacking( new_prev_repack );
							modified_place_stub->default_fold_tree_ = default_fold_tree_;
							modified_place_stub->previous_coordinate_constraints_ = curr_coordinate_constraints_;
							modified_place_stub->coord_cst_func_ = coord_cst_func_;
							modified_place_stub->apply( pose );
							if( modified_place_stub->get_last_move_status() == protocols::moves::FAIL_RETRY ){
								stats_TR<<"Subsequent stub placement failed. Trying another stub placement"<<std::endl;
								subsequent_stub_placement_failure = true;
								break;
							}
						} //a place stub mover
						else{ // other design movers
							TR_debug<<"setting coordinate constraint weights to 1.0 in design movers"<<std::endl;
							core::scoring::ScoreFunctionOP scorefxn_rep( it->first->scorefxn_repack() );
							core::scoring::ScoreFunctionOP scorefxn_min( it->first->scorefxn_minimize() );
							core::Real const coord_cst_weight( use_constraints ? 1 : 0 );
							if( scorefxn_rep ) scorefxn_rep->set_weight( coordinate_constraint, coord_cst_weight );
							if( scorefxn_min ) scorefxn_min->set_weight( coordinate_constraint, coord_cst_weight );
							it->first->apply( pose );
						}
						++it_sdev;
						core::Size const after_apply_design_mover( pose.constraint_set()->get_all_constraints().size() );
						TR_debug<<"before applying design mover "<<before_apply_design_mover<<" constraints. After: "<<after_apply_design_mover<<std::endl;
						if( before_apply_design_mover != after_apply_design_mover ){
							TR<<"ERROR: This design mover changed the number of constraints on the pose. Before: "<< before_apply_design_mover<<" after: "<<after_apply_design_mover<<". This behaviour is unsupported."<<std::endl;
							runtime_assert( before_apply_design_mover != after_apply_design_mover );
						}
//						TR<<"removing coordinate constraints\n";
//						remove_coordinate_constraints_from_pose( pose );
					}//Design movers
					if( subsequent_stub_placement_failure ){
						cst_cleanup( pose );
						break;
					}
					core::Size const before_final_filter( pose.constraint_set()->get_all_constraints().size() );
					TR_debug<<"before final filter "<<before_final_filter<<" constraints"<<std::endl;
					/// modify final filter to contain the stub's residue information
					protocols::filters::FilterOP modified_final_filter( final_filter_->clone() );
					protocols::moves::modify_ResId_based_object( modified_final_filter, res );
					bool const final_filter_pass( modified_final_filter->apply( pose ) );
					if( final_filter_pass ) {
						stats_TR<<"SUCCESS: Stub passed final filter and placed at position "<<res<<". after "<<count_two_sided<<" two sided evaluations."<< std::endl;
						if( residue_level_tasks_for_placed_hotspots_ ){
							TR_debug<<"Adding taskawareness to taskaware design movers\n";
							core::Size const size_before_adding( residue_level_tasks_for_placed_hotspots_->size() );
							for( utility::vector1< std::pair< core::Size, bool > >::const_iterator placed_stub( placed_stubs_.begin()); placed_stub!=placed_stubs_.end(); ++placed_stub ){
								using namespace core::pack::task::operation;
								PreventRepackingOP pr( new PreventRepacking );
								pr->include_residue( placed_stub->first );
								residue_level_tasks_for_placed_hotspots_->push_back( pr );
							}//for placed_stubs_
							core::Size const size_after_adding( residue_level_tasks_for_placed_hotspots_->size() );
							TR_debug<<"Before adding we had "<<size_before_adding<<" task operations. Now we have "<<size_after_adding<<std::endl;
						}//fi residue_level_tasks
						set_last_move_status( protocols::moves::MS_SUCCESS );
						final_cleanup( pose );
						core::Size const after_success( pose.constraint_set()->get_all_constraints().size() );
						TR_debug<<"after success "<<after_success<<" constraints"<<std::endl;
						TR.flush();

						// Experimental: add residue position of successful stub placement to the score.sc file
						protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
						//std::string column_header = this->get_user_defined_name();

						//convert residue number into a string
						std::ostringstream convert;
						convert << res;

						// set column name
						//std::string column_header = "hotspot_" + convert.str();
						std::string column_header = user_defined_name_;

						//job->add_string_real_pair(column_header, res);
						job->add_string_string_pair(column_header, convert.str());

						return;
					}//final_filter_pass
					else{
						pose = saved_pose;
						cst_cleanup( pose );
						stats_TR<<"Stub at position "<<res<<" failed final filter"<<std::endl;
					}
				}//interface energy filter
				else
					stats_TR<<"Stub at position "<<res<<" failed interface energy filter"<<std::endl;
			}//self_energy pass
			else
				stats_TR<<"Stub failed self-energy filter"<<std::endl;
			pose = saved_pose; // revert to beginning
			cst_cleanup( pose );
			TR.flush();
		} //loop over host chain residues
	} while( count_two_sided <= stub_set_->size() ); // this condition will never be violated here b/c SelectStubIteratively above has the same condition. Nevertheless, better have this check to ensure we don't launch an infinite loop
	stats_TR<<"Stub placement failed after "<<count_two_sided<<" two-sided_trials."<<std::endl;
	TR.flush();
	pose = saved_pose;
	final_cleanup( pose );
	set_last_move_status( protocols::moves::FAIL_RETRY );
}

std::string
PlaceStubMover::get_name() const {
	return PlaceStubMoverCreator::mover_name();
}

/// @details reapply saved coord constraints and refresh placed_stubs and prevent_repacking
/// Nov09 changed logic: curr_coordinate constraints are removed only
void
PlaceStubMover::cst_cleanup( core::pose::Pose & pose )
{
	core::Size const crd_cst_size( curr_coordinate_constraints_.size() );
	if( crd_cst_size > 0 ){
		TR<<"removed "<<crd_cst_size<<" coordinate constraints"<<std::endl;
		pose.remove_constraints( curr_coordinate_constraints_ );
//		remove_coordinate_constraints_from_pose( pose );
		curr_coordinate_constraints_.clear();
//		previous_coordinate_constraints_ = pose.add_constraints( previous_coordinate_constraints_ );
	}
	placed_stubs_ = saved_placed_stubs_;
	prevent_repacking_ = saved_prevent_repacking_;
	pose.update_residue_neighbors();
}

/// @details This should be called before exiting placestub.
void
PlaceStubMover::final_cleanup( core::pose::Pose & pose )
{
	if( get_last_move_status() == protocols::moves::MS_SUCCESS &&
			leave_coord_csts_after_placement_ ){
		runtime_assert( post_placement_sdev_ >= 0.0 );
		TR<<"Cleaning up after successful stubplacement, but NOT removing coordinate constraints\n";
		refresh_coordinate_constraints( pose, post_placement_sdev_ );
	}
	else
		cst_cleanup( pose );
	if( !saved_bb_constraints_.empty() ){
		core::Size const bb_constraints_size( saved_bb_constraints_.size() );
		if( bb_constraints_size > 0 ){
			saved_bb_constraints_ = pose.add_constraints( saved_bb_constraints_ );
			TR<<"reapplied "<<bb_constraints_size<<" bb constraints to pose"<<std::endl;
		}
	}

	//don't accumulate state between ntrials
	placed_stubs_.clear();
	prevent_repacking_.clear();
	restrict_to_repacking_.clear();
	previous_coordinate_constraints_.clear();
	curr_coordinate_constraints_.clear();
// commenting out fold tree reseting b/c downstream movers would want to use it
//	pose.fold_tree( *default_fold_tree_ );//this is important for allowing further centroid-based moves/filters on the pose
}

/// @details changes the pose's fold-tree to connect
/// the nearest residue on the target to the stub on the scaffold. The foldtree connects the last carbon atom before the
/// stub's functional group. That can be useful in minimization of the pose, b/c the stub's interaction with the target will
/// not be lost due to minimization. Notice that the pose's fold tree changes in the process and so it is usually a good idea to
/// save the old foldtree and reinstate it after minimization.
/// Can work with multiple stubs. In that case, rb jumps are introduced between the target and each of the stubs
/// and a cutpoint is introduced in the putative circle that has just been formed in the fold tree.
/// Alternatively, if the current stub is to be added via a constraint, then a coordinate
/// constraint is set up for this stub
void
PlaceStubMover::stub_based_atom_tree( core::pose::Pose & pose, core::conformation::Residue const res_stub, core::Real const cst_sdev )
{
//	protocols::loops::remove_cutpoint_variants( pose, true/*force*/ );
//	core::Size const before_removing_coord_cst( pose.constraint_set()->get_all_constraints().size() );
//	remove_coordinate_constraints_from_pose( pose );
//	core::Size const after_removing_coord_cst( pose.constraint_set()->get_all_constraints().size() );
//	TR_debug<<"before removing coordinate cst "<<before_removing_coord_cst<<" constraints. after "<<after_removing_coord_cst<<std::endl;

/// Add constraints pertaining to the last placed stub
	core::Size const curr_stub_resnum( placed_stubs_[ placed_stubs_.size() ].first );
	bool const add_constraints( placed_stubs_[ placed_stubs_.size() ].second );
	if( add_constraints )
		curr_coordinate_constraints_ = add_coordinate_constraints( pose, res_stub, host_chain_, curr_stub_resnum, cst_sdev, coord_cst_func_ );
/// Set atom tree according to the first placed stub.
	protocols::protein_interface_design::movers::SetAtomTree sat;
	core::Size const host_residue_num( placed_stubs_.begin()->first );
	pose.fold_tree( *sat.create_atom_tree( pose, host_chain_, host_residue_num ) );
	TR<<"Stub based atom tree: "<<pose.fold_tree()<<std::endl;
}

void
PlaceStubMover::parse_my_tag( TagCOP const tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		Movers_map const &movers,
		core::pose::Pose const & pose )
{
	using namespace protocols::hotspot_hashing;

	TR<<"Parsing PlaceStubMover----"<<std::endl;

  if( tag->hasOption( "name" ) ){
		// Experimental
		user_defined_name_ = tag->getOption< std::string >( "name" );
	}

	if( tag->hasOption( "residue_numbers_setter" ) ){
		std::string const residue_numbers_name( tag->getOption( "residue_numbers", tag->getOption< std::string >( "residue_numbers_setter" ) ) );
		residue_numbers_ = basic::datacache::get_set_from_datamap< basic::datacache::DataMapObj< utility::vector1< core::Size > > >( "residue_numbers", residue_numbers_name, data );
	}

	if( tag->hasOption( "task_operations" ) ){
		if( task_factory() )
			TR<<"*****WARNING: ERASING existing task_factory, b/c of specifications for new task operations in\n"<<tag<<std::endl;
		task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	}
	host_chain_ = tag->getOption<core::Size>( "chain_to_design", 2 );
	score_threshold_ = tag->getOption<core::Real>( "score_threshold", 0 );
	hurry_ = tag->getOption<bool>( "hurry", 0 );
	triage_positions_ = tag-> getOption<bool>( "triage_positions", 1 );
	stub_energy_threshold_ = tag-> getOption<core::Real>( "stub_energy_threshold", 1.0 );

	if( tag->hasOption( "stubfile" ) ) { //assigning a unique stubset
		std::string const stub_fname = tag->getOption<std::string>( "stubfile" ); // Experimental: stub file name
		stub_set_ = protocols::hotspot_hashing::HotspotStubSetOP( new protocols::hotspot_hashing::HotspotStubSet );
		if( data.has( "hotspot_library", stub_fname ) ){
			stub_set_ = data.get_ptr<protocols::hotspot_hashing::HotspotStubSet>( "hotspot_library", stub_fname );
			TR<<"Associated mover with an already read stubset named "<<stub_fname<<std::endl;
		}
		else
			stub_set_->read_data( stub_fname );
	}
	else { //assigning the stubset taken from the constraints setup
		std::string const hs( "hotspot_stubset" );
		stub_set_ = data.get_ptr<HotspotStubSet>( "constraints", hs );
	}

	bool const minimize_rb( tag->getOption< bool >( "minimize_rb", 0 ) );
	min_rb( minimize_rb );

	std::string const after_placement_filter_name( tag->getOption<std::string>( "after_placement_filter", "true_filter" ) );
	protocols::filters::Filters_map::const_iterator find_ap_filter( filters.find( after_placement_filter_name ));

	bool const ap_filter_found( find_ap_filter != filters.end() );
	if( ap_filter_found )
		after_placement_filter_ = find_ap_filter->second->clone();
	else {
		if( after_placement_filter_name != "true_filter" ){
			TR<<"***WARNING WARNING! Filter defined for PlaceStubMover not found in filter_list!!!! Defaulting to truefilter***"<<std::endl;
			runtime_assert( ap_filter_found );
		}
		else
			after_placement_filter_ = protocols::filters::FilterOP( new protocols::filters::TrueFilter );
	}

	std::string const final_filter_name( tag->getOption<std::string>( "final_filter", "true_filter" ) );
	protocols::filters::Filters_map::const_iterator find_filter( filters.find( final_filter_name ));

	bool const filter_found( find_filter != filters.end() );
	if( filter_found )
		final_filter_ = find_filter->second->clone();
	else {
		if( final_filter_name != "true_filter" ){
			TR<<"***WARNING WARNING! Filter defined for PlaceStubMover not found in filter_list!!!! Defaulting to truefilter***"<<std::endl;
			runtime_assert( filter_found );
		}
		else
			final_filter_ = protocols::filters::FilterOP( new protocols::filters::TrueFilter );
	}

	//parse allowed residues
	disallowed_host_pos_.clear();
	if( tag->hasOption("allowed_host_res") ) {
		utility::vector1<Size> allowed_host_pos(
			core::pose::get_resnum_list(
				tag, "allowed_host_res", pose ));
		//disallowed is the set complement of allowed

		core::Size const host_chain_begin( pose.conformation().chain_begin( host_chain_ ));
		core::Size const host_chain_end( pose.conformation().chain_end( host_chain_ ));

		for(Size host_pos( host_chain_begin); host_pos <= host_chain_end; ++host_pos )
		{
			if( std::find( allowed_host_pos.begin(), allowed_host_pos.end(), host_pos)
					== allowed_host_pos.end() )
			{
				// not allowed
				disallowed_host_pos_.push_back(host_pos);
			}
		}
		TR.Debug<<"Disallowing "<<disallowed_host_pos_.size()
			<<" of "<<host_chain_end - host_chain_begin + 1
			<<" host positions."<<std::endl;
	}

	//parsing stub minimize movers and design movers for place stub
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );
	for( utility::vector0< TagCOP >::const_iterator btag=branch_tags.begin(); btag!=branch_tags.end(); ++btag ) {
		if( (*btag)->getName() == "StubMinimize" ){
			utility::vector0< TagCOP > const & stub_min_tags( (*btag)->getTags() );
			for( utility::vector0< TagCOP >::const_iterator stub_m_tag=stub_min_tags.begin(); stub_m_tag!=stub_min_tags.end(); ++stub_m_tag ) {
				std::string const stub_mover_name( (*stub_m_tag)->getOption<std::string>( "mover_name" ) );
				core::Real  const bb_stub_constraint_weight( (*stub_m_tag)->getOption< core::Real > ( "bb_cst_weight", 10.0 ) );
				std::map< std::string const, MoverOP >::const_iterator find_mover( movers.find( stub_mover_name ));
				bool const stub_mover_found( find_mover != movers.end() );
				if( stub_mover_found ){
					simple_moves::DesignRepackMoverOP drSOP = utility::pointer::dynamic_pointer_cast< simple_moves::DesignRepackMover > ( find_mover->second->clone() );
					if( !drSOP ){
						TR<<"dynamic cast failed in tag "<<tag<<". Make sure that the mover is derived from DesignRepackMover"<<std::endl;
						runtime_assert( drSOP != 0 );
					}//done cast check
					stub_minimize_movers_.push_back( std::make_pair( drSOP, bb_stub_constraint_weight) );
					TR<<"added stub minimize mover "<<stub_mover_name<<" to minimize towards the stub. Using this weight for the bb stub constraints: "<< bb_stub_constraint_weight<<'\n';
				}
			}
		}
		else if( (*btag)->getName() == "DesignMovers" ){
			utility::vector0< TagCOP > const & design_tags( (*btag)->getTags() );
			for( utility::vector0< TagCOP >::const_iterator m_it=design_tags.begin(); m_it!=design_tags.end(); ++m_it ) {
				TagCOP const m_tag_ptr = *m_it;
				std::string const mover_name( m_tag_ptr->getOption< std::string >( "mover_name" ) );
				bool const apply_coord_constraints( m_tag_ptr->getOption< bool >( "use_constraints", 1 ) );
				core::Real const coord_cst_std( m_tag_ptr->getOption< core::Real >( "coord_cst_std", 0.5 ) );

				std::map< std::string const, MoverOP >::const_iterator find_mover( movers.find( mover_name ));
				bool const mover_found( find_mover != movers.end() );
				if( mover_found ){
					simple_moves::DesignRepackMoverOP drOP = utility::pointer::dynamic_pointer_cast< simple_moves::DesignRepackMover > ( find_mover->second->clone() );
					if( !drOP ){
						TR<<"dynamic cast failed in tag "<<tag<<". Make sure that the mover is derived from DesignRepackMover"<<std::endl;
						runtime_assert( drOP != 0 );
					}
					design_movers_.push_back( std::make_pair( drOP, apply_coord_constraints ) );
					coord_cst_std_.push_back( coord_cst_std );
					TR<<"added design mover "<<mover_name<<" to place stub with apply_coord_constraints switched to "<< apply_coord_constraints<<" with std "<< coord_cst_std<< " and hurry=" << hurry_ << '\n';
				}
				else{
					TR<<"***WARNING WARNING! Mover defined for PlaceStubMover not found in mover_list. EXITING ***"<<std::endl;
					runtime_assert( mover_found );
				}
			}
		}
		else if( (*btag)->getName() != "NotifyMovers" )
			utility_exit_with_message( "ERROR: tag in PlaceStub not defined\n" );
		generate_taskfactory_and_add_task_awareness( *btag, movers, data, residue_level_tasks_for_placed_hotspots_ );
	}
	if( stub_minimize_movers_.size() == 0 )
		TR<<"No StubMinimize movers defined by user, defaulting to minimize_rb and _sc of stubs only\n";

	add_constraints_ = tag->getOption< bool >( "add_constraints", 1 );


	{ // pair stubset with an ala_pose scope
		//we want the stubset to be aware of the host chain. This way, when we
		//ask a stub set whether a stub can be placed on a particular position
		//in terms of its self-energy, it can cache this information and return
		//it to us at a later point if we ask for that stub again. Since we
		//expect the host_chain to be redesignable, we switch all residues
		//(other than gly/pro) to ala before pairing.
		core::pose::PoseOP ala_pose( new core::pose::Pose( pose ) );
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( *ala_pose ));
		task->initialize_from_command_line().or_include_current( true );

		utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, false );
		allowed_aas[ chemical::aa_ala ] = true;

		core::Size const chain_begin( ala_pose->conformation().chain_begin( host_chain_ ) );
		core::Size const chain_end( ala_pose->conformation().chain_end( host_chain_ ) );

		for ( core::Size i = 1; i <= pose.total_residue(); i++) {
			if ( !pose.residue(i).is_protein() ) continue;
			if( i >= chain_begin && i <=chain_end ) {
				core::Size const restype( ala_pose->residue(i).aa() );
				if( (restype == chemical::aa_pro && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline] ) || restype == chemical::aa_gly )
					task->nonconst_residue_task(i).prevent_repacking();
				else
					task->nonconst_residue_task(i).restrict_absent_canonical_aas( allowed_aas );
			}
			else {
				task->nonconst_residue_task( i ).prevent_repacking();
			}
		}
		if( basic::options::option[basic::options::OptionKeys::packing::resfile].user() )
			core::pack::task::parse_resfile(pose, *task);

		core::scoring::ScoreFunctionOP scorefxn( get_score_function() );
		pack::pack_rotamers( *ala_pose, *scorefxn, task);
		(*scorefxn)( *ala_pose );
		stub_set_->pair_with_scaffold( *ala_pose, host_chain_, protocols::filters::FilterCOP( protocols::filters::FilterOP( new protocols::filters::TrueFilter ) ) );
	}// pair stubset scope

	max_cb_cb_dist_ = tag->getOption< core::Real >( "max_cb_dist", 4.0 );
	leave_coord_csts_after_placement_ = tag->getOption< bool >( "leave_coord_csts", 0 );
	if( leave_coord_csts_after_placement_ ){
		post_placement_sdev_ = tag->getOption< core::Real >("post_placement_sdev", 1.0 );
		TR<<"leaving constraints on after successful placement\n";
	}
	else
		post_placement_sdev_ = -1.0; // this will assert later on
	TR<<"max cb cb distance set to "<<max_cb_cb_dist_<<'\n';
	place_scaffold_ = tag->getOption< bool >( "place_scaffold", 0 );
	TR<<"place stub mover on chain "<<host_chain_<<" with score threshold of "<<score_threshold_<<" minimize_rb to "<<minimize_rb<<" final filter "<<final_filter_name<<" and place scaffold="<<place_scaffold_<<std::endl;
}

} //movers
} //protein_interface_design
} //protocols

