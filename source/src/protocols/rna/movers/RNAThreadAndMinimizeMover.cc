// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/rna/movers/RNAThreadAndMinimizeMover.cc
/// @brief Thread a new sequence over a given RNA scaffold and do a little optimization
/// @author Andy Watkins (amw579@nyu.edu)

// Unit headers
#include <protocols/rna/movers/RNAThreadAndMinimizeMover.hh>
#include <protocols/rna/movers/RNAThreadAndMinimizeMoverCreator.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>


// Core headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/annotated_sequence.hh>
//#include <core/import_pose/import_pose.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/io/NomenclatureManager.hh>
#include <core/pose/rna/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Protocols
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/stepwise/sampler/rna/RNA_KIC_Sampler.hh>


// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <numeric/random/random.hh>


#include <iostream>
#include <fstream>

// TODO:
// Should deprecate the non-long strategy
// stress test mutation list
// write iupac sequence support
// loop modeling should work eventually


using namespace core;

static basic::Tracer TR( "protocols.rna.movers.RNAThreadAndMinimizeMover" );

namespace protocols {
namespace rna {
namespace movers {

RNAThreadAndMinimizeMover::RNAThreadAndMinimizeMover( std::string const & seq,
	std::string const & template_sequence_from_alignment,
	bool long_strategy/* = false*/,
	std::string const & input_sequence_type/* = ""*/,
	std::string const & mutation_list/* = ""*/,
	std::string const & insertion_list/* = ""*/
):
	seq_( seq ),
	template_sequence_from_alignment_( template_sequence_from_alignment ),
	long_strategy_( long_strategy ),
	input_sequence_type_( input_sequence_type ),
	scorefxn_( core::scoring::get_score_function() ),
	target_sequence_( "" ),
	working_res_( utility::vector1< Size >() ),
	mutation_list_( mutation_list ),
	insertion_list_( insertion_list )
{}

void
RNAThreadAndMinimizeMover::obtain_rtypes_for_target_sequence() {
	core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	rtypes_ = core::pose::residue_types_from_sequence( target_sequence_, *rsd_set, true /*auto_termini*/ );

	// Go through mutlist and modify as needed
	if ( mutation_list_.size() == 0 ) return;

	auto const mutation_list = utility::string_split( mutation_list_, ' ' );
	for ( auto const & mutation : mutation_list ) {
		// expect: seqpos,identity
		std::istringstream iss( mutation );
		std::string seqpos_str, new_identity;
		std::getline( iss, seqpos_str, ',');
		Size const seqpos = stoi( seqpos_str );
		iss >> new_identity;//std::getline( iss, identity, ',')
		if ( input_sequence_type_ == "MODOMICS" ) {
			new_identity = core::io::NomenclatureManager::get_instance()
				->annotated_sequence_from_modomics_oneletter_sequence( new_identity );
		} else if ( input_sequence_type_ == "IUPAC" ) {
			new_identity = core::io::NomenclatureManager::get_instance()
				->annotated_sequence_from_IUPAC_sequence( new_identity );
		}
		rtypes_[ seqpos ] = core::pose::residue_types_from_sequence( new_identity, *rsd_set, false /*auto_termini*/ )[ 1 ];
		//change_sequence_accordingly( target_sequence_from_alignment, seqpos, new_identity );
	}
}

void RNAThreadAndMinimizeMover::set_up_target_sequence( core::pose::Pose & pose ) {

	using namespace core::pose;

	Size const nres = pose.total_residue();
	Pose const template_pose = pose;

	utility::vector1< Size > residues_missing_from_intended_sequence;

	// AMW TODO: consider processing mutlist first if provided
	// so that you can "mutlist" to a gap

	// AMW TODO:
	if ( input_sequence_type_ == "MODOMICS" ) {
		seq_ = core::io::NomenclatureManager::annotated_sequence_from_modomics_oneletter_sequence( seq_ );
	} else if ( input_sequence_type_ == "IUPAC" ) {
		seq_ = core::io::NomenclatureManager::annotated_sequence_from_IUPAC_sequence( seq_ );
	}

	Size const alignment_length = core::pose::rna::remove_bracketed( seq_ ).size();
	std::map <Size, Size> template_alignment2sequence, target_alignment2sequence;
	//if ( sequence_mask_.size() == 0 ) sequence_mask_ = ObjexxFCL::FArray1D_bool( alignment_length, true );


	std::cout << "Current length of pose is " << nres << std::endl;
	std::cout << "Target sequence: " << seq_ << std::endl;
	std::cout << "...has true length " << alignment_length << std::endl;


	target_sequence_ = seq_;
	// AMW: to handle all of this, we need a method that takes an annotated sequence and pushes back
	// residues worth of that. My thought is actually creating a second vector of
	// rtypes that only has the desired rtypes
	for ( Size i = 1; i <= alignment_length; i++ ) {
		//target_sequence += seq_[i-1];
		working_res_.push_back( target_alignment2sequence[i] /*+ offset_*/); // will be used for PDB numbering
	}
}


void
RNAThreadAndMinimizeMover::long_mutate_strategy( core::pose::Pose & pose ) {

	core::pose::set_reasonable_fold_tree( pose );

	utility::vector1< Size > changed_pos;
	utility::vector1< int> changed_pos_working;
	//Size cycles = 0;
	//core::kinematics::MoveMapOP mm;

	bool missing_residues = false;

	std::set< Size > remaining_res_in_pose;
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) remaining_res_in_pose.insert( ii );
	//while ( ! identical_save_for_deletions( pose.annotated_sequence(), target_sequence_ ) ) {
	while ( remaining_res_in_pose.size() != 0 ) {
		Size const i = numeric::random::rg().uniform() * remaining_res_in_pose.size();
		auto it = remaining_res_in_pose.begin();
		std::advance( it, i ); if ( it == remaining_res_in_pose.end() ) continue;

		if ( !rtypes_[*it] ) {
			// residue is marked for deletion
			remaining_res_in_pose.erase( it );
			missing_residues = true;
			continue;
		}

		if ( pose::rna::mutate_position( pose, *it, *rtypes_[*it] ) ) {
			changed_pos.push_back( *it );
			changed_pos_working.push_back( working_res_[*it] );
			utility::vector1< Size > just_one;
			just_one.push_back( *it );

			protocols::simple_moves::MinMoverOP minm(
				new protocols::simple_moves::MinMover(
				mm_from_residues( pose, just_one, true ),
				scorefxn_,
				"lbfgs_armijo_nonmonotone",
				0.001,
				true ) );
			// mm = mm_from_residues( pose, just_one, true );
			minm->apply( pose );

		}

		remaining_res_in_pose.erase( it );
	}

	std::cout << "Changed residues (without offset applied): " << make_tag_with_dashes( changed_pos ) << std::endl;
	std::cout << "Changed residues (in residue numbering): "    << make_tag_with_dashes( changed_pos_working ) << std::endl;

	if ( missing_residues ) {
		process_deletions( pose );
	}

	if ( insertion_list_.size() > 0 ) {
		process_insertions( pose );
	}
}

void RNAThreadAndMinimizeMover::process_deletions( core::pose::Pose & pose ) {
	std::set< Size > seqpos_for_deletion;
	for ( Size ii = 1; ii <= rtypes_.size(); ++ii ) {
		if ( ! rtypes_[ii] ) seqpos_for_deletion.insert( ii );
	}

	for ( auto const seqpos : seqpos_for_deletion ) {

		//TR << "pre:  " << pose.fold_tree() << std::endl;
		pose.delete_residue_slow( seqpos );

		//TR << "post: " << pose.fold_tree() << std::endl;

		// process deletions here, merge long deletions, etc.
		// oh, delete pose residue first.
		// use some kind of loopmodel around residue seq method...

		// Simple minimization to close
		accomodate_length_change( pose, seqpos, seqpos );
	}
}

void RNAThreadAndMinimizeMover::accomodate_length_change( Pose & pose, Size const insertion_begin, Size const seqpos ) {
	core::chemical::ResidueTypeSetCOP rsd_set = pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t );
	core::pose::Pose ref_pose = pose;

	Size const width = 3;
	Size begin = insertion_begin > width ? insertion_begin - width : 1;
	Size end   = seqpos < pose.total_residue() - width ? seqpos + width : pose.total_residue();
	core::pose::addVirtualResAsRoot( pose );
	auto old_ft = pose.fold_tree();


	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {

		if ( ! pose.residue_type( ii ).has( "P" ) ) continue;
		if ( ii >= insertion_begin && ii <= seqpos ) continue;

		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		Real const stdev = ( ii > begin || ii < end ) ? 6 : 2;
		pose.add_constraint( ConstraintOP( new CoordinateConstraint(
			id::AtomID( pose.residue_type( ii ).atom_index( "P" ), ii ),
			id::AtomID( 1, pose.total_residue() ),
			ref_pose.residue( ii ).xyz( pose.residue_type( ii ).atom_index( "P" ) ),
			FuncOP( new HarmonicFunc( 0.0, stdev ) ) ) ) );
	}

	///*
	// Real prior_rep = scorefxn_->get_weight( core::scoring::fa_rep );
	scorefxn_->set_weight( core::scoring::coordinate_constraint, 5 );
	// scorefxn_->set_weight( core::scoring::fa_rep, 0 );



	// AMW: this reflects a bug in this scoring term, but let's deal with that later
	core::scoring::ScoreFunctionOP loop_close_score = scorefxn_->clone();
	loop_close_score->set_weight( core::scoring::suiteness_bonus, 0 );
	loop_close_score->set_weight( core::scoring::rna_torsion, 0 );

	Pose recovered_low;
	//Real const temp = 1.0;
	Size const nrep = 3;
	for ( Size ii = 1; ii <= nrep; ++ii ) {
		Size moving = 0;
		while ( moving == 0 || moving == insertion_begin ) {
			moving = Size( numeric::random::rg().uniform() * ( end - begin ) + begin );
		}

		if ( pose::rna::mutate_position( pose, moving+1,
				rsd_set->get_residue_type_with_variant_added(
				pose.residue_type( moving+1 ), core::chemical::CUTPOINT_LOWER ) ) ) {
			std::cout << "Successfully applied cutpoint variant" << std::endl;
		}
		if ( pose::rna::mutate_position( pose, moving+2,
				rsd_set->get_residue_type_with_variant_added(
				pose.residue_type( moving+2 ), core::chemical::CUTPOINT_UPPER ) ) ) {
			std::cout << "Successfully applied cutpoint variant" << std::endl;
		}

		// Randomly perturb all bb torsions from begin to end
		Size const npert = 1000;
		Real pert_size = 40;
		for ( Size rr = 1; rr <= npert; ++rr ) {
			if ( rr > npert/2 ) pert_size = 20;
			if ( rr > 3*npert/4 ) pert_size = 10;
			if ( rr > 7*npert/8 ) pert_size = 5;

			Real old_score = ( *loop_close_score )( pose );

			std::map< Size, std::map< Size, Real > > old_torsions;
			for ( Size jj = begin; jj <= end; ++jj ) {
				old_torsions[ jj ] = std::map< Size, Real >();
				for ( Size kk = 1; kk <= pose.residue_type( jj ).mainchain_atoms().size(); ++kk ) {
					auto const tid = id::TorsionID( jj, id::BB, kk );
					old_torsions[ jj ][ kk ] = pose.torsion( tid );
					if ( jj == moving && kk == 6 ) continue;
					if ( jj == moving+1 && kk == 1 ) continue;

					pose.set_torsion( tid, pose.torsion( tid ) + pert_size * numeric::random::rg().gaussian() );
				}
			}

			Real new_score = ( *loop_close_score )( pose );

			if ( new_score < old_score ) {
				TR << "Perturbation " << rr << ": moving " << moving << " move accepted absolutely " << old_score << " to " << new_score << "." << std::endl;
				recovered_low = pose;
				// Temp is 10 - rr/10
			} else if ( numeric::random::rg().uniform() < std::exp( -1 * ( new_score - old_score ) / ((npert-rr)/10) ) ) {
				TR << "Perturbation " << rr << ": moving " << moving << " move accepted thermally  " << old_score << " to " << new_score << "." << std::endl;
			} else {
				TR << "Perturbation " << rr << ": moving " << moving << " move rejected! " << old_score << " to " << new_score << "." << std::endl;

				// Reset
				for ( auto const & res : old_torsions ) {
					for ( auto const & torsion : res.second ) {
						auto const tid = id::TorsionID( res.first, id::BB, torsion.first );
						pose.set_torsion( tid, torsion.second );
					}
				}
			}
		}


		core::kinematics::MoveMapOP ins_mm( new core::kinematics::MoveMap );
		ins_mm->set_chi( false );
		ins_mm->set_bb( false );
		// add additional changed pos - also add base pairing, maybe use slice mechanism
		for ( Size ii = begin; ii <= end; ++ii ) {
			ins_mm->set_bb( ii, true );
			ins_mm->set_chi( ii, true );
			for ( Size jj = 1; jj <= pose.residue_type( ii ).natoms(); ++jj ) {
				ins_mm->set( id::DOF_ID( id::AtomID( jj, ii ), id::D ), true );
				ins_mm->set( id::DOF_ID( id::AtomID( jj, ii ), id::THETA ), true );
			}
		}

		core::scoring::ScoreFunctionOP cart_scorefxn = loop_close_score->clone();
		cart_scorefxn->set_weight( core::scoring::cart_bonded, 1 );
		cart_scorefxn->set_weight( core::scoring::pro_close, 0 );
		for ( Real ii = 0.01; ii <= 500; ii *= 5 ) {
			//if ( ii > 0.01 ) cart_scorefxn->set_weight( core::scoring::fa_rep, prior_rep );
			cart_scorefxn->set_weight( core::scoring::chainbreak, ii );
			cart_scorefxn->set_weight( core::scoring::linear_chainbreak, ii );
			if ( ii > 0.2 ) scorefxn_->set_weight( core::scoring::coordinate_constraint, 1/ii );
			protocols::simple_moves::MinMoverOP ins_minm( new protocols::simple_moves::MinMover( ins_mm, cart_scorefxn, "lbfgs_armijo_nonmonotone", 1, true ) );
			ins_minm->apply( pose );
			ins_minm->apply( recovered_low );

		}



		if ( pose::rna::mutate_position( pose, moving+1,
				rsd_set->get_residue_type_with_variant_removed(
				pose.residue_type( moving+1 ), core::chemical::CUTPOINT_LOWER ) ) ) {
			std::cout << "Successfully removed cutpoint variant" << std::endl;
		}
		if ( pose::rna::mutate_position( pose, moving+2,
				rsd_set->get_residue_type_with_variant_removed(
				pose.residue_type( moving+2 ), core::chemical::CUTPOINT_UPPER ) ) ) {
			std::cout << "Successfully removed cutpoint variant" << std::endl;
		}

		if ( pose::rna::mutate_position( recovered_low, moving+1,
				rsd_set->get_residue_type_with_variant_removed(
				recovered_low.residue_type( moving+1 ), core::chemical::CUTPOINT_LOWER ) ) ) {
			std::cout << "Successfully removed cutpoint variant" << std::endl;
		}
		if ( pose::rna::mutate_position( pose, moving+2,
				rsd_set->get_residue_type_with_variant_removed(
				recovered_low.residue_type( moving+2 ), core::chemical::CUTPOINT_UPPER ) ) ) {
			std::cout << "Successfully removed cutpoint variant" << std::endl;
		}
	}

	// Restore old FT -- suppose loop is basically closed -- and:
	// 1. cartmin once
	// 2. optimize internal coords.

	pose.fold_tree( old_ft );

	core::kinematics::MoveMapOP ins_mm( new core::kinematics::MoveMap );
	ins_mm->set_chi( false );
	ins_mm->set_bb( false );
	// add additional changed pos - also add base pairing, maybe use slice mechanism
	for ( Size ii = begin; ii <= end; ++ii ) {
		ins_mm->set_bb( ii, true );
		ins_mm->set_chi( ii, true );
		for ( Size jj = 1; jj <= pose.residue_type( ii ).natoms(); ++jj ) {
			ins_mm->set( id::DOF_ID( id::AtomID( jj, ii ), id::D ), true );
			ins_mm->set( id::DOF_ID( id::AtomID( jj, ii ), id::THETA ), true );
		}
	}
	core::scoring::ScoreFunctionOP cart_scorefxn = scorefxn_->clone();
	cart_scorefxn->set_weight( core::scoring::cart_bonded, 1 );
	cart_scorefxn->set_weight( core::scoring::pro_close, 0 );

	protocols::simple_moves::MinMoverOP ins_minm( new protocols::simple_moves::MinMover( ins_mm, cart_scorefxn, "lbfgs_armijo_nonmonotone", 1, true ) );
	//ins_minm->cartesian( true );
	ins_minm->apply( pose );
	ins_minm->apply( recovered_low );


	// By the end, there shouldn't be much strain left in the structure--we can
	// remove coordinate constraints.
	pose.constraint_set( core::scoring::constraints::ConstraintSetOP( new core::scoring::constraints::ConstraintSet ) );

	ins_minm->apply( pose );
	ins_minm->apply( recovered_low );

	begin = insertion_begin > 3*width ? insertion_begin - 3*width : 1;
	end   = seqpos < pose.total_residue() - 3*width ? seqpos + 3*width : pose.total_residue();
	core::kinematics::MoveMapOP final_mm( new core::kinematics::MoveMap );
	final_mm->set_chi( false );
	final_mm->set_bb( false );
	for ( Size ii = begin; ii <= end; ++ii ) {
		final_mm->set_bb( ii, true );
		final_mm->set_chi( ii, true );
	}

	protocols::simple_moves::MinMoverOP final_minm( new protocols::simple_moves::MinMover( final_mm, scorefxn_, "lbfgs_armijo_nonmonotone", 0.001, true ) );
	final_minm->apply( pose );
	final_minm->apply( recovered_low );

	Real low = ( *scorefxn_ )( recovered_low );
	Real score = ( *scorefxn_ )( pose );
	if ( low < score ) {
		pose = recovered_low;
	}
}

void RNAThreadAndMinimizeMover::process_insertions( core::pose::Pose & pose ) {
	core::chemical::ResidueTypeSetCOP rsd_set = pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t );

	auto const insertions = utility::string_split( insertion_list_, ' ' );
	for ( auto const & insertion : insertions ) {

		std::istringstream iss( insertion );
		std::string to_be_added, seqpos_str;
		std::getline( iss, seqpos_str, ',');
		iss >> to_be_added;
		Size seqpos = stoi( seqpos_str );
		Size const insertion_begin = seqpos;
		if ( input_sequence_type_ == "MODOMICS" ) {
			to_be_added = core::io::NomenclatureManager::annotated_sequence_from_modomics_oneletter_sequence( to_be_added );
		} else if ( input_sequence_type_ == "IUPAC" ) {
			to_be_added = core::io::NomenclatureManager::annotated_sequence_from_IUPAC_sequence( to_be_added );
		}

		TR << "Appending " << to_be_added << " to seqpos " << seqpos << std::endl;

		// New rtypes to add
		auto const add_rtypes = core::pose::residue_types_from_sequence( to_be_added, *rsd_set, false);

		TR.Debug << "Now sequence " << pose.annotated_sequence() << std::endl;

		for ( auto const & add_rtype : add_rtypes ) {
			pose.conformation().append_polymer_residue_after_seqpos(
				*(new core::conformation::Residue( *add_rtype, true /*dummy*/ ) ), seqpos++, true );
		}

		TR.Debug << "Now sequence " << pose.annotated_sequence() << std::endl;

		accomodate_length_change( pose, insertion_begin, seqpos );
	}
}

void
RNAThreadAndMinimizeMover::mutate_all_at_once( core::pose::Pose & pose ) {
	utility::vector1< Size > changed_pos;
	utility::vector1< int> changed_pos_working;
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( pose::rna::mutate_position( pose, i, *rtypes_[i] ) ) {
			//if ( pose::rna::mutate_position( pose, i, special_res[i-1] ) ) {
			changed_pos.push_back( i );
			changed_pos_working.push_back( working_res_[i] );
		}
	}

	std::cout << "Changed residues (without offset applied): " << make_tag_with_dashes( changed_pos ) << std::endl;
	std::cout << "Changed residues (in residue numbering): "    << make_tag_with_dashes( changed_pos_working ) << std::endl;

	core::pose::set_reasonable_fold_tree( pose );

	// OK, we have a vector of changed residues. Setup a MoveMap accordingly.
	core::kinematics::MoveMapOP mm( mm_from_residues( pose, changed_pos, true ) );

	protocols::simple_moves::MinMoverOP minm( new protocols::simple_moves::MinMover( mm, scorefxn_, "lbfgs_armijo_nonmonotone", 0.001, true ) );

	minm->apply( pose );
}

core::kinematics::MoveMapOP
RNAThreadAndMinimizeMover::mm_from_residues( core::pose::Pose const & pose, utility::vector1< Size > const & changed_pos, bool add_nearby/*=false*/ ) {

	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_chi( false );
	mm->set_bb( false );
	if ( add_nearby ) {
		// add additional changed pos - also add base pairing, maybe use slice mechanism
		for ( Size ii = 1; ii <= changed_pos.size(); ++ii ) {
			if ( changed_pos[ ii ] > 1 ) {
				mm->set_bb( changed_pos[ ii ]-1, true );
				mm->set_chi( changed_pos[ ii ]-1, true );
			}
			if ( changed_pos[ ii ] < pose.total_residue() ) {
				mm->set_bb( changed_pos[ ii ]+1, true );
				mm->set_chi( changed_pos[ ii ]+1, true );
			}
		}
	}
	for ( Size ii = 1; ii <= changed_pos.size(); ++ii ) {
		mm->set_bb( changed_pos[ ii ], true );
		mm->set_chi( changed_pos[ ii ], true );
	}
	return mm;
}

void
RNAThreadAndMinimizeMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

}

protocols::moves::MoverOP
RNAThreadAndMinimizeMover::clone() const
{
	return protocols::moves::MoverOP( new RNAThreadAndMinimizeMover( *this ) );
}

protocols::moves::MoverOP
RNAThreadAndMinimizeMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new RNAThreadAndMinimizeMover );
}

std::string
RNAThreadAndMinimizeMover::get_name() const
{
	return RNAThreadAndMinimizeMover::class_name();
}

std::string
RNAThreadAndMinimizeMover::class_name()
{
	return "RNAThreadAndMinimizeMover";
}

void
RNAThreadAndMinimizeMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, RNAThreadAndMinimizeMover const & mover )
{
	mover.show(os);
	return os;
}

void
RNAThreadAndMinimizeMover::apply( core::pose::Pose & pose )
{
	core::pose::Pose const starting_pose( pose );

	// The target sequence needs a considerable amount of setup
	// to reconcile input formats and potential sequence alignment
	set_up_target_sequence( pose );

	// Obtain a vector of target residue types.
	obtain_rtypes_for_target_sequence();

	std::string current_sequence = pose.annotated_sequence();
	if ( target_sequence_ != current_sequence ) {
		std::cout << "TARGET:  " << target_sequence_ << std::endl;
	} else {
		std::cout << "MUTLIST: " << mutation_list_;
	}
	std::cout << "CURRENT: " << current_sequence << std::endl;

	std::string const clean_target_sequence = core::pose::rna::remove_bracketed( target_sequence_ );

	// Mutate sequence
	if ( long_strategy_ || insertion_list_ != "" || std::count( seq_.begin(), seq_.end(), '-' ) != 0 ) {
		long_mutate_strategy( pose );
	} else {
		mutate_all_at_once( pose );
	}
	core::pose::setPoseExtraScore( pose, "rms_from_starting", core::scoring::all_atom_rmsd( starting_pose, pose ) );
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
RNAThreadAndMinimizeMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new RNAThreadAndMinimizeMover );
}

std::string
RNAThreadAndMinimizeMoverCreator::keyname() const
{
	return RNAThreadAndMinimizeMover::class_name();
}

} //movers
} //rna
} //protocols

