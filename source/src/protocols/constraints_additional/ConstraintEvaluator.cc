// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PoseEvaluator
/// @brief PoseEvaluator
/// @details
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/constraints_additional/ConstraintEvaluator.hh>

// Package Headers
#include <protocols/constraints_additional/MaxSeqSepConstraintSet.hh>

#include <protocols/constraints_additional/AdditionalConstraintCreators.hh>

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/id/Exceptions.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>


#include <core/kinematics/FoldTree.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>
#include <utility>
#include <utility/file/file_sys_util.hh>

#include <core/id/SequenceMapping.hh>
#include <utility/vector1.hh>


// C++ headers

static basic::Tracer tr( "protocols.constraints_additional.ConstraintEvaluator" );

namespace protocols {
namespace constraints_additional {

using namespace core;
using namespace scoring;
using namespace constraints;

ConstraintEvaluator::ConstraintEvaluator( std::string tag, ConstraintSet const& cst, Size /*viol_level*/, Real /*threshold*/,  Size max_seq_sep )
: name_(std::move( tag )),
	constraints_( core::scoring::constraints::ConstraintSetOP( new ConstraintSet( cst ) ) ),
	tried_fullatom_pose_( false ),
	tried_centroid_pose_( false ),
	file_name_( "" ),
	// viol_level_ ( viol_level  ),
	// threshold_( threshold ),
	max_seq_sep_( max_seq_sep ),
	constraints_combine_ratio_( 1 ),
	cst_source_( "n/a" )
{}

ConstraintEvaluator::ConstraintEvaluator( std::string tag, ConstraintCOPs const& csts, Size /*viol_level*/, Real /*threshold*/, Size max_seq_sep )
: name_(std::move( tag )),
	constraints_( core::scoring::constraints::ConstraintSetOP( new ConstraintSet() ) ),
	tried_fullatom_pose_( false ),
	tried_centroid_pose_( false ),
	file_name_( "" ),
	// viol_level_ ( viol_level  ),
	// threshold_( threshold ),
	max_seq_sep_( max_seq_sep ),
	constraints_combine_ratio_( 1 ),
	cst_source_( "n/a" )
{
	constraints_->add_constraints( csts );
}

ConstraintEvaluator::ConstraintEvaluator( std::string tag, std::string file_name, Size /*viol_level*/, Real /*threshold*/, Size max_seq_sep )
: name_(std::move( tag )),
	constraints_( /* NULL */ ),
	tried_fullatom_pose_( false ),
	tried_centroid_pose_( false ),
	file_name_(std::move( file_name )),
	// viol_level_ ( viol_level  ),
	// threshold_( threshold ),
	max_seq_sep_( max_seq_sep ),
	constraints_combine_ratio_( 1 ),
	cst_source_( "n/a" )
{
	//check file exists
	if ( !utility::file::file_exists( file_name_ ) ) {
		utility_exit_with_message(" could not find file " + file_name_ );
	}
}


void ConstraintEvaluator::prepare_pose( core::pose::Pose const& pose_in, core::pose::Pose& pose ) const {

	toolbox::pose_manipulation::remove_chainbreaks_according_to_jumps( pose );

	ConstraintSetOP now_cst = pose.is_fullatom() ? fa_constraints_ : constraints_;
	if ( !now_cst ) {

		if ( pose.is_fullatom() && tried_fullatom_pose_ ) {
			pose.constraint_set( nullptr );
			return;
		}
		if ( !pose.is_fullatom() && tried_centroid_pose_ ) {
			pose.constraint_set( nullptr );
			return;
		}

		runtime_assert( utility::file::file_exists( file_name_) ); //it has already been checked... here it is an assertion

		using namespace core::scoring::constraints;
		//  ConstraintCreatorCOP new_atom_pair_creator( new constraints_additional::NamedAtomPairConstraintCreator );
		//  ConstraintCreatorCOP orig_atom_pair_creator( ConstraintFactory::get_instance()->get_creator( "AtomPair" ) ); // <-- this may actually be a NamedAtomPairConstraintCreator, we don't know; restore it, when done.
		//  ConstraintFactory::get_instance()->replace_creator( new_atom_pair_creator );
		try{
			now_cst = ConstraintIO::get_instance()->read_constraints( file_name_, ConstraintSetOP( new ConstraintSet ), pose );
			scoring::constraints::ConstraintCOPs added_constraints = now_cst->get_all_constraints();
			kinematics::ShortestPathInFoldTree sp( pose.fold_tree() );
			scoring::constraints::choose_effective_sequence_separation( sp, added_constraints );
			utility::vector1< bool > combine_exclude_res;
			scoring::constraints::combine_constraints( added_constraints, constraints_combine_ratio_, combine_exclude_res, sp ); // if combine_ratio_ > 1 this wil
			now_cst = ConstraintSetOP( new ConstraintSet() );
			now_cst->add_constraints( added_constraints );
		} catch ( core::id::EXCN_AtomNotFound& excn ) {
			tr.Warning << " cannot use constraint file " << file_name_ << " on " << ( pose.is_fullatom() ? " fullatom " : " centroid " ) << " pose " << std::endl;
			tr.Warning << " because: " << excn << std::endl;
			now_cst = nullptr;
			pose.constraint_set( nullptr );
			if ( pose.is_fullatom() ) tried_fullatom_pose_ = true;
			else tried_centroid_pose_ = true;
			return;
		}

		//restore original cst-type
		//ConstraintIO::get_cst_factory().add_type( orig_atom_pair_type );
		//ConstraintFactory::get_instance()->replace_creator( orig_atom_pair_creator );

	}

	runtime_assert( now_cst != nullptr );

	constraints_additional::MaxSeqSepConstraintSetOP new_cst( nullptr );
	if ( max_seq_sep_ > 0 ) {
		new_cst = constraints_additional::MaxSeqSepConstraintSetOP( new constraints_additional::MaxSeqSepConstraintSet( *now_cst, pose.fold_tree() ) );
		new_cst->set_max_seq_sep( max_seq_sep_ );
	} else {
		//if not specified  we copy the max_seq_separation if present
		constraints_additional::MaxSeqSepConstraintSetCOP ms_set = utility::pointer::dynamic_pointer_cast< constraints_additional::MaxSeqSepConstraintSet const > ( pose_in.constraint_set() );
		if ( ms_set ) {
			new_cst = constraints_additional::MaxSeqSepConstraintSetOP( new constraints_additional::MaxSeqSepConstraintSet( *now_cst, pose.fold_tree()  ) );
			new_cst->set_max_seq_sep( ms_set->max_seq_sep() );
		}
	}

	if ( new_cst ) {
		pose.constraint_set( new_cst );
	} else {
		pose.constraint_set( now_cst );
	}

	if ( pose.is_fullatom() ) fa_constraints_ = now_cst;
	else constraints_ = now_cst;

}

Real ConstraintEvaluator::apply( core::pose::Pose& pose_in ) const {
	// PROF_START( basic::TEST1 );
	pose::Pose pose( pose_in );
	// PROF_STOP( basic::TEST1 );

	// PROF_START( basic::TEST2 );
	prepare_pose( pose_in, pose );
	// PROF_STOP( basic::TEST2 );

	// PROF_START( basic::TEST3 );
	ScoreFunction scfxn;
	scfxn.set_weight( atom_pair_constraint, 1.0 );
	core::Real score( scfxn( pose ) );
	// PROF_STOP( basic::TEST3 );

	return score;
}

std::string ConstraintEvaluator::name( core::Size i ) const {
	if ( i == 1 ) { return name_; }
	// if ( i == 2 ) { return name_ + "_viol"; }
	runtime_assert( i <= 1 && i > 0 );
	return ""; //make compiler happy
}


void ConstraintEvaluator::apply( core::pose::Pose& pose_in, std::string, core::io::silent::SilentStruct &pss ) const {
	// PROF_START( basic::TEST1 );
	pose::Pose pose( pose_in );
	// PROF_STOP( basic::TEST1 );

	// PROF_START( basic::TEST2 );
	prepare_pose( pose_in, pose );
	// PROF_STOP( basic::TEST2 );

	// PROF_START( basic::TEST3 );
	ScoreFunction scfxn;
	scfxn.set_weight( atom_pair_constraint, 1.0 );
	core::Real score( scfxn( pose ) );
	// PROF_STOP( basic::TEST3 );

	pss.add_energy( name( 1 ), score );
	if ( cst_source_ != "n/a" ) {
		pss.add_string_value( "cst_source_"+name(1), cst_source_ );
	}
	// pss.add_energy( name( 2 ), pose.constraint_set()->show_violations( tr.Info, pose, viol_level_, threshold_ ) );
	// return scfxn( pose );

}

}
}
