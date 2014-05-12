// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/simple_filters/ScoreEvaluator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.hh>

// Project Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <basic/options/option.hh>

#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
// AUTO-REMOVED #include <iterator>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>
static basic::Tracer tr("protocols.evaluation.Score");


// option key includes
// option key includes

#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <protocols/simple_filters/RDC_Evaluator.hh>
#include <protocols/evaluation/util.hh>
#include <utility/vector1.hh>



// ObjexxFCL Headers

// Utility headers

//// C++ headers

namespace protocols {
namespace simple_filters {

using namespace core;

ScoreEvaluator::ScoreEvaluator( std::string tag, core::scoring::ScoreFunctionOP scorefxn, bool fullname ) :
	evaluation::SingleValuePoseEvaluator< Real >(  ( !fullname ? ("score"+tag) : tag ) ),
  scorefxn_ ( NULL )

{
	if( !scorefxn ) {
		utility_exit_with_message("Score evaluator score function is NULL");
	} else {
		scorefxn_ = scorefxn->clone();
	}
	tr.Info << "ScoreEvaluator: " << "score" << tag << " " << std::endl;
	scorefxn->show( tr.Info );
	tr.Info << std::endl;
}

ScoreEvaluator::~ScoreEvaluator() {}

core::Real
ScoreEvaluator::apply(
  core::pose::Pose& pose
) const {
	core::scoring::ScoreFunctionOP scorefxn( scorefxn_->clone() );
	core::pose::Pose chainbreak_pose( pose );
	core::pose::Pose nochainbreak_pose( pose );
	toolbox::pose_manipulation::add_chainbreaks_according_to_jumps( chainbreak_pose );
	//	js.remove_chainbreaks( nochainbreak_pose );  --- whatever comes in should be consistent with constraints..this might make it worse
	scoring::ScoreFunction chainbreaks_scfxn;
	chainbreaks_scfxn.set_weight(  scoring::linear_chainbreak, scorefxn->get_weight(  scoring::linear_chainbreak ) );
	chainbreaks_scfxn.set_weight(  scoring::overlap_chainbreak, scorefxn->get_weight(  scoring::overlap_chainbreak ) );
	chainbreaks_scfxn.set_weight(  scoring::chainbreak, scorefxn->get_weight(  scoring::chainbreak ) );
	scorefxn->set_weight(  scoring::linear_chainbreak, 0);
	scorefxn->set_weight(  scoring::overlap_chainbreak, 0);
	scorefxn->set_weight(  scoring::chainbreak, 0);
	core::Real val = scorefxn->score( nochainbreak_pose );
	return val + chainbreaks_scfxn( chainbreak_pose );
}

bool ScoreEvaluator::applicable( core::pose::Pose const&pose ) const {
	bool centroid( scorefxn_->get_weight( scoring::vdw ) > 0.0 );
	bool fa( scorefxn_->get_weight( scoring::fa_atr ) > 0.0 );
	bool pose_fa( pose.is_fullatom() );
	return (!fa && !centroid) || ( pose_fa && fa && !centroid ) || ( !pose_fa && centroid && !fa );
}

TruncatedScoreEvaluator::TruncatedScoreEvaluator(
	 std::string tag,
   core::scoring::ResidueSelectionVector const& selection,
 	 core::scoring::ScoreFunctionOP scorefxn,
	 bool fullname
) :
	ScoreEvaluator( tag, scorefxn, fullname ),
	selection_( selection ),
	rdcs_( NULL )
{
	nres_ = 400; //no worries: if the pose turns out to be larger we repeat that inversion in apply().
	evaluation::invert_include_residues( 400, selection_, exclude_list_ );
	if ( tr.Trace.visible() ) {
		for ( core::scoring::ResidueSelectionVector::const_iterator it = selection.begin(); it != selection.end(); ++it ) {
			tr.Trace << *it << " ";
		}
		tr.Trace << std::endl;
		scorefxn->show( tr.Trace );
		tr.Trace << std::endl;
	}
	if ( !rdcs_ && scorefxn->get_weight( scoring::rdc ) > 0.0 && basic::options::option[ basic::options::OptionKeys::in::file::rdc ].user() ) {
		rdcs_ = new SelectRDC_Evaluator( selection_, "none");
	}
}

TruncatedScoreEvaluator::~TruncatedScoreEvaluator() {}

Real TruncatedScoreEvaluator::apply( core::pose::Pose& pose ) const {
	//	core::Pose my_pose( pose );
	PROF_START( basic::TRUNCATED_SCORE_EVALUATOR );
	scoring::ScoreFunctionOP scorefxn( scorefxn_ );
	if ( !scorefxn ) {
		tr.Trace << "no scorefunction specified in TruncatedScoreEvaluator... make appropriate standard score " << std::endl;
		if ( pose.is_fullatom() ) {
			scorefxn = core::scoring::getScoreFunction();
		} else {
			scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
		}
		scorefxn->set_weight( scoring::linear_chainbreak, 1.0 );
		scorefxn->set_weight( scoring::overlap_chainbreak, 1.0 );
		if ( pose.constraint_set()->has_residue_pair_constraints()  ) {
			scorefxn->set_weight( scoring::atom_pair_constraint, basic::options::option[ basic::options::OptionKeys::constraints::cst_weight ]() );
			scorefxn->set_weight( scoring::angle_constraint, basic::options::option[ basic::options::OptionKeys::constraints::cst_weight ]() );
			scorefxn->set_weight( scoring::dihedral_constraint, basic::options::option[ basic::options::OptionKeys::constraints::cst_weight ]() );
		}
	}

	if ( pose.total_residue() != nres_ ) {
		nres_ = pose.total_residue();
		evaluation::invert_include_residues( nres_, selection_, exclude_list_ );
	}

	core::pose::Pose chainbreak_pose( pose );
	core::pose::Pose nochainbreak_pose( pose );
	toolbox::pose_manipulation::add_chainbreaks_according_to_jumps( chainbreak_pose );
	//	js.remove_chainbreaks( nochainbreak_pose );
	scoring::ScoreFunction chainbreaks_scfxn;
	chainbreaks_scfxn.set_weight(  scoring::linear_chainbreak, scorefxn->get_weight(  scoring::linear_chainbreak ) );
	chainbreaks_scfxn.set_weight(  scoring::overlap_chainbreak, scorefxn->get_weight(  scoring::overlap_chainbreak ) );
	chainbreaks_scfxn.set_weight(  scoring::chainbreak, scorefxn->get_weight(  scoring::chainbreak ) );
	scorefxn->set_weight(  scoring::linear_chainbreak, 0);
	scorefxn->set_weight(  scoring::overlap_chainbreak, 0);
	scorefxn->set_weight(  scoring::chainbreak, 0);

	tr.Debug << "compute score without using these residues: ";
	std::copy( exclude_list_.begin(), exclude_list_.end(), std::ostream_iterator< Size >( tr.Debug, " " ));
	tr.Debug << std::endl;
	(*scorefxn)( nochainbreak_pose );
	core::Real val = scorefxn->get_sub_score_exclude_res( nochainbreak_pose, exclude_list_ );
	PROF_STOP( basic::TRUNCATED_SCORE_EVALUATOR );
	if ( rdcs_ ) val += scorefxn->get_weight( scoring::rdc )*rdcs_->apply( nochainbreak_pose );
	val += chainbreaks_scfxn( chainbreak_pose );
	return val;
}

}
}
