// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/eval/EnergyEval.cc
/// @brief  scores a fragment by inserting its backbone angles into a Pose
///         and evaluating its energy using a given ScoreFunction
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/fragment/picking_old/vall/eval/EnergyEval.hh>

// project headers
#include <basic/Tracer.hh>

// numeric headers
#include <numeric/random/random.hh>

//Auto Headers
#include <core/fragment/picking_old/concepts/Extent.hh>
#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace eval {


// static initialization

static thread_local basic::Tracer TR( "core.fragment.picking_old.vall.eval.EnergyEval" );


/// @brief default constructor
EnergyEval::EnergyEval() :
	Super(),
	insert_position_( 0 ),
	score_function_( ScoreFunctionOP( new core::scoring::ScoreFunction ) ),
	randomize_( false )
{}


/// @brief constructor
/// @param[in] pose insert backbone angles using a copy of this Pose
/// @param[in] insert_position insert backbone angles starting from this
///  position in the Pose
/// @param[in] score_function evaluate the Pose using a copy of this
///  ScoreFunction
/// @param[in] randomize flags that indicates whether a small amount
///  of noise between [0, 0.000001) will be added to the energy
EnergyEval::EnergyEval(
	Pose const & pose,
	Size const insert_position,
	ScoreFunction const & score_function,
	bool const randomize
) :
	Super(),
	pose_( pose ),
	insert_position_( insert_position ),
	score_function_( score_function.clone() ),
	randomize_( randomize )
{}


/// @brief default copy constructor
EnergyEval::EnergyEval( EnergyEval const & rval ) :
	Super( rval ),
	pose_( rval.pose_ ),
	insert_position_( rval.insert_position_ ),
	score_function_( rval.score_function_->clone() ),
	randomize_( rval.randomize_ )
{}


/// @brief default destructor
EnergyEval::~EnergyEval() {}


/// @brief copy assignment
EnergyEval & EnergyEval::operator =( EnergyEval const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );

		pose_ = rval.pose_;
		insert_position_ = rval.insert_position_;
		score_function_ = rval.score_function_->clone();
		randomize_ = rval.randomize_;
	}
	return *this;
}


/// @brief clone this object
VallFragmentEvalOP EnergyEval::clone() const {
	return VallFragmentEvalOP( new EnergyEval( *this ) );
}


/// @brief for a fragment extent, evaluate and store results in a VallFragmentScore
/// @return true, so score is always stored during VallLibrarian::catalog()
bool EnergyEval::eval_impl(
	Extent const & extent,
	VallFragmentScore & fs
)
{
	// insert backbone angles
	Size position = insert_position_;
	for ( VallResidueConstIterator i = extent.begin; i != extent.end; ++i, ++position ) {
		pose_.set_phi( position, i->phi() );
		pose_.set_psi( position, i->psi() );
		pose_.set_omega( position, i->omega() );
	}

	// evaluate the energy
	fs.score += score_function_->score( pose_ );

	if ( randomize_ ) {
		fs.score += ( numeric::random::rg().uniform() * 0.000001 );
	}

	return true;
}


/// @brief operation to be perform before catalog() starts
void EnergyEval::pre_catalog_op( VallLibrary const & ) {
	score_function_->show_line_headers( TR );
	TR << std::endl;
}


} // namespace eval
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core

