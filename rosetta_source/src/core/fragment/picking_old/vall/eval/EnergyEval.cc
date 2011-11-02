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

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/fragment/BBTorsionSRFD.fwd.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/Frame.fwd.hh>
#include <core/fragment/SecstructSRFD.hh>
#include <core/fragment/SingleResidueFragData.fwd.hh>
#include <core/fragment/SingleResidueFragData.hh>
#include <core/fragment/picking_old/concepts/Book.fwd.hh>
#include <core/fragment/picking_old/concepts/Book.hh>
#include <core/fragment/picking_old/concepts/Extent.fwd.hh>
#include <core/fragment/picking_old/concepts/Extent.hh>
#include <core/fragment/picking_old/vall/VallLibrary.fwd.hh>
#include <core/fragment/picking_old/vall/VallResidue.fwd.hh>
#include <core/fragment/picking_old/vall/VallResidue.hh>
#include <core/fragment/picking_old/vall/VallSection.fwd.hh>
#include <core/fragment/picking_old/vall/VallSection.hh>
#include <core/fragment/picking_old/vall/eval/EnergyEval.fwd.hh>
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.fwd.hh>
#include <core/fragment/picking_old/vall/eval/VallFragmentEval.hh>
#include <core/fragment/picking_old/vall/scores/VallFragmentScore.fwd.hh>
#include <core/fragment/picking_old/vall/scores/VallFragmentScore.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/uniform.fwd.hh>
#include <numeric/random/uniform.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>



namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace eval {


// static initialization
static numeric::random::RandomGenerator RG( 167872 ); // magic number, don't change

static basic::Tracer TR( "core.fragment.picking_old.vall.eval.EnergyEval" );


/// @brief default constructor
EnergyEval::EnergyEval() :
	Super(),
	insert_position_( 0 ),
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
	score_function_( score_function ),
	randomize_( randomize )
{}


/// @brief default copy constructor
EnergyEval::EnergyEval( EnergyEval const & rval ) :
	Super( rval ),
	pose_( rval.pose_ ),
	insert_position_( rval.insert_position_ ),
	score_function_( rval.score_function_ ),
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
		score_function_ = rval.score_function_;
		randomize_ = rval.randomize_;
	}
	return *this;
}


/// @brief clone this object
VallFragmentEvalOP EnergyEval::clone() const {
	return new EnergyEval( *this );
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
	fs.score += score_function_( pose_ );

	if ( randomize_ ) {
		fs.score += ( RG.uniform() * 0.000001 );
	}

	return true;
}


/// @brief operation to be perform before catalog() starts
void EnergyEval::pre_catalog_op( VallLibrary const & ) {
	score_function_.show_line_headers( TR );
	TR << std::endl;
}


} // namespace eval
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core

