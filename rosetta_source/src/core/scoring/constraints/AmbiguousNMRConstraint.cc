// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief contains declarations for a type of constraint that holds a number of constraints
/// @brief where only the lowest one is evaluated
/// @author Florian Richter (floric@u.washington.edu, march 2008)


#include <core/scoring/constraints/AmbiguousNMRConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/FuncFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/prof.hh>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/AmbiguousNMRConstraint.fwd.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintCreator.fwd.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/MultiConstraint.fwd.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
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
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.fwd.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <time.h>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>


namespace core {
namespace scoring {
namespace constraints {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
AmbiguousNMRConstraint::AmbiguousNMRConstraint( FuncOP func ):
	MultiConstraint( atom_pair_constraint ),
	func_( func )
{
	//	init_cst_score_types();
	assert ( member_constraints().size() == 0 );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
AmbiguousNMRConstraint::AmbiguousNMRConstraint( ConstraintCOPs const& cst_in, FuncOP func ):
	MultiConstraint( cst_in, atom_pair_constraint ),
	func_( func )
{
	//	init_cst_score_types();
	assert ( member_constraints().size() > 0 );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
// void
// AmbiguousNMRConstraint::init_cst_score_types()
// {
// 	cst_score_types_.clear();
// 	cst_score_types_.push_back(atom_pair_constraint);
// }
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief ScoreFunction, scores all member constraints but only reports the lowest one
void
AmbiguousNMRConstraint::score( XYZ_Func const & xyz_func, EnergyMap const & /*weights*/, EnergyMap & emap ) const
{

	core::Real cum_invdist6 = 0;

	for( ConstraintCOPs::const_iterator member_it = member_constraints().begin(); member_it != member_constraints().end(); member_it++){
		AmbiguousNMRDistanceConstraint const* cst_in_casted;
		cst_in_casted = dynamic_cast< AmbiguousNMRDistanceConstraint const* >( (*member_it).get() );
		if ( cst_in_casted ) cum_invdist6 += cst_in_casted->inv_dist6( xyz_func );
		if ( !cst_in_casted ) {
			AtomPairConstraint const* cst_in_casted;
			cst_in_casted = dynamic_cast< AtomPairConstraint const* >( (*member_it).get() );
			if ( cst_in_casted ) {
				Real dist = cst_in_casted->dist( xyz_func );
				Real inv_dist = 1.0/dist;
				Real inv_dist2 = inv_dist*inv_dist;
				cum_invdist6 += inv_dist2*inv_dist2*inv_dist2;
			} else {
				runtime_assert( 0 == 1 );
			}
		}
	}
	//add lowest score to the actual emap
	Real eff_dist = pow( cum_invdist6, -1.0/6 );
	emap[ score_type() ] +=  get_func().func( eff_dist );
} //score


core::Real
AmbiguousNMRConstraint::dist( core::pose::Pose const& pose ) const
{
	return dist( ConformationXYZ( pose.conformation() ) );
}

core::Real
AmbiguousNMRConstraint::dist( XYZ_Func const & xyz ) const
{
	core::Real cum_invdist6 = 0;

	for( ConstraintCOPs::const_iterator member_it = member_constraints().begin(); member_it != member_constraints().end(); member_it++){
		AmbiguousNMRDistanceConstraint const* cst_in_casted;
		PROF_START( basic::NOESY_ASSIGN_DIST_CST_CAST );
		cst_in_casted = dynamic_cast< AmbiguousNMRDistanceConstraint const* >( (*member_it).get() );
		PROF_STOP( basic::NOESY_ASSIGN_DIST_CST_CAST );
		if ( cst_in_casted ) cum_invdist6 += cst_in_casted->inv_dist6( xyz );
		if ( !cst_in_casted ) {
			AtomPairConstraint const* cst_in_casted;
			PROF_START( basic::NOESY_ASSIGN_DIST_CST_CAST );
			cst_in_casted = dynamic_cast< AtomPairConstraint const* >( (*member_it).get() );
			PROF_STOP( basic::NOESY_ASSIGN_DIST_CST_CAST );
			if ( cst_in_casted ) {
				Real dist = cst_in_casted->dist( xyz );
				Real inv_dist = 1.0/dist;
				Real inv_dist2 = inv_dist*inv_dist;
				cum_invdist6 += inv_dist2*inv_dist2*inv_dist2;
			} else {
				runtime_assert( 0 == 1 );
			}
		}
	}
	//add lowest score to the actual emap
	Real eff_dist = pow( cum_invdist6, -1.0/6 );
	return eff_dist;
}

ConstraintOP
AmbiguousNMRConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	ConstraintCOPs new_csts;
	for( ConstraintCOPs::const_iterator cst_it = member_constraints_.begin(); cst_it != member_constraints_.end(); ++cst_it ){
		ConstraintOP new_cst = (*cst_it)->remap_resid( seqmap );
		if( new_cst ) new_csts.push_back( new_cst );
	}
	if( new_csts.size() > 0 ){
		return ConstraintOP( new AmbiguousNMRConstraint( new_csts, get_func().clone() ) );
	}
	else return NULL;
}


/// @brief function to minimize lowest scoring member constraint
void
AmbiguousNMRConstraint::fill_f1_f2(
	AtomID const & atom,
	XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{
	//active_constraint_->fill_f1_f2(atom, conformation, F1, F2, weights);
	Real eff_dist = dist( xyz );
	Real out_wderiv( weights[ score_type() ] * get_func().dfunc( eff_dist ));
	Real in_deriv = -1.0/6.0 * pow( eff_dist, 7.0 );
	//	tr.Trace << "deriv for atom " << atom << eff_dist << " " << out_wderiv << " " << in_deriv << std::endl;

	//	tr.Trace << "the_other_atoms: " << the_other_atoms.size() << " " << the_other_atoms.front() << std::endl;
	//Vector f1(0.0), f2(0.0);
	for( ConstraintCOPs::const_iterator member_it = member_constraints().begin(); member_it != member_constraints().end(); member_it++){
              Vector f1(0.0), f2(0.0);

               //fpd hack to get at vector from atom i->j
                (*member_it)->fill_f1_f2( atom, xyz, f1, f2, weights );
    core::Real member_cst_is_scaled_by = weights[ (*member_it)->score_type() ] * (*member_it)->get_func().dfunc( eff_dist
);

               core::Real scale_i = -6.0*pow((*member_it)->dist( xyz ),-7.0);

               if (std::fabs(member_cst_is_scaled_by) > 1e-14) scale_i /= member_cst_is_scaled_by;

               F1 += 1.0 * scale_i * out_wderiv * in_deriv * f1;
               F2 += 1.0 * scale_i * out_wderiv * in_deriv * f2;
        }
        //              tr.Trace << "wderiv " << wderiv << std::endl;	


	//		tr.Trace << "wderiv " << wderiv << std::endl;
//	F1 += out_wderiv * in_deriv * f1;
//	F2 += out_wderiv * in_deriv * f2;
}


// void
// AmbiguousNMRConstraint::show( std::ostream& out) const
// {
// 	out << "AmbiguousNMRConstraint Active constraint:" << std::endl;
// 	out << "AmbiguousNMRConstraint containing the following " << member_constraints().size() << " constraints: " << std::endl;
// 	for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
// 		(*cst_it)->show(out);
// 	}

// 	out << " ...all member constraints of this AmbiguousNMRConstraint shown." << std::endl;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////
///@details read definition of a multiconstraint. Since a MultiConstraint is essentially a vector of
void
AmbiguousNMRConstraint::read_def(
	std::istream& data,
	core::pose::Pose const& pose,
	FuncFactory const & func_factory
)
{
	std::string func_type;
	data >> func_type;

	func_ = func_factory.new_func( func_type );
	func_->read_data( data );

	if ( data.good() ) {
	//chu skip the rest of line since this is a single line defintion.
		while( data.good() && (data.get() != '\n') ) {}
		if ( !data.good() ) data.setstate( std::ios_base::eofbit );
	}
	MultiConstraint::read_def( data, pose, func_factory );
}

void AmbiguousNMRConstraint::show_def( std::ostream& out, pose::Pose const& pose ) const {
	out << type() << " ";
	if ( func_ ) func_->show_definition( out );
	else out << std::endl;
  for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
    (*cst_it)->show_def( out, pose );
		//		out<<std::endl;
  }
	out << "End_"<< type() << std::endl;
}

Size
AmbiguousNMRConstraint::show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold ) const
{
	Size total_viol = 0;
	bool passed( false );
	for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
		Size viol = (*cst_it)->show_violations( out, pose, verbose_level, threshold);
		if ( viol == 0 && verbose_level > 70 ) {
			passed = true;
			out << "\nResiduePairConstraints (" <<  (*cst_it)->atom(1).rsd() << ", " << (*cst_it)->atom((*cst_it)->natoms()).rsd() << " ) . of total: 1  0 violated" << std::endl;
		}
		total_viol += viol;
	}
	if ( !passed && verbose_level > 70 ) {
		for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
			out << "\nResiduePairConstraints (" <<  (*cst_it)->atom(1).rsd() << ", " << (*cst_it)->atom((*cst_it)->natoms()).rsd() << " ) + of total: 1  1 violated" << std::endl;
		}
	}
	return func_->show_violations( out, dist( pose ), verbose_level, threshold );
}

void
AmbiguousNMRConstraint::add_individual_constraint( ConstraintCOP cst_in )
{
	Constraint const* cst_in_casted;
	cst_in_casted = dynamic_cast< AmbiguousNMRDistanceConstraint const* >( cst_in.get() );
	if ( !cst_in_casted ) {
		cst_in_casted = dynamic_cast< AtomPairConstraint const* >( cst_in.get() );
	}
	if ( !cst_in_casted ) {
		throw utility::excn::EXCN_BadInput( "failed attempt to add " + cst_in->type() + " to AmbiguousNMRConstraint. Can only add AmbiguousNMRDistanceConstraint and AtomPairConstraint");
	}

	MultiConstraint::add_individual_constraint( cst_in );
}


} //constraints
} //scoring
} //core
