// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/scoring/methods/ResidualDipolarCouplingEnergy_Rohl.cc
/// @brief  RDC energy - comparing experimental RDC values to calculated values
/// @author Srivatsan Raman


//Unit headers
#include <core/scoring/methods/ResidualDipolarCouplingEnergy_Rohl.hh>
#include <core/scoring/methods/ResidualDipolarCouplingEnergy_RohlCreator.hh>
#include <core/scoring/ResidualDipolarCoupling_Rohl.hh>
#include <core/scoring/ResidualDipolarCoupling_Rohl.fwd.hh>

//Package headers

#include <core/conformation/Residue.hh>
//#include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
#include <core/pose/Pose.hh>
//#include <core/pose/datacache/CacheableDataType.hh>

//numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

//utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/exit.hh>

//Objexx headers
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
// AUTO-REMOVED #include <ObjexxFCL/char.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/Fmath.hh>

//C++ headers
#include <iostream>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
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
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
// AUTO-REMOVED #include <core/scoring/DerivVectorPair.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
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
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>

//Auto Headers



namespace core {
namespace scoring {
namespace methods {

using namespace ObjexxFCL;

/// @details This must return a fresh instance of the ResidualDipolarCouplingEnergy_Rohl class,
/// never an instance already in use
methods::EnergyMethodOP
ResidualDipolarCouplingEnergy_RohlCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new ResidualDipolarCouplingEnergy_Rohl;
}

ScoreTypes
ResidualDipolarCouplingEnergy_RohlCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rdc_rohl );
	return sts;
}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
ResidualDipolarCouplingEnergy_Rohl::ResidualDipolarCouplingEnergy_Rohl() :
	parent( new ResidualDipolarCouplingEnergy_RohlCreator )
{}


//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
EnergyMethodOP
ResidualDipolarCouplingEnergy_Rohl::clone() const
{
  return new ResidualDipolarCouplingEnergy_Rohl();
}

//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
void ResidualDipolarCouplingEnergy_Rohl::finalize_total_energy(
  pose::Pose & pose,
  ScoreFunction const &,
  EnergyMap & totals
) const
{

	Real dipolar_score = eval_dipolar( pose );
	totals[ rdc_rohl ] = dipolar_score;

}

//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
ResidualDipolarCoupling_Rohl const &
ResidualDipolarCouplingEnergy_Rohl::rdc_from_pose(
	pose::Pose & pose
) const
{
// 	//using core::pose::datacache::CacheableDataType::RESIDUAL_DIPOLAR_COUPLING_DATA;

// 	if( pose.data().has( RESIDUAL_DIPOLAR_COUPLING_DATA ) )
// 		return *( static_cast< ResidualDipolarCoupling const * >( pose.data().get_const_ptr( RESIDUAL_DIPOLAR_COUPLING_DATA )() ) );

// 	ResidualDipolarCouplingOP rdc_info = new ResidualDipolarCoupling;
// 	pose.data().set( RESIDUAL_DIPOLAR_COUPLING_DATA, rdc_info );
 	ResidualDipolarCoupling_RohlOP rdc_info( retrieve_RDC_ROHL_from_pose( pose ) );
	if ( !rdc_info ) {
		rdc_info = new ResidualDipolarCoupling_Rohl;
		store_RDC_ROHL_in_pose( rdc_info, pose );
	}
	return *rdc_info;
}

//////////////////////////////////////////////////////
//@brief main computation routine for RDC energy... everything is happening here right now.
// this has to be spread out over different routines to make this energy yield derivatives
//////////////////////////////////////////////////////
Real ResidualDipolarCouplingEnergy_Rohl::eval_dipolar(
  pose::Pose & pose
) const
{

	ResidualDipolarCoupling_Rohl const & rdc_data( rdc_from_pose( pose ) );
	utility::vector1< core::scoring::RDC_Rohl > All_RDC_lines( rdc_data.get_RDC_data() );

	Size const nrow( All_RDC_lines.size() ); //number of experimental couplins
	Size const ORDERSIZE = { 5 }; //Syy,Szz,Sxy,Sxz,Syz

	ObjexxFCL::FArray2D< Real > A( nrow, ORDERSIZE ); // N x 5, the 5 assymetric tensor elements per relevant vector in pose
	ObjexxFCL::FArray1D< Real > b( nrow );            // experimental values
	ObjexxFCL::FArray1D< Real > x( ORDERSIZE ); //previously dimensioned to nrow, which is wrong.
	ObjexxFCL::FArray1D< Real > weights( nrow ); //previously dimensioned to nrow, which is wrong.
	ObjexxFCL::FArray2D< Real > vec( 3, 3 );
	Real Azz, eta;
	bool reject = false;


	assemble_datamatrix( pose, All_RDC_lines, A, b, weights);
 	/*	for ( core::Size i = 1; i <= nrow; ++i ) {
		std::cout << "Matrices A & b  " << A(i,1) << " " << A(i,2) << " " << A(i,3) << " " << A(i,4) << " " << A(i,5) << " " << b(i) << std::endl;
	}
	*/
	calc_ordermatrix( nrow, ORDERSIZE, A, b, x, weights, reject );
	if( reject ) std::cout << "SET SCORE VALUE, FIX THIS LATER " << std::endl;
	calc_orderparam( x, vec, Azz, eta );
	Real score( calc_dipscore( A, x, b, All_RDC_lines, ORDERSIZE, Azz )*nrow );

	return score;

}

//////////////////////////////////////////////////////
//@brief
//////////////////////////////////////////////////////
void ResidualDipolarCouplingEnergy_Rohl::assemble_datamatrix(
	pose::Pose const & pose,
	utility::vector1< core::scoring::RDC_Rohl> const & All_RDC_lines,
	ObjexxFCL::FArray2D< Real > & A,
	ObjexxFCL::FArray1D< Real > & b,
	ObjexxFCL::FArray1D< Real > & weights
) const
{


	numeric::xyzVector< Real > umn;
	utility::vector1< core::scoring::RDC_Rohl >::const_iterator it;
	Size nrow( 0 );

	for( it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it) {

		++nrow;
			umn = pose.residue(it->res()).atom("N").xyz() - pose.residue(it->res()).atom("H").xyz();

			Real umn_x = umn.x()/it->fixed_dist();
			Real umn_y = umn.y()/it->fixed_dist();
			Real umn_z = umn.z()/it->fixed_dist();

			//filling matrix A
			A( nrow, 1 ) = umn_y*umn_y - umn_x*umn_x;
			A( nrow, 2 ) = umn_z*umn_z - umn_x*umn_x;
			A( nrow, 3 ) = 2.0*umn_x*umn_y;
			A( nrow, 4 ) = 2.0*umn_x*umn_z;
			A( nrow, 5 ) = 2.0*umn_z*umn_y;

			//filling matrix b
			b( nrow ) = it->Reduced_Jdipolar();
			weights( nrow ) = it->weight();
			//			std::cout << "nrow " << nrow << std::endl;
			//		std::cout << "Matrix A " << A(nrow,1) << " " << A(nrow,2) << " " << A(nrow,3) << " " << A(nrow,4) << " " << A(nrow,5) << std::endl;
	}


}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void ResidualDipolarCouplingEnergy_Rohl::calc_ordermatrix(
		Size const & nrow,
		Size const & ORDERSIZE,
		ObjexxFCL::FArray2D< Real > & A,
		ObjexxFCL::FArray1D< Real > & b,
		ObjexxFCL::FArray1D< Real > & x,
		ObjexxFCL::FArray1D< Real > & weights,
		bool & reject
) const
{


	core::Real factor = { 1e-6 }; // cutoff factor for singular values in svd

	ObjexxFCL::FArray2D< Real > U( nrow, ORDERSIZE );
	ObjexxFCL::FArray1D< Real > w( ORDERSIZE ); // singular values
	ObjexxFCL::FArray2D< Real > v( ORDERSIZE, ORDERSIZE );
	ObjexxFCL::FArray1D< Real > bweighted( nrow ); // singular values
	Real wmin, wmax; // min and max values of w
	Size sing; // number of singular values in w
	Real Sxx;

	// why not	U=A; ? should even be faster!
	Size ct_align = 0;
 	for( core::Size i = 1; i <= nrow; ++i ) { // copy A
		if ( weights( i ) > 0.000001 ) {
			++ct_align;
			U( ct_align, 1 ) = A(i,1) * weights( i ); // copy into U
			U( ct_align, 2 ) = A(i,2) * weights( i );
			U( ct_align, 3 ) = A(i,3) * weights( i );
			U( ct_align, 4 ) = A(i,4) * weights( i );
			U( ct_align, 5 ) = A(i,5) * weights( i );
			bweighted( ct_align  ) = b( i ) * weights( i );
		}
 	}

	svdcmp( U, ct_align, ORDERSIZE, w, v );

	wmax = 0.0;
	for ( core::Size j = 1; j <= ORDERSIZE; ++j ) {
		if ( w(j) > wmax ) wmax = w(j);
	}
	wmin = wmax * factor;
	sing = 0;
	for ( core::Size j = 1; j <= ORDERSIZE; ++j ) {
		if ( w(j) < wmin ) {
			w(j) = 0.0;
			++sing;
		}
	}
	if ( (int)sing > std::abs( int( ct_align ) - int( ORDERSIZE ) ) )
	 std::cout << "SVD yielded a matrix singular above expectation " <<
	 "in get_ordermatrix" << std::endl;

	// find solution for exact dipolar values

	svbksb( U, w, v, ct_align, ORDERSIZE, bweighted, x );

// x components: (Syy,Szz,Sxy,Sxz,Syz)
// check for acceptable values

	reject = false;
	Sxx = -x(1) - x(2);
	if ( Sxx < -0.5 || Sxx > 1.0 ) reject = true; // Sxx
	if ( x(1) < -0.5 || x(1) > 1.0 ) reject = true; // Syy
	if ( x(2) < -0.5 || x(2) > 1.0 ) reject = true; // Szz
	if ( x(3) < -0.75 || x(3) > 0.75 ) reject = true; // Sxy
	if ( x(4) < -0.75 || x(4) > 0.75 ) reject = true; // Sxz
	if ( x(5) < -0.75 || x(5) > 0.75 ) reject = true; // Syz

	if ( reject ) {
		std::cout << "order matrix not physically meaningful" << std::endl;
//		try with errors on dipolar values? map error?
//		score = 0.0;
		return;
	}



}


////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void ResidualDipolarCouplingEnergy_Rohl::svdcmp(
	ObjexxFCL::FArray2D< Real > & a,
	Size const & m,
	Size const & n,
	ObjexxFCL::FArray1D< Real > & w,
	ObjexxFCL::FArray2D< Real > & v
) const
{

//U    USES pythag
	Size i,its,j,jj,k,l,nm;
	ObjexxFCL::FArray1D< core::Real > rv1( n );
	Real anorm, c, f, g, h, s, scale, x, y, z;
	g = 0.0;
	scale = 0.0;
	anorm = 0.0;
	l = 0;
	nm = 0;
	for ( i = 1; i <= n; ++i ) {
		l = i+1;
		rv1(i) = scale*g;
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if ( i <= m ) {
			for ( k = i; k <= m; ++k ) {
				scale += std::abs(a(k,i));
			}
			if ( scale != 0.0 ) {
				for ( k = i; k <= m; ++k ) {
					a(k,i) /= scale;
					s += a(k,i)*a(k,i);
				}
				f = a(i,i);
				g = -sign(std::sqrt(s),f);
				h = f*g-s;
				a(i,i) = f-g;
				for ( j = l; j <= n; ++j ) {
					s = 0.0;
					for ( k = i; k <= m; ++k ) {
						s += a(k,i)*a(k,j);
					}
					f = s/h;
					for ( k = i; k <= m; ++k ) {
						a(k,j) += f*a(k,i);
					}
				}
				for ( k = i; k <= m; ++k ) {
					a(k,i) *= scale;
				}
			}
		}
		w(i) = scale *g;
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		if ( (i <= m) && (i != n) ) {
			for ( k = l; k <= n; ++k ) {
				scale += std::abs(a(i,k));
			}
			if ( scale != 0.0 ) {
				for ( k = l; k <= n; ++k ) {
					a(i,k) /= scale;
					s += a(i,k)*a(i,k);
				}
				f = a(i,l);
				g = -sign(std::sqrt(s),f);
				h = f*g-s;
				a(i,l) = f-g;
				for ( k = l; k <= n; ++k ) {
					rv1(k) = a(i,k)/h;
				}
				for ( j = l; j <= m; ++j ) {
					s = 0.0;
					for ( k = l; k <= n; ++k ) {
						s += a(j,k)*a(i,k);
					}
					for ( k = l; k <= n; ++k ) {
						a(j,k) += s*rv1(k);
					}
				}
				for ( k = l; k <= n; ++k ) {
					a(i,k) *= scale;
				}
			}
		}
		anorm = std::max(anorm,(std::abs(w(i))+std::abs(rv1(i))));
	}
	for ( i = n; i >= 1; --i ) {
		if ( i < n ) {
			if ( g != 0.0 ) {
				for ( j = l; j <= n; ++j ) {
					v(j,i) = (a(i,j)/a(i,l))/g;
				}
				for ( j = l; j <= n; ++j ) {
					s = 0.0;
					for ( k = l; k <= n; ++k ) {
						s += a(i,k)*v(k,j);
					}
					for ( k = l; k <= n; ++k ) {
						v(k,j) += s*v(k,i);
					}
				}
			}
			for ( j = l; j <= n; ++j ) {
				v(i,j) = 0.0;
				v(j,i) = 0.0;
			}
		}
		v(i,i) = 1.0;
		g = rv1(i);
		l = i;
	}
	for ( i = std::min(m,n); i >= 1; --i ) {
		l = i+1;
		g = w(i);
		for ( j = l; j <= n; ++j ) {
			a(i,j) = 0.0;
		}
		if ( g != 0.0 ) {
			g = 1.0/g;
			for ( j = l; j <= n; ++j ) {
				s = 0.0;
				for ( k = l; k <= m; ++k ) {
					s += a(k,i)*a(k,j);
				}
				f = (s/a(i,i))*g;
				for ( k = i; k <= m; ++k ) {
					a(k,j) += f*a(k,i);
				}
			}
			for ( j = i; j <= m; ++j ) {
				a(j,i) *= g;
			}
		} else {
			for ( j = i; j <= m; ++j ) {
				a(j,i) = 0.0;
			}
		}
		a(i,i) += 1.0;
	}
	for ( k = n; k >= 1; --k ) {
		for ( its = 1; its <= 30; ++its ) {
			for ( l = k; l >= 1; --l ) {
				nm = l-1;
				if ( (std::abs(rv1(l))+anorm) == anorm ) goto L2;
				if ( (std::abs(w(nm))+anorm) == anorm ) goto L1;
			}
L1:
			c = 0.0;
			s = 1.0;
			for ( i = l; i <= k; ++i ) {
				f = s*rv1(i);
				rv1(i) *= c;
				if ( (std::abs(f)+anorm) == anorm ) goto L2;
				g = w(i);
				h = pythag(f,g);
				w(i) = h;
				h = 1.0/h;
				c = (g*h);
				s = -(f*h);
				for ( j = 1; j <= m; ++j ) {
					y = a(j,nm);
					z = a(j,i);
					a(j,nm) = (y*c)+(z*s);
					a(j,i) = -(y*s)+(z*c);
				}
			}
L2:
			z = w(k);
			if ( l == k ) {
				if ( z < 0.0 ) {
					w(k) = -z;
					for ( j = 1; j <= n; ++j ) {
						v(j,k) = -v(j,k);
					}
				}
				goto L3;
			}
			if ( its == 30) utility_exit_with_message("no convergence in svdcmp \n" );
			x = w(l);
			nm = k-1;
			y = w(nm);
			g = rv1(nm);
			h = rv1(k);
			f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g = pythag(f,1.0);
			f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x;
			c = 1.0;
			s = 1.0;
			for ( j = l; j <= nm; ++j ) {
				i = j+1;
				g = rv1(i);
				y = w(i);
				h = s*g;
				g *= c;
				z = pythag(f,h);
				rv1(j) = z;
				c = f/z;
				s = h/z;
				f = (x*c)+(g*s);
				g = -(x*s)+(g*c);
				h = y*s;
				y *= c;
				for ( jj = 1; jj <= n; ++jj ) {
					x = v(jj,j);
					z = v(jj,i);
					v(jj,j) = (x*c)+(z*s);
					v(jj,i) = -(x*s)+(z*c);
				}
				z = pythag(f,h);
				w(j) = z;
				if ( z != 0.0 ) {
					z = 1.0/z;
					c = f*z;
					s = h*z;
				}
				f = (c*g)+(s*y);
				x = -(s*g)+(c*y);
				for ( jj = 1; jj <= m; ++jj ) {
					y = a(jj,j);
					z = a(jj,i);
					a(jj,j) = (y*c)+(z*s);
					a(jj,i) = -(y*s)+(z*c);
				}
			}
			rv1(l) = 0.0;
			rv1(k) = f;
			w(k) = x;
		}
L3:;
	}


}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
Real ResidualDipolarCouplingEnergy_Rohl::pythag(
	Real const & a,
	Real const & b
) const
{

	Real pythag;

	Real absa = std::abs(a);
	Real absb = std::abs(b);
	if ( absa > absb ) {
		Real const ratio = absb/absa;
		pythag = absa * std::sqrt( 1.0 + ( ratio * ratio ) );
	} else {
		if ( absb == 0.0 ) {
			pythag = 0.0;
		} else {
			Real const ratio = absa/absb;
			pythag = absb * std::sqrt( 1.0 + ( ratio * ratio ) );
		}
	}
	return pythag;


}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void ResidualDipolarCouplingEnergy_Rohl::svbksb(
	ObjexxFCL::FArray2D< Real > const & u,
	ObjexxFCL::FArray1D< Real > const & w,
	ObjexxFCL::FArray2D< Real > const & v,
	Size const & m,
	Size const & n,
	ObjexxFCL::FArray1D< Real > const & b,
	ObjexxFCL::FArray1D< Real > & x
) const
{


	ObjexxFCL::FArray1D< core::Real > tmp( n );
	Real s;

	for ( Size j = 1; j <= n; ++j ) {
		s = 0.0;
		if ( w(j) != 0.0 ) {
			for ( Size i = 1; i <= m; ++i ) {
				s += u(i,j) * b(i);
			}
			s /= w(j);
		}
		tmp(j) = s;
	}
	for ( Size j = 1; j <= n; ++j ) {
		s = 0.0;
		for ( Size jj = 1; jj <= n; ++jj ) {
			s += v(j,jj) * tmp(jj);
		}
		x(j) = s;
	}


}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
void ResidualDipolarCouplingEnergy_Rohl::calc_orderparam(
	ObjexxFCL::FArray1D< Real > x,
	ObjexxFCL::FArray2D< Real > vec,
	Real & Azz,
	Real & eta
) const
{

	using numeric::xyzMatrix;

	ObjexxFCL::FArray1D< Size > sort( 3 ); // sorted index to val, val(sort(1)) = largest abs val.
	Real temp1, temp2;

	// Assemble order matrix
	numeric::xyzMatrix< Real > S = numeric::xyzMatrix< core::Real >::rows(	-x(1) - x(2), x(3), x(4),x(3), x(1), x(5), x(4), x(5), x(2) );


	numeric::xyzVector< Real > val; // Eigenvalues
	numeric::xyzMatrix< Real > xyz_vec; // Eigenvectors

	// Find eigenvalues, eigenvectors of symmetric matrix S
	val = eigenvector_jacobi( S, 1E-9, xyz_vec );

	// sort eigenvalues
	sort(1) = 1;
	sort(2) = 2;
	sort(3) = 3;

	if ( std::abs(val(1)) < std::abs(val(2)) ) {     // largest absolute value
		sort(2) = 1;
		sort(1) = 2;
	}
	if ( std::abs(val(sort(2))) < std::abs(val(3)) ) {
		sort(3) = sort(2);
		sort(2) = 3;
		if ( std::abs(val(sort(1))) < std::abs(val(3)) ) {
			sort(2) = sort(1);
			sort(1) = 3;
		}
	}


	Azz = val(sort(1));
	eta = (2.0/3.0) * std::abs(val(sort(2))-val(sort(3))/Azz);

// sort eigen values      // largest to smallest : Azz,Ayy,Axx
	temp1 = val(sort(1));
	temp2 = val(sort(2));
	val(3) = val(sort(3));
	val(2) = temp2;
	val(1) = temp1;

// sort eigen vectors
	for ( Size i = 1; i <= 3; ++i ) {
		temp1 = xyz_vec(i,sort(3));
		temp2 = xyz_vec(i,sort(2));
		vec(i,3) = xyz_vec(i,sort(1));
		vec(i,1) = temp1;
		vec(i,2) = temp2;
	}

	Azz = val(1);
	eta = (2.0/3.0) * (val(3)-val(2))/val(1);



}
////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
Real ResidualDipolarCouplingEnergy_Rohl::calc_dipscore(
	ObjexxFCL::FArray2D< Real > const & A,
	ObjexxFCL::FArray1D< Real > const & x,
	ObjexxFCL::FArray1D< Real > const & b,
	utility::vector1< core::scoring::RDC_Rohl > const & All_RDC_lines,
	Size const & ORDERSIZE,
	Real const & Azz
) const
{

	assert( Azz != 0 );

	Real score( 0.0 );

	utility::vector1< core::scoring::RDC_Rohl >::const_iterator it;
	Size nrow( 0 );
	for( it = All_RDC_lines.begin(); it != All_RDC_lines.end(); ++it) {
		Real Jcalc( 0.0 );
		++nrow;
		for( Size j = 1; j <= ORDERSIZE; ++j ) {
			Jcalc += A( nrow, j )*x(j);
		}
		score += ( b( nrow ) -Jcalc )*( b( nrow ) - Jcalc ); //these are reduced Jd

		//v		std::cout << b(nrow)/it->invDcnst() << " " << Jcalc/it->invDcnst() << " " << numeric::square( b(nrow)/it->invDcnst() - Jcalc/it->invDcnst() ) << " " << it->res() << std::endl;
	}

	//std::cout << "score, nrow, Azz " << score << " " << nrow << " " << Azz << std::endl;

	score /= ( nrow*( Azz*Azz ) );

	//std::cout << "Total score " << score << std::endl;

	return score;

}
core::Size
ResidualDipolarCouplingEnergy_Rohl::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core
