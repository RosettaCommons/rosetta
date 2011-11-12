// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/DihedralConstraint.cc
/// @brief

#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/FuncFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/deriv/dihedral_deriv.hh>

#include <utility/exit.hh>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
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
#include <core/chemical/types.hh>
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
#include <core/conformation/Residue.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/MiniPose.fwd.hh>
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
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintCreator.fwd.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <utility/down_cast.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/factory/WidgetRegistrator.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
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
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/sphericalVector.hh>
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <numeric/random/random.fwd.hh>
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
#include <ObjexxFCL/proxy_const_assert.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
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
#include <typeinfo>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Auto using namespaces
namespace std { } using namespace std; // AUTO USING NS
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end

namespace core {
namespace scoring {
namespace constraints {

static basic::Tracer TR("core.io.constraints");

/////////////////////////////////////////////////////////////////////////////
// contributions for a term that looks like
//
//   d(M-F)
// ( ------ x v ) dot w
//   d phi
//
// F1 collects terms f1 that look like: (-u_phi) dot f1
// F2 collects terms f2 that look like: (-u_phi x R_phi) dot f2
//
//
// where u_phi is the torsion axis of rotation
// and R_phi is a terminal atom of this axis
//
// static class function
void
DihedralConstraint::helper(
	Vector const & M,
	Vector const & v,
	Vector const & w,
	Vector & F1,
	Vector & F2
)
{
	Vector const f2 = cross(v,w);
	F2 += f2;
	F1 += cross(f2,M);
}

/////////////////////////////////////////////////////////////////////////////
//
// static class function
void
DihedralConstraint::p1_cosine_deriv(
	Vector const & p1,
	Vector const & p2,
	Vector const & p3,
	Vector const & p4,
	Real & x,
	Vector & F1,
	Vector & F2
)
{
	F1 = 0.0;
	F2 = 0.0;

	Vector v1( p1-p2 );
	Vector v2( p2-p3 );
	Vector v3( p3-p4 );

	Vector v12( cross( v1, v2 ));
	Vector v23( cross( v2, v3 ));

	Real const n12( v12.length() );
	Real const n23( v23.length() );

	if ( n12 < 1e-9 || n23 < 1e-9 ) return;

	x = dot( v12, v23) / ( n12 * n23 );

	// first term:
	{
		Real const f( 1.0f/ ( n12 * n23 ) );
		helper( p1, f * v2, v23 , F1, F2);
	}

	// second term
	{
		Real const f( -1.0f * x / ( n12 * n12 ) );
		helper( p1, f * v2, v12, F1, F2 );
	}


	// debugging
	// translation of p1 in the place spanned by v1 and v2
	// does not affect the torsion angle
	// ==> rotation of p1 about an axis perpendicular to this plane
	// also does not change the torsion angle, ie deriv should be 0
	assert( std::abs( dot( F2, v1 ) ) < 1e-3 );
	assert( std::abs( dot( F2, v2 ) ) < 1e-3 );
	assert( std::abs( dot( F1, cross( v1, v2 ) ) ) < 1e-3 );
}

/////////////////////////////////////////////////////////////////////////////
// static class function
//
void
DihedralConstraint::p2_cosine_deriv(
	Vector const & p1,
	Vector const & p2,
	Vector const & p3,
	Vector const & p4,
	Real & x,
	Vector & F1,
	Vector & F2
)
{
	//std::cout << "p2_cosine_deriv!" << std::endl;

	F1 = 0.0;
	F2 = 0.0;

	typedef Vector Vec;
	Vector v1( p1-p2 );
	Vector v2( p2-p3 );
	Vector v3( p3-p4 );

	Vector v12( cross( v1, v2 ));
	Vector v23( cross( v2, v3 ));

	Real const n12( v12.length() );
	Real const n23( v23.length() );

	if ( n12 < 1e-9 || n23 < 1e-9 ) return;

	x = dot( v12, v23) / ( n12 * n23 );

	// here we are taking derivatives of an expression that looks like
	//
	//                   dot( v12, v23 )
	// x = cos theta =  -----------------
	//                      n12 * n23
	//
	// where theta is our dihedral angle
	//

	{ // derivatives of the numerator
		// v1 and v2 both depend on position of p2

		{ // first term
			Real const f( -1.0f/ ( n12 * n23 ) );
			helper( p2, f * v2, v23 , F1, F2);
		}

		{ // second term
			Real const f( -1.0f/ ( n12 * n23 ) );
			helper( p2, f * v1, v23 , F1, F2);
		}

		{ // third term
			Real const f( 1.0f/ ( n12 * n23 ) );
			helper( p2, f * v3, v12 , F1, F2);
		}
	}

	{ // derivatives of the denominator
		// v1 and v2 both depend on position of p2

		{ // first term
			Real const f( x/ ( n12 * n12 ) );
			helper( p2, f * v2, v12, F1, F2 );
		}

		{ // second term
			Real const f( x/ ( n12 * n12 ) );
			helper( p2, f * v1, v12, F1, F2 );
		}

		{ // third term
			Real const f( -1.0f * x/ ( n23 * n23 ) );
			helper( p2, f * v3, v23, F1, F2 );
		}
	}

	// debugging
	// translation of p2 along v2 does not change the torsion angle
	assert( std::abs( dot( F2, v2 ) ) < 1e-3 );

}

/////////////////////////////////////////////////////////////////////////////
Real
DihedralConstraint::score(
	conformation::Conformation const & conformation
) const
{
	return score(  conformation.xyz( atom1_ ),conformation.xyz( atom2_ ),
								 conformation.xyz( atom3_ ),conformation.xyz( atom4_ ) );
}


/////////////////////////////////////////////////////////////////////////////
void
DihedralConstraint::fill_f1_f2(
	AtomID const & atom,
	XYZ_Func const & xyz,
	Vector & F1,
 	Vector & F2,
 	EnergyMap const & weights
) const
{


	/*using numeric::conversions::radians;

	//static bool dihedral_deriv( true );
	//if ( dihedral_deriv ) { std::cout << "DIHEDRAL DERIV" << std::endl; }
	//dihedral_deriv = false;

	// to avoid problems with dtheta/dx around 0 and 180 degrees
	// truncate x a bit in the calculation of the derivative
	static Real const small_angle( radians( Real(0.1) ) );
	static Real const big_angle( radians( Real(179.9) ) );
	static Real const max_x( std::cos( small_angle ));
	static Real const min_x( std::cos( big_angle ));
	// dtheta_dx has a value of ~ 572.96 for min_x and max_x
	// this goes to infinity as x goes to -1 or 1

	Vector f1(0.0) ,f2(0.0);

	Real x(0.0); // cos theta
	if ( atom == atom1_ ) {
		p1_cosine_deriv( conformation.xyz( atom1_ ),conformation.xyz( atom2_ ),
										 conformation.xyz( atom3_ ),conformation.xyz( atom4_ ), x, f1, f2 );
	} else if ( atom == atom2_ ) {
		p2_cosine_deriv( conformation.xyz( atom1_ ),conformation.xyz( atom2_ ),
										 conformation.xyz( atom3_ ),conformation.xyz( atom4_ ), x, f1, f2 );
	} else if ( atom == atom3_ ) {
		p2_cosine_deriv( conformation.xyz( atom4_ ),conformation.xyz( atom3_ ),
										 conformation.xyz( atom2_ ),conformation.xyz( atom1_ ), x, f1, f2 );
	} else if ( atom == atom4_ ) {
		p1_cosine_deriv( conformation.xyz( atom4_ ),conformation.xyz( atom3_ ),
										 conformation.xyz( atom2_ ),conformation.xyz( atom1_ ), x, f1, f2 );
	} else {
		return;
	}

	Real const thetaU( numeric::arccos( x )); // unsigned version of theta

	Real const theta ( dihedral_radians( conformation.xyz( atom1_ ),conformation.xyz( atom2_ ),
																			 conformation.xyz( atom3_ ),conformation.xyz( atom4_ ) ) );

	assert( std::abs( std::abs( theta ) - thetaU ) < 1e-2 );

	Real const dE_dtheta( dfunc( theta ) );

	x = std::min( std::max( min_x, x ), max_x );
	Real const dthetaU_dx = -1 / sqrt( 1- x*x );
	Real const dtheta_dthetaU( theta < 0 ? -1 : 1 );

	F1 += ( dE_dtheta * dtheta_dthetaU * dthetaU_dx * f1 ) * weights[ this->score_type() ];

	F2 += ( dE_dtheta * dtheta_dthetaU * dthetaU_dx * f2 ) * weights[ this->score_type() ];

*//**/

	using namespace numeric::deriv;

	Vector f1(0.0) ,f2(0.0);

	Real theta(0.0);
	if ( atom == atom1_ ) {
		dihedral_p1_cosine_deriv( xyz( atom1_ ),xyz( atom2_ ), xyz( atom3_ ), xyz( atom4_ ), theta, f1, f2 );
	} else if ( atom == atom2_ ) {
		dihedral_p2_cosine_deriv( xyz( atom1_ ), xyz( atom2_ ), xyz( atom3_ ), xyz( atom4_ ), theta, f1, f2 );
	} else if ( atom == atom3_ ) {
		dihedral_p2_cosine_deriv( xyz( atom4_ ), xyz( atom3_ ), xyz( atom2_ ), xyz( atom1_ ), theta, f1, f2 );
	} else if ( atom == atom4_ ) {
		dihedral_p1_cosine_deriv( xyz( atom4_ ), xyz( atom3_ ), xyz( atom2_ ), xyz( atom1_ ), theta, f1, f2 );
	} else {
		return;
	}

	Real const dE_dtheta( dfunc( theta ) );

	F1 += dE_dtheta * weights[ this->score_type() ] * f1;
	F2 += dE_dtheta * weights[ this->score_type() ] * f2;
/**/
}

ConstraintOP
DihedralConstraint::remap_resid(
	core::id::SequenceMapping const & seqmap
) const {
  if ( seqmap[atom1_.rsd()] != 0 && seqmap[atom2_.rsd()] != 0 && seqmap[atom3_.rsd()] != 0 && seqmap[atom4_.rsd()] != 0 ) {
    AtomID remap_a1( atom1_.atomno(), seqmap[atom1_.rsd()] ),
      remap_a2( atom2_.atomno(), seqmap[atom2_.rsd()] ),
			remap_a3( atom3_.atomno(), seqmap[atom3_.rsd()] ),
			remap_a4( atom4_.atomno(), seqmap[atom4_.rsd()] );
    return ConstraintOP( new DihedralConstraint( remap_a1, remap_a2, remap_a3, remap_a4, this->func_ ) );
  } else {
    return NULL;
  }
}


/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP DihedralConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
  id::NamedAtomID atom1( core::pose::atom_id_to_named_atom_id(atom(1), src ) );
  id::NamedAtomID atom2( core::pose::atom_id_to_named_atom_id(atom(2), src ) );
  id::NamedAtomID atom3( core::pose::atom_id_to_named_atom_id(atom(3), src ) );
  id::NamedAtomID atom4( core::pose::atom_id_to_named_atom_id(atom(4), src ) );
//  std::cout << "making a clone";
  if ( smap ) {
    atom1.rsd() = (*smap)[ atom1_.rsd() ];
    atom2.rsd() = (*smap)[ atom2_.rsd() ];
    atom3.rsd() = (*smap)[ atom3_.rsd() ];
    atom4.rsd() = (*smap)[ atom4_.rsd() ];
  }

  //get AtomIDs for target pose
  id::AtomID id1( core::pose::named_atom_id_to_atom_id(atom1, dest ));
  id::AtomID id2( core::pose::named_atom_id_to_atom_id(atom2, dest ));
  id::AtomID id3( core::pose::named_atom_id_to_atom_id(atom3, dest ));
  id::AtomID id4( core::pose::named_atom_id_to_atom_id(atom4, dest ));
  if ( id1.valid() && id2.valid() &&  id3.valid() && id4.valid()  ) {
    return new DihedralConstraint( id1, id2, id3, id4, func_, score_type() );
  } else {
    return NULL;
  }
}


id::AtomID const &
DihedralConstraint::atom( Size const n ) const {
	switch( n ) {
	case 1:
		return atom1_;
	case 2:
		return atom2_;
	case 3:
		return atom3_;
	case 4:
		return atom4_;
	default:
		utility_exit_with_message( "DihedralConstraint::atom() bad argument" );
	}
	return atom1_;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
///@details one line definition "Dihedral atom1 res1 atom2 res2 atom3 res3 atom4 res4 function_type function_definition"
void
DihedralConstraint::read_def(
	std::istream & in,
	pose::Pose const & pose,
	FuncFactory const & func_factory
) {
	Size res1, res2, res3, res4;
	std::string tempres1, tempres2, tempres3, tempres4;
	std::string name1, name2, name3, name4;
	std::string func_type;
	std::string type;

	in
		>> name1 >> tempres1
		>> name2 >> tempres2
		>> name3 >> tempres3
		>> name4 >> tempres4
		>> func_type;

	ConstraintIO::parse_residue( pose, tempres1, res1 );
	ConstraintIO::parse_residue( pose, tempres2, res2 );
	ConstraintIO::parse_residue( pose, tempres3, res3 );
	ConstraintIO::parse_residue( pose, tempres4, res4 );

	TR.Debug 	<< "read: " << name1 << " " << name2 << " " << name3 << " " << name4 << " "
						<< res1 << " " << res2 << " " << res3 << " " << res4 << " func: " << func_type
						<< std::endl;
	if ( res1 > pose.total_residue() || res2 > pose.total_residue() || res3 > pose.total_residue() || res4 > pose.total_residue() ) {
		TR.Warning 	<< "ignored constraint (no such atom in pose!)"
								<< name1 << " " << name2 << " " << name3 << " " << name4 << " "
													<< res1 << " " << res2 << " " << res3 << " " << res4 << " func: " << func_type
													<< std::endl;
		in.setstate( std::ios_base::failbit );
		return;
	}

	atom1_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name1, res1 ), pose ) );
	atom2_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name2, res2 ), pose ) );
	atom3_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name3, res3 ), pose ) );
	atom4_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name4, res4 ), pose ) );

	if ( atom1_.atomno() == 0 || atom2_.atomno() == 0 || atom3_.atomno() == 0 || atom4_.atomno() == 0 ) {
		TR.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "," << name3 << "," << name4 << "), "
			<< "and found AtomIDs (" << atom1_ << "," << atom2_ << "," << atom3_ << "," << atom4_ << ")"
			<< std::endl;
			in.setstate( std::ios_base::failbit );
			return;
	}

	func_ = func_factory.new_func( func_type );
	func_->read_data( in );

	//chu skip the rest of line since this is a single line defintion.
	while( in.good() && (in.get() != '\n') ) {}

	if ( TR.Debug.visible() ) {
		func_->show_definition( std::cout );
		std::cout << std::endl;
	}
}



/////////////////////////////////////////////////////////////////////////////
Real
DihedralConstraint::score(
	Vector const & p1,
	Vector const & p2,
	Vector const & p3,
	Vector const & p4
) const
{
	return func( dihedral_radians( p1, p2, p3, p4 ) );
}


bool
DihedralConstraint::operator == ( Constraint const & other_cst ) const
{
	if( !dynamic_cast< DihedralConstraint const * > ( &other_cst ) ) return false;

	DihedralConstraint const & other( static_cast< DihedralConstraint const & > (other_cst) );

	if( atom1_ != other.atom1_ ) return false;
	if( atom2_ != other.atom2_ ) return false;
	if( atom3_ != other.atom3_ ) return false;
	if( atom4_ != other.atom4_ ) return false;
	if( func_ != other.func_ ) return false;
	if( this->score_type() != other.score_type() ) return false;

	return true;
}

void
DihedralConstraint::score(
	XYZ_Func const & xyz,
	EnergyMap const &,
	EnergyMap & emap
) const {
	emap[ this->score_type() ] += score(
		xyz( atom1_ ), xyz( atom2_ ), xyz( atom3_ ), xyz( atom4_ )
	);
}

void DihedralConstraint::show( std::ostream & out ) const {
	out << "DihedralConstraint";
	for ( Size i = 1; i <= natoms(); ++i ) {
		AtomID const & id = atom(i);
		out << ' ' << id.rsd() << ' ' << id.atomno();
	}
	out << ' ';
	func_->show_definition(out);
}


Size DihedralConstraint::show_violations(
        std::ostream& out,
        pose::Pose const& pose,
        Size verbose_level,
        Real threshold
) const {

	     if ( verbose_level > 80 ) {
                out << "Dihedral ("
                        << pose.residue_type(atom1_.rsd() ).atom_name( atom1_.atomno() ) << ":"
                        << atom1_.atomno() << "," << atom1_.rsd() << "-"
                        << pose.residue_type(atom2_.rsd() ).atom_name( atom2_.atomno() ) << ":"
                        << atom2_.atomno() << "," << atom2_.rsd() << ") ";
        }
        if ( verbose_level > 120 ) { //don't ask but I had a really weird bug to track down!
                conformation::Conformation const & conformation( pose.conformation() );
                Vector const & xyz1( conformation.xyz( atom1_ ) ), xyz2( conformation.xyz( atom2_ ) );
                out << "\ncoords1: " << xyz1[ 1 ] << " " << xyz1[ 2 ] << " " << xyz1[ 3 ] << " --- ";
                out << "coords1: " << xyz2[ 1 ] << " " << xyz2[ 2 ] << " " << xyz2[ 3 ] << "\n";
        }

        return func_->show_violations( out, 0.0, verbose_level, threshold );
}


Real
DihedralConstraint::func( Real const theta ) const {
	return func_->func( theta );
}

Real
DihedralConstraint::dfunc( Real const theta ) const {
	return func_->dfunc( theta );
}

} // constraints
} // scoring
} // core
