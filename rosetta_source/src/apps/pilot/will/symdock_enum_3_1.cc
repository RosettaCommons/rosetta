// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file /src/apps/pilat/will/symdock_enum.cc
/// @brief docks trimer center and pentamer center together

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmData.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/sasa.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/string_util.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

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
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueSelector.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
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
// AUTO-REMOVED #include <core/conformation/symmetry/SymDof.fwd.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymSlideInfo.fwd.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymSlideInfo.hh>
#include <core/conformation/symmetry/SymmData.fwd.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/VirtualCoordinate.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
// AUTO-REMOVED #include <core/id/DOF_ID_Map.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/io/pdb/file_data.fwd.hh>
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
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/vector0_bool.hh>
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
#include <utility/file/gzip_util.hh>
// AUTO-REMOVED #include <utility/io/irstream.fwd.hh>
// AUTO-REMOVED #include <utility/io/irstream.hh>
// AUTO-REMOVED #include <utility/io/izstream.fwd.hh>
#include <utility/io/mpistream.hh>
#include <utility/io/orstream.fwd.hh>
#include <utility/io/orstream.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/io/zipstream.hpp>
#include <utility/io/zipstream.ipp>
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
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/IOTraits.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/sphericalVector.hh>
#include <numeric/trig.functions.hh>
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzMatrix.io.hh>
#include <numeric/xyzTriple.fwd.hh>
#include <numeric/xyzTriple.hh>
#include <numeric/xyzTriple.io.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.fwd.hh>
#include <numeric/random/uniform.hh>
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
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/Fstring.fwd.hh>
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
#include <ObjexxFCL/byte.fwd.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
// AUTO-REMOVED #include <ObjexxFCL/ubyte.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <execinfo.h>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/BasicDataCache.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableData.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/CacheableData.hh>
// AUTO-REMOVED #include <basic/datacache/DataCache.fwd.hh>
// AUTO-REMOVED #include <basic/datacache/DataCache.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <boost/algorithm/string/erase.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>
#include <zlib/zlib.h>
#include <zlib/zutil.h>


static basic::Tracer TR("symdock_enum");

using core::Size;
using core::Real;
using core::pose::Pose;
using utility::vector1;
using ObjexxFCL::fmt::I;
using ObjexxFCL::fmt::F;
using numeric::min;
using numeric::max;

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;
typedef numeric::xyzVector<double> Vecf;
typedef numeric::xyzMatrix<double> Matf;

void trans_pose( Pose & pose, Vec const & trans ) {
  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, pose.xyz(aid) + trans );
    }
  }
}

void rot_pose( Pose & pose, Mat const & rot ) {
  for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
    for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, rot * pose.xyz(aid) );
    }
  }
}

void rot_pose( Pose & pose, Mat const & rot, Vec const & cen ) {
  trans_pose(pose,-cen);
  rot_pose(pose,rot);
  trans_pose(pose,cen);
}

void rot_pose( Pose & pose, Vec const & axis, Real const & ang ) {
  rot_pose(pose,rotation_matrix_degrees(axis,ang));
}

void rot_pose( Pose & pose, Vec const & axis, Real const & ang, Vec const & cen ) {
  rot_pose(pose,rotation_matrix_degrees(axis,ang),cen);
}


void alignaxis(Pose & pose, Vec newaxis, Vec oldaxis, Vec cen = Vec(0,0,0) ) {
  newaxis.normalize();
  oldaxis.normalize();
  Vec axis = newaxis.cross(oldaxis).normalized();
  Real ang = -acos(numeric::max(-1.0,numeric::min(1.0,newaxis.dot(oldaxis))))*180/numeric::constants::d::pi;
  rot_pose(pose,axis,ang,cen);
}

int cbcount_vec(vector1<Vecf> & cba, vector1<Vecf> & cbb) {
  int cbcount = 0;
  for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia)
    for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib)
      if( ib->distance_squared(*ia) < 100.0 ) cbcount++;
  return cbcount;
}

void set_cb_pairs(vector1<Vecf> & cba, vector1<Vecf> & cbb) {
  vector1<Vecf> a,b;
  for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
    for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
      if( ib->distance_squared(*ia) < 100.0 ) {
        a.push_back(*ia);
        b.push_back(*ib);
      }
    }
  }
  cba = a;
  cbb = b;
}

int pose_cbcount(Pose const & a, Pose const & b) {
  int count = 0;
  for(Size i = 1; i <= a.n_residue(); ++i) {
    for(Size j = 1; j <= b.n_residue(); ++j) {
      if(a.residue(i).xyz(2).distance_squared(b.residue(j).xyz(2)) < 100.0) {
        count++;
      }
    }
  }
  return count;
}

double
sicfast(
        vector1<Vecf>  pa,
        vector1<Vecf>  pb,
        vector1<Vecf> & cba,
        vector1<Vecf> & cbb,
        Vecf ori,
        int & cbcount,
        bool debug = false
        ) {
  double BIN = 2.0;

  // get points, rotated ro ori is 0,0,1, might already be done
  Matf rot = Matf::identity();
  if     ( ori.dot(Vec(0,0,1)) < -0.99999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
  else if( ori.dot(Vec(0,0,1)) <  0.99999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
  if( rot != Matf::identity() ) {
    for(vector1<Vecf>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
    for(vector1<Vecf>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
  }

  // get bounds for plane hashes
  double xmx1=-9e9,xmn1=9e9,ymx1=-9e9,ymn1=9e9,xmx=-9e9,xmn=9e9,ymx=-9e9,ymn=9e9;
  for(vector1<Vecf>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
    xmx1 = max(xmx1,ia->x()); xmn1 = min(xmn1,ia->x());
    ymx1 = max(ymx1,ia->y()); ymn1 = min(ymn1,ia->y());
  }
  for(vector1<Vecf>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
    xmx = max(xmx,ib->x()); xmn = min(xmn,ib->x());
    ymx = max(ymx,ib->y()); ymn = min(ymn,ib->y());
  }
  xmx = min(xmx,xmx1); xmn = max(xmn,xmn1);
  ymx = min(ymx,ymx1); ymn = max(ymn,ymn1);


  int xlb = (int)floor(xmn/BIN)-2; int xub = (int)ceil(xmx/BIN)+2; // one extra on each side for correctness,
  int ylb = (int)floor(ymn/BIN)-2; int yub = (int)ceil(ymx/BIN)+2; // and one extra for outside atoms

  // TR << "BOUNDS " << xmn << " " << xmx << " " << ymn << " " << ymx << std::endl;
  // TR << "BOUNDS " << xlb << " " << xub << " " << ylb << " " << yub << std::endl;

  // insert points into hashes
  int const xsize = xub-xlb+1;
  int const ysize = yub-ylb+1;
  ObjexxFCL::FArray2D<Vecf> ha(xsize,ysize,Vecf(0,0,-9e9)),hb(xsize,ysize,Vecf(0,0,9e9));
  for(vector1<Vecf>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
    // int const ix = min(xsize,max(1,(int)ceil(ia->x()/BIN)-xlb));
    // int const iy = min(ysize,max(1,(int)ceil(ia->y()/BIN)-ylb));
    int const ix = (int)ceil(ia->x()/BIN)-xlb;
    int const iy = (int)ceil(ia->y()/BIN)-ylb;
    if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
    if( ha(ix,iy).z() < ia->z() ) ha(ix,iy) = *ia;
  }
  for(vector1<Vecf>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
    // int const ix = min(xsize,max(1,(int)ceil(ib->x()/BIN)-xlb));
    // int const iy = min(ysize,max(1,(int)ceil(ib->y()/BIN)-ylb));
    int const ix = (int)ceil(ib->x()/BIN)-xlb;
    int const iy = (int)ceil(ib->y()/BIN)-ylb;
    if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
    if( hb(ix,iy).z() > ib->z() ) hb(ix,iy) = *ib;
  }

  // check hashes for min dis
  int imna=0,jmna=0,imnb=0,jmnb=0;
  double mindis = 9e9;
  for(int i = 1; i <= xsize; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
    for(int j = 1; j <= ysize; ++j) {
      for(int k = -2; k <= 2; ++k) {
        if(i+k < 1 || i+k > xsize) continue;
        for(int l = -2; l <= 2; ++l) {
          if(j+l < 1 || j+l > ysize) continue;
          double const xa = ha(i  ,j  ).x();
          double const ya = ha(i  ,j  ).y();
          double const xb = hb(i+k,j+l).x();
          double const yb = hb(i+k,j+l).y();
          double const d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);

          if( d2 < 16.0 ) {
            double dz = hb(i+k,j+l).z() - ha(i,j).z() - sqrt(16.0-d2);
            if( dz < mindis ) {
              mindis = dz;
              imna = i;
              jmna = j;
              imnb = i+k;
              jmnb = j+l;
            }
          }
        }
      }
    }
  }

  // {
  //  utility::io::ozstream out("cba.pdb");
  //  for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
  //    Vec viz = (*ia) + (mindis*ori);
  //    out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"A"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
  //  }
  //  out.close();
  // }
  // {
  //  utility::io::ozstream out("cbb.pdb");
  //  for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
  //    Vec viz = (*ib);
  //    out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"B"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
  //  }
  //  out.close();
  // }

  cbcount = 0;
  // utility::io::ozstream out("cb8.pdb");
  // TR << "CB0 " << cbcount << std::endl;
  for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
    for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
      if( ib->distance_squared( (*ia) + (mindis*ori) ) < 100.0 ) {
        cbcount++;
        // Vec viz = (*ia) + (mindis*ori);
        // out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"A"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        // viz = *ib;
        // out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<  "VIZ"<<' '<<"B"<<I(4,100)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      }
    }
  }
  // out.close();
  // TR << "CB1 " << cbcount << std::endl;

  // // rotate points back -- needed iff pa/pb come by reference
  rot = rot.transposed();
  // if( rot != Matf::identity() ) {
  //  for(vector1<Vecf>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
  //  for(vector1<Vecf>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
  // }

  // uncomment this to get hashes in local space
  // rot = Matf::identity();
  // ori = Vec(0,0,1);

  if(debug){
    {
      utility::io::ozstream out("hasha.pdb");
      for(int i = 2; i <= xsize-1; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
        for(int j = 2; j <= ysize-1; ++j) {
          Vecf viz = rot*ha(i,j) + mindis*ori;
          if(viz.z() < -9e8 || 9e8 < viz.z()) continue;
          out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"B"<<I(4,100+j)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        }
      }
      Vecf viz = rot*ha(imna,jmna) + mindis*ori;
      out<<"HETATM"<<I(5,1000+imna)<<' '<<"MIN "<<' ' <<  "MIN"<<' '<<"B"<<I(4,100+jmna)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      out.close();
    }
    {
      utility::io::ozstream out("hashb.pdb");
      for(int i = 2; i <= xsize-1; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
        for(int j = 2; j <= ysize-1; ++j) {
          Vecf viz = rot*hb(i,j);
          if(viz.z() < -9e8 || 9e8 < viz.z()) continue;
          out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"C"<<I(4,100+j)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
        }
      }
      Vecf viz = rot*hb(imnb,jmnb);
      out<<"HETATM"<<I(5,1000+imnb)<<' '<<"MIN "<<' ' <<  "MIN"<<' '<<"C"<<I(4,100+jmnb)<<"    "<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
      out.close();
    }
  }

  return mindis;
}

double sicfast(
               Pose const & a,
               Pose const & b,
               Vecf ori_in,
               int & cbcount
               ) {
  // get points, rotated ro ori is 0,0,1
  vector1<Vecf> pa,pb;
  vector1<Vecf> cba,cbb;
  Vecf ori = ori_in.normalized();
  Matf rot = Matf::identity();
  if     ( ori.dot(Vec(0,0,1)) < -0.999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
  else if( ori.dot(Vec(0,0,1)) <  0.999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
  for(int i = 1; i <= (int)a.n_residue(); ++i) {
    cba.push_back(rot*Vecf(a.residue(i).xyz(2)));
    int const natom = (a.residue(i).name3()=="GLY") ? 4 : 5;
    for(int j = 1; j <= natom; ++j) pa.push_back(rot*Vecf(a.residue(i).xyz(j)));
  }
  for(int i = 1; i <= (int)b.n_residue(); ++i) {
    cbb.push_back(rot*Vecf(b.residue(i).xyz(2)));
    int const natom = (b.residue(i).name3()=="GLY") ? 4 : 5;
    for(int j = 1; j <= natom; ++j) pb.push_back(rot*Vecf(b.residue(i).xyz(j)));
  }
  return sicfast( pa, pb, cba, cbb, Vec(0,0,1), cbcount );
}

void make_trimer(core::pose::Pose & pose) {
	core::pose::Pose t2(pose),t3(pose);
	rot_pose(t2,Vec(0,0,1),120.0);
	rot_pose(t3,Vec(0,0,1),240.0);
	for(Size i = 1; i <= t2.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
	for(Size i = 1; i <= t3.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t3.residue(i),1); else pose.append_residue_by_bond(t3.residue(i));
}

void make_pentamer(core::pose::Pose & pose) {
	core::pose::Pose t2(pose),t3(pose),t4(pose),t5(pose);
	rot_pose(t2,Vec(0,0,1), 72.0);
	rot_pose(t3,Vec(0,0,1),144.0);
	rot_pose(t4,Vec(0,0,1),216.0);
	rot_pose(t5,Vec(0,0,1),288.0);
	for(Size i = 1; i <= t2.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
	for(Size i = 1; i <= t3.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t3.residue(i),1); else pose.append_residue_by_bond(t3.residue(i));
	for(Size i = 1; i <= t4.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t4.residue(i),1); else pose.append_residue_by_bond(t4.residue(i));
	for(Size i = 1; i <= t5.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t5.residue(i),1); else pose.append_residue_by_bond(t5.residue(i));
}

void run( Size itrifile, Size ipntfile ) {
  using basic::options::option;
  using namespace basic::options::OptionKeys;
  using namespace core::id;
  using numeric::conversions::radians;

  core::chemical::ResidueTypeSetCAP crs=core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID);

  Pose t_in,p_in;
  core::import_pose::pose_from_pdb(t_in,*crs,option[in::file::s]()[itrifile]);
  core::import_pose::pose_from_pdb(p_in,*crs,option[in::file::s]()[ipntfile]);
	std::string trifile = utility::file_basename(option[in::file::s]()[itrifile]);
	//TR << "make trimer and pentamer" << std::endl;
	make_trimer(t_in);
	//t_in.dump_pdb("trimer.pdb");
	make_pentamer(p_in);
	//p_in.dump_pdb("pentamer.pdb");
	//std::exit(-1);

  // set up geometry
  Vecf taxs = Vec( 0.000000, 0.000000,1.000000).normalized();
  Vecf tax2 = Vec(-0.333333,-0.577350,0.745356).normalized(); // 33.4458470159, 10.42594
  Vecf paxs = Vec(-0.607226, 0.000000,0.794529).normalized();
  Vecf pax2 = Vec(-0.491123,-0.850651,0.187593).normalized(); // 63.4311873349, 5.706642
  Real alpha = angle_degrees(taxs,Vec(0,0,0),paxs);
  rot_pose(p_in,Vec(0,1,0),-alpha,Vec(0,0,0));

  // Vecf tcom(0,0,0),pcom(0,0,0);
  // for(Size i = 1; i <= t_in.n_residue()/3; ++i) for(Size j = 1; j <= 5; ++j) tcom += t_in.residue(i).xyz(j);
  // tcom /= double(5*t_in.n_residue()/3);
  // for(Size i = 1; i <= p_in.n_residue()/5; ++i) for(Size j = 1; j <= 5; ++j) pcom += p_in.residue(i).xyz(j);
  // pcom /= double(5*t_in.n_residue()/5);
  // rot_pose(t_in,taxs,dihedral_degrees(paxs,Vec(0,0,0),taxs,tcom));
  // rot_pose(p_in,paxs,dihedral_degrees(taxs,Vec(0,0,0),paxs,pcom));
  // t_in.dump_pdb("t0.pdb");
  // p_in.dump_pdb("p0.pdb");
  // utility_exit_with_message("debug");

  Pose const tinit(t_in);
  Pose const pinit(p_in);


  int ANGLE_INCR = 3;
  if(ANGLE_INCR != 3) utility_exit_with_message("first ANGLE_INCR must be 3!!!");
  ObjexxFCL::FArray3D<int>   cbcount3((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0);
  {
    // compute high/low min dis for pent and tri here, input to sicfast and don't allow any below
    TR << "precomputing pent-pent and tri-tri interactions every 3°" << std::endl;
    vector1<double>            pntmnpos(72/ANGLE_INCR,0),   trimnpos(120/ANGLE_INCR,0);
    vector1<double>            pntmnneg(72/ANGLE_INCR,0),   trimnneg(120/ANGLE_INCR,0);
    ObjexxFCL::FArray2D<int>   pntcbpos(72/ANGLE_INCR,97,0),tricbpos(120/ANGLE_INCR,145,0);
    ObjexxFCL::FArray2D<int>   pntcbneg(72/ANGLE_INCR,97,0),tricbneg(120/ANGLE_INCR,145,0);
    {
      Pose p = pinit;
      for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
        rot_pose(p,paxs,ANGLE_INCR);
        vector1<Vecf> ppnt,cbp; // precompute these
        for(int i = 1; i <= (int)p.n_residue(); ++i) {
          cbp.push_back(Vecf(p.residue(i).xyz(2)));
          int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
          for(int j = 1; j <= natom; ++j) ppnt.push_back(Vecf(p.residue(i).xyz(j)));
        }
        vector1<Vecf> ppn2(ppnt),cb2(cbp);
        Matf r = rotation_matrix_degrees( paxs.cross(pax2), angle_degrees(paxs,Vec(0,0,0),pax2) );
        for(vector1<Vecf>::iterator i = ppn2.begin(); i != ppn2.end(); ++i) *i = r*(*i);
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
        int cbcount = 0;
        double const d = sicfast(ppnt,ppn2,cbp,cb2,(pax2-paxs).normalized(),cbcount,false);
        if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntpos! "+ObjexxFCL::string_of(ipnt));
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pax2-paxs).normalized();
        pntmnpos[ipnt/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(pax2,Vec(0,0,0),paxs)/2.0 );
        pntcbpos(ipnt/ANGLE_INCR+1,1) = cbcount;
        // std::cerr << "compute CB" << std::endl;
        set_cb_pairs(cbp,cb2);
        for(int i = 2; i <= 97; ++i) {
          for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) + 0.1*paxs;
          for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1*pax2;
          int cbc = 0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100) cbc++;
          pntcbpos(ipnt/ANGLE_INCR+1,i) = cbc;
          if(cbc==0) break;
        }
      }
      for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
        // Pose p = pinit;
        rot_pose(p,paxs,ANGLE_INCR);
        vector1<Vecf> ppnt,cbp; // precompute these
        for(int i = 1; i <= (int)p.n_residue(); ++i) {
          cbp.push_back(Vecf(p.residue(i).xyz(2)));
          int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
          for(int j = 1; j <= natom; ++j) ppnt.push_back(Vecf(p.residue(i).xyz(j)));
        }
        vector1<Vecf> ppn2(ppnt),cb2(cbp);
        Matf r = rotation_matrix_degrees( paxs.cross(pax2), angle_degrees(paxs,Vec(0,0,0),pax2) );
        for(vector1<Vecf>::iterator i = ppn2.begin(); i != ppn2.end(); ++i) *i = r*(*i);
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
        int cbcount = 0;
        double const d = sicfast(ppnt,ppn2,cbp,cb2,(paxs-pax2).normalized(),cbcount,false);
        if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntneg! "+ObjexxFCL::string_of(ipnt));
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(paxs-pax2).normalized();
        pntmnneg[ipnt/ANGLE_INCR+1] = d/2.0/sin( angle_radians(pax2,Vec(0,0,0),paxs)/2.0 );
        pntcbneg(ipnt/ANGLE_INCR+1,1) = cbcount;
        // std::cerr << "compute CB" << std::endl;
        set_cb_pairs(cbp,cb2);
        for(int i = 2; i <= 97; ++i) {
          for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) - 0.1*paxs;
          for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1*pax2;
          int cbc = 0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100) cbc++;
          pntcbneg(ipnt/ANGLE_INCR+1,i) = cbc;
          if(cbc==0) break;
        }
      }


      for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
        Pose t = tinit;
        rot_pose(t,taxs,(Real)itri);
        vector1<Vecf> ptri,cbt; // precompute these
        for(int i = 1; i <= (int)t.n_residue(); ++i) {
          cbt.push_back(Vecf(t.residue(i).xyz(2)));
          int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
          for(int j = 1; j <= natom; ++j) ptri.push_back(Vecf(t.residue(i).xyz(j)));
        }
        vector1<Vecf> ptr2(ptri),cb2(cbt);
        Matf r = rotation_matrix_degrees( taxs.cross(tax2), angle_degrees(taxs,Vec(0,0,0),tax2) );
        for(vector1<Vecf>::iterator i = ptr2.begin(); i != ptr2.end(); ++i) *i = r*(*i);
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
        int cbcount = 0;
        double const d = sicfast(ptri,ptr2,cbt,cb2,(tax2-taxs).normalized(),cbcount,false);
        if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for tripos! "+ObjexxFCL::string_of(itri));
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(tax2-taxs).normalized();
        trimnpos[itri/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 );
        tricbpos(itri/ANGLE_INCR+1,1) = cbcount;
        // std::cerr << "compute CB" << std::endl;
        set_cb_pairs(cbt,cb2);
        for(int i = 2; i <= 145; ++i) {
          for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) + 0.1*taxs;
          for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1*tax2;
          int cbc = 0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100) cbc++;
          tricbpos(itri/ANGLE_INCR+1,i) = cbc;
          if(cbc==0) break;
        }
        // if(itri != 60) continue;
        // Pose p1(t),p2(t);
        // rot_pose(p2,r);
        // trans_pose(p1,taxs*trimnpos[itri/ANGLE_INCR+1]);
        // trans_pose(p2,tax2*trimnpos[itri/ANGLE_INCR+1]);
        // p1.dump_pdb("tritri1.pdb");
        // p2.dump_pdb("tritri2.pdb");
        // TR << "PNT NEG D " << -d/2.0/sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 ) << std::endl;
        // for(int i = 1; i <= 145; ++i) {
        //  TR << "CB " << i << " " << pose_cbcount(p1,p2) << " " << tricbpos(itri/ANGLE_INCR+1,i) << std::endl;
        //  trans_pose(p1,taxs*0.1);
        //  trans_pose(p2,tax2*0.1);
        // }
        // TR << "D " << trimnpos[itri/ANGLE_INCR+1] << std::endl;
        // // utility_exit_with_message("test");
      }
      for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
        Pose t = tinit;
        rot_pose(t,taxs,(Real)itri);
        vector1<Vecf> ptri,cbt; // precompute these
        for(int i = 1; i <= (int)t.n_residue(); ++i) {
          cbt.push_back(Vecf(t.residue(i).xyz(2)));
          int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
          for(int j = 1; j <= natom; ++j) ptri.push_back(Vecf(t.residue(i).xyz(j)));
        }
        vector1<Vecf> ptr2(ptri),cb2(cbt);
        Matf r = rotation_matrix_degrees( taxs.cross(tax2), angle_degrees(taxs,Vec(0,0,0),tax2) );
        for(vector1<Vecf>::iterator i = ptr2.begin(); i != ptr2.end(); ++i) *i = r*(*i);
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
        int cbcount = 0;
        double const d = sicfast(ptri,ptr2,cbt,cb2,(taxs-tax2).normalized(),cbcount,false);
        if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for trineg! "+ObjexxFCL::string_of(itri));
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(taxs-tax2).normalized();
        trimnneg[itri/ANGLE_INCR+1] = d/2.0/sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 );
        tricbneg(itri/ANGLE_INCR+1,1) = cbcount;
        // std::cerr << "compute CB" << std::endl;
        set_cb_pairs(cbt,cb2);
        for(int i = 2; i <= 145; ++i) {
          for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) - 0.1*taxs;
          for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1*tax2;
          int cbc = 0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100) cbc++;
          tricbneg(itri/ANGLE_INCR+1,i) = cbc;
          if(cbc==0) break;
        }
        // if(itri < 70) continue;
        // Pose p1(t),p2(t);
        // rot_pose(p2,r);
        // trans_pose(p1,taxs*trimnneg[itri/ANGLE_INCR+1]);
        // trans_pose(p2,tax2*trimnneg[itri/ANGLE_INCR+1]);
        // p1.dump_pdb("test1.pdb");
        // p2.dump_pdb("test2.pdb");
        // TR << "PNT NEG D " << d << " " << sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 ) << std::endl;
        // for(int i = 1; i <= 145; ++i) {
        //  TR << "CB " << i << " " << pose_cbcount(p1,p2) << " " << tricbneg(itri/ANGLE_INCR+1,i) << std::endl;
        //  trans_pose(p1,-taxs*0.1);
        //  trans_pose(p2,-tax2*0.1);
        // }
        // TR << "D " << trimnneg[itri/ANGLE_INCR+1] << std::endl;
        // utility_exit_with_message("test");
      }
    }

    TR << "looping over ipnt, itri, iori every 3 degrees" << std::endl;
    ObjexxFCL::FArray3D<int>   cbcount((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0);
    ObjexxFCL::FArray3D<float> surfdis3((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0.0);
    {
      double mxd = 0;
      int cbmax = 0, mxiori = 0, mxipnt = 0, mxitri = 0;
      for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
        if(ipnt%3==0 && ipnt!=0) TR << "   loop1 " << trifile << " " << (100*ipnt)/72 << "%% done, cbmax:" << cbmax << std::endl;
        Pose p = pinit;
        rot_pose(p,paxs,(Real)ipnt);
        vector1<Vecf> pb,cbb; // precompute these
        for(int i = 1; i <= (int)p.n_residue(); ++i) {
          cbb.push_back(Vecf(p.residue(i).xyz(2)));
          int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
          for(int j = 1; j <= natom; ++j) pb.push_back(Vecf(p.residue(i).xyz(j)));
        }
        for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
          Pose t = tinit;
          rot_pose(t,taxs,(Real)itri);
          vector1<Vecf> pa,cba; // precompute these
          for(int i = 1; i <= (int)t.n_residue(); ++i) {
            cba.push_back(Vecf(t.residue(i).xyz(2)));
            int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
            for(int j = 1; j <= natom; ++j) pa.push_back(Vecf(t.residue(i).xyz(j)));
          }

          int iori = -1, ori_stage = 1;
          bool newstage = true;
          // for(iori = 0; iori < 360; iori+=ANGLE_INCR)
          while(ori_stage < 5) {
            if(newstage) {
              if( ori_stage == 1 || ori_stage == 2 ) iori = (int)( 90.0+double(ANGLE_INCR)/2.0+angle_degrees(taxs,Vecf(0,0,0),paxs));
              if( ori_stage == 3 || ori_stage == 4 ) iori = (int)(270.0+double(ANGLE_INCR)/2.0+angle_degrees(taxs,Vecf(0,0,0),paxs));
              iori = (iori / ANGLE_INCR) * ANGLE_INCR; // round to closest multiple of angle incr
              if( ori_stage == 2 || ori_stage == 4 ) iori -= ANGLE_INCR;
              newstage = false;
            } else {
              if( ori_stage == 1 || ori_stage == 3 ) iori += ANGLE_INCR;
              if( ori_stage == 2 || ori_stage == 4 ) iori -= ANGLE_INCR;
            }
            // TR << "IORI " << iori << std::endl;
            Vecf sicaxis = (rotation_matrix_degrees(-Vec(0,1,0),(Real)iori) * Vec(0,0,1)).normalized();
            int tmpcbc;
            double d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);
            double theta = iori;
            double gamma = theta-alpha;
            double x = d * sin(numeric::conversions::radians(gamma));
            double y = d * cos(numeric::conversions::radians(gamma));
            double w = x / sin(numeric::conversions::radians(alpha));
            double z = x / tan(numeric::conversions::radians(alpha));
            double dpnt = y+z;
            double dtri = w;
            double pntmn,trimn;
            if( w > 0 ) {
              pntmn = pntmnpos[ipnt/ANGLE_INCR+1];
              trimn = trimnpos[itri/ANGLE_INCR+1];
              int dp = (int)(dpnt-pntmn)*10+1;
              int dt = (int)(dtri-trimn)*10+1;
              if( dp < 1 ) { ori_stage++; newstage=true; continue; };
              if( dt < 1 ) { ori_stage++; newstage=true; continue; };
              // if(ipnt==18 && itri==72 && iori==276)
              //  TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbpos(ipnt/ANGLE_INCR+1,dp)
              //    << "    " << dt << " " << dtri << " " << trimn << " " << tricbpos(itri/ANGLE_INCR+1,dt) << std::endl;
              if( dp <= 97  ) tmpcbc += pntcbpos(ipnt/ANGLE_INCR+1,dp);
              if( dt <= 145 ) tmpcbc += tricbpos(itri/ANGLE_INCR+1,dt);
              // TR << "CHK " << dpnt << " " << pntmn << "    " << dtri << " " << trimn << std::endl;
            } else {
              pntmn = pntmnneg[ipnt/ANGLE_INCR+1];
              trimn = trimnneg[itri/ANGLE_INCR+1];
              int dp = (int)(-dpnt+pntmn)*10+1;
              int dt = (int)(-dtri+trimn)*10+1;
              if( dp < 1 ) { ori_stage++; newstage=true; continue; };
              if( dt < 1 ) { ori_stage++; newstage=true; continue; };
              // if(ipnt==18 && itri==72 && iori==276)
              //  TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbneg(ipnt/ANGLE_INCR+1,dp)
              //    << "    " << dt << " " << dtri << " " << trimn << " " << pntcbneg(itri/ANGLE_INCR+1,dt) << std::endl;
              if( dp <= 97  ) tmpcbc += pntcbneg(ipnt/ANGLE_INCR+1,dp);
              if( dt <= 145 ) tmpcbc += tricbneg(itri/ANGLE_INCR+1,dt);
            }

            surfdis3(ipnt/ANGLE_INCR+1,itri/ANGLE_INCR+1,iori/ANGLE_INCR+1) = d;
            cbcount(ipnt/ANGLE_INCR+1,itri/ANGLE_INCR+1,iori/ANGLE_INCR+1) = tmpcbc;
            // d = sicfast(t,p,sicaxis,cbcount);
            // TR << "trial " << ipnt << " " << itri << " " << iori << " " << d << " " << cbcount << std::endl;
            if(tmpcbc > cbmax) {
              cbmax = tmpcbc;
              mxiori = iori;
              mxipnt = ipnt;
              mxitri = itri;
              mxd = d;
            }
          }
          // utility_exit_with_message("testing");
        }
      }
      TR << "MAX3 " << mxipnt << " " << mxitri << " " << mxiori << " " << cbmax << " " << mxd << std::endl;
      // TR << "MAX " << mxipnt << " " << mxitri << " " << mxiori << " " << cbcount(mxipnt/ANGLE_INCR+1,mxitri/ANGLE_INCR+1,mxiori/ANGLE_INCR+1) << " " << surfdis3(mxipnt/ANGLE_INCR+1,mxitri/ANGLE_INCR+1,mxiori/ANGLE_INCR+1) << std::endl;
    }

    cbcount3 = cbcount;
    // int delta = ceil(3.0/(Real)ANGLE_INCR);
    // TR << "scanning results for local maxima +- " << delta*ANGLE_INCR << "°" << std::endl;
    // for(int ipnt = 1; ipnt <=  72/ANGLE_INCR; ++ipnt) {
    //  for(int itri = 1; itri <= 120/ANGLE_INCR; ++itri) {
    //    for(int iori = 1; iori <= 360/ANGLE_INCR; ++iori) {
    //      // std::cerr << ipnt << " " << itri << " " << iori << " " << cbcount(ipnt,itri,iori) << std::endl;
    //      int lmaxcb = 0;
    //      for(int dipnt = -delta; dipnt <= delta; ++dipnt) {
    //        for(int ditri = -delta; ditri <= delta; ++ditri) {
    //          for(int diori = -delta; diori <= delta; ++diori) {
    //            int i = (ipnt+dipnt-1+( 72/ANGLE_INCR)) % ( 72/ANGLE_INCR) + 1;
    //            int j = (itri+ditri-1+(120/ANGLE_INCR)) % (120/ANGLE_INCR) + 1;
    //            int k = (iori+diori-1+(360/ANGLE_INCR)) % (360/ANGLE_INCR) + 1;
    //            if( cbcount(i,j,k) > lmaxcb ) {
    //              lmaxcb = cbcount(i,j,k);
    //            }
    //          }
    //        }
    //      }
    //      if( cbcount(ipnt,itri,iori) >= lmaxcb ) cbcount3(ipnt,itri,iori) = cbcount(ipnt,itri,iori);
    //      if( ipnt==2 && itri==15 && iori==99 ) TR << "NEARMAX: " << cbcount(ipnt,itri,iori) << " " << lmaxcb << std::endl;
    //    }
    //  }
    // }
  }

  vector1<int> hitpnt,hittri,hitori,hitcbc;
  ANGLE_INCR = 1;
  if(ANGLE_INCR != 1) utility_exit_with_message("second ANGLE_INCR must be 1!!!");
  ObjexxFCL::FArray3D<float> surfdis1((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0.0);
  ObjexxFCL::FArray3D<int>   cbcount1((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0);
  {
    // compute high/low min dis for pent and tri here, input to sicfast and don't allow any below
    TR << "precomputing pent-pent and tri-tri interactions every 1°" << std::endl;
    vector1<double>            pntmnpos(72/ANGLE_INCR,0),   trimnpos(120/ANGLE_INCR,0);
    vector1<double>            pntmnneg(72/ANGLE_INCR,0),   trimnneg(120/ANGLE_INCR,0);
    ObjexxFCL::FArray2D<int>   pntcbpos(72/ANGLE_INCR,97,0),tricbpos(120/ANGLE_INCR,145,0);
    ObjexxFCL::FArray2D<int>   pntcbneg(72/ANGLE_INCR,97,0),tricbneg(120/ANGLE_INCR,145,0);
    {
      Pose p = pinit;
      for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
        rot_pose(p,paxs,ANGLE_INCR);
        vector1<Vecf> ppnt,cbp; // precompute these
        for(int i = 1; i <= (int)p.n_residue(); ++i) {
          cbp.push_back(Vecf(p.residue(i).xyz(2)));
          int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
          for(int j = 1; j <= natom; ++j) ppnt.push_back(Vecf(p.residue(i).xyz(j)));
        }
        vector1<Vecf> ppn2(ppnt),cb2(cbp);
        Matf r = rotation_matrix_degrees( paxs.cross(pax2), angle_degrees(paxs,Vec(0,0,0),pax2) );
        for(vector1<Vecf>::iterator i = ppn2.begin(); i != ppn2.end(); ++i) *i = r*(*i);
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
        int cbcount = 0;
        double const d = sicfast(ppnt,ppn2,cbp,cb2,(pax2-paxs).normalized(),cbcount,false);
        if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntpos! "+ObjexxFCL::string_of(ipnt));
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pax2-paxs).normalized();
        pntmnpos[ipnt/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(pax2,Vec(0,0,0),paxs)/2.0 );
        pntcbpos(ipnt/ANGLE_INCR+1,1) = cbcount;
        // std::cerr << "compute CB" << std::endl;
        set_cb_pairs(cbp,cb2);
        for(int i = 2; i <= 97; ++i) {
          for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) + 0.1*paxs;
          for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1*pax2;
          int cbc = 0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100) cbc++;
          pntcbpos(ipnt/ANGLE_INCR+1,i) = cbc;
          if(cbc==0) break;
        }
      }
      for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
        // Pose p = pinit;
        rot_pose(p,paxs,ANGLE_INCR);
        vector1<Vecf> ppnt,cbp; // precompute these
        for(int i = 1; i <= (int)p.n_residue(); ++i) {
          cbp.push_back(Vecf(p.residue(i).xyz(2)));
          int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
          for(int j = 1; j <= natom; ++j) ppnt.push_back(Vecf(p.residue(i).xyz(j)));
        }
        vector1<Vecf> ppn2(ppnt),cb2(cbp);
        Matf r = rotation_matrix_degrees( paxs.cross(pax2), angle_degrees(paxs,Vec(0,0,0),pax2) );
        for(vector1<Vecf>::iterator i = ppn2.begin(); i != ppn2.end(); ++i) *i = r*(*i);
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
        int cbcount = 0;
        double const d = sicfast(ppnt,ppn2,cbp,cb2,(paxs-pax2).normalized(),cbcount,false);
        if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntneg! "+ObjexxFCL::string_of(ipnt));
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(paxs-pax2).normalized();
        pntmnneg[ipnt/ANGLE_INCR+1] = d/2.0/sin( angle_radians(pax2,Vec(0,0,0),paxs)/2.0 );
        pntcbneg(ipnt/ANGLE_INCR+1,1) = cbcount;
        // std::cerr << "compute CB" << std::endl;
        set_cb_pairs(cbp,cb2);
        for(int i = 2; i <= 97; ++i) {
          for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) - 0.1*paxs;
          for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1*pax2;
          int cbc = 0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100) cbc++;
          pntcbneg(ipnt/ANGLE_INCR+1,i) = cbc;
          if(cbc==0) break;
        }
      }


      for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
        Pose t = tinit;
        rot_pose(t,taxs,(Real)itri);
        vector1<Vecf> ptri,cbt; // precompute these
        for(int i = 1; i <= (int)t.n_residue(); ++i) {
          cbt.push_back(Vecf(t.residue(i).xyz(2)));
          int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
          for(int j = 1; j <= natom; ++j) ptri.push_back(Vecf(t.residue(i).xyz(j)));
        }
        vector1<Vecf> ptr2(ptri),cb2(cbt);
        Matf r = rotation_matrix_degrees( taxs.cross(tax2), angle_degrees(taxs,Vec(0,0,0),tax2) );
        for(vector1<Vecf>::iterator i = ptr2.begin(); i != ptr2.end(); ++i) *i = r*(*i);
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
        int cbcount = 0;
        double const d = sicfast(ptri,ptr2,cbt,cb2,(tax2-taxs).normalized(),cbcount,false);
        if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for tripos! "+ObjexxFCL::string_of(itri));
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(tax2-taxs).normalized();
        trimnpos[itri/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 );
        tricbpos(itri/ANGLE_INCR+1,1) = cbcount;
        // std::cerr << "compute CB" << std::endl;
        set_cb_pairs(cbt,cb2);
        for(int i = 2; i <= 145; ++i) {
          for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) + 0.1*taxs;
          for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1*tax2;
          int cbc = 0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100) cbc++;
          tricbpos(itri/ANGLE_INCR+1,i) = cbc;
          if(cbc==0) break;
        }
        // if(itri != 60) continue;
        // Pose p1(t),p2(t);
        // rot_pose(p2,r);
        // trans_pose(p1,taxs*trimnpos[itri/ANGLE_INCR+1]);
        // trans_pose(p2,tax2*trimnpos[itri/ANGLE_INCR+1]);
        // p1.dump_pdb("tritri1.pdb");
        // p2.dump_pdb("tritri2.pdb");
        // TR << "PNT NEG D " << -d/2.0/sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 ) << std::endl;
        // for(int i = 1; i <= 145; ++i) {
        //  TR << "CB " << i << " " << pose_cbcount(p1,p2) << " " << tricbpos(itri/ANGLE_INCR+1,i) << std::endl;
        //  trans_pose(p1,taxs*0.1);
        //  trans_pose(p2,tax2*0.1);
        // }
        // TR << "D " << trimnpos[itri/ANGLE_INCR+1] << std::endl;
        // // utility_exit_with_message("test");
      }
      for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
        Pose t = tinit;
        rot_pose(t,taxs,(Real)itri);
        vector1<Vecf> ptri,cbt; // precompute these
        for(int i = 1; i <= (int)t.n_residue(); ++i) {
          cbt.push_back(Vecf(t.residue(i).xyz(2)));
          int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
          for(int j = 1; j <= natom; ++j) ptri.push_back(Vecf(t.residue(i).xyz(j)));
        }
        vector1<Vecf> ptr2(ptri),cb2(cbt);
        Matf r = rotation_matrix_degrees( taxs.cross(tax2), angle_degrees(taxs,Vec(0,0,0),tax2) );
        for(vector1<Vecf>::iterator i = ptr2.begin(); i != ptr2.end(); ++i) *i = r*(*i);
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = r*(*i);
        int cbcount = 0;
        double const d = sicfast(ptri,ptr2,cbt,cb2,(taxs-tax2).normalized(),cbcount,false);
        if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for trineg! "+ObjexxFCL::string_of(itri));
        for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(taxs-tax2).normalized();
        trimnneg[itri/ANGLE_INCR+1] = d/2.0/sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 );
        tricbneg(itri/ANGLE_INCR+1,1) = cbcount;
        // std::cerr << "compute CB" << std::endl;
        set_cb_pairs(cbt,cb2);
        for(int i = 2; i <= 145; ++i) {
          for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) - 0.1*taxs;
          for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1*tax2;
          int cbc = 0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100) cbc++;
          tricbneg(itri/ANGLE_INCR+1,i) = cbc;
          if(cbc==0) break;
        }
        // if(itri < 70) continue;
        // Pose p1(t),p2(t);
        // rot_pose(p2,r);
        // trans_pose(p1,taxs*trimnneg[itri/ANGLE_INCR+1]);
        // trans_pose(p2,tax2*trimnneg[itri/ANGLE_INCR+1]);
        // p1.dump_pdb("test1.pdb");
        // p2.dump_pdb("test2.pdb");
        // TR << "PNT NEG D " << d << " " << sin( angle_radians(tax2,Vec(0,0,0),taxs)/2.0 ) << std::endl;
        // for(int i = 1; i <= 145; ++i) {
        //  TR << "CB " << i << " " << pose_cbcount(p1,p2) << " " << tricbneg(itri/ANGLE_INCR+1,i) << std::endl;
        //  trans_pose(p1,-taxs*0.1);
        //  trans_pose(p2,-tax2*0.1);
        // }
        // TR << "D " << trimnneg[itri/ANGLE_INCR+1] << std::endl;
        // utility_exit_with_message("test");
      }
    }

    int top100_3 = 0;
    {
      vector1<int> cbtmp;
      for(Size i = 0; i < cbcount3.size(); ++i) {
        if(cbcount3[i] > 0) cbtmp.push_back(cbcount3[i]);
      }
      std::sort(cbtmp.begin(),cbtmp.end());
      top100_3 = cbtmp[max(1,(int)cbtmp.size()-999)];
      TR << "scanning top100 with cbcount3 >= " << top100_3 << std::endl;
    }
    // assuming ANGLE_INCR = 1!!!!!!!!!!!
    utility::vector1<vector1<int> > pntlmx,trilmx,orilmx;
    for(int ipnt = 0; ipnt < 72; ipnt+=3) {
      for(int itri = 0; itri < 120; itri+=3) {
        for(int iori = 0; iori < 360; iori+=3) {
          if( cbcount3(ipnt/3+1,itri/3+1,iori/3+1) >= top100_3) {
            // TR << "checking around " << ipnt << " " << itri << " " << iori << " " << cbcount3(ipnt/3+1,itri/3+1,iori/3+1) << std::endl;
            vector1<int> pnt,tri,ori;
            for(int i = -1; i <= 1; ++i) pnt.push_back( (ipnt+i+ 72)% 72 );
            for(int j = -1; j <= 1; ++j) tri.push_back( (itri+j+120)%120 );
            for(int k = -1; k <= 1; ++k) ori.push_back( (iori+k+360)%360 );
            pntlmx.push_back(pnt);
            trilmx.push_back(tri);
            orilmx.push_back(ori);
          }
        }
      }
    }

    TR << "looping over top " << pntlmx.size() << " 3degree hits, +-1 degree " << std::endl;
    {
      int max1 = 0;
      for(Size ilmx = 1; ilmx <= pntlmx.size(); ++ilmx) {
        if( (ilmx-1)%100==0 && ilmx!=1) TR << "  loop2 " << trifile << " " << (double(ilmx-1)/10) << "%% done, cbmax: " << max1 << std::endl;
        // TR << "checking around local max # " << ilmx << std::endl;
        int cbmax = 0, mxiori = 0, mxipnt = 0, mxitri = 0;
        for(vector1<int>::const_iterator pipnt = pntlmx[ilmx].begin(); pipnt != pntlmx[ilmx].end(); ++pipnt) {
          int ipnt = *pipnt;
          Pose p = pinit;
          rot_pose(p,paxs,(Real)ipnt);
          vector1<Vecf> pb,cbb; // precompute these
          for(int i = 1; i <= (int)p.n_residue(); ++i) {
            cbb.push_back(Vecf(p.residue(i).xyz(2)));
            int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
            for(int j = 1; j <= natom; ++j) pb.push_back(Vecf(p.residue(i).xyz(j)));
          }
          for(vector1<int>::const_iterator pitri = trilmx[ilmx].begin(); pitri != trilmx[ilmx].end(); ++pitri) {
            int itri = *pitri;
            Pose t = tinit;
            rot_pose(t,taxs,(Real)itri);
            vector1<Vecf> pa,cba; // precompute these
            for(int i = 1; i <= (int)t.n_residue(); ++i) {
              cba.push_back(Vecf(t.residue(i).xyz(2)));
              int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
              for(int j = 1; j <= natom; ++j) pa.push_back(Vecf(t.residue(i).xyz(j)));
            }
            for(vector1<int>::const_iterator piori = orilmx[ilmx].begin(); piori != orilmx[ilmx].end(); ++piori) {
              int iori = *piori;
              Vecf sicaxis = (rotation_matrix_degrees(-Vec(0,1,0),(Real)iori) * Vec(0,0,1)).normalized();
              int tmpcbc;
              double d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);
              double theta = iori;
              double gamma = theta-alpha;
              double x = d * sin(numeric::conversions::radians(gamma));
              double y = d * cos(numeric::conversions::radians(gamma));
              double w = x / sin(numeric::conversions::radians(alpha));
              double z = x / tan(numeric::conversions::radians(alpha));
              double dpnt = y+z;
              double dtri = w;
              double pntmn,trimn;
              if( w > 0 ) {
                pntmn = pntmnpos[ipnt/ANGLE_INCR+1];
                trimn = trimnpos[itri/ANGLE_INCR+1];
                int dp = (int)(dpnt-pntmn)*10+1;
                int dt = (int)(dtri-trimn)*10+1;
                if( dp < 1 ) { continue; };
                if( dt < 1 ) { continue; };
                // if(ipnt==18 && itri==72 && iori==276)
                //  TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbpos(ipnt/ANGLE_INCR+1,dp)
                //    << "    " << dt << " " << dtri << " " << trimn << " " << tricbpos(itri/ANGLE_INCR+1,dt) << std::endl;
                if( dp <= 97  ) tmpcbc += pntcbpos(ipnt/ANGLE_INCR+1,dp);
                if( dt <= 145 ) tmpcbc += tricbpos(itri/ANGLE_INCR+1,dt);
                // TR << "CHK " << dpnt << " " << pntmn << "    " << dtri << " " << trimn << std::endl;
              } else {
                pntmn = pntmnneg[ipnt/ANGLE_INCR+1];
                trimn = trimnneg[itri/ANGLE_INCR+1];
                int dp = (int)(-dpnt+pntmn)*10+1;
                int dt = (int)(-dtri+trimn)*10+1;
                if( dp < 1 ) { continue; };
                if( dt < 1 ) { continue; };
                // if(ipnt==18 && itri==72 && iori==276)
                //  TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbneg(ipnt/ANGLE_INCR+1,dp)
                //    << "    " << dt << " " << dtri << " " << trimn << " " << pntcbneg(itri/ANGLE_INCR+1,dt) << std::endl;
                if( dp <= 97  ) tmpcbc += pntcbneg(ipnt/ANGLE_INCR+1,dp);
                if( dt <= 145 ) tmpcbc += tricbneg(itri/ANGLE_INCR+1,dt);
              }
              if(tmpcbc > cbmax) {
                cbmax = tmpcbc;
                mxiori = iori;
                mxipnt = ipnt;
                mxitri = itri;
              }
            }
          }
        }
        hitpnt.push_back(mxipnt);
        hittri.push_back(mxitri);
        hitori.push_back(mxiori);
        hitcbc.push_back(cbmax);
        if(cbmax > max1) max1 = cbmax;
        // TR << "HIT " << ilmx << " " << mxipnt << " " << mxitri << " " << mxiori << " " << cbmax << std::endl;
      }
    }

    int top10,max1;
    {
      vector1<int> hittmp = hitcbc;
      std::sort(hittmp.begin(),hittmp.end());
      max1  = hittmp[hittmp.size()];
      top10 = hittmp[hittmp.size()-9];
      TR << "top10 " << top10 << std::endl;
    }
    for(Size ihit = 1; ihit <= hitcbc.size(); ++ihit) {
      if(hitcbc[ihit]==max1)  TR << "MAX1 " << hitpnt[ihit] << " " << hittri[ihit] << " " << hitori[ihit] << " " << hitcbc[ihit] << std::endl;
    }

    for(Size ihit = 1; ihit <= hitcbc.size(); ++ihit) {
      if(hitcbc[ihit] >= top10 ) {
        int ipnt = hitpnt[ihit];
        int itri = hittri[ihit];
        int iori = hitori[ihit];

        Pose p = pinit;
        Pose t = tinit;
        rot_pose(p,paxs,(Real)ipnt);
        rot_pose(t,taxs,(Real)itri);

        vector1<Vecf> pb,cbb; // precompute these
        for(int i = 1; i <= (int)p.n_residue(); ++i) {
          cbb.push_back(Vecf(p.residue(i).xyz(2)));
          int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
          for(int j = 1; j <= natom; ++j) pb.push_back(Vecf(p.residue(i).xyz(j)));
        }
        vector1<Vecf> pa,cba; // precompute these
        for(int i = 1; i <= (int)t.n_residue(); ++i) {
          cba.push_back(Vecf(t.residue(i).xyz(2)));
          int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
          for(int j = 1; j <= natom; ++j) pa.push_back(Vecf(t.residue(i).xyz(j)));
        }
        Vecf sicaxis = (rotation_matrix_degrees(-Vec(0,1,0),(Real)iori) * Vec(0,0,1)).normalized();
        int tmpcbc;
        double d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);

        double theta = iori;
        double gamma = theta-alpha;
        double x = d * sin(numeric::conversions::radians(gamma));
        double y = d * cos(numeric::conversions::radians(gamma));
        double w = x / sin(numeric::conversions::radians(alpha));
        double z = x / tan(numeric::conversions::radians(alpha));
        double dpnt = y+z;
        double dtri = w;
        double pntmn,trimn;
        trans_pose(p,dpnt*paxs);
        trans_pose(t,dtri*taxs);
        if( w > 0 ) {
          pntmn = pntmnpos[ipnt/ANGLE_INCR+1];
          trimn = trimnpos[itri/ANGLE_INCR+1];
          int dp = (int)(dpnt-pntmn)*10+1;
          int dt = (int)(dtri-trimn)*10+1;
          if( dp < 1 ) { continue; };
          if( dt < 1 ) { continue; };
          // if(ipnt==18 && itri==72 && iori==276)
          //  TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbpos(ipnt/ANGLE_INCR+1,dp)
          //    << "    " << dt << " " << dtri << " " << trimn << " " << tricbpos(itri/ANGLE_INCR+1,dt) << std::endl;
          if( dp <= 97  ) tmpcbc += pntcbpos(ipnt/ANGLE_INCR+1,dp);
          if( dt <= 145 ) tmpcbc += tricbpos(itri/ANGLE_INCR+1,dt);
          // TR << "CHK " << dpnt << " " << pntmn << "    " << dtri << " " << trimn << std::endl;
        } else {
          pntmn = pntmnneg[ipnt/ANGLE_INCR+1];
          trimn = trimnneg[itri/ANGLE_INCR+1];
          int dp = (int)(-dpnt+pntmn)*10+1;
          int dt = (int)(-dtri+trimn)*10+1;
          if( dp < 1 ) { continue; };
          if( dt < 1 ) { continue; };
          // if(ipnt==18 && itri==72 && iori==276)
          //  TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbneg(ipnt/ANGLE_INCR+1,dp)
          //    << "    " << dt << " " << dtri << " " << trimn << " " << pntcbneg(itri/ANGLE_INCR+1,dt) << std::endl;
          if( dp <= 97  ) tmpcbc += pntcbneg(ipnt/ANGLE_INCR+1,dp);
          if( dt <= 145 ) tmpcbc += tricbneg(itri/ANGLE_INCR+1,dt);
        }
				Real sizefac = pow(t_in.n_residue(),4.0/3.0) + pow(p_in.n_residue(),4.0/3.0);

        std::string fname;
				fname  = ObjexxFCL::lead_zero_string_of(Size(Real(tmpcbc)/sizefac*1000000),10);
				fname += "_" + utility::file_basename(option[in::file::s]()[1]);
        fname += "_" + utility::file_basename(option[in::file::s]()[2]);
        fname += "_" + ObjexxFCL::lead_zero_string_of(tmpcbc,4);
        fname += "_" + ObjexxFCL::lead_zero_string_of(ipnt,3);
        fname += "_" + ObjexxFCL::lead_zero_string_of(itri,3);
        fname += "_" + ObjexxFCL::lead_zero_string_of(iori,3);
        fname += ".pdb.gz";
        // continue;
        TR << "dumping file " << fname << std::endl;

        Pose symm;
        symm.append_residue_by_jump(t.residue(1),1);
        for(Size i = 2; i <= t.n_residue()/3; ++i) {
          if(symm.residue(i-1).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(symm,i-1);
          if(   t.residue(i  ).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(t,i);
          symm.append_residue_by_bond(t.residue(i));
        }
        symm.append_residue_by_jump(p.residue(1),1);
        for(Size i = 2; i <= p.n_residue()/5; ++i) symm.append_residue_by_bond(p.residue(i));

				TR << "making symm" << std::endl;
				core::util::switch_to_residue_type_set(symm,"fa_standard");
				core::chemical::ResidueTypeSetCAP rs = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
				for(Size i = 1; i <= symm.n_residue(); ++i) core::pose::replace_pose_residue_copying_existing_coordinates(symm,i,rs->name_map("ALA"));
        core::pose::symmetry::make_symmetric_pose(symm);
        core::io::pdb::dump_pdb(symm,option[out::file::o]()+"/"+fname);

        TR << ihit << " " << ipnt << " " << itri << " " << iori << " " << hitcbc[ihit] << " " << tmpcbc << " " << Real(tmpcbc)/sizefac << std::endl;

			}
    }
  }
}


int main (int argc, char *argv[]) {
  core::init(argc,argv);
	using namespace basic::options::OptionKeys;
	for(Size i = 2; i <= basic::options::option[in::file::s]().size(); ++i) {
		run(i,1);
	}
}




//
//








