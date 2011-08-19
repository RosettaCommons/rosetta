// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/rms/rms_util.tmpl.hh
/// @brief  RMS stuff from rosetta++ that requires templates
/// @author James Thompson
/// @author Ian Davis
/// @date   Wed Aug 22 12:10:37 2007

#ifndef core_scoring_rms_util_templ_HH
#define core_scoring_rms_util_templ_HH

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/rms_util.hh>

#include <numeric/model_quality/rms.hh>
// AUTO-REMOVED #include <numeric/model_quality/maxsub.hh>
#include <numeric/xyzVector.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <list>

//Auto Headers
#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
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
//XRW_B_T1
//#include <core/coarse/Translator.fwd.hh>
//XRW_E_T1
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
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
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <basic/MetricValue.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <utility/assert.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
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
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/types.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1A.fwd.hh>
#include <ObjexxFCL/FArray1P.fwd.hh>
#include <ObjexxFCL/FArray1P.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/StaticIndexRange.fwd.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
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
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>


namespace core {
namespace scoring {

template< class T >
core::Real
rmsd_with_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	std::list< core::Size > const & subset_residues,
	T* predicate
)
{
	ASSERT_ONLY(core::Size const nres1 = pose1.total_residue();)
	ASSERT_ONLY(core::Size const nres2 = pose2.total_residue();)
	assert( nres1 == nres2 );

	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for(std::list<Size>::const_iterator list_iter = subset_residues.begin();
	    list_iter != subset_residues.end(); list_iter++)
	  {
	    Size const i( *list_iter );
	    Size num_atoms ( pose1.residue(i).natoms() );
	    if ( predicate == is_ligand_heavyatom ||
		 predicate == is_polymer_heavyatom ||
		 predicate == is_heavyatom ) {
	      //			assert( pose1.residue(i).nheavyatoms() == pose2.residue(i).nheavyatoms() );
	    } else if ( num_atoms > pose2.residue(i).natoms() ){
	      num_atoms = pose2.residue(i).natoms();
	    }
	    for ( core::Size j = 1; j <= num_atoms; ++j ) {
	      if( j <= pose1.residue(i).natoms() &&
		  j <= pose2.residue(i).natoms() ) {
		if ( predicate( pose1, pose2, i, j ) ) {
		  p1_coords.push_back( pose1.residue(i).xyz(j) );
		  p2_coords.push_back( pose2.residue(i).xyz(j) );
		}
	      }
	    }
	  }
	assert( p1_coords.size() == p2_coords.size() );

	int const natoms = p1_coords.size();
	ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
	ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
	for ( int i = 0; i < natoms; ++i ) {
	  for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
	    p1a(k+1,i+1) = p1_coords[i][k];
	    p2a(k+1,i+1) = p2_coords[i][k];
	  }
	}

	return numeric::model_quality::rms_wrapper( natoms, p1a, p2a );
}

/// @brief Select atoms for RMS via a predicate function/functor.
/// @details Calculates minimal rms, allowing rotation/translation for best fit.
/// Parameter "predicate" should be a function pointer
/// [ or class that defines operator() ] with the following signature:
///
///   bool my_pred_func(Pose const & pose1, Pose const & pose2, Size resno, Size atomno);
///
/// It should return true if the atom should be included and false otherwise.
///
/// Example of use, to calculate C-alpha RMSD:
///   rmsd_with_super(pose1, pose2, is_protein_CA);
template< class T >
core::Real
rmsd_with_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	T* predicate
)
{
  std::list< core::Size > all_residues;
  for (core::Size n = 1; n <= pose1.total_residue(); n++ ){
    all_residues.push_back( n );
  }
  return rmsd_with_super( pose1, pose2, all_residues, predicate );
}



/// @brief Select a subset atoms for RMS via a predicate function/functor.
/// @details Calculates minimal rms, allowing rotation/translation for best fit.
///		Same as above function, but allows a subset of residues over which to
///		superposition to be passed in
/// Parameter "predicate" should be a function pointer
/// [ or class that defines operator() ] with the following signature:
///
///   bool my_pred_func(Pose const & pose1, Pose const & pose2, Size resno, Size atomno);
///
/// It should return true if the atom should be included and false otherwise.
///
/// Example of use, to calculate C-alpha RMSD:
///   rmsd_with_super(pose1, pose2, is_protein_CA, subset);
template< class T >
core::Real
rmsd_with_super_subset(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & subset,
	T* predicate
)
{
	core::Size const nres1 = pose1.total_residue();
	//core::Size const nres2 = pose2.total_residue();
	//assert( nres1 == nres2 );

	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for ( core::Size i = 1; i <= nres1; ++i ) {
		if ( subset( i ) ) {
			Size num_atoms ( pose1.residue(i).natoms() );
			if ( predicate == is_ligand_heavyatom ||
				predicate == is_polymer_heavyatom ||
				predicate == is_heavyatom ) {
				assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
			} else if ( num_atoms > pose2.residue(i).natoms() ){
				num_atoms = pose2.residue(i).natoms();
			}
			for ( core::Size j = 1; j <= num_atoms; ++j ) {
				if ( predicate( pose1, pose2, i, j ) ) {
					p1_coords.push_back( pose1.residue(i).xyz(j) );
					p2_coords.push_back( pose2.residue(i).xyz(j) );
				}
			}
		}
	}
	assert( p1_coords.size() == p2_coords.size() );

	int const natoms = p1_coords.size();
	ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
	ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
	for ( int i = 0; i < natoms; ++i ) {
		for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
			p1a(k+1,i+1) = p1_coords[i][k];
			p2a(k+1,i+1) = p2_coords[i][k];
		}
	}

	return numeric::model_quality::rms_wrapper( natoms, p1a, p2a );
}


/// @brief Select atoms for RMS via a predicate function/functor.
/// @details Calculates minimal rms, NOT allowing rotation/translation --
/// uses current coordinates as-is.
/// Parameter "predicate" should be a function pointer
/// [ or class that defines operator() ] with the following signature:
///
///   bool my_pred_func(Pose const & pose1, Pose const & pose2, Size resno, Size atomno);
///
/// It should return true if the atom should be included and false otherwise.
///
/// Example of use, to calculate C-alpha RMSD:
///   rmsd_no_super(pose1, pose2, is_protein_CA);
template< class T >
core::Real
rmsd_no_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	T* predicate
)
{
	core::Size const nres1 = pose1.total_residue();
	ASSERT_ONLY(core::Size const nres2 = pose2.total_residue();)
	assert( nres1 == nres2 );

	core::Real sum2( 0.0 );
	core::Size natoms( 0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {
		Size num_atoms ( pose1.residue(i).natoms() );
		if ( predicate == is_ligand_heavyatom ||
			predicate == is_polymer_heavyatom ||
			predicate == is_heavyatom ) {
			assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
		} else if ( num_atoms > pose2.residue(i).natoms() ){
			num_atoms = pose2.residue(i).natoms();
		}
		for ( core::Size j = 1; j <= num_atoms; ++j ) {
			if ( predicate( pose1, pose2, i, j ) ) {
				core::Vector diff = pose1.residue(i).xyz(j) - pose2.residue(i).xyz(j);
				sum2 += diff.length_squared();
				natoms += 1;
			}
		}
	}
	return std::sqrt(sum2 / natoms);
}

template< class T >
core::Real
rmsd_no_super(
	core::conformation::ResidueCOPs const & residues1,
	core::conformation::ResidueCOPs const & residues2,
	T* predicate
)
{
	core::Size const nres1 = residues1.size();
	ASSERT_ONLY(core::Size const nres2 = residues2.size();)
	assert( nres1 == nres2 );

	core::Real sum2( 0.0 );
	core::Size natoms( 0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {
		Size num_atoms ( residues1[i]->natoms() );
		if (  predicate == is_ligand_heavyatom_residues
			/* predicate == is_polymer_heavyatom || */
			/* predicate == is_heavyatom */) { // for Rotate.cc
			assert( residues1[i]->natoms() == residues2[i]->natoms() );
		} else if ( num_atoms > residues2[i]->natoms() ){
			num_atoms = residues2[i]->natoms();
		}
		for ( core::Size j = 1; j <= num_atoms; ++j ) {
			if ( predicate( *residues1[i], *residues2[i], j ) ) {
				core::Vector diff = residues1[i]->xyz(j) - residues2[i]->xyz(j);
				sum2 += diff.length_squared();
				natoms += 1;
			}
		}
	}
	return std::sqrt(sum2 / natoms);
}

/// @brief Select atoms for RMS via a predicate function/functor.
/// @details Calculates minimal rms over a subset of residues,
/// NOT allowing rotation/translation --
/// uses current coordinates as-is.
/// Parameter "predicate" should be a function pointer
/// [ or class that defines operator() ] with the following signature:
///
///   bool my_pred_func(Pose const & pose1, Pose const & pose2, Size resno, Size atomno);
/// the subset is a vector of all the residues, with true for those
/// over which to calculate the rms
///
/// It should return true if the atom should be included and false otherwise.
///
/// Example of use, to calculate C-alpha RMSD:
///   rmsd_no_super_subset(pose1, pose2, subset, is_protein_CA);
template< class T >
core::Real
rmsd_no_super_subset(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & subset,
	T* predicate
)
{
	core::Size const nres1 = pose1.total_residue();
	ASSERT_ONLY(core::Size const nres2 = pose2.total_residue();)
	assert( nres1 == nres2 );

	core::Real sum2( 0.0 );
	core::Size natoms( 0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {
	if ( subset( i ) ) {
		Size num_atoms ( pose1.residue(i).natoms() );
			if ( predicate == is_ligand_heavyatom ||
				predicate == is_polymer_heavyatom ||
				predicate == is_heavyatom ) {
				assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
			} else if ( num_atoms > pose2.residue(i).natoms() ){
				num_atoms = pose2.residue(i).natoms();
			}
			for ( core::Size j = 1; j <= num_atoms; ++j ) {
				if ( predicate( pose1, pose2, i, j ) ) {
					core::Vector diff = pose1.residue(i).xyz(j) - pose2.residue(i).xyz(j);
					sum2 += diff.length_squared();
					natoms += 1;
				}
			}
		}
	}
	return std::sqrt(sum2 / natoms);
}

/// @brief function to return the biggest deviation between an atom in a pair of poses,
template< class T >
core::Real
biggest_residue_deviation_no_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	T* predicate
){

	core::Size const nres1 = pose1.total_residue();
	ASSERT_ONLY(core::Size const nres2 = pose2.total_residue();)
	assert( nres1 == nres2 );

	core::Real biggest_dev2( 0.0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {

		Real res_dev2(0.0);
		Size res_count_atoms(0);
		Size num_atoms ( pose1.residue(i).natoms() );
		if ( predicate == is_ligand_heavyatom ||
			predicate == is_polymer_heavyatom ||
			predicate == is_heavyatom ) {
			assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
		} else if ( num_atoms > pose2.residue(i).natoms() ){
			num_atoms = pose2.residue(i).natoms();
		}
		for ( core::Size j = 1; j <= num_atoms; ++j ) {
			if ( predicate( pose1, pose2, i, j ) ) {
				core::Vector diff = pose1.residue(i).xyz(j) - pose2.residue(i).xyz(j);
				res_dev2 += diff.length_squared();
				res_count_atoms++;
			}
		}
		res_dev2 /= res_count_atoms;
		if( res_dev2 > biggest_dev2 ) biggest_dev2 = res_dev2;
	}
	return std::sqrt( biggest_dev2 );
}


/// @brief function to return the biggest deviation between an atom in a pair of poses,
/// @brief as specified by the predicate and the subset
template< class T >
core::Real
biggest_residue_deviation_no_super_subset(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & subset,
	T* predicate
){

	core::Size const nres1 = pose1.total_residue();
	core::Size const nres2 = pose2.total_residue();
	assert( nres1 == nres2 );

	core::Real biggest_dev2( 0.0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {
	if ( subset( i ) ) {
		Real res_dev2(0.0);
		Size res_count_atoms(0);
		Size num_atoms ( pose1.residue(i).natoms() );
			if ( predicate == is_ligand_heavyatom ||
				predicate == is_polymer_heavyatom ||
				predicate == is_heavyatom ) {
				assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
			} else if ( num_atoms > pose2.residue(i).natoms() ){
				num_atoms = pose2.residue(i).natoms();
			}
			for ( core::Size j = 1; j <= num_atoms; ++j ) {
				if ( predicate( pose1, pose2, i, j ) ) {
					core::Vector diff = pose1.residue(i).xyz(j) - pose2.residue(i).xyz(j);
					res_dev2 += diff.length_squared();
					res_count_atoms++;
				}
			}
			res_dev2 /= res_count_atoms;
			if( res_dev2 > biggest_dev2 ) biggest_dev2 = res_dev2;
		}
	}
	return std::sqrt( biggest_dev2 );
}

template< class T >
void
fill_rmsd_coordinates(
	int & natoms,
	ObjexxFCL::FArray2D< core::Real > & p1a,
	ObjexxFCL::FArray2D< core::Real > & p2a,
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	T* predicate
)
{
	core::Size const nres1 = pose1.total_residue();
	core::Size const nres2 = pose2.total_residue();
	//assert( nres1 == nres2 );

	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for ( core::Size i = 1; i <= std::min( nres1, nres2 ); ++i ) {
		//assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
		for ( core::Size j = 1; j <= pose1.residue(i).natoms(); ++j ) {
			if ( (*predicate)( pose1, pose2, i, j ) ) {
				p1_coords.push_back( pose1.residue(i).xyz(j) );
				p2_coords.push_back( pose2.residue(i).xyz(j) );
			}
		}
	}
	assert( p1_coords.size() == p2_coords.size() );

	natoms = p1_coords.size();
	p1a.dimension( 3, natoms );
	p2a.dimension( 3, natoms );
	for ( int i = 0; i < natoms; ++i ) {
		for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
			p1a(k+1,i+1) = p1_coords[i][k];
			p2a(k+1,i+1) = p2_coords[i][k];
		}
	}
}

/// @brief Select a subset atoms for Symmetric RMS via a predicate function/functor.
/// Example of use, to calculate C-alpha RMSD:
///   rmsd_with_super(pose1, pose2, is_protein_CA, subset);
template< class T >
core::Real
sym_rmsd_with_super_subset(
	core::pose::Pose const & native_pose,
	core::pose::Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & subset,
	T* predicate
)
{
	core::Size const nres1 = native_pose.total_residue();
	//core::Size const nres2 = pose2.total_residue();
	//assert( nres1 == nres2 );

	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for ( core::Size i = 1; i <= nres1; ++i ) {
		if ( subset( i ) ) {
			Size num_atoms ( native_pose.residue(i).natoms() );
			if ( predicate == is_ligand_heavyatom ||
				predicate == is_polymer_heavyatom ||
				predicate == is_heavyatom ) {
				assert( native_pose.residue(i).natoms() == pose2.residue(i).natoms() );
			} else if ( num_atoms > pose2.residue(i).natoms() ){
				num_atoms = pose2.residue(i).natoms();
			}
			for ( core::Size j = 1; j <= num_atoms; ++j ) {
				if ( predicate( native_pose, pose2, i, j ) ) {
					p1_coords.push_back( native_pose.residue(i).xyz(j) );
					p2_coords.push_back( pose2.residue(i).xyz(j) );
				}
			}
		}
	}
	assert( p1_coords.size() == p2_coords.size() );

	int const natoms = p1_coords.size();
	ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
	ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
	for ( int i = 0; i < natoms; ++i ) {
		for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
			p1a(k+1,i+1) = p1_coords[i][k];
			p2a(k+1,i+1) = p2_coords[i][k];
		}
	}


	assert( dynamic_cast< core::conformation::symmetry::SymmetricConformation const * >( &pose2.conformation() ) );

	core::conformation::symmetry::SymmetricConformation const & symm_conf (
		dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose2.conformation()) );
	core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

    int const N ( symm_info->subunits() );
    int const nres ( symm_info->num_total_residues_without_pseudo() );
    FArray2D< core::Real > p1a_shuffle( 3, nres );

  core::Real rms = 1e3; //Since fast_rms has not been evaluated yet
	int const num_atoms_subunit (natoms);
	std::vector< std::vector<int> > shuffle_map;
	create_shuffle_map_recursive_rms(std::vector<int>(), N,shuffle_map);
	for (int j=1; j < int (shuffle_map.size()); j++ ){
		for (int i=0; i < N; ++i ) {
			int const begin ( shuffle_map.at(j).at(i)*num_atoms_subunit*3);
			for ( int k = 0; k < num_atoms_subunit*3; ++k ) {
				int const begin_shuffled (i*num_atoms_subunit*3);
					p1a_shuffle[begin_shuffled+k] = p1a[begin+k];
			}
		}
		Real rms_shuffle = numeric::model_quality::rms_wrapper( natoms, p1a_shuffle, p2a );
		if ( rms_shuffle < rms ) {
			rms = rms_shuffle;
		}
	}
	return rms;
}

} // end namespace scoring
} // end namespace core

#endif
