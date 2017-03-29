// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/rms/rms_util.tmpl.hh
/// @brief  RMS stuff from rosetta++ that requires templates
/// @author James Thompson
/// @author Ian Davis
/// @date   Wed Aug 22 12:10:37 2007

#ifndef INCLUDED_core_scoring_rms_util_tmpl_HH
#define INCLUDED_core_scoring_rms_util_tmpl_HH

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/SequenceMapping.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <list>

namespace core {
namespace scoring {

/// @brief Select atoms for RMS via a predicate function/functor.
/// @details Calculates minimal rms, allowing rotation/translation for best fit.
/// The two vectors should be the same length, and taken pairwise should be a mapping.
/// Parameter "predicate" should be a function pointer
/// [ or class that defines operator() ] with the following signature:
///
///   bool my_pred_func(Pose const & pose1, Pose const & pose2, Size resno, Size atomno);
///
/// It should return true if the atom should be included and false otherwise.
///
/// The "resno" and "atomno" are calculated based on pose1 (which is what most predicates expect)
///
template< class T >
core::Real
rmsd_with_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	utility::vector1< core::Size > const & pose1_residues,
	utility::vector1< core::Size > const & pose2_residues,
	T* predicate
) {
	debug_assert( pose1_residues.size() == pose2_residues.size() );

	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for ( core::Size ii(1); ii <= pose1_residues.size(); ++ii ) {
		core::Size res1( pose1_residues[ii] );
		core::Size res2( pose2_residues[ii] );

		Size num_atoms ( pose1.residue(res1).natoms() );
		if ( predicate == is_ligand_heavyatom ||
				predicate == is_polymer_heavyatom ||
				predicate == is_heavyatom ) {
			//  debug_assert( pose1.residue(res1).nheavyatoms() == pose2.residue(res2).nheavyatoms() );
		} else if ( num_atoms > pose2.residue(res2).natoms() ) {
			num_atoms = pose2.residue(res2).natoms();
		}
		for ( core::Size jj = 1; jj <= num_atoms; ++jj ) {
			if ( jj <= pose1.residue(res1).natoms() &&
					jj <= pose2.residue(res2).natoms() ) {
				if ( predicate( pose1, pose2, res1, jj ) ) {
					p1_coords.push_back( pose1.residue(res1).xyz(jj) );
					p2_coords.push_back( pose2.residue(res2).xyz(jj) );
				}
			}
		}
	}
	debug_assert( p1_coords.size() == p2_coords.size() );

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

template< class T >
core::Real
rmsd_with_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	std::list< core::Size > const & subset_residues,
	T* predicate
) {

	utility::vector1< core::Size > res_vector( subset_residues.begin(), subset_residues.end() );

	return rmsd_with_super( pose1, pose2, res_vector, res_vector, predicate );
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
	debug_assert( pose1.size() == pose2.size() );
	utility::vector1< core::Size > all_residues;
	for ( core::Size n = 1; n <= pose1.size(); n++ ) {
		all_residues.push_back( n );
	}
	return rmsd_with_super( pose1, pose2, all_residues, all_residues, predicate );
}


/// @brief Select a subset atoms for RMS via a predicate function/functor.
/// @details Calculates minimal rms, allowing rotation/translation for best fit.
///  Same as above function, but allows a subset of residues over which to
///  superposition to be passed in
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
	debug_assert( subset.size() == pose1.size() );
	debug_assert( subset.size() == pose2.size() );

	utility::vector1< core::Size > selected_residues;
	for ( core::Size ii = 1; ii <= pose1.size(); ++ii ) {
		if ( subset( ii ) ) {
			selected_residues.push_back( ii );
		}
	}

	return rmsd_with_super( pose1, pose2, selected_residues, selected_residues, predicate );
}

/// @brief like function above, but uses sequence mapping,
/// i.e. sections of poses of different lengths can be compared
/// at the moment sorta assumes that residues at corresponding
/// positions have the same identity, mainly becaue the predicates
/// are structured that way...
template< class T >
core::Real
rmsd_with_super_subset(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & subset,
	core::id::SequenceMapping const & seqmap,
	T* predicate
)
{
	utility::vector1< core::Size > selected_residues1;
	utility::vector1< core::Size > selected_residues2;
	for ( core::Size ii = 1; ii <= pose1.size(); ++ii ) {
		if ( subset( ii ) ) {
			selected_residues1.push_back( ii );
			selected_residues2.push_back( seqmap[ ii ] );
		}
	}

	return rmsd_with_super( pose1, pose2, selected_residues1, selected_residues2, predicate );
}

/// @brief Select atoms for RMS via a predicate function/functor at given residues
/// @details Calculates rms, NOT allowing rotation/translation --
/// uses current coordinates as-is.
/// The two vectors should be the same length, and taken pairwise should be a mapping.
/// Parameter "predicate" should be a function pointer
/// [ or class that defines operator() ] with the following signature:
///
///   bool my_pred_func(Pose const & pose1, Pose const & pose2, Size resno, Size atomno);
///
/// It should return true if the atom should be included and false otherwise.
///
/// The "resno" and "atomno" are calculated based on pose1 (which is what most predicates expect)
///
template< class T >
core::Real
rmsd_no_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	utility::vector1< core::Size > const & pose1_residues,
	utility::vector1< core::Size > const & pose2_residues,
	T* predicate
)
{
	debug_assert( pose1_residues.size() == pose2_residues.size() );

	core::Real sum2( 0.0 );
	core::Size natoms( 0 );
	for ( core::Size ii = 1; ii <= pose1_residues.size(); ++ii ) {
		core::Size res1( pose1_residues[ii] );
		core::Size res2( pose2_residues[ii] );

		Size num_atoms ( pose1.residue(res1).natoms() );
		if ( predicate == is_ligand_heavyatom ||
				predicate == is_polymer_heavyatom ||
				predicate == is_heavyatom ) {
			// Should we just be checking heavy atoms here?
			debug_assert( pose1.residue(res1).natoms() == pose2.residue(res2).natoms() );
		} else if ( num_atoms > pose2.residue(res2).natoms() ) {
			// Should we been warning here?
			num_atoms = pose2.residue(res2).natoms();
		}
		for ( core::Size jj = 1; jj <= num_atoms; ++jj ) {
			if ( predicate( pose1, pose2, res1, jj ) ) {
				core::Vector diff = pose1.residue(res1).xyz(jj) - pose2.residue(res2).xyz(jj);
				sum2 += diff.length_squared();
				natoms += 1;
			}
		}
	}
	return std::sqrt(sum2 / natoms);
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
	debug_assert( pose1.size() == pose2.size() );

	utility::vector1<core::Size> all_res;
	all_res.reserve( pose1.size() );
	for ( core::Size ii(1); ii <= pose1.size(); ++ii ) {
		all_res.push_back( ii );
	}

	return rmsd_no_super( pose1, pose2, all_res, all_res, predicate );
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
	debug_assert( nres1 == nres2 );

	core::Real sum2( 0.0 );
	core::Size natoms( 0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {
		Size num_atoms ( residues1[i]->natoms() );
		if (  predicate == is_ligand_heavyatom_residues
				/* predicate == is_polymer_heavyatom || */
				/* predicate == is_heavyatom */ ) { // for Rotate.cc
			debug_assert( residues1[i]->natoms() == residues2[i]->natoms() );
		} else if ( num_atoms > residues2[i]->natoms() ) {
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
	debug_assert( subset.size() == pose1.size() );
	debug_assert( subset.size() == pose2.size() );

	utility::vector1< core::Size > selected_residues;
	for ( core::Size ii = 1; ii <= pose1.size(); ++ii ) {
		if ( subset( ii ) ) {
			selected_residues.push_back( ii );
		}
	}

	return rmsd_no_super( pose1, pose2, selected_residues, selected_residues, predicate );
}

/// @brief like function above, but uses sequence mapping,
/// i.e. sections of poses of different lengths can be compared
/// at the moment sorta assumes that residues at corresponding
/// positions have the same identity, mainly becaue the predicates
/// are structured that way...
template< class T >
core::Real
rmsd_no_super_subset(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	ObjexxFCL::FArray1D_bool const & subset,
	core::id::SequenceMapping const & seqmap,
	T* predicate
)
{
	utility::vector1< core::Size > selected_residues1;
	utility::vector1< core::Size > selected_residues2;
	for ( core::Size ii = 1; ii <= pose1.size(); ++ii ) {
		if ( subset( ii ) ) {
			selected_residues1.push_back( ii );
			selected_residues2.push_back( seqmap[ ii ] );
		}
	}

	return rmsd_no_super( pose1, pose2, selected_residues1, selected_residues2, predicate );
}


/// @brief function to return the biggest deviation between an atom in a pair of poses,
template< class T >
core::Real
biggest_residue_deviation_no_super(
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	T* predicate
){

	core::Size const nres1 = pose1.size();
	ASSERT_ONLY(core::Size const nres2 = pose2.size();)
	debug_assert( nres1 == nres2 );

	core::Real biggest_dev2( 0.0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {

		Real res_dev2(0.0);
		Size res_count_atoms(0);
		Size num_atoms ( pose1.residue(i).natoms() );
		if ( predicate == is_ligand_heavyatom ||
				predicate == is_polymer_heavyatom ||
				predicate == is_heavyatom ) {
			debug_assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
		} else if ( num_atoms > pose2.residue(i).natoms() ) {
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
		if ( res_dev2 > biggest_dev2 ) biggest_dev2 = res_dev2;
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

	core::Size const nres1 = pose1.size();
	ASSERT_ONLY(core::Size const nres2 = pose2.size());
	debug_assert( nres1 == nres2 );

	core::Real biggest_dev2( 0.0 );
	for ( core::Size i = 1; i <= nres1; ++i ) {
		if ( subset( i ) ) {
			Real res_dev2(0.0);
			Size res_count_atoms(0);
			Size num_atoms ( pose1.residue(i).natoms() );
			if ( predicate == is_ligand_heavyatom ||
					predicate == is_polymer_heavyatom ||
					predicate == is_heavyatom ) {
				debug_assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
			} else if ( num_atoms > pose2.residue(i).natoms() ) {
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
			if ( res_dev2 > biggest_dev2 ) biggest_dev2 = res_dev2;
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
	core::Size const nres1 = pose1.size();
	core::Size const nres2 = pose2.size();
	//debug_assert( nres1 == nres2 );

	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for ( core::Size i = 1; i <= std::min( nres1, nres2 ); ++i ) {
		if ( pose1.residue(i).is_virtual_residue() || pose2.residue(i).is_virtual_residue() ) continue;
		//debug_assert( pose1.residue(i).natoms() == pose2.residue(i).natoms() );
		for ( core::Size j = 1; j <= pose1.residue(i).natoms(); ++j ) {
			if ( (*predicate)( pose1, pose2, i, j ) ) {
				p1_coords.push_back( pose1.residue(i).xyz(j) );
				p2_coords.push_back( pose2.residue(i).xyz(j) );
			}
		}
	}
	debug_assert( p1_coords.size() == p2_coords.size() );

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
	core::Size const nres1 = native_pose.size();
	//core::Size const nres2 = pose2.size();
	//debug_assert( nres1 == nres2 );

	std::vector< core::Vector > p1_coords;
	std::vector< core::Vector > p2_coords;

	for ( core::Size i = 1; i <= nres1; ++i ) {
		if ( subset( i ) ) {
			Size num_atoms ( native_pose.residue(i).natoms() );
			if ( predicate == is_ligand_heavyatom ||
					predicate == is_polymer_heavyatom ||
					predicate == is_heavyatom ) {
				debug_assert( native_pose.residue(i).natoms() == pose2.residue(i).natoms() );
			} else if ( num_atoms > pose2.residue(i).natoms() ) {
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
	debug_assert( p1_coords.size() == p2_coords.size() );

	int const natoms = p1_coords.size();
	ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
	ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
	for ( int i = 0; i < natoms; ++i ) {
		for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
			p1a(k+1,i+1) = p1_coords[i][k];
			p2a(k+1,i+1) = p2_coords[i][k];
		}
	}


	debug_assert( dynamic_cast< core::conformation::symmetry::SymmetricConformation const * >( &pose2.conformation() ) );

	core::conformation::symmetry::SymmetricConformation const & symm_conf (
		dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose2.conformation()) );
	core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	int const N ( symm_info->subunits() );
	int const nres ( symm_info->num_total_residues_without_pseudo() );
	ObjexxFCL::FArray2D< core::Real > p1a_shuffle( 3, nres );

	core::Real rms = 1e3; //Since fast_rms has not been evaluated yet
	int const num_atoms_subunit (natoms);
	std::vector< std::vector<int> > shuffle_map;
	create_shuffle_map_recursive_rms(std::vector<int>(), N,shuffle_map);
	for ( int j=1; j < int (shuffle_map.size()); j++ ) {
		for ( int i=0; i < N; ++i ) {
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
