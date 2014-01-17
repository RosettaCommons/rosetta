// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rms_util.cc
/// @brief  RMS stuff from rosetta++
/// @author Christopher Miles (cmiles@uw.edu)
/// @author James Thompson
/// @author Ian Davis
/// @date   Wed Aug 22 12:10:37 2007
///

// Unit headers
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

// C/C++ headers
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <utility>

// External headers
#include <boost/assign.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/unordered/unordered_map.hpp>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/util.hh>
#include <core/id/types.hh>
#include <core/chemical/automorphism.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/TMscore.hh>

// Utility headers
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <numeric/util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/maxsub.hh>
#include <numeric/model_quality/rms.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <utility/exit.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

using namespace ObjexxFCL;

namespace core {
namespace scoring {

static basic::Tracer tr("core.scoring.rms_util");

core::Real gdtsc(const core::pose::Pose& ref,
                 const core::pose::Pose& mod,
                 const std::map<core::Size, core::Size>& residues) {
  using core::Real;
  using core::Size;
  using core::id::NamedAtomID;
  using numeric::xyzVector;
  using std::map;
  using std::string;

  static map<char, string> gdtsc_atom = boost::assign::map_list_of
      ('A',  "CA")
      ('C',  "SG")
      ('D', "OD2")
      ('E', "OE2")
      ('F',  "CZ")
      ('G',  "CA")
      ('H', "NE2")
      ('I', "CD1")
      ('K',  "NZ")
      ('L', "CD1")
      ('M',  "CE")
      ('N', "OD1")
      ('P',  "CG")
      ('Q', "OE1")
      ('R', "NH2")
      ('S',  "OG")
      ('T', "OG1")
      ('V', "CG1")
      ('W', "CH2")
      ('Y',  "OH");

  if (!ref.is_fullatom() || !mod.is_fullatom()) {
    tr.Warning << "Reference and model must be fullatom for gdtsc()" << std::endl;
    return -1;
  }

  // Retrieve ref, mod coordinates
  int expected_num_atoms = residues.size();
  int actual_num_atoms = 0;
  FArray2D<Real> coords_ref(3, expected_num_atoms);
  FArray2D<Real> coords_mod(3, expected_num_atoms);

  int count = 1;
  for (map<Size, Size>::const_iterator i = residues.begin(); i != residues.end(); ++i, ++count) {
    const Size ref_idx = i->first;
    const Size mod_idx = i->second;
    const char ref_residue = ref.residue(ref_idx).name1();
    const char mod_residue = mod.residue(mod_idx).name1();

    if (ref_residue != mod_residue) {
      tr.Warning << "Reference and model must have identical sequences for gdtha-- "
                 << ref_residue << " != " << mod_residue << std::endl;
      continue;
    }

    ++actual_num_atoms;
    const NamedAtomID ref_atom(gdtsc_atom[ref_residue], ref_idx);
    const xyzVector<Real>& xyz_ref = ref.xyz(ref_atom);
    coords_ref(1, count) = xyz_ref.x();
    coords_ref(2, count) = xyz_ref.y();
    coords_ref(3, count) = xyz_ref.z();

    const NamedAtomID mod_atom(gdtsc_atom[mod_residue], mod_idx);
    const xyzVector<Real>& xyz_mod = mod.xyz(mod_atom);
    coords_mod(1, count) = xyz_mod.x();
    coords_mod(2, count) = xyz_mod.y();
    coords_mod(3, count) = xyz_mod.z();
  }

  // Calculate maxsub over several distance thresholds
  Real sum = 0;
  Size num_dists = 10;

  for (Size i = 1; i <= num_dists; ++i) {
    Real dist_threshold = 0.5 * i;  // 0.5, 1.0, ...

    int nali;
    double mxrms, mxpsi, mxzscore, mxscore, mxeval;
    numeric::model_quality::maxsub(
        actual_num_atoms, coords_ref, coords_mod,
        mxrms, mxpsi, nali, mxzscore, mxeval, mxscore,
        dist_threshold, dist_threshold);

    Real pct_residues = static_cast<Real>(nali) / static_cast<Real>(actual_num_atoms);
    sum += pct_residues;
  }

  return sum / num_dists;
}

core::Real gdtha(const core::pose::Pose& ref,
                 const core::pose::Pose& mod,
                 const std::map<core::Size, core::Size>& residues) {
  using core::Real;
  using core::Size;
  using core::id::NamedAtomID;
  using numeric::xyzVector;

  // Retrieve ref, mod coordinates
  int expected_num_atoms = residues.size();
  int actual_num_atoms = 0;
  FArray2D<Real> coords_ref(3, expected_num_atoms);
  FArray2D<Real> coords_mod(3, expected_num_atoms);

  int count = 1;
  for (std::map<Size, Size>::const_iterator i = residues.begin(); i != residues.end(); ++i, ++count) {
    const Size ref_idx = i->first;
    const Size mod_idx = i->second;
    const char ref_residue = ref.residue(ref_idx).name1();
    const char mod_residue = mod.residue(mod_idx).name1();

    if (ref_residue != mod_residue) {
      tr.Warning << "Reference and model must have identical sequences for gdtha-- "
                 << ref_residue << " != " << mod_residue << std::endl;
      continue;
    }

    ++actual_num_atoms;
    const NamedAtomID ref_atom("CA", ref_idx);
    const xyzVector<Real>& xyz_ref = ref.xyz(ref_atom);
    coords_ref(1, count) = xyz_ref.x();
    coords_ref(2, count) = xyz_ref.y();
    coords_ref(3, count) = xyz_ref.z();

    const NamedAtomID mod_atom("CA", mod_idx);
    const xyzVector<Real>& xyz_mod = mod.xyz(mod_atom);
    coords_mod(1, count) = xyz_mod.x();
    coords_mod(2, count) = xyz_mod.y();
    coords_mod(3, count) = xyz_mod.z();
  }

  // Calculate maxsub over several distance thresholds
  Real dists[] = {0.5, 1.0, 2.0, 4.0};
  Size num_dists = 4;
  Real sum = 0;

  for (Size i = 0; i < num_dists; ++i) {
    Real dist_threshold = dists[i];

    int nali;
    double mxrms, mxpsi, mxzscore, mxscore, mxeval;
    numeric::model_quality::maxsub(
        actual_num_atoms, coords_ref, coords_mod,
        mxrms, mxpsi, nali, mxzscore, mxeval, mxscore,
        dist_threshold, dist_threshold);

    Real pct_residues = static_cast<Real>(nali) / static_cast<Real>(actual_num_atoms);
    sum += pct_residues;
  }

  return sum / num_dists;
}

void invert_exclude_residues( Size nres, utility::vector1<int> const& exclude_list, ResidueSelection& residue_selection ) {
	residue_selection.clear();

	for( Size ir = 1; ir <= nres; ++ir ) {
		bool exclude_residue = false;
		for( Size ex = 1; ex <= exclude_list.size(); ex ++ ){
			if( int(exclude_list[ex]) == int(ir) ) {
				exclude_residue = true;
				break;
			}
		}

		if ( !exclude_residue ) {
			residue_selection.push_back( ir );
		}
	} // for ( Size ir = 1; ir <= native_pose.total_residue(); ++ir )
}

ResidueSelection invert_exclude_residues( core::Size nres, utility::vector1<int> const& exclude_list ) {
	ResidueSelection tmp;
	invert_exclude_residues( nres, exclude_list, tmp );
	return tmp;
}

Real native_CA_rmsd( const core::pose::Pose & native_pose, const core::pose::Pose & pose ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ in::file::native_exclude_res ].user() ) {
		ResidueSelection residues;
		invert_exclude_residues( native_pose.total_residue(), option[ in::file::native_exclude_res ](), residues );
		return core::scoring::CA_rmsd( native_pose, pose, residues );
	} else {
		// no residues excluded from the native.
		return core::scoring::CA_rmsd( native_pose, pose );
	}
} // native_CA_rmsd

Real native_CA_gdtmm( const core::pose::Pose & native_pose, const core::pose::Pose & pose ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ in::file::native_exclude_res ].user() ) {
		ResidueSelection residues;
		invert_exclude_residues( native_pose.total_residue(), option[ in::file::native_exclude_res ](), residues );
		return core::scoring::CA_gdtmm( native_pose, pose, residues );
	} else {
		// no residues excluded from the native.
		return core::scoring::CA_gdtmm( native_pose, pose );
	}
} // native_CA_gdtmm

/// @details Just iterates over all automorphisms for this residue type and chooses the minimum RMS.
core::Real
automorphic_rmsd(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	bool superimpose
)
{
	using namespace core;
	using namespace core::chemical;
	using namespace core::conformation;
	// name() and total number of atoms may actually be different, if we're comparing e.g. tautomers
	if( rsd1.type().name3() != rsd2.type().name3() ) utility_exit_with_message("Residue type name3 mismatch");
	if( rsd1.nheavyatoms()  != rsd2.nheavyatoms()  ) utility_exit_with_message("Residue number-of-heavy-atoms mismatch");
	core::Real best_rms = 1e99;
	int counter = 0;
	// Make atom-number translation table
	ResidueTypeCOP rsd1_type( &(rsd1.type()) );
	AutomorphismIterator ai( rsd1_type );
	AtomIndices old2new( ai.next() );
	// For each permutation of automorphisms...
	while( old2new.size() > 0 ) {
		counter++;
		if( counter%10000 == 0 ) tr.Info << counter << " so far..." << std::endl;

		// Print out translation table for debugging
		//std::cout << "[";
		//for(Size i = 1; i <= old2new.size(); ++i) std::cout << " " << old2new[i];
		//std::cout << " ]\n";
		//for(Size j = 1; j <= old2new.size(); ++j) std::cout << "  " << j << " --> " << old2new[j] << "  /  " << rsd1.type().atom_name(j) << " --> " << rsd1.type().atom_name(old2new[j]) << "\n";

		// Compute rmsd
		if( superimpose ) {
			utility::vector0< core::Vector > p1_coords;
			utility::vector0< core::Vector > p2_coords;
			for ( core::Size j = 1; j <= rsd1.type().natoms(); ++j ) {
				if ( !rsd1.atom_type(j).is_hydrogen() ) {
					// This is the step where we effectively re-assign atom names
					// in hopes of reducing RMS (e.g. by "flipping" a phenyl ring).
					p1_coords.push_back( rsd1.xyz( j ) );
					p2_coords.push_back( rsd2.xyz( old2new[j] ) );
				}
			}
			runtime_assert( p1_coords.size() == p2_coords.size() );
			int const natoms = p1_coords.size();
			ObjexxFCL::FArray2D< core::Real > p1a( 3, natoms );
			ObjexxFCL::FArray2D< core::Real > p2a( 3, natoms );
			for ( int j = 0; j < natoms; ++j ) {
				for ( int k = 0; k < 3; ++k ) { // k = X, Y and Z
					p1a(k+1,j+1) = p1_coords[j][k];
					p2a(k+1,j+1) = p2_coords[j][k];
				}
			}
			core::Real const curr_rms = numeric::model_quality::rms_wrapper( natoms, p1a, p2a );
			// Check vs. minimum rmsd
			if( curr_rms < best_rms ) {
				//tr.Debug << "New rms of " << curr_rms << " beats previous best of " << best_rms << std::endl;
				best_rms = curr_rms;
			}
		} else { // don't superimpose
			core::Real sum2( 0.0 );
			core::Size natoms( 0 );
			for ( core::Size j = 1; j <= rsd1.type().natoms(); ++j ) {
				if ( !rsd1.atom_type(j).is_hydrogen() ) {
					// This is the step where we effectively re-assign atom names
					// in hopes of reducing RMS (e.g. by "flipping" a phenyl ring).
					core::Vector diff = rsd1.xyz( j ) - rsd2.xyz( old2new[j] );
					sum2 += diff.length_squared();
					natoms += 1;
				}
			}
			core::Real const curr_rms = std::sqrt(sum2 / natoms);
			// Check vs. minimum rmsd
			if( curr_rms < best_rms ) {
				//tr.Debug << "New rms of " << curr_rms << " beats previous best of " << best_rms << std::endl;
				best_rms = curr_rms;
			}
		}
		old2new = ai.next();
	} // done checking all automorphisms
	tr.Debug << counter << " automorphisms from iterator; best rms is " << best_rms << std::endl;
	return best_rms;
}


//////////////////////////////////////////////////////////////////////////////
// Predicate functions to use with rmsd_no_super() and rmsd_with_super()

bool
is_protein_CA(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	return rsd.is_protein() && rsd.has("CA") && rsd.atom_index("CA") == atomno;
}

bool
is_protein_CA_or_CB(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	return rsd.is_protein() && ( ( rsd.has("CA") && rsd.atom_index("CA") == atomno ) || ( rsd.has("CB") && rsd.atom_index("CB") == atomno ) );
}

bool
is_protein_backbone(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	return  (rsd.is_protein() && ( rsd.has("CA") && rsd.atom_index("CA") == atomno )) ||
		( rsd.has("N") && rsd.atom_index("N") == atomno ) ||
		( rsd.has("C") && rsd.atom_index("C") == atomno );
}

bool
is_protein_backbone_including_O(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	return  (rsd.is_protein() && ( rsd.has("CA") && rsd.atom_index("CA") == atomno )) ||
		( rsd.has("N") && rsd.atom_index("N") == atomno ) ||
		( rsd.has("C") && rsd.atom_index("C") == atomno ) ||
		( rsd.has("O") && rsd.atom_index("O") == atomno );
}

bool
is_protein_sidechain_heavyatom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::chemical::ResidueType const & rsd = pose1.residue_type(resno);
	return rsd.is_protein() && (rsd.first_sidechain_atom() <= atomno ) && ( !rsd.atom_is_hydrogen( atomno ) );
}

bool
is_ligand_heavyatom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	return !rsd.is_polymer() && !rsd.atom_is_hydrogen(atomno);
}

bool
is_ligand_heavyatom_residues(
		core::conformation::Residue const & residue1,
		core::conformation::Residue const &, // residue2
		core::Size atomno
){
	return !residue1.is_polymer() && !residue1.atom_is_hydrogen(atomno);
}

bool
is_polymer_heavyatom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	return rsd.is_polymer() && !rsd.atom_is_hydrogen(atomno);
}

bool
is_heavyatom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	return !rsd.atom_is_hydrogen(atomno);
}

bool
is_scatom(
  core::pose::Pose const & pose1,
  core::pose::Pose const & ,//pose2,
  core::Size resno,
  core::Size atomno
)
{
  core::conformation::Residue const & rsd = pose1.residue(resno);
  return !rsd.atom_is_backbone(atomno);
}

bool
is_nbr_atom(
	core::pose::Pose const & pose1,
	core::pose::Pose const & ,//pose2,
	core::Size resno,
	core::Size atomno
)
{
	core::conformation::Residue const & rsd = pose1.residue(resno);
	return rsd.nbr_atom() == atomno;
}

/////////////////////////////////////////////
// Predicate classes for more complex control

bool ResRangePredicate::operator()(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::Size resno,
		core::Size atomno) const {
	if ( resno < start_ || resno > end_ ) { return false; }
	else { return (*pred_)(pose1, pose2, resno, atomno); }
}

bool SelectedResPredicate::operator()(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::Size resno,
		core::Size atomno) const {
	if ( std::find( selected_.begin(), selected_.end(), resno ) == selected_.end() ) { return false; }
	else { return (*pred_)(pose1, pose2, resno, atomno); }
}

bool ExcludedResPredicate::operator()(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::Size resno,
		core::Size atomno) const {
	if ( std::find( excluded_.begin(), excluded_.end(), resno ) != excluded_.end() ) { return false; }
	else { return (*pred_)(pose1, pose2, resno, atomno); }
}


/////////////////////////////////////////////////////////////////////////

core::Real
CA_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
) {
	using namespace core;
	Size start ( 1 );
	Size end ( std::min( pose1.total_residue(), pose2.total_residue() ) );
	return CA_rmsd( pose1, pose2, start, end );
} // CA_rmsd

core::Real
CA_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	Size start,
	Size end
) {
	PROF_START( basic::CA_RMSD_EVALUATION );
	using namespace core;
	// copy coords into Real arrays
	int natoms;
	FArray2D< core::Real > p1a;
	FArray2D< core::Real > p2a;
	PredicateOP pred( new ResRangePredicate( start, end, new IsProteinCAPredicate ) );
	fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, pred() );

	if ( (int) (end - start + 1) > natoms ) { tr.Warning << "WARNING: In CA_rmsd, residue range " << start << " to " << end
			<< " requested but only " << natoms << " protein CA atoms found." << std::endl; }

	Real rms = numeric::model_quality::rms_wrapper( natoms, p1a, p2a );
	if(rms < 0.00001) rms = 0.0;
	PROF_STOP( basic::CA_RMSD_EVALUATION );
	return rms;
} // CA_rmsd

/// @detail Populates the output parameter with the xyz coordinates of
/// a subset of <pose>'s CA atoms, which are specified in <residues>
void retrieve_coordinates(const core::pose::Pose& pose,
													const utility::vector1<core::Size>& residues,
													FArray2D<core::Real>* coords) {
	using core::Real;
	using core::Size;
	using core::id::NamedAtomID;
	using numeric::xyzVector;
	using utility::vector1;

	coords->dimension(3, residues.size());

	for (Size i = 1; i <= residues.size(); ++i) {
		const NamedAtomID id("CA", residues[i]);
		const xyzVector<Real>& xyz = pose.xyz(id);
		(*coords)(1, i) = xyz.x();
		(*coords)(2, i) = xyz.y();
		(*coords)(3, i) = xyz.z();
	}
}

core::Real CA_rmsd(const core::pose::Pose& pose1,
									 const core::pose::Pose& pose2,
									 const std::map<core::Size, core::Size>& residues) {
	using core::Real;
	using core::Size;
	using utility::vector1;

	vector1<Size> residues_1;  // residues in pose1
	vector1<Size> residues_2;  // residues in pose2
	for (std::map<Size, Size>::const_iterator i = residues.begin(); i != residues.end(); ++i) {
		Size res_1 = i->first;
		Size res_2 = i->second;
		residues_1.push_back(res_1);
		residues_2.push_back(res_2);
	}

	FArray2D<Real> p1;  // coordinates of CA atoms of selected residues in pose1
	FArray2D<Real> p2;  // coordinates of CA atoms of selected residues in pose2
	retrieve_coordinates(pose1, residues_1, &p1);
	retrieve_coordinates(pose2, residues_2, &p2);

	return numeric::model_quality::rms_wrapper(residues.size(), p1, p2);
}

core::Real CA_gdtmm(const core::pose::Pose& pose1,
										const core::pose::Pose& pose2,
										const std::map<core::Size, core::Size>& residues) {
  using core::Real;
	using core::Size;
	using utility::vector1;

	vector1<Size> residues_1;  // residues in pose1
	vector1<Size> residues_2;  // residues in pose2
	for (std::map<Size, Size>::const_iterator i = residues.begin(); i != residues.end(); ++i) {
		Size res_1 = i->first;
		Size res_2 = i->second;
		residues_1.push_back(res_1);
		residues_2.push_back(res_2);
	}

	FArray2D<Real> p1;  // coordinates of CA atoms of selected residues in pose1
	FArray2D<Real> p2;  // coordinates of CA atoms of selected residues in pose2
	retrieve_coordinates(pose1, residues_1, &p1);
	retrieve_coordinates(pose2, residues_2, &p2);

	Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
	return xyz_gdtmm(p1, p2, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4);
}

core::Real
CA_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	Size start,
	Size end,
	utility::vector1< Size > const& exclude
) {
	PROF_START( basic::CA_RMSD_EVALUATION );
	using namespace core;
	// copy coords into Real arrays
	int natoms;
	FArray2D< core::Real > p1a;//( 3, pose1.total_residue() );
	FArray2D< core::Real > p2a;//( 3, pose2.total_residue() );
	PredicateOP pred( new ResRangePredicate( start, end, new ExcludedResPredicate( exclude, new IsProteinCAPredicate )) );
	fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, pred() );

	// Calc rms
	Real rms = numeric::model_quality::rms_wrapper( natoms, p1a, p2a );
	if(rms < 0.00001) rms = 0.0;
	PROF_STOP( basic::CA_RMSD_EVALUATION );
	return rms;
} // CA_rmsd

core::Real
bb_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
) {

	using namespace core;
	Size start ( 1 );
	Size end ( std::min( pose1.total_residue(), pose2.total_residue() ) );
	utility::vector1<Size> blank;
	return bb_rmsd( pose1, pose2, start, end, blank );
} // bb_rmsd

core::Real
bb_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	Size,
	Size,
	utility::vector1< Size > const&
) {
	using namespace core;
	Real rms = rmsd_with_super( pose1, pose2, is_protein_backbone );
	return rms;
} // bb_rmsd

core::Real
bb_rmsd_including_O(
        const core::pose::Pose & pose1,
        const core::pose::Pose & pose2
) {

        using namespace core;
        Size start ( 1 );
        Size end ( std::min( pose1.total_residue(), pose2.total_residue() ) );
        utility::vector1<Size> blank;
        return bb_rmsd_including_O( pose1, pose2, start, end, blank );
} // bb_rmsd_including_O

core::Real
bb_rmsd_including_O(
        const core::pose::Pose & pose1,
        const core::pose::Pose & pose2,
        Size,
        Size,
        utility::vector1< Size > const&
) {
        using namespace core;
        Real rms = rmsd_with_super( pose1, pose2, is_protein_backbone_including_O );
        return rms;
} // bb_rmsd_including_O

core::Real
CA_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	std::list< Size > residue_selection //the std::list can be sorted!
) {
	PROF_START( basic::CA_RMSD_EVALUATION );
	using namespace core;
	// copy coords into Real arrays
	int natoms;
	FArray2D< core::Real > p1a;//( 3, pose1.total_residue() );
	FArray2D< core::Real > p2a;//( 3, pose2.total_residue() );
	PredicateOP pred( new SelectedResPredicate( residue_selection, new IsProteinCAPredicate ) );
	fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, pred() );

	if ( (int) residue_selection.size() > natoms ) { tr.Warning << "WARNING: In CA_rmsd " << residue_selection.size()
				<< " residues selected but only " << natoms << " protein CA atoms found." << std::endl; }

	// Calc rms
	PROF_STOP( basic::CA_RMSD_EVALUATION );
	return numeric::model_quality::rms_wrapper( natoms, p1a, p2a );
} // CA_rmsd

core::Real
all_atom_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
) {

	using namespace core;
	Real rms = rmsd_with_super( pose1, pose2, is_heavyatom );
	return rms;
} // all atom rmsd

core::Real
all_scatom_rmsd_nosuper(
  const core::pose::Pose & pose1,
  const core::pose::Pose & pose2
) {

  using namespace core;
  Real rms = rmsd_no_super( pose1, pose2, is_scatom );
  return rms;
} // sc all atom rmsd no super

core::Real
all_atom_rmsd_nosuper(
  const core::pose::Pose & pose1,
  const core::pose::Pose & pose2
) {

  using namespace core;
  Real rms = rmsd_no_super( pose1, pose2, is_heavyatom);
  return rms;
} // all atom rmsd no super

core::Real
all_atom_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	std::list< Size > residue_selection //the std::list can be sorted!
)
{
	using namespace core;
	Real rms = rmsd_with_super( pose1, pose2, residue_selection, is_heavyatom );
	return rms;
} // CA_rmsd

core::Real
nbr_atom_rmsd(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2
) {

	using namespace core;
	Real rms = rmsd_with_super( pose1, pose2, is_nbr_atom );
	return rms;
} // nbr_atom_rmsd

// fill_rmsd_coordinates() is in rms_util.tmpl.hh

int
CA_maxsub(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	Real  rms
)
{
	const int nres1( pose1.total_residue() );

	static std::string atom_name = "CA";
	int natoms(0);
	FArray2D< double > p1a( 3, nres1 );
	FArray2D< double > p2a( 3, nres1 );
	// fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, "CA" );
	fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, is_protein_CA );

	double mxrms, mxpsi, mxzscore, mxscore, mxeval;
	int nali;
	numeric::model_quality::maxsub( natoms, p1a, p2a, mxrms, mxpsi, nali, mxzscore, mxeval, mxscore, rms );
	//logeval = std::log(mxeval);
	return nali;
}

int
CA_maxsub(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	std::list< Size > residue_selection, //the std::list can be sorted!
	Real rms
) {

	using namespace core;
	// copy coords into Real arrays
	int natoms;
	FArray2D< core::Real > p1a;//( 3, pose1.total_residue() );
	FArray2D< core::Real > p2a;//( 3, pose2.total_residue() );
	PredicateOP pred( new SelectedResPredicate( residue_selection, new IsProteinCAPredicate ) );
	fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, pred() );

	if ( (int) residue_selection.size() > natoms ) { tr.Warning << "WARNING: In CA_maxsub " << residue_selection.size()
			<< " residues selected but only " << natoms << " protein CA atoms found." << std::endl; }

	double mxrms, mxpsi, mxzscore, mxscore, mxeval;
	int nali;
	numeric::model_quality::maxsub( natoms, p1a, p2a, mxrms, mxpsi, nali, mxzscore, mxeval, mxscore, rms );
	//logeval = std::log(mxeval);
	return nali;
}

int xyz_maxsub(
	FArray2D< core::Real > p1a,
	FArray2D< core::Real > p2a,
	int natoms
) {
	double mxrms, mxpsi, mxzscore, mxscore, mxeval;
	int nali;
	numeric::model_quality::maxsub( natoms, p1a, p2a, mxrms, mxpsi, nali, mxzscore, mxeval, mxscore );
	return nali;
}

int
CA_maxsub_by_subset(
	const core::pose::Pose & pose1,
	const core::pose::Pose & pose2,
	utility::vector1< bool > //subset
)
{
	const int nres1( pose1.total_residue() );

	static std::string atom_name = "CA";
	int natoms(0);
	FArray2D< double > p1a( 3, nres1 );
	FArray2D< double > p2a( 3, nres1 );
	// fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, "CA" );
	fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, is_protein_CA );

	double mxrms, mxpsi, mxzscore, mxscore, mxeval;
	int nali;
	numeric::model_quality::maxsub( natoms, p1a, p2a, mxrms, mxpsi, nali, mxzscore, mxeval, mxscore );
	//logeval = std::log(mxeval);
	return nali;
} // CA_maxsub_by_subset

core::Real
CA_gdtmm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2,
	std::list< Size > residue_selection, //the std::list can be sorted! -- note std::sort can be applied to vectors
	core::Real& m_1_1,
	core::Real& m_2_2,
	core::Real& m_3_3,
	core::Real& m_4_3,
	core::Real& m_7_4
) {
	int natoms;
	FArray2D< core::Real > p1a( 3, pose1.total_residue() );
	FArray2D< core::Real > p2a( 3, pose2.total_residue() );
	PredicateOP pred( new SelectedResPredicate( residue_selection, new IsProteinCAPredicate ) );
	fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, pred() );

	if ( (int) residue_selection.size() > natoms ) { tr.Warning << "WARNING: In CA_gdtmm " << residue_selection.size()
			<< " residues selected but only " << natoms << " protein CA atoms found." << std::endl; }

	core::Real gdtmm = xyz_gdtmm( p1a, p2a, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
	return gdtmm;
}

core::Real
CA_gdtmm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2,
	core::Real& m_1_1,
	core::Real& m_2_2,
	core::Real& m_3_3,
	core::Real& m_4_3,
	core::Real& m_7_4
) {
	int natoms;
	FArray2D< core::Real > p1a( 3, pose1.total_residue() );
	FArray2D< core::Real > p2a( 3, pose2.total_residue() );
	fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, is_protein_CA );

	core::Real gdtmm = xyz_gdtmm( p1a, p2a, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
	return gdtmm;
}

core::Real
xyz_gdtmm(
	FArray2D< core::Real > p1a,
	FArray2D< core::Real > p2a
) {
	core::Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
	core::Real gdtmm = xyz_gdtmm( p1a, p2a, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
	return gdtmm;
}

core::Real
xyz_gdtmm(
	FArray2D< core::Real > p1a,
	FArray2D< core::Real > p2a,
	core::Real& m_1_1,
	core::Real& m_2_2,
	core::Real& m_3_3,
	core::Real& m_4_3,
	core::Real& m_7_4
) {
	int natoms = p1a.size2();
	runtime_assert( (p1a.size() == p2a.size()) && (p1a.size2() == p2a.size2()) );
	double mxrms, mxpsi, mxzscore, mxscore, mxeval;
	int nali;

	core::Real rmstol, disttol;
	rmstol = 1.0;
	disttol = 1.0;
	tr.Trace << "call maxsub with rmstol " << rmstol << " and disttol " << disttol << std::endl;
	numeric::model_quality::maxsub( natoms, p1a, p2a, mxrms, mxpsi, nali, mxzscore, mxeval, mxscore, rmstol, disttol );
	m_1_1 = core::Real( nali ) / core::Real( natoms );

	rmstol = 2.0;
	disttol = 2.0;
	tr.Trace << "call maxsub with rmstol " << rmstol << " and disttol " << disttol << std::endl;
	numeric::model_quality::maxsub( natoms, p1a, p2a, mxrms, mxpsi, nali, mxzscore, mxeval, mxscore, rmstol, disttol );
	m_2_2 = core::Real( nali ) / core::Real( natoms );

	rmstol = 3.0;
	disttol = 3.0;
	tr.Trace << "call maxsub with rmstol " << rmstol << " and disttol " << disttol << std::endl;
	numeric::model_quality::maxsub( natoms, p1a, p2a, mxrms, mxpsi, nali, mxzscore, mxeval, mxscore, rmstol, disttol );
	m_3_3 = core::Real( nali ) / core::Real( natoms );

	rmstol = 3.0;
	disttol = 4.0;
	tr.Trace << "call maxsub with rmstol " << rmstol << " and disttol " << disttol << std::endl;
	numeric::model_quality::maxsub( natoms, p1a, p2a, mxrms, mxpsi, nali, mxzscore, mxeval, mxscore, rmstol, disttol );
	m_4_3 = core::Real( nali ) / core::Real( natoms );

	rmstol = 4.0;
	disttol = 7.0;
	tr.Trace << "call maxsub with rmstol " << rmstol << " and disttol " << disttol << std::endl;
	numeric::model_quality::maxsub( natoms, p1a, p2a, mxrms, mxpsi, nali, mxzscore, mxeval, mxscore, rmstol, disttol );
	m_7_4 = core::Real( nali ) / core::Real( natoms );

	return (m_1_1 + m_2_2 + m_3_3 + m_4_3 + m_7_4 )/5.0;
}

core::Real
CA_gdtmm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2
) {
	core::Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
	return CA_gdtmm( pose1, pose2, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
}

core::Real
CA_gdtmm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2,
	std::list< Size> residue_selection
) {
	core::Real m_1_1, m_2_2, m_3_3, m_4_3, m_7_4;
	return CA_gdtmm( pose1, pose2, residue_selection, m_1_1, m_2_2, m_3_3, m_4_3, m_7_4 );
}

///////// GDT values by TMscore
void
CA_gdttm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2,
	core::Real &gdttm_score,
	core::Real &gdtha_score,
	std::list< Size > residue_selection //the std::list can be sorted! -- note std::sort can be applied to vectors
) {
	int natoms;
	FArray2D< core::Real > p1a( 3, pose1.total_residue() );
	FArray2D< core::Real > p2a( 3, pose2.total_residue() );
	PredicateOP pred( new SelectedResPredicate( residue_selection, new IsProteinCAPredicate ) );
	fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, pred() );

	if ( (int) residue_selection.size() > natoms ) { tr.Warning << "WARNING: In CA_gdtmm " << residue_selection.size()
			<< " residues selected but only " << natoms << " protein CA atoms found." << std::endl; }

	TMscore tm( p1a );
	tm.apply( p2a );
	gdttm_score = tm.get_GDTTS();
	gdtha_score = tm.get_GDTHA();
}

void
CA_gdttm(
	core::pose::Pose const& pose1,
	core::pose::Pose const& pose2,
	core::Real &gdttm_score,
	core::Real &gdtha_score
) {
	int natoms;
	FArray2D< core::Real > p1a( 3, pose1.total_residue() );
	FArray2D< core::Real > p2a( 3, pose2.total_residue() );
	fill_rmsd_coordinates( natoms, p1a, p2a, pose1, pose2, is_protein_CA );

	TMscore tm( p1a );
	tm.apply( p2a );
	gdttm_score = tm.get_GDTTS();
	gdtha_score = tm.get_GDTHA();
}

void
xyz_gdttm(
	FArray2D< core::Real > p1a,
	FArray2D< core::Real > p2a,
	core::Real &gdttm_score,
	core::Real &gdtha_score
) {

	TMscore tm( p1a );
	tm.apply( p2a );
	gdttm_score = tm.get_GDTTS();
	gdtha_score = tm.get_GDTHA();
}

/// @details  Superimpose mod_pose onto ref_pose using the AtomID mapping, which maps atoms in mod_pose onto
/// atoms in ref_pose. Returns rmsd over alignment.
///
/// @note  Atoms in mod_pose whose ids map to bogus atom ids will not be used in the fitting
///
/// Usage example: superimpose pose1 onto pose2 by mapping C-alphas of residue 10-30 onto residues 20-40
///
///		id::AtomID_Map< id::AtomID > atom_map;
///		id::initialize( atom_map, pose1, id::BOGUS_ATOM_ID ); // maps every atomid to bogus atom
///
///   for ( Size i=10; i<=30; ++i ) {
///     id::AtomID const id1( pose1.residue(i).atom_index("CA"), i );
///     id::AtomID const id2( pose2.residue(i+10).atom_index("CA"), i+10 );
///     atom_map[ id1 ] = id2;
///   }
///   superimpose_pose( pose1, pose2, atom_map );
///

Real
superimpose_pose(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose,
	id::AtomID_Map< id::AtomID > const & atom_map // from mod_pose to ref_pose
)
{
	pose::MiniPose mini_ref_pose( ref_pose );  //minipose is a lightweight pose with just xyz positions (& fold tree)
	return superimpose_pose(mod_pose, mini_ref_pose, atom_map );
}

Real
superimpose_pose(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose,
	std::map< id::AtomID, id::AtomID > const & atom_map // from mod_pose to ref_pose
)
{
	id::AtomID_Map< id::AtomID > const atom_id_map = convert_from_std_map( atom_map, mod_pose );
	pose::MiniPose mini_ref_pose( ref_pose );  //minipose is a lightweight pose with just xyz positions (& fold tree)
	return superimpose_pose(mod_pose, mini_ref_pose, atom_id_map );
}

Real
superimpose_pose(
	pose::Pose & mod_pose,
	pose::MiniPose const & ref_pose,
	id::AtomID_Map< id::AtomID > const & atom_map // from mod_pose to ref_pose
)
{
	using namespace numeric::model_quality;
	using namespace id;

	/// how many atoms in mod_pose?
	Size natoms(0);
	for ( Size i=1; i<= mod_pose.total_residue(); ++i ) natoms += mod_pose.residue(i).natoms();

	// pack coords into local arrays for passing into the rms fitting routine
	FArray2D_double xx1(3,natoms);
	FArray2D_double xx2(3,natoms);
	FArray1D_double wt(natoms);


	Size nsup(0);
	{ // pack coordinates into the local arrays
		Size atomno(0);
		Vector const zero_vector(0.0);
		for ( Size i=1; i<= mod_pose.total_residue(); ++i ) {
			for ( Size j=1; j<= mod_pose.residue(i).natoms(); ++j ) {
				++atomno;
				AtomID const & aid( atom_map[ id::AtomID( j,i) ] );
				Vector const & x1( aid.valid() ? ref_pose.xyz( aid ) : zero_vector );
				Vector const & x2( mod_pose.residue(i).xyz(j) );
				wt( atomno ) = ( aid.valid() ? 1.0 : 0.0 );
				if ( aid.valid() ) ++nsup;
				for ( Size k=1; k<= 3; ++k ) {
					xx1(k,atomno) = x1(k);
					xx2(k,atomno) = x2(k);
				}
			}
		}
		runtime_assert( atomno == natoms );
	}

	// calculate starting center of mass (COM):
	FArray1D_double COM(3);
	COMAS(xx1,wt,natoms,COM(1),COM(2),COM(3));

	// superimpose:: shifts xx1, shifts and transforms xx2;
	double rms;
	rmsfitca2(natoms,xx1,xx2,wt,nsup,rms);


	if ( true ) { // debug:
		double tmp1,tmp2,tmp3;
		COMAS(xx1,wt,natoms,tmp1,tmp2,tmp3); // store xcen,ycen,zcen vals for later
		//std::cout << "zero??: " << std::abs(tmp1) + std::abs(tmp2) + std::abs(tmp3)
		//					<< std::endl;
		runtime_assert( std::abs(tmp1) + std::abs(tmp2) + std::abs(tmp3) < 1e-3 );
	}

	{ // translate xx2 by COM and fill in the new mod_pose coordinates
		Size atomno(0);
		Vector x2;
		for ( Size i=1; i<= mod_pose.total_residue(); ++i ) {
			for ( Size j=1; j<= mod_pose.residue_type(i).natoms(); ++j ) { // use residue_type to prevent internal coord update
				++atomno;
				for ( Size k=1; k<= 3; ++k ) x2(k) = xx2(k,atomno) + COM(k);
				mod_pose.set_xyz( id::AtomID( j,i), x2 );
			}
		}
		runtime_assert( atomno == natoms );
	}


	return ( static_cast < Real > ( rms ) );
}

/// @details both poses must have the same length.
Real
calpha_superimpose_pose(
	pose::Pose & mod_pose,
	pose::Pose const & ref_pose
)
{
	runtime_assert( mod_pose.total_residue() == ref_pose.total_residue() );
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, mod_pose, id::BOGUS_ATOM_ID );
	for ( Size ii = 1; ii <= mod_pose.total_residue(); ++ii ) {
		if ( ii > ref_pose.total_residue() ) break;
		if ( ! mod_pose.residue(ii).has("CA") ) continue;
		if ( ! ref_pose.residue(ii).has("CA") ) continue;

		id::AtomID const id1( mod_pose.residue(ii).atom_index("CA"), ii );
		id::AtomID const id2( ref_pose.residue(ii).atom_index("CA"), ii );
		atom_map.set( id1, id2 );

	}
	return superimpose_pose( mod_pose, ref_pose, atom_map );
}

core::Real
CA_rmsd_symmetric(
  const core::pose::Pose & native_pose,
  const core::pose::Pose & pose
)
{
	using namespace core;
	using namespace conformation::symmetry;

	runtime_assert( core::pose::symmetry::is_symmetric( native_pose ) || core::pose::symmetry::is_symmetric( pose ) );

	SymmetricConformation const & symm_conf (
		dynamic_cast<SymmetricConformation const & > ( pose.conformation()) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );


		int const nres_monomer ( symm_info->num_independent_residues() );
    int const N ( symm_info->subunits() );
    int const nres ( symm_info->num_total_residues_without_pseudo() );
    FArray2D< core::Real > p1a_shuffle( 3, nres );

    core::Real rms = 1e3; //Since fast_rms has not been evaluated yet

  // copy coords into Real arrays
  int natoms;
  FArray2D< core::Real > p1a;//( 3, pose1.total_residue() );
  FArray2D< core::Real > p2a;//( 3, pose2.total_residue() );
  fill_rmsd_coordinates( natoms, p1a, p2a, native_pose, pose, is_protein_CA );
	if (natoms%nres_monomer != 0 ) {
		tr.Warning << "CA atoms in fill_rmsd " << natoms << "is not a multiple of number of residues per subunit " << nres_monomer << std::endl;
	}

	// Calc rms
	std::vector< std::vector<int> > shuffle_map;
	create_shuffle_map_recursive_rms(std::vector<int>(), N,shuffle_map);
	for (int j=1; j < int (shuffle_map.size()); j++ ){
		for (int i=0; i < N; ++i ) {
			int const begin ( shuffle_map.at(j).at(i)*nres_monomer*3);
			for ( int k = 0; k < nres_monomer*3; ++k ) {
				int const begin_shuffled (i*nres_monomer*3);
					p1a_shuffle[begin_shuffled+k] = p1a[begin+k];
			}
		}
		Real rms_shuffle = numeric::model_quality::rms_wrapper( natoms, p1a_shuffle, p2a );
		if ( rms_shuffle < rms ) {
			rms = rms_shuffle;
		}
	}

	if(rms < 0.00001) rms = 0.0;
  return rms;
}

// @details This is a recursive algorithm to generate all combinations of
// n digits where a number can only occur once in the sequence.
// The size scales as N! so don't use this for large values of N!!!
void
create_shuffle_map_recursive_rms(
	std::vector<int> sequence,
	int const N,
	std::vector< std::vector<int> > & map
)
{
	if ( int(sequence.size()) == N ){
		map.push_back(sequence);
		return;
	}
	for (int i=0; i< N; i++) {
		bool exist (false);
		for (int j=0; j < int(sequence.size()); j++) {
			if (sequence.at(j) == i )
				exist = true;
		}
		if (!exist) {
			std::vector<int> sequence_tmp (sequence);
			sequence_tmp.push_back(i);
			create_shuffle_map_recursive_rms(sequence_tmp,N,map);
		}
	}
}

/////////////////////////////////////////////////////////////
/// @details Should be more robust to crazy variant type mismatches. Both poses must have the same length.
Real
rms_at_corresponding_atoms(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	std::map< core::id::AtomID, core::id::AtomID > atom_id_map
)
{
	utility::vector1< Size > calc_rms_res;
	for ( Size n = 1; n <= mod_pose.total_residue(); n++ ) calc_rms_res.push_back( n );

	return rms_at_corresponding_atoms( mod_pose, ref_pose, atom_id_map, calc_rms_res );
}

/////////////////////////////////////////////////////////////
/// @details Should be more robust to crazy variant type mismatches. Both poses must have the same length.
Real
rms_at_corresponding_atoms(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	std::map< core::id::AtomID, core::id::AtomID > atom_id_map,
	utility::vector1< Size > const & calc_rms_res
)
{
	utility::vector1< bool > is_calc_rms( mod_pose.total_residue(), false );
	for ( Size n = 1; n <= calc_rms_res.size(); n++ ) is_calc_rms[ calc_rms_res[ n ] ] = true;

	utility::vector1< Vector > p1_coords, p2_coords;

	for ( std::map< core::id::AtomID, core::id::AtomID >::const_iterator iter = atom_id_map.begin();
				iter != atom_id_map.end(); iter++ ) {

		assert ( mod_pose.residue( (iter->first).rsd() ).atom_name(  (iter->first).atomno() ) ==
						 ref_pose.residue( (iter->second).rsd() ).atom_name(  (iter->second).atomno() ) );

		if ( !is_calc_rms[ (iter->first).rsd() ] ) continue;

		Vector const & p1(  mod_pose.xyz( iter->first ));
		Vector const & p2(  ref_pose.xyz( iter->second ));
		p1_coords.push_back(  p1 );
		p2_coords.push_back(  p2 );
	}
	return numeric::model_quality::calc_rms( p1_coords, p2_coords );
}

/// @details Calculates RMSD of all atoms in AtomID map, no need for the poses to be the same length.
Real
rms_at_all_corresponding_atoms(
        pose::Pose const & mod_pose,
        pose::Pose const & ref_pose,
        std::map< core::id::AtomID, core::id::AtomID > atom_id_map
)
{
	utility::vector1< Vector > p1_coords, p2_coords;

	for ( std::map< core::id::AtomID, core::id::AtomID >::const_iterator iter = atom_id_map.begin();
				iter != atom_id_map.end(); iter++ ) {

		runtime_assert ( mod_pose.residue( (iter->first).rsd() ).atom_name(  (iter->first).atomno() ) ==
										 ref_pose.residue( (iter->second).rsd() ).atom_name(  (iter->second).atomno() ) );

		Vector const & p1(  mod_pose.xyz( iter->first ));
		Vector const & p2(  ref_pose.xyz( iter->second ));
		p1_coords.push_back(  p1 );
		p2_coords.push_back(  p2 );
	}
	return numeric::model_quality::calc_rms( p1_coords, p2_coords );
}

Real
rms_at_corresponding_atoms_no_super(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	std::map< core::id::AtomID, core::id::AtomID > atom_id_map ){

	utility::vector1< Size > calc_rms_res;
	for ( Size n = 1; n <= mod_pose.total_residue(); n++ ) calc_rms_res.push_back( n );

	return rms_at_corresponding_atoms_no_super( mod_pose, ref_pose, atom_id_map, calc_rms_res );
}

Real
rms_at_corresponding_atoms_no_super(
	pose::Pose const & mod_pose,
	pose::Pose const & ref_pose,
	std::map< core::id::AtomID, core::id::AtomID > atom_id_map,
	utility::vector1< Size > const & calc_rms_res
)
{
	utility::vector1< bool > is_calc_rms( mod_pose.total_residue(), false );
	for ( Size n = 1; n <= calc_rms_res.size(); n++ ) is_calc_rms[ calc_rms_res[ n ] ] = true;

	Size natoms( 0 );
	Real sum( 0.0 );
	for ( std::map< core::id::AtomID, core::id::AtomID >::const_iterator iter = atom_id_map.begin();
				iter != atom_id_map.end(); iter++ ) {

		assert ( mod_pose.residue( (iter->first).rsd() ).atom_name(  (iter->first).atomno() ) ==
						 ref_pose.residue( (iter->second).rsd() ).atom_name(  (iter->second).atomno() ) );

		if ( !is_calc_rms[ (iter->first).rsd() ] ) continue;

		Vector const & p1(  mod_pose.xyz( iter->first ));
		Vector const & p2(  ref_pose.xyz( iter->second ));

		sum += (p1 - p2).length_squared();
		natoms++;
	}
	return std::sqrt( sum / natoms );
}

Real
rms_at_corresponding_heavy_atoms(
													 pose::Pose const & mod_pose,
													 pose::Pose const & ref_pose
													 )
{
 	std::map< core::id::AtomID, core::id::AtomID > atom_id_map;
 	setup_matching_heavy_atoms( mod_pose, ref_pose, atom_id_map );
 	return rms_at_corresponding_atoms( mod_pose, ref_pose, atom_id_map );
}

void
setup_matching_heavy_atoms( core::pose::Pose const & pose1, core::pose::Pose const & pose2, 	std::map< core::id::AtomID, core::id::AtomID > & atom_id_map ){

	using namespace core::id;
	using namespace core::conformation;

	atom_id_map.clear();
	assert( pose1.sequence() == pose2.sequence() );

	for (Size i = 1; i <= pose1.total_residue(); i++ ) {
		Residue const & rsd1 = pose1.residue( i );
		Residue const & rsd2 = pose2.residue( i );

		for ( Size j = 1; j <= rsd1.nheavyatoms(); j++ ) {
			std::string name( rsd1.atom_name( j )  );
			if ( !rsd2.has( name ) ) continue;

			if( rsd1.is_virtual(j)) continue;

			Size const j2( rsd2.atom_index( name ) );
			if( rsd2.is_virtual(j2) ) continue;

			atom_id_map[ AtomID( j, i ) ] = AtomID(  j2, i ) ;
		}
	}
}

/// @details Iterates over all non-hydrogen sidechain atoms of two residues and returns their rmsd without superposition.
core::Real
residue_sc_rmsd_no_super( core::conformation::ResidueCOP res1, core::conformation::ResidueCOP res2, bool const fxnal_group_only /*false*/ )
{
	runtime_assert( res1->name3() == res2->name3() );
	std::vector< core::Size > compare_atoms;
	core::Real sum2( 0.0 );

	if( fxnal_group_only ) {
		std::string const name1 = res1->name3();
		if( name1 == "GLY" )
			compare_atoms.push_back( res1->atom_index("CA"));
		if( name1 == "ALA" )
			compare_atoms.push_back( res1->atom_index("CB"));
		if( name1 == "SER" )
			compare_atoms.push_back( res1->atom_index("OG"));
		if( name1 == "THR" )
			compare_atoms.push_back( res1->atom_index("OG1"));
		if( name1 == "CYS" )
			compare_atoms.push_back( res1->atom_index("SG"));
		if( name1 == "VAL" ) {
			compare_atoms.push_back( res1->atom_index("CG1"));
			compare_atoms.push_back( res1->atom_index("CG2"));
		}
		if( name1 == "LEU" ) {
			compare_atoms.push_back( res1->atom_index("CD1"));
			compare_atoms.push_back( res1->atom_index("CD2"));
		}
		if( name1 == "ILE" ) {
			compare_atoms.push_back( res1->atom_index("CD1"));
			compare_atoms.push_back( res1->atom_index("CG2"));
		}
		if( name1 == "MET" ) {
			compare_atoms.push_back( res1->atom_index("CE"));
			compare_atoms.push_back( res1->atom_index("SD"));
		}
		if( name1 == "PRO" )
			compare_atoms.push_back( res1->atom_index("CG"));
		if( name1 == "PHE" ) {
			compare_atoms.push_back( res1->atom_index("CZ"));
			compare_atoms.push_back( res1->atom_index("CE1"));
			compare_atoms.push_back( res1->atom_index("CE2"));
		}
		if( name1 == "TYR" )
			compare_atoms.push_back( res1->atom_index("OH"));
		if( name1 == "TRP" )
			compare_atoms.push_back( res1->atom_index("NE1"));
		if( name1 == "ASP" ) {
			compare_atoms.push_back( res1->atom_index("OD1"));
			compare_atoms.push_back( res1->atom_index("OD2"));
		}
		if( name1 == "GLU" ) {
			compare_atoms.push_back( res1->atom_index("OE1"));
			compare_atoms.push_back( res1->atom_index("OE2"));
		}
		if( name1 == "ASN" ) {
			compare_atoms.push_back( res1->atom_index("OD1"));
			compare_atoms.push_back( res1->atom_index("ND2"));
		}
		if( name1 == "GLN" ) {
			compare_atoms.push_back( res1->atom_index("OE1"));
			compare_atoms.push_back( res1->atom_index("NE2"));
		}
		if( name1 == "HIS" ) {
			compare_atoms.push_back( res1->atom_index("ND1"));
			compare_atoms.push_back( res1->atom_index("NE2"));
		}
		if( name1 == "LYS" )
			compare_atoms.push_back( res1->atom_index("NZ"));
		if( name1 == "ARG" ) {
			compare_atoms.push_back( res1->atom_index("NH1"));
			compare_atoms.push_back( res1->atom_index("NH2"));
		}

	}
	else {
		core::Size num_atoms ( res1->natoms() );
		if ( num_atoms > res2->natoms() ){
			num_atoms = res2->natoms();
		}
		for ( core::Size j = res1->first_sidechain_atom(); j <= num_atoms; ++j ) {
			if( !res1->atom_type(j).is_heavyatom() ) continue;
			compare_atoms.push_back( j );
		}
	}

	// compare over sidechain heavy atoms only, ignoring hydrogens
	for ( std::vector< core::Size >::const_iterator it = compare_atoms.begin(); it != compare_atoms.end(); ++it ) {
		core::Size const atomno = *it;
		core::Vector const diff = res1->xyz(atomno) - res2->xyz(atomno);
		sum2 += diff.length_squared();
	}

	return std::sqrt(sum2 / compare_atoms.size() );
}

//////////////////////////////////////////////////////////////////////////
void
setup_matching_CA_atoms( core::pose::Pose const & pose1, core::pose::Pose const & pose2, 	std::map< core::id::AtomID, core::id::AtomID > & atom_id_map ){

	utility::vector1< std::string > protein_backbone_heavy_atoms;
	protein_backbone_heavy_atoms.push_back( " CA " );
	setup_matching_atoms_with_given_names( pose1, pose2, protein_backbone_heavy_atoms, atom_id_map );

}

//////////////////////////////////////////////////////////////////////////
void
setup_matching_protein_backbone_heavy_atoms( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
																						 std::map< core::id::AtomID, core::id::AtomID > & atom_id_map ){

	utility::vector1< std::string > protein_backbone_heavy_atoms;
	protein_backbone_heavy_atoms.push_back( " N  ");
	protein_backbone_heavy_atoms.push_back( " CA ");
	protein_backbone_heavy_atoms.push_back( " C  ");
	protein_backbone_heavy_atoms.push_back( " O  ");
	setup_matching_atoms_with_given_names( pose1, pose2, protein_backbone_heavy_atoms, atom_id_map );
}

//////////////////////////////////////////////////////////////////////////
void
setup_matching_atoms_with_given_names( core::pose::Pose const & pose1, core::pose::Pose const & pose2,
																			 utility::vector1< std::string > const & atom_names_to_find,
																			 std::map< core::id::AtomID, core::id::AtomID > & atom_id_map ){
	using namespace core::id;
	using namespace core::conformation;
	using namespace core::chemical;

	atom_id_map.clear();

	for (Size i = 1; i <= pose1.total_residue(); i++ ) {
		Residue const & rsd1 = pose1.residue( i );
		Residue const & rsd2 = pose2.residue( i );
		ResidueType const & rsd_type1 = pose1.residue_type( i );
		ResidueType const & rsd_type2 = pose2.residue_type( i );

		for ( Size n = 1; n <= atom_names_to_find.size(); n++ ) {
			std::string const & atom_name = atom_names_to_find[ n ];

			if ( !rsd_type1.has( atom_name )  ) continue;
			if ( !rsd_type2.has( atom_name )  ) continue;

			Size const j1 = rsd1.atom_index( atom_name );
			if( rsd1.is_virtual( j1 )) continue;

			Size const j2 = rsd2.atom_index( atom_name );
			if( rsd2.is_virtual( j2 )) continue;

			atom_id_map[ AtomID( j1, i ) ] = AtomID(  j2, i ) ;
		}
	}
}

/// @detail Computes the RMSD of the jump residues (jump point +/- 1 residue) of
/// <model> and <reference>. Jump residues are identified by scanning <reference>'s
/// FoldTree. Results are stored in the output parameter <rmsds>, keyed by the index
/// of the jump point. For example,
///
/// Jump 100 => 10
/// rmsds[10] = rmsd(residues 9-11 in reference, residues 9-11 in model)
void compute_jump_rmsd(const core::pose::Pose& reference,
                       const core::pose::Pose& model,
                       boost::unordered_map<core::Size, core::Real>* rmsds) {
  using core::Size;
  using core::kinematics::FoldTree;
  assert(rmsds);

  // Identify jump residues
  const FoldTree& tree = model.fold_tree();
  for (Size i = 1; i <= tree.nres(); ++i) {
    if (!tree.is_jump_point(i))
      continue;

    // Edge cases at pose start/stop
    if ((i - 1) < 1 || (i + 1) > model.total_residue())
      continue;

    (*rmsds)[i] = CA_rmsd(reference, model, i - 1, i + 1);
  }
}

} // namespace core
} // namespace scoring
