// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @details
///
///
///
/// @author olange: ported from original bblum-rosetta++ version $

// Unit Headers
#include <core/scoring/dssp/Dssp.hh>

// Package Headers

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/id/NamedAtomID.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>


// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

// #include <basic/options/option.hh>
// #include <basic/options/keys/OptionKeys.hh>

// numeric headers
#include <numeric/xyzVector.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>

#include <core/id/AtomID.hh>
#include <core/scoring/dssp/StrandPairing.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "core.scoring.dssp" );

using namespace core;
using namespace basic;
using namespace ObjexxFCL;
//using namespace basic::options;

namespace core {
namespace scoring {
namespace dssp {


Dssp::Dssp( core::pose::Pose const& pose ) {
	pair_set_ = NULL;
	compute( pose );
}

Dssp::~Dssp() {}

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Populates the hbond_bb_pair_score_ array with dssp energies
///
/// @details
/// Uses hydrogen bond energies computed a la dssp to fill the
/// hbond_bb_pair_score_ array.  Entry (i,j) is the backbone
/// hydrogen bond energy between residues i (acceptor) and j (donor).
///
/// @global_read
/// misc::current_pose::Eposition
/// misc::current_pose::full_coord
/// misc::Sizes::res
/// misc::current_pose::res_variant
/// aaproperties_pack::properties_per_aa_aav::HNpos
/// aaproperties_pack::properties_per_aa_aav::HNpos
/// param_aa::aa_pro
/// termini_ns::is_N_terminus
/// termini_ns::is_C_terminus
///
/// @global_write
/// dssp_ns::hbond_bb_pair_score_
///
/// @remarks
///
/// @references
///
/// @author olange: ported from original bblum-rosetta++ version $
///
//////////////////////////////////////////////////////////////////////////////
void
fill_hbond_bb_pair_score_dssp( pose::Pose const& pose, ObjexxFCL::FArray2D_float &hbond_bb_pair_score ) {
	Size const total_residue( pose.total_residue() );

	//initialize FArray
	hbond_bb_pair_score.dimension(total_residue, total_residue);
	hbond_bb_pair_score = 0.0f;

	// FArray3D_float full_coord_backup( 3, param::MAX_ATOM(), total_residue );
	// copy_the_relevant_atoms( full_coord_backup, full_coord );

	// copy_position_to_fullcoord(Eposition,full_coord,total_residue);
	// put_nh(full_coord, 1, total_residue, res, res_variant, is_N_terminus, is_C_terminus, 0.0f);
	// initialize_fullcoord_array(Eposition, full_coord, total_residue, res, res_variant);

	// float C[3], O[3], N[3], H[3], rON, rCH, rOH, rCN, E, dist, total;
	for ( Size i = 1; i <= total_residue; ++i ) {
		if ( pose.residue(i).aa() == core::chemical::aa_vrt ) continue;
		for ( Size j = 1; j <= total_residue; ++j ) {
			if ( !pose.residue(i).is_protein() ) continue;
			if ( !pose.residue(j).is_protein() ) continue;
			if ( i == j || i - j == 1 || j - i == 1 ) continue;
			//chu skip non-protein residues
			if ( !(pose.residue(i).is_protein() && pose.residue(j).is_protein()) ) continue;

			if ( pose.residue(j).aa() == core::chemical::aa_vrt ) continue;

			//ignore if CA-CA > 10A
			id::NamedAtomID CA1("CA",i );
			id::NamedAtomID CA2("CA",j );
			if ( pose.xyz( CA1 ).distance( pose.xyz( CA2 ) ) > 10.0 ) continue;

			//get C at i
			PointPosition pC;
			pC = pose.xyz( id::NamedAtomID( "C", i ) );

			//get O at i
			PointPosition pO;
			pO = pose.xyz( id::NamedAtomID( "O", i ) );

			//get N at j
			PointPosition pN;
			pN = pose.xyz( id::NamedAtomID( "N", j ) );

			Real rON, rCN;
			rON = pO.distance(  pN );
			rCN = pC.distance(  pN );

			//   Size HN1_pos = HNpos( res(j), res_variant(j) );

			chemical::ResidueType const& rtj ( pose.residue_type ( j ) );
			runtime_assert( rtj.has( "N" ) );
			for ( Size hatom = rtj.attached_H_begin( rtj.atom_index( "N" ) ), ehatom = rtj.attached_H_end( rtj.atom_index( "N" ) );
					hatom <= ehatom; ++hatom ) {

				Real rOH, rCH;
				PointPosition pH;
				pH = pose.xyz( id::AtomID( hatom, j ) );
				rOH = pH.distance( pO );
				rCH = pH.distance( pC );
				Real E = 27.888f * ( 1.0f / rON + 1.0f / rCH - 1.0f / rOH - 1.0f / rCN);
				//tr.Trace << "rOH " << rOH << " rCH " << rCH <<  " rCN " << rCN << " rON " << rON << " rNH " << distance( pN, pH ) << std::endl;
				//tr.Trace << "hbond energy for H " << hatom << " respair " << i << " " << j << " " << E << std::endl;
				if ( E < hbond_bb_pair_score(i, j) ) {
					hbond_bb_pair_score(i, j) = E;
					break;
				}
			}
		}
	}
	//copy_the_relevant_atoms( full_coord, full_coord_backup);
}

void
Dssp::dssp( FArray1_char &secstruct ) {
	for ( Size i = 1; i <= dssp_secstruct_.size(); i++ ) {
		secstruct(i) = dssp_secstruct_(i);
	}
}

void
Dssp::dssp_featurizer( FArray1_char &secstruct ) {
	for ( Size i = 1; i <= dssp_secstruct_.size(); i++ ) {
		if ( dssp_secstruct_(i) == 'H' || dssp_secstruct_(i) == 'G' || dssp_secstruct_(i) == 'I' ) {
			secstruct(i) = 'H';
		} else if ( dssp_secstruct_(i) == 'B' || dssp_secstruct_(i) == 'E' ) {
			// see if paired on one or both sides
			secstruct(i) = pair_set_->featurizer_state(i);
		} else secstruct(i) = 'L';
	}
}

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Reduces to E/H/L secondary structure alphabet
///
/// @details
/// This function simply reduces dssp's secondary structure alphabet (which
/// includes 3- and 5-turn helices and various kinds of loop, and
/// differentiates lone beta bridges from extended beta strand pairings)
/// Sizeo the standard E/H/L alphabet, as follows:
/// G,H,I --> H
/// E,B --> E
/// S,T,blank --> L
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author olange: ported from original bblum-rosetta++ version $
///
//////////////////////////////////////////////////////////////////////////////
void
Dssp::dssp_reduced( FArray1_char &secstruct ) {
	for ( Size i = 1; i <= dssp_secstruct_.size(); i++ ) {
		if ( dssp_secstruct_(i) == 'H' || dssp_secstruct_(i) == 'G' || dssp_secstruct_(i) == 'I' ) {
			secstruct(i) = 'H';
		} else if ( dssp_secstruct_(i) == 'B' || dssp_secstruct_(i) == 'E' ) {
			secstruct(i) = 'E';
		} else secstruct(i) = 'L';
	}
}

void
Dssp::dssp_reduced_IG_as_L_if_adjcent_H( FArray1_char &secstruct ) {
	for ( Size i = 1; i <= dssp_secstruct_.size(); i++ ) {
		if ( dssp_secstruct_(i) == 'H' || dssp_secstruct_(i) == 'G' || dssp_secstruct_(i) == 'I' ) {
			if ( ( dssp_secstruct_(i)=='G'|| dssp_secstruct_(i)=='I' ) &&
					( (i>1 && dssp_secstruct_(i-1)=='H') || (i<dssp_secstruct_.size() && dssp_secstruct_(i+1)=='H'))
					) {
				secstruct(i) = 'L';
			} else {
				secstruct(i) = 'H';
			}
		} else if ( dssp_secstruct_(i) == 'B' || dssp_secstruct_(i) == 'E' ) {
			secstruct(i) = 'E';
		} else secstruct(i) = 'L';
	}
}

void
Dssp::dssp_reduced_IG_as_L( FArray1_char &secstruct ) {
	for ( Size i = 1; i <= dssp_secstruct_.size(); i++ ) {
		if ( dssp_secstruct_(i) == 'H' ) {
			secstruct(i) = 'H';
		} else if ( dssp_secstruct_(i) == 'B' || dssp_secstruct_(i) == 'E' ) {
			secstruct(i) = 'E';
		} else secstruct(i) = 'L';
	}
}

void
Dssp::dssp_reduced() {
	dssp_reduced( dssp_secstruct_ );
}

//////////////////////////////////////////////////////////////////////////////
///
/// @brief Runs dssp, calculating per-residue secondary structure.
///
/// @details
/// dssp is a standard algorithm for per-residue secondary structure analysis.
/// It has the following alphabet:
/// H: 4-turn helix
/// B: beta bridge
/// E: extended beta strand
/// G: 3-turn helix
/// I: 5-turn helix
/// T: helix-like loop (some nearby hbonds)
/// S: loop with high curvature
///  : loop
/// Most of these determinations are made on the basis of hydrogen bonds,
/// with torsion angles disregarded.  The designations are reported with priority
/// given to categories higher in the list (e.g. if a residue has both high-
/// curvature (S) and helical hydrogen bonds (H), it will be reported as H).
/// Experiments indicate that this function agrees almost completely with
/// the true dssp results (run on dumped pdbs).  Slight differences can perhaps
/// be chalked up to unequal placement of hydrogens.
///
/// @global_read
/// dssp_ns::hbond_bb_pair_score_
///
/// @global_write
/// dssp_ns::hbond_bb_pair_score_
///
/// @remarks
///
/// @references
///
/// @author olange: ported from original bblum-rosetta++ version $
///
//////////////////////////////////////////////////////////////////////////////
void
Dssp::compute( pose::Pose const& pose ) {

	float dssp_hbond_threshold = -0.5;

	core::Size total_residue( pose.total_residue() );
	ObjexxFCL::FArray1D_bool invalid; //Should we omit this residue when doing DSSP?
	dssp_secstruct_.dimension(total_residue);
	invalid.dimension(total_residue);

	fill_hbond_bb_pair_score_dssp( pose, hbond_bb_pair_score_ ); // fills hbond_bb_pair_score_ array

	// Initialize to all loops
	for ( Size i = 1; i <= total_residue; i++ ) {
		dssp_secstruct_(i) = ' ';
		invalid(i) = (! pose.residue(i).is_protein()) || ( pose.residue(i).aa() == core::chemical::aa_vrt);
	}

	bool helix;

	// Record all 5-turn helices (I)
	if ( total_residue > 5 ) {
		for ( Size i=1; i <= total_residue-5; i++ ) {
			if ( invalid(i) || invalid(i+1) || invalid(i+2) || invalid(i+3) || invalid(i+4) || invalid(i+5) ) continue;
			if ( hbond_bb_pair_score_(i, i + 5) < dssp_hbond_threshold ) {
				helix = i < total_residue - 5 &&
					hbond_bb_pair_score_(i + 1, i + 6) < dssp_hbond_threshold;
				for ( Size j = 1; j < 6; j++ ) {
					if ( helix ) {
						dssp_secstruct_(i + j) = 'I';
					} else if ( j < 5 && dssp_secstruct_(i + j) == ' ' ) {
						dssp_secstruct_(i + j) = 'T';
					}
				}
			}
		}
	}
	// Record all 3-turn helices (G)
	if ( total_residue > 3 ) {
		for ( Size i = 1; i <= total_residue - 3; i++ ) {
			if ( invalid(i) || invalid(i+1) || invalid(i+2) || invalid(i+3) ) continue;
			if ( hbond_bb_pair_score_(i, i + 3) < dssp_hbond_threshold ) {
				helix = i < total_residue - 3 &&
					hbond_bb_pair_score_(i + 1, i + 4) < dssp_hbond_threshold;
				for ( Size j = 1; j < 4; j++ ) {
					if ( helix ) {
						dssp_secstruct_(i+j) = 'G';
					} else if ( j < 3 && dssp_secstruct_(i+j) == ' ' ) {
						dssp_secstruct_(i+j) = 'T';
					}
				}
			}
		}
	}

	// Record all strands (B and E)
	pair_set_ = StrandPairingSetOP( new StrandPairingSet( hbond_bb_pair_score_, dssp_hbond_threshold, pose ) );

	for ( Size i = 1; i <= total_residue; i++ ) {
		char state = pair_set_->dssp_state(i);
		if ( state != ' ' ) dssp_secstruct_(i) = state;
	}

	// Record all 4-turn helices (H)
	if ( total_residue > 4 ) {
		for ( Size i = 1; i <= total_residue - 4; i++ ) {
			if ( invalid(i) || invalid(i+1) || invalid(i+2) || invalid(i+3) || invalid(i+4) ) continue;
			if ( hbond_bb_pair_score_(i, i + 4) < dssp_hbond_threshold ) {
				helix = i < total_residue - 4 &&
					hbond_bb_pair_score_(i + 1, i + 5) < dssp_hbond_threshold;
				for ( Size j = 1; j < 5; j++ ) {
					if ( helix ) {
						/*
						&& torsion_bin(phi(i+j), psi(i+j), omega(i+j)) != 'A')
						// only allow helix chunks in which every residue has bin A
						*/
						dssp_secstruct_(i + j) = 'H';
					} else if ( j < 4 && dssp_secstruct_(i + j) == ' ' ) {
						dssp_secstruct_(i + j) = 'T';
					}
				}
			}
		}
	}

	// Record all tight turns (S), only if still a loop
	if ( total_residue > 2 ) {
		for ( Size i = 3; i <= total_residue - 2; i++ ) {
			if ( invalid(i-2) || invalid(i) || invalid(i+2) ) continue;
			if ( dssp_secstruct_(i) == ' ' ) {
				Vector const v1 ( pose.xyz( id::NamedAtomID("CA",i ) ) - pose.xyz(  id::NamedAtomID("CA",i-2 ) ) );
				Vector const v2 ( pose.xyz( id::NamedAtomID("CA",i+2 ) ) - pose.xyz(  id::NamedAtomID("CA",i ) ) );
				Real dot = angle_of( v1, v2 );
				if ( dot < .34202014 ) {
					dssp_secstruct_(i) = 'S';
				}
			}
		}
	}
}

bool
Dssp::paired(Size res1, Size res2, bool antiparallel) {
	return pair_set_->paired(res1, res2, antiparallel);
}


void
Dssp::insert_ss_into_pose( core::pose::Pose & pose, bool recompute ) {
	if ( recompute ) compute( pose );
	dssp_reduced( dssp_secstruct_ );
	for ( core::Size i = 1; i <= dssp_secstruct_.size(); ++i ) {
		pose.set_secstruct( i, dssp_secstruct_(i) );
	}
}

void
Dssp::insert_ss_into_pose_no_IG_helix( core::pose::Pose & pose, bool recompute ) {
	if ( recompute ) compute( pose );
	dssp_reduced_IG_as_L( dssp_secstruct_ );
	for ( core::Size i = 1; i <= dssp_secstruct_.size(); ++i ) {
		pose.set_secstruct( i, dssp_secstruct_(i) );
	}
}

void
Dssp::insert_dssp_ss_into_pose( core::pose::Pose & pose, bool recompute ) {
	if ( recompute ) compute( pose );
	for ( core::Size i = 1; i <= dssp_secstruct_.size(); ++i ) {
		char ss = dssp_secstruct_(i);
		if ( ss==' ' ) ss = '_';
		pose.set_secstruct( i, ss );
	}
}
void
Dssp::insert_edge_ss_into_pose( core::pose::Pose & pose, bool recompute ) {
	if ( recompute ) compute( pose );
	for ( core::Size i = 1; i <= dssp_secstruct_.size(); ++i ) {
		char ss = dssp_secstruct_(i);
		if ( ss==' ' ) ss = 'L';
		if ( ss=='E' && num_pairings(i)<2 ) ss = 'U';
		pose.set_secstruct( i, ss );
	}
}


char
Dssp::get_dssp_secstruct( core::Size resid ) {
	return dssp_secstruct_( resid );
}

std::string
Dssp::get_dssp_secstruct() {
	dssp_reduced( dssp_secstruct_ );
	std::string sequence;
	for ( core::Size i = 1; i <= dssp_secstruct_.size(); ++i ) {
		sequence += dssp_secstruct_( i );
	}
	return sequence;
}

std::string
Dssp::get_dssp_reduced_IG_as_L_secstruct() {
	dssp_reduced_IG_as_L( dssp_secstruct_ );
	std::string sequence;
	for ( core::Size i = 1; i <= dssp_secstruct_.size(); ++i ) {
		sequence += dssp_secstruct_( i );
	}
	return sequence;
}

float
Dssp::bb_pair_score( Size res1, Size res2 )
{
	return hbond_bb_pair_score_( res1, res2 );
}

Size
Dssp::num_pairings(Size resi ) const {
	Size npair = 0;
	for ( StrandPairingSet::const_iterator i = pair_set_->begin(); i != pair_set_->end(); ++i ) {
		if ( i->contains(resi) ) ++npair;
	}
	return npair;
}
bool
Dssp::in_paired_strands(Size res1, Size res2 ) const {
	for ( StrandPairingSet::const_iterator i = pair_set_->begin(); i != pair_set_->end(); ++i ) {
		if ( i->contains(res1) && i->contains(res2) ) return true;
	}
	return false;
}


} //dssp
} //scoring
} //core
