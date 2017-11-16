// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/flxbb/utility.cc
/// @brief
/// @author Nobuyasu Koga (nobuyasua@uw.edu)

// Unit headers
#include <protocols/flxbb/utility.hh>

// Package headers
#include <protocols/parser/BluePrint.hh>

// Project headers
#include <basic/options/option.hh>
#include <basic/options/keys/flxbb.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/SS_Info.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

/// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

static basic::Tracer TR( "protocols.flxbb.FlxbbDesign.utility" );

namespace protocols {
namespace flxbb {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef core::Size Size;
typedef core::Real Real;
typedef std::string String;
typedef core::pose::Pose Pose;
typedef protocols::parser::BluePrintOP BluePrintOP;
typedef core::scoring::constraints::ConstraintOPs ConstraintOPs;

ConstraintOPs
constraints_sheet( Pose const & pose, BluePrintOP const & blueprint, Real const coef, Real const condist )
{
	using core::scoring::constraints::AtomPairConstraint;
	using core::scoring::constraints::BoundFunc;
	using core::scoring::constraints::ConstraintOPs;
	using core::scoring::func::ScalarWeightedFunc;
	using core::scoring::func::ScalarWeightedFuncOP;
	using protocols::fldsgn::topology::StrandPairing;
	using protocols::fldsgn::topology::StrandPairings;
	using protocols::fldsgn::topology::StrandPairingOP;
	using protocols::fldsgn::topology::StrandPairingSet;
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;

	// returned data
	ConstraintOPs csts;

	// set constraint func
	Real lb( 0.0 );
	Real ub( condist );
	Real sd( 1.0 );
	String tag( "constraints_in_beta_sheet" );
	ScalarWeightedFuncOP cstfunc( new ScalarWeightedFunc( coef, core::scoring::func::FuncOP( new BoundFunc( lb, ub, sd, tag ) ) ) );

	//flo sep '12 add more accurate constraints by also constraining proper angles along paired residues
	core::scoring::func::FuncOP cacb_dihedral_func( new core::scoring::constraints::OffsetPeriodicBoundFunc(-0.9,0.9, sqrt(1.0/42.0), "dihed_cacb", 6.28, 0.0 ) );
	core::scoring::func::FuncOP bb_dihedral_func( new core::scoring::constraints::OffsetPeriodicBoundFunc(-0.52,0.52, sqrt(1.0/42.0), "dihed_bb", 3.14, 0.0 ) );
	core::scoring::func::FuncOP bb_angle_func( new BoundFunc(1.22,1.92, sqrt(1.0/42.0), "angle_bb") );

	// set constraints to csts
	Size nres( pose.size() );
	//flo sep '12 in case we have ligands in the pose, don't count them
	for ( core::Size i = nres; i != 0; i-- ) {
		if ( pose.residue_type(i).is_ligand() ) nres--;
		else break;
	}
	runtime_assert( nres == blueprint->total_residue() );

	TR << "Blueprint file is used for determining constrained residue pairs.  " << std::endl;
	TR << "Constrains between CA-CA atoms in sheet are applied for the following residues. " << std::endl;
	TR << "dist=" << condist << ", coef=" << coef << std::endl;

	SS_Info2_OP ssinfo( new SS_Info2( pose, blueprint->secstruct() ) );
	StrandPairingSet spairset( blueprint->strand_pairings(), ssinfo );
	StrandPairings spairs = spairset.strand_pairings();
	for ( utility::vector1< StrandPairingOP >::const_iterator it=spairs.begin(); it!=spairs.end(); ++it ) {

		StrandPairing & spair=**it;
		for ( Size iaa=spair.begin1(); iaa<=spair.end1(); iaa++ ) {
			Size jaa( spair.residue_pair( iaa ) );
			TR << iaa << ' ' << jaa << std::endl;
			core::id::AtomID atom1( pose.residue_type( iaa ).atom_index( "CA" ), iaa );
			core::id::AtomID atom2( pose.residue_type( jaa ).atom_index( "CA" ), jaa );
			csts.push_back( core::scoring::constraints::ConstraintOP( new AtomPairConstraint( atom1, atom2, cstfunc ) ) );
			//flo sep '12: constrain dihedral, might be more accurate
			if ( basic::options::option[ basic::options::OptionKeys::flxbb::constraints_sheet_include_cacb_pseudotorsion ].value() ) {
				core::id::AtomID resi_n( pose.residue_type( iaa ).atom_index( "N" ), iaa );
				core::id::AtomID resi_c( pose.residue_type( iaa ).atom_index( "C" ), iaa );
				core::id::AtomID resi_o( pose.residue_type( iaa ).atom_index( "O" ), iaa );
				core::id::AtomID resj_n( pose.residue_type( jaa ).atom_index( "N" ), jaa );
				core::id::AtomID resj_c( pose.residue_type( jaa ).atom_index( "C" ), jaa );
				core::id::AtomID resj_o( pose.residue_type( jaa ).atom_index( "O" ), jaa );
				csts.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::DihedralConstraint( resi_o, resi_n, resi_c, resj_c, bb_dihedral_func ) ) );
				csts.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::DihedralConstraint( resj_o, resj_n, resj_c, resi_c, bb_dihedral_func ) ) );
				if ( spair.orient() == 'P' ) {
					csts.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint( resi_n, resi_c, resj_c, bb_angle_func ) ) );
					csts.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint( resj_n, resj_c, resi_c, bb_angle_func ) ) );
				} else if ( spair.orient() == 'A' ) {
					csts.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint( resi_n, resi_c, resj_n, bb_angle_func ) ) );
					csts.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AngleConstraint( resj_n, resj_c, resi_n, bb_angle_func ) ) );
				}
				if ( (pose.residue_type( iaa ).name3() == "GLY") || (pose.residue_type( jaa ).name3() == "GLY" ) ) continue; // don't bother restraining cacb dihedral with gly
				core::id::AtomID resi_cb( pose.residue_type( iaa ).atom_index( "CB" ), iaa );
				core::id::AtomID resj_cb( pose.residue_type( jaa ).atom_index( "CB" ), jaa );
				csts.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::DihedralConstraint( resi_cb, atom1, atom2, resj_cb, cacb_dihedral_func ) ) );
			}
			// flo sep '12 over
		} // for( Size i=1 )

	} // StrandPairingOP

	return csts;
} // constraint_sheet


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ConstraintOPs
constraints_NtoC( Pose const & pose, Real const coef, Real const condist )
{
	using core::scoring::constraints::AtomPairConstraint;
	using core::scoring::constraints::BoundFunc;
	using core::scoring::constraints::ConstraintOPs;
	using core::scoring::func::ScalarWeightedFunc;
	using core::scoring::func::ScalarWeightedFuncOP;

	// returned data
	ConstraintOPs csts;

	// set constraint func
	Real lb( 0.0 );
	Real ub( condist );
	Real sd( 1.0 );
	String tag( "constraint_between_N_&_C_terminal_Calpha" );
	ScalarWeightedFuncOP cstfunc( new ScalarWeightedFunc( coef, core::scoring::func::FuncOP( new BoundFunc( lb, ub, sd, tag ) ) ) );

	Size nres( pose.size() );
	core::id::AtomID atom1( pose.residue_type( 1 ).atom_index( "CA" ), 1 );
	core::id::AtomID atom2( pose.residue_type( nres ).atom_index( "CA" ), nres );
	csts.push_back( core::scoring::constraints::ConstraintOP( new AtomPairConstraint( atom1, atom2, cstfunc ) ) );

	TR << "Constraints between N- and C- terminal: 1-" << nres << ", dist=" << condist << ", coef=" << coef << std::endl;

	return csts;
} // constraints_NtoC

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
ConstraintOPs
constraints_sheet( Pose const & pose, Real const coef, Real const condist )
{
	using core::scoring::constraints::AtomPairConstraint;
	using core::scoring::constraints::BoundFunc;
	using core::scoring::constraints::ConstraintOPs;
	using core::scoring::func::ScalarWeightedFunc;
	using core::scoring::func::ScalarWeightedFuncOP;
	using core::scoring::dssp::Dssp;
	using core::scoring::Strands;

	// returned data
	ConstraintOPs csts;

	// set constraint func
	Real lb( 0.0 );
	Real ub( condist );
	Real sd( 1.0 );
	std::string tag( "constraints_in_beta_sheet" );
	ScalarWeightedFuncOP cstfunc( new ScalarWeightedFunc( coef, core::scoring::func::FuncOP( new BoundFunc( lb, ub, sd, tag ) ) ) );

	// set secondary structure
	Dssp dssp( pose );

	// set strands
	Size nres( pose.size() );
	bool flag( false );
	Size istrand ( 0 );
	Strands strands( nres );

	for ( Size i=1; i<=nres; ++i ) {
		char ss( dssp.get_dssp_secstruct( i ) );
		if ( ss =='E' && flag == false ) {
			istrand ++;
			strands.SS_strand_end( 1, istrand ) = i;
			flag = true;
		}
		if ( ss !='E' && flag == true ) {
			strands.SS_strand_end( 2, istrand ) = i - 1;
			flag = false;
		}
	}
	strands.total_strands = istrand;

	TR << "# Constrains between CA-CA atoms in sheet are applied for the following residues " << std::endl;
	TR << "dist=" << condist << ", coef=" << coef << std::endl;

	Real condist2 = condist*condist;
	for ( int i=1; i<=strands.total_strands-1; ++i ) {
		for ( int j=i+1; j<=strands.total_strands; ++j ) {

			for ( int iresid=strands.SS_strand_end( 1, i ); iresid<=strands.SS_strand_end( 2, i ); iresid++ ) {
				for ( int jresid=strands.SS_strand_end( 1, j ); jresid<=strands.SS_strand_end( 2, j ); jresid++ ) {

					core::conformation::Residue const & ires( pose.residue( iresid ) );
					core::conformation::Residue const & jres( pose.residue( jresid ) );
					Size ica = ires.atom_index( "CA" );
					Size jca = jres.atom_index( "CA" );
					Real const dsq( ires.xyz( ica ).distance_squared( jres.xyz( jca ) ));
					if ( dsq<=condist2 ) {
						TR << iresid << ' ' << jresid << std::endl;
						core::id::AtomID atom1( pose.residue_type( iresid ).atom_index( "CA" ), iresid );
						core::id::AtomID atom2( pose.residue_type( jresid ).atom_index( "CA" ), jresid );
						csts.push_back( core::scoring::constraints::ConstraintOP( new AtomPairConstraint( atom1, atom2, cstfunc ) ) );
					}

				}// jresid
			} // iresid

		} // j
	} // i

	return csts;

} // constraint_sheet

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Looks for unknown amino acids in the pose and returns their indices
utility::vector1<Size>
find_ligands( Pose const & pose )
{
	utility::vector1<Size> retval;
	// look at each amino acid to see if it is of unknown type
	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( ! pose.residue( i ).is_protein() ) {
			TR << "Residue " << i << " (type=" << pose.residue(i).name3() << ") is probably a ligand" << std::endl;
			retval.push_back( i );
		}
	}
	// WARNING: This doesn't rearrange residue numbers to put the ligands at the end
	// However, this behavior IS important for many functions
	return retval;
}

} //namespace flxbb
} //namespace protocols

