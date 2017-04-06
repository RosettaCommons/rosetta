// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file util.cc
/// @brief util functions for general rna use
/// @details See also core/pose/rna/util.hh, core/chemical/rna/util.hh
/// @author Rhiju Das

#include <protocols/rna/util.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/rna/BasePair.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <basic/Tracer.hh>

using namespace core;
using namespace core::chemical::rna;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

static THREAD_LOCAL basic::Tracer TR( "protocols.rna.util" );

namespace protocols {
namespace rna {

/////////////////////////////////////////////////////////
/// @details
///
///   Adds 2-3 constraints between donor hydrogen and acceptor
///    across C-G, A-U, or G-U pairs, *and* asks that C1'-C1'
///    distance be near 10.5 Angstroms (prevents stacking!)
///
///   "suppress_factor" reduces strength, in use by FARNA/RNA_DeNovo
///
///   use_harmonic is standard 'spring' like function
///    that can pull in bases from very far away to optimal
///    distance of 1.9 Angstroms:
///
///     H-bond   ( x - 1.9 A )/( 0.25 / scale_factor )^2
///
///     C1'-C1'  ( x - 10.5 A )/( 1/0 / scale_factor )^2
///
///   If use_flat_harmonic is set to true, functions are zero out to
///    tolerance of 1.0 A and 2.0 A, and then go up quadratically.
///    [used in RECCES/free-energy estimation where we want constraint = 0
///     when base pairs are formed rather than some small number].
///
void
setup_base_pair_constraints(
	pose::Pose & pose,
	utility::vector1< std::pair< Size, Size > > const &  pairings,
	Real const scale_factor /* = 1.0 */,
	bool const use_flat_harmonic /* = false */ )
{

	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	Distance const WC_distance( 1.9 );
	Distance const distance_stddev( 0.25 / scale_factor ); //Hmm. Maybe try linear instead?
	FuncOP distance_func( new HarmonicFunc( WC_distance, distance_stddev ) );

	// Need to force base pairing -- not base stacking!
	Distance const C1prime_distance( 10.5 );
	Distance const C1prime_distance_stddev( 1.0 / scale_factor ); //Hmm. Maybe try linear instead?
	FuncOP C1prime_distance_func( new HarmonicFunc( C1prime_distance, C1prime_distance_stddev ) );

	if ( use_flat_harmonic ) {
		Distance const distance_tolerance( 1.0 );
		distance_func = FuncOP( new FlatHarmonicFunc( WC_distance, distance_stddev, distance_tolerance ) );
		Distance const C1prime_distance_tolerance( 2.0 );
		C1prime_distance_func = FuncOP( new FlatHarmonicFunc( C1prime_distance, C1prime_distance_stddev, C1prime_distance_tolerance  ) );
	}

	for ( auto const & pairing : pairings ) {

		Size const i = pairing.first;
		Size const j = pairing.second;

		if ( !pose.residue_type(i).is_RNA() ) continue;
		if ( !pose.residue_type(j).is_RNA() ) continue;

		if ( !pose.residue_type(i).is_coarse() ) { //fullatom

			Size const atom1 = pose.residue_type(i).RNA_info().c1prime_atom_index();
			Size const atom2 = pose.residue_type(j).RNA_info().c1prime_atom_index();
			pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(
				id::AtomID(atom1,i),
				id::AtomID(atom2,j),
				C1prime_distance_func, core::scoring::base_pair_constraint ) ) ) );

			utility::vector1< std::string > atom_ids1, atom_ids2;
			get_watson_crick_base_pair_atoms( pose.residue_type(i), pose.residue_type(j), atom_ids1, atom_ids2 );

			for ( Size p = 1; p <= atom_ids1.size(); p++ ) {

				Size const atom1 = pose.residue_type(i).atom_index( atom_ids1[p] ) ;
				Size const atom2 = pose.residue_type(j).atom_index( atom_ids2[p] ) ;

				TR << "BASEPAIR: Adding rna_force_atom_pair constraint: " << pose.residue_type(i).name1() << I(3,i) << " <-->  " <<
					pose.residue_type(j).name1() << I(3,j) << "   " <<
					atom_ids1[p] << " <--> " <<
					atom_ids2[p] << ".  [ " << atom1 << "-" << atom2 << "]" << std::endl;

				pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(
					id::AtomID(atom1,i),
					id::AtomID(atom2,j),
					distance_func, core::scoring::base_pair_constraint ) ) ) );
			}

		} else { //coarse-grained RNA

			static Real const coarse_WC_SUG_distance_min( 12.3 );
			static Real const coarse_WC_SUG_distance_max( 14.3 );
			static Real const coarse_WC_SUG_distance_stddev( 2.0 );
			static Real const coarse_WC_SUG_bonus( -2.5 );
			static core::scoring::func::FuncOP const coarse_SUG_distance_func( new core::scoring::func::FadeFunc( coarse_WC_SUG_distance_min - coarse_WC_SUG_distance_stddev,
				coarse_WC_SUG_distance_max + coarse_WC_SUG_distance_stddev,
				coarse_WC_SUG_distance_stddev,
				coarse_WC_SUG_bonus ) );

			{
				Size const atom1 = pose.residue_type(i).atom_index( " S  " );
				Size const atom2 = pose.residue_type(j).atom_index( " S  " );
				TR << "BASEPAIR: Adding rna_force_atom_pair constraint: " << pose.residue_type(i).name1() << I(3,i) << " <-->  " <<
					pose.residue_type(j).name1() << I(3,j) << "   " <<
					" S  " << " <--> " <<
					" S  " << ".  [ " << atom1 << "-" << atom2 << "]" << std::endl;
				pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(
					id::AtomID(atom1,i),
					id::AtomID(atom2,j),
					coarse_SUG_distance_func, core::scoring::base_pair_constraint ) ) ) );
			}

			static Real const coarse_WC_CEN_distance( 5.5 );
			static Real const coarse_WC_CEN_distance_stddev( 3.0 );
			static Real const coarse_WC_CEN_bonus( -5.0 );
			static core::scoring::func::FuncOP const coarse_CEN_distance_func( new core::scoring::func::FadeFunc( coarse_WC_CEN_distance - coarse_WC_CEN_distance_stddev,
				coarse_WC_CEN_distance + coarse_WC_CEN_distance_stddev,
				coarse_WC_CEN_distance_stddev,
				coarse_WC_CEN_bonus ) );

			{
				Size const atom1 = pose.residue_type(i).atom_index( " CEN" );
				Size const atom2 = pose.residue_type(j).atom_index( " CEN" );
				TR << "BASEPAIR: Adding rna_force_atom_pair constraint: " << pose.residue_type(i).name1() << I(3,i) << " <-->  " <<
					pose.residue_type(j).name1() << I(3,j) << "   " <<
					" CEN" << " <--> " <<
					" CEN" << ".  [ " << atom1 << "-" << atom2 << "]" << std::endl;
				pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(
					id::AtomID(atom1,i),
					id::AtomID(atom2,j),
					coarse_CEN_distance_func, core::scoring::base_pair_constraint ) ) ) );
			}

			static Real const coarse_WC_X_distance( 3.5 );
			static Real const coarse_WC_X_distance_stddev( 2.0 );
			static Real const coarse_WC_X_bonus( -5.0 );
			static core::scoring::func::FuncOP const coarse_X_distance_func( new core::scoring::func::FadeFunc( coarse_WC_X_distance - coarse_WC_X_distance_stddev,
				coarse_WC_X_distance + coarse_WC_X_distance_stddev,
				coarse_WC_X_distance_stddev, coarse_WC_X_bonus ) );

			{
				Size const atom1 = pose.residue_type(i).atom_index( " X  " );
				Size const atom2 = pose.residue_type(j).atom_index( " X  " );
				TR << "BASEPAIR: Adding rna_force_atom_pair constraint: " << pose.residue_type(i).name1() << I(3,i) << " <-->  " <<
					pose.residue_type(j).name1() << I(3,j) << "   " <<
					" X  " << " <--> " <<
					" X  " << ".  [ " << atom1 << "-" << atom2 << "]" << std::endl;
				pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(
					id::AtomID(atom1,i),
					id::AtomID(atom2,j),
					coarse_X_distance_func, core::scoring::base_pair_constraint ) ) ) );
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
// @brief score pose & get base pairing list.
// @details copies code from get_base_pairing_info [?]. Should be in core/pose/rna/util? Shouldn't have to score?
void
get_base_pairing_list(
	pose::Pose & pose,
	utility::vector1< std::pair<Size, Size> > & base_pairing_list )
{

	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::chemical::rna;
	using namespace core::conformation;
	using namespace core::scoring::rna;
	using namespace core::pose::rna;

	// Need some stuff to figure out which residues are base paired. First score.
	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( rna_base_pair, 1.0 );
	(*scorefxn)( pose );

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	EnergyBasePairList const & scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	// bool forms_canonical_base_pair( false );

	BaseEdge k( ANY_BASE_EDGE ), m( ANY_BASE_EDGE );

	std::list < std::pair< Size,Size > > base_pair_list0;

	for ( auto const & it : scored_base_pair_list ) {

		BasePair const & base_pair = it.second;

		Size const i = base_pair.res1();
		Size const j = base_pair.res2();

		if ( i > j ) continue;

		k = base_pair.edge1();
		m = base_pair.edge2();

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK
				&& base_pair.orientation() == ANTIPARALLEL )  &&
				possibly_canonical( rsd_i.aa(), rsd_j.aa() ) ) {
			std::string atom1, atom2;

			if ( !rsd_i.is_coarse() ) { // doesn't work for coarse-grained RNA
				get_watson_crick_base_pair_atoms( rsd_i.type(), rsd_j.type(), atom1, atom2 );
				if ( ( rsd_i.xyz( atom1 ) - rsd_j.xyz( atom2 ) ).length() > 3.5 ) continue;
			}

			base_pair_list0.emplace_back( i, j );
		}
	}


	base_pair_list0.sort();
	base_pairing_list.clear();
	for ( auto const & elem : base_pair_list0 ) {
		base_pairing_list.push_back( elem );
	}
}

// @brief   check that pose has OP1/OP2 with chirality matching Rosetta standard
// @details an ancient function; see also ensure_phosphate_nomenclature_matches_mini()
void
assert_phosphate_nomenclature_matches_mini( pose::Pose const & pose)
{
	runtime_assert( pose.size() > 1 );

	for ( Size res_num=1; res_num<=pose.size(); res_num++ ) {

		Real sign1 = core::pose::rna::get_op2_op1_sign( pose,  res_num);

		pose::Pose mini_pose; //Could move this part outside the for loop
		make_pose_from_sequence( mini_pose, "aa", pose.residue_type_set_for_pose( pose.residue_type(res_num).mode() ) );
		Real const sign2 = core::pose::rna::get_op2_op1_sign( mini_pose);

		if ( sign1 * sign2 < 0 ) {
			utility_exit_with_message("In the assert_phosphate_nomenclature_matches_mini function: phosphate_nomenclature_matches does not match mini!");
			if ( !pose.residue_type(res_num).is_RNA() ) { //Consistency check!
				utility_exit_with_message("residue # " + ObjexxFCL::string_of(res_num)+ " should be a RNA nucleotide!");
			}
		}
	}
}


} //rna
} //protocols
