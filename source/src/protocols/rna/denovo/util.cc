// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file util.cc
/// @brief util functions for use in rna/denovo, i.e. fragment assembly of RNA
/// @details
/// @author Rhiju Das

#include <protocols/rna/denovo/util.hh>
#include <protocols/rna/denovo/secstruct_legacy/RNA_SecStructLegacyInfo.hh>
#include <protocols/rna/denovo/base_pairs/BasePairStep.hh>
#include <protocols/rna/denovo/libraries/RNA_ChunkLibrary.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/forge/methods/fold_tree_functions.hh>
#include <protocols/toolbox/AtomLevelDomainMap.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/rna/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>

// External library headers
#include <utility/io/ozstream.hh>

//C++ headers
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include <protocols/rna/denovo/fragments/RNA_MatchType.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>


using namespace core;
using namespace core::chemical::rna;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
using basic::T;

static THREAD_LOCAL basic::Tracer TR( "protocols.rna.denovo.util" );

using namespace protocols::rna::denovo::secstruct_legacy;
using namespace protocols::rna::denovo::fragments;

namespace protocols {
namespace rna {
namespace denovo {

///////////////////////////////////////////////////////////////////////////////
// A bunch of helper functions used in rna apps...
///////////////////////////////////////////////////////////////////////////////

// @details copies code from get_base_pairing_info [?]. Should be in core/pose/rna/util?
void
get_base_pairing_info( pose::Pose const & pose,
	Size const & seqpos,
	char & secstruct,
	FArray1D <bool> & edge_is_base_pairing ){

	using namespace core::scoring::rna;
	using namespace core::pose::rna;
	using namespace core::chemical;
	using namespace core::conformation;

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	EnergyBasePairList const & scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	edge_is_base_pairing.dimension( 3 );
	edge_is_base_pairing = false;

	bool forms_canonical_base_pair( false ), forms_base_pair( false );

	BaseEdge k( ANY_BASE_EDGE ), m( ANY_BASE_EDGE );
	for ( auto const & it : scored_base_pair_list ) {

		BasePair const & base_pair = it.second;

		Size const i = base_pair.res1();
		Size const j = base_pair.res2();

		if ( i == seqpos ) {
			k = base_pair.edge1();
			m = base_pair.edge2();
		} else if ( j == seqpos ) {
			k = base_pair.edge2();
			m = base_pair.edge1();
		} else {
			continue;
		}

		edge_is_base_pairing( k ) = true;
		forms_base_pair = true;

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

			forms_canonical_base_pair = true;
		}
	}

	secstruct = 'N';
	if ( forms_canonical_base_pair ) {
		secstruct = 'H';
	} else if ( forms_base_pair ) {
		secstruct = 'P';
	}
}



///////////////////////////////////////////////////////////////////////////////
void
figure_out_secstruct( pose::Pose & pose ){
	using namespace core::scoring;
	using namespace ObjexxFCL;

	// Need some stuff to figure out which residues are base paired. First score.
	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( rna_base_pair, 1.0 );
	(*scorefxn)( pose );

	std::string secstruct = "";
	FArray1D < bool > edge_is_base_pairing( 3, false );
	char secstruct1( 'X' );
	for ( Size i=1; i <= pose.size() ; ++i ) {
		get_base_pairing_info( pose, i, secstruct1, edge_is_base_pairing );
		secstruct += secstruct1;
	}

	TR << "SECSTRUCT: " << secstruct << std::endl;

	set_rna_secstruct_legacy( pose, secstruct );
}


///////////////////////////////////////////////////////////////////////////////
void
create_rna_vall_torsions( pose::Pose & pose,
	utility::io::ozstream & torsions_out,
	utility::vector1 <Size> const & exclude_res_list )
{

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace core::pose::rna;
	using namespace protocols::rna::denovo;

	Size const total_residue = pose.size();

	core::pose::rna::figure_out_reasonable_rna_fold_tree( pose );

	// Need some stuff to figure out which residues are base paired. First score.
	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( rna_base_pair, 1.0 );
	(*scorefxn)( pose );
	scorefxn->show( TR, pose );

	bool const idealize_frag( false );
	protocols::idealize::IdealizeMover idealizer;
	idealizer.fast( false /* option[ fast ] */ );

	// EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	for ( Size i=1; i <= total_residue; ++i ) {
		if ( is_num_in_list(i, exclude_res_list) ) continue;

		torsions_out << pose.residue_type( i ).name1() << " " ;

		if ( idealize_frag ) {

			////////////////////////////////////////////////
			// NEW: can we idealize this thing?
			pose::Pose mini_pose;
			Size offset = 1;
			if ( i > 1 ) {
				mini_pose.append_residue_by_bond( pose.residue( i-1 ) );
				offset = 2;
			}
			mini_pose.append_residue_by_bond( pose.residue( i ) );
			if ( i < total_residue ) {
				mini_pose.append_residue_by_bond( pose.residue( i+1 ) );
			}

			core::pose::rna::figure_out_reasonable_rna_fold_tree( mini_pose );
			idealizer.apply( mini_pose );
			idealizer.apply( mini_pose );

			for ( Size j=1; j <= chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
				id::TorsionID my_ID( offset, id::BB, j );
				torsions_out << F( 12, 6, mini_pose.torsion( my_ID ) );
			}

			for ( Size j=1; j <= chemical::rna::NUM_RNA_CHI_TORSIONS; ++j ) {
				id::TorsionID my_ID( offset, id::CHI, j );
				torsions_out << F( 12, 6, mini_pose.torsion( my_ID ) ) << " ";
			}

		} else {
			for ( Size j=1; j <= chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
				id::TorsionID my_ID( i, id::BB, j );
				torsions_out << F( 12, 6, pose.torsion( my_ID ) );
			}

			for ( Size j=1; j <= chemical::rna::NUM_RNA_CHI_TORSIONS; ++j ) {
				id::TorsionID my_ID( i, id::CHI, j );
				torsions_out << F( 12, 6, pose.torsion( my_ID ) ) << " ";
			}

			//New (Feb., 2009) ...
			// x-y-z of coordinates of C2', C1', and O4', in a local coordiante system defined
			// by C3', C4', and C5' (as "stub" atoms).
			conformation::Residue rsd = pose.residue( i );
			chemical::rna::RNA_Info rna_type = rsd.type().RNA_info();
			kinematics::Stub const input_stub( rsd.xyz( rna_type.c3prime_atom_index() ), rsd.xyz( rna_type.c3prime_atom_index() ), rsd.xyz( rna_type.c4prime_atom_index() ), rsd.xyz( rna_type.c5prime_atom_index() ) );

			torsions_out << " S  " ;
			for ( Size n = 1; n <= non_main_chain_sugar_atoms.size(); n++  ) {
				Vector v = input_stub.global2local( rsd.xyz( non_main_chain_sugar_atoms[ n ] ) );
				torsions_out << F( 12, 6, v.x() ) << " " ;
				torsions_out << F( 12, 6, v.y() ) << " " ;
				torsions_out << F( 12, 6, v.z() ) << "  " ;
			}

		}

		// "Secondary structure" ==>
		//    H = forming canonical Watson/Crick base pair
		//    P = forming some kind of base pair
		//    N = no base pairs.
		FArray1D < bool > edge_is_base_pairing( 3, false );
		char secstruct( 'X' );
		get_base_pairing_info( pose, i, secstruct, edge_is_base_pairing );
		torsions_out << secstruct << " " <<
			edge_is_base_pairing( 1 ) << " " <<
			edge_is_base_pairing( 2 ) << " " <<
			edge_is_base_pairing( 3 ) << " ";

		bool is_cutpoint = false;
		if ( pose.fold_tree().is_cutpoint( i ) ||
				is_num_in_list(i + 1, exclude_res_list) ) {
			is_cutpoint = true;
		}

		torsions_out << is_cutpoint << I(6, i)  << std::endl;
	}
}

///////////////////////////////////////////////////////////////////////////////
void
create_rna_vall_torsions( pose::Pose & pose,
	std::string const & outfile,
	utility::vector1 <Size> const & exclude_res_list )
{
	utility::io::ozstream torsions_out ( outfile );
	create_rna_vall_torsions( pose, torsions_out, exclude_res_list );
}


////////////////////////////////////////////////////////////////
//Jan 08, 2012 Parin S:
//To Rhiju, I think this version is buggy and should be deprecated (use make_phosphate_nomenclature_matches_mini instead!)
void
ensure_phosphate_nomenclature_matches_mini( pose::Pose & pose )
{
	runtime_assert( pose.size() > 1 );
	Real sign1 = core::pose::rna::get_op2_op1_sign( pose );

	pose::Pose mini_pose;
	make_pose_from_sequence( mini_pose, "aa", pose.residue_type_set_for_pose() );
	Real sign2 = core::pose::rna::get_op2_op1_sign( mini_pose );

	if ( sign1 * sign2 > 0 ) return;

	TR << "*************************************************************" << std::endl;
	TR << " Warning ... flipping OP2 <--> OP1 to match mini convention  " << std::endl;
	TR << "*************************************************************" << std::endl;

	for ( Size i = 1; i <= pose.size(); i++ ) {

		conformation::Residue const & rsd( pose.residue(i) );
		if ( !rsd.is_RNA() ) continue;

		if ( !rsd.type().has( " OP2") ) continue;
		if ( !rsd.type().has( " OP1") ) continue;

		Vector const temp1 = rsd.xyz( " OP2" );
		Vector const temp2 = rsd.xyz( " OP1" );
		pose.set_xyz( id::AtomID( rsd.atom_index( " OP2" ), i ), temp2 );
		pose.set_xyz( id::AtomID( rsd.atom_index( " OP1" ), i ), temp1 );
	}
}

///////////////////////////////////////////////////////////////////////////////
void
export_packer_results(  utility::vector1< std::pair< Real, std::string > > & results,
	utility::vector1< pose::PoseOP > pose_list,
	scoring::ScoreFunctionOP & scorefxn,
	std::string const & outfile,
	bool const dump )
{
	utility::io::ozstream out( outfile );
	for ( Size n = 1; n <= results.size() ; n++ ) {
		out << F(8,3,results[n].first) << " " << results[n].second << std::endl;
	}
	out.close();

	using namespace core::io::silent;
	SilentFileOptions opts;
	SilentFileData silent_file_data( opts );

	std::string silent_file( outfile );
	Size pos( silent_file.find( ".txt" ) );
	silent_file.replace( pos, 4, ".out" );

	for ( Size n = 1; n <= results.size() ; n++ ) {
		pose::Pose & pose( *pose_list[n] );
		(*scorefxn)( pose );
		std::string const tag( "S_"+lead_zero_string_of( n, 4 ) );
		RNA_SilentStruct s( opts, pose, tag );
		if ( dump ) pose.dump_pdb( tag+".pdb");
		silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );
	}
}


///////////////////////////////////////////////////////////////////////////////
void
check_base_pair( pose::Pose & pose, FArray1D_int & struct_type )
{
	using namespace core::scoring::rna;
	using namespace core::pose::rna;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::conformation;

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );
	(*scorefxn)( pose );

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	EnergyBasePairList const & scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	Size const nres( pose.size() );
	FArray1D_bool forms_noncanonical_base_pair( nres, false );
	FArray1D_bool forms_canonical_base_pair( nres, false );
	FArray1D_int  WC_base_pair_partner( nres, 0 );

	for ( auto const & it : scored_base_pair_list ) {

		BasePair const & base_pair = it.second;

		Size const i = base_pair.res1();
		Size const j = base_pair.res2();

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		bool WC_base_pair( false );
		if ( ( base_pair.edge1() == WATSON_CRICK && base_pair.edge2() == WATSON_CRICK
				&& base_pair.orientation() == ANTIPARALLEL )  &&
				possibly_canonical( rsd_i.aa(), rsd_j.aa() ) )    {
			//std::string atom1, atom2;
			//  get_watson_crick_base_pair_atoms( rsd_i.aa(), rsd_j.aa(), atom1, atom2 );
			//   if ( ( rsd_i.xyz( atom1 ) - rsd_j.xyz( atom2 ) ).length() < 3.5 ) {
			WC_base_pair = true;
			//   }
		}

		TR << rsd_i.name1() << I(3,i) << "--" << rsd_j.name1() << I(3,j) << "  WC:" << WC_base_pair << "  score: " << F(10,6,it.first) << std::endl;

		if ( WC_base_pair ) {
			forms_canonical_base_pair( i ) = true;
			forms_canonical_base_pair( j ) = true;
			WC_base_pair_partner( i ) = j;
			WC_base_pair_partner( j ) = i;
		} else {
			forms_noncanonical_base_pair( i ) = true;
			forms_noncanonical_base_pair( j ) = true;
		}
	}

	struct_type.dimension( nres );
	for ( Size i = 1; i <= nres; i++ ) {
		if ( !forms_noncanonical_base_pair(i) && !forms_canonical_base_pair(i) ) {
			struct_type( i ) = 0;
		} else if ( forms_canonical_base_pair(i) && !forms_noncanonical_base_pair( WC_base_pair_partner(i) )  ) {
			struct_type( i ) = 1;
		} else {
			struct_type( i ) = 2;
		}
	}
}

/////////////////////////////////////////////////////////
// Following only works for coarse-grained RNA poses --
//  perhaps should set up something similar for regular RNA
//  as an alternative to linear_chainbreak.
void
setup_coarse_chainbreak_constraints( pose::Pose & pose, Size const & n )
{
	using namespace core::scoring::constraints;
	using namespace core::scoring::rna;

	if ( !pose.residue_type( n ).has( " S  " ) ) return;
	if ( !pose.residue_type( n ).is_RNA() ) return;

	TR << "Adding CHAINBREAK constraints to " << n << " " << n+1 << std::endl;

	static Real const S_P_distance_min( 3 );
	static Real const S_P_distance_max( 4 );
	static Real const S_P_distance_fade( 1.0 );
	static Real const S_P_distance_bonus( -10.0 );
	static core::scoring::func::FuncOP const S_P_distance_func( new core::scoring::func::FadeFunc( S_P_distance_min - S_P_distance_fade,
		S_P_distance_max + S_P_distance_fade,
		S_P_distance_fade,
		S_P_distance_bonus ) );

	// loose long-range constraint too
	static Real const S_P_distance( 4.0 );
	static Real const S_P_distance_stdev( 10.0 );
	static core::scoring::func::FuncOP const S_P_harmonic_func( new core::scoring::func::HarmonicFunc( S_P_distance,
		S_P_distance_stdev ) );

	static Real const S_S_distance_min( 4.5 );
	static Real const S_S_distance_max( 7.5 );
	static Real const S_S_distance_fade( 2.0 );
	static Real const S_S_distance_bonus( -5.0 );
	static core::scoring::func::FuncOP const S_S_distance_func( new core::scoring::func::FadeFunc( S_S_distance_min - S_S_distance_fade,
		S_S_distance_max + S_S_distance_fade,
		S_S_distance_fade,
		S_S_distance_bonus ) );

	static Real const P_P_distance_min( 4.5 );
	static Real const P_P_distance_max( 7.5 );
	static Real const P_P_distance_fade( 2.0 );
	static Real const P_P_distance_bonus( -5.0 );
	static core::scoring::func::FuncOP const P_P_distance_func( new core::scoring::func::FadeFunc( P_P_distance_min - P_P_distance_fade,
		P_P_distance_max + P_P_distance_fade,
		P_P_distance_fade,
		P_P_distance_bonus ) );

	Size const atom_S1 = pose.residue_type( n ).atom_index( " S  " );
	Size const atom_P1 = pose.residue_type( n ).atom_index( " P  " );
	Size const atom_S2 = pose.residue_type( n+1 ).atom_index( " S  " );
	Size const atom_P2 = pose.residue_type( n+1 ).atom_index( " P  " );

	pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(
		id::AtomID(atom_S1, n),
		id::AtomID(atom_P2, n+1),
		S_P_distance_func, core::scoring::coarse_chainbreak_constraint ) ) ) );

	pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(
		id::AtomID(atom_S1, n),
		id::AtomID(atom_P2, n+1),
		S_P_harmonic_func, core::scoring::coarse_chainbreak_constraint ) ) ) );

	pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(
		id::AtomID(atom_P1, n),
		id::AtomID(atom_P2, n+1),
		P_P_distance_func, core::scoring::coarse_chainbreak_constraint ) ) ) );

	pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new AtomPairConstraint(
		id::AtomID(atom_S1, n),
		id::AtomID(atom_S2, n+1),
		S_S_distance_func, core::scoring::coarse_chainbreak_constraint ) ) ) );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
std::string const
convert_based_on_match_type( std::string const & RNA_string, Size const type ){

	// AMW: this can probably be made more concise using a functional approach.

	std::string RNA_string_local = RNA_string;

	Size const size = RNA_string.length();

	static bool print_warning( false );

	//Obey orders to match exactly, match pyrimidine/purine, or match all.
	if ( type == MATCH_ALL ) {
		for ( Size i = 0; i < size; i++ )  RNA_string_local[ i ] = 'n';
	} else if ( type == MATCH_YR ) {
		for ( Size i = 0; i < size; i++ ) {
			if ( RNA_string[ i ] == 'g' || RNA_string[ i ] == 'a' ) {
				RNA_string_local[ i ] = 'r';
			} else {
				runtime_assert( RNA_string[ i ] == 'u' || RNA_string[ i ] == 'c' || RNA_string[ i ] == 't' );
				RNA_string_local[ i ] = 'y';
			}
		}
	} else {
		for ( Size i = 0; i < size; i++ )  {
			if ( RNA_string[ i ] == 't' ) {
				if ( !print_warning ) {
					TR.Warning << TR.Red << "Requesting an RNA fragment for t. Instead choosing fragment based on u!" << std::endl;
					print_warning = true;
				}
				RNA_string_local[ i ] = 'u';
			}
		}
	}

	return RNA_string_local;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
bool
compare_RNA_char( char const char1, char const char2 ) {
	//Man this is silly, there must be a more elegant way to do this.
	if ( char1 == char2 ) return true;
	if ( char1 == 'n' || char2 == 'n' ) return true;
	if ( char1 == 'r' && (char2 == 'a' || char2 == 'g') ) return true;
	if ( char1 == 'y' && (char2 == 'c' || char2 == 'u') ) return true;
	if ( char2 == 'r' && (char1 == 'a' || char1 == 'g') ) return true;
	if ( char2 == 'y' && (char1 == 'c' || char1 == 'u') ) return true;
	return false;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
bool
compare_RNA_secstruct( char const char1, char const char2 ) {
	if ( char1 == char2 ) return true;
	if ( char1 == 'X' || char2 == 'X' ) return true;
	if ( char1 == 'L' && ( char2 == 'N' || char2 == 'P') ) return true;
	if ( char2 == 'L' && ( char1 == 'N' || char1 == 'P') ) return true;
	return false;
}

///////////////////////////////////////////////////////////////////////////
Vector
get_sugar_centroid( core::conformation::Residue const & rsd ){
	Vector cen( 0.0 );
	core::chemical::rna::RNA_Info rna_type = rsd.type().RNA_info();
	cen += rsd.xyz( rna_type.c1prime_atom_index() );
	cen += rsd.xyz( rna_type.c2prime_atom_index() );
	cen += rsd.xyz( rna_type.c3prime_atom_index() );
	cen += rsd.xyz( rna_type.c4prime_atom_index() );
	cen += rsd.xyz( rna_type.o4prime_atom_index() );
	cen /= 5.0;
	return cen;
}

/////////////////////////////////////////////////
void
make_extended_coarse_pose( pose::Pose & coarse_pose, std::string const & full_sequence ){

	using namespace core::chemical;
	using namespace core::id;

	ResidueTypeSetCOP rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );

	make_pose_from_sequence( coarse_pose, full_sequence, *rsd_set_coarse );

	for ( Size n = 1; n <= coarse_pose.size(); n++ ) {
		coarse_pose.set_torsion( TorsionID( n, BB, 1 ), -150.0 );
		coarse_pose.set_torsion( TorsionID( n, BB, 2 ),  180.0 );
	}
}


///////////////////////////////////////////////////////////////////////////
void
make_coarse_pose( pose::Pose const & pose, pose::Pose & coarse_pose ){

	using namespace core::chemical;
	using namespace core::id;
	using namespace core::scoring::rna;

	ResidueTypeSetCAP rsd_set, rsd_set_coarse;
	RNA_CentroidInfo rna_centroid_info;

	rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );
	make_extended_coarse_pose( coarse_pose, pose.sequence() );

	for ( Size n = 1; n <= pose.size(); n++ ) {
		coarse_pose.set_xyz(  NamedAtomID( " P  ", n ),  pose.xyz( NamedAtomID( " P  ", n )) );

		//coarse_pose.set_xyz(  NamedAtomID( " S  ", n ),  pose.xyz( NamedAtomID( " C4'", n )) );
		Vector sugar_centroid = get_sugar_centroid( pose.residue( n ) );
		coarse_pose.set_xyz(  NamedAtomID( " S  ", n ),  sugar_centroid );

		Vector base_centroid = rna_centroid_info.get_base_centroid( pose.residue( n ) );
		kinematics::Stub stub = rna_centroid_info.get_base_coordinate_system( pose.residue( n ), base_centroid );

		coarse_pose.set_xyz(  NamedAtomID( " CEN", n ),  base_centroid );
		coarse_pose.set_xyz(  NamedAtomID( " X  ", n ),  base_centroid + stub.M.col_x() );
		coarse_pose.set_xyz(  NamedAtomID( " Y  ", n ),  base_centroid + stub.M.col_y() );
	}

	core::kinematics::FoldTree f( pose.fold_tree() );
	for ( Size n = 1; n <= pose.num_jump(); n++ ) {
		f.set_jump_atoms( n, " Y  ", " Y  " );
	}
	coarse_pose.fold_tree( f );
}

////////////////////////////////////////////////////////
void
remove_cutpoint_closed( pose::Pose & pose, Size const i ) {

	using namespace core::chemical;
	using namespace core::kinematics;

	remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, i );
	remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, i+1 );

	// using namespace protocols::forge::methods;
	//  FoldTree f( pose.fold_tree() );
	//  remove_cutpoint( i, f );
	//  pose.fold_tree( f );

	utility::vector1< Size > const & cutpoints = pose.fold_tree().cutpoints();

	Size const num_jump =  pose.fold_tree().num_jump();
	utility::vector1< Size > upstream_pos, downstream_pos;
	for ( Size n = 1; n <= num_jump; n++ ) {
		upstream_pos.push_back( pose.fold_tree().upstream_jump_residue( n ) );
		downstream_pos.push_back( pose.fold_tree().downstream_jump_residue( n ) );
	}

	ObjexxFCL::FArray1D< Size > cuts( num_jump-1 );
	Size count( 0 );
	for ( Size n = 1; n <= num_jump; n++ ) {
		if ( cutpoints[n] == i ) continue;
		count++;
		cuts( count ) = cutpoints[ n ];
	}

	Size const nres( pose.size() );

	// Just brute-force iterate through to find a jump we can remove.
	Size jump_to_remove( 0 );
	for ( Size k = 1; k <= num_jump; k++ ) {
		FArray1D< bool > partition_definition( nres, false );
		pose.fold_tree().partition_by_jump( k, partition_definition );
		if ( partition_definition( i ) != partition_definition( i+1 ) ) {
			jump_to_remove = k; break;
		}
	}

	bool success( false );

	ObjexxFCL::FArray2D< Size > jump_point( 2, num_jump-1 );

	count = 0;
	for ( Size n = 1; n <= num_jump; n++ ) {
		if ( n == jump_to_remove ) continue;
		count++;
		if ( upstream_pos[ n ] < downstream_pos[ n ] ) {
			jump_point( 1, count ) = upstream_pos[ n ];
			jump_point( 2, count ) = downstream_pos[ n ];
		} else {
			jump_point( 1, count ) = downstream_pos[ n ];
			jump_point( 2, count ) = upstream_pos[ n ];
		}
	}

	FoldTree f( nres );

	success = f.tree_from_jumps_and_cuts( nres, num_jump-1, jump_point, cuts, 1, false /*verbose*/ );

	if ( !success ) utility_exit_with_message( "FAIL to remove cutpoint "+string_of( i ) );

	pose.fold_tree( f );
}

////////////////////////////////////////////////////////
void
remove_cutpoints_closed( pose::Pose & pose ){
	// Make a list of each cutpoint_closed.
	for ( Size i = 1; i < pose.size(); i++ ) {
		if ( pose.fold_tree().is_cutpoint( i ) &&
				pose.residue_type( i   ).has_variant_type( chemical::CUTPOINT_LOWER ) &&
				pose.residue_type( i+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ) {
			remove_cutpoint_closed( pose, i ); // this will cycle through to find a jump that is removable.
		}
	}
}

///////////////////////////////////////////////////////////////////////
// One of these days I'll put in residue 1 as well.
void
print_internal_coords( core::pose::Pose const & pose ) {

	using namespace core::id;
	using numeric::conversions::degrees;

	for ( Size i = 2;  i <= pose.size(); i++ ) {

		chemical::ResidueType const & rsd( pose.residue_type( i ) ) ;

		TR << "----------------------------------------------------------------------" << std::endl;
		TR << "RESIDUE: " << rsd.name3() << " " << i << std::endl;

		for ( Size j = 1; j <= rsd.natoms(); j++ ) {
			core::kinematics::tree::AtomCOP current_atom( pose.atom_tree().atom( AtomID(j,i) ).get_self_ptr() );
			core::kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );
			core::kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
			core::kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() );
			if ( !(current_atom && input_stub_atom1 && input_stub_atom2 && input_stub_atom3) )  continue;

			if ( current_atom->is_jump() ) {

				core::kinematics::tree::AtomCOP jump_stub_atom1( current_atom->stub_atom1() );
				core::kinematics::tree::AtomCOP jump_stub_atom2( current_atom->stub_atom2() );
				core::kinematics::tree::AtomCOP jump_stub_atom3( current_atom->stub_atom3() );

				// Now try to reproduce this jump based on coordinates of atoms.
				TR <<   "JUMP! " <<
					//     A( 5, rsd.atom_name( j )) << " " <<
					"STUB ==> " <<
					"FROM " << input_stub_atom1->id().rsd() << " " <<
					pose.residue_type( (input_stub_atom1->id()).rsd() ).atom_name( (input_stub_atom1->id()).atomno() ) << "  " <<
					pose.residue_type( (input_stub_atom2->id()).rsd() ).atom_name( (input_stub_atom2->id()).atomno() ) << "  " <<
					pose.residue_type( (input_stub_atom3->id()).rsd() ).atom_name( (input_stub_atom3->id()).atomno() );
				TR << "TO " << current_atom->id().rsd() << " " <<
					pose.residue_type( (jump_stub_atom1->id()).rsd() ).atom_name( (jump_stub_atom1->id()).atomno() ) << "  " <<
					pose.residue_type( (jump_stub_atom2->id()).rsd() ).atom_name( (jump_stub_atom2->id()).atomno() ) << "  " <<
					pose.residue_type( (jump_stub_atom3->id()).rsd() ).atom_name( (jump_stub_atom3->id()).atomno() ) << std::endl;
				TR << " MY JUMP: " << current_atom->jump() << std::endl;

				kinematics::Stub const input_stub( input_stub_atom1->xyz(), input_stub_atom1->xyz(), input_stub_atom2->xyz(), input_stub_atom3->xyz());
				kinematics::Stub const jump_stub ( jump_stub_atom1->xyz(), jump_stub_atom1->xyz(), jump_stub_atom2->xyz(), jump_stub_atom3->xyz());

				TR << " MY JUMP: " <<  kinematics::Jump( input_stub, jump_stub ) << std::endl;

			} else {
				TR << "ICOOR_INTERNAL  " <<
					ObjexxFCL::format::A( 5, rsd.atom_name( j )) << " " <<
					ObjexxFCL::format::F(11,6, degrees( pose.atom_tree().dof( DOF_ID( current_atom->id(), id::PHI ) ) ) )  << " " <<
					ObjexxFCL::format::F(11,6, degrees(  pose.atom_tree().dof( DOF_ID( current_atom->id(), id::THETA ) ) ) ) << " " <<
					ObjexxFCL::format::F(11,6, pose.atom_tree().dof( DOF_ID( current_atom->id(), id::D ) ) )    << "  " <<
					pose.residue_type( (input_stub_atom1->id()).rsd() ).atom_name( (input_stub_atom1->id()).atomno() ) << "  " <<
					pose.residue_type( (input_stub_atom2->id()).rsd() ).atom_name( (input_stub_atom2->id()).atomno() ) << "  " <<
					pose.residue_type( (input_stub_atom3->id()).rsd() ).atom_name( (input_stub_atom3->id()).atomno() ) << "  " <<
					" [" <<
					" " << (input_stub_atom1->id()).rsd() <<
					" " << (input_stub_atom2->id()).rsd() <<
					" " << (input_stub_atom3->id()).rsd() <<
					"]" <<
					std::endl;
			}
		}

		//   /////////////////////////////////
		//   if ( option[ rsd_type_set]() == "rna" ){
		//    std::cout << "Give me LOWER-P-O5' " << 180.0 - numeric::conversions::degrees( numeric::angle_radians( pose.residue(i-1).xyz( "O3'" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ) ) ) << std::endl;
		//    std::cout << "Give me OP2-P-O5' " << 180.0 - numeric::conversions::degrees( numeric::angle_radians( rsd.xyz( "OP2" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ) ) ) << std::endl;
		//    std::cout << "Give me OP1-P-O5' " << 180.0 - numeric::conversions::degrees( numeric::angle_radians( rsd.xyz( "OP1" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ) ) ) << std::endl;

		//    Real main_torsion = numeric::conversions::degrees(  numeric::dihedral_radians( pose.residue(i-1).xyz( "O3'" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ), rsd.xyz( "C5'" ) ) );
		//    Real main_torsion1 = numeric::conversions::degrees(  numeric::dihedral_radians( pose.residue(i).xyz( "OP2" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ), rsd.xyz( "C5'" ) ) );
		//    Real main_torsion2 = numeric::conversions::degrees(  numeric::dihedral_radians( pose.residue(i).xyz( "OP1" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ), rsd.xyz( "C5'" ) ) );

		//    std::cout << "OFFSETS: " << main_torsion1 - main_torsion << " " << main_torsion2 - main_torsion1 << std::endl;
		//   }

	}
}

////////////////////////////////////////////////////////////////////////////////////
// Apparently can only reroot the tree at a "vertex", i.e. beginning or end of an "edge".
bool
possible_root( core::kinematics::FoldTree const & f, core::Size const & n ){
	return f.possible_root( n );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_rigid_body_jumps( core::pose::Pose const & pose ) {

	utility::vector1< Size > rigid_body_jumps;

	TR.Debug << "Initialize RigidBodyMover: Is last residue virtual? " <<  pose.residue_type( pose.size() ).name3()  << std::endl;

	if ( pose.residue_type( pose.size() ).name3() != "XXX" ) return rigid_body_jumps; // must have a virtual anchor residue.

	for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ) {
		TR.Debug << "checking jump: " <<  pose.fold_tree().upstream_jump_residue( n ) << " to " <<  pose.fold_tree().downstream_jump_residue( n ) << std::endl;
		if ( pose.fold_tree().upstream_jump_residue( n ) == pose.size()  ||
				pose.fold_tree().downstream_jump_residue( n ) == pose.size()  ) {
			TR.Debug << "found jump to virtual anchor at: " << n << std::endl;

			rigid_body_jumps.push_back( n );
			if ( rigid_body_jumps.size() > 1 ) {
				TR.Debug << "found moveable jump!" << std::endl;
			}
		}
	}

	return rigid_body_jumps;
}


////////////////////////////////////////////////////////////////////////////////////
bool
let_rigid_body_jumps_move( core::kinematics::MoveMap & movemap,
	pose::Pose const & pose,
	bool const move_first_rigid_body /* = false */ ){

	utility::vector1< Size > const rigid_body_jumps = get_rigid_body_jumps( pose );
	Size const found_jumps = rigid_body_jumps.size();
	if ( found_jumps <= 1 )  return false; // nothing to rotate/translate relative to another object.

	Size start = ( move_first_rigid_body ) ? 1 : 2;
	for ( Size n = start; n <= rigid_body_jumps.size(); n++ ) movemap.set_jump( rigid_body_jumps[n], true );

	return true;
}

/////////////////////////////////////////////////////////////////////////////////////
Size
get_anchor_rsd( pose::Pose const & pose )
{
	utility::vector1< Size > const rigid_body_jumps = get_rigid_body_jumps( pose );
	if ( rigid_body_jumps.size() == 0 ) return 0;

	Size const nres = pose.size(); // This better be a virtual residue -- checked in get_rigid_body_jumps() above.

	Size anchor_rsd = pose.fold_tree().downstream_jump_residue( rigid_body_jumps[1] );
	if ( anchor_rsd == nres ) anchor_rsd = pose.fold_tree().upstream_jump_residue( rigid_body_jumps[1] );
	return anchor_rsd;
}

/////////////////////////////////////////////////////////////////////////////////////
void
translate_virtual_anchor_to_first_rigid_body( pose::Pose & pose ){

	Size anchor_rsd = get_anchor_rsd( pose );
	if ( anchor_rsd == 0 ) return;

	Size const nres = pose.size(); // This better be a virtual residue -- checked in get_rigid_body_jumps() above.
	Vector anchor1 = pose.xyz( id::AtomID( 1, anchor_rsd ) );
	Vector root1   = pose.xyz( id::AtomID( 1, nres ) );
	Vector const offset = anchor1 - root1;

	for ( Size j = 1; j <= pose.residue_type( nres ).natoms(); j++ ) {
		id::AtomID atom_id( j, nres );
		pose.set_xyz( atom_id, pose.xyz( atom_id ) + offset );
	}
}

////////////////////////////////////////////////////////////////////////////////////////
bool
involved_in_phosphate_torsion( std::string atomname )
{
	utility::vector1< std::string > const & atoms_involved = core::chemical::rna::atoms_involved_in_phosphate_torsion;

	for ( Size n = 1; n <= atoms_involved.size(); n++ ) {
		if (  atomname == atoms_involved[ n ] ) return true;
	}
	return false;
}

///////////////////////////////////////////////////////////////////////////////
void
set_output_res_and_chain( pose::Pose & extended_pose,
	std::pair< utility::vector1< int >, utility::vector1< char > > const & output_resnum_and_chain )
{
	using namespace pose;
	utility::vector1< int > const & output_res_num = output_resnum_and_chain.first;
	utility::vector1< char > const & output_chain = output_resnum_and_chain.second;
	if ( output_res_num.size() == 0 ) return;
	runtime_assert( output_res_num.size() == extended_pose.size() );

	PDBInfoOP pdb_info( new PDBInfo( extended_pose ) );
	pdb_info->set_numbering( output_res_num );

	bool chain_interesting( false );
	for ( Size n = 1; n <= output_chain.size(); n++ ) {
		if ( output_chain[ n ] != ' ' ) chain_interesting = true;
	}
	if ( chain_interesting ) pdb_info->set_chains( output_chain );

	extended_pose.pdb_info( pdb_info );
}


//////////////////////////////////////////////////////////////////////////////////////
void
figure_out_base_pair_partner( pose::Pose & pose, std::map< Size, Size > & partner,
	bool const strict /* = true */ )
{
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace core::pose::rna;
	using namespace core::chemical::rna;
	using namespace core::chemical;
	using namespace core::conformation;

	partner.clear();

	utility::vector1< core::pose::rna::BasePair > base_pair_list = classify_base_pairs_lores( pose );

	for ( Size n = 1; n <= base_pair_list.size(); n++ ) {

		BasePair const & base_pair = base_pair_list[ n ];

		Size const i = base_pair.res1();
		Size const j = base_pair.res2();

		BaseEdge const & k = base_pair.edge1();
		BaseEdge const & m = base_pair.edge2();

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK
				&& base_pair.orientation() == 1 )  &&
				possibly_canonical( rsd_i.aa(), rsd_j.aa() ) &&
				pose.torsion( id::TorsionID( i, id::CHI, 1 ) ) > 0  && //Need to check syn/anti
				pose.torsion( id::TorsionID( j, id::CHI, 1 ) ) > 0     //Need to check syn/anti
				) {

			if ( strict && !possibly_canonical_strict( rsd_i.aa(), rsd_j.aa() ) ) continue;

			std::string atom1, atom2;
			get_watson_crick_base_pair_atoms( rsd_i.type(), rsd_j.type(), atom1, atom2 );
			if ( ( rsd_i.xyz( atom1 ) - rsd_j.xyz( atom2 ) ).length() < 3.5 ) {
				partner[ i ] = j;
				partner[ j ] = i;
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////
utility::vector1< core::pose::rna::BasePair >
classify_base_pairs_lores( pose::Pose const & pose_input )
{

	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace core::pose::rna;
	using namespace core::pose;

	Pose pose = pose_input; // need a non-const copy to apply scorefunction
	ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );
	(*lores_scorefxn)( pose );
	lores_scorefxn->show( TR, pose );

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	EnergyBasePairList scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	utility::vector1< core::pose::rna::BasePair > base_pair_list;
	for ( EnergyBasePairList::const_iterator it = scored_base_pair_list.begin();
			it != scored_base_pair_list.end(); ++it ) base_pair_list.push_back( it->second );

	return base_pair_list;
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
process_input_file( std::string const & input_file,
	utility::vector1< pose::PoseOP > & pose_list,
	bool is_pdb /*= false*/,
	bool coarse_rna /* = false */)
{
	using namespace core::io::silent;
	using namespace protocols::rna::denovo;

	core::chemical::ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	if ( is_pdb ) {

		pose::PoseOP pose_op( new pose::Pose );
		core::import_pose::pose_from_file( *pose_op, *rsd_set, input_file , core::import_pose::PDB_file);
		//   ensure_phosphate_nomenclature_matches_mini( *pose_op );
		core::pose::rna::figure_out_reasonable_rna_fold_tree( *pose_op );
		pose_list.push_back( pose_op );

	} else { //its a silent file.
		SilentFileOptions opts;
		SilentFileData silent_file_data( opts );
		silent_file_data.read_file( input_file );
		for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(),
				end = silent_file_data.end(); iter != end; ++iter ) {
			pose::PoseOP pose_op( new pose::Pose );
			iter->fill_pose( *pose_op );
			pose_list.push_back( pose_op );
		}

	}

	// further cleanup.
	for ( Size n = 1; n <= pose_list.size(); n++ ) {

		pose::PoseOP pose_op = pose_list[ n ];

		remove_cutpoints_closed( *pose_op );

		if ( coarse_rna && !pose_op->residue(1).is_coarse() ) {
			pose::Pose coarse_pose;
			make_coarse_pose( *pose_op, coarse_pose );
			*pose_op = coarse_pose;
		}

		core::pose::rna::virtualize_5prime_phosphates( *pose_op );
	}

	if ( pose_list.size() < 1 )  {
		utility_exit_with_message(  "No structure found in input file  " + input_file );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
void
print_hbonds( pose::Pose & pose ){

	using namespace core::scoring;

	hbonds::HBondOptionsOP hbond_options( new hbonds::HBondOptions() );
	hbond_options->use_hb_env_dep( false );
	hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( *hbond_options ) );

	hbonds::fill_hbond_set( pose, false /*calc deriv*/, *hbond_set );

	for ( Size i = 1; i <= hbond_set->nhbonds(); i++ ) {
		hbonds::HBond const & hbond( hbond_set->hbond( i ) );

		Size const don_res_num = hbond.don_res();
		Size const don_hatm = hbond.don_hatm();

		Size const acc_res_num = hbond.acc_res();
		Size const acc_atm = hbond.acc_atm();

		TR << "HBOND: " << pose.residue_type( don_res_num ).name1() << don_res_num <<
			" " << pose.residue_type( don_res_num ).atom_name( don_hatm ) << " --- " <<
			pose.residue_type( acc_res_num).name1() << acc_res_num << " " << pose.residue_type( acc_res_num ).atom_name( acc_atm ) << " ==> " << hbond.energy()
			<< std::endl;
	}
}


/////////////////////////////////////////////////////////
bool
moveable_jump( id::AtomID const & jump_atom_id1,
	id::AtomID const & jump_atom_id2,
	protocols::toolbox::AtomLevelDomainMap const & atom_level_domain_map)  {
	if ( atom_level_domain_map.get( jump_atom_id1 ) ) return true;
	if ( atom_level_domain_map.get( jump_atom_id2 ) ) return true;
	if ( atom_level_domain_map.get_domain( jump_atom_id1 ) != atom_level_domain_map.get_domain( jump_atom_id2 ) ) return true;
	return false;
}

/////////////////////////////////////////////////////////
// should this really be separate from AtomID-base moveable_jump?
bool
moveable_jump( Size const jump_pos1,
	Size const jump_pos2,
	protocols::toolbox::AtomLevelDomainMap const & atom_level_domain_map)  {
	if ( atom_level_domain_map.get( jump_pos1 ) ) return true;
	if ( atom_level_domain_map.get( jump_pos2 ) ) return true;
	if ( atom_level_domain_map.get_domain( jump_pos1 ) != atom_level_domain_map.get_domain( jump_pos2 ) ) return true;
	return false;
}

///////////////////////////////////////////////////////////////
void
fill_in_default_jump_atoms( kinematics::FoldTree & f, pose::Pose const & pose )
{

	// default atoms for jump connections.
	for ( Size i = 1; i <= f.num_jump(); i++ ) {
		Size const jump_pos1( f.upstream_jump_residue( i ) );
		Size const jump_pos2( f.downstream_jump_residue( i ) );
		f.set_jump_atoms( i,
			chemical::rna::default_jump_atom( pose.residue_type( jump_pos1 ) ),
			chemical::rna::default_jump_atom( pose.residue_type( jump_pos2) ) );
	}
}

///////////////////////////////////////////////////////////////
void
fill_in_default_jump_atoms( pose::Pose & pose )
{
	core::kinematics::FoldTree f( pose.fold_tree() );
	fill_in_default_jump_atoms( f, pose );
	pose.fold_tree( f );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
get_default_allowed_bulge_res(
	utility::vector1< core::Size > & allow_bulge_res_list,
	core::pose::Pose const & pose,
	bool const verbose ){

	if ( allow_bulge_res_list.size() != 0 ) {
		utility_exit_with_message( "allow_bulge_res_list.size() != 0" );
	}

	if ( verbose ) {
		TR << "allow_bulge_res_list.size() == 0, ";
		TR << "Getting default_allowed_bulge_res!" << std::endl;
	}

	for ( Size seq_num = 1; seq_num <= pose.size(); seq_num++ ) {

		//exclude edge residues:
		if ( seq_num == 1 ) continue;
		if ( seq_num == pose.size() ) continue;

		//bool is_cutpoint_closed=false;

		bool is_cutpoint_lower = pose.residue_type( seq_num ).has_variant_type(
			chemical::CUTPOINT_LOWER );
		bool is_cutpoint_upper = pose.residue_type( seq_num ).has_variant_type(
			chemical::CUTPOINT_UPPER );

		bool near_cutpoint_closed = is_cutpoint_lower || is_cutpoint_upper;

		bool near_cutpoint = pose.fold_tree().is_cutpoint( seq_num ) ||
			pose.fold_tree().is_cutpoint( seq_num - 1 );

		bool near_cutpoint_open = near_cutpoint && !near_cutpoint_closed;

		if ( near_cutpoint_open ) continue;

		allow_bulge_res_list.push_back( seq_num );
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
core::Size
virtualize_bulges( core::pose::Pose & input_pose,
	utility::vector1< core::Size > const & in_allow_bulge_res_list,
	core::scoring::ScoreFunctionCOP const & scorefxn,
	std::string const & tag,
	bool const allow_pre_virtualize,
	bool const allow_consecutive_bulges,
	bool const verbose ){

	using namespace core::pose;
	using namespace core::scoring;
	using namespace ObjexxFCL;


	Size const total_res = input_pose.size();

	Real const rna_bulge_bonus = ( scorefxn->get_weight( rna_bulge ) )*10;

	utility::vector1< core::Size > allow_bulge_res_list = in_allow_bulge_res_list;

	if ( allow_bulge_res_list.size() == 0 ) {
		get_default_allowed_bulge_res( allow_bulge_res_list, input_pose, verbose );
	}

	if ( verbose ) {
		TR << "Enter virtualize_bulges() " << std::endl;
		TR << "rna_bulge_bonus = " << rna_bulge_bonus << std::endl;
		//  output_boolean( "allow_pre_virtualize = ", allow_pre_virtualize, TR ); TR << std::endl;
		//  output_boolean( "allow_consecutive_bulges = ", allow_consecutive_bulges, TR ); TR << std::endl;
		//  output_seq_num_list( "allow_bulge_res_list = ", allow_bulge_res_list, TR );

		//Testing to see if checking in core::pose::rna::apply_virtual_rna_residue_variant_type can be violated!//////////
		pose::Pose testing_pose = input_pose;

		for ( Size seq_num = 1; seq_num <= total_res; seq_num++ ) {
			if ( !testing_pose.residue_type( seq_num ).is_RNA() ) continue; //Fang's electron density code
			if ( allow_bulge_res_list.has_value( seq_num ) == false ) continue;
			core::pose::rna::apply_virtual_rna_residue_variant_type( testing_pose, seq_num, true /*apply_check*/ );
		}
		//////////////////////////////////////////////////////////////////////////////////////////////////
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	if ( !allow_pre_virtualize ) {
		for ( Size seq_num = 1; seq_num <= total_res; seq_num++ ) {
			if ( !input_pose.residue_type( seq_num ).is_RNA() ) continue; //Fang's electron density code
			if ( input_pose.residue_type( seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
				utility_exit_with_message( "allow_pre_virtualize == false but seq_num = " + string_of( seq_num ) +
					"  is already virtualized!!" );
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	pose::Pose working_pose = input_pose;
	Real const start_score = ( *scorefxn )( working_pose );


	Size num_res_virtualized = 0;
	Size round_num = 0;

	while ( true ) { //Need multiple round, since virtualizing a particular res will reduce the energy score of neighoring res and might effect wether neighoring res should be virtualized.
		pose::Pose base_pose = input_pose;
		Real const base_score = ( *scorefxn )( base_pose );

		round_num++;

		Size num_res_virtualized_in_this_round = 0;

		for ( Size seq_num = 1; seq_num <= total_res; seq_num++ ) {
			if ( !input_pose.residue_type( seq_num ).is_RNA() ) continue; //Fang's electron density code
			if ( !allow_bulge_res_list.has_value( seq_num ) ) continue;

			if ( input_pose.residue_type( seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
				if ( ! input_pose.residue_type( seq_num + 1 ).has_variant_type( core::chemical::VIRTUAL_PHOSPHATE ) ) { //consistency_check
					utility_exit_with_message( "seq_num = " + string_of( seq_num ) + "  is a virtual res but seq_num + 1 is not a virtual_res_upper!" );
				}

				if ( ! base_pose.residue_type( seq_num ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) { //consistency check
					utility_exit_with_message( "input_pose have virtual at seq_num = " + string_of( seq_num ) + "  but input_pose doesn't!" );
				}

				continue;
			}

			if ( allow_consecutive_bulges == false ) {
				if ( ( seq_num + 1 ) <= total_res ) {
					if ( input_pose.residue_type( seq_num + 1 ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
						continue;
					}
				}

				if ( ( seq_num - 1 ) >= 1 ) {
					if ( input_pose.residue_type( seq_num - 1 ).has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
						continue;
					}
				}
			}


			working_pose = base_pose; //reset working_pose to base_pose
			core::pose::rna::apply_virtual_rna_residue_variant_type( working_pose, seq_num, true ) ;
			Real const new_score = ( *scorefxn )( working_pose );

			if ( new_score < base_score ) {
				num_res_virtualized++;
				num_res_virtualized_in_this_round++;

				TR << "tag = " << tag << " round_num = " << round_num << " seq_num = " << seq_num << ". new_score ( " << new_score << " ) is lesser than base_score ( " << base_score << " ). " << std::endl;

				core::pose::rna::apply_virtual_rna_residue_variant_type( input_pose, seq_num, true ) ;

			}
		}

		if ( num_res_virtualized_in_this_round == 0 ) break;
	}


	working_pose = input_pose;
	Real const final_score = ( *scorefxn )( working_pose );

	if ( num_res_virtualized > 0 ) {
		TR << "----------------------------------------------------------" << std::endl;
		TR << "Inside virtualize_bulges() " << std::endl;
		TR << "TOTAL_NUM_ROUND = " << round_num << std::endl;
		TR << "tag = " << tag << std::endl;
		TR << "num_res_virtualized = " << num_res_virtualized << std::endl;
		TR << "start_score = " << start_score << std::endl;
		TR << "final_score = " << final_score << std::endl;
		TR << "----------------------------------------------------------" << std::endl;

	}

	return num_res_virtualized;
}

//////////////////////////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionOP
get_rna_hires_scorefxn() {
	core::scoring::ScoreFunctionOP scorefxn_;
	if ( basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
		scorefxn_ = core::scoring::get_score_function();
	} else {
		scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::RNA_HIRES_WTS );
	}
	return scorefxn_;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_moving_res( core::pose::Pose const & pose,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map ) {

	utility::vector1< Size > moving_res;

	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( atom_level_domain_map->get( n ) ) {
			moving_res.push_back( n );
		}
	}

	return moving_res;
}


////////////////////////////////////////////////////////////////////////////////////////
bool
base_pair_step_moving( base_pairs::BasePairStep const & base_pair_step,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
	pose::Pose const & pose )
{
	// look that at least one base in this base pair step is moveable.
	bool pair_moving( false );
	utility::vector1< Size > bps_res = utility::tools::make_vector1( base_pair_step.i(),
		base_pair_step.i_next(),
		base_pair_step.j(),
		base_pair_step.j_next() );
	for ( Size i = 1; i <= 4; i++ ) {
		Size domain(  atom_level_domain_map->get_domain( core::id::NamedAtomID( " C1'", bps_res[i] ), pose  ) );
		// TR << "DOMAIN for " << bps_res[i] << " --> " << domain << std::endl;
		if ( domain == 0 || domain == libraries::ROSETTA_LIBRARY_DOMAIN ) {
			pair_moving = true; //break;
		}
	}
	// TR << "MOVING? " << pair_moving << std::endl;
	return pair_moving;
}

/////////////////////////////////////////////////////////////////////////////////
bool
base_pair_moving( core::pose::rna::BasePair const & base_pair,
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map,
	core::pose::Pose const & pose )
{
	// look that at least one base in this base pair step is moveable.
	bool pair_moving( false );
	utility::vector1< Size > bp_res = utility::tools::make_vector1( base_pair.res1(),
		base_pair.res2() );
	for ( Size i = 1; i <= bp_res.size(); i++ ) {
		Size domain(  atom_level_domain_map->get_domain( core::id::NamedAtomID( " C1'", bp_res[i] ), pose  ) );
		//  TR << "DOMAIN for " << bp_res[i] << " --> " << domain << std::endl;
		if ( domain == 0 || domain == libraries::ROSETTA_LIBRARY_DOMAIN ) {
			pair_moving = true; //break;
		}
	}
	// TR << "MOVING? " << pair_moving << std::endl;
	return pair_moving;
}

/////////////////////////////////////////////////////////////////////////////////
void
delete_non_protein_from_pose( core::pose::Pose & pose ) {

	// Delete all residues that aren't protein from a pose

	utility::vector1< std::pair< core::Size, core::Size > > rna_start_and_end_residues;
	core::Size start = 0;
	core::Size end = 0;
	// Figure out the non-protein start and end residues
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		if ( !pose.residue( i ).is_protein() && start == 0 ) {
			start = i;
		} else if ( pose.residue( i ).is_protein() ) {
			if ( start != 0 ) {
				end = i-1;
				rna_start_and_end_residues.push_back( std::make_pair( start, end ) );
				start = 0;
			}
		}
		if ( (i == pose.total_residue()) && (start > end) ) {
			end = i;
			rna_start_and_end_residues.push_back( std::make_pair( start, end ) );
		}
	}

	// Now delete all the non-protein residues
	core::Size offset = 0;
	for ( core::Size i=1; i<=rna_start_and_end_residues.size(); ++i ) {
		pose.delete_residue_range_slow( rna_start_and_end_residues[i].first - offset, rna_start_and_end_residues[i].second - offset );
		offset += rna_start_and_end_residues[i].second - rna_start_and_end_residues[i].first + 1;
	}

}

} //denovo
} //rna
} //protocols
