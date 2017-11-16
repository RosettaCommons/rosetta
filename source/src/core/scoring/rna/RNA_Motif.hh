// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/rna/RNA_Motif.hh
/// @brief find RNA motifs (u-turns, T-loops, etc.) in a pose
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_rna_RNA_MotifFinder_HH
#define INCLUDED_core_scoring_rna_RNA_MotifFinder_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/scoring/rna/RNA_Motif.fwd.hh>

#include <core/chemical/rna/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/rna/RNA_BasePairClassifier.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/RNA_FilteredBaseBaseInfo.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>

#include <boost/iterator/indirect_iterator.hpp>


////////////////////////////////////////////////////////////////////////////////////////////////
/// @detailed
///
/// Defines in detail the classic motifs & submotifs that
///   recur in RNA folding, collated from numerous sources.
///
/// TODO:
///
///   Move this into core::pose::rna. To do that, remove
///    dependency on core::scoring::rna as following:
///    - Need to take as input RNA_BaseStackList & RNA_BasePairList instead
///         of scoring::rna::FilteredBaseBaseInfo
///    - Need to write a helper function to look for 2'-OH or O1P contacts
///        with nucleobases instead of using scoring::rna::RNA_LowResolutionPotential
///
///   Maybe: Devise 'bonus' system for nice features (e.g., "canonical" sequence), but
///     use mainly geometry to define obligate features. Depends on results of benchmarks.
///
////////////////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {

enum RNA_MotifType
{ U_TURN,
	UA_HANDLE,
	T_LOOP,
	INTERCALATED_T_LOOP,
	LOOP_E_SUBMOTIF,
	BULGED_G,
	GNRA_TETRALOOP,
	STRICT_WC_STACKED_PAIR,
	WC_STACKED_PAIR,
	A_MINOR,
	PLATFORM,
	TL_RECEPTOR,
	TETRALOOP_TL_RECEPTOR
};

std::map< RNA_MotifType, core::Real > const rna_motif_bonus =
{ {U_TURN, -1.0},
{UA_HANDLE, -1.0},
{T_LOOP, -1.0 },
{INTERCALATED_T_LOOP, -1.0 },
{LOOP_E_SUBMOTIF, -1.0 },
{BULGED_G, -1.0 },
{GNRA_TETRALOOP, -0.1},  // GNRA is already pretty favorable.
{STRICT_WC_STACKED_PAIR, 0.0},
{WC_STACKED_PAIR, 0.0},
{A_MINOR, -1.0},
{PLATFORM, -1.0},
{TL_RECEPTOR, -1.0},
{TETRALOOP_TL_RECEPTOR, -1.0}
};


std::map< RNA_MotifType, std::string > const motif_color =
{ {U_TURN, "lightblue"},
{UA_HANDLE, "marine" },
{T_LOOP, "tv_blue" },
{INTERCALATED_T_LOOP, "deepblue" },
{LOOP_E_SUBMOTIF, "salmon" },
{BULGED_G, "red" },
{GNRA_TETRALOOP, "ruby"},  // GNRA is already pretty favorable.
{STRICT_WC_STACKED_PAIR, "gray20"},
{WC_STACKED_PAIR, "gray50"},
{A_MINOR, "gold"},
{PLATFORM, "sand"},
{TL_RECEPTOR, "limon"},
{TETRALOOP_TL_RECEPTOR, "orange"}
};

inline
std::string
to_string( RNA_MotifType type ) {
	switch ( type ) {
	case U_TURN : return "U_TURN";
	case UA_HANDLE : return "UA_HANDLE";
	case T_LOOP : return "T_LOOP";
	case INTERCALATED_T_LOOP : return "INTERCALATED_T_LOOP";
	case LOOP_E_SUBMOTIF : return "LOOP_E_SUBMOTIF";
	case BULGED_G : return "BULGED_G";
	case GNRA_TETRALOOP : return "GNRA_TETRALOOP";
	case STRICT_WC_STACKED_PAIR : return "STRICT_WC_STACKED_PAIR";
	case WC_STACKED_PAIR : return "WC_STACKED_PAIR";
	case A_MINOR : return "A_MINOR";
	case PLATFORM : return "PLATFORM";
	case TL_RECEPTOR : return "TL_RECEPTOR";
	case TETRALOOP_TL_RECEPTOR : return "TETRALOOP_TL_RECEPTOR";
	default :
		utility_exit_with_message( "Unrecognized type" );
	}
	return "";
}


// @brief RNA_Motif has a type (e.g., U_TURN) and residues (in pose numbering).
class RNA_Motif: public utility::pointer::ReferenceCount {

public:

	//constructor
	RNA_Motif( RNA_MotifType const & rna_motif_type,
		utility::vector1< Size > const & residues ) :
		type_( rna_motif_type ),
		residues_( residues )
	{}

	//destructor
	~RNA_Motif()
	{}

public:

	void set_type( RNA_MotifType const & setting ){ type_ = setting; }
	RNA_MotifType type() const { return type_; }

	void set_residues( utility::vector1< core::Size > const & setting ){ residues_ = setting; }
	utility::vector1< core::Size > residues() const { return residues_; }

	core::Size
	operator[]( core::Size const index ) const { return residues_[ index ]; }

	friend
	std::ostream &
	operator << ( std::ostream & out, RNA_Motif const & s );

	friend
	bool operator < ( RNA_Motif const & lhs, RNA_Motif const & rhs );

	utility::vector1< core::Size >::const_iterator begin() const { return residues_.begin(); }
	utility::vector1< core::Size >::const_iterator end() const { return residues_.end(); }

private:

	RNA_MotifType type_;
	utility::vector1< core::Size > residues_;

};

inline
std::ostream &
operator << ( std::ostream & out, RNA_Motif const & s ) {
	out << to_string( s.type_ ) << ' ' << s.residues_;
	return out;
}

///////////////////////////////////////////////////////////////////
inline
bool operator < ( RNA_Motif const & lhs, RNA_Motif const & rhs )
{
	return ( lhs[1] < rhs[1] );
}

// @brief Collection of RNA motifs. get_motifs() can get sub-list of just particular type.
class RNA_Motifs: public utility::pointer::ReferenceCount {

public:

	//constructor
	RNA_Motifs()
	{}

	//destructor
	~RNA_Motifs()
	{}

public:

	void
	push_back( RNA_Motif const & rna_motif ) {
		rna_motifs_.push_back( rna_motif );
		rna_motif_map_[ rna_motif.type() ].push_back( rna_motif );
	}

	utility::vector1< RNA_Motif > const &
	get_motifs() const { return rna_motifs_; }

	utility::vector1< RNA_Motif >::const_iterator begin() const { return rna_motifs_.begin(); }
	utility::vector1< RNA_Motif >::const_iterator end() const { return rna_motifs_.end(); }

	utility::vector1< RNA_Motif > const &
	get_motifs( RNA_MotifType const & type ) const {
		return rna_motif_map_[ type ];
	}

	friend
	std::ostream &
	operator << ( std::ostream & out, RNA_Motifs const & s );

private:

	utility::vector1< RNA_Motif > rna_motifs_;
	mutable std::map< RNA_MotifType, utility::vector1< RNA_Motif> > rna_motif_map_;

};

inline
std::ostream &
operator << ( std::ostream & out, RNA_Motifs const & s ) {
	for ( auto const & motif : s ) out << motif << std::endl;
	return out;
}


// @brief check that n, n+1, ... n+nloop-1 is a continuous RNA loop.
inline
bool
check_rna_loop( core::pose::Pose const & pose,
	core::Size const & n,
	core::Size const & nloop )
{
	using namespace core::chemical;
	if ( n < 1 ) return false;
	for ( Size i = n; i <= n + nloop - 1; i++ ) {
		if ( i > pose.size() ) return false;
		if ( !pose.residue( i ).is_RNA() ) return false;
	}
	for ( Size i = n; i < n + nloop - 1; i++ ) {
		if ( pose.fold_tree().is_cutpoint( i ) &&
				!( pose.residue( i ).has_variant_type( CUTPOINT_LOWER ) &&
				pose.residue( i+1 ).has_variant_type( CUTPOINT_UPPER ) ) ) return false;
	}
	return true;
}

// @brief check that i stacks on j (i<j)
inline
bool
check_stack( Size const & i, Size const & j,
	core::pose::rna::RNA_BaseStackList const & base_stack_list )
{
	for ( auto const & base_stack : base_stack_list ) {
		if ( base_stack.res1() == i && base_stack.res2() == j ) {
			return true;
		}
	}
	return false;
}

// @brief output for RNA_Motif
inline
void
output_rna_motif( pose::Pose const & pose,
	RNA_Motif const & motif )
{
	if ( rna_motif_bonus.count( motif.type() ) == 0 ) return;
	if ( rna_motif_bonus.find( motif.type() )->second >= 0.0  ) return;
	std::cout << ObjexxFCL::right_string_of(to_string(motif.type()),20) << ": ";
	for ( auto const & res : motif ) {
		std::cout << ' ' << ObjexxFCL::right_string_of( pose.residue( res ).annotated_name(), 6 ) <<
			"-" <<pose.pdb_info()->chain( res ) << ":";
		if ( pose.pdb_info()->segmentID( res ).size() > 0 &&
				pose.pdb_info()->segmentID( res ) != "    " ) std::cout << pose.pdb_info()->segmentID( res ) << ":";
		std::cout << ObjexxFCL::left_string_of(pose.pdb_info()->number( res ),4);
	}
	std::cout << std::endl;
}

// @brief output for RNA_Motifs (with detailed residue names, etc.)
inline
void
output_rna_motifs_detailed( pose::Pose const & pose,
	RNA_Motifs const & motifs )
{
	for ( auto const & motif : motifs.get_motifs() ) output_rna_motif( pose, motif );
}

// @brief check if sequence is compatible with strict W/C (G-C, A-U, *not* G-U).
inline
bool
check_watson_crick_sequence( pose::Pose const & pose, Size const & res1, Size const & res2,
	bool const strict = false ) {
	using namespace core::chemical;
	AA aa1( pose.residue( res1 ).aa() );
	if ( aa1 == aa_unp || aa1 == aa_unk ) aa1 = pose.residue( res1 ).na_analogue();
	AA aa2( pose.residue( res2 ).aa() );
	if ( aa2 == aa_unp || aa2 == aa_unk ) aa2 = pose.residue( res2 ).na_analogue();

	if ( aa1 == na_rcy && aa2 == na_rgu ) return true;
	if ( aa2 == na_rcy && aa1 == na_rgu ) return true;

	if ( aa1 == na_ura && aa2 == na_rad ) return true;
	if ( aa2 == na_ura && aa1 == na_rad ) return true;

	if ( !strict ) {
		if ( aa1 == na_ura && aa2 == na_rgu ) return true;
		if ( aa2 == na_ura && aa1 == na_rgu ) return true;
	}

	return false;
}

/// @brief identify RNA motifs inside pose
/// @details this version allows user to supply preinstantiated potential
///          and precalculated filtered_base_base_info for speed.
inline
RNA_Motifs
get_rna_motifs( pose::Pose const & pose,
	core::scoring::rna::RNA_LowResolutionPotential const & potential,
	core::pose::rna::RNA_FilteredBaseBaseInfo const & filtered_base_base_info )
{
	using namespace core::pose::rna;
	using namespace core::scoring::rna;
	using namespace core::chemical::rna;
	using namespace core::chemical;
	using namespace utility::tools;

	RNA_Motifs rna_motifs;

	// we should pre-cache this information.
	utility::vector1< core::pose::rna::BasePair> base_pairs;
	RNA_BasePairList  base_pair_list  = filtered_base_base_info.base_pair_list();
	for ( auto const & base_pair : base_pair_list ) {
		base_pairs.push_back( base_pair );
		base_pairs.push_back( base_pair.flipped() );
	}
	utility::vector1< core::pose::rna::BaseStack> base_stacks;
	RNA_BaseStackList base_stack_list = filtered_base_base_info.base_stack_list();
	for ( auto const & base_stack : base_stack_list ) {
		base_stacks.push_back( base_stack );
		base_stacks.push_back( base_stack.flipped() );
	}

	/////////////////////////////////////////
	// Straight up W/C base pair step?
	/////////////////////////////////////////
	for ( auto const & base_pair1  : base_pairs ) {
		if ( base_pair1.edge1() != WATSON_CRICK )  continue;
		if ( base_pair1.edge2() != WATSON_CRICK )  continue;
		if ( base_pair1.orientation() != ANTIPARALLEL )  continue;
		if ( !check_watson_crick_sequence( pose, base_pair1.res1(), base_pair1.res2() ) ) continue;
		for ( auto const & base_pair2  : base_pairs ) {
			if ( base_pair2.res1() != base_pair1.res1() + 1 ) continue;
			if ( base_pair1.res2() != base_pair2.res2() + 1 ) continue;
			if ( pose.pdb_info() && pose.pdb_info()->number( base_pair2.res1() ) != pose.pdb_info()->number( base_pair1.res1() ) + 1 ) continue;
			if ( pose.pdb_info() && pose.pdb_info()->chain(  base_pair2.res1() ) != pose.pdb_info()->chain(  base_pair1.res1() ) ) continue;
			if ( pose.pdb_info() && pose.pdb_info()->number( base_pair1.res2() ) != pose.pdb_info()->number( base_pair2.res2() ) + 1 ) continue;
			if ( pose.pdb_info() && pose.pdb_info()->chain(  base_pair1.res2() ) != pose.pdb_info()->chain(  base_pair2.res2() ) ) continue;
			if ( base_pair2.edge1() != WATSON_CRICK )  continue;
			if ( base_pair2.edge2() != WATSON_CRICK )  continue;
			if ( base_pair2.orientation() != ANTIPARALLEL )  continue;
			if ( !check_watson_crick_sequence( pose, base_pair2.res1(), base_pair2.res2() ) ) continue;
			// Following ends up being too stringent:
			//  if ( !check_stack( base_pair2.res1(), base_pair1.res1(), base_stacks ) ) continue;
			// if ( !check_stack( base_pair2.res2(), base_pair1.res2(), base_stacks ) ) continue;
			rna_motifs.push_back( RNA_Motif( WC_STACKED_PAIR,
				{ base_pair1.res1(), base_pair2.res1(), base_pair2.res2(), base_pair1.res2() } ) );
		}
	}

	/////////////////////////
	//  U-turns
	//
	//     2
	//     N
	//    /  R 3
	// 1 U .. p
	//
	// See, e.g.,
	//
	//   Moore, Annual Review Biochem., 1999.
	//
	/////////////////////////
	for ( Size n = 1; n <= pose.size(); n++ ) {
		utility::vector1< Size > u_turn;
		if ( n > pose.size() - 3 ) continue;
		if ( !check_rna_loop( pose, n, 4 ) ) continue;
		if ( potential.get_base_backbone( pose.residue( n ), pose.residue( n + 3 ), 1 /* O1P */) >= 0.0 &&
				potential.get_base_backbone( pose.residue( n ), pose.residue( n + 3 ), 2 /* O2P */) >= 0.0 ) continue;
		if ( potential.get_base_backbone( pose.residue( n+2 ), pose.residue( n ), 6 /* O2' */) >= 0.0 )  continue;
		// should we also enforce stacking n+1, n+2?
		for ( Size i = n; i <= n+2; i++ ) u_turn.push_back( i );
		rna_motifs.push_back( RNA_Motif( U_TURN, u_turn ) );
	}


	/////////////////////////////////////////
	//  UA-handles
	//
	// 2 U o-@ A 3
	//   |      N <-- [could be lots of N's]
	// 1 G -o- C 4
	//
	/////////////////////////////////////////
	for ( Size n = 2; n <= pose.size(); n++ ) {
		for ( auto const & base_pair : base_pairs ) {
			// Look for the U-A Watson-Crick/Hoogsteen pair
			if ( base_pair.res1() == n && base_pair.edge1() == WATSON_CRICK && base_pair.edge2() == HOOGSTEEN ) {
				utility::vector1< Size > ua_handle;
				// Look for W/C pair stacked immediately previous, connected to U
				if ( !check_stack( n-1, n, base_stacks ) ) continue;
				for ( auto const & base_pair_prev : base_pairs ) {
					if ( base_pair_prev.res1() == n-1 && base_pair_prev.edge1() == WATSON_CRICK && base_pair_prev.edge2() == WATSON_CRICK ) {
						if ( base_pair_prev.res2() == base_pair.res2()+1 ) continue;
						if ( !check_stack( base_pair.res2(), base_pair_prev.res2(), base_stacks ) ) continue;
						ua_handle.push_back( base_pair_prev.res1() );
						ua_handle.push_back( base_pair.res1() /* n */);
						ua_handle.push_back( base_pair.res2() );
						ua_handle.push_back( base_pair_prev.res2() );
						rna_motifs.push_back( RNA_Motif( UA_HANDLE, ua_handle ) );
						break;
					}
				}
			}
		}
	}

	/////////////////////////////////////////
	//  T-loop
	//
	//     4
	//     N
	//    /  R 5
	// 3 U .. p
	//   |     \                           //
	// 2 U o-@ A 6
	//   |      N <-- [could be lots of N's]
	// 1 G -o- C 7
	//
	/////////////////////////////////////////
	for ( auto const & u_turn : rna_motifs.get_motifs( U_TURN ) ) {
		for ( auto const & ua_handle : rna_motifs.get_motifs( UA_HANDLE ) ) {
			// right on top:
			if ( u_turn[1] != ua_handle[2]+1 ) continue;
			if ( u_turn[3] != ua_handle[3]-1 ) continue;
			// check for stacking
			if ( !check_stack( ua_handle[2], u_turn[1], base_stacks ) ) continue;
			// check for pocket (lack of stacking!) that would allow intercalation
			// (this should probably be a bonus feature, not obligate)
			if ( check_stack(  ua_handle[3], u_turn[3], base_stacks ) ) continue;
			utility::vector1< Size > t_loop = make_vector1( ua_handle[1],ua_handle[2],u_turn[1],u_turn[2],u_turn[3],
				ua_handle[3],ua_handle[4]);
			rna_motifs.push_back( RNA_Motif( T_LOOP, t_loop ) );
		}
	}

	/////////////////////////////////////////
	//  Intercalated T-loops
	//
	//     4
	//     N
	//    /  R 5
	// 3 U    X 8 <--- intercalated base from *outside* T-loop
	//   |     \                           //
	// 2 U o-@ A 6
	//   |      N <-- [could be lots of N's]
	// 1 G -o- C 7
	//
	/////////////////////////////////////////
	for ( auto const & t_loop : rna_motifs.get_motifs( T_LOOP ) ) {
		for ( auto const & base_stack : base_stacks ) {
			if ( base_stack.res1() == t_loop[ 5 ] ) {
				Size const & intercalator = base_stack.res2();
				for ( auto const & base_stack2 : base_stacks ) {
					if ( base_stack2.res1() == intercalator &&
							base_stack2.res2() == t_loop[ 6 ] ) {
						utility::vector1< Size > intercalated_t_loop = t_loop.residues();
						intercalated_t_loop.push_back( intercalator );
						rna_motifs.push_back( RNA_Motif( INTERCALATED_T_LOOP, intercalated_t_loop ) );
						break;
					}
				}
			}
		}
	}


	/////////////////////////////////////////
	// GNRA tetraloop
	//
	//     2
	//     N
	//    /  R 3
	// 1 G .. p A 4
	//   G -o- C
	//
	// Following generalizes to GNRA with bulges, e.g.,
	//  pentaloop in SARS virus, GAGUA. Bulge must be
	//  after U-turn (4) and before WC (6).
	/////////////////////////////////////////
	for ( auto const & u_turn : rna_motifs.get_motifs( U_TURN ) ) {
		for ( auto const & base_pair : base_pairs ) {
			if ( base_pair.edge1() != WATSON_CRICK ) continue;
			if ( base_pair.edge2() != WATSON_CRICK ) continue;
			if ( base_pair.res1() + 1 != u_turn[ 1 ] ) continue;
			if ( !check_rna_loop( pose, base_pair.res1(), base_pair.res2() - base_pair.res1() + 1 ) ) continue;
			if ( !check_stack( base_pair.res1(), u_turn[1], base_stacks ) ) continue;
			if ( u_turn[3] >= base_pair.res2() ) continue;
			// Now look for 'A' in GNRA:
			for ( Size i = u_turn[3]+1; i < base_pair.res2(); i++ ) {
				if ( !check_stack( u_turn[3],        i, base_stacks ) ) continue;
				if ( !check_stack( i, base_pair.res2(), base_stacks ) ) continue;
				rna_motifs.push_back( RNA_Motif( GNRA_TETRALOOP,
					{ u_turn[1], u_turn[2], u_turn[3], i } ) );
				break;
			}
		}
	}

	/////////////////////////////////////////
	// A-minor motif
	//
	//   1 A o-@ 2 C -o- G 3
	//     A       G -o- C
	//
	//  Note that 1 & 2 are point of WC-Sugar contact.
	//    The other A may be 5' or 3' of the A
	//    The stacked pair may be 5' or 3' of the other stacked pair.
	//
	// Could also set up Type 0, I, II, III ?
	//  [note this supercedes 'ribose zipper']
	//
	//  An alternative would be to use
	//   Grabow et al., 2013 WIRES nomenclature --
	//   A-planar & A-twisted.
	//  Note that they force stereotyped sugar-sugar base pair as anchor:
	//
	//      A > G -o- C
	//      A   X -o- X
	//      5'  5'    3'
	//   Not implemented yet.
	/////////////////////////////////////////
	for ( Size i = 1; i < pose.size(); i++ ) {
		if ( pose.residue( i ).aa() != na_rad &&  pose.residue( i ).na_analogue() != na_rad ) continue;
		if ( pose.residue( i+1 ).aa() != na_rad && pose.residue( i+1 ).na_analogue() != na_rad ) continue;
		if ( !check_rna_loop( pose, i, 2 ) ) continue;
		if ( !check_stack( i, i+1, base_stacks ) ) continue;
		bool found_a_minor( false );
		for ( Size a_offset = 0; a_offset <= 1; a_offset++ ) {
			for ( auto const & base_pair : base_pairs ) {
				if ( base_pair.res1() != i + a_offset ) continue;
				// if ( base_pair.edge1() != WATSON_CRICK ) continue; // should be SUGAR in Grabow.
				if ( base_pair.edge2() != SUGAR ) continue;
				for ( auto const & stacked_pair : rna_motifs.get_motifs( WC_STACKED_PAIR ) ) {
					for ( Size stacked_pair_offset = 0; stacked_pair_offset <= 1; stacked_pair_offset++ ) {
						if ( stacked_pair[ 1+stacked_pair_offset ] == base_pair.res2() ) {
							// utility::vector1< Size > a_minor_res = a_offset ? make_vector1(i+1,i) : make_vector1(i+1,i);
							// if ( stacked_pair_offset ) {
							//  a_minor_res.append( {stacked_pair[2],stacked_pair[1],stacked_pair[4],stacked_pair[3]} );
							// } else {
							//  a_minor_res.append( stacked_pair.residues() );
							// }
							utility::vector1< Size > a_minor_res = { i + a_offset };
							if ( stacked_pair_offset ) {
								a_minor_res.append( {stacked_pair[2],stacked_pair[3]} );
							} else {
								a_minor_res.append( {stacked_pair[1],stacked_pair[4]} );
							}
							rna_motifs.push_back( RNA_Motif( A_MINOR, a_minor_res ) );
							found_a_minor = true; break;
						}
					}
					if ( found_a_minor ) break;
				} // stacked_pair
				if ( found_a_minor ) break;
			} // base_pair
			if ( found_a_minor ) {
				if ( a_offset ) i += 1; // already covered i,i+1
				break;
			}
		} // a_offset
	}

	/////////////////////////////////////////
	// Platform
	//   1     2
	//   A <-@ A
	//   U -o- G 3
	//   4
	//
	//  1 -> 2 -> 3 should be sequence contiguous!
	//
	// See Cate et al., P4-P6 RNA
	//
	/////////////////////////////////////////
	// "platform"
	for ( auto const & base_pair : base_pairs ) {
		// sequence adjacent!
		if ( base_pair.res1()+1 != base_pair.res2() /* canonical */ &&
				base_pair.res1()+2 != base_pair.res2() /* C7.2 */ ) continue;
		if ( !check_rna_loop( pose, base_pair.res2(), base_pair.res2() - base_pair.res1() + 1 ) ) continue;
		if ( base_pair.edge1() != SUGAR ) continue;
		if ( base_pair.edge2() != HOOGSTEEN /* canonical */ &&  base_pair.edge2() != WATSON_CRICK /* C7.2 */ ) continue;
		for ( auto const & base_pair2 : base_pairs ) {
			if ( base_pair2.res1() != base_pair.res2() + 1 ) continue; // sequence adjacent
			if ( !check_rna_loop( pose, base_pair.res2(), 2 ) ) continue;
			if ( !check_stack( base_pair.res1(), base_pair2.res2(), base_stacks ) ) continue;
			if ( !check_stack( base_pair.res2(), base_pair2.res1(), base_stacks ) ) continue;
			// perhaps could relax this -- then we could include bulged-G motif
			// if ( base_pair2.edge1() != WATSON_CRICK ) continue;
			// if ( base_pair2.edge2() != WATSON_CRICK ) continue;
			rna_motifs.push_back( RNA_Motif( PLATFORM,
				{ base_pair.res1(), base_pair.res2(),
				base_pair2.res1(), base_pair2.res2() } ) );
		}
	}


	/////////////////////////////////////////
	// Tetraloop receptor
	//
	//  classic 11-nt:
	//    1 C -o- G 10 <-- stacked pair
	//    2 C -o- G 9 -[U] <-- bulge is optional.
	//    3 U o-@ A 8 /
	//    4 A  <-@  A 5   <-- PLATFORM! Note strand crossover
	//    7 U  -o-  G 6
	//
	//  R(1):
	//    1 C -o- G 10 <-- stacked pair
	//    2 C -o- G 9
	//    3 U -o- U 8  <-- does not have to be UA handle.
	//    4 G  <-@  U 5   <-- PLATFORM! Note strand crossover
	//    7 A  -o-  G 6
	//
	// Based on similarities across
	// models of C7.2, C7.10, R1 &
	// classic 11-nt tetraloop receptors
	//
	/////////////////////////////////////////
	for ( auto const & platform : rna_motifs.get_motifs( PLATFORM ) ) {
		for ( auto const & base_pair : base_pairs ) {
			if ( base_pair.res1() == platform[4] ) continue;
			if ( base_pair.res2() == platform[3] ) continue;

			// sequence adjacent on 5' strand
			if ( base_pair.res1() >= platform[1] ) continue;
			if ( platform[1] - base_pair.res1() > 3 ) continue;

			// sequence adjacent on 3' strand
			if ( base_pair.res2() <= platform[4] ) continue;
			if ( base_pair.res2() - platform[4] > 3 ) continue;

			if ( !check_stack( base_pair.res1(), platform[1], base_stacks ) ) continue;

			for ( auto const & stacked_pair : rna_motifs.get_motifs( WC_STACKED_PAIR ) ) {

				// sequence adjacent on 5' strand
				if ( stacked_pair[2] >= base_pair.res1() ) continue;
				if ( base_pair.res1() - stacked_pair[2] > 3 ) continue;

				// sequence adjacent on 3' strand
				if ( stacked_pair[3] <= base_pair.res2() ) continue;
				if ( stacked_pair[3] - base_pair.res2() > 3 ) continue;

				if ( !check_stack( stacked_pair[2], base_pair.res1(), base_stacks ) &&
						!check_stack( stacked_pair[3], base_pair.res2(), base_stacks ) ) continue;
				rna_motifs.push_back( RNA_Motif( TL_RECEPTOR,
					{ stacked_pair[1], stacked_pair[2], base_pair.res1(), platform[1], platform[2], platform[3], platform[4], base_pair.res2(), stacked_pair[3], stacked_pair[4] } ) );
			}
		}
	}



	///////////////////////////////////////////
	// Tetraloop/tetraloop-receptor (docked!)
	//
	//    1 C -o- G 10  _ 4A  G 1
	//    2 C -o- G 9 _/  3A /
	//    3 U o-@ A 8     2A
	//    4 A  <-@  A 5
	//    7 U  -o-  G 6
	//
	///////////////////////////////////////////
	for ( auto const & tetraloop : rna_motifs.get_motifs( GNRA_TETRALOOP ) ) {
		for ( auto const & receptor : rna_motifs.get_motifs( TL_RECEPTOR ) ) {
			if ( !check_stack( tetraloop[2], receptor[5], base_stacks ) ) continue;
			for ( auto const & a_minor : rna_motifs.get_motifs( A_MINOR ) ) {
				if ( a_minor[1] != tetraloop[ 4 ] ) continue;
				if ( a_minor[2] != receptor[ 9 ] ) continue;
				utility::vector1< Size > res( tetraloop.residues() );
				res.append( receptor.residues() );
				rna_motifs.push_back( RNA_Motif( TETRALOOP_TL_RECEPTOR, res ) );
			}
		}
	}


	/////////////////////////////////////////////////////////////////////////////////
	// loopE-submotif
	//
	// 3  A @-o U 4
	// 2  G <-@ A 5
	// 1  G -o- C 6
	//
	//  occurs in loopE, SRL, SRP -- note that
	//  sequence can shift in, e.g., SRP
	//
	// See, e.g., Leontis & Westhof, Comp Funct Genom 2002; 3: 518–524.
	//  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2448414/pdf/CFG-03-518.pdf
	//
	// Specific case of "cross-strand purine stack" from Moore, Ann. Rev. Biochem. 1999.
	// Also called GA_UA motif [see, e.g., Grabow et al., WIREs RNA 2013]
	//
	/////////////////////////////////////////////////////////////////////////////////
	for ( auto const & base_pair1 : base_pairs ) {
		if ( base_pair1.edge1() != WATSON_CRICK ) continue;
		if ( base_pair1.edge2() != WATSON_CRICK ) continue;
		if ( !check_rna_loop( pose, base_pair1.res1(),   3 ) ) continue;
		if ( !check_rna_loop( pose, base_pair1.res2()-2, 3 ) ) continue;
		for ( auto const & base_pair2 : base_pairs ) {
			// direct sequence adjacent
			if ( base_pair2.res1() != base_pair1.res1()+1 ) continue;
			if ( base_pair1.res2() != base_pair2.res2()+1 ) continue;
			// 'sheared' base pair on top of W/C pair.
			if ( base_pair2.edge1() != SUGAR &&
					potential.get_base_backbone( pose.residue( base_pair2.res2() ),
					pose.residue( base_pair2.res1() ), 6 /* O2' */) >= 0.0  /* for 1lnt edge case*/ )  continue;
			if ( base_pair2.edge2() != HOOGSTEEN ) continue;
			// directly stacking of 'sheared' base pair on W/C
			if ( !check_stack( base_pair1.res1(), base_pair2.res1(), base_stacks ) ) continue;
			if ( !check_stack( base_pair1.res2(), base_pair2.res2(), base_stacks ) ) continue;
			for ( auto const & base_pair3 : base_pairs ) {
				if ( base_pair3.edge1() != HOOGSTEEN ) continue;
				if ( base_pair3.edge2() != WATSON_CRICK ) continue;
				// direct sequence adjacent
				if ( base_pair3.res1() != base_pair2.res1()+1 ) continue;
				if ( base_pair2.res2() != base_pair3.res2()+1 ) continue;
				// actually, these two will not stack if in bulged-G
				// if ( !check_stack( base_pair2.res1(), base_pair3.res1(), base_stacks ) ) continue;
				// this is the "cross-strand" purine stack:
				if ( !check_stack( base_pair2.res2(), base_pair3.res1(), base_stacks ) ) continue;
				rna_motifs.push_back( RNA_Motif( LOOP_E_SUBMOTIF,
					{base_pair1.res1(), base_pair2.res1(), base_pair3.res1(),
					base_pair3.res2(), base_pair2.res2(), base_pair1.res2() } ) );
			}
		}
	}

	/////////////////////////////////////////////////////
	// bulged-G motif
	//           5     4
	// 3  A @-o  U @-> G
	// 2  G <-@  A 6
	// 1  G -o-  C 7
	//
	//  bulged-G motif from, e.g., Sarcin-Ricin loop
	//
	// Note: should we consider the G4-U5 a 'platform'?
	// TODO: These motifs usually have an additional
	//         N and A 5' of G4 with a cool turn -- include here?
	//         Or define a more extended motif to include the S-turn?
	//
	/////////////////////////////////////////////////////
	for ( auto const & loop_e_submotif : rna_motifs.get_motifs( LOOP_E_SUBMOTIF ) ) {
		for ( auto const & base_pair : base_pairs ) {
			if ( base_pair.res2() == loop_e_submotif[ 4 ] ) {
				if ( base_pair.edge1() != SUGAR ) continue; /*bulged G*/
				if ( base_pair.edge2() != HOOGSTEEN ) continue;
				// "platform"
				if ( base_pair.res1()+1 == base_pair.res2() /* SRL*/ ||
						base_pair.res1()+2 == base_pair.res2() /* FMN aptamer, PDB 1FMN*/ ) {
					if ( !check_rna_loop( pose, base_pair.res1(), base_pair.res2() - base_pair.res1() ) ) continue;
					// check that bulged G makes contact with phosphate
					if ( potential.get_base_backbone( pose.residue(base_pair.res1()), pose.residue(loop_e_submotif[3]), 1 /* O1P */) >= 0.0 &&
							potential.get_base_backbone( pose.residue(base_pair.res1()), pose.residue(loop_e_submotif[3]), 2 /* O2P */) >= 0.0 ) continue;
					rna_motifs.push_back( RNA_Motif( BULGED_G,
						{ loop_e_submotif[1], loop_e_submotif[2], loop_e_submotif[3],
						base_pair.res1(),
						loop_e_submotif[4], loop_e_submotif[5], loop_e_submotif[6] } ) );
				}
			}
		}
	}


	/////////////////////////////////////////
	// Z-turn
	//  see Auffinger
	/////////////////////////////////////////

	/////////////////////////////////////////
	// U-turn loop
	//  generalization of GNRA
	//  see Auffinger
	/////////////////////////////////////////

	/////////////////////////////////////////
	// G-quartet
	/////////////////////////////////////////

	/////////////////////////////////////////
	// G-quartet-stack
	/////////////////////////////////////////

	/////////////////////////////////////////
	// Tandem G/A (sheared)
	// a la Cruz/Westhof.
	// or cross-strand purine stack [Moore, 1999]
	/////////////////////////////////////////

	/////////////////////////////////////////
	// Kink-turn
	// contains tandem G/A
	/////////////////////////////////////////

	/////////////////////////////////////////
	// Tandem W/C
	//  A -o- A
	//  G -o- A
	// should this include 'normal' base pair
	//  steps?
	// may be needed to counter-balance
	//  sheared G/A
	// Happens in 1FMN (FMN aptamer)
	/////////////////////////////////////////

	/////////////////////////////////////////
	// what is P4-P6 j5/5a? --
	//  W/C stacked pair with minor groove
	//  docked into non-W/C stacked-pair
	// May be describable as an 'off-label'
	// UA_h turn, since it has a UA_handle
	//  made with a CC instead of a UA.
	/////////////////////////////////////////

	/////////////////////////////////////////
	// G-ribo
	/////////////////////////////////////////

	/////////////////////////////////////////
	// Double T-loop
	//  a la T-box, RNase P
	/////////////////////////////////////////

	/////////////////////////////////////////
	// A-minor junctions
	//  a la Jaeger
	/////////////////////////////////////////

	/////////////////////////////////////////
	// T-loop/PK
	//  a la Jaeger
	/////////////////////////////////////////

	return rna_motifs;
}

/// @brief wrapper around get_rna_motifs -- outputs score.
inline
Real
get_rna_motif_score(
	pose::Pose const & pose,
	core::scoring::rna::RNA_LowResolutionPotential const & potential,
	core::pose::rna::RNA_FilteredBaseBaseInfo const & filtered_base_base_info )
{
	RNA_Motifs const rna_motifs = get_rna_motifs( pose, potential, filtered_base_base_info );
	Real score( 0.0 );
	for ( auto const & rna_motif : rna_motifs.get_motifs() ) {
		if ( rna_motif_bonus.count( rna_motif.type() ) ) continue;
		score += rna_motif_bonus.find( rna_motif.type() )->second;
	}
	return score;
}

// @brief output for RNA_Motifs
void
output_rna_motifs(
	std::ostream & out,
	core::pose::Pose const & pose,
	RNA_Motifs const & motifs,
	bool const output_WC_stacked_pair = false);

/// @brief sets up .pml file that will color motifs.
void
output_motifs_to_pymol(
	std::ostream & out,
	core::pose::Pose const & pose,
	RNA_Motifs const & rna_motifs );


} //rna
} //scoring
} //core

#endif
