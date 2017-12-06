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

enum RNA_MotifType {
	U_TURN,
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
	TETRALOOP_TL_RECEPTOR,
	Z_TURN,
	TANDEM_GA_SHEARED,
	TANDEM_GA_WATSON_CRICK,
	GA_MINOR,
	DOUBLE_T_LOOP
};

std::map< RNA_MotifType, core::Real > const rna_motif_bonus = {
{ U_TURN, -1.0 },
{ UA_HANDLE, -1.0 },
{ T_LOOP, -1.0 },
{ INTERCALATED_T_LOOP, -1.0 },
{ LOOP_E_SUBMOTIF, -1.0 },
{ BULGED_G, -1.0 },
{ GNRA_TETRALOOP, -0.1 },  // GNRA is already pretty favorable.
{ STRICT_WC_STACKED_PAIR, 0.0 },
{ WC_STACKED_PAIR, 0.0 },
{ A_MINOR, -1.0 },
{ PLATFORM, -1.0 },
{ TL_RECEPTOR, -1.0 },
{ TETRALOOP_TL_RECEPTOR, -1.0 },
{ Z_TURN, -1.0 },
{ TANDEM_GA_SHEARED, -1.0 },
{ TANDEM_GA_WATSON_CRICK, -1.0 },
{ GA_MINOR, -1.0 },
{ DOUBLE_T_LOOP, -1.0 }
};


std::map< RNA_MotifType, std::string > const motif_color = {
{ U_TURN, "lightblue" },
{ UA_HANDLE, "marine" },
{ T_LOOP, "tv_blue" },
{ INTERCALATED_T_LOOP, "deepblue" },
{ LOOP_E_SUBMOTIF, "salmon" },
{ BULGED_G, "red" },
{ GNRA_TETRALOOP, "ruby" },
{ STRICT_WC_STACKED_PAIR, "gray20" },
{ WC_STACKED_PAIR, "gray50" },
{ A_MINOR, "gold" },
{ PLATFORM, "sand" },
{ TL_RECEPTOR, "limon" },
{ TETRALOOP_TL_RECEPTOR, "orange" },
{ Z_TURN, "violet" },
{ TANDEM_GA_SHEARED, "lime" },
{ TANDEM_GA_WATSON_CRICK, "lime" },
{ GA_MINOR, "purple" },
{ DOUBLE_T_LOOP, "forest" }
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
	case Z_TURN : return "Z_TURN";
	case TANDEM_GA_SHEARED : return "TANDEM_GA_SHEARED";
	case TANDEM_GA_WATSON_CRICK : return "TANDEM_GA_WATSON_CRICK";
	case GA_MINOR : return "GA_MINOR";
	case DOUBLE_T_LOOP : return "DOUBLE_T_LOOP";
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

	void set_type( RNA_MotifType const & setting ) { type_ = setting; }
	RNA_MotifType type() const { return type_; }

	void set_residues( utility::vector1< core::Size > const & setting ) { residues_ = setting; }
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
	core::Size const n,
	core::Size const nloop )
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
check_stack( Size const i, Size const j,
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
	std::cout << ObjexxFCL::right_string_of( to_string( motif.type() ), 20 ) << ": ";
	for ( auto const & res : motif ) {
		std::cout << ' ' << ObjexxFCL::right_string_of( pose.residue( res ).annotated_name(), 6 ) <<
			"-" << pose.pdb_info()->chain( res ) << ":";
		if ( pose.pdb_info()->segmentID( res ).size() > 0 &&
				pose.pdb_info()->segmentID( res ) != "    " ) std::cout << pose.pdb_info()->segmentID( res ) << ":";
		std::cout << ObjexxFCL::left_string_of( pose.pdb_info()->number( res ), 4 );
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
check_watson_crick_sequence( pose::Pose const & pose, Size const res1, Size const res2,
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
///          Used to be inlined -- since this function is long, this is probably
///          to no one's advantage.
RNA_Motifs
get_rna_motifs( pose::Pose const & pose,
	core::scoring::rna::RNA_LowResolutionPotential const & potential,
	core::pose::rna::RNA_FilteredBaseBaseInfo const & filtered_base_base_info );

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
	bool const output_WC_stacked_pair = false );

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
