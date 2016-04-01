// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/EnvPairPotential.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/EnvPairPotential.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/prof.hh>
// Utility headers
#include <utility/io/izstream.hh>

#include <utility/vector1.hh>


// just for debugging
//#include <ObjexxFCL/format.hh>

// C++


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

/// @details Copy constructors must copy all data, not just some...
CenListInfo::CenListInfo( CenListInfo const & src ) :
	CacheableData()
{
	fcen6_ = src.fcen6_;
	fcen10_ = src.fcen10_;
	fcen12_ = src.fcen12_;
	calculated_ = src.calculated_;
}

void
CenListInfo::initialize( pose::Pose const & pose )
{
	Size const nres( pose.total_residue() );

	fcen6_.resize( nres, 0.0 );
	fcen10_.resize( nres, 0.0 );
	fcen12_.resize( nres, 0.0 );

	std::fill( fcen6_.begin(), fcen6_.end(), 1.0 );   // 1 because a residue is w/i 6 A of itself
	std::fill( fcen12_.begin(), fcen12_.end(), 0.0 ); // 0 because a residue is not between 6 and 12 A of itself
	std::fill( fcen10_.begin(), fcen10_.end(), 1.0 ); // 1 because a residue is w/i 10A of itself
}

EnvPairPotential::EnvPairPotential():
	cen_dist_cutoff2( 144.0 ),

	//cems transition regions between environment bins
	//cems transition is from +/- sqrt(36+pad6) +/- sqrt(100+pad10) etc
	cen_dist5_pad( 0.5 ),
	cen_dist6_pad( 0.6 ),
	cen_dist7_pad( 0.65 ),
	cen_dist10_pad( 1.0 ),
	cen_dist12_pad( 1.2 ),

	cen_dist5_pad_plus ( cen_dist5_pad  + 25.0 ),
	cen_dist6_pad_plus( cen_dist6_pad + 36.0 ),
	cen_dist7_pad_plus ( cen_dist7_pad  + 56.25 ),
	cen_dist10_pad_plus( cen_dist10_pad + 100.0 ),
	cen_dist12_pad_plus( cen_dist12_pad + 144.0 ),

	cen_dist5_pad_minus ( cen_dist5_pad  - 25.0 ),
	cen_dist7_pad_minus ( cen_dist7_pad  - 56.25 ),
	cen_dist10_pad_minus( cen_dist10_pad - 100.0 ),
	cen_dist12_pad_minus( cen_dist12_pad - 144.0 ),

	cen_dist5_pad_hinv ( 0.5 / cen_dist5_pad ),
	cen_dist6_pad_hinv ( 0.5 / cen_dist6_pad ),
	cen_dist7_pad_hinv ( 0.5 / cen_dist7_pad ),
	cen_dist10_pad_hinv( 0.5 / cen_dist10_pad ),
	cen_dist12_pad_hinv( 0.5 / cen_dist12_pad ),

	cen_dist_cutoff_12_pad( cen_dist_cutoff2 + cen_dist12_pad )
{
	// load the data
	Size const max_aa( 20 ); // just the standard aa's for now
	Size const env_log_table_size( 40 );
	Size const pair_log_table_size( 5 );
	Size const cbeta_den_table_size( 45 );
	Size const cenpack_log_table_size( 120 );

	std::string tag,line;
	chemical::AA aa;

	{ // env_log
		env_log_.dimension( max_aa, env_log_table_size );

		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/EnvPairPotential/env_log.txt" );
		while ( getline( stream, line ) ) {
			std::istringstream l(line);
			l >> tag >> aa;
			for ( Size i=1; i<= env_log_table_size; ++i ) {
				l >> env_log_(aa,i);
			}
			if ( l.fail() || tag != "ENV_LOG:"  ) utility_exit_with_message("bad format for scoring/score_functions/EnvPairPotential/env_log.txt");
		}
	}

	{ // cebeta_den_6/12
		cbeta_den6_.dimension( cbeta_den_table_size );
		cbeta_den12_.dimension( cbeta_den_table_size );

		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/EnvPairPotential/cbeta_den.txt" );

		{ // den6
			getline( stream, line );
			std::istringstream l(line);
			l >> tag;
			for ( Size i=1; i<= cbeta_den_table_size; ++i ) {
				l >> cbeta_den6_(i);
			}
			if ( l.fail() || tag != "CBETA_DEN6:"  ) utility_exit_with_message("bad format for scoring/score_functions/EnvPairPotential/cbeta_den.txt");
		}

		{ // den12
			getline( stream, line );
			std::istringstream l(line);
			l >> tag;
			for ( Size i=1; i<= cbeta_den_table_size; ++i ) {
				l >> cbeta_den12_(i);
			}
			if ( l.fail() || tag != "CBETA_DEN12:"  ) utility_exit_with_message("bad format for scoring/score_functions/EnvPairPotential/cbeta_den.txt");
		}
	}

	{ // pair_log
		pair_log_.dimension( pair_log_table_size, max_aa, max_aa );

		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/EnvPairPotential/pair_log.txt" );
		for ( Size j=1; j<= pair_log_table_size; ++j ) {
			for ( Size k=1; k<= max_aa; ++k ) {
				getline( stream, line );
				std::istringstream l(line);
				Size jj;
				l >> tag >> jj >> aa;
				debug_assert( Size(aa) == k );
				for ( Size i=1; i<= max_aa; ++i ) {
					l >> pair_log_(j,aa,i);
				}
				if ( l.fail() || jj != j || tag != "PAIR_LOG:"  ) utility_exit_with_message("bad format for scoring/score_functions/EnvPairPotential/pair_log.txt");
			}
		}
	}

	{ // cenpack_log
		cenpack_log_.dimension( cenpack_log_table_size ); //sequence independent

		utility::io::izstream stream;
		basic::database::open( stream, "scoring/score_functions/EnvPairPotential/cenpack_log.txt" );
		for ( Size j=1; j<= cenpack_log_table_size; ++j ) {
			getline( stream, line );
			std::istringstream l(line);
			Size jj;
			l >> tag >> jj;
			l >> cenpack_log_(j);
			if ( l.fail() || jj != j || tag != "CENPACK_LOG:"  ) utility_exit_with_message("bad format for scoring/score_functions/EnvPairPotential/cenpack_log.txt");
		}
	}
}

/// @brief fill the cenlist using interpolation
/// @details
///cems--------------------------------------------------------------------------
///  interpolation notes --Historically we have broken the
///  centroid density statistics into three bins: i) pairs
///  less than 6 angstroms ii) pairs less than 10 angstroms ems
///  iii) and pairs between 6 and 12 angstroms the resulting
///  abruptness in the scoring functions due to the arbitrary radius
///  cutoffs has caused some problems during gradient minimization.
///  therefore this was replaced with an interpolated binning
///  schema as follows: When a pairwise distance lies within "+/-
///  dr" of the bin boundary (6,10,12) then partial credit is given
///  to the enclosing bins.  So for example, if fgap=0.5 angstroms, and
///  a pair radius were 6.4 angstroms, then a fractional count is
///  given to BOTH the "less-than-6" bin AND to the
///  "between-6-and-10" bin.  The sum of these fractions always adds to
///  one.  So that we dont have to re-do the statistics we
///  currently use we want to keep "fgap" small.  ideally fgap
///  should be large compared to the search algorithm step size, and
///  larger than the expected roundoff error in any refold
///  operation, and otherwise as small as possible.  Also we want
///  to cleverly choose the interpolation function so that the average
///  number of counts getting into the bins is the same as under
///  the old schema.  As long as dr is small then we can use either
///  r+/-fgap or alternatively r^2+/-fgap^2 and this will be
///  approximately satsified.  since the squared from is easier to work
///  we will use this.  in the code below the frag^2 term is called
///  a _pad, and we allow for different pad_sizes on the three radii.
///cems--------------------------------------------------------------------------
void
EnvPairPotential::fill_cenlist(
	CenListInfo & cenlist,
	Size const res1,
	Size const res2,
	Real const cendist
) const {
	debug_assert( cendist <= cen_dist12_pad_plus );

	//  compute arrays needed for C-beta  energy function
	Real const one( 1.0 );

	// NOTE: *_hinv is negative of vdw.cc version and positive in structure.cc version.
	// We are using postive hinvs.
	if ( cendist <= cen_dist10_pad_plus ) {
		Real interp = std::min( ( cen_dist10_pad_plus - cendist ) * cen_dist10_pad_hinv, one );
		cenlist.fcen10(res1) += interp;
		cenlist.fcen10(res2) += interp;
	}

	if ( cendist <= cen_dist6_pad_plus ) { // its sort of a "6" and not so much a "12"
		Real interp = std::min( ( cen_dist6_pad_plus - cendist ) * cen_dist6_pad_hinv, one );

		cenlist.fcen6(res1) += interp;
		cenlist.fcen6(res2) += interp;
		cenlist.fcen12(res1) += 1.0 - interp;
		cenlist.fcen12(res2) += 1.0 - interp;

	} else { // then its sort of a "12" but definitely not a "6"

		Real interp = std::min( ( cen_dist12_pad_plus - cendist ) * cen_dist12_pad_hinv, one );

		cenlist.fcen12(res1) += interp;
		cenlist.fcen12(res2) += interp;
	}
}

void
EnvPairPotential::truncate_cenlist_values( CenListInfo & cenlist ) const
{
	for ( Size ii = 1; ii <= cenlist.size(); ++ii ) {
		if ( cenlist.fcen6(ii) >= 45.0 ) cenlist.fcen6(ii) = 44.9999;
		if ( cenlist.fcen10(ii) >= 31.0 ) cenlist.fcen10(ii) = 30.9999;
		if ( cenlist.fcen12(ii) < 1 ) {
			cenlist.fcen12(ii) = 1;
		} else if ( cenlist.fcen12(ii) >= 45.0 ) {
			cenlist.fcen12(ii) = 44.9999;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////
void
EnvPairPotential::compute_centroid_environment(
	pose::Pose & pose
) const
{
	// basic::ProfileThis doit( basic::ENERGY_ENVPAIR_POTENTIAL );

	CenListInfo & cenlist( nonconst_cenlist_from_pose( pose ));

	/// Energy graph contains edges for all residue pairs with
	/// centroids w/i cen_dist_cutoff_12_pad
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	Size const nres( energy_graph.num_nodes() );

	/// calculate the cenlist info only if it has not been calculated since the last score evaluation
	if ( cenlist.calculated() ) return;

	// ensure that cenlist has pose.total_residue() elements in case the pose has
	// changed its sequence lenght since the last cenlist update
	cenlist.initialize( pose );
	
	for ( Size i = 1; i < nres; ++i ) {
		conformation::Residue const & rsd1 ( pose.residue(i) );
		if ( !rsd1.is_protein() ) continue;
		for ( graph::Graph::EdgeListConstIter
			 iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
			 irue = energy_graph.get_node(i)->const_upper_edge_list_end();
			 iru != irue; ++iru ) {
			EnergyEdge const * edge( static_cast< EnergyEdge const *> (*iru) );
			Size const j( edge->get_second_node_ind() );
			conformation::Residue const & rsd2 ( pose.residue(j) );
			if ( !rsd2.is_protein() ) continue;
			
			Real const cendist = edge->square_distance();
			
			//  compute arrays needed for C-beta  energy function
			//  first do a coarse grain reality check on centroid separations
			if ( cendist <= cen_dist_cutoff_12_pad ) {
				fill_cenlist( cenlist, i, j, cendist );
			}
		}
	}
	
	truncate_cenlist_values( cenlist );
	cenlist.calculated() = true;
}

void
EnvPairPotential::finalize( pose::Pose & pose ) const
{
	CenListInfo & cenlist( nonconst_cenlist_from_pose( pose ));
	cenlist.calculated() = false;
}

////////////////////////////////////////////////////////////////////////////////////
void
EnvPairPotential::evaluate_env_and_cbeta_scores(
	pose::Pose const & pose,
	conformation::Residue const & rsd,
	Real & env_score,
	Real & cb_score6,
	Real & cb_score12
) const {
	//using ObjexxFCL::format::F; // debugging
	//using ObjexxFCL::format::I;
	// basic::ProfileThis doit( basic::ENERGY_ENVPAIR_POTENTIAL );

	CenListInfo const & cenlist( cenlist_from_pose( pose ));

	int const position ( rsd.seqpos() );

	Real const fcen6  ( cenlist.fcen6( position) );
	Real const fcen10 ( cenlist.fcen10(position) );
	Real const fcen12 ( cenlist.fcen12(position) );

	if ( ! rsd.is_protein() ) { // amino acid check
		env_score = 0.0;
		cb_score6  = 0.0;
		cb_score12 = 0.0;
	} else {

		env_score = env_log_( rsd.aa(), static_cast< int >( fcen10 ) );

		// interp1 rounds down to nearest (non-negative) integer.
		int interp1 = static_cast< int >( fcen6 );
		// note cen6 is always at least 1.0

		// fraction remainder after nearest lower integer is removed
		Real interp2 = fcen6 - interp1;

		//  use interp2 to linearly interpolate the two nearest bin values
		cb_score6 =
			( ( 1.0 - interp2 ) * cbeta_den6_( interp1 ) +
			(   interp2 ) * cbeta_den6_( interp1+1 ) );

		interp1 = static_cast< int >( fcen12 );
		// note cen12 is always at least 1.0 -- this is in fact false for fcen12
		interp2 = fcen12 - interp1;
		cb_score12 =
			( ( 1.0 - interp2 ) * cbeta_den12_( interp1   ) +
			(   interp2 ) * cbeta_den12_( interp1+1 ) );

		//std::cout << "eval_env_cbeta: " << I(4,rsd.seqpos()) << F(9,3,fcen6) << F(9,3,fcen10) << F(9,3,fcen12) <<
		// F(9,3,env_score) << F(9,3,cb_score6) << F(9,3,cb_score12) << ' ' << rsd.name() << std::endl;
		//std::cout << "fcen6( " << position << " ) = " << fcen6 << " fcen10( " <<  position << " ) " << fcen10 << " fcen12( " << position << " ) = ";
		//std::cout << fcen12 << " "; //<< std::endl;
		// " interp1: " << interp1 << " interp2: " << interp2 << std::endl;
	}
}


///////////////////////////////////////////////////////////////////////////////////////////////

void
EnvPairPotential::evaluate_pair_and_cenpack_score(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Real const cendist,
	Real & pair_contribution,
	Real & cenpack_contribution
) const {
	// basic::ProfileThis doit( basic::ENERGY_ENVPAIR_POTENTIAL );

	pair_contribution = 0.0;
	cenpack_contribution = 0.0;

	if ( !rsd1.is_protein() || !rsd2.is_protein() ) return;

	chemical::AA const aa1( rsd1.aa() );
	chemical::AA const aa2( rsd2.aa() );

	//CAR  no pair score if a disulfide
	if ( rsd1.type().is_disulfide_bonded()
			&& rsd2.type().is_disulfide_bonded()
			&& rsd1.is_bonded( rsd2 )
			&& rsd1.polymeric_sequence_distance( rsd2 ) > 1
			&& rsd1.has_variant_type( chemical::DISULFIDE )
			&& rsd2.has_variant_type( chemical::DISULFIDE ) ) return;

	// no pair score for residues closer than 9 in sequence
	if ( rsd1.polymeric_sequence_distance( rsd2 ) /* j - i */ <= 8 ) return;

	//$$$  we now try to find which bin the pair distance lies in
	//$$$  I note this could in principle be calculated and updatded
	//$$$  just like cen_dist is if there is a need for speed.
	//$$$  this function interpolates between bins.
	//$$$  An important(!) requirement on pair_log is that the
	//$$$  value should approach zero as the radius increases.
	//$$$  this fact permits us not to have to compute and score pairs are larger
	//$$$  than cen_dist > cutoff.

	int icon = 5;
	Real interp2( 0.0 );

	if ( cendist > cen_dist10_pad_plus ) {
		icon = 4;
		interp2 = ( cendist + cen_dist12_pad_minus ) * cen_dist12_pad_hinv;
	} else {
		if ( cendist > cen_dist7_pad_plus ) {
			icon = 3;
			interp2 = ( cendist + cen_dist10_pad_minus ) * cen_dist10_pad_hinv;
		} else {
			if ( cendist > cen_dist5_pad_plus ) {
				icon = 2;
				interp2 = ( cendist + cen_dist7_pad_minus ) * cen_dist7_pad_hinv;
			} else {
				icon = 1;
				interp2 = ( cendist + cen_dist5_pad_minus ) * cen_dist5_pad_hinv;
			}
		}
	}
	if ( interp2 < 0.0 ) interp2 = 0.0;

	// note in theory this will never happen but in practice round off
	// error can cause problem
	if ( interp2 > 1.0 ) interp2 = 1.0;

	// handle last bin specially since icon+1 would be past array end
	// pb note -- I don't think icon will ever be 5 here, wonder if it has always been this way?
	if ( icon != 5 ) {
		pair_contribution =
			( ( 1.0f - interp2 ) * pair_log_( icon  , aa1, aa2 ) +
			(     interp2 ) * pair_log_( icon+1, aa1, aa2 ) );
	} else {
		pair_contribution =   ( 1.0f - interp2 ) * pair_log_( icon  , aa1, aa2 );
	}


	// Adding a term that should help reproduce pairwise correlation function between centroids
	//   as observed in the PDB.
	int cendist_bin = static_cast <int> ( sqrt( cendist ) * 10 + 1); //Binned with 0.1 A width.

	if ( cendist_bin > 120 )   cendist_bin = 120;
	if ( cendist_bin <   1 )   cendist_bin = 1;

	cenpack_contribution = cenpack_log_( cendist_bin );
}

/// @details Pose must already contain a cenlist object or this method will fail.
CenListInfo const &
EnvPairPotential::cenlist_from_pose( pose::Pose const & pose )
{
	using namespace core::pose::datacache;
	return static_cast< core::scoring::CenListInfo const & > ( pose.data().get( CacheableDataType::CEN_LIST_INFO ) );
}

/// @details Either returns a non-const reference to the cenlist object already stored
/// in the pose, or creates a new cenist object, places it in the pose, and returns
/// a non-const reference to it.
CenListInfo &
EnvPairPotential::nonconst_cenlist_from_pose( pose::Pose & pose )
{
	if ( pose.data().has( core::pose::datacache::CacheableDataType::CEN_LIST_INFO ) ) {
		return static_cast< core::scoring::CenListInfo & > ( pose.data().get( core::pose::datacache::CacheableDataType::CEN_LIST_INFO ));
	}
	// else
	CenListInfoOP cenlist( new CenListInfo );
	pose.data().set( core::pose::datacache::CacheableDataType::CEN_LIST_INFO, cenlist );
	return *cenlist;
}


}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::CenListInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( fcen6_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( fcen10_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( fcen12_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( calculated_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::CenListInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( fcen6_ ); // utility::vector1<Real>
	arc( fcen10_ ); // utility::vector1<Real>
	arc( fcen12_ ); // utility::vector1<Real>
	arc( calculated_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::CenListInfo );
CEREAL_REGISTER_TYPE( core::scoring::CenListInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_EnvPairPotential )
#endif // SERIALIZATION
