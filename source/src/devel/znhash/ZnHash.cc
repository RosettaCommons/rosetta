// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/znhash/ZnHash.cc
/// @brief  Implementation of zinc-match hash for use in optimizing zinc coordination
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Bryan Der (bder@email.unc.edu)

// Unit headers
#include <devel/znhash/ZnHash.hh>

// Core headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <numeric/polynomial.hh>
#include <core/scoring/hbonds/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>

// Protocols headers
#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/enzdes/EnzdesMovers.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/FixedSizeLexicographicalIterator.hh>
#include <utility/FixedSizeLexicographicalIterator.tmpl.hh>


namespace devel {
namespace znhash {

/// @details Auto-generated virtual destructor
ZnHash::~ZnHash() {}

static thread_local basic::Tracer TR( "devel.znhash.ZnHash" );

ZnMatchData::ZnMatchData() :
	res1_(0),
	res2_(0)
{}

ZnMatchData::~ZnMatchData() {}

ZnMatchData::ZnMatchData( ZnMatchData const & src ) :
	index_( src.index_ ),
	res1_( src.res1_ ),
	res2_( src.res2_ ),
	zn_and_orbitals_( src.zn_and_orbitals_ ),
	match_pdb_file_( src.match_pdb_file_ ),
	match_cst_file_( src.match_cst_file_ ),
	res1conf_( src.res1conf_ ),
	res2conf_( src.res2conf_ ),
	znconf_( src.znconf_ )
{}

ZnMatchData const &
ZnMatchData::operator = ( ZnMatchData const & src )
{
	if ( this != & src ) {
		index_ = src.index_;
		res1_ = src.res1_;
		res2_ = src.res2_;
		zn_and_orbitals_ = src.zn_and_orbitals_;
		match_pdb_file_ = src.match_pdb_file_;
		match_cst_file_ = src.match_cst_file_;
		res1conf_ = src.res1conf_;
		res2conf_ = src.res2conf_;
		znconf_ = src.znconf_;
	}
	return *this;
}

void ZnMatchData::res1conf( core::conformation::Residue const & r1 )
{
	res1conf_ = r1.clone();
}

void ZnMatchData::res2conf( core::conformation::Residue const & r2 )
{
	res2conf_ = r2.clone();
}

void ZnMatchData::znconf( core::conformation::Residue const & zn )
{
	znconf_ = zn.clone();
}

core::conformation::Residue const &
ZnMatchData::res1conf() const
{
	return *res1conf_;
}

core::conformation::Residue const &
ZnMatchData::res2conf() const
{
	return *res2conf_;
}

core::conformation::Residue const &
ZnMatchData::znconf() const
{
	return *znconf_;
}


ZnHash::ZnHash() :
	grid_size_( 0.5 ),
	inv_grid_size_( 2.0 ),
	ngrid_cells_( /* 0 */ ),
	hash_has_been_built_( false )
{}

/// @brief must be called before build_hash is called
void ZnHash::set_uniform_bin_width( Real width )
{
	assert( ! hash_has_been_built_ );
	grid_size_ = width;
	inv_grid_size_ = 1 / grid_size_;
}

/// @brief First, add all zinc coordinates to the hash, then invoke build_hash.
void ZnHash::add_zn_coordinate( ZnCoord const & zn )
{
	assert( ! hash_has_been_built_ );
	zn_coords_.push_back( zn );
}

// @brief Build the hash after all ZnCoordinates have been added.  Must be called before query_hash().
void ZnHash::build_hash()
{
	//1. compute a bounding box for all zn coords
	//2. determine the number of grid cells in each dimension
	//2b. extend the bounding box to fit tightly over all cells, even if it ends up a little too big
	//3. setup the hash data
	//4. place all the coordinates into the hash.

	for ( std::list< ZnCoord >::const_iterator iter = zn_coords_.begin(), iter_end = zn_coords_.end();
			iter != iter_end; ++iter ) {
		if ( iter == zn_coords_.begin() ) { bb_.set_lower((*iter)[1]); bb_.set_upper((*iter)[1]); }
		for ( Size ii = 1; ii <= ZnCoord::natoms(); ++ii ) bb_.add( (*iter)[ii] );
	}

	core::Vector dim_ranges = bb_.upper() - bb_.lower();
	for ( Size ii = 1; ii <= 3; ++ii ) {
		ngrid_cells_[ ii ] = static_cast<core::Size>( std::ceil( dim_ranges(ii) / (grid_size_) ) );
		if ( ngrid_cells_[ ii ] == 0 ) ngrid_cells_[ ii ] = 1; // imagine there's only one point
	}
	core::Vector new_upper(
		bb_.lower().x() + grid_size_ * ngrid_cells_[ 1 ],
		bb_.lower().y() + grid_size_ * ngrid_cells_[ 2 ],
		bb_.lower().z() + grid_size_ * ngrid_cells_[ 3 ] );
	bb_.set_upper( new_upper );
	bb_ext_.set_lower( bb_.lower() - grid_size_ );
	bb_ext_.set_upper( bb_.upper() + grid_size_ );

	ndim_prods_[3] = 1;
	ndim_prods_[2] = ngrid_cells_[3];
	ndim_prods_[1] = ndim_prods_[2]*ngrid_cells_[2];

	ZnCoordinateHash::iterator hash_end = zn_hash_.end();
	for ( std::list< ZnCoord >::const_iterator iter = zn_coords_.begin(), iter_end = zn_coords_.end();
			iter != iter_end; ++iter ) {
		boost::uint64_t iter_bin = bin_for_point( (*iter)[1] );
		ZnCoordinateHash::iterator hash_iter = zn_hash_.find( iter_bin );
		if ( hash_iter == hash_end ) {
			zn_hash_[ iter_bin ].push_back( *iter );
		} else {
			(*hash_iter).second.push_back( *iter );
		}

	}
	hash_has_been_built_ = true;
}

boost::uint64_t ZnHash::bin_for_point( Vector const & query_point ) const
{
	assert( bb_.contains( query_point ) );
	Vector shifted = query_point - bb_.lower();
	boost::uint64_t sum(0);
	for ( Size ii = 1; ii <= 3; ++ii ) {
		sum += static_cast< core::Size > (shifted(ii) * inv_grid_size_) * ndim_prods_[ii];
	}
	return sum;
}

ZnHash::ZnCoordinateHash::const_iterator ZnHash::query_hash( Vector const & query_point ) const
{
	if ( ! bb_.contains( query_point ) ) return hash_end();
	return zn_hash_.find( bin_for_point( query_point ) );
}

ZnHash::ZnCoordinateHash::const_iterator ZnHash::hash_begin() const
{
	return zn_hash_.begin();
}

ZnHash::ZnCoordinateHash::const_iterator ZnHash::hash_end() const
{
	return zn_hash_.end();
}

ZnCoordinationScorer::ZnCoordinationScorer() : utility::pointer::ReferenceCount(),
	well_depth_( 1.0 ),
	idealize_input_virtual_atoms_( true ),
	max_natoms_(0),
	clash_weight_( 1.0 ),
	require_3H_( false ),
	znreach_( 3.0 ),
	orbital_dist_( 1.0 ),
	orbital_reach_( 1.0 ),
	znwelldepth_( 3.0 ),
	asymm_atids_( core::id::AtomID( 0, 0 ) ),
	focused_clone_atids_( core::id::AtomID( 0, 0 ) ),
	third_resid_( 0u, 0u )
{
	reset_reach(); // the reach for the hash is deteremined by the znreach, the orbital distances from the zn and the orbital reach.

	// by default, the asymm_atids_ are initialized to atoms 1,2, and 3 on residue 1.
	// the focused_clone_atids are not set and must be set explicitly.
	asymm_atids_[1] = core::id::AtomID( 1, 1 );
	asymm_atids_[2] = core::id::AtomID( 2, 1 );
	asymm_atids_[3] = core::id::AtomID( 3, 1 );

	if ( ! core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->has_name( "ZNX" ) ) {
		utility_exit_with_message( "Error in construction of ZnCoordinateConstraint: ResidueType ZNX has not been defined" );
	}

	// Measure the virtual atom coordinates in an idealized frame
	core::chemical::ResidueType const & znx_restype =
		core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->
		name_map( "ZNX" );

	core::Vector znpos = znx_restype.atom( znx_restype.atom_index( "ZN" )).ideal_xyz();
	core::Vector v1pos = znx_restype.atom( znx_restype.atom_index( "V1" )).ideal_xyz();
	core::Vector v2pos = znx_restype.atom( znx_restype.atom_index( "V2" )).ideal_xyz();
	core::Vector v12_midpoint = v1pos + v2pos / 2;

	HTReal znframe( v1pos, v12_midpoint, znpos );

	/// Now measure the coordinates of the four virtual atoms in this idealized frame
	/// IE convert them to their local coordinates.
	HTReal invznframe = znframe.inverse();
	znx_ideal_coords_.resize( znx_restype.natoms() );
	for ( core::Size ii = 1; ii <= znx_restype.natoms(); ++ii ) {
		znx_ideal_coords_[ ii ] = invznframe * znx_restype.atom( ii ).ideal_xyz();
	}

	reset_znx_orbital_coords();

	// Simple polynomial: derivative of zero at 0 and 1; values of
	// -1 and 0 at 0 and 1.
	utility::vector1< core::Real > coeffs( 4 );
	coeffs[ 1 ] = -2; coeffs[ 2 ] =  3; coeffs[ 3 ] =  0; coeffs[ 4 ] = -1;

	ramp_to_zero_poly_ = numeric::Polynomial_1dOP( new numeric::Polynomial_1d(
		"ramp", 0, 1, -1, 0, 1, 12, 4, coeffs ) );
}

ZnCoordinationScorer::~ZnCoordinationScorer() {}

ZnCoordinationScorer::ZnCoordinationScorer( ZnCoordinationScorer const & src ) : utility::pointer::ReferenceCount(),
	znx_ideal_coords_( src.znx_ideal_coords_ ),
	ramp_to_zero_poly_( src.ramp_to_zero_poly_ ),
	well_depth_( src.well_depth_ ),
	reach_( src.reach_ ),
	reach2_( src.reach2_ ),
	reference_frame_( src.reference_frame_ ),
	inv_reference_frame_( src.inv_reference_frame_ ),
	idealize_input_virtual_atoms_( src.idealize_input_virtual_atoms_ ),
	zn_matches_( src.zn_matches_ ),
	max_natoms_( src.max_natoms_ ),
	clash_weight_( src.clash_weight_ ),
	require_3H_( src.require_3H_ ),
	znreach_( src.znreach_ ),
	orbital_dist_( src.orbital_dist_ ),
	orbital_reach_( src.orbital_reach_ ),
	znwelldepth_( src.znwelldepth_ ),
	asymm_chain_( src.asymm_chain_ ),
	focsed_clone_chain_( src.focsed_clone_chain_ ),
	asymm_atids_( src.asymm_atids_ ),
	focused_clone_atids_( src.focused_clone_atids_ ),
	third_resid_( src.third_resid_ ),
	hash_( src.hash_ )
{}

/// @brief set an individual atom ID to use for determining the refrence frame for the asymmetric unit
/// There are three atoms total which are used to define the reference frame.  By default, these are
/// initialized to atoms 1,2,&3 on reisude 1.  All asymm ids
void ZnCoordinationScorer::set_asymm_atid( Size which_atid, core::id::AtomID atid )
{
	asymm_atids_[ which_atid ] = atid;
}

/// @brief set an individual atom ID to use for determining the reference frame for the symmetric clone.
void ZnCoordinationScorer::set_symm_atid( Size which_atid, core::id::AtomID atid )
{
	focused_clone_atids_[ which_atid ] = atid;
}

/// @brief set the residue on the asymmetric unit and use atoms 1,2,&3 on that residue to define the reference
/// frame.
void ZnCoordinationScorer::set_asymm_resid( Size resid )
{
	for ( Size ii = 1; ii <= 3; ++ii ) asymm_atids_[ii] = core::id::AtomID( ii, resid );
}

/// @brief set the residue on the symmetric clone and use atoms 1,2,&3 on that residue to define the reference
/// frame.
void ZnCoordinationScorer::set_symm_resid( Size resid )
{
	for ( Size ii = 1; ii <= 3; ++ii ) focused_clone_atids_[ii] = core::id::AtomID( ii, resid );
}

void ZnCoordinationScorer::set_third_resid( Size resid ) {
	core::id::AtomID new_third_resid( 1, resid );
	third_resid_ = new_third_resid;
}

void ZnCoordinationScorer::set_zn_reach( Real reach ) {
	znreach_ = reach;
	reset_reach();
}

void ZnCoordinationScorer::set_orbital_dist( Real dist ) {
	orbital_dist_ = dist;
	reset_znx_orbital_coords();
	reset_reach();
}

void ZnCoordinationScorer::set_orbital_reach( Real reach )
{
	orbital_reach_ = reach;
	reset_reach();
}

void ZnCoordinationScorer::set_zn_well_depth( Real depth )
{
	znwelldepth_ = depth;
}

void ZnCoordinationScorer::set_reference_pdb(
	std::string const & start_pdb
)
{
	core::pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, start_pdb );
	reference_frame_ = HTReal( pose.xyz( asymm_atids_[1] ),pose.xyz( asymm_atids_[2] ),pose.xyz( asymm_atids_[3] ) );
	inv_reference_frame_ = reference_frame_.inverse();
}

void ZnCoordinationScorer::set_matcher_constraint_file_name( std::string const & fname )
{
	matchcst_file_name_ = fname;
}


/// @details Since this step requires converting the match coordinates to the local coordinates
/// in the reference frame, the reference PDB must already have been set.
void ZnCoordinationScorer::add_match_from_file(
	std::string const & match_file_name
)
{
	assert( matchcst_file_name_ != "" );

	// Load the match pdb.
	// Place the ZN virtual atoms.
	// Transform their coordinates to the local coordinates in the reference frame
	// Save details about this match: file name.

	assert( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->has_name( "ZNX" ));


	/*core::chemical::ResidueType const & znx_restype =
		core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->
		name_map( "ZNX" ); */ // Unused variable causes warning.

	core::pose::Pose match_pose;
	core::import_pose::pose_from_pdb( match_pose, match_file_name );

	/// Input validation.
	if ( match_pose.total_residue() != 3 || match_pose.residue(3).name() != "ZNX" ) {
		std::cerr << "Error, match file " << match_file_name << " does not conform to the standard" << std::endl;
		if ( match_pose.total_residue() != 3 ) {
			std::cerr << "Total number of residues is not equal to 3 (read " << match_pose.total_residue()  << " residues)" << std::endl;
		} else {
			std::cerr << "Residue 3 is not a ZNX residue!  Instead, its residue type " << match_pose.residue(3).name() << std::endl;
		}
		utility_exit_with_message( "Bad input to ZnCoordinateConstraint add_match_from_file()" );
	}

	//
	Size nhistadines( 0 );
	for ( Size ii = 1; ii <= 2; ++ii ) {
		if ( match_pose.residue(ii).aa() == core::chemical::aa_his ) ++nhistadines;
	}

	// The next block of code either idealizes the input match geometry,
	// or uses the input match geometry. Both paths store their computations
	// in the next two variables.
	ZnCoord zn_match_coords;
	core::conformation::ResidueOP znres;

	if ( idealize_input_virtual_atoms_ ) {
		// determine which are the coordinating atoms by asking the
		// Enzyme design constraint addition code to identify those atoms.
		protocols::enzdes::AddOrRemoveMatchCsts cst_adder;
		cst_adder.cstfile( matchcst_file_name_ );
		cst_adder.set_cst_action( protocols::enzdes::ADD_NEW );

		cst_adder.apply( match_pose );

		protocols::enzdes::EnzdesConstraintReporter reporter;
		reporter.ligand_resno( 3 );
		reporter.find_constraints_to_ligand( match_pose );

		Size const num_geom_csts = 2;
		utility::vector1< std::pair< core::id::AtomID, core::Real > > closest_atoms( num_geom_csts, std::make_pair( core::id::AtomID(), -1 ));

		core::id::AtomID zn_atom_id( match_pose.residue(3).atom_index( "ZN" ), 3);

		for ( core::Size ii = 1; ii <= reporter.constrained_nonligand_atoms().size(); ++ii ) {
			core::id::AtomID iiatid = reporter.constrained_nonligand_atoms()[ ii ];
			assert( iiatid.rsd() == 1 || iiatid.rsd() == 2 );
			//std::cout << "Constraint to ligand #" << ii << " Res " <<
			//	iiatid.rsd() << " atom " <<
			//	match_pose.residue( iiatid.rsd() ).atom_name( iiatid.atomno() ) <<
			//	std::endl;
			core::Real const d2 = match_pose.xyz(zn_atom_id).distance_squared( match_pose.xyz( iiatid ) );
			if ( closest_atoms[ iiatid.rsd() ].second < 0 || d2 < closest_atoms[ iiatid.rsd() ].second ) {
				closest_atoms[ iiatid.rsd() ].first = iiatid;
				closest_atoms[ iiatid.rsd() ].second = d2;
			}
		}

		//for ( core::Size ii = 1; ii <= 2; ++ii ) {
		//	std::cout << "Closest atom rsd: " << closest_atoms[ii].first.rsd() << " " <<
		//		match_pose.residue( closest_atoms[ii].first.rsd() ).atom_name( closest_atoms[ii].first.atomno() ) <<
		//		" with distance " << std::sqrt( closest_atoms[ii].second ) << std::endl;
		//}

		core::Vector r1coord_atom = match_pose.xyz( closest_atoms[1].first );
		core::Vector r2coord_atom = match_pose.xyz( closest_atoms[2].first );
		//core::Vector r12_halfway = r1coord_atom + r2coord_atom / 2;

		core::Vector zn_atom = match_pose.xyz( core::id::AtomID( 1, 3 ) );
		core::Vector zn_vect = ((r1coord_atom-zn_atom).normalize() + (r2coord_atom-zn_atom).normalize()) / 2;

		core::Vector walk_from_zn = zn_atom + zn_vect;

		HTReal actual_frame( r1coord_atom, walk_from_zn, zn_atom );
		HTReal transform_frame = inv_reference_frame_ * actual_frame;
		ZnCoord idealized_zn_match_coords(
			zn_matches_.size() + 1,
			nhistadines,
			transform_frame * znx_ideal_coords_[ 1 ],
			transform_frame * znx_ideal_coords_[ 2 ],
			transform_frame * znx_ideal_coords_[ 3 ],
			transform_frame * znx_ideal_coords_[ 4 ],
			transform_frame * znx_ideal_coords_[ 5 ] );

		zn_match_coords = idealized_zn_match_coords;

		// also save the coordinates of the idealized virtual atoms
		znres = match_pose.residue(3).clone();
		for ( Size ii = 2; ii <= 5; ++ii ) {
			znres->set_xyz( ii, actual_frame * znx_ideal_coords_[ ii ] );
		}
	} else {
		znres = match_pose.residue(3).clone();

		// take the coordinates of the input virtual atoms and pull them
		// toward the zinc so that the magnitude of the distance between the
		// virtual atoms and the zinc matches the desired length.
		core::conformation::ResidueOP tmp_coords = match_pose.residue(3).clone();
		for ( Size ii = 2; ii <= tmp_coords->natoms(); ++ii ) {
			core::Vector zn2vrt = tmp_coords->xyz(ii) - tmp_coords->xyz(1);
			core::Vector zn2vrt_scaled = zn2vrt.normalize() *  orbital_dist_;
			tmp_coords->set_xyz( ii, tmp_coords->xyz(1) + zn2vrt_scaled );
		}
		ZnCoord raw_zn_match_coords(
			zn_matches_.size() + 1,
			nhistadines,
			inv_reference_frame_ * tmp_coords->xyz( 1 ),
			inv_reference_frame_ * tmp_coords->xyz( 2 ),
			inv_reference_frame_ * tmp_coords->xyz( 3 ),
			inv_reference_frame_ * tmp_coords->xyz( 4 ),
			inv_reference_frame_ * tmp_coords->xyz( 5 ));
		zn_match_coords = raw_zn_match_coords;
	}

	ZnMatchData zndat;
	zndat.index( zn_matches_.size() + 1 );
	zndat.res1( match_pose.pdb_info()->number( 1 ) );
	zndat.res2( match_pose.pdb_info()->number( 2 ) );
	zndat.zn_and_orbitals( zn_match_coords );
	zndat.match_pdb_file( match_file_name );
	zndat.match_cst_file( matchcst_file_name_ );
	zndat.res1conf( match_pose.residue(1) );
	zndat.res2conf( match_pose.residue(2) );
	zndat.znconf( *znres );

	if ( max_natoms_ < match_pose.residue(1).natoms() ) max_natoms_ = match_pose.residue(1).natoms();
	if ( max_natoms_ < match_pose.residue(2).natoms() ) max_natoms_ = match_pose.residue(2).natoms();

	zn_matches_.push_back( zndat );
}

void ZnCoordinationScorer::add_matches_from_files(
	std::list< std::string > const &
)
{}

void ZnCoordinationScorer::add_match_from_istream( std::istream & )
{}

void ZnCoordinationScorer::finalize_after_all_matches_added()
{
	hash_ = ZnHashOP( new ZnHash );
	hash_->set_uniform_bin_width( reach_ );
	for ( Size ii = 1; ii <= zn_matches_.size(); ++ii ) {
		hash_->add_zn_coordinate( zn_matches_[ ii ].zn_and_orbitals() );
	}
	hash_->build_hash();
}

bool
ZnCoordinationScorer::optimal_coordination_is_reversed(
	core::pose::Pose const & p
)
{
	CoordinationData result = score_and_index_for_best_match( p.residue(r1()), p.residue(r2()) );
	return result.second.second;
}

core::Real
ZnCoordinationScorer::score( core::pose::Pose const & p ) const
{
	return score( p.residue(r1()), p.residue(r2()) );
}

core::Real
ZnCoordinationScorer::score(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2
) const
{
	CoordinationData result = score_and_index_for_best_match( res1, res2 );
	return result.second.first;

}


/// @details zn1 coordinates chain A with virtual atoms 1 and 2; zn2 coordinates chain B also with virtual atoms 1 and 2.
/// So the correspondence we want to find is virtual atoms 1 and 2 on zn1 with virtual atoms 3 and 4 on zn2, and
/// virtual atoms 3 and 4 on zn1 with virtual atoms 1 and 2 on zn2.  These can be in any order, so examine all.
/// This function returns a pair: the score, and a boolean representing whether zn2 is being coordinated "in reverse"
///
ZnCoordinationScorer::ZnScoreAndFlipState
ZnCoordinationScorer::score_zn_pair( ZnCoord const & zn1, ZnCoord const & zn2 ) const
{
	Real zn_d2 = zn1[1].distance_squared( zn2[1] );
	if ( zn_d2 < reach2_ ) {
		Real const znscore = znwelldepth_*ramp_to_zero_poly_->eval( std::sqrt( zn_d2 ) / znreach_ );

		//Real A12score(0), B12score(0), A34score(0), B34score(0);
		//Real d2_v1v1A = zn1[2].distance_squared( zn2[4] );
		//Real d2_v1v1B = zn1[2].distance_squared( zn2[5] );

		//Real d2_v2v2A = zn1[3].distance_squared( zn2[5] );
		//Real d2_v2v2B = zn1[3].distance_squared( zn2[4] );


		//Real d2_v3v3A = zn1[4].distance_squared( zn2[2] );
		//Real d2_v3v3B = zn1[4].distance_squared( zn2[3] );

		//Real d2_v4v4A = zn1[5].distance_squared( zn2[5] );
		//Real d2_v4v4B = zn1[5].distance_squared( zn2[4] );

		//if ( d2_v1v1A < 1 ) { A12score += ramp_to_zero_poly_->eval( std::sqrt( d2_v1v1A )); }
		//if ( d2_v2v2A < 1 ) { A12score += ramp_to_zero_poly_->eval( std::sqrt( d2_v2v2A )); }
		//if ( d2_v1v1B < 1 ) { B12score += ramp_to_zero_poly_->eval( std::sqrt( d2_v1v1B )); }
		//if ( d2_v2v2B < 1 ) { B12score += ramp_to_zero_poly_->eval( std::sqrt( d2_v2v2B )); }

		//if ( d2_v3v3A < 1 ) { A34score += ramp_to_zero_poly_->eval( std::sqrt( d2_v3v3A )); }
		//if ( d2_v4v4A < 1 ) { A34score += ramp_to_zero_poly_->eval( std::sqrt( d2_v4v4A )); }
		//if ( d2_v3v3B < 1 ) { B34score += ramp_to_zero_poly_->eval( std::sqrt( d2_v3v3B )); }
		//if ( d2_v4v4B < 1 ) { B34score += ramp_to_zero_poly_->eval( std::sqrt( d2_v4v4B )); }

		//return znscore + numeric::min( A12score, B12score ) + numeric::min( A34score, B34score );

		// Take 2: only two forms of chirality
		// i)  v1a -- v1b, v2a -- v2b, v3a -- v3b, v4a -- v4b, and
		// ii) v1a -- v1b, v2a -- v2b, v3a -- v4b, v4a -- v3b, and

		/// Which, because Zn2 is coordinated at v1 and v2 gives us this mapping:
		// i)  v1zn1 -- v3zn1, v2zn1 -- v4zn2, v3zn1 -- v1zn2, v4zn1 -- v2zn2, and
		// ii) v1zn1 -- v3zn1, v2zn1 -- v4zn2, v3zn1 -- v2zn2, v4zn1 -- v1zn2, and

		Real i_score(0), ii_score(0);
		Real d2_v1 = zn1.orbital(1).distance_squared( zn2.orbital(3) );
		Real d2_v2 = zn1.orbital(2).distance_squared( zn2.orbital(4) );

		Real d2_v3_i  = zn1.orbital(3).distance_squared( zn2.orbital(1) );
		Real d2_v4_i  = zn1.orbital(4).distance_squared( zn2.orbital(2) );
		Real d2_v3_ii = zn1.orbital(3).distance_squared( zn2.orbital(2) );
		Real d2_v4_ii = zn1.orbital(4).distance_squared( zn2.orbital(1) );

		if ( d2_v1 < 1 ) { ii_score = i_score += ramp_to_zero_poly_->eval( std::sqrt( d2_v1 ) / orbital_reach_ ); }
		if ( d2_v2 < 1 ) { ii_score = i_score += ramp_to_zero_poly_->eval( std::sqrt( d2_v2 ) / orbital_reach_ ); }

		if ( d2_v3_i  < 1 ) { i_score  += ramp_to_zero_poly_->eval( std::sqrt( d2_v3_i  ) / orbital_reach_ ); }
		if ( d2_v4_i  < 1 ) { i_score  += ramp_to_zero_poly_->eval( std::sqrt( d2_v4_i  ) / orbital_reach_ ); }
		if ( d2_v3_ii < 1 ) { ii_score += ramp_to_zero_poly_->eval( std::sqrt( d2_v3_ii ) / orbital_reach_ ); }
		if ( d2_v4_ii < 1 ) { ii_score += ramp_to_zero_poly_->eval( std::sqrt( d2_v4_ii ) / orbital_reach_ ); }

		return std::make_pair( znscore + numeric::min( i_score, ii_score ), i_score < ii_score );
	}
	return std::make_pair( 0.0, true );
}

ZnCoordinationScorer::ZnIndexPair
ZnCoordinationScorer::best_match( core::pose::Pose const & p ) const
{
	return best_match( p.residue(r1()), p.residue(r2()) );
}

ZnCoordinationScorer::ZnIndexPair
ZnCoordinationScorer::best_match(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2
) const
{
	CoordinationData result = score_and_index_for_best_match( res1, res2 );
	return result.first;
}

ZnCoordinationScorer::CoordinationData
ZnCoordinationScorer::score_and_index_for_best_match(
	core::conformation::Residue const & res1,
	core::conformation::Residue const & res2
) const
{
	CoordinationData result;
	result.first = std::make_pair( core::Size(0), core::Size(0) );
	result.second.first = 0;
	result.second.second = false;

	assert( hash_ );

	HTReal chAframe( res1.xyz( asymm_atids_[1].atomno() ), res1.xyz( asymm_atids_[2].atomno() ), res1.xyz( asymm_atids_[3].atomno() ) );
	HTReal chBframe( res2.xyz( focused_clone_atids_[1].atomno() ), res2.xyz( focused_clone_atids_[2].atomno() ), res2.xyz( focused_clone_atids_[3].atomno() ) );

	//HTReal to_hash_frame_transform = query_frame_to_original_frame( res1, res2 );
	HTReal to_hash_frame_transform = chAframe.inverse() * chBframe;
	HTReal clash_frame_transform   = reference_frame_ * chAframe.inverse() * chBframe *  inv_reference_frame_;

	utility::fixedsizearray1< Size, 3 > threes(3);
	utility::FixedSizeLexicographicalIterator< 3 > lex( threes );
	utility::fixedsizearray1< Real, 3 > step_sizes(0);
	step_sizes[1] = -reach_;
	step_sizes[3] =  reach_;

	numeric::geometry::BoundingBox< numeric::xyzVector< core::Real > > const & bb( hash_->bb());
	numeric::geometry::BoundingBox< numeric::xyzVector< core::Real > > const & bb_ext( hash_->bb_ext());

	/// for calculating the transformed coordinates of the
	utility::vector1< core::Vector > symmclone_coords_res1( max_natoms_, core::Vector( 0.0 ) );
	utility::vector1< core::Vector > symmclone_coords_res2( max_natoms_, core::Vector( 0.0 ) );

	for ( ZnHash::ZnCoordinateHash::const_iterator
			hash_iter = hash_->hash_begin(), hash_iter_end = hash_->hash_end();
			hash_iter != hash_iter_end; ++hash_iter ) {
		utility::vector1< ZnCoord > const & zncoords( hash_iter->second );
		for ( Size ii = 1, iiend = zncoords.size(); ii <= iiend; ++ii ) {
			ZnCoord const & iicoord = zncoords[ ii ];
			core::Vector ii_zn_local_coord = to_hash_frame_transform * iicoord[1];
			if ( ! bb_ext.contains( ii_zn_local_coord ) ) continue;
			ZnCoord ii_local_coords(
				iicoord.index(),
				iicoord.nhis(),
				ii_zn_local_coord,
				to_hash_frame_transform * iicoord[2],
				to_hash_frame_transform * iicoord[3],
				to_hash_frame_transform * iicoord[4],
				to_hash_frame_transform * iicoord[5] );
			ZnMatchData const & ii_match_data = zn_matches_[ iicoord.index() ];
			for ( lex.begin(); ! lex.at_end(); ++lex ) {
				core::Vector query_coord;
				for ( Size jj = 1; jj <= 3; ++jj ) query_coord(jj) = ii_zn_local_coord(jj)+step_sizes[lex[jj]];
				if ( ! bb.contains( query_coord ) ) continue;
				ZnHash::ZnCoordinateHash::const_iterator query_hits = hash_->query_hash( query_coord );
				if ( query_hits == hash_iter_end ) continue;
				utility::vector1< ZnCoord > const & query_hit_coords( query_hits->second );
				for ( Size jj = 1, jjend = query_hit_coords.size(); jj <= jjend; ++jj ) {
					ZnMatchData const & jj_match_data = zn_matches_[ query_hit_coords[ jj ].index() ];
					if ( ii_match_data.res1() == jj_match_data.res1() ||
							ii_match_data.res1() == jj_match_data.res2() ||
							ii_match_data.res2() == jj_match_data.res1() ||
							ii_match_data.res2() == jj_match_data.res2() ) {
						continue; // don't consider a match where one residue is being used twice
					}

					if ( require_3H_ && iicoord.nhis() + query_hit_coords[ jj ].nhis() != 3 ) {
						//std::cout << "ii nhis: " << iicoord.nhis() << " query nhis: " << query_hit_coords[ jj ].nhis() << std::endl;
						continue;
					}

					ZnScoreAndFlipState score_and_flip = score_zn_pair( query_hit_coords[ jj ], ii_local_coords );
					if ( score_and_flip.first < 0 ) {
						Real clash = clash_score( query_hit_coords[ jj ].index(), iicoord.index(),
							clash_frame_transform,
							symmclone_coords_res1, symmclone_coords_res2 );
						score_and_flip.first += clash_weight_ * clash;
					}
					if ( score_and_flip.first < result.second.first ) {
						result.second = score_and_flip;
						result.first = std::make_pair( query_hit_coords[ jj ].index(), iicoord.index() );
						//std::cout << "new best pair: " << result.first.first << " and " << result.first.second << " score: " << result.second.first << std::endl;
					}
				}
			}
		}
	}
	//std::cout << "Zn Hash Score: " << result.second << std::endl;
	return result;
}


core::Real
ZnCoordinationScorer::clash_score(
	Size m1,
	Size m2,
	HTReal const & to_hash_frame_transform,
	utility::vector1< core::Vector > & symmclone_coords_res1,
	utility::vector1< core::Vector > & symmclone_coords_res2
) const
{
	ZnMatchData const & match1 = zn_matches_[ m1 ];
	ZnMatchData const & match2 = zn_matches_[ m2 ];
	core::conformation::Residue const & r1m1 = match1.res1conf();
	core::conformation::Residue const & r2m1 = match1.res2conf();

	core::conformation::Residue const & r1m2 = match2.res1conf();
	core::conformation::Residue const & r2m2 = match2.res2conf();

	assert( symmclone_coords_res1.size() >= r1m2.natoms() );
	assert( symmclone_coords_res2.size() >= r2m2.natoms() );


	/*for ( Size ii = 1; ii <= r1m1.natoms(); ++ii ) {
		std::cout << " r1m1 atom " << ii
			<< " " << r1m1.xyz(ii).x()
			<< " " << r1m1.xyz(ii).y()
			<< " " << r1m1.xyz(ii).z() << std::endl;
	}
	for ( Size ii = 1; ii <= r2m1.natoms(); ++ii ) {
		std::cout << " r2m1 atom " << ii
			<< " " << r1m1.xyz(ii).x()
			<< " " << r1m1.xyz(ii).y()
			<< " " << r1m1.xyz(ii).z() << std::endl;
	}*/

	for ( Size ii = 1; ii <= r1m2.natoms(); ++ii ) {
		symmclone_coords_res1[ ii ] = to_hash_frame_transform*r1m2.xyz(ii);
		//std::cout << " r1m2 atom " << ii
		//	<< " " << symmclone_coords_res1[ ii ].x()
		//	<< " " << symmclone_coords_res1[ ii ].y()
		//	<< " " << symmclone_coords_res1[ ii ].z() << std::endl;

	}
	for ( Size ii = 1; ii <= r2m2.natoms(); ++ii ) {
		symmclone_coords_res2[ ii ] = to_hash_frame_transform*r2m2.xyz(ii);
		//std::cout << " r2m2 atom " << ii
		//	<< " " << symmclone_coords_res2[ ii ].x()
		//	<< " " << symmclone_coords_res2[ ii ].y()
		//	<< " " << symmclone_coords_res2[ ii ].z() << std::endl;

	}

	return clash_score_residue_pair( r1m1, r1m2,  symmclone_coords_res1 ) +
		clash_score_residue_pair( r1m1, r2m2,  symmclone_coords_res2 ) +
		clash_score_residue_pair( r2m1, r1m2,  symmclone_coords_res1 ) +
		clash_score_residue_pair( r2m1, r2m2,  symmclone_coords_res2 );

}

core::Real
ZnCoordinationScorer::clash_score_residue_pair(
	core::conformation::Residue const & r1,
	core::conformation::Residue const & r2,
	utility::vector1< core::Vector > const & r2coords
) const
{
	Real score = 0;
	Real const hvd2 = 2.8*2.8;
	Real const hvhd2 = 2.4*2.4;
	Real const hd2 = 2.0*2.0;
	for ( Size ii = 1, iiend = r1.nheavyatoms(), jjend = r2.nheavyatoms(); ii <= iiend; ++ii ) {
		for ( Size jj = 1; jj <= jjend; ++jj ) {
			//Real d2 = r1.xyz(ii).distance_squared( r2coords[ jj ] );
			//if ( d2 < 100 ) {
			//	std::cout << " d2 : " << ii << " " << jj << " " << d2 << std::endl;
			//}
			Real d2_minus_hvd2 = r1.xyz(ii).distance_squared( r2coords[ jj ] ) - hvd2;
			if ( d2_minus_hvd2 < 0 ) {
				score += std::sqrt( -1 * d2_minus_hvd2); // linear distance penalty
			}
		}
	}
	// hydrogen / heavyatom collision
	for ( Size ii = 1, iiend = r1.nheavyatoms(), jjstart = r2.nheavyatoms()+1, jjend = r2.natoms(); ii <= iiend; ++ii ) {
		for ( Size jj = jjstart; jj <= jjend; ++jj ) {
			Real d2_minus_hvhd2 = r1.xyz(ii).distance_squared( r2coords[ jj ] ) - hvhd2;
			if ( d2_minus_hvhd2 < 0 ) {
				score += std::sqrt( -1 * d2_minus_hvhd2); // linear distance penalty
			}
		}
	}

	// hydrogen / heavyatom collision
	for ( Size ii = r1.nheavyatoms()+1, iiend = r1.natoms(), jjstart = 1, jjend = r2.nheavyatoms(); ii <= iiend; ++ii ) {
		for ( Size jj = jjstart; jj <= jjend; ++jj ) {
			Real d2_minus_hvhd2 = r1.xyz(ii).distance_squared( r2coords[ jj ] ) - hvhd2;
			if ( d2_minus_hvhd2 < 0 ) {
				score += std::sqrt( -1 * d2_minus_hvhd2); // linear distance penalty
			}
		}
	}

	for ( Size ii = r1.nheavyatoms()+1, iiend = r1.natoms(), jjstart = r2.nheavyatoms()+1, jjend = r2.natoms(); ii <= iiend; ++ii ) {
		for ( Size jj = jjstart; jj <= jjend; ++jj ) {
			Real d2_minus_hd2 = r1.xyz(ii).distance_squared( r2coords[ jj ] ) - hd2;
			if ( d2_minus_hd2 < 0 ) {
				score += std::sqrt( -1 * d2_minus_hd2); // linear distance penalty
			}
		}
	}
	return score;
}

void
ZnCoordinationScorer::insert_match_onto_pose(
	core::pose::Pose & p, // <-- should be a symmetric pose if this match is to be applied to all clones
	Size match_index,
	Size chain_insertion_id // 1 for the asymm unit, 2 for the symm clone.
) const
{
	ZnMatchData const & match = zn_matches_[ match_index ];

	utility::fixedsizearray1< Size, 2 > destresids;
	destresids[1] = chain_insertion_id == 1 ? match.res1() : match.res1() + r2() - r1();
	destresids[2] = chain_insertion_id == 1 ? match.res2() : match.res2() + r2() - r1();

	core::chemical::ResidueTypeSet const & restypeset( p.residue_type(r1()).residue_type_set() );

	for ( Size ii = 1; ii <= 2; ++ii ) {
		core::conformation::Residue const & matchres = ii == 1 ? match.res1conf() : match.res2conf();
		core::conformation::Residue const & dstres = p.residue( destresids[ ii ] );
		assert( &restypeset == & dstres.type().residue_type_set() );
		assert( &restypeset == & matchres.type().residue_type_set() );

		// find the appropriate residue type for this position
		core::chemical::ResidueTypeCOP newrestype( matchres.type().get_self_ptr() );
		utility::vector1< std::string > const & matchres_variants =
				matchres.type().properties().get_list_of_variants();
		for ( Size jj = 1; jj <= matchres_variants.size(); ++jj ) {
			if ( ! dstres.type().has_variant_type( matchres_variants[ jj ]  )) {
				core::chemical::ResidueTypeCOP variantfree_newrestype;
				// TODO: Refactor this to avoid working with strings.
				variantfree_newrestype =
						restypeset.get_residue_type_with_variant_removed(
								*newrestype,
								core::chemical::ResidueProperties::get_variant_from_string( matchres_variants[ jj ] )
						).get_self_ptr();
				if ( ! variantfree_newrestype  ) {
					std::cerr << "Error could not remove variant " << matchres_variants[ jj ] << " from restype " <<
						newrestype->name() << std::endl;
					utility_exit_with_message( "get_residue_type_with_variant_removed failed" );
				} else {
					newrestype = variantfree_newrestype;
				}
			}
		}
		utility::vector1< std::string > const & dstres_variants =
				dstres.type().properties().get_list_of_variants();
		for ( Size jj = 1; jj <= dstres_variants.size(); ++jj ) {
			if ( ! newrestype->has_variant_type( dstres_variants[ jj ]  )) {
				core::chemical::ResidueTypeCOP variantful_newrestype;
				// TODO: Refactor this to avoid working with strings.
				variantful_newrestype =
						restypeset.get_residue_type_with_variant_added(
								*newrestype,
								core::chemical::ResidueProperties::get_variant_from_string( dstres_variants[ jj ] )
						).get_self_ptr();
				if ( ! variantful_newrestype  ) {
					std::cerr << "Error could not add variant " << dstres_variants[ jj ] << " to restype " <<
						newrestype->name() << std::endl;
					utility_exit_with_message( "get_residue_type_with_variant_added failed" );
				} else {
					newrestype = variantful_newrestype;
				}
			}
		}

		// This should be unreachable.
		if ( ! variants_match( *newrestype, dstres.type() ) ) {
			std::cerr << "ERROR: could not get newrestype and dstres to agree on variant types.\n";
			std::cerr << "Dstres (" << dstres.name() << ") variants:\n";
			for ( Size jj = 1; jj <= dstres_variants.size(); ++jj ) {
				std::cerr << "  " << dstres_variants[ jj ] << "\n";
			}
			std::cerr << "Newrestype variants:\n" << std::endl;
			for ( Size jj = 1; jj <= newrestype->properties().get_list_of_variants().size(); ++jj ) {
				std::cerr << "  " << newrestype->properties().get_list_of_variants()[ jj ] << "\n";
			}

			utility_exit_with_message( "Could not match the variant types for newres and dstres" );
		}

		TR << "Inserting zinc coordinating residue " << newrestype->name() << " at " << dstres.seqpos() << std::endl;

		core::conformation::ResidueOP newres = core::conformation::ResidueFactory::create_residue(
			*newrestype, dstres, p.conformation(), false );

		// now transform the match coordinates onto the new backbone position.
		HTReal dstframe( dstres.xyz(1), dstres.xyz(2), dstres.xyz(3));
		HTReal matchframe( matchres.xyz(1), matchres.xyz(2), matchres.xyz(3));
		HTReal transform = dstframe * matchframe.inverse();

		utility::vector1< Size > newcoords_calculated( newres->natoms(), 0 );
		for ( Size jj = 1; jj <= matchres.natoms(); ++jj ) {
			if ( newres->has( matchres.atom_name( jj ) )) {
				Size newres_atind = newres->atom_index( matchres.atom_name( jj ) );
				newcoords_calculated[ newres_atind ] = 1;
				newres->set_xyz( newres_atind, transform * matchres.xyz(jj) );
			}
		}
		// make sure we got everything -- probably shouldn't be an assert, but an actual if
		for ( Size jj = 1; jj <= newcoords_calculated.size(); ++jj ) {
			assert( newres->atom_is_backbone(jj) || newcoords_calculated[ jj ] == 1 );
		}

		// now insert this residue into the pose
		p.replace_residue( destresids[ ii ], *newres, false );
	}

}



// @brief Compute the homogeneous transform that will take Zn coordinates from the reference frame,
// translate them onto the the location of the coordinates on the slave clone, and then
// translate them further into their effective location in the original reference frame given
// the new location of the master clone.  Pretty simple, really.  Public for testing purposes only.
ZnCoordinationScorer::HTReal
ZnCoordinationScorer::query_frame_to_original_frame(
	core::conformation::Residue const & r1,
	core::conformation::Residue const & r2
) const
{
	HTReal chAframe( r1.xyz( asymm_atids_[1].atomno() ), r1.xyz( asymm_atids_[2].atomno() ), r1.xyz( asymm_atids_[3].atomno() ) );
	HTReal chBframe( r2.xyz( focused_clone_atids_[1].atomno() ), r2.xyz( focused_clone_atids_[2].atomno() ), r2.xyz( focused_clone_atids_[3].atomno() ) );

	return chAframe.inverse() * chBframe;
}

ZnCoordinationScorer::HTReal
ZnCoordinationScorer::query_frame_to_original_frame( core::pose::Pose const & p ) const
{
	//HTReal chAframe( p.xyz( asymm_atids_[1] ), p.xyz( asymm_atids_[2] ), p.xyz( asymm_atids_[3] ) );
	//HTReal chBframe( p.xyz( focused_clone_atids_[1] ), p.xyz( focused_clone_atids_[2] ), p.xyz( focused_clone_atids_[3] ) );

	//return chAframe.inverse() * chBframe;
	return query_frame_to_original_frame( p.residue( asymm_atids_[1].rsd() ), p.residue( focused_clone_atids_[1].rsd() ));
}

ZnCoord
ZnCoordinationScorer::original_frame_coordinate_for_match(
	Size index,
	core::pose::Pose const & p
) const
{

	HTReal transform = query_frame_to_original_frame( p );
	ZnCoord const & zn( zn_matches_[ index ].zn_and_orbitals() );
	return ZnCoord(
		index,
		zn.nhis(),
		transform * zn[ 1 ],
		transform * zn[ 2 ],
		transform * zn[ 3 ],
		transform * zn[ 4 ],
		transform * zn[ 5 ] );
}

void
ZnCoordinationScorer::reset_reach()
{
	assert( ! hash_ ); // must be called before finalize gets called if you're going to call it.
	reach_ = std::max( znreach_, orbital_dist_ * 2 + orbital_reach_ );
	reach2_ = reach_ * reach_;
}


void
ZnCoordinationScorer::reset_znx_orbital_coords()
{
	assert( zn_matches_.size() == 0 ); // no matches should have been added yet.
	// Scale the coordinates so that the virtual atoms are unit length from the center
	for ( core::Size ii = 1; ii <= znx_ideal_coords_.size(); ++ii ) {
		Real d = znx_ideal_coords_[ii].length();
		if ( d > 1e-6 ) {
			// don't move the zinc; it should be right at the origin.
			znx_ideal_coords_[ ii ] /= d;
			znx_ideal_coords_[ ii ] *= orbital_dist_;
		}
	}
}


ZnCoordinationConstraint::ZnCoordinationConstraint( ZnCoordinationScorerCOP zn_score
):
	core::scoring::constraints::Constraint( core::scoring::metalhash_constraint ),
	zn_score_( zn_score )
{
	assert( zn_score_ ); // do not create this with a null pointer
}

ZnCoordinationConstraint::~ZnCoordinationConstraint() {}

core::scoring::constraints::ConstraintOP
ZnCoordinationConstraint::clone() const {
	return core::scoring::constraints::ConstraintOP( new ZnCoordinationConstraint( zn_score_ ) );
}

/// @brief Returns the number of atoms involved in defining this constraint.
core::Size
ZnCoordinationConstraint::natoms() const {
	return 7;
}

/// @brief Returns the AtomID referred to by index.
core::id::AtomID const &
ZnCoordinationConstraint::atom( Size const index ) const
{
	assert( index <= natoms() );
	if ( index <= 3 ) {
		return zn_score_->asymm_atids()[ index ];
	} else if (index <= 6 ) {
		return zn_score_->focused_clone_atids()[ index-3 ];
	} else {
		// The seventh atom qualifies this as a non-pairwise decomposable score term
		return zn_score_->atom_one_on_third_residue();
	}
}

void
ZnCoordinationConstraint::score(
	core::scoring::func::XYZ_Func const & xyz_func,
	core::scoring::EnergyMap const & /*weights*/,
	core::scoring::EnergyMap & emap
) const
{
	emap[ core::scoring::metalhash_constraint ] += zn_score_->score( xyz_func.residue(r1()),  xyz_func.residue(r2()) );
}

/// @brief Fill the f1 and f2 vectors, necessary for considering the
/// derivative this constraint during minimization. (someone please reference
/// Bill Wedermeyer's paper here, as I'm in an airport and can't fill it in
/// myself!)
void
ZnCoordinationConstraint::fill_f1_f2(
	core::id::AtomID const & /*atom*/,
	core::scoring::func::XYZ_Func const & /*xyz_func*/,
	core::Vector & /*F1*/,
	core::Vector & /*F2*/,
	core::scoring::EnergyMap const & /*weights*/
) const
{
}


}
}
