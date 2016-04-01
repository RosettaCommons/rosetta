// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/VdWTinkerPotential.cc
/// @brief
/// @author Jim Havranek

// Project headers
#include <core/scoring/VdWTinkerPotential.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/prof.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>

#include <numeric/constants.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.io.hh>
#include <numeric/xyz.functions.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <boost/lexical_cast.hpp>

#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <map>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.VdWTinkerPotential" );

typedef numeric::xyzVector< core::Real > Vector;
typedef numeric::xyzMatrix< core::Real > Matrix;

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {

VdWTinkerResidueInfo::~VdWTinkerResidueInfo() {}

bool
VdWShouldItCount(
	conformation::Residue const & rsd,
	Size const & atm
) {
	if ( rsd.seqpos() == 2 ) return false;

	if ( rsd.atom_name( atm ) == " CB " ) {
		return true;
	} else {
		return false;
	}
}

///
void
VdWTinkerResidueInfo::initialize( conformation::Residue const & rsd )
{
	Size const natoms( rsd.natoms() );
	vdw_type_.resize( natoms );

	// Look up values by Amoeba type
	for ( Size i=1; i<= natoms; ++i ) {
		vdw_type_[ i ] = 0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
VdWTinkerPoseInfo::VdWTinkerPoseInfo( VdWTinkerPoseInfo const & src ):
	CacheableData()
{
	Size const src_size( src.size() );

	residue_info_.resize( src_size );
	placeholder_residue_.resize( src_size );
	placeholder_info_.resize( src_size );

	for ( Size i=1; i<= src_size; ++i ) {
		residue_info_[i] = src.residue_info_[i]->copy_clone();
		if ( src.placeholder_residue_[i] ) {
			placeholder_residue_[i] = src.placeholder_residue_[i]->clone();
			placeholder_info_[i] = src.placeholder_info_[i]->copy_clone();
		} else {
			placeholder_residue_[i] = 0;
			placeholder_info_[i] = 0;
		}
	}
	being_packed_ = src.being_packed_;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
void
VdWTinkerPoseInfo::initialize( pose::Pose const & pose )
{
	Size const nres( pose.total_residue() );

	residue_info_.resize( nres, 0 );
	placeholder_residue_.resize( nres, 0 );
	placeholder_info_.resize( nres, 0 );

	for ( Size i=1; i<= nres; ++i ) {
		if ( !residue_info_[i] ) residue_info_[i] = VdWTinkerResidueInfoOP( new VdWTinkerResidueInfo( pose.residue(i) ) );
		else  residue_info_[i]->initialize( pose.residue(i) );
	}

	being_packed_.clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
VdWTinkerPoseInfo::set_placeholder( Size const i, ResidueOP rsd, VdWTinkerResidueInfoOP info )
{
	placeholder_residue_[ i ] = rsd;
	placeholder_info_[ i ] = info;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
VdWTinkerPoseInfo::set_repack_list( utility::vector1< bool > const & repacking_residues )
{
	being_packed_.resize( size(), false );
	for ( Size i=1; i<= size(); ++i ) {
		being_packed_[i] = repacking_residues[ i ];
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
VdWTinkerRotamerSetInfo::initialize( RotamerSetBase const & rotamer_set )
{
	Size const nrot( rotamer_set.num_rotamers() );
	residue_info_.resize( nrot );
	for ( Size i=1; i<= nrot; ++i ) {
		residue_info_[i] = VdWTinkerResidueInfoOP( new VdWTinkerResidueInfo( *rotamer_set.rotamer(i) ) );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
VdWTinkerPotential::read_in_amoeba_parameters() {

	/////////////////////////////////////
	// First read in vdw type parameters
	/////////////////////////////////////
	std::string type_file_name = basic::database::full_name( "chemical/amoeba/amoeba09_types.txt" );
	std::ifstream types_file( type_file_name.c_str() );
	std::string input_line;
	while ( getline( types_file, input_line ) ) {
		//TR << "Processing line: " << input_line << std::endl;
		utility::vector1< std::string > tokens( utility::split_whitespace( input_line ) );
		Size const num_tokens( tokens.size() );
		//TR << "Found " << num_tokens << " tokens " << std::endl;
		std::string atomname( tokens[2] );
		utility::trim( atomname, " " );
		//TR << "Atom name " << atomname << " Residue name " << tokens[3] << std::endl;
		std::string new_key = atomname + tokens[3];
		type_lookup_[ new_key ] = static_cast<core::Size>( boost::lexical_cast< core::Size >( tokens[ num_tokens ] ) );
		//TR << "Stashing key X" << new_key << "X" << std::endl;
	}
	types_file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
VdWTinkerPotential::read_in_vdw_tinker_parameters() {

	std::string vdw_file_name = basic::database::full_name( "chemical/amoeba/amoeba09_vdw_params.txt" );
	std::ifstream vdw_file( vdw_file_name.c_str() );
	std::string input_line;

	// Find the maximum vdw type number
	while ( getline( vdw_file, input_line ) ) {

		//TR << "Processing line: " << input_line << std::endl;

		utility::vector1< std::string > tokens( utility::split_whitespace( input_line ) );
		Size const num_tokens( tokens.size() );

		//  Size const vdw_type( static_cast<core::Size>( boost::lexical_cast< core::Size >( tokens[ 2 ] ) ) );
		Real const vdw_radius( static_cast<core::Real>( boost::lexical_cast< core::Real >( tokens[ 3 ] ) ) );
		Real const vdw_depth( static_cast<core::Real>( boost::lexical_cast< core::Real >( tokens[ 4 ] ) ) );
		Real const vdw_reduce( num_tokens == 5 ? static_cast<core::Real>( boost::lexical_cast< core::Real >( tokens[ 5 ] ) ) : 1.0 );

		//TR << "Processed:  type " << vdw_type << " has radius " << vdw_radius << " with well depth " << vdw_depth << " and reduction factor " << vdw_reduce << std::endl;

		vdw_radius_.push_back( vdw_radius );
		vdw_depth_.push_back( vdw_depth );
		vdw_reduce_.push_back( vdw_reduce );
	}

	vdw_file.close();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

core::Size
VdWTinkerPotential::amoeba_type_lookup(
	std::string const & atomname,
	std::string const & resname,
	std::string const & variantname
) const {

	std::string const not_variant( "NONE" );
	core::Size type( 0 );

	// Lookup default value w/o variant info
	std::string tmp_key = atomname + resname;
	tmp_key.erase( std::remove_if( tmp_key.begin(), tmp_key.end(), ::isspace ), tmp_key.end() );
	std::map< std::string, core::Size >::const_iterator map_it;

	//TR << "Querying key X" << tmp_key << "X" << std::endl;
	map_it = type_lookup_.find( tmp_key );
	if ( map_it != type_lookup_.end() ) {
		type = map_it->second;
	}

	if ( !utility::trimmed_compare( variantname, not_variant ) ) {
		std::string tmp_key2 = atomname + resname + variantname;
		tmp_key2.erase( std::remove_if( tmp_key2.begin(), tmp_key2.end(), ::isspace ), tmp_key2.end() );
		// Only overwrite default if an entry is found

		//TR << "Querying key X" << tmp_key2 << "X" << std::endl;
		map_it = type_lookup_.find( tmp_key2 );
		if ( map_it != type_lookup_.end() ) {
			type = map_it->second;
		}
	}

	// Give a holler if nothing has been found

	if ( type == 0 ) {
		TR << "PROBLEM - NO VDW TYPE FOUND FOR " << atomname << " " << resname << " " << variantname << std::endl;
	}

	return type;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
VdWTinkerPotential::assign_residue_amoeba_type(
	Residue const & rsd,
	VdWTinkerResidueInfo & vdw
) const {
	utility::vector1< std::string > parsed_resname( utility::string_split_simple( rsd.name(), ':' ) );
	std::string resname( parsed_resname[1] );

	// Handle variants
	std::string variantname( "NONE" );
	if ( parsed_resname.size() == 2 ) {
		//TR << "Using variant specialization " << parsed_resname[2] << std::endl;
		variantname = parsed_resname[2];
	} else if ( parsed_resname.size() > 2 ) {
		//   TR << "MULTIPLE VARIANT SITUATION -> NOT YET CODED!!!" << std::endl;
		//   TR << "Full resname is " << rsd.name() <<  std::endl;
		//   TR << "Using parsed resname " << parsed_resname[2] <<  std::endl;
		//   TR << "Using first variant name only!!!" << std::endl;
		variantname = parsed_resname[2];
		//for( Size ivar = 1 ; ivar <= parsed_resname.size() ; ++ivar ) {
		//TR << "Residue info " << parsed_resname[ ivar ] << std::endl;
		//}
	}

	for ( Size j = 1 ; j <= rsd.natoms() ; j++ ) {
		//TR << "Residue " << resname << " Atom " << rsd.atom_name( j )   << std::endl;
		// Lookup amoeba type
		core::Size this_type( amoeba_type_lookup( rsd.atom_name(j), resname, variantname ) );
		//TR << "Found Amoeba VdW type " << this_type << std::endl;
		vdw.nonconst_vdw_type( j ) = this_type;
	}
}

void
VdWTinkerPotential::assign_all_amoeba_types(
	pose::Pose & pose
) const {
	//TR << "Assigning amoeba vdw types!" << std::endl;

	Size const nres( pose.total_residue() );

	// Look up the Amoeba vdw type information
	VdWTinkerPoseInfoOP vdw_info;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) ) {
		vdw_info = utility::pointer::static_pointer_cast< VdWTinkerPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) );
	} else {
		vdw_info = VdWTinkerPoseInfoOP( new VdWTinkerPoseInfo() );
	}

	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd( pose.residue( res1 ) );
		VdWTinkerResidueInfo & vdw1( vdw_info->residue_info( res1 ) );
		assign_residue_amoeba_type( rsd, vdw1 );
	}

	pose.data().set( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO, vdw_info );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
VdWTinkerPotential::setup_for_scoring(
	pose::Pose & pose
) const {
	// TR << "in setup_for_scoring use_nblist is " << pose.energies().use_nblist() << std::endl;

	VdWTinkerPoseInfoOP vdw_info;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) ) {
		vdw_info = utility::pointer::static_pointer_cast< VdWTinkerPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) );
	} else {
		//TR << "Allocating fresh vdw info objects" << std::endl;
		vdw_info = VdWTinkerPoseInfoOP( new VdWTinkerPoseInfo() );
		vdw_info->initialize( pose );
		pose.data().set( pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO, vdw_info );
		assign_all_amoeba_types( pose );
	}

	// TR << "In setup_for_scoring, past initialization" << std::endl;

	// These shouldn't be necessary, but I don't know how to
	// transfer from RotamerSetInfo to ResidueInfo
	// The lines with the extra indent are only necessary because
	// RotamerSetInfo isn't automatically transferred to ResidueInfo
	// after a repack.

	vdw_info->initialize( pose );
	assign_all_amoeba_types( pose );
	// TR << "Exiting setup_for_scoring" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief
/// Note: when called at the beginning of rotamer_trials, task.being_packed(i) will be false for all i
/// this ensures that we use all the information we have to compute the current set of radii

void
VdWTinkerPotential::setup_for_packing(
	pose::Pose & , // ellided
	utility::vector1< bool > const & // repacking_residues
) const {
	// jjh Commenting out for now.  First need to get regular scoring done.  Worry
	// about packing later.
#ifdef NOTDEF
	PROF_START( basic::GB_SETUP_FOR_PACKING );

	VdWTinkerPoseInfoOP vdw_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) ) {
		vdw_info = utility::pointer::static_pointer_cast< VdWTinkerPoseInfo >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) );
	} else {
		vdw_info = new VdWTinkerPoseInfo();
	}

	//jjh zero out arrays
	vdw_info->initialize( pose );

	/// store info about which positions are moving
	vdw_info->set_repack_list( repacking_residues );
	build_placeholders( pose, *vdw_info );
	get_template_born_radii( pose, *vdw_info );

	pose.data().set( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO, vdw_info );

	PROF_STOP( basic::GB_SETUP_FOR_PACKING );
#endif
}



void
VdWTinkerPotential::get_rotamers_vdw_info(
	core::pose::Pose const &,
	conformation::RotamerSetBase & rotamer_set
) const {
	VdWTinkerRotamerSetInfoOP vdw_info_rotamers( new VdWTinkerRotamerSetInfo( rotamer_set ) );
	for ( Size n=1; n<= rotamer_set.num_rotamers(); ++n ) {
		assign_residue_amoeba_type( *rotamer_set.rotamer(n), vdw_info_rotamers->residue_info( n ) );
	}

	rotamer_set.data().set( core::conformation::RotamerSetCacheableDataType::VDWTINKER_ROTAMER_SET_INFO, vdw_info_rotamers );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// called eg after a rotamer substitution is accepted during rotamer trials
void
VdWTinkerPotential::update_residue_for_packing(
	pose::Pose & pose,
	Size const seqpos
) const {
	VdWTinkerPoseInfo & vdw_info( static_cast< VdWTinkerPoseInfo & >( pose.data().get( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO ) ) );
	VdWTinkerResidueInfo & vdw_residue_info( vdw_info.residue_info( seqpos ) );

	Residue const & rsd( pose.residue( seqpos ) );

	vdw_residue_info.initialize( rsd );
	assign_residue_amoeba_type( rsd, vdw_residue_info );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
VdWTinkerPotential::get_res_res_vdw(
	Residue const & rsd1,
	VdWTinkerResidueInfo const & vdw1,
	Residue const & rsd2,
	VdWTinkerResidueInfo const & vdw2
) const {
	using namespace etable::count_pair;

	Size natoms1 = rsd1.natoms();
	Size natoms2 = rsd2.natoms();

	Size path_dist( 1 );
	Real weight( 0.0 );
	Real vdwE( 0.0 );

	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );

	etable::count_pair::CountPairFunctionOP cpfxn( 0 );

	if ( same_res ) {
		cpfxn = etable::count_pair::CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
	} else {
		cpfxn = etable::count_pair::CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	}

	for ( Size atm1 = 1 ; atm1 <= natoms1 ; ++atm1 ) {
		if ( rsd1.is_virtual( atm1 ) ) continue;

		//  if( !VdWShouldItCount( rsd1, atm1 ) ) continue;

		Real const rad1( vdw_radius_[ vdw1.vdw_type( atm1 ) ] );
		Real const dep1( vdw_depth_[ vdw1.vdw_type( atm1 ) ] );
		Real const red1( vdw_reduce_[ vdw1.vdw_type( atm1 ) ] );

		// Hydrogens are shifted towards their heavy atom base
		Vector p1;
		if ( rsd1.atom_is_hydrogen( atm1 ) ) {
			Vector delta( rsd1.xyz( atm1 ) - rsd1.xyz( rsd1.atom_base( atm1 ) ) );
			delta *= red1;
			p1 = rsd1.xyz( rsd1.atom_base( atm1 ) ) + delta;
		} else {
			p1 = rsd1.xyz( atm1 );
		}

		for ( Size atm2 = (same_res ? atm1 : 1 ) ; atm2 <= natoms2 ; ++atm2 ) {
			if ( rsd2.is_virtual( atm2 ) ) continue;

			Real atom_atomE(0.0);

			//   if( !VdWShouldItCount( rsd2, atm2 ) ) continue;

			Real const rad2( vdw_radius_[ vdw2.vdw_type( atm2 ) ] );
			Real const dep2( vdw_depth_[ vdw2.vdw_type( atm2 ) ] );
			Real const red2( vdw_reduce_[ vdw2.vdw_type( atm2 ) ] );

			// Hydrogens are shifted towards their heavy atom base
			Vector p2;
			if ( rsd2.atom_is_hydrogen( atm2 ) ) {
				Vector delta( rsd2.xyz( atm2 ) - rsd2.xyz( rsd2.atom_base( atm2 ) ) );
				delta *= red2;
				p2 = rsd2.xyz( rsd2.atom_base( atm2 ) ) + delta;
			} else {
				p2 = rsd2.xyz( atm2 );
			}

			// Assuming all the checks tell us to calculate the interaction

			if ( !cpfxn->count( atm1, atm2, weight, path_dist ) &&
				(!same_res || (atm1 != atm2) ) ) continue;
			
			Real const dist( p1.distance( p2 ) );
			Real const eff_dep( 4.0*dep1*dep2 / std::pow( std::sqrt( dep1 ) + std::sqrt( dep2 ), 2.0 ) );
			Real const eff_rad( ( rad1*rad1*rad1 + rad2*rad2*rad2 ) / ( rad1*rad1 + rad2*rad2 ) );
			Real const rad_ratio( dist / eff_rad );
			
			Real const factor1( std::pow( 1.07/( rad_ratio + 0.07 ), 7.0 ) );
			Real const factor2( 1.12/( std::pow( rad_ratio, 7.0 ) + 0.12 ) - 2.0 );
			
			
			atom_atomE = eff_dep*factor1*factor2;

			vdwE += atom_atomE;
		}
	}
	// TR << "res-res energy between " << rsd1.seqpos() << " and " << rsd2.seqpos() << " is " << vdwE << std::endl;
	// TR << "res-res interaction between residue " << rsd1.seqpos() << " and " << rsd2.seqpos() << std::endl;
	// TR << "Returning vdwE of " << vdwE << std::endl;

	return vdwE;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
VdWTinkerPotential::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	VdWTinkerResidueInfo const &,
	VdWTinkerResidueInfo const &,
	pose::Pose const & pose, // provides context
	Real const & factor,
	utility::vector1< DerivVectorPair > & r1_at_derivs,
	utility::vector1< DerivVectorPair > & r2_at_derivs
) const {
	using namespace etable::count_pair;

	Size natoms1 = rsd1.natoms();
	Size natoms2 = rsd2.natoms();

	Size path_dist( 1 );
	Real weight( 0.0 );
	
	Size const resi( rsd1.seqpos() );
	Size const resj( rsd2.seqpos() );
	bool const same_res( resi == resj );

	VdWTinkerPoseInfo const & vdw_info( static_cast< VdWTinkerPoseInfo const & >( pose.data().get( core::pose::datacache::CacheableDataType::VDWTINKER_POSE_INFO)));

	VdWTinkerResidueInfo const & vdw1( vdw_info.residue_info( resi ) );
	VdWTinkerResidueInfo const & vdw2( vdw_info.residue_info( resj ) );

	assert( pose.energies().use_nblist() );
	assert( r1_at_derivs.size() >= rsd1.natoms() );
	assert( r2_at_derivs.size() >= rsd2.natoms() );

	CountPairFunctionOP cpfxn( 0 );
	if ( same_res ) {
		cpfxn = etable::count_pair::CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
	} else {
		cpfxn = etable::count_pair::CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	}

	for ( Size atm1 = 1 ; atm1 <= natoms1 ; ++atm1 ) {
		if ( rsd1.is_virtual( atm1 ) ) continue;

		//  if( !VdWShouldItCount( rsd1, atm1 ) ) continue;

		Real const rad1( vdw_radius_[ vdw1.vdw_type( atm1 ) ] );
		Real const dep1( vdw_depth_[ vdw1.vdw_type( atm1 ) ] );
		Real const red1( vdw_reduce_[ vdw1.vdw_type( atm1 ) ] );

		// Hydrogens are shifted towards their heavy atom base
		Vector p1;
		if ( rsd1.atom_is_hydrogen( atm1 ) ) {
			Vector delta( rsd1.xyz( atm1 ) - rsd1.xyz( rsd1.atom_base( atm1 ) ) );
			delta *= red1;
			p1 = rsd1.xyz( rsd1.atom_base( atm1 ) ) + delta;
		} else {
			p1 = rsd1.xyz( atm1 );
		}

		for ( Size atm2 = (same_res ? atm1 : 1 ) ; atm2 <= natoms2 ; ++atm2 ) {
			if ( rsd2.is_virtual( atm2 ) ) continue;

			//   if( !VdWShouldItCount( rsd2, atm2 ) ) continue;

			Real const rad2( vdw_radius_[ vdw2.vdw_type( atm2 ) ] );
			Real const dep2( vdw_depth_[ vdw2.vdw_type( atm2 ) ] );
			Real const red2( vdw_reduce_[ vdw2.vdw_type( atm2 ) ] );

			// Hydrogens are shifted towards their heavy atom base
			Vector p2;
			if ( rsd2.atom_is_hydrogen( atm2 ) ) {
				Vector delta( rsd2.xyz( atm2 ) - rsd2.xyz( rsd2.atom_base( atm2 ) ) );
				delta *= red2;
				p2 = rsd2.xyz( rsd2.atom_base( atm2 ) ) + delta;
			} else {
				p2 = rsd2.xyz( atm2 );
			}

			// Assuming all the checks tell us to calculate the interaction

			if ( !cpfxn->count( atm1, atm2, weight, path_dist ) &&
				(!same_res || (atm1 != atm2) ) ) continue;
			
			Real const dist( p1.distance( p2 ) );
			Real const eff_dep( 4.0*dep1*dep2 / std::pow( std::sqrt( dep1 ) + std::sqrt( dep2 ), 2.0 ) );
			Real const eff_rad( ( rad1*rad1*rad1 + rad2*rad2*rad2 ) / ( rad1*rad1 + rad2*rad2 ) );
			Real const rad_ratio( dist / eff_rad );
			
			Real const factor1( std::pow( 1.07/( rad_ratio + 0.07 ), 7.0 ) );
			Real const factor2( 1.12/( std::pow( rad_ratio, 7.0 ) + 0.12 ) - 2.0 );
			
			Real const dfactor1( -6.542*std::pow( 1.07/(rad_ratio+0.07), 8.0 ) );
			Real const dfactor2( -7.84*std::pow(rad_ratio,6.0)/( std::pow( ( std::pow(rad_ratio,7.0) + 0.12 ) , 2.0 ) ) );
			
			Real const dEdr( -1.0*eff_dep*factor*( factor1*dfactor2 + factor2*dfactor1 )/eff_rad );
			
			Vector const deriv_dr_f1( p2.cross_product( p1 ) );
			Vector const deriv_dr_f2( p2 - p1 );
			
			Vector const f1( dEdr * deriv_dr_f1 / dist );
			Vector const f2( dEdr * deriv_dr_f2 / dist );
			
			r1_at_derivs[ atm1 ].f1() += f1;
			r1_at_derivs[ atm1 ].f2() += f2;
			r2_at_derivs[ atm2 ].f1() -= f1;
			r2_at_derivs[ atm2 ].f2() -= f2;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace scoring
} // namespace core
