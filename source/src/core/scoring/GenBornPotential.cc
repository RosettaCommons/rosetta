// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/GenBornPotential.cc
/// @brief
/// @author JIM HAVRANEK (originally)

// Project headers
#include <core/scoring/GenBornPotential.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/types.hh>

#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/conformation/RotamerSetCacheableDataType.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/prof.hh>

// Utility headers
#include <utility/exit.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/id/AtomID.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.hh>


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
/**

This is a reimplementation of Jim Havranek's original rosetta++ Gen Born code.
source files: rosetta++/gb_elec*

**/
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

GenBornResidueInfo::~GenBornResidueInfo() = default;


void
GenBornResidueInfo::initialize( conformation::Residue const & rsd )
{
	Size const natoms( rsd.natoms() );
	born_radius_.resize( natoms );
	atomic_radius_.resize( natoms );
	scale_factor_.resize( natoms );

	// important: initialize with 0.0 prior to burial calculations
	std::fill( born_radius_.begin(), born_radius_.end(), 0.0 );

	// atom index for looking up an extra data type stored in the AtomTypes
	Size const       GB_RADIUS_INDEX( rsd.atom_type_set().extra_parameter_index( "GB_RADIUS"       ) );
	Size const GB_SCALE_FACTOR_INDEX( rsd.atom_type_set().extra_parameter_index( "GB_SCALE_FACTOR" ) );

	for ( Size i=1; i<= natoms; ++i ) {
		atomic_radius_[i] = rsd.atom_type(i).extra_parameter( GB_RADIUS_INDEX );
		scale_factor_ [i] = rsd.atom_type(i).extra_parameter( GB_SCALE_FACTOR_INDEX );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
GenBornPoseInfo::GenBornPoseInfo( GenBornPoseInfo const & src ):
	CacheableData()
{
	Size const src_size( src.size() );

	residue_info_.resize( src_size );
	placeholder_residue_.resize( src_size );
	placeholder_info_.resize( src_size );

	for ( Size i=1; i<= src_size; ++i ) {
		residue_info_[i] = src.residue_info_[i]->clone();
		if ( src.placeholder_residue_[i] ) {
			placeholder_residue_[i] = src.placeholder_residue_[i]->clone();
			placeholder_info_[i] = src.placeholder_info_[i]->clone();
		} else {
			placeholder_residue_[i] = nullptr;
			placeholder_info_[i] = nullptr;
		}
	}
	being_packed_ = src.being_packed_;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
void
GenBornPoseInfo::initialize( pose::Pose const & pose )
{
	Size const nres( pose.total_residue() );

	residue_info_.resize( nres, nullptr );
	placeholder_residue_.resize( nres, nullptr );
	placeholder_info_.resize( nres, nullptr );

	for ( Size i=1; i<= nres; ++i ) {
		if ( !residue_info_[i] ) residue_info_[i] = GenBornResidueInfoOP( new GenBornResidueInfo( pose.residue(i) ) );
		else  residue_info_[i]->initialize( pose.residue(i) );
	}

	being_packed_.clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
GenBornPoseInfo::set_placeholder( Size const i, ResidueOP rsd, GenBornResidueInfoOP info )
{
	placeholder_residue_[ i ] = rsd;
	placeholder_info_[ i ] = info;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
GenBornPoseInfo::set_repack_list( utility::vector1< bool > const & repacking_residues )
{
	being_packed_.resize( size(), false );
	for ( Size i=1; i<= size(); ++i ) {
		being_packed_[i] = repacking_residues[ i ];
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// dont forget to 0 the born_radii
void
GenBornRotamerSetInfo::initialize( RotamerSetBase const & rotamer_set )
{
	Size const nrot( rotamer_set.num_rotamers() );
	residue_info_.resize( nrot );
	for ( Size i=1; i<= nrot; ++i ) {
		residue_info_[i] = GenBornResidueInfoOP( new GenBornResidueInfo( *rotamer_set.rotamer(i) ) );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
GenBornPotential::res_res_burial(
	Residue const & rsd1,
	GenBornResidueInfo & gb1,
	Residue const & rsd2,
	GenBornResidueInfo const & gb2
) const {
	bool const same_res( rsd1.seqpos() == rsd2.seqpos() );

	//jjh local
	utility::vector1< Real > atm_radii1( rsd1.natoms() );
	utility::vector1< Real > atm_radii2( rsd2.natoms() );

	//jjh Fill the radii array - later will allow for overrides
	for ( int atm1 = 1, atm1e = rsd1.natoms(); atm1 <= atm1e; ++atm1 ) {
		atm_radii1[ atm1 ] = gb1.atomic_radius( atm1 );
	}
	for ( int atm2 = 1, atm2e = rsd2.natoms(); atm2 <= atm2e; ++atm2 ) {
		atm_radii2[ atm2]  = gb2.scale_factor( atm2 ) * ( gb2.atomic_radius( atm2 ) - ParamS );
	}

	for ( int atm1 = 1, atm1e = rsd1.natoms(); atm1 <= atm1e; ++atm1 ) {
		Real const rad1 = atm_radii1[ atm1 ];
		Real const rwork1 = rad1 - ParamS;
		Real const inv_r1 = 1.0/rwork1;
		Vector const & xyz1( rsd1.xyz( atm1 ) );

		for ( int atm2 = 1, atm2e = rsd2.natoms(); atm2 <= atm2e; ++atm2 ) {
			Real const rwork2 = atm_radii2[ atm2 ];
			Real const dis2( xyz1.distance_squared( rsd2.xyz( atm2 ) ) );
			Real const dis = std::sqrt(dis2);
			Real const inv_dis = 1.0/(dis+1.0e-6);
			Real const rwork2_sq = rwork2*rwork2;

			if ( same_res && (atm1 == atm2) ) continue;

			if ( dis > (3.5 * rwork2) ) {
				Real const inv_dis2 = inv_dis * inv_dis;
				Real const tmpsd = rwork2_sq * inv_dis2;
				Real const dumbo = Param_TA+tmpsd*(Param_TB+tmpsd*(Param_TC+tmpsd*(Param_TD+tmpsd*Param_TDD)));
				gb1.born_radius( atm1 ) -= rwork2_sq*rwork2*inv_dis2*inv_dis2*dumbo;

			} else if ( dis > (rwork1 + rwork2) ) {
				gb1.born_radius( atm1 ) -= (0.5*(rwork2/(dis2 - rwork2_sq) +
					0.5*inv_dis*std::log((dis-rwork2)/(dis+rwork2))));

			} else if ( dis > std::abs( rwork1 - rwork2 ) ) {
				Real const theta = 0.5*inv_r1*inv_dis*(dis2+rwork1*rwork1 - rwork2_sq);
				Real const U12 = 1.0/(dis+rwork2);
				gb1.born_radius( atm1 ) -= 0.25*(inv_r1 * (2.0-theta) - U12 + inv_dis*std::log(rwork1*U12));

			} else if ( rwork1 < rwork2 ) {
				gb1.born_radius( atm1 ) -= 0.5*(rwork2/(dis2-rwork2_sq) + 2.0*inv_r1 +
					0.5*inv_dis*std::log( (rwork2-dis)/(rwork2+dis)));

			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// helper function

void
GenBornPotential::finalize_born_radii( GenBornResidueInfo & gb_info ) const
{
	for ( Size atm=1; atm<= gb_info.size(); ++atm ) {
		Real const radius = gb_info.atomic_radius( atm );
		Real const factor = radius - ParamS;
		Real const br_save = gb_info.born_radius( atm );
		Real const integral = (-1.0)*factor*br_save;
		Real const inv_brad = (1.0/factor) - std::tanh( ( ParamD - ParamB*integral +
			ParamG*integral*integral)*integral)/radius;
		gb_info.born_radius( atm ) = 1.0 / inv_brad;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void
GenBornPotential::get_all_born_radii(
	pose::Pose & pose
) const {
	PROF_START( basic::GB_GET_ALL_BORN_RADII );

	Size const nres( pose.total_residue() );

	// initialize gb info:
	//  * zeros born radii for start of calculations
	//  * fills the atomic radii  (from the Pose/AtomTypeSet)
	//  * fills the scale factors (from the Pose/AtomTypeSet)

	GenBornPoseInfoOP gb_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO ) ) {
		gb_info = utility::pointer::static_pointer_cast< core::scoring::GenBornPoseInfo > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO ) );
	} else {
		gb_info = GenBornPoseInfoOP( new GenBornPoseInfo() );
	}

	//jjh zero out arrays
	gb_info->initialize( pose );

	//jjh get the residue-residue burial
	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		Residue const & rsd1( pose.residue( res1 ) );
		GenBornResidueInfo & gb1( gb_info->residue_info( res1 ) );
		debug_assert( std::abs( gb1.born_radius(1) ) < 1e-3 ); // ensure that this has been zeroed out
		for ( Size res2 = 1; res2 <= nres; ++res2 ) {
			res_res_burial( rsd1, gb1, pose.residue( res2 ), gb_info->residue_info( res2 ) );
		}
	}

	//jjh convert to Born radius

	for ( Size res1 = 1; res1 <= nres; ++res1 ) {
		finalize_born_radii( gb_info->residue_info( res1 ) );
	}

	pose.data().set( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO, gb_info );

	PROF_STOP( basic::GB_GET_ALL_BORN_RADII );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief
/// Note: when called at the beginning of rotamer_trials, task.being_packed(i) will be false for all i
/// this ensures that we use all the information we have to compute the current set of radii

void
GenBornPotential::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const & repacking_residues
) const {
	PROF_START( basic::GB_SETUP_FOR_PACKING );

	GenBornPoseInfoOP gb_info;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO ) ) {
		gb_info = utility::pointer::static_pointer_cast< core::scoring::GenBornPoseInfo > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO ) );
	} else {
		gb_info = GenBornPoseInfoOP( new GenBornPoseInfo() );
	}

	//jjh zero out arrays
	gb_info->initialize( pose );

	/// store info about which positions are moving
	gb_info->set_repack_list( repacking_residues );

	build_placeholders( pose, *gb_info );
	get_template_born_radii( pose, *gb_info );

	pose.data().set( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO, gb_info );

	PROF_STOP( basic::GB_SETUP_FOR_PACKING );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// private
void
GenBornPotential::build_placeholders(
	pose::Pose const & pose,
	GenBornPoseInfo & gb_info
) const {
	Size const nres( pose.total_residue() );

	chemical::ResidueTypeSetCOP residue_set( pose.residue(1).residue_type_set() );

	for ( Size i=1; i<= nres; ++i ) {
		if ( !gb_info.being_packed(i) ) continue;
		Residue const & existing_rsd( pose.residue(i) );

		if ( !existing_rsd.is_protein() ) {
			std::cout << "WARNING: no mechanism for building genborn placeholders at non-protein positions\n" <<
				"Using existing residue coords" << std::endl;
			gb_info.set_placeholder( i, existing_rsd.clone(), GenBornResidueInfoOP( new GenBornResidueInfo( existing_rsd ) ) );
		} else {
			// build a placeholder at this position
			chemical::ResidueTypeCOP protein_placeholder_residue_type( residue_set->name_map("GB_AA_PLACEHOLDER").get_self_ptr() );
			// use appropriate termini variants if necessary:
			if ( existing_rsd.is_lower_terminus() ) {
				protein_placeholder_residue_type = chemical::ResidueTypeCOP(
					residue_set->get_residue_type_with_variant_added( *protein_placeholder_residue_type,
					chemical::LOWER_TERMINUS_VARIANT ).get_self_ptr() );
			}
			if ( existing_rsd.is_upper_terminus() ) {
				protein_placeholder_residue_type = chemical::ResidueTypeCOP(
					residue_set->get_residue_type_with_variant_added( *protein_placeholder_residue_type,
					chemical::UPPER_TERMINUS_VARIANT ).get_self_ptr() );
			}

			conformation::ResidueOP rsd
				( conformation::ResidueFactory::create_residue( *protein_placeholder_residue_type, existing_rsd,
				pose.conformation() ) );
			GenBornResidueInfoOP rsd_info( new GenBornResidueInfo( *rsd ) );

			Size const dummy_index( rsd->atom_index("DUMM") );
			rsd_info->atomic_radius( dummy_index ) = dummy_radius;
			rsd_info->scale_factor ( dummy_index ) = dummy_scale;

			// debug placement of dummy atom
			runtime_assert( std::abs( rsd->xyz("CA").distance( rsd->xyz( dummy_index )) - dummy_distance ) < 1e-2 );
			gb_info.set_placeholder( i, rsd, rsd_info );
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// private
void
GenBornPotential::get_template_born_radii(
	pose::Pose const & pose,
	//pack::task::PackerTask const & task,
	GenBornPoseInfo & gb_info
) const
{
	Size const nres( pose.total_residue() );

	debug_assert( gb_info.size() == nres );

	for ( Size i=1; i<= nres; ++i ) {
		if ( gb_info.being_packed( i ) ) continue;
		Residue const & rsd1( pose.residue( i ) );
		GenBornResidueInfo & gb1( gb_info.residue_info( i ) );
		debug_assert( rsd1.natoms()<1 || std::abs( gb1.born_radius(1) ) < 1e-3 ); // ensure radii have been initialized to 0
		for ( Size j=1; j<= nres; ++j ) {

			if ( gb_info.being_packed( j ) ) {
				// use placeholder
				res_res_burial( rsd1, gb1, gb_info.placeholder_residue( j ), gb_info.placeholder_info( j ) );
			} else {
				res_res_burial( rsd1, gb1, pose.residue(j), gb_info.residue_info( j ) );
			}
		}
	}

	for ( Size i=1; i<= nres; ++i ) {
		if ( gb_info.being_packed( i ) ) continue;
		finalize_born_radii( gb_info.residue_info(i) );
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// called eg after a rotamer substitution is accepted during rotamer trials
void
GenBornPotential::update_residue_for_packing(
	pose::Pose & pose,
	Size const seqpos
) const {
	GenBornPoseInfo & gb_info( static_cast< GenBornPoseInfo & >( pose.data().get( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO ) ) );
	GenBornResidueInfo & gb( gb_info.residue_info( seqpos ) );

	Residue const & rsd( pose.residue( seqpos ) );
	gb.initialize( rsd );
	get_single_rotamer_born_radii( rsd, pose, gb_info, gb );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// uses placeholder info at positions i with gb_info.being_packed(i) == true


void
GenBornPotential::get_single_rotamer_born_radii(
	Residue const & rsd1,
	pose::Pose const & pose,
	GenBornPoseInfo const & gb_info,
	GenBornResidueInfo & gb1
) const {
	debug_assert( rsd1.natoms()<1 || std::abs( gb1.born_radius(1) ) < 1e-3 ); // ensure radii have been initialized to 0

	for ( Size res2=1; res2<= pose.total_residue(); ++res2 ) {
		if ( gb_info.being_packed( res2 ) ) {
			// use placeholder
			res_res_burial( rsd1, gb1, gb_info.placeholder_residue( res2 ), gb_info.placeholder_info( res2));
		} else {
			res_res_burial( rsd1, gb1, pose.residue( res2 ), gb_info.residue_info( res2 ) );
		}
	}

	finalize_born_radii( gb1 );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
GenBornPotential::get_rotamers_born_radii(
	pose::Pose const & pose,
	conformation::RotamerSetBase & rotamer_set
) const {
	using core::conformation::RotamerSetCacheableDataType::GEN_BORN_ROTAMER_SET_INFO;

	// this holds placeholders, info for non-packed residues
	GenBornPoseInfo const & gb_info_pose
		( static_cast< GenBornPoseInfo const & >(pose.data().get( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO )));

	// this will get cached in the rotamer set
	// this call should initialize the residue_info objects with the appropriate Residue info
	GenBornRotamerSetInfoOP gb_info_rotamers( new GenBornRotamerSetInfo( rotamer_set ) );

	for ( Size n=1; n<= rotamer_set.num_rotamers(); ++n ) {
		get_single_rotamer_born_radii( *rotamer_set.rotamer(n), pose, gb_info_pose, gb_info_rotamers->residue_info( n ) );
	}

	rotamer_set.data().set( GEN_BORN_ROTAMER_SET_INFO, gb_info_rotamers );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
GenBornPotential::get_res_res_elecE(
	Residue const & rsd1,
	GenBornResidueInfo const & gb1,
	Residue const & rsd2,
	GenBornResidueInfo const & gb2
) const {
	using namespace etable::count_pair;

	int natoms1 = rsd1.natoms();
	int natoms2 = rsd2.natoms();

	Real const inv_Ep = 1.0/Ep;
	Real const tau = (1.0/Ep - 1.0/Ew);
	bool const same_res = ( rsd1.seqpos() == rsd2.seqpos() );

	etable::count_pair::CountPairFunctionOP cpfxn( nullptr );
	if ( same_res ) {
		cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
	} else if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
	}

	Real elecE = 0.0;
	for ( int atm1 = 1; atm1 <= natoms1; ++atm1 ) {
		Real const q1 = rsd1.atomic_charge( atm1 );

		if ( std::fabs( q1 ) < 0.00001 ) continue;

		Real const brad1( gb1.born_radius( atm1 ) );
		Real const r1 = gb1.atomic_radius( atm1 );
		Vector const & xyz1( rsd1.xyz( atm1 ) );

		for ( int atm2 = 1; atm2 <= natoms2; ++atm2 ) {
			Real const q2 = rsd2.atomic_charge( atm2 );
			if ( std::fabs( q2 ) < 0.0001 ) continue;

			Real const brad2( gb2.born_radius( atm2 ) );

			Real const dis2( xyz1.distance_squared( rsd2.xyz( atm2 ) ) );
			Real const exparg = (-dis2)/(4.0 * brad1 * brad2 );
			Real const denom = std::sqrt( dis2 + brad1 * brad2 * std::exp( exparg ) );
			Real this_intx;
			if ( same_res ) {
				// tau = (1.0/Ep - 1.0/Ew)    Ep = 4.0  Ew=80.0
				this_intx = (-166.0)*tau*q1*q2/denom;
			} else {
				this_intx = (-332.0)*tau*q1*q2/denom;
			}
			elecE += this_intx;
			//      total_gb += this_intx;

			this_intx = 0.0;
			Real weight; // unused
			Size path_dist( 0 );
			if ( !cpfxn || cpfxn->count( atm1, atm2, weight, path_dist ) ) {
				// 3 or more bonds away
				Real const dis = std::sqrt( dis2 );
				Real const r2 = gb2.atomic_radius( atm2 );
				this_intx = 166.0 * gb_shell_intxn( inv_Ep * q1, r1, q2, r2, dis);
				if ( !same_res ) this_intx *= 2.0;
			}
			elecE += this_intx;
		}
	}

	return elecE;
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Calculates the interaction energy of two shells of charge.
///      Doesn't blow up as shells pass through each other
///
/// @author jjh 5/17/2004
///
/////////////////////////////////////////////////////////////////////////////////
Real
GenBornPotential::gb_shell_intxn(
	Real const qai,
	Real const rai,
	Real const qbi,
	Real const rbi,
	Real const dist
) const {
	if ( dist >= (rai+rbi) ) return (qai * qbi /dist);

	Real qa;
	Real ra;
	Real qb;
	Real rb;

	// Make sure rb is larger than ra
	if ( rai > rbi ) {
		qa = qbi;
		ra = rbi;
		qb = qai;
		rb = rai;
	} else {
		qa = qai;
		ra = rai;
		qb = qbi;
		rb = rbi;
	}

	if ( ((ra+rb) > dist) && (dist > (rb - ra)) ) {
		Real fout;
		if ( ra == rb ) {
			fout = 0.5 * (1.0+0.5*dist/ra);
		} else {
			fout = 0.5*(1.0 + 0.5*(ra*ra - rb*rb)/(ra*dist) + 0.5*dist/ra);
		}
		Real fin = 1.0 - fout;
		return qa*qb*(fin/rb + fout*2.0/(dist+ra+rb));
	} else {
		return qa*qb/rb;
	}
}

///////////////////////////////////////////////////////////////////////////////
Real
GenBornPotential::gb_shell_intxn_deriv(
	Real const qai,
	Real const rai,
	Real const qbi,
	Real const rbi,
	Real const dist
) const {
	if ( dist >= (rai+rbi) ) return (-1.0 * qai * qbi / ( dist * dist ) );

	Real qa;
	Real ra;
	Real qb;
	Real rb;

	// Make sure rb is larger than ra
	if ( rai > rbi ) {
		qa = qbi;
		ra = rbi;
		qb = qai;
		rb = rai;
	} else {
		qa = qai;
		ra = rai;
		qb = qbi;
		rb = rbi;
	}

	if ( ((ra+rb) > dist) && (dist > (rb - ra)) ) {
		Real dfout;
		Real fout;
		if ( ra == rb ) {
			fout = 0.5 * (1.0+0.5*dist/ra);
			dfout = 0.25/ra;
		} else {
			fout = 0.5*(1.0 + 0.5*(ra*ra - rb*rb)/(ra*dist) + 0.5*dist/ra);
			dfout = 0.25*(rb*rb - ra*ra)/(ra*dist*dist) + 0.25/ra;
		}
		Real dfin = -1.0*dfout;
		return qa*qb*(dfin/rb -
			2.0*fout/( (dist+ra+rb)*(dist+ra+rb) ) +
			2.0*dfout/( dist + ra + rb ) );
	} else {
		///jjh already fully buried
		return 0.0;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
GenBornPotential::eval_atom_derivative(
	id::AtomID const & id,
	Real const weight,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	bool const exclude_DNA_DNA,
	Vector & F1,
	Vector & F2
) const {
	using namespace etable::count_pair;

	Size const i( id.rsd() ), ii( id.atomno() );

	GenBornPoseInfo const & gb_info( static_cast< GenBornPoseInfo const & >( pose.data().get( core::pose::datacache::CacheableDataType::GEN_BORN_POSE_INFO)));

	Residue const & rsd1( pose.residue( i ) );
	GenBornResidueInfo const & gb1( gb_info.residue_info( i ) );


	// pose stuff
	Size const nres( pose.total_residue() );
	int const i_map( domain_map( i ) );
	bool const i_fixed( i_map != 0 );

	Vector const & xyzi( rsd1.xyz( ii ) );

	// gb stuff
	Real const tau = (1.0/Ep - 1.0/Ew);
	Real const q1 = rsd1.atomic_charge( ii );
	Real const b1 = gb1.born_radius( ii );

	for ( Size j=1; j<= nres; ++j ) {
		Residue const & rsd2( pose.residue( j ) );
		if ( i_fixed && domain_map( j ) == i_map ) continue; // DANGER
		if ( exclude_DNA_DNA && rsd1.is_DNA() && rsd2.is_DNA() ) continue;

		GenBornResidueInfo const & gb2( gb_info.residue_info( j ) );

		CountPairFunctionOP cpfxn( nullptr );
		if ( i == j ) {
			cpfxn = CountPairFactory::create_intrares_count_pair_function( rsd1, CP_CROSSOVER_3 );
		} else if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
			cpfxn = CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_3 );
		}

		for ( Size jj=1, jj_end=rsd2.natoms(); jj<= jj_end; ++jj ) {
			Real cp_weight;
			Size path_dist( 0 );
			if ( cpfxn && !cpfxn->count( ii, jj, cp_weight, path_dist ) ) continue; // less than 3 bonds away

			//bool const same_atom( j == i && jj == ii );
			debug_assert( j != i || jj != ii );

			Vector const & xyzj( rsd2.xyz(jj) );
			Vector const f2( xyzi - xyzj );

			Real const dis2( f2.length_squared() );
			Real const dis = std::sqrt( dis2 );

			///jjh First the dielectric screening term
			Real const q2 = rsd2.atomic_charge( jj );
			Real const b2 = gb2.born_radius( jj );

			Real const exparg = (-dis2)/(4.0*b1*b2);
			Real const expon = std::exp( exparg );
			Real const denom  = 1.0 / ( std::sqrt( dis2 + b1*b2*expon ) );
			Real const deriv_denom = denom * denom * denom;

			Real dE_dR
				( +166.0*tau*q1*q2*dis * ( 2.0 - 0.5*expon ) * deriv_denom );

			Real const r1( gb1.atomic_radius(ii) );
			Real const r2( gb2.atomic_radius(jj) );

			Real const dE_dR_coul( gb_shell_intxn_deriv( q1, r1, q2, r2, dis ) );

			dE_dR += 332.0*dE_dR_coul/Ep;

			//if ( same_atom ) dE_dR *= 0.5;

			if ( dis > 0.0 ) {
				Real const factor( weight * dE_dR / dis );
				Vector const f1( xyzi.cross( xyzj ) );
				F1 += factor * f1;
				F2 += factor * f2;
			}
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


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::GenBornResidueInfo::GenBornResidueInfo() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::GenBornResidueInfo::save( Archive & arc ) const {
	arc( CEREAL_NVP( atomic_radius_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( born_radius_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( scale_factor_ ) ); // utility::vector1<Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::GenBornResidueInfo::load( Archive & arc ) {
	arc( atomic_radius_ ); // utility::vector1<Real>
	arc( born_radius_ ); // utility::vector1<Real>
	arc( scale_factor_ ); // utility::vector1<Real>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::GenBornResidueInfo );
CEREAL_REGISTER_TYPE( core::scoring::GenBornResidueInfo )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::GenBornPoseInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( residue_info_ ) ); // utility::vector1<GenBornResidueInfoOP>
	arc( CEREAL_NVP( placeholder_residue_ ) ); // utility::vector1<ResidueOP>
	arc( CEREAL_NVP( placeholder_info_ ) ); // utility::vector1<GenBornResidueInfoOP>
	arc( CEREAL_NVP( being_packed_ ) ); // utility::vector1<_Bool>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::GenBornPoseInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( residue_info_ ); // utility::vector1<GenBornResidueInfoOP>
	arc( placeholder_residue_ ); // utility::vector1<ResidueOP>
	arc( placeholder_info_ ); // utility::vector1<GenBornResidueInfoOP>
	arc( being_packed_ ); // utility::vector1<_Bool>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::GenBornPoseInfo );
CEREAL_REGISTER_TYPE( core::scoring::GenBornPoseInfo )


/// @brief Default constructor required by cereal to deserialize this class
core::scoring::GenBornRotamerSetInfo::GenBornRotamerSetInfo() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::GenBornRotamerSetInfo::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( residue_info_ ) ); // utility::vector1<GenBornResidueInfoOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::GenBornRotamerSetInfo::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( residue_info_ ); // utility::vector1<GenBornResidueInfoOP>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::GenBornRotamerSetInfo );
CEREAL_REGISTER_TYPE( core::scoring::GenBornRotamerSetInfo )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_GenBornPotential )
#endif // SERIALIZATION

