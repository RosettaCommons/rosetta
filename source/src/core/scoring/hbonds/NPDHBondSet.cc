// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/NPDHBondSet.hh
/// @brief  Hydrogen bond set class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/hbonds/NPDHBondSet.hh>


// Package headers
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <map>

#ifdef    SERIALIZATION
// Project serialization headers
#include <core/id/AtomID_Map.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static THREAD_LOCAL basic::Tracer TR( "core.scoring.hbonds.NPDHBondSet" );


namespace core {
namespace scoring {
namespace hbonds {

using namespace ObjexxFCL::format;


NPDHBondSet::NPDHBondSet():
	parent()
{}

NPDHBondSet::~NPDHBondSet() {}

NPDHBondSet::NPDHBondSet( Size const nres ):
	parent( nres )
{
}


// Note that the following constructors are duplicated to prevent
// automatic casting.  Since owning-ptr is castible to Size and the
// above constructor takes a Size, automatic casting would
// ambiguous.
NPDHBondSet::NPDHBondSet(
	HBondOptions const & opts
) :
	parent( opts )
{}

NPDHBondSet::NPDHBondSet( HBondOptions const & opts, Size const nres ) :
	parent( opts, nres )
{
}

/// @brief convenience constructor. If you need more controlled
///construction please use one of the other constructors
///
/// The pose must be non-const because the neighbor graph may need to
/// be initialized.
NPDHBondSet::NPDHBondSet(
	pose::Pose & pose,
	EnergyMap const & weights
) :
	parent( pose, false /*bb_only*/ ),
	sfxn_weights_( weights )
{
	derive_per_hbond_donor_and_acceptor_weights( pose );
}

/// @brief convenience constructor. If you need more controlled
///construction please use one of the other constructors
///
/// The pose must be non-const because the neighbor graph may need to
/// be initialized.
NPDHBondSet::NPDHBondSet(
	HBondOptions const & opts,
	pose::Pose & pose,
	EnergyMap const & weights
) :
	parent( opts, pose, false /*bb_only*/ ),
	sfxn_weights_( weights )
{
	// verbose = true;
	derive_per_hbond_donor_and_acceptor_weights( pose );

	// verbose = false;
	//TR << "Beginning numeric deriv check" << std::endl;
	//
	//utility::vector1< Size > hbs_to_alter( 4, 0 );
	//hbs_to_alter[ 1 ] = 2;
	//hbs_to_alter[ 2 ] = 3;
	//hbs_to_alter[ 3 ] = 5;
	//hbs_to_alter[ 4 ] = 7;
	//
	//Real e_step = 1e-4;
	//for ( Size ii = 1; ii <= 4; ++ii ) {
	//	HBond const & ii_hb_const = hbond( hbs_to_alter[ ii ] );
	//	HBond & ii_hb( const_cast< HBond & > ( ii_hb_const ) );
	//	Real ii_orig_energy = ii_hb.energy();
	//	Real ii_acc_wt = ii_hb.acc_npd_weight();
	//	Real ii_don_wt = ii_hb.don_npd_weight();
	//	utility::vector1< Real > energies( 2 );
	//	for ( Size jj = 1; jj <= 2; ++jj ) {
	//		Real deltaE = ( jj == 1 ? -1 : 1 ) * e_step;
	//		Real jj_energy = ii_orig_energy + deltaE;
	//		ii_hb.energy_ = jj_energy;
	//		derive_per_hbond_donor_and_acceptor_weights( pose );
	//		for ( Size kk = 1; kk <= nhbonds(); ++kk ) {
	//			HBond const & kk_hb( hbond( kk ));
	//			energies[ jj ] += kk_hb.energy() * kk_hb.don_npd_weight() * kk_hb.acc_npd_weight();
	//		}
	//	}
	//	TR << "HBond dTotE / dE " << hbs_to_alter[ ii ] << ":    " <<
	//		( energies[ 2 ] - energies[ 1 ] ) / ( 2 * e_step ) - ii_acc_wt * ii_don_wt  << std::endl;
	//}

}


/// copy ctor
NPDHBondSet::NPDHBondSet( NPDHBondSet const & src ):
	parent( src ),
	sfxn_weights_( src.sfxn_weights_ ),
	atom_hbond_energies_( src.atom_hbond_energies_ ),
	atom_hbond_weights_(  src.atom_hbond_weights_ ),
	dwt_dE_for_hbond_(    src.dwt_dE_for_hbond_   )
{
}

/// @brief NOTE: this is going to be wrong...
NPDHBondSet::NPDHBondSet(
	NPDHBondSet const & src,
	utility::vector1< core::Size > const & exclude_list
) :
	parent( src, exclude_list ),
	sfxn_weights_( src.sfxn_weights_ ),
	atom_hbond_energies_( src.atom_hbond_energies_ ),
	atom_hbond_weights_(  src.atom_hbond_weights_ ),
	dwt_dE_for_hbond_(    src.dwt_dE_for_hbond_   )
{
}

NPDHBondSet::NPDHBondSet(
	NPDHBondSet const & src,
	utility::vector1< bool > const & residue_mask
) :
	parent( src, residue_mask ),
	sfxn_weights_( src.sfxn_weights_ ),
	atom_hbond_energies_( src.atom_hbond_energies_ ),
	atom_hbond_weights_(  src.atom_hbond_weights_ ),
	dwt_dE_for_hbond_(    src.dwt_dE_for_hbond_   )
{
}


/// copy ctor
NPDHBondSet::NPDHBondSet(
	NPDHBondSet const & src,
	Size seqpos
) :
	parent( src, seqpos ),
	sfxn_weights_( src.sfxn_weights_ ),
	atom_hbond_energies_( src.atom_hbond_energies_ ),
	atom_hbond_weights_(  src.atom_hbond_weights_ ),
	dwt_dE_for_hbond_(    src.dwt_dE_for_hbond_   )
{}


/// @brief clone this object
basic::datacache::CacheableDataOP
NPDHBondSet::clone() const
{
	return basic::datacache::CacheableDataOP( new NPDHBondSet( *this ) );
}

void
NPDHBondSet::setup_for_residue_pair_energies(
	pose::Pose const & pose,
	bool const calculate_derivative /*= false*/
) {
	parent::setup_for_residue_pair_energies( pose, calculate_derivative, false /*bb/bb only*/ );
	derive_per_hbond_donor_and_acceptor_weights( pose );

}

std::ostream &
operator<< ( std::ostream & out, const NPDHBondSet & hbond_set ){
	hbond_set.show( out );
	return out;
}

void
NPDHBondSet::show(
	std::ostream & out/*=std::cout*/
) const {
	for ( core::Size i=1; i<=nhbonds(); ++i ) {
		out << hbond(i);
	}
}


bool
operator==(NPDHBondSet const & a, NPDHBondSet const & b)
{
	if ( operator==( static_cast< HBondSet const & > (a), static_cast< HBondSet const & > (b) ) ) {
		return true;
	} else {
		return false;
	}
}

/// @detail Optionally print a header, and then a row for each hydrogen bond in the set using the iterprable version of the hbond show format formatted for easy parsing by R, Excel, etc.
void
NPDHBondSet::show(
	pose::Pose const & pose,
	bool print_header/*=true*/,
	std::ostream & out/*=std::cout*/
) const {

	for ( Size i=1; i<= nhbonds(); ++i ) {
		hbond(i).show(pose, i==1 && print_header, out);
	}

}


/// @detail Optionally print a header, and then a row for each hydrogen
///bond in the set using the iterprable version of the hbond show
///format formatted for easy parsing by R, Excel, etc.
void
NPDHBondSet::show(
	pose::Pose const & pose,
	Size const residue_number,
	bool const print_header /*=true*/,
	std::ostream & out /*=std::cout*/
) const {

	bool first_found(true);
	for ( Size i=1; i<=nhbonds(); ++i ) {
		if ( (hbond(i).don_res() == residue_number)
				|| (hbond(i).acc_res() == residue_number) ) {

			hbond(i).show(pose, print_header && first_found, out);
			first_found = false;
		}
	}
}

Real NPDHBondSet::hbond_weight( id::AtomID const & id, Size hbond_index ) const
{
	return atom_hbond_weights_[ id ][ hbond_index ];
}

Real NPDHBondSet::d_hbond_weight_dE( id::AtomID const & id, Size hbond_fixed_index, Size hbond_changing_index ) const
{
	return dwt_dE_for_hbond_[ id ][ hbond_fixed_index ][ hbond_changing_index ];
}

Real NPDHBondSet::dEtot_dEhb( Size hbond_index ) const
{
	return dEtot_dEhb_[ hbond_index ];
}

void
NPDHBondSet::derive_per_hbond_donor_and_acceptor_weights( pose::Pose const & pose )
{

	//atom_hbond_weights_.default_value( 0.0 );
	pose::initialize_atomid_map( atom_hbond_energies_, pose );
	pose::initialize_atomid_map( atom_hbond_sfxn_wtd_energies_, pose );
	pose::initialize_atomid_map( atom_hbond_weights_, pose );
	pose::initialize_atomid_map( dwt_dE_for_hbond_, pose );

	dEtot_dEhb_.resize( nhbonds() ); std::fill( dEtot_dEhb_.begin(), dEtot_dEhb_.end(), 0.0 );
	//don_dEtot_dEhb_.resize( nhbonds() ); std::fill( don_dEtot_dEhb_.begin(), don_dEtot_dEhb_.end(), 0.0 );


	//Real invKT = 1; // kT of 1?

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		conformation::Residue const & ii_rsd( pose.residue(ii));
		if ( ii_rsd.is_virtual_residue() ) continue;

		// process the acceptors first
		for ( Size jj : ii_rsd.accpt_pos() ) {

			id::AtomID jj_id( jj, ii );
			auto const & jj_hbonds = atom_hbonds_all( jj_id );
			if ( jj_hbonds.size() == 0 ) continue;

			HBAccChemType jj_type = get_hb_acc_chem_type( jj, pose.residue( ii ) );

			// Almost every acceptor type can accept two hydrogen bonds, so look for the small
			// set that only want a single hbond.
			if ( jj_type == hbacc_IMD || jj_type == hbacc_IME || jj_type == hbacc_GENERIC_RINGBB || jj_type == hbacc_GENERIC_RINGSC ) {
				get_weights_for_one_partner_hbonder( jj_id );
			} else {
				get_weights_for_two_partner_hbonder( jj_id );
			}
			for ( Size kk = 1; kk <= atom_hbonds_all( jj_id ).size(); ++kk ) {
				auto const & kk_hbond = atom_hbonds_all( jj_id )[ kk ];
				if ( kk_hbond->energy() > 0 ) continue;
				//TR << "Setting acceptor weight for hbond " << kk_hbond->index() << " with E: " << kk_hbond->energy() << " to " << atom_hbond_weights_( jj_id )[ kk ] << std::endl;
				hbond_acc_npd_weight( kk_hbond->index(), atom_hbond_weights_( jj_id )[ kk ] );
			}
		}

		// process the donors next
		for ( Size jj : ii_rsd.Hpos_polar() ) {
			id::AtomID jj_id( jj, ii );
			get_weights_for_one_partner_hbonder( jj_id );
			for ( Size kk = 1; kk <= atom_hbonds_all( jj_id ).size(); ++kk ) {
				auto const & kk_hbond = atom_hbonds_all( jj_id )[ kk ];
				if ( kk_hbond->energy() > 0 ) continue;
				hbond_don_npd_weight( kk_hbond->index(), atom_hbond_weights_( jj_id )[ kk ] );
				//TR << "Setting donor weight for hbond " << kk_hbond->index() << " with E: " << kk_hbond->energy() << " to " << atom_hbond_weights_( jj_id )[ kk ] << std::endl;
			}
		}

	}


	// Finally, iterate across all hbonds, and for each, calculate for their donors and acceptors
	// the change in the total energy of the Pose wrt the change in the energy of that hydrogen bond
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		conformation::Residue const & ii_rsd( pose.residue(ii));
		if ( ii_rsd.is_virtual_residue() ) continue;
		// process the acceptors first
		for ( Size jj : ii_rsd.accpt_pos() ) {

			id::AtomID jj_id( jj, ii );
			auto const & jj_hbonds = atom_hbonds_all( jj_id );
			if ( jj_hbonds.size() == 0 ) continue;

			sum_dEtot_dEhb_for_atom( jj_id );
		}
		// then the donors
		for ( Size jj : ii_rsd.Hpos_polar() ) {

			id::AtomID jj_id( jj, ii );
			auto const & jj_hbonds = atom_hbonds_all( jj_id );
			if ( jj_hbonds.size() == 0 ) continue;

			sum_dEtot_dEhb_for_atom( jj_id );
		}
	}
}

void
NPDHBondSet::get_weights_for_one_partner_hbonder(
	id::AtomID jj_id
)
{
	auto const & jj_hbonds = atom_hbonds_all( jj_id );
	if ( jj_hbonds.size() == 0 ) return;

	auto & jj_hbond_energies = atom_hbond_energies_[ jj_id ];
	jj_hbond_energies.resize( jj_hbonds.size(), 0.0 );
	std::transform( jj_hbonds.begin(), jj_hbonds.end(), jj_hbond_energies.begin(), []( HBondCOP const & hbond ) { return hbond->energy(); } );

	auto & jj_hbond_sfxn_wtd_energies = atom_hbond_sfxn_wtd_energies_[ jj_id ];
	jj_hbond_sfxn_wtd_energies.resize( jj_hbonds.size(), 0.0 );
	std::transform( jj_hbonds.begin(), jj_hbonds.end(), jj_hbond_sfxn_wtd_energies.begin(), [this]( HBondCOP const & hb )
		{ return npd_hb_eval_type_weight( hb->eval_type(), sfxn_weights_, hb->don_res() == hb->acc_res() ); } );
	for ( Size kk = 1; kk <= jj_hbond_energies.size(); ++kk ) {
		jj_hbond_sfxn_wtd_energies[ kk ] *= jj_hbond_energies[ kk ];
	}

	auto & jj_hbond_weights = atom_hbond_weights_[ jj_id ];
	jj_hbond_weights.resize( jj_hbonds.size(), 0.0 );

	auto & jj_hbond_dwt_dE = dwt_dE_for_hbond_[ jj_id ];
	jj_hbond_dwt_dE.resize( jj_hbonds.size() );
	for ( Size kk = 1; kk <= jj_hbonds.size(); ++kk ) {
		jj_hbond_dwt_dE[ kk ].resize( jj_hbonds.size(), 0.0 );
	}

	hbonds::get_weights_for_one_partner_hbonder( jj_hbond_sfxn_wtd_energies, jj_hbond_weights, jj_hbond_dwt_dE );

}

void
NPDHBondSet::get_weights_for_two_partner_hbonder(
	id::AtomID jj_id
)
{
	auto const & jj_hbonds = atom_hbonds_all( jj_id );
	if ( jj_hbonds.size() == 0 ) return;

	auto & jj_hbond_energies = atom_hbond_energies_[ jj_id ];
	jj_hbond_energies.resize( jj_hbonds.size(), 0.0 );
	std::transform( jj_hbonds.begin(), jj_hbonds.end(), jj_hbond_energies.begin(), []( HBondCOP const & hbond ) { return hbond->energy(); } );

	auto & jj_hbond_sfxn_wtd_energies = atom_hbond_sfxn_wtd_energies_[ jj_id ];
	jj_hbond_sfxn_wtd_energies.resize( jj_hbonds.size(), 0.0 );
	std::transform( jj_hbonds.begin(), jj_hbonds.end(), jj_hbond_sfxn_wtd_energies.begin(), [this]( HBondCOP const & hb )
		{ return npd_hb_eval_type_weight( hb->eval_type(), sfxn_weights_, hb->don_res() == hb->acc_res() ); } );
	for ( Size kk = 1; kk <= jj_hbond_energies.size(); ++kk ) {
		jj_hbond_sfxn_wtd_energies[ kk ] *= jj_hbond_energies[ kk ];
	}

	auto & jj_hbond_weights = atom_hbond_weights_[ jj_id ];
	jj_hbond_weights.resize( jj_hbonds.size(), 0.0 );

	auto & jj_hbond_dwt_dE = dwt_dE_for_hbond_[ jj_id ];
	jj_hbond_dwt_dE.resize( jj_hbonds.size() );
	for ( Size kk = 1; kk <= jj_hbonds.size(); ++kk ) {
		jj_hbond_dwt_dE[ kk ].resize( jj_hbonds.size(), 0.0 );
	}

	hbonds::get_weights_for_two_partner_hbonder( jj_hbond_sfxn_wtd_energies, jj_hbond_weights, jj_hbond_dwt_dE );
}

void
NPDHBondSet::sum_dEtot_dEhb_for_atom( id::AtomID const & jj_id )
{
	auto const & jj_hbonds = atom_hbonds_all( jj_id );
	if ( jj_hbonds.size() == 0 ) return;

	Size const nhbonds = jj_hbonds.size();
	bool is_don = jj_hbonds[ 1 ]->don_res() == jj_id.rsd();
	//utility::vector1< Real > & dEtot_dEhb( is_don ? don_dEtot_dEhb_ : acc_dEtot_dEhb_ );

	auto const & jj_hbond_dwt_dE = dwt_dE_for_hbond_[ jj_id ];

	for ( Size kk = 1; kk <= nhbonds; ++kk ) {
		HBond const & kkhb = *jj_hbonds[ kk ];
		Real kkE = kkhb.energy();
		if ( kkE > 0 ) continue;
		for ( Size ll = 1; ll <= nhbonds; ++ll ) {
			HBond const & llhb = *jj_hbonds[ ll ];
			Real llE = llhb.energy();
			if ( llE > 0 ) continue;

			Real other_wt( is_don ? kkhb.acc_npd_weight() : kkhb.don_npd_weight() );
			dEtot_dEhb_[ llhb.index() ] += jj_hbond_dwt_dE[ kk ][ ll ] * other_wt * kkE;
		}
	}
}


/// @brief depending on the hybridization type of the atom, compute the weights
/// for its hydrogen bonds. No derivative evaluation.
void
weights_for_hbonds(
	conformation::Residue const & res,
	Size atom,
	utility::vector1< Real > const & energies,
	utility::vector1< Real > & weights
)
{
	debug_assert( energies.size() == weights.size() );
	if ( energies.size() == 1 ) { weights[ 1 ] = 1; return; }

	if ( res.atom_is_polar_hydrogen( atom ) ) {
		get_weights_for_one_partner_hbonder( energies, weights );
	} else {

		HBAccChemType jj_type = get_hb_acc_chem_type( atom, res );

		// Almost every acceptor type can accept two hydrogen bonds, so look for the small
		// set that only want a single hbond.
		if ( jj_type == hbacc_IMD || jj_type == hbacc_IME || jj_type == hbacc_GENERIC_RINGBB || jj_type == hbacc_GENERIC_RINGSC ) {
			get_weights_for_one_partner_hbonder( energies, weights );
		} else {
			get_weights_for_two_partner_hbonder( energies, weights );
		}
	}
}

/// @brief depending on the hybridization type of the atom, compute the weights
/// for its hydrogen bonds. No derivative evaluation.
void
weights_and_derivs_for_hbonds(
	conformation::Residue const & res,
	Size atom,
	utility::vector1< Real > const & energies,
	utility::vector1< Real > & weights,
	utility::vector1< utility::vector1< Real > > & dwt_dE
)
{
	debug_assert( energies.size() == weights.size() );
	if ( energies.size() == 1 ) { weights[ 1 ] = 1; return; }

	if ( res.atom_is_polar_hydrogen( atom ) ) {
		get_weights_for_one_partner_hbonder( energies, weights, dwt_dE );
	} else {

		HBAccChemType jj_type = get_hb_acc_chem_type( atom, res );

		// Almost every acceptor type can accept two hydrogen bonds, so look for the small
		// set that only want a single hbond.
		if ( jj_type == hbacc_IMD || jj_type == hbacc_IME || jj_type == hbacc_GENERIC_RINGBB || jj_type == hbacc_GENERIC_RINGSC ) {
			get_weights_for_one_partner_hbonder( energies, weights, dwt_dE );
		} else {
			get_weights_for_two_partner_hbonder( energies, weights, dwt_dE );
		}
	}
}

void
get_weights_for_one_partner_hbonder(
 	utility::vector1< Real > const & energies,
	utility::vector1< Real > & weights,
	utility::vector1< utility::vector1< Real > > & dwt_dE
)
{
	debug_assert( energies.size() == weights.size() );
	debug_assert( dwt_dE.size() == energies.size() );

	Size const nhbonds = energies.size();
	if ( nhbonds == 1 ) {
		weights[ 1 ] = 1;
		dwt_dE[1][1] = 0;
		return;
	}


	Real energy_sum = 0;
	for ( Size ii = 1; ii <= nhbonds; ++ii ) {
		if ( energies[ ii ] > 0 ) continue;
		energy_sum += energies[ ii ];
	}
	if ( energy_sum == 0 ) {
		// no hydrogen bonds; don't compute 1/0
		for ( Size ii = 1; ii <= nhbonds; ++ii ) {
			weights[ ii ] = 0;
			std::fill( dwt_dE[ ii ].begin(), dwt_dE[ ii ].end(), 0.0 );
		}
	} else {

		// energy_sum better not be zero!
		Real inv_energy_sum = 1 / energy_sum;

		for ( Size ii = 1; ii <= nhbonds; ++ii ) {
			if ( energies[ ii ] > 0 ) continue;
			weights[ ii ] = energies[ ii ] * inv_energy_sum;
		}
		for ( Size ii = 1; ii <= nhbonds; ++ii ) {
			debug_assert( dwt_dE[ ii ].size() == nhbonds );
			if ( energies[ ii ] > 0 ) continue;
			for ( Size jj = 1; jj <= nhbonds; ++jj ) {
				if ( energies[ jj ] > 0 ) continue;
				if ( ii == jj ) {
					dwt_dE[ ii ][ jj ] = inv_energy_sum - weights[ ii ] * inv_energy_sum;
				} else {
					dwt_dE[ ii ][ jj ] = -1 * weights[ ii ] * inv_energy_sum;
				}
			}
		}
	}
}

Real
get_weights_for_one_partner_hbonder(
 	utility::vector1< Real > const & energies, // the sfxn-weighted energies
	utility::vector1< Real > & weights
)
{
	debug_assert( energies.size() == weights.size() );

	Size const nhbonds = energies.size();
	if ( nhbonds == 1 ) {
		weights[ 1 ] = 1;
		return energies[ 1 ];
	}

	Real energy_sum = 0;
	for ( Size ii = 1; ii <= nhbonds; ++ii ) {
		if ( energies[ ii ] > 0 ) continue;
		energy_sum += energies[ ii ];
	}
	if ( energy_sum == 0 ) {
		// no hydrogen bonds; don't compute 1/0
		for ( Size ii = 1; ii <= nhbonds; ++ii ) {
			weights[ ii ] = 0;
		}
	} else {

		// energy_sum better not be zero!
		Real inv_energy_sum = 1 / energy_sum;

		for ( Size ii = 1; ii <= nhbonds; ++ii ) {
			weights[ ii ] = energies[ ii ] < 0 ? energies[ ii ] * inv_energy_sum : 0;
		}
	}
	return energy_sum;
}


void
get_weights_for_two_partner_hbonder(
 	utility::vector1< Real > const & energies, // the sfxn-weighted energies
	utility::vector1< Real > & weights,
	utility::vector1< utility::vector1< Real > > & dwt_dE
)
{
	debug_assert( energies.size() == weights.size() );
	debug_assert( dwt_dE.size() == energies.size() );

	Size const nhbonds = energies.size();
	if ( nhbonds == 1 ) {
		weights[ 1 ] = 1;
		dwt_dE[1][1] = 0;
	}

	Real energy_sum = get_weights_for_one_partner_hbonder( energies, weights );
	if ( energy_sum < 0 ) {

		// energy_sum better not be zero!
		Real inv_energy_sum = 1 / energy_sum;

		for ( Size ii = 1; ii <= nhbonds; ++ii ) {
			debug_assert( dwt_dE[ ii ].size() == nhbonds );
			if ( energies[ ii ] > 0 ) continue;
			for ( Size jj = 1; jj <= nhbonds; ++jj ) {
				if ( energies[ jj ] > 0 ) continue;
				if ( ii == jj ) {
					dwt_dE[ ii ][ jj ] = inv_energy_sum - weights[ ii ] * inv_energy_sum;
				} else {
					dwt_dE[ ii ][ jj ] = -1 * weights[ ii ] * inv_energy_sum;
				}
			}
		}

		for ( Size ii = 1; ii <= nhbonds; ++ii ) {
			weights[ ii ] = std::sqrt( weights[ ii ] );
		}
		for ( Size ii = 1; ii <= nhbonds; ++ii ) {
			if ( weights[ ii ] == 0 ) continue;
			Real const half_ii_inv_weight = 0.5 / weights[ ii ];
			for ( Size jj = 1; jj <= nhbonds; ++jj ) {
				dwt_dE[ ii ][ jj ] *= half_ii_inv_weight;
			}
		}
	}
}

void
get_weights_for_two_partner_hbonder(
 	utility::vector1< Real > const & energies,
	utility::vector1< Real > & weights
)
{
	get_weights_for_one_partner_hbonder( energies, weights );
	for ( Size ii = 1; ii <= energies.size(); ++ii ) {
		weights[ ii ] = std::sqrt( weights[ ii ] );
	}
}


}
}
}


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::hbonds::NPDHBondSet::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::hbonds::HBondSet >( this ) );
	arc( CEREAL_NVP( sfxn_weights_ ) ); // EnergyMap
	arc( CEREAL_NVP( atom_hbond_energies_ ) ); // id::AtomID_Map<utility::vector1<Real> >
	arc( CEREAL_NVP( atom_hbond_sfxn_wtd_energies_ ) ); // id::AtomID_Map<utility::vector1<Real> >
	arc( CEREAL_NVP( atom_hbond_weights_ ) ); // id::AtomID_Map<utility::vector1<Real> >
	arc( CEREAL_NVP( dwt_dE_for_hbond_ ) ); // id::AtomID_Map<utility::vector1<utility::vector1<Real> > >
	arc( CEREAL_NVP( dEtot_dEhb_ ) ); // utility::vector1<Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::hbonds::NPDHBondSet::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::hbonds::HBondSet >( this ) );
	arc( sfxn_weights_ ); // EnergyMap
	arc( atom_hbond_energies_ ); // id::AtomID_Map<utility::vector1<Real> >
	arc( atom_hbond_sfxn_wtd_energies_ ); // id::AtomID_Map<utility::vector1<Real> >
	arc( atom_hbond_weights_ ); // id::AtomID_Map<utility::vector1<Real> >
	arc( dwt_dE_for_hbond_ ); // id::AtomID_Map<utility::vector1<utility::vector1<Real> > >
	arc( dEtot_dEhb_ ); // utility::vector1<Real>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::hbonds::NPDHBondSet );
CEREAL_REGISTER_TYPE( core::scoring::hbonds::NPDHBondSet )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_hbonds_NPDHBondSet )
#endif // SERIALIZATION
