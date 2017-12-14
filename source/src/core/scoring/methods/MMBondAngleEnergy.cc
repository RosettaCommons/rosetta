// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MMBondAngleAnergy.cc
/// @brief  molecular mechanics torsion energy
/// @author Colin A. Smith (colin.smith@ucsf.edu)

// Unit headers
#include <core/scoring/methods/MMBondAngleEnergy.hh>
#include <core/scoring/methods/MMBondAngleEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>  // for is_virtual()
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParam.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/string_util.hh>

// Numeric headers
#include <numeric/xyz.functions.hh>
#include <numeric/deriv/angle_deriv.hh>

// C++ headers
#include <iostream>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the MMBondAngleEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
MMBondAngleEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new MMBondAngleEnergy( options ) );
}

ScoreTypes
MMBondAngleEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( mm_bend );
	return sts;
}


static basic::Tracer TR( "core.mm.MMBondAngleEnergy" );

MMBondAngleEnergy::MMBondAngleEnergy( methods::EnergyMethodOptions const & options ):
	parent( methods::EnergyMethodCreatorOP( new MMBondAngleEnergyCreator ) ),
	param_set_( /* NULL */ ),
	central_atoms_to_score_( options.bond_angle_central_atoms_to_score() )
{
	if ( options.bond_angle_residue_type_param_set() ) {
		param_set_ = core::scoring::mm::MMBondAngleResidueTypeParamSetOP( new core::scoring::mm::MMBondAngleResidueTypeParamSet( *(options.bond_angle_residue_type_param_set()) ) );
	}
}

MMBondAngleEnergy::MMBondAngleEnergy( MMBondAngleEnergy const & src ):
	parent( src ),
	param_set_( /* NULL */ ),
	central_atoms_to_score_( src.central_atoms_to_score_ )
{
	if ( src.param_set_ ) {
		param_set_ = core::scoring::mm::MMBondAngleResidueTypeParamSetOP( new core::scoring::mm::MMBondAngleResidueTypeParamSet( *src.param_set_ ) );
	}
}

MMBondAngleEnergy::~MMBondAngleEnergy() = default;


/// clone
EnergyMethodOP
MMBondAngleEnergy::clone() const
{
	return EnergyMethodOP( new MMBondAngleEnergy( *this ) );
}


void
MMBondAngleEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const & ) const
{
	pose.update_residue_neighbors();
	pose.update_actcoords();
}


void
MMBondAngleEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	pose.update_actcoords();

	// make sure the parameters for all residues are loaded
	if ( param_set_ ) {
		for ( Size i = 1; i <= pose.size(); ++i ) {
			param_set_->get(pose.residue_type(i));
		}
	}
}


void
MMBondAngleEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
	pose.update_actcoords();
}


bool
MMBondAngleEnergy::defines_intrares_energy( EnergyMap const & ) const
{
	return true;
}

void
MMBondAngleEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & ,
	EnergyMap & emap
) const
{
	// bail out if the residues aren't bonded
	if ( !rsd1.is_bonded(rsd2) ) return;

	Real energy = 0;

	// get residue types
	chemical::ResidueType const & rsd1_type = rsd1.type();
	chemical::ResidueType const & rsd2_type = rsd2.type();

	TR(basic::t_trace) << "residue_pair_energy: processing residues "
		<< rsd1.seqpos() << "." << rsd1_type.name() << "-"
		<< rsd2.seqpos() << "." << rsd2_type.name() << std::endl;

	utility::vector1< Size > const & r1_resconn_ids( rsd1.connections_to_residue( rsd2 ) );

	for ( Size ii = 1; ii <= r1_resconn_ids.size(); ++ii ) {

		Size const resconn_id1( r1_resconn_ids[ii] );
		Size const resconn_id2( rsd1.residue_connection_conn_id( resconn_id1 ) );

		Size const resconn_atomno1( rsd1.residue_connection( resconn_id1 ).atomno() );
		Size const resconn_atomno2( rsd2.residue_connection( resconn_id2 ).atomno() );

		TR(basic::t_trace) << "Found residue connection id " << resconn_id1 << "-" << resconn_id2 << ": "
			<< rsd1.atom_name( resconn_atomno1 ) << "-" << rsd2.atom_name( resconn_atomno2 ) << std::endl;

		Size const resconn_mmat1 = rsd1_type.atom( resconn_atomno1 ).mm_atom_type_index();
		Size const resconn_mmat2 = rsd2_type.atom( resconn_atomno2 ).mm_atom_type_index();

		/// Score only a subset of the bong angles using the central_atoms_to_score_ list
		/// of atom names.  If this list is empty, then score all atom pairs.
		/// determine if residue connection atoms should be scored
		bool should_score1( true );
		bool should_score2( true );

		if ( !param_set_ && central_atoms_to_score_.size() != 0 ) {
			should_score1 = score_atom_centrally( rsd1_type, resconn_atomno1 );
			should_score2 = score_atom_centrally( rsd2_type, resconn_atomno2 );
		}

		if ( should_score1 ) {
			/// compute the bond-angle energies from pairs of atoms within-1 bond on rsd1 with
			/// the the connection atom on rsd2.
			utility::vector1< chemical::two_atom_set > const & rsd1_atoms_wi1_bond_of_ii(
				rsd1_type.atoms_within_one_bond_of_a_residue_connection( resconn_id1 ));
			for ( Size jj = 1; jj <= rsd1_atoms_wi1_bond_of_ii.size(); ++jj ) {
				debug_assert( rsd1_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno1 );
				Size const res1_lower_atomno = rsd1_atoms_wi1_bond_of_ii[ jj ].key2();
				Size const res1_lower_mmtype = rsd1_type.atom( res1_lower_atomno ).mm_atom_type_index();

				Real const angle = numeric::angle_radians(
					rsd1.atom( res1_lower_atomno ).xyz(),
					rsd1.atom( resconn_atomno1 ).xyz(),
					rsd2.atom( resconn_atomno2 ).xyz() );

				TR(basic::t_trace)
					<< "r1 " << res1_lower_atomno << " " << rsd1.atom_name( res1_lower_atomno ) << "("
					<< rsd1_type.atom( res1_lower_atomno ).mm_name() << ") - "
					<< "r1 " << resconn_atomno1 << " " << rsd1.atom_name( resconn_atomno1 ) << "("
					<< rsd1_type.atom( resconn_atomno1 ).mm_name() << ") - "
					<< "r2 " << resconn_atomno2 << " " << rsd2.atom_name( resconn_atomno2 ) << "("
					<< rsd2_type.atom( resconn_atomno2 ).mm_name() << ")" << std::endl;

				if ( param_set_ ) {

					// lookup Ktheta and theta0
					Real Ktheta;
					Real theta0;
					param_set_->lookup(
						pose.conformation(),
						id::AtomID(res1_lower_atomno, rsd1.seqpos()),
						id::AtomID(resconn_atomno1, rsd1.seqpos()),
						id::AtomID(resconn_atomno2, rsd2.seqpos()),
						Ktheta,
						theta0
					);

					// accumulate the energy
					energy += potential_.score( Ktheta, theta0, angle );

					continue;
				}

				energy += potential_.mm::MMBondAngleScore::score
					(  mm::mm_bondangle_atom_tri( res1_lower_mmtype, resconn_mmat1, resconn_mmat2 ), angle );

				//TR << rsd1.atom_name( res1_lower_atomno ) << "(" << rsd1_type.atom( res1_lower_atomno ).mm_name() << ") - "
				//<< rsd1.atom_name( resconn_atomno1 ) << "(" << rsd1_type.atom( resconn_atomno1 ).mm_name() << ") - "
				//<< rsd2.atom_name( resconn_atomno2 ) << "(" << rsd2_type.atom( resconn_atomno2 ).mm_name() << ")   \t"
				//<< potential_.mm::MMBondAngleScore::score
				//(  mm::mm_bondangle_atom_tri( res1_lower_mmtype, resconn_mmat1, resconn_mmat2 ), angle ) << std::endl;

			}
		}

		if ( should_score2 ) {
			/// compute the bond-angle energies from pairs of atoms within-1 bond on rsd2 with
			/// the the connection atom on rsd1.
			utility::vector1< chemical::two_atom_set > const & rsd2_atoms_wi1_bond_of_ii(
				rsd2_type.atoms_within_one_bond_of_a_residue_connection( resconn_id2 ));
			for ( Size jj = 1; jj <= rsd2_atoms_wi1_bond_of_ii.size(); ++jj ) {
				debug_assert( rsd2_atoms_wi1_bond_of_ii[ jj ].key1() == resconn_atomno2 );
				Size const res2_lower_atomno = rsd2_atoms_wi1_bond_of_ii[ jj ].key2();
				Size const res2_lower_mmtype = rsd2_type.atom( res2_lower_atomno ).mm_atom_type_index();

				Real const angle = numeric::angle_radians(
					rsd2.atom( res2_lower_atomno ).xyz(),
					rsd2.atom( resconn_atomno2 ).xyz(),
					rsd1.atom( resconn_atomno1 ).xyz() );

				TR(basic::t_trace) << "r2 " << res2_lower_atomno << " "
					<< rsd2.atom_name( res2_lower_atomno ) << "("
					<< rsd2_type.atom( res2_lower_atomno ).mm_name() << ") - " << "r2 " << resconn_atomno2
					<< " " << rsd2.atom_name( resconn_atomno2 ) << "("
					<< rsd2_type.atom( resconn_atomno2 ).mm_name() << ") - " << "r1 " << resconn_atomno1
					<< " " << rsd1.atom_name( resconn_atomno1 ) << "("
					<< rsd1_type.atom( resconn_atomno1 ).mm_name() << ")" << std::endl;

				if ( param_set_ ) {

					// lookup Ktheta and theta0
					Real Ktheta;
					Real theta0;
					param_set_->lookup(
						pose.conformation(),
						id::AtomID(res2_lower_atomno, rsd2.seqpos()),
						id::AtomID(resconn_atomno2, rsd2.seqpos()),
						id::AtomID(resconn_atomno1, rsd1.seqpos()),
						Ktheta,
						theta0
					);

					// accumulate the energy
					energy += potential_.score( Ktheta, theta0, angle );

					continue;
				}

				energy += potential_.mm::MMBondAngleScore::score
					(  mm::mm_bondangle_atom_tri( res2_lower_mmtype, resconn_mmat2, resconn_mmat1 ), angle );
				//TR << rsd2.atom_name( res2_lower_atomno ) << "(" << rsd2_type.atom( res2_lower_atomno ).mm_name() << ") - "
				//<< rsd2.atom_name( resconn_atomno2 ) << "(" << rsd2_type.atom( resconn_atomno2 ).mm_name() << ") - "
				//<< rsd1.atom_name( resconn_atomno1 ) << "(" << rsd1_type.atom( resconn_atomno1 ).mm_name() << ")   \t"
				//<< potential_.mm::MMBondAngleScore::score
				//(  mm::mm_bondangle_atom_tri( res2_lower_mmtype, resconn_mmat2, resconn_mmat1 ), angle ) << std::endl;
			}
		}
	}

	emap[ mm_bend ] += energy;
}

void
MMBondAngleEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap & emap
) const
{
	Real energy = 0;

	// get residue type
	chemical::ResidueType const & rsd_type = rsd.type();

	if ( param_set_ ) {

		mm::MMBondAngleResidueTypeParam const & rsd_param(param_set_->get(rsd_type));

		for ( Size bondang = 1; bondang <= rsd_param.num_bondangles(); ++bondang ) {

			if ( rsd_param.Ktheta( bondang ) ) {

				// get ResidueType ints
				Size rt1 = ( rsd_param.bondangle( bondang ) ).key1();
				Size rt2 = ( rsd_param.bondangle( bondang ) ).key2();
				Size rt3 = ( rsd_param.bondangle( bondang ) ).key3();

				// get angle
				Real angle = numeric::angle_radians
					( rsd.atom( rt1 ).xyz(), rsd.atom( rt2 ).xyz(), rsd.atom( rt3 ).xyz() );

				// accumulate the energy
				energy += potential_.score( rsd_param.Ktheta( bondang ), rsd_param.theta0( bondang ), angle );
			}
		}

		// add energy to emap
		emap[ mm_bend ] += energy;

		return;
	}

	TR(basic::t_trace) << "MMIntrares Processing residue " << rsd.seqpos() << " " << rsd_type.name()
		<< std::endl;

	// for each dihedral angle in the residue type
	for ( Size bondang = 1; bondang <= rsd_type.num_bondangles(); ++bondang ) {
		// get ResidueType ints
		Size rt1 = ( rsd_type.bondangle( bondang ) ).key1();
		Size rt2 = ( rsd_type.bondangle( bondang ) ).key2();
		Size rt3 = ( rsd_type.bondangle( bondang ) ).key3();

		if ( central_atoms_to_score_.size() ) {

			bool should_score( score_atom_centrally( rsd_type, rt2 ));
			if ( !should_score ) continue;
		}

		// check for vrt
		if ( rsd_type.atom_type(rt1).is_virtual()
				|| rsd_type.atom_type(rt2).is_virtual()
				|| rsd_type.atom_type(rt3).is_virtual() ) {
			continue;
		}

		// get mm atom type index
		int mmat1 = rsd_type.atom( rt1 ).mm_atom_type_index();
		int mmat2 = rsd_type.atom( rt2 ).mm_atom_type_index();
		int mmat3 = rsd_type.atom( rt3 ).mm_atom_type_index();

		// get angle
		Real angle = numeric::angle_radians
			( rsd.atom( rt1 ).xyz(), rsd.atom( rt2 ).xyz(), rsd.atom( rt3 ).xyz() );

		TR(basic::t_trace) << rt1 << " " << rsd.atom_name( rt1 ) << "(" << rsd_type.atom( rt1 ).mm_name()
			<< ") - " << rt2 << " " << rsd.atom_name( rt2 ) << "(" << rsd_type.atom( rt2 ).mm_name()
			<< ") - " << rt3 << " " << rsd.atom_name( rt3 ) << "(" << rsd_type.atom( rt3 ).mm_name()
			<< ") angle " << angle << std::endl;

		// score bond angle
		energy += potential_.mm::MMBondAngleScore::score
			(  mm::mm_bondangle_atom_tri( mmat1, mmat2, mmat3 ), angle );

		// debug code
		//TR << rsd.atom_name( rt1 ) << "(" << rsd_type.atom( rt1 ).mm_name() << ") - "
		//   << rsd.atom_name( rt2 ) << "(" << rsd_type.atom( rt2 ).mm_name() << ") - "
		//   << rsd.atom_name( rt3 ) << "(" << rsd_type.atom( rt3 ).mm_name() << ")   \t"
		//   << potential_.mm::MMBondAngleScore::score
		//              (  mm::mm_bondangle_atom_tri( mmat1, mmat2, mmat3 ), angle ) << std::endl;
	}

	// add energy to emap
	emap[ mm_bend ] += energy;
}

bool
MMBondAngleEnergy::score_atom_centrally(
	core::chemical::ResidueType const & restype,
	Size atomno
) const
{
	if ( central_atoms_to_score_.size() == 0 ) return true;

	for ( Size ii = 1; ii <= central_atoms_to_score_.size(); ++ii ) {
		if ( utility::same_ignoring_spaces( restype.atom_name( atomno ), central_atoms_to_score_[ ii ] ) ) {
			//std::cout << "Counting atom " << atomno << " " << restype.name() << " " << restype.atom_name( atomno ) << std::endl;
			return true;
		}
	}
	return false;
}

void
MMBondAngleEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	//TR << "Evaluating derivatives for rsd# " << id.rsd() << " atom# " << id.atomno() << std::endl;

	Vector LF1( 0.0 );
	Vector LF2( 0.0 );

	//Vector f1(0.0), f2(0.0);
	/// 0. If we're scoring only a subset of atom centers,
	///    Ask if this atom should be scored as a center atom and
	///    then for each neighbor of this atom, ask if
	///    it should be scored.  If this atom forms an inter-residue
	///    connection or is bonded to an atom that does, figure out
	///    if the interresidue bond angles should be scored.
	/// 1. Evaluate intraresidue bond angle derivatives.
	///    for bonds where this atom is scored
	/// 2. Evaluate interresidue bond angle derivatives.
	///    KNOWN BUG: This does not evaluate derivatives for bond angles
	///    spanning three or more residues.
	///    a.  Calculate derivatives when exactly two of the atoms come from id.rsd()
	///    b.  Calculate derivatives when exactly one of the atoms comes from id.rsd()

	core::chemical::ResidueType const & restype( pose.residue_type( id.rsd() ) );
	core::conformation::Residue const & res( pose.residue( id.rsd() ) );
	Size const atomno( id.atomno());

	bool const score_everything( central_atoms_to_score_.size() == 0 || param_set_ );
	bool const score_this_atom_centrally( score_everything || score_atom_centrally( restype, atomno ));

	/// 1.
	utility::vector1< Size > const & angs( restype.bondangles_for_atom( id.atomno() ) );

	for ( Size ii = 1, ii_end = angs.size(); ii <= ii_end; ++ii ) {
		chemical::bondangle_atom_set const & ii_bangle( restype.bondangle( angs[ ii ] ) );
		debug_assert( ii_bangle.key1() == atomno || ii_bangle.key2() == atomno || ii_bangle.key3() == atomno );
		if ( score_everything ||
				( score_this_atom_centrally && ii_bangle.key2() == atomno ) ||
				( score_atom_centrally( restype, ii_bangle.key2() ) ) ) {    // could be faster here

			//TR << "  1. Deriv for angle# " << angs[ ii ] << " between " << ii_bangle.key1() << " " << ii_bangle.key2()  << " " << ii_bangle.key3() << std::endl;

			Real Ktheta;
			Real theta0;
			if ( param_set_ ) {
				param_set_->lookup(
					pose.conformation(),
					id::AtomID(ii_bangle.key1(), id.rsd()),
					id::AtomID(ii_bangle.key2(), id.rsd()),
					id::AtomID(ii_bangle.key3(), id.rsd()),
					Ktheta,
					theta0
				);
				// skip this bond angle if Ktheta is 0
				if ( !Ktheta ) continue;
			}

			// check for vrt
			if ( restype.atom_type(ii_bangle.key1()).is_virtual()
					|| restype.atom_type(ii_bangle.key2()).is_virtual()
					|| restype.atom_type(ii_bangle.key3()).is_virtual() ) {
				continue;
			}

			Vector f1(0.0), f2(0.0);
			Real theta(0.0);
			if ( ii_bangle.key1() == atomno ) {
				numeric::deriv::angle_p1_deriv(
					res.xyz( ii_bangle.key1() ),
					res.xyz( ii_bangle.key2() ),
					res.xyz( ii_bangle.key3() ),
					theta, f1, f2 );
			} else if ( ii_bangle.key2() == atomno ) {
				numeric::deriv::angle_p2_deriv(
					res.xyz( ii_bangle.key1() ),
					res.xyz( ii_bangle.key2() ),
					res.xyz( ii_bangle.key3() ),
					theta, f1, f2 );
			} else {
				numeric::deriv::angle_p1_deriv(
					res.xyz( ii_bangle.key3() ),
					res.xyz( ii_bangle.key2() ),
					res.xyz( ii_bangle.key1() ),
					theta, f1, f2 );
			}

			Real dE_dtheta;
			if ( param_set_ ) {
				dE_dtheta = weights[ mm_bend ] * potential_.dscore( Ktheta, theta0, theta );
			} else {
				Size const mmat1 = restype.atom( ii_bangle.key1() ).mm_atom_type_index();
				Size const mmat2 = restype.atom( ii_bangle.key2() ).mm_atom_type_index();
				Size const mmat3 = restype.atom( ii_bangle.key3() ).mm_atom_type_index();

				dE_dtheta = weights[ mm_bend ] * potential_.mm::MMBondAngleScore::dscore(
					mm::mm_bondangle_atom_tri( mmat1, mmat2, mmat3 ), theta );
			}
			//TR << "    theta " << numeric::conversions::degrees( theta ) << " dE_dtheta " << dE_dtheta << " f1: " << f1.x() << " " << f1.y() << " " << f1.z() << " f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;

			LF1 += dE_dtheta * f1;
			LF2 += dE_dtheta * f2;
		}
	}

	/// 2.
	/// a. Bond angles involving two atoms on this residue.
	utility::vector1< std::pair< Size, Size > > const & interres_wi1_for_this_atom(
		restype.within1bonds_sets_for_atom( atomno ));
	for ( Size ii = 1; ii <= interres_wi1_for_this_atom.size(); ++ii ) {
		Size const ii_resconn   = interres_wi1_for_this_atom[ ii ].first;
		Size const ii_whichpair = interres_wi1_for_this_atom[ ii ].second;
		chemical::two_atom_set const & ii_pair = restype.
			atoms_within_one_bond_of_a_residue_connection( ii_resconn )[ ii_whichpair ];
		debug_assert( ii_pair.key1() == atomno || ii_pair.key2() == atomno );

		/// Find the neighbor residue and atom
		Size const ii_neighb      = res.residue_connection_partner( ii_resconn );
		Size const neighb_resconn = res.residue_connection_conn_id( ii_resconn );
		conformation::Residue const & neighb_res( pose.residue( ii_neighb ));
		chemical::ResidueType const & neighb_restype( pose.residue_type( ii_neighb ) );
		Size const neighb_atom    = neighb_restype.residue_connection( neighb_resconn ).atomno();

		if ( score_everything ||
				( score_this_atom_centrally && ii_pair.key1() == atomno ) ||
				( score_atom_centrally( restype, ii_pair.key1() )) ) {

			//TR << "  2a. Deriv for interres angle between (" << ii_neighb << "," << neighb_atom << ") " << ii_pair.key1()  << " " << ii_pair.key2() << std::endl;

			Real Ktheta;
			Real theta0;
			if ( param_set_ ) {
				param_set_->lookup(
					pose.conformation(),
					id::AtomID(neighb_atom, neighb_res.seqpos()),
					id::AtomID(ii_pair.key1(), id.rsd()),
					id::AtomID(ii_pair.key2(), id.rsd()),
					Ktheta,
					theta0
				);
				// skip this bond angle if Ktheta is 0
				if ( !Ktheta ) continue;
			}

			Vector f1(0.0), f2(0.0);
			Real theta(0.0);
			if ( ii_pair.key1() == atomno ) {
				numeric::deriv::angle_p2_deriv(
					neighb_res.xyz( neighb_atom ),
					res.xyz( ii_pair.key1() ),
					res.xyz( ii_pair.key2() ),
					theta, f1, f2 );
			} else { // ( ii_pair.key2() == atomno )
				numeric::deriv::angle_p1_deriv(
					res.xyz( ii_pair.key2() ),
					res.xyz( ii_pair.key1() ),
					neighb_res.xyz( neighb_atom ),
					theta, f1, f2 );
			}

			Real dE_dtheta;
			if ( param_set_ ) {
				dE_dtheta = weights[ mm_bend ] * potential_.dscore( Ktheta, theta0, theta );
			} else {
				Size const mmat1 = neighb_restype.atom( neighb_atom ).mm_atom_type_index();
				Size const mmat2 = restype.atom( ii_pair.key1() ).mm_atom_type_index();
				Size const mmat3 = restype.atom( ii_pair.key2() ).mm_atom_type_index();

				dE_dtheta = weights[ mm_bend ] * potential_.mm::MMBondAngleScore::dscore(
					mm::mm_bondangle_atom_tri( mmat1, mmat2, mmat3 ), theta );
			}

			//TR << "    theta " << numeric::conversions::degrees( theta ) << " dE_dtheta " << dE_dtheta << " f1: " << f1.x() << " " << f1.y() << " " << f1.z() << " f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;

			LF1 += dE_dtheta * f1;
			LF2 += dE_dtheta * f2;
		}
	}

	/// b. Bond angles involving two atoms on the neighbor residues.
	/// Iterate across residue connections for this atom -- for each residue connection,
	/// find the residue connection partner, and the list of within-1-bond atom pairs
	/// for the atom forming the residue connection.  Then iterate across this list.
	utility::vector1< Size > const & connections( restype.residue_connections_for_atom( atomno ) );
	for ( Size ii = 1; ii <= connections.size(); ++ii ) {
		Size const ii_resconn   = connections[ ii ];

		/// Find the neighbor residue and connection atom
		Size const ii_neighb      = res.residue_connection_partner( ii_resconn );
		Size const neighb_resconn = res.residue_connection_conn_id( ii_resconn );
		conformation::Residue const & neighb_res( pose.residue( ii_neighb ));
		chemical::ResidueType const & neighb_restype( pose.residue_type( ii_neighb ) );
		Size const neighb_atom1    = neighb_restype.residue_connection( neighb_resconn ).atomno();

		if ( ! score_everything && ! score_atom_centrally( neighb_restype, neighb_atom1 ) ) continue;
		// Beyond this point, we will score all atom triples

		Size const mmat1 = restype.atom( atomno ).mm_atom_type_index();
		Size const mmat2 = neighb_restype.atom( neighb_atom1 ).mm_atom_type_index();


		utility::vector1< chemical::two_atom_set > const & neighb_atoms_wi1(
			neighb_restype.atoms_within_one_bond_of_a_residue_connection( neighb_resconn ));

		for ( Size jj = 1; jj <= neighb_atoms_wi1.size(); ++jj ) {

			chemical::two_atom_set const & neighb_pair = neighb_atoms_wi1[ jj ];
			debug_assert( neighb_pair.key1() == neighb_atom1 );

			Size const neighb_atom2 = neighb_pair.key2();
			Size const mmat3 = neighb_restype.atom( neighb_atom2 ).mm_atom_type_index();

			//TR << "  2b. Deriv for interres angle between (" << ii_neighb << "," << neighb_atom2 << ") ("<< ii_neighb << "," << neighb_atom1 << ") " << atomno << std::endl;

			Real Ktheta;
			Real theta0;
			if ( param_set_ ) {
				param_set_->lookup(
					pose.conformation(),
					id::AtomID(atomno, id.rsd()),
					id::AtomID(neighb_atom1, neighb_res.seqpos()),
					id::AtomID(neighb_atom2, neighb_res.seqpos()),
					Ktheta,
					theta0
				);
				// skip this bond angle if Ktheta is 0
				if ( !Ktheta ) continue;
			}

			Vector f1(0.0), f2(0.0);
			Real theta(0.0);

			numeric::deriv::angle_p1_deriv(
				res.xyz( atomno ),
				neighb_res.xyz( neighb_atom1 ),
				neighb_res.xyz( neighb_atom2 ),
				theta, f1, f2 );

			Real dE_dtheta;
			if ( param_set_ ) {
				dE_dtheta = weights[ mm_bend ] * potential_.dscore( Ktheta, theta0, theta );
			} else {
				dE_dtheta = weights[ mm_bend ] * potential_.mm::MMBondAngleScore::dscore(
					mm::mm_bondangle_atom_tri( mmat1, mmat2, mmat3 ), theta );
			}

			//TR << "    theta " << numeric::conversions::degrees( theta ) << " dE_dtheta " << dE_dtheta << " f1: " << f1.x() << " " << f1.y() << " " << f1.z() << " f2: " << f2.x() << " " << f2.y() << " " << f2.z() << std::endl;

			LF1 += dE_dtheta * f1;
			LF2 += dE_dtheta * f2;
		}
	}
	//TR << "Finished derivatives for rsd# " << id.rsd() << " atom# " << id.atomno() <<  " F1: " << LF1.x() << " " << LF1.y() << " " << LF1.z() << " F2: " << LF2.x() << " " << LF2.y() << " " << LF2.z() << std::endl;
	F1 += LF1;
	F2 += LF2;
}


/// @brief MMBondAngleEnergy does not have an atomic interation threshold
Distance
MMBondAngleEnergy::atomic_interaction_cutoff() const
{
	return 0.0;
}

/// @brief MMBondAngleEnergy is context independent; indicates that no context graphs are required
void
MMBondAngleEnergy::indicate_required_context_graphs(utility::vector1< bool > & ) const
{}
core::Size
MMBondAngleEnergy::version() const
{
	return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core
