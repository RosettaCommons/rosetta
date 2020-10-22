// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/dunbrack/PeptoidDOFReporters.hh
/// @brief   Class to measure the DOFs used by a RotamerLibrary
/// @author  Andrew Leaver-Fay

// Unit Headers
#include <core/pack/dunbrack/PeptoidDOFReporters.hh>

// Package headers
#include <core/id/PartialAtomID.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/pose/Pose.hh>

namespace core {
namespace pack {
namespace dunbrack {

PeptoidOmegaReporter::PeptoidOmegaReporter(Real const neutral_omega) :
	peptoid_neutral_omega_(neutral_omega)
{}

PeptoidOmegaReporter::~PeptoidOmegaReporter() = default;

/// @note Ugh, it's real ugly to require that the atoms that define
/// the omega dof live in the Pose
Real
PeptoidOmegaReporter::get_dof(
	conformation::Residue const & rsd,
	pose::Pose const & pose
) const
{
	if ( pose.conformation().num_chains() == 2 || rsd.chain() == pose.conformation().num_chains()-1 ) {
		// chain_end won't be the last residue, because the last residue of a conformation isn't the chain ending for some reason...
		if ( rsd.has_variant_type( chemical::NTERM_CONNECT ) && pose.residue( pose.size() ).has_variant_type( chemical::CTERM_CONNECT ) ) {
			debug_assert( pose.residue( pose.size() ).is_protein() || pose.residue( pose.size() ).is_peptoid() );
			return pose.residue( pose.size() ).mainchain_torsion( rsd.mainchain_torsions().size() );
		} else if ( rsd.is_lower_terminus() ) {
			return peptoid_neutral_omega_;
		} else {
			// Don't look for seqpos-1 but rather what's bonded at lower.
			Size const lower_seqpos = rsd.residue_connection_partner( rsd.type().lower_connect_id() );
			debug_assert( pose.residue( lower_seqpos ).is_protein() || pose.residue( lower_seqpos ).is_peptoid() );
			return pose.residue( lower_seqpos ).mainchain_torsion( rsd.mainchain_torsions().size() ); // last torsion
		}
	} else {
		if ( rsd.has_variant_type( chemical::NTERM_CONNECT ) && pose.residue( pose.conformation().chain_end( rsd.chain() ) ).has_variant_type( chemical::CTERM_CONNECT ) ) {
			debug_assert( pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_protein() ||
				pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_peptoid() );
			return pose.residue( pose.conformation().chain_end( rsd.chain() ) ).mainchain_torsion( rsd.mainchain_torsions().size() );
		} else if ( rsd.is_lower_terminus() ) {
			return peptoid_neutral_omega_;
		} else {
			Size const lower_seqpos = rsd.residue_connection_partner( rsd.type().lower_connect_id() );
			debug_assert( pose.residue( lower_seqpos ).is_protein() || pose.residue( lower_seqpos ).is_peptoid() );
			return pose.residue( lower_seqpos ).mainchain_torsion( rsd.mainchain_torsions().size() );
		}
	}
}

void
PeptoidOmegaReporter::insert_atoms_defining_dof(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	std::set< id::PartialAtomID > & atoms
) const
{
	auto insert_for_residue = [&]( Size other_resid ) {
		Size const lower_connect_ind = rsd.lower_connect().index();
		Size const conn_id_on_lower_residue = rsd.residue_connection_conn_id(lower_connect_ind);
		for ( Size ii = 0; ii <= 1; ++ii ) {
			atoms.insert( id::PartialAtomID( conn_id_on_lower_residue, other_resid, ii ));
		}
		atoms.insert( id::PartialAtomID( rsd.mainchain_atoms()[1], rsd.seqpos()) );
		atoms.insert( id::PartialAtomID( rsd.mainchain_atoms()[2], rsd.seqpos()) );
	};

	if ( pose.conformation().num_chains() == 2 || rsd.chain() == pose.conformation().num_chains()-1 ) {
		// chain_end won't be the last residue, because the last residue of a conformation isn't the chain ending for some reason...
		if ( rsd.has_variant_type( chemical::NTERM_CONNECT ) && pose.residue( pose.size() ).has_variant_type( chemical::CTERM_CONNECT ) ) {
			debug_assert( pose.residue( pose.size() ).is_protein() || pose.residue( pose.size() ).is_peptoid() );
			insert_for_residue( pose.size() );
		} else if ( rsd.is_lower_terminus() ) {
			// no op
		} else {
			insert_for_residue( rsd.residue_connection_partner( rsd.type().lower_connect_id() ));
		}
	} else {
		if ( rsd.has_variant_type( chemical::NTERM_CONNECT ) && pose.residue( pose.conformation().chain_end( rsd.chain() ) ).has_variant_type( chemical::CTERM_CONNECT ) ) {
			debug_assert( pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_protein() ||
				pose.residue( pose.conformation().chain_end( rsd.chain() ) ).is_peptoid() );
			insert_for_residue( pose.conformation().chain_end( rsd.chain() ) );
		} else if ( rsd.is_lower_terminus() ) {
			// no op
		} else {
			insert_for_residue( rsd.residue_connection_partner( rsd.type().lower_connect_id()));
		}
	}
}



PeptoidGeneralDOFReporter::PeptoidGeneralDOFReporter(Size const mainchain_torsion, Size const n_tors, Real const neutral_val):
	mainchain_torsion_( mainchain_torsion ),
	n_tors_(  n_tors ),
	neutral_val_( neutral_val )
{}

PeptoidGeneralDOFReporter::~PeptoidGeneralDOFReporter() = default;

Real
PeptoidGeneralDOFReporter::get_dof( conformation::Residue const & rsd, pose::Pose const & ) const
{
	if      ( rsd.is_lower_terminus() && mainchain_torsion_ == 2 ) return neutral_val_;
	else if ( rsd.is_upper_terminus() && mainchain_torsion_ == n_tors_ ) return neutral_val_;
	else return rsd.mainchain_torsion( mainchain_torsion_ ); // 2 == phi, 3 == psi
}

void
PeptoidGeneralDOFReporter::insert_atoms_defining_dof(
	conformation::Residue const & rsd,
	pose::Pose const &,
	std::set< id::PartialAtomID > & atoms
) const
{
	if (
			! ( rsd.is_lower_terminus() && mainchain_torsion_ == 2 ) &&
			! ( rsd.is_upper_terminus() && mainchain_torsion_ == n_tors_ ) ) {
		conformation::insert_partial_atom_ids_for_mainchain_torsion(
			rsd, mainchain_torsion_, atoms );
	}
}


}
}
}
