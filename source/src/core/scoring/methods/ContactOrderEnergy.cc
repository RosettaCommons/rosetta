// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/methods/ContactOrderEnergy.cc
/// @brief calculates the contact order of a given conformation average sequence.
/// @detailed contact order is defined as the average sequence separation of
/// residues that are in contact.
///
/// @author James Thompson

// Unit Headers
#include <core/scoring/methods/ContactOrderEnergy.hh>
#include <core/scoring/methods/ContactOrderEnergyCreator.hh>

// Package Headers

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/scoring/EnergyMap.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.hh>

#include <basic/Tracer.hh>

//// C++ headers

static thread_local basic::Tracer tr( "core.scoring.methods.ContactOrderEnergy" );

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the ContactOrderEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ContactOrderEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new ContactOrderEnergy );
}

ScoreTypes
ContactOrderEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( co );
	return sts;
}


/// c-tor
ContactOrderEnergy::ContactOrderEnergy() :
	parent( methods::EnergyMethodCreatorOP( new ContactOrderEnergyCreator ) )
{}

/// clone
EnergyMethodOP
ContactOrderEnergy::clone() const
{
    return EnergyMethodOP( new ContactOrderEnergy() );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void ContactOrderEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {

  totals[ co ] = calculate_contact_order( pose );

}

core::Real
ContactOrderEnergy::calculate_contact_order( pose::Pose const & pose ) const
{

// tex: below is old code for calculating contact order from rosetta++
// 	int nco = 0;
// 	for ( int i = 1; i <= pose.total_residue(); ++i ) {
// 		if ( is_protein(res(i)) || is_nonnatural(res(i)) ) { /// <---
// 			for ( int kk = 1, kke = cen12up(i); kk <= kke; ++kk ) { /// <---
// 				int j = cen_list(kk,i); /// <---
// 				if ( cendist(i,j) < 64.0 && std::abs(j-i) > 2 ) { /// <---
// 					co += std::abs(j-i);
// 					++nco;
// 				}
// 			}
// 		}
// 	}
// notes on what things mean:
// cenlist(*,i) contains the list of residues within 12A of residue i
// cen12up(i) is the number of residues in cen_list(*,i)

	using core::Real;
	using core::Size;
	using core::Vector;

	Real co_score = 0.0;
	Size n_in_contact = 0;
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		Vector const v1( pose.residue(i).nbr_atom_xyz() );

		for ( Size j = i + 3; j <= pose.total_residue(); ++j ) {
			Vector const v2( pose.residue(j).nbr_atom_xyz() );
			if ( v1.distance_squared( v2 ) < 64.0 ) {
				co_score += j - i;
				++n_in_contact;
			}
		}
	}

	if ( n_in_contact > 0 ) {
		co_score /= static_cast< Real >( n_in_contact );
	}

	return co_score;
}
core::Size
ContactOrderEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core
