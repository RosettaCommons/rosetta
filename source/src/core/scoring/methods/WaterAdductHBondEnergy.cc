// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/WaterAdductHBondEnergy.cc
/// @brief
/// @author Jim Havranek
/// @author Phil Bradley
/// @author Andrew Leaver-Fay


// Unit headers
#include <core/scoring/methods/WaterAdductHBondEnergy.hh>
#include <core/scoring/methods/WaterAdductHBondEnergyCreator.hh>

// Package headers
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/WaterAdductHBondPotential.hh>
#include <core/scoring/hbonds/HBondSet.hh>
// AUTO-REMOVED #include <core/scoring/hbonds/hbonds.hh>

// Project headers
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


// Utility headers


// C++


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the WaterAdductHBondEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
WaterAdductHBondEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new WaterAdductHBondEnergy;
}

ScoreTypes
WaterAdductHBondEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( h2o_hbond );
	return sts;
}


WaterAdductHBondEnergy::WaterAdductHBondEnergy() :
	parent( methods::EnergyMethodCreatorOP( new WaterAdductHBondEnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_WaterAdductHBondPotential() )
{}


/// clone
EnergyMethodOP
WaterAdductHBondEnergy::clone() const
{
	return new WaterAdductHBondEnergy();
}

///
void
WaterAdductHBondEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	using EnergiesCacheableDataType::H2O_HBOND_SET;
	pose.update_residue_neighbors();
	hbonds::HBondSetOP h2o_hbond_set( new hbonds::HBondSet( pose.total_residue() ) );
	potential_.fill_h2o_hbond_set( pose, *h2o_hbond_set );
	pose.energies().data().set( H2O_HBOND_SET, h2o_hbond_set );
}

///
void
WaterAdductHBondEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}

void
WaterAdductHBondEnergy::setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const
{
	pose.update_residue_neighbors();
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///
void
WaterAdductHBondEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & ,  //pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{

	emap[ h2o_hbond ] += potential_.water_adduct_hbond_score( rsd1, rsd2 );

}


void
WaterAdductHBondEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
 	) const
{
	using EnergiesCacheableDataType::H2O_HBOND_SET;
	/// f1 and f2 are zeroed
	hbonds::HBondSet const & hbond_set
	( static_cast< hbonds::HBondSet const & >
	( pose.energies().data().get( H2O_HBOND_SET ) ) );
	Vector f1,f2;
	get_atom_h2o_hbond_derivative( atom_id, hbond_set, weights, f1, f2 );
	F1 += f1;
	F2 += f2;
}


///////////////////////////////////////////////////////////////////////////////
/// \brief  Get the f1 and f2 contributions from all hbonds involving this atom

void
WaterAdductHBondEnergy::get_atom_h2o_hbond_derivative(
	id::AtomID const & atom,
	hbonds::HBondSet const & hbond_set,
	EnergyMap const & weights,
	Vector & f1,
	Vector & f2
) const
{
	f1 = Vector(0.0);
	f2 = Vector(0.0);

	utility::vector1< hbonds::HBondCOP > const & hbonds
		( hbond_set.atom_hbonds( atom ) );

	for ( Size i=1; i<= hbonds.size(); ++i ) {
		hbonds::HBond const & hbond( *hbonds[ i ] );
		Real sign_factor( 0.0 );

		// This part is different from the straight hbond version
		// since there is no real hydrogen when water adducts donate
		if ( hbond.atom_is_acceptor( atom ) ) {
			sign_factor = -1.0;
		} else {
			sign_factor = 1.0;
		}

//		std::cout << "Processing h2o hbond with energy" << hbond.energy() << std::endl;

		// get the appropriate type of hbond weight
		Real const weight ( sign_factor * hbond.weight() * weights[ h2o_hbond ] );
//		std::cout << "Applying weight " << weight << std::endl;
//		std::cout << "sign_factor " << sign_factor << std::endl;
//		std::cout << "hbond stored weight " << hbond.weight() << std::endl;
//		std::cout << "stupid type weight " << weights[ h2o_hbond ]  << std::endl;
		f1 += weight * hbond.derivs().h_deriv.f1();
		f2 += weight * hbond.derivs().h_deriv.f2();
//		std::cout << "F1 is " <<
//		f1[0] << " " <<
//		f1[1] << " " <<
//		f1[2] << " " <<
//		std::endl;
//		std::cout << "F2 is " <<
//		f2[0] << " " <<
//		f2[1] << " " <<
//		f2[2] << " " <<
//		std::endl;
	}
}



/// @brief distance cutoff
Distance
WaterAdductHBondEnergy::atomic_interaction_cutoff() const
{
	return 5.5; // -- temporary hack to allow us to use the standard neighbor array
}

/// @brief WaterAdductHBondEnergy
void
WaterAdductHBondEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}
core::Size
WaterAdductHBondEnergy::version() const
{
	return 1; // Initial versioning
}



} // methods
} // scoring
} // core
