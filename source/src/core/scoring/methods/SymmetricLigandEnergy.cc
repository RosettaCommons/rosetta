// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SymmetricLigandEnergy.cc
/// @brief  Dunbrack energy method implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <core/scoring/methods/SymmetricLigandEnergy.hh>
#include <core/scoring/methods/SymmetricLigandEnergyCreator.hh>

// Package Headers
#include <core/scoring/EnergyMap.hh>
//#include <core/scoring/ScoringManager.hh>

#include <core/scoring/ScoreType.hh>

// Project headers
#include <core/pose/Pose.hh>

// Utility headers
#include <core/conformation/Residue.hh>

#include <core/id/AtomID.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the SymmetricLigandEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
SymmetricLigandEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new SymmetricLigandEnergy );
}

ScoreTypes
SymmetricLigandEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( sym_lig );
	return sts;
}

////////////////////////////////////////////////////////////////////////////
// SymmetricLigandEnergy
////////////////////////////////////////////////////////////////////////////

SymmetricLigandEnergy::SymmetricLigandEnergy() :
	parent( methods::EnergyMethodCreatorOP( new SymmetricLigandEnergyCreator ) )
{}


SymmetricLigandEnergy::~SymmetricLigandEnergy() {}

/// clone
EnergyMethodOP
SymmetricLigandEnergy::clone() const
{
	return EnergyMethodOP( new SymmetricLigandEnergy );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////


void
SymmetricLigandEnergy::residue_energy(
	conformation::Residue const & rsd,
	pose::Pose const &,
	EnergyMap & emap
) const
{
	numeric::xyzVector<core::Real> target(0,0,4);
	if( rsd.has("CEN") ) { // centroid
		if( /*"HIS" == rsd.name3() &&*/ rsd.xyz("CEN").z() > 0.0 ) {
			numeric::xyzVector<core::Real> cen  = rsd.xyz("CEN");
			numeric::xyzVector<core::Real> base = rsd.xyz("CA" );
			Real score = numeric::min( 0.0, (-10.0 / (cen.distance(target) + 1.0) + 1.0) );
			score *= dot( (cen-base).normalized(),(target-cen).normalized());
			emap[sym_lig] += score;

			Real penalty = rsd.xyz("CA" ).y();
			penalty *= penalty;
			emap[sym_lig] += penalty * 0.01;

		}
		emap[sym_lig] += numeric::min( 0.0, rsd.xyz("CEN").distance(target) - 20.0 ) / 20;
	} else { // full atom
		if( "HIS" == rsd.name3() && rsd.xyz("NE2").z() > 0.0 ) {
			numeric::xyzVector<core::Real> cen  = rsd.xyz("NE2");
			numeric::xyzVector<core::Real> base = (rsd.xyz("CG")+rsd.xyz("ND1"))/2.0;
			Real score = numeric::min( 0.0, (-8.0 / (cen.distance(target) + 1.0) + 1.0) );
			// no orientation here yet because I don't wanna do derives for it
			// if( score != 0.0 ) {
			// 	std::cerr << "FA HIS BONUS " << score << std::endl;
			// }
			emap[sym_lig] += score;
		}
	}


}


void
SymmetricLigandEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const & /*domain_map*/,
	ScoreFunction const & /*sfxn*/,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const {
	if( "HIS" != pose.residue(id.rsd()).name3() ) return;

	if( "NE2" != pose.residue(id.rsd()).atom_name(id.atomno()) ) return;

	// std::cerr << "SymmetricLigandEnergy deriv " << id << std::endl;

	numeric::xyzVector<core::Real> target(0,0,4.0);
	//if( pose.xyz(id).z() < 0.0 ) target *= -1.0;

	if( pose.xyz(id).distance(target) > 5.0 ) return;

  	numeric::xyzVector<core::Real> atom_x = pose.xyz(id);
	core::Real mag = pose.xyz(id).distance(target) + 1.0;
	mag = 6.0 / (mag*mag);
	numeric::xyzVector<core::Real> const f2( mag * ( pose.xyz(id) - target ).normalized() );
	numeric::xyzVector<core::Real> const atom_y = atom_x - f2;   // a "fake" atom in the direcion of the gradient
	numeric::xyzVector<core::Real> const f1( atom_x.cross( atom_y ) );

	F1 += weights[ sym_lig ] * f1;
	F2 += weights[ sym_lig ] * f2;

}


/// @brief SymmetricLigandEnergy is context independent; indicates that no context graphs are required
void
SymmetricLigandEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const
{}
core::Size
SymmetricLigandEnergy::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

