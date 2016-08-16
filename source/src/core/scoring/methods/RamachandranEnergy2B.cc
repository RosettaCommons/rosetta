// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/RamachandranEnergy.cc
/// @brief  Ramachandran energy method class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/methods/RamachandranEnergy2B.hh>
#include <core/scoring/methods/RamachandranEnergy2BCreator.hh>

// Package Headers
#include <core/scoring/Ramachandran2B.hh>
#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <basic/options/option.hh>

// Utility headers
#include <numeric/conversions.hh>

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/vector1.hh>


// C++


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the RamachandranEnergy2B class,
/// never an instance already in use
methods::EnergyMethodOP
RamachandranEnergy2BCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new RamachandranEnergy2B );
}

ScoreTypes
RamachandranEnergy2BCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rama2b );
	return sts;
}


/// ctor
RamachandranEnergy2B::RamachandranEnergy2B() :
	parent( EnergyMethodCreatorOP( new RamachandranEnergy2BCreator ) ),
	potential_( ScoringManager::get_instance()->get_Ramachandran2B() )
{}

/// clone
EnergyMethodOP
RamachandranEnergy2B::clone() const
{
	return EnergyMethodOP( new RamachandranEnergy2B );
}

/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentOneBodyEnergies
/////////////////////////////////////////////////////////////////////////////

void
RamachandranEnergy2B::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &, // unnneeded
	ScoreFunction const &, // unneeded
	EnergyMap & emap
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( ! option[ score::ramaneighbors ] ) return;

	/// This is called for all nearby residue pairs, so first check to make sure that we've got an i, i+1 pair
	if ( rsd1.seqpos() + 1 != rsd2.seqpos() || rsd1.seqpos() != rsd2.seqpos() + 1 ) return;
	if ( rsd1.chain() != rsd2.chain() ) return;
	if ( ! rsd1.is_protein() || ! rsd2.is_protein() ) return;

	conformation::Residue const & lower_residue( rsd1.seqpos() < rsd2.seqpos() ? rsd1 : rsd2 );
	conformation::Residue const & upper_residue( rsd1.seqpos() < rsd2.seqpos() ? rsd2 : rsd1 );

	//// also need to treat cutpoints correctly.

	if ( ! lower_residue.is_lower_terminus() ) {
		emap[ rama2b ] += potential_.RamaE_Upper( lower_residue, upper_residue.aa() );
	}

	if ( ! upper_residue.is_upper_terminus() ) {
		emap[ rama2b ] += potential_.RamaE_Lower( upper_residue, lower_residue.aa() );
	}
}

bool
RamachandranEnergy2B::defines_intrares_energy( EnergyMap const & /*weights*/ ) const {
	return true;
}

/// @details fictional Cprev-Nnext distance + fudge.
Real
RamachandranEnergy2B::atomic_interaction_cutoff() const
{
	return 2.0;
}

void
RamachandranEnergy2B::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const &, // unused,
	ScoreFunction const &, // unused,
	EnergyMap & emap
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !rsd.is_protein() ) return;

	if ( rsd.is_terminus() ) {
		/// -- no op --  Rama does not have a defined score for termini
		/// emap[ rama ] += potential_.Rama_E( rsd ); // add in neighbor-independent rama scores for termini
	}  else if ( option[ score::ramaneighbors ] ) {
		emap[ rama2b ] -= potential_.RamaE( rsd ); // subtract double-counted neighbor-independent rama scores for mid residues
	} else {
		emap[ rama2b ] += potential_.RamaE( rsd ); // add the neighbor-independent rama score, since there's no double counting.
	}
}


Real
RamachandranEnergy2B::eval_dof_derivative(
	id::DOF_ID const &,// dof_id,
	id::TorsionID const & tor_id,
	pose::Pose const & pose,
	ScoreFunction const &,// sfxn,
	EnergyMap const & weights
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// ignore invalid or non-BB torsions
	if ( !tor_id.valid() || tor_id.type() != id::BB ) return 0.0;

	Real deriv(0.0);
	conformation::Residue const & rsd( pose.residue( tor_id.rsd() ) );
	if ( rsd.is_protein() && tor_id.torsion() <= 2 && ! rsd.is_terminus() ) {
		Real rama_score, drama_dphi, drama_dpsi;
		if ( option[ score::ramaneighbors ] ) {
			// Neighbor dependent rama score + derivatives.
			Size const seqpos( rsd.seqpos() );
			potential_.eval_rama_score_residue( rsd,
				pose.residue_type( seqpos - 1 ).aa(),
				pose.residue_type( seqpos + 1 ).aa(),
				rama_score, drama_dphi, drama_dpsi );
		} else {
			/// Neighbor independent rama score + derivatives
			potential_.eval_rama_score_residue( rsd,
				rama_score, drama_dphi, drama_dpsi );
		}
		deriv = ( tor_id.torsion() == 1 ? drama_dphi : drama_dpsi );
	}
	// note that the atomtree PHI dofs are in radians
	// use degrees since dE/dangle has angle in denominator
	return numeric::conversions::degrees( weights[ rama ] * deriv );
}

/// @brief Ramachandran Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
RamachandranEnergy2B::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
)
const
{}
core::Size
RamachandranEnergy2B::version() const
{
	return 1; // Initial versioning
}


} // methods
} // scoring
} // core

