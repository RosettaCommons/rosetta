// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/PeptideBondEnergy.cc
/// @author James Thompson <tex@u.washington.edu>

// Unit headers
#include <core/scoring/methods/PeptideBondEnergy.hh>
#include <core/scoring/methods/PeptideBondEnergyCreator.hh>

// Package headers
//#include <core/scoring/ScoringManager.hh>

// Project headers
// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/scoring/func/GaussianFunc.hh>

#include <utility/vector1.hh>


// Numeric headers


namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the PeptideBondEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
PeptideBondEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new PeptideBondEnergy;
}

ScoreTypes
PeptideBondEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( peptide_bond );
	return sts;
}


PeptideBondEnergy::PeptideBondEnergy() : parent( methods::EnergyMethodCreatorOP( new PeptideBondEnergyCreator ) ) {}

/// called at the end of energy evaluation
/// In this case (PeptideBondEnergy), all the calculation is done here

//void
//PeptideBondEnergy::finalize_total_energy(
// pose::Pose & pose,
// ScoreFunction const &,
// EnergyMap & totals
//) const
//{
//	using conformation::Residue;
//
//	Real const mean( 1.325883 );
//	Real const sdev( 0.012547 );
//	Real total_dev(0.0);
//	core::scoring::constraints::GaussianFunc gfunc( mean, sdev );
//
//	std::string const bbN_( "N" );
//	std::string const bbC_( "C" );
//
//	typedef core::conformation::ResidueOPs ResidueOPs;
//	for ( ResidueOPs::iterator this_res = pose.res_begin(),
//				next_res = this_res + 1,
//				end = pose.res_end();
//				next_res != end; ++this_res, ++next_res
//	) {
//		total_dev += gfunc.func(
//			(*this_res)->xyz( bbC_ ).distance( (*next_res)->xyz( bbN_ ) )
//		);
//	}
//	totals[ peptide_bond ] = total_dev;
//}

void
PeptideBondEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & /* pose */,
	ScoreFunction const &,
	EnergyMap & emap
) const {
	using conformation::Residue;

	Real const mean( 1.325883 );
	Real const sdev( 0.012547 );
	Real total_dev(0.0);
	core::scoring::func::GaussianFunc gfunc( mean, sdev );

	std::string const bbN_( "N" );
	std::string const bbC_( "C" );

	// make certain that both require atoms are present.
	if ( !rsd1.type().has(bbN_) || !rsd2.type().has(bbC_) )
		return;
	// make sure we're bonded and in relatively the right sequence orientation
	if ( !rsd1.is_bonded(rsd2) || ( rsd1.seqpos() > rsd2.seqpos() ) )
		return;

	total_dev += gfunc.func(
		rsd1.xyz( bbC_ ).distance( rsd2.xyz( bbN_ ) )
	);
	emap[ peptide_bond ] += total_dev;
} // residue_pair_energy

	/// called during gradient-based minimization inside dfunc
	/**
		 F1 and F2 are not zeroed -- contributions from this atom are
		 just summed in
	**/

void
PeptideBondEnergy::eval_atom_derivative(
	id::AtomID const & /* id */,
	pose::Pose const & /* pose */,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const &, // sfxn,
	EnergyMap const & /* weights */,
	Vector & /* F1 */,
	Vector & /* F2 */
) const
{} // need to fill this in if we ever want to minimize with it!

/// @brief Chainbreak Energy is context independent and thus indicates that no context graphs need to
/// be maintained by class Energies
void
PeptideBondEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & /*context_graphs_required*/
) const
{}

Distance
PeptideBondEnergy::atomic_interaction_cutoff() const {
	return 3.0;
}
core::Size
PeptideBondEnergy::version() const
{
	return 1; // Initial versioning
}


} // namespace methods
} // namespace scoring
} // namespace core
