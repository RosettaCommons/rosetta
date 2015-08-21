// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/chemical_shift/RNA_ChemicalShiftEnergy.cc
/// @brief  Energy score based on the agreement between experimentally determined and theoretically calculated NMR chemical_shift
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)


// Unit headers
#include <core/scoring/rna/chemical_shift/RNA_ChemicalShiftEnergy.hh>
#include <core/scoring/rna/chemical_shift/RNA_ChemicalShiftEnergyCreator.hh>

#include <core/scoring/ScoringManager.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>


// Utility headers

//Auto Headers
#include <core/id/AtomID.hh>
#include <core/scoring/EnergyMap.hh>

////////////////////////////////////////////////////////


// C++

namespace core {
namespace scoring {
namespace rna {
namespace chemical_shift {


/// @details This must return a fresh instance of the RNA_ChemicalShiftEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
RNA_ChemicalShiftEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new RNA_ChemicalShiftEnergy );
}

ScoreTypes
RNA_ChemicalShiftEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( rna_chem_shift );
	return sts;
}


/// c-tor
RNA_ChemicalShiftEnergy::RNA_ChemicalShiftEnergy():
	parent( methods::EnergyMethodCreatorOP( methods::EnergyMethodCreatorOP( new RNA_ChemicalShiftEnergyCreator ) ) ),
	rna_chemical_shift_potential_( ScoringManager::get_instance()->get_RNA_ChemicalShiftPotential() )
{}

/// clone
methods::EnergyMethodOP
RNA_ChemicalShiftEnergy::clone() const
{
	return methods::EnergyMethodOP( new RNA_ChemicalShiftEnergy );
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
void
RNA_ChemicalShiftEnergy::setup_for_scoring( pose::Pose &, ScoreFunction const & ) const
{


	//MIGHT BENEFIT FROM SOME BOOKKEEPING such as MAPPING atom_namerealatomdata_index/ to atom_index for fast look-up;

}


/////////////////////////////////////////////////////////////////////////////
void
RNA_ChemicalShiftEnergy::setup_for_derivatives( pose::Pose &, ScoreFunction const & ) const
{


	//MIGHT BENEFIT FROM SOME BOOKKEEPING such as MAPPING atom_namerealatomdata_index/ to atom_index for fast look-up;
	//or precalculate the calc_chem_shift values of the atoms with CS_data
	//this assumes that the pose conformation does change between the call to setup_for_derivative() and the calls to eval_atom_derivative()

}


///////////////////////////////////////////////////////////////////////////////
void
RNA_ChemicalShiftEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {

	rna_chemical_shift_potential_.finalize_total_energy( pose, totals );

} // finalize_total_energy

///////////////////////////////////////////////////////////////////////////////
void
RNA_ChemicalShiftEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{

	rna_chemical_shift_potential_.eval_atom_derivative( atom_id, pose, domain_map, weights, F1, F2 );


} // eval atom derivative

core::Size
RNA_ChemicalShiftEnergy::version() const
{
	return 1; // Initial versioning
}


} //chemical_shift
} //rna
} //scoring
} //core
