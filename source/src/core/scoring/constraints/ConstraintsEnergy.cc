// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/ConstraintsEnergy.hh
/// @brief  Constraints Energy Method declaration
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/scoring/constraints/ConstraintsEnergy.hh>
#include <core/scoring/constraints/ConstraintsEnergyCreator.hh>

// Package headers
#include <core/scoring/constraints/ConstraintEnergyContainer.hh>
#include <core/scoring/constraints/CstMinimizationData.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/XYZ_Func.hh>

#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/ScoreFunction.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
//#include <core/scoring/ScoringManager.hh>

// utility headers
#include <basic/Tracer.hh>
#include <basic/prof.hh>

#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace constraints {

/// @details This must return a fresh instance of the ConstraintsEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ConstraintsEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new ConstraintsEnergy;
}

ScoreTypes
ConstraintsEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( atom_pair_constraint );
	sts.push_back( angle_constraint );
	sts.push_back( dihedral_constraint );
	sts.push_back( constant_constraint );
	sts.push_back( coordinate_constraint );
	sts.push_back( dof_constraint );
	sts.push_back( res_type_constraint );
	sts.push_back( res_type_linking_constraint );
	sts.push_back( backbone_stub_constraint );
	sts.push_back( backbone_stub_linear_constraint );
	sts.push_back( big_bin_constraint );
	sts.push_back( dunbrack_constraint );
	sts.push_back( pocket_constraint );
	sts.push_back( rna_bond_geometry );
	sts.push_back( metalhash_constraint );
	return sts;
}


static thread_local basic::Tracer tr( "core.scoring.ConstraintsEnergy" );

ConstraintsEnergy::ConstraintsEnergy() : parent( methods::EnergyMethodCreatorOP( new ConstraintsEnergyCreator ) ) {}

ConstraintsEnergy::~ConstraintsEnergy() {}

methods::EnergyMethodOP
ConstraintsEnergy::clone() const
{
	return new ConstraintsEnergy;
}

methods::LongRangeEnergyType
ConstraintsEnergy::long_range_type() const { return methods::constraints_lr; }


bool
ConstraintsEnergy::defines_residue_pair_energy(
	pose::Pose const & pose,
	Size res1,
	Size res2
) const
{
	return pose.constraint_set()->residue_pair_constraint_exists( res1, res2 );
}

/// @brief Evaluate constraint residue_pair_energy.  Defers entirely to cst_set.
void
ConstraintsEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	if ( pose.constraint_set() == 0 ) return;
	//	PROF_START( basic::CONSTRAINT_SCORE );
	pose.constraint_set()->residue_pair_energy( rsd1, rsd2, pose, sfxn, emap );
	//	PROF_STOP( basic::CONSTRAINT_SCORE );
}

bool
ConstraintsEnergy::minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

bool
ConstraintsEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

void
ConstraintsEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const & min_data,
	pose::Pose const & , //pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	CstMinimizationDataCOP res_cst = static_cast< CstMinimizationData const * > ( min_data.get_data( cst_respair_data )() );
	if (!res_cst) return;
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	res_cst->constraints().residue_pair_energy( rsd1, rsd2, sfxn.weights(), emap );
}

bool
ConstraintsEnergy::defines_score_for_residue_pair(
	conformation::Residue const &,
	conformation::Residue const &,
	bool res_moving_wrt_eachother
) const
{
	return res_moving_wrt_eachother;
}

void
ConstraintsEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & minmap,
	ResSingleMinimizationData & res_data_cache
) const
{
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	pose.constraint_set()->setup_for_minimizing_for_residue( rsd, pose, sfxn, minmap, res_data_cache );
}

void
ConstraintsEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & minmap,
	ResSingleMinimizationData const & res1_data_cache,
	ResSingleMinimizationData const & res2_data_cache,
	ResPairMinimizationData & respair_data_cache
) const
{
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	pose.constraint_set()->setup_for_minimizing_for_residue_pair(
		rsd1, rsd2, pose, sfxn, minmap,
		res1_data_cache, res2_data_cache, respair_data_cache );
}

bool
ConstraintsEnergy::requires_a_setup_for_scoring_for_residue_pair_opportunity( pose::Pose const & ) const
{
	return true;
}

void
ConstraintsEnergy::setup_for_scoring_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	pose::Pose const &,
	ScoreFunction const & sfxn,
	ResPairMinimizationData & data_cache
) const
{
	CstMinimizationDataOP res_cst = static_cast< CstMinimizationData * > ( data_cache.get_data( cst_respair_data )() );
	if (!res_cst) return;
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	res_cst->constraints().setup_for_scoring( func::ResiduePairXYZ( rsd1, rsd2 ), sfxn );

}

bool
ConstraintsEnergy::requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const
{
	return true;
}

void
ConstraintsEnergy::setup_for_derivatives_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	pose::Pose const &,
	ScoreFunction const & sfxn,
	ResPairMinimizationData & data_cache
) const
{
	CstMinimizationDataOP res_cst = static_cast< CstMinimizationData * > ( data_cache.get_data( cst_respair_data )() );
	if (!res_cst) return;
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	res_cst->constraints().setup_for_derivatives( func::ResiduePairXYZ( rsd1, rsd2 ), sfxn );
}

/*void
ConstraintsEnergy::eval_atom_derivative_for_residue_pair(
	Size const atom_index,
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,// minsingle_data1,
	ResSingleMinimizationData const &,// minsingle_data2,
	ResPairMinimizationData const & minpair_data,
	pose::Pose const & pose, // provides context
	kinematics::DomainMap const &,// domain_map,
	ScoreFunction const & ,//sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	CstMinimizationDataCOP res_cst = static_cast< CstMinimizationData const * > ( minpair_data.get_data( cst_respair_data )() );
	if (!res_cst) return;
	res_cst->constraints().eval_respair_atom_derivative(
		id::AtomID( atom_index, rsd1.seqpos() ), rsd1, rsd2, weights, F1, F2 );
}*/

void
ConstraintsEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	pose::Pose const &,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	CstMinimizationDataCOP res_cst = static_cast< CstMinimizationData const * > ( min_data.get_data( cst_respair_data )() );
	if (!res_cst) return;
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	for ( Size ii = 1; ii <= rsd1.natoms(); ++ii ) {
		res_cst->constraints().eval_respair_atom_derivative(
			id::AtomID( ii, rsd1.seqpos() ), rsd1, rsd2, weights, r1_atom_derivs[ ii ].f1(), r1_atom_derivs[ ii ].f2() );
	}
	for ( Size ii = 1; ii <= rsd2.natoms(); ++ii ) {
		res_cst->constraints().eval_respair_atom_derivative(
			id::AtomID( ii, rsd2.seqpos() ), rsd2, rsd1, weights, r2_atom_derivs[ ii ].f1(), r2_atom_derivs[ ii ].f2() );
	}

}



/////////////////////////////////////////////////////////////////////////////
// methods for ContextIndependentTwoBodyEnergies
/////////////////////////////////////////////////////////////////////////////

/// @brief Intraresidue constraints can exist; the ConstraintsEnergy class
/// cannot tell whether or not any intraresidue constraints exist when only
/// given the weight set, so for safety it says "yes, I define intraresidue
/// energies" and guarantees they will be evaluated properly if they do exist.
bool
ConstraintsEnergy::defines_intrares_energy( EnergyMap const & ) const
{
	return true;
}

/// @brief Evaluate the intra-residue constraint energy for a given residue
void
ConstraintsEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	if ( pose.constraint_set() == 0 ) return;
	//	PROF_START( basic::CONSTRAINT_SCORE );
	pose.constraint_set()->eval_intrares_energy( rsd, pose, sfxn, emap );
	//	PROF_STOP( basic::CONSTRAINT_SCORE );
}

/// @brief request of minimization routines that they use the extended intraresidue energy
/// interface
bool
ConstraintsEnergy::use_extended_intrares_energy_interface() const
{
	return true;
}

/// @brief Evaluate the intra-residue energies using ConstraintCOPs cached in the data_cache object
void
ConstraintsEnergy::eval_intrares_energy_ext(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & data_cache,
	pose::Pose const &,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) const
{
	CstMinimizationDataCOP res_cst = static_cast< CstMinimizationData const * > ( data_cache.get_data( cst_res_data )() );
	if (!res_cst) return;
	res_cst->constraints().intra_residue_energy( rsd, sfxn.weights(), emap );
}

bool
ConstraintsEnergy::requires_a_setup_for_scoring_for_residue_opportunity( pose::Pose const & ) const
{
	return true;
}


void
ConstraintsEnergy::setup_for_scoring_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const &,
	scoring::ScoreFunction const & sfxn,
	ResSingleMinimizationData & min_data
) const
{
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	CstMinimizationDataOP res_cst = static_cast< CstMinimizationData * > ( min_data.get_data( cst_res_data )() );
	if (!res_cst) return;
	res_cst->constraints().setup_for_scoring( func::ResidueXYZ( rsd ), sfxn );
}


bool
ConstraintsEnergy::requires_a_setup_for_derivatives_for_residue_opportunity( pose::Pose const & ) const
{
	return true;
}


void
ConstraintsEnergy::setup_for_derivatives_for_residue(
	conformation::Residue const & rsd,
	pose::Pose const &,
	scoring::ScoreFunction const & sfxn,
	ResSingleMinimizationData & min_data
) const
{
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	CstMinimizationDataOP res_cst = static_cast< CstMinimizationData * > ( min_data.get_data( cst_res_data )() );
	if (!res_cst) return;
	res_cst->constraints().setup_for_derivatives( func::ResidueXYZ( rsd ), sfxn );
}

/*void
ConstraintsEnergy::eval_intrares_atom_derivative(
	Size const atom_index,
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & min_data,
	pose::Pose const &,
	kinematics::DomainMap const &,
	ScoreFunction const & ,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	CstMinimizationDataCOP res_cst = static_cast< CstMinimizationData const * > ( min_data.get_data( cst_res_data )() );
	if (!res_cst) return;
	res_cst->constraints().eval_intrares_atom_derivative(
		id::AtomID( atom_index, rsd.seqpos() ), rsd, weights, F1, F2 );
}*/

void
ConstraintsEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const & min_data,
	pose::Pose const &,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	CstMinimizationDataCOP res_cst = static_cast< CstMinimizationData const * > ( min_data.get_data( cst_res_data )() );
	if (!res_cst) return;
	for ( Size ii = 1; ii <= rsd.natoms(); ++ii ) {
		res_cst->constraints().eval_intrares_atom_derivative(
			id::AtomID( ii, rsd.seqpos() ), rsd, weights,
			atom_derivs[ ii ].f1(), atom_derivs[ ii ].f2() );
	}

}

bool
ConstraintsEnergy::defines_intrares_dof_derivatives( pose::Pose const & ) const
{
	// Only a handful of constraints actually define dof derivatives, but it's hard to know
	// in advance if any of those constraints are in the ConstraintSet; maybe if dof-deriv
	// defining constraints had to be held in a special place in the constraint set, they'd
	// be easier to find and access?
	return true;
}

/// @brief Evaluate the DOF derivative for a particular residue.  The Pose merely serves as context,
/// and the input residue is not required to be a member of the Pose.
Real
ConstraintsEnergy::eval_intraresidue_dof_derivative(
	conformation::Residue const &,// rsd,
	ResSingleMinimizationData const &,// min_data,
	id::DOF_ID const &,// dof_id,
	id::TorsionID const &,// torsion_id,
	pose::Pose const &,// pose,
	ScoreFunction const &,// sfxn,
	EnergyMap const &// weights
) const
{
	return 0.0;
}

// @brief store a constraint graph in the pose.energies object, or if one is
// already present, make sure that the constraint graph accurately reflects
// the constraint set.  The constraint energy container keeps the revision id
// for a particular constraint set.  If a pose constraint set changes, either
// in identity (i.e. a new object) or in composition (i.e. an addition/removal
// of a constraint) then thene ConstraintsEnergy object throws away the old
// ConstraintsEnergyContainer and creats a new one... somewhat wasteful.
//
// In the future, either the CE or the ConstraintSet will have to intelligently
// decide between a constraint graph or a constraint table.  In the case of
// dense residue pair constraints (i.e. all pairs), a table would be more efficient.
void
ConstraintsEnergy::prepare_constraints_energy_container( pose::Pose & pose ) const
{
	using namespace methods;

	//do nothing if there are no constraints
	if ( pose.constraint_set() == 0 ) return; //this will never happen: since empty constraint_set is created as soon as method constraint_set() is called
	// however, if one jumps out here after asking !has_constraints() it breaks in ScoringFunction.cc:604 because the
	// long_rang_container is not setup for constraints.
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	Energies & energies( pose.energies() );
	bool create_new_cstcontainer( false );
	if ( energies.long_range_container( constraints_lr ) == 0 ) {
		create_new_cstcontainer = true; // pbmod
	} else {
		LREnergyContainerOP lrc = energies.nonconst_long_range_container( constraints_lr );
		CstEnergyContainerOP cec( static_cast< CstEnergyContainer * > ( lrc.get() ) );
		if ( ! cec->matches( pose.constraint_set()) ) {
			create_new_cstcontainer = true;
		}
	}

	if ( create_new_cstcontainer ) {
		CstEnergyContainerOP new_cec = new CstEnergyContainer( pose );
		energies.set_long_range_container( constraints_lr, new_cec );
	}

}

void
ConstraintsEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & sfxn ) const
{
	prepare_constraints_energy_container( pose );
	if ( pose.constraint_set() ) {
		pose.constraint_set()->setup_for_scoring( pose, sfxn );  //fpd
	}
}

void
ConstraintsEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	prepare_constraints_energy_container( pose );
}

///@brief Setup constraint-set specific derivatives
void
ConstraintsEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sfxn ) const
{
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	pose.constraint_set()->setup_for_derivatives( pose, sfxn );  //fpd
}

bool
ConstraintsEnergy::defines_high_order_terms( pose::Pose const & pose ) const
{
	return pose.constraint_set()->has_non_residue_pair_constraints();
}


/// @brief called at the end of energy evaluation; allows
/// for the evaluation of constraints that cannot be
/// decomposed into residue pairs.  Defers entirely to cst_set.
void
ConstraintsEnergy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	EnergyMap & totals
) const
{
	if ( pose.constraint_set() == 0 ) return;
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	pose.constraint_set()->eval_non_residue_pair_energy( pose, sfxn, totals );
}


/// called during gradient-based minimization inside dfunc
/**
	 F1 and F2 are not zeroed -- contributions from this atom are
	 just summed in
**/
void
ConstraintsEnergy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &, // domain_map,
	ScoreFunction const & sfxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	if ( !pose.constraint_set()->has_constraints() ) return;
	//	basic::ProfileThis doit( basic::CONSTRAINT_SCORE );
	pose.constraint_set()->eval_multibody_atom_derivative( id, pose, sfxn, weights, F1, F2 );
}

/// uses the dof constraints -- depricated interface
Real
ConstraintsEnergy::eval_dof_derivative(
	id::DOF_ID const &,// id,
	id::TorsionID const &,// tor,
	pose::Pose const &,// pose,
	ScoreFunction const &,// scorefxn,
	EnergyMap const & // weights
) const
{
	// if ( pose.constraint_set() == 0 ) return 0.0;
	//this -- oddly -- never happens, since the pose will make an empty set as soon as you ask.

	// instead...
	//if ( ! pose.constraint_set()->has_constraints() ) return 0.0;
	//PROF_START( basic::CONSTRAINT_SCORE );
	//Real result ( pose.constraint_set()->eval_dof_derivative( id, tor, pose, scorefxn, weights ) );
	//PROF_STOP( basic::CONSTRAINT_SCORE );
	return 0.0;
}


/// @brief constraints are context independent
void
ConstraintsEnergy::indicate_required_context_graphs(
	utility::vector1< bool > &
) const
{}
core::Size
ConstraintsEnergy::version() const
{
	return 1; // Initial versioning
}





} // constraints
} // scoring
} // core
