// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/FA_GrpElecEnergy.cc
/// @brief  Electrostatic energy with a distance-dependant dielectric
/// @author Hahnbeom Park

// Unit headers
#include <core/scoring/elec/FA_GrpElecEnergy.hh>
#include <core/scoring/elec/FA_GrpElecEnergyCreator.hh>
#include <core/scoring/elec/GroupElec.hh>

// Package headers
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/types.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/hbonds/HBondSet.hh>

#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Project headers
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/conformation/RotamerSetBase.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Basic headers
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

static basic::Tracer TR("core.scoring.elec.FA_GrpElecEnergy");

namespace core {
namespace scoring {
namespace elec {


/// @details This must return a fresh instance of the FA_GrpElecEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
FA_GrpElecEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new FA_GrpElecEnergy( options ) );
}

ScoreTypes
FA_GrpElecEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( fa_grpelec );
	return sts;
}

FAElecContextData::FAElecContextData(){}
FAElecContextData::~FAElecContextData(){}

void 
FAElecContextData::initialize( Size const nres )
{
	n_.resize( nres );
	for( Size i = 1; i <= nres; ++i ) n_[i] = 0;
	dw_dr_.resize( nres );
	for( Size i = 1; i <= nres; ++i ) dw_dr_[i] = Vector( 0.0 );
	boundary_neighs_.resize( nres );
}


////////////////////////////////////////////////////////////////////////////
FA_GrpElecEnergy::FA_GrpElecEnergy( methods::EnergyMethodOptions const & options ):
	parent( methods::EnergyMethodCreatorOP( new FA_GrpElecEnergyCreator ) ),
	coulomb_( options ),
	groupelec_( options ),
	exclude_protein_protein_( options.exclude_protein_protein_fa_elec() ),
	exclude_monomer_( options.exclude_monomer_fa_elec() ),
	exclude_DNA_DNA_( options.exclude_DNA_DNA() ),
	intrares_scale_( options.intrares_elec_correction_scale() )
{
	coulomb_.initialize();
	groupelec_.initialize( coulomb() );
}


////////////////////////////////////////////////////////////////////////////
FA_GrpElecEnergy::FA_GrpElecEnergy( FA_GrpElecEnergy const & src ):
	parent( src ),
	coulomb_( src.coulomb() ),
	groupelec_( src.groupelec_ ),
	exclude_protein_protein_( src.exclude_protein_protein_ ),
	exclude_monomer_( src.exclude_monomer_ ),
	exclude_DNA_DNA_( src.exclude_DNA_DNA_ ),
	intrares_scale_( src.intrares_scale_ )
{
	coulomb_.initialize();
	groupelec_.initialize( coulomb() );
}


void
FA_GrpElecEnergy::initialize() {
	coulomb_.initialize();
	groupelec_.initialize( coulomb() );
}


/// clone
methods::EnergyMethodOP
FA_GrpElecEnergy::clone() const
{
	return methods::EnergyMethodOP( new FA_GrpElecEnergy( *this ) );
}

void
FA_GrpElecEnergy::setup_for_minimizing(
	pose::Pose & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	TR.Debug << "setup_for_minimization" << std::endl;

	set_nres_mono(pose);

	if ( pose.energies().use_nblist() ) {
		// stash our nblist inside the pose's energies object
		Energies & energies( pose.energies() );

		// setup the atom-atom nblist
		NeighborListOP nblist;
		Real const tolerated_motion = pose.energies().use_nblist_auto_update() ? option[ run::nblist_autoupdate_narrow ] : 1.5;
		Real const XX = coulomb().max_dis() + 2 * tolerated_motion;
		nblist = NeighborListOP( new NeighborList( min_map.domain_map(), XX*XX, XX*XX, XX*XX) );
		if ( pose.energies().use_nblist_auto_update() ) {
			nblist->set_auto_update( tolerated_motion );
		}
		// this partially becomes the EtableEnergy classes's responsibility
		nblist->setup( pose, sfxn, *this);
		energies.set_nblist( EnergiesCacheableDataType::ELEC_NBLIST, nblist );
	}

	TR.Debug << "done setup_for_minimization" << std::endl;

}

//
void
FA_GrpElecEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn ) const
{
	using namespace core::pose::datacache;

	TR.Debug << "setup_for_scoring" << std::endl;

  Eres_.resize( pose.total_residue(), 0.0 ); // temporary array

	set_nres_mono(pose);
	pose.update_residue_neighbors();

	if ( pose.energies().use_nblist() ) {
		NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::ELEC_NBLIST ) );
		nblist.prepare_for_scoring( pose, scfxn, *this );
	}

	FAElecContextDataOP data;
	if( pose.data().has( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA ) ){
		//data = static_cast< FAElecContextData * >
		//( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA )() );
		data = utility::pointer::static_pointer_cast< FAElecContextData >
			( pose.data().get_ptr( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA ) );
	} else {
		data = FAElecContextDataOP( new FAElecContextData() );
	}

	precalc_context( pose, data );

  pose.data().set( pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA, data );
}

void
FA_GrpElecEnergy::setup_for_derivatives( pose::Pose & pose,
																			 ScoreFunction const &scfxn
																			 ) const
{
	TR.Debug << "setup_for_deriv" << std::endl;

	set_nres_mono(pose);
	pose.update_residue_neighbors();

  //Eres_.resize( pose.total_residue(), 0.0 ); // temporary array
	setup_for_scoring( pose, scfxn );

	TR.Debug << "done: setup_for_deriv" << std::endl;
}


// The FA_ElectEnergy method stores a vector of rotamer trie objects in the Energies
// object for use in rapid rotamer/background energy calculations.  Overrides default
// do-nothing behavior.
void
FA_GrpElecEnergy::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	pose.update_residue_neighbors();
}

// @brief Creates a rotamer trie for the input set of rotamers and stores the trie
// in the rotamer set.
void
FA_GrpElecEnergy::prepare_rotamers_for_packing(
	pose::Pose const &,
	conformation::RotamerSetBase &
) const
{}

// @brief Updates the cached rotamer trie for a residue if it has changed during the course of
// a repacking
void
FA_GrpElecEnergy::update_residue_for_packing(
	pose::Pose &,
	Size
) const
{}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

///
void
FA_GrpElecEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{
	if ( pose.energies().use_nblist() ) return;

	TR.Debug << "res pair energy: " << rsd1.seqpos() << " " << rsd2.seqpos() << std::endl;

	using namespace etable::count_pair;

	Real score(0.0);
	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return;

	// use to avoid double counting with hbond_set
	/*
	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( EnergiesCacheableDataType::HBOND_SET )));
	*/

	FAElecContextDataCOP data = utility::pointer::static_pointer_cast< FAElecContextData const >
		( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA ) );

	Size res1( rsd1.seqpos()  );
	Size res2( rsd2.seqpos()  );
  Real const nb1 = data->get_n( res1 );
  Real const nb2 = data->get_n( res2 );

	score = groupelec().eval_respair_group_coulomb( rsd1, rsd2 );
	Real w( burial_weight( nb1 ) +  burial_weight( nb2 ) );
	score *= w;

	/*
	TR << rsd1.seqpos() << " " << rsd2.seqpos() << " "
		 << burial_weight( nb1 ) << " " << burial_weight( nb2 ) << " " << score << std::endl;
	*/

	emap[ fa_grpelec ] += score;

	TR.Debug << "done res pair energy: " << score << std::endl;
}

bool
FA_GrpElecEnergy::minimize_in_whole_structure_context( pose::Pose const & pose ) const
{
	return pose.energies().use_nblist_auto_update();
}

bool
FA_GrpElecEnergy::defines_score_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	bool res_moving_wrt_eachother
) const
{
	if ( rsd1.seqpos() == rsd2.seqpos() ) {
		return false; // turn it off for now; TODO: let intrares_elec_correction_scale control it
	} else if ( exclude_protein_protein_ && rsd1.is_protein() && rsd2.is_protein() ) {
		return false;
	} else if ( exclude_monomer_ && monomer_test( rsd1.seqpos(), rsd2.seqpos()) ) {
		return false;
	} else if ( exclude_DNA_DNA_ && rsd1.is_DNA() && rsd2.is_DNA() ) {
		return false;
	}

	return res_moving_wrt_eachother;
}


bool
FA_GrpElecEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

// this function is being used instead of residue_pair_energy in min_test;
void
FA_GrpElecEnergy::residue_pair_energy_ext(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResPairMinimizationData const &,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const
{

	TR.Debug << "residue_pair_energy_ext" << std::endl;

	if ( pose.energies().use_nblist_auto_update() ) return;

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return;
	bool intrares( rsd1.seqpos() == rsd2.seqpos() );
	Real score( 0.0 );

	// use to avoid double counting with hbond_set
	/*
	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( EnergiesCacheableDataType::HBOND_SET )));
	*/

	FAElecContextDataCOP data = utility::pointer::static_pointer_cast< FAElecContextData const >
		( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA ) );

	Size res1 = rsd1.seqpos();
	Size res2 = rsd2.seqpos();
  Real const &nb1 = data->get_n( res1 );
  Real const &nb2 = data->get_n( res2 );

  Real w = burial_weight( nb1 ) + burial_weight( nb2 );
	if( intrares ) w *= intrares_scale_;

	score = groupelec().eval_respair_group_coulomb( rsd1, rsd2 );

	/*
	TR << rsd1.seqpos() << " " << rsd2.seqpos() << " "
		 << burial_weight( nb1 ) << " " << burial_weight( nb2 ) << " " << score << std::endl;
	*/

	emap[ fa_grpelec ] += w*score;

	TR.Debug << "done residue_pair_energy_ext" << std::endl;
}

void
FA_GrpElecEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData &
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	etable::count_pair::CountPairFunctionCOP count_pair = get_count_pair_function( rsd1, rsd2 );
	assert( rsd1.seqpos() < rsd2.seqpos() );

	/*
	// update the existing nblist if it's already present in the min_data object
	ResiduePairNeighborListOP nblist( static_cast< ResiduePairNeighborList * > (pair_data.get_data( elec_pair_nblist )() ));
	if ( ! nblist ) nblist = new ResiduePairNeighborList;

	/// STOLEN CODE!
	Real const tolerated_narrow_nblist_motion = 0.75; //option[ run::nblist_autoupdate_narrow ];
	Real const XX2 = std::pow( coulomb().max_dis() + 2*tolerated_narrow_nblist_motion, 2 );

	nblist->initialize_from_residues( XX2, XX2, XX2, rsd1, rsd2, count_pair );

	pair_data.set_data( elec_pair_nblist, nblist );
	*/
}


void
FA_GrpElecEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;
	bool intrares( rsd1.seqpos() == rsd2.seqpos() );

	TR.Debug << "eval residue pair deriv" << std::endl;

	assert( utility::pointer::static_pointer_cast< FAElecContextData const > 
					( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA ) ));

	FAElecContextDataCOP data = utility::pointer::static_pointer_cast< FAElecContextData const >
	( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA ) );

	Size res1( rsd1.seqpos() );
	Size res2( rsd2.seqpos() );

	// use to avoid double counting with hbond_set
	/*
	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( EnergiesCacheableDataType::HBOND_SET )));
	*/

	Real elec_weight = weights[ fa_grpelec ];
	Real Erespair( 0.0 );
	Real w = burial_weight( data->get_n( res1 ) ) + burial_weight( data->get_n( res2 ) );
	if( intrares ) w *= intrares_scale_;

	Real total_weight = elec_weight*w;

	groupelec().eval_respair_group_derivatives( rsd1, rsd2,
																							r1_atom_derivs, r2_atom_derivs,
																							total_weight, Erespair );

	// mutable 
	Eres_[ rsd1.seqpos() ] += 0.5*Erespair;
	Eres_[ rsd2.seqpos() ] += 0.5*Erespair;

}

void
FA_GrpElecEnergy::eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
 ) const
{
	if ( pose.energies().use_nblist() ) return;

	TR.Debug << "intrares energy: " << rsd.seqpos() << std::endl;

	using namespace etable::count_pair;

	Real score(0.0);
	if ( ! defines_score_for_residue_pair(rsd, rsd, true) ) return;

	// use to avoid double counting with hbond_set
	/*
	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( EnergiesCacheableDataType::HBOND_SET )));
	*/

	FAElecContextDataCOP data = utility::pointer::static_pointer_cast< FAElecContextData const >
		( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA ) );

	Size res( rsd.seqpos() );
  Real const nb = data->get_n( res );

	score = groupelec().eval_respair_group_coulomb( rsd, rsd );
	Real w( 2.0*burial_weight( nb ) );
	w *= intrares_scale_;
	score *= w;

	/*
	TR << rsd.seqpos() << " "
		 << 2.0*burial_weight( nb ) << " " << intrares_scale_ << " " << score << std::endl;
	*/

	emap[ fa_grpelec ] += score;

	TR.Debug << "done intrares energy: " << score << std::endl;
}

void
FA_GrpElecEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;

	TR.Debug << "eval residue pair deriv" << std::endl;

	assert( utility::pointer::static_pointer_cast< FAElecContextData const > 
					( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA ) ));

	FAElecContextDataCOP data = utility::pointer::static_pointer_cast< FAElecContextData const >
	( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA ) );

	Size res( rsd.seqpos() );

  // get deriv on context 
	eval_context_derivatives( rsd, data, weights, r1_atom_derivs );

	// use to avoid double counting with hbond_set
	/*
	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( EnergiesCacheableDataType::HBOND_SET )));
	*/

	Real elec_weight = weights[ fa_grpelec ];
	Real Erespair( 0.0 );
	Real w = 2.0*burial_weight( data->get_n( res ) );
	w *= intrares_scale_;
	Real total_weight = elec_weight*w;

	// TODO
	groupelec().eval_respair_group_derivatives( rsd, rsd,
																							r1_atom_derivs, r1_atom_derivs,
																							total_weight, Erespair );

	// mutable
	Eres_[ rsd.seqpos() ] += Erespair;

	TR.Debug << "done eval residue pair deriv" << std::endl;
}

/// @details for use only with the nblist auto-update algorithm
void
FA_GrpElecEnergy::eval_atom_derivative(
	id::AtomID const &,
	pose::Pose const &,
	kinematics::DomainMap const &,// domain_map,
	ScoreFunction const &,
	EnergyMap const &,
	Vector &,
	Vector &
) const
{}

void
FA_GrpElecEnergy::finalize_total_energy(
	pose::Pose & ,
	ScoreFunction const &,
	EnergyMap & 
) const
{}

void
FA_GrpElecEnergy::evaluate_rotamer_pair_energies(
	conformation::RotamerSetBase const & set1,
	conformation::RotamerSetBase const & set2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap const & /*weights*/,
	ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
) const
{
	assert( set1.resid() != set2.resid() );

	// Since a rotamer set may include multiple residue types,
	// we'll make our decision based on what's currently in the Pose.
	if ( !defines_score_for_residue_pair( pose.residue(set1.resid()), 
																				pose.residue(set2.resid()), true) ) return;

	FAElecContextDataCOP data = utility::pointer::static_pointer_cast< FAElecContextData const >
	( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA ) );

	/*
	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( EnergiesCacheableDataType::HBOND_SET )));
	*/

	for( Size ii = 1; ii <= set1.num_rotamers(); ++ii ){
		conformation::Residue const &rot1 = *set1.rotamer(ii);
		for( Size jj = 1; jj <= set2.num_rotamers(); ++jj ){
			conformation::Residue const &rot2 = *set2.rotamer(jj);

			Real w = burial_weight( data->get_n( rot1.seqpos() ) ) 
				+ burial_weight( data->get_n( rot2.seqpos() ) );

			Real res_energy = groupelec().eval_respair_group_coulomb( rot1, rot2 );

			energy_table( jj, ii ) += w*res_energy;
		}
	}
	return;
}

void
FA_GrpElecEnergy::evaluate_rotamer_background_energies(
	conformation::RotamerSetBase const & set,
	conformation::Residue const & residue,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap const & /*weights*/,
	utility::vector1< core::PackerEnergy > & energy_vector
) const
{
	// Since a rotamer set may include multiple residue types,
	// we'll make our decision based on what's currently in the Pose.
	if ( ! defines_score_for_residue_pair( pose.residue(set.resid()), residue, true) ) return;

	FAElecContextDataCOP data = utility::pointer::static_pointer_cast< FAElecContextData const >
	( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::FAELEC_CONTEXT_DATA ) );

	/*
	hbonds::HBondSet const & hbond_set
		( static_cast< hbonds::HBondSet const & >
			( pose.energies().data().get( EnergiesCacheableDataType::HBOND_SET )));
	*/

	for( Size ii = 1; ii <= set.get_n_residue_types(); ++ii ){ // restype

		Size const ii_offset = set.get_residue_type_begin( ii );
		conformation::Residue const & ii_example_rotamer( *set.rotamer( ii_offset ));
		Vector const & ii_coord( ii_example_rotamer.nbr_atom_xyz() );
		Real const ii_radius( ii_example_rotamer.nbr_radius() );
		if ( exclude_DNA_DNA_ && ii_example_rotamer.is_DNA() && residue.is_DNA() ) continue;
		//bool intrares( rsd1.seqpos() == rsd2.seqpos() );

		Vector const & jj_coord( residue.nbr_atom_xyz() );
		Real const jj_radius( residue.nbr_radius() );

		Real dcut = ii_radius+jj_radius+5.5;

		if ( ii_coord.distance_squared( jj_coord ) < dcut*dcut ){

			Size const n = set.get_n_rotamers_for_residue_type( ii );
			for ( Size kk = 1; kk <= n; ++kk ) {
				Size const kk_rot_id = ii_offset + kk - 1;

				Size res2 = set.rotamer( kk_rot_id)->seqpos();
				Real w = burial_weight( data->get_n( residue.seqpos() ) ) 
					+ burial_weight( data->get_n( res2 ) );
				//if( intrares ) w *= intrares_scale_;

				Real const res_energy = 
					groupelec().eval_respair_group_coulomb( *set.rotamer( kk_rot_id ),
																									residue );
				energy_vector[ kk_rot_id ] += w*res_energy;
			}
		}
	}
	return;
}

/// @brief FA_GrpElecEnergy distance cutoff
 ///
 /// Reports the maximum heavy atom/heavy atom distance at which two residues have a non-zero fa_elec interaction energy.
Distance
FA_GrpElecEnergy::atomic_interaction_cutoff() const
{
	return coulomb().max_dis() + 2*core::chemical::MAX_CHEMICAL_BOND_TO_HYDROGEN_LENGTH;
}

etable::count_pair::CountPairFunctionCOP
FA_GrpElecEnergy::get_intrares_countpair(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &
) const
{
	utility_exit_with_message( "FA_GrpElecEnergy does not define intra-residue pair energies; do not call get_intrares_countpair()" );
	return 0;
}

etable::count_pair::CountPairFunctionCOP
FA_GrpElecEnergy::get_count_pair_function(
	Size const res1,
	Size const res2,
	pose::Pose const & pose,
	ScoreFunction const &
) const
{
	using namespace etable::count_pair;
	if ( res1 == res2 ) {
		return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairNone ) );
	}

	conformation::Residue const & rsd1( pose.residue( res1 ) );
	conformation::Residue const & rsd2( pose.residue( res2 ) );
	return get_count_pair_function( rsd1, rsd2 );
}

etable::count_pair::CountPairFunctionCOP
FA_GrpElecEnergy::get_count_pair_function(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
) const
{
	using namespace etable::count_pair;

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairNone ) );

	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		return CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	}
	return etable::count_pair::CountPairFunctionCOP( etable::count_pair::CountPairFunctionOP( new CountPairAll ) );

}

/// @brief FA_GrpElecEnergy
void
FA_GrpElecEnergy::indicate_required_context_graphs( utility::vector1< bool > & /* context_graphs_required */ ) const
{
}

void
FA_GrpElecEnergy::precalc_context( pose::Pose & pose,
																	FAElecContextDataOP data
																	) const
{

	TR.Debug << "precalc_context" << std::endl;

	TenANeighborGraph const & tenA_neighbor_graph( pose.energies().tenA_neighbor_graph() );

  Real cen_dist_cut2 = 10.0*10.0;

	//if( data->) TR.Debug << "?" << data->get_n( 1 ) << std::endl;

	data->initialize( pose.total_residue() );
	//TR.Debug << "??" << data->get_n( 1 ) << std::endl;

  // First get boundary neighs feeling context-dependent derivatives
	utility::vector1< Vector > dn_dr( pose.total_residue() );
	for ( Size i = 1; i <= pose.total_residue(); ++i ) dn_dr[i] = Vector( 0.0 );
	
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		conformation::Residue const & rsd1 ( pose.residue(i) );
		if ( !rsd1.is_protein() ) continue;

		for ( graph::Graph::EdgeListConstIter
						iru  = tenA_neighbor_graph.get_node(i)->const_upper_edge_list_begin(),
						irue = tenA_neighbor_graph.get_node(i)->const_upper_edge_list_end();
					iru != irue; ++iru ) {

			Size const j = (*iru)->get_second_node_ind();

			conformation::Residue const & rsd2 ( pose.residue(j) );
			if ( !rsd2.is_protein() ) continue;

			Vector cenvec = rsd2.atom( rsd2.nbr_atom() ).xyz() - rsd1.atom( rsd1.nbr_atom() ).xyz();
			Real const cendist2 = cenvec.length_squared();

			if ( cendist2 > cen_dist_cut2 ) continue;

			Real dn_drij( 0.0 );
			Real cendist = std::sqrt( cendist2 );
			Real const n_ij = eval_n( cendist, dn_drij, true );

			Real &ni = data->n( i );
			Real &nj = data->n( j );
			ni += n_ij;
			nj += n_ij;

			//TR << "Context: " << i << " " << j << " " << n_ij << " " << ni << " " << nj << std::endl;

			if( dn_drij > 0.0 ){
				data->boundary_neighs( i ).push_back( j );
				data->boundary_neighs( j ).push_back( i );
				dn_dr[ i ] += dn_drij*cenvec/cendist;
				dn_dr[ j ] -= dn_drij*cenvec/cendist;
			}
		}
	}

  for ( Size res1 = 1; res1 <= pose.total_residue(); ++res1 ){
		Real dw_dn = burial_deriv( data->get_n( res1 ) );
		data->dw_dr( res1 ) = dw_dn*dn_dr[ res1 ];
	}

	TR.Debug << "done: precalc_context" << std::endl;
}

Real
FA_GrpElecEnergy::eval_n( Real const cendist,
												 Real &dn_drij,
												 bool const eval_deriv ) const 
{
	Real interp( 0.0 );
	dn_drij = 0.0;
	if ( cendist < 9.0 ){ // < 9ang
		interp = 1.0;
	}	else if( cendist > 10.0 ){
		interp = 0;
	} else {
		interp = 0.5 + cos( M_PI*(cendist-9.0) );

		if( eval_deriv ){
			dn_drij = 0.5*std::sqrt(1.0 - interp*interp ); // just sin
		}
	}

	return interp;
}


void
FA_GrpElecEnergy::eval_context_derivatives(
	conformation::Residue const & rsd1,
	FAElecContextDataCOP data,
	EnergyMap const &,
	utility::vector1< DerivVectorPair > & r1_atom_derivs
) const
{
	TR.Debug << "eval context deriv" << std::endl;

	Vector f1( 0.0 ), f2( 0.0 );
	Size const atm1 = rsd1.nbr_atom();

	utility::vector1< Size > neighs = data->get_boundary_neighs( rsd1.seqpos() );
	for( Size jres = 1; jres <= neighs.size(); ++jres ){
		Size res2 = neighs[ jres ];
		f2 += Eres_[ res2 ]*(data->get_dw_dr( res2 ));
	}
	r1_atom_derivs[ atm1 ].f1() += f1;

	f2 = rsd1.xyz( atm1 ).cross( f1 );
	r1_atom_derivs[ atm1 ].f2() += f2;

	TR.Debug << "done eval context deriv" << std::endl;

}

// same as in hbond; 
inline core::Real
FA_GrpElecEnergy::burial_weight( core::Real const nb ) const
{
	if ( nb < 7.0 ) return 0.1;
	if ( nb > 24.0 ) return 0.5;
	return (nb-2.75)*(0.5/21.25);
}

// same as in hbond; 
inline core::Real
FA_GrpElecEnergy::burial_deriv( core::Real const nb ) const
{
	if ( nb < 7.0 ) return 0.0;
	if ( nb > 24.0 ) return 0.0;
	return 0.5/21.25;
}

core::Size
FA_GrpElecEnergy::version() const
{
	return 1; // Initial versioning
}

bool
FA_GrpElecEnergy::monomer_test(
	Size irsd,
	Size jrsd
) const {
	return (irsd <= nres_monomer_ && jrsd <= nres_monomer_ ) ||
	       (irsd >  nres_monomer_ && jrsd >  nres_monomer_ );

}

void
FA_GrpElecEnergy::set_nres_mono(
	core::pose::Pose const & pose
) const {
	for(Size i = 1; i <= pose.n_residue(); ++i) {
		if(pose.residue(i).is_upper_terminus()) {
			nres_monomer_ = i;
			//std::cerr << "nres_monomer_ " << i << std::endl;
			return;
		}
	}
}


} // namespace elec
} // namespace scoring
} // namespace core
