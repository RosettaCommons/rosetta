// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/scoring/geometric_solvation/OccludedHbondSolEnergy.cc
/// @brief  Solvation model based on penalizing potential for Hbonding to solvent
/// @author John Karanicolas




// NOTES FOR IMPROVED PERFORMANCE / POTENTIAL IMPROVEMENTS....

// The GeometricSolvation implementation was context-dependent,
// because backbone groups participating in secondary structure
// were considered exempt. Here, though, we'll compute their solvation
// energy as for any other group (partly since it's not obvious how
// else they *should* be treated). This in turn allows this energy
// term to be context-independent. The real question is....
// what should be the solvation energy for CO in secondary structure??

// Probably the best alternative would be to *NOT* compute solvation energies
// here for backbone groups in secondary structure, and instead assign
// them a fixed solvation penalty. Not clear what this penalty should be though...

// It might make sense to precompute and cache scores and derivatives
// in eg. a "geometric solvation potential" object,
// so that they don't need to be computed over and over again


// Unit Headers
#include <core/scoring/geometric_solvation/OccludedHbondSolEnergy.hh>
#include <core/scoring/geometric_solvation/OccludedHbondSolEnergyCreator.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/DerivVectorPair.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/geometric_solvation/DatabaseOccSolEne.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/prof.hh>

// Package headers

// Project headers
#include <numeric/trig.functions.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/angle_deriv.hh>

// Utility headers
#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "core.scoring.geometric_solvation.OccludedHbondSolEnergy" );

namespace core {
namespace scoring {
namespace geometric_solvation {

using namespace ObjexxFCL::format;

/// @details This must return a fresh instance of the OccludedHbondSolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
OccludedHbondSolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new geometric_solvation::OccludedHbondSolEnergy( options ) );
}

ScoreTypes
OccludedHbondSolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( occ_sol_fitted );
	return sts;
}



Vector dummy_deriv_vector_;

// jumpouts will apply if this is the best possible energy.
// this value corresponds to the discontinuity we'll deem acceptable.
// deriv_check starts to give bad results with a value of 0.05 (dist=6.6), but is mostly acceptable with 0.01 (dist=7.5)
core::Real const MIN_OCC_ENERGY = { 0.01 };


OccludedHbondSolEnergy::OccludedHbondSolEnergy(
	methods::EnergyMethodOptions const & options,
	bool const verbose )
:
	parent( methods::EnergyMethodCreatorOP( new OccludedHbondSolEnergyCreator ) ),
	occ_hbond_sol_database_( ScoringManager::get_instance()->get_DatabaseOccSolEne( options.etable_type(), MIN_OCC_ENERGY ) ),
	verbose_( verbose )
{
	if ( verbose_ ) tr <<"OccludedHbondSolEnergy constructor" << std::endl;
}

OccludedHbondSolEnergy::OccludedHbondSolEnergy( OccludedHbondSolEnergy const & src ):
	parent( src ),
	occ_hbond_sol_database_( src.occ_hbond_sol_database_ ),
	verbose_( src.verbose_ )
{
	if ( verbose_ ) tr <<"OccludedHbondSolEnergy constructor" << std::endl;
}

methods::EnergyMethodOP
OccludedHbondSolEnergy::clone() const
{
	return methods::EnergyMethodOP( new OccludedHbondSolEnergy( *this ) );
}

void
OccludedHbondSolEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}

void
OccludedHbondSolEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
 	pose.update_residue_neighbors();
}


void
OccludedHbondSolEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & ,
	ScoreFunction const &,
	EnergyMap & emap
) const
{

	if ( verbose_ ) tr << "jk evaluating residue pair energy" << std::endl;

	// Note: no count-pair stuff, these will just be computed normally
	// jk is there double-counting with other stuff, eg. backbone-dependent Dunbrack when we include self-terms like this?

	assert ( rsd1.seqpos() != rsd2.seqpos() ); // this should be computed via eval_intrares_energy

	Real occ_solE =
		res_res_occ_sol_one_way( rsd1, rsd2 ) +
		res_res_occ_sol_one_way( rsd2, rsd1 ) ;

	// store the energies
	emap[ occ_sol_fitted ] += occ_solE;

}

void
OccludedHbondSolEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & ,
	ScoreFunction const & ,
	EnergyMap & emap ) const {

	if ( verbose_ ) tr << "jk evaluating intrares energy" << std::endl;

	// behaves as not-same-residue, except that we only do the calculation once
	Real occ_solE = res_res_occ_sol_one_way( rsd, rsd );
	emap[ occ_sol_fitted ] += occ_solE;
}

/// @details return true if the two residues are moving with respect to each other.
bool
OccludedHbondSolEnergy::defines_score_for_residue_pair(
	conformation::Residue const &,
	conformation::Residue const &,
	bool res_moving_wrt_eachother
) const
{
	return res_moving_wrt_eachother;
}

void
OccludedHbondSolEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const &,
	pose::Pose const &, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	eval_residue_pair_derivatives_one_way( rsd1, rsd2, weights, r1_atom_derivs, r2_atom_derivs );
	eval_residue_pair_derivatives_one_way( rsd2, rsd1, weights, r2_atom_derivs, r1_atom_derivs );
}

void
OccludedHbondSolEnergy::eval_residue_pair_derivatives_one_way(
	conformation::Residue const & rsd1, // polar residue
	conformation::Residue const & rsd2, // occluding residue
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{

	Real energy(0); // dummy variable
	Real const atom_pair_cutoff( occ_hbond_sol_database_.atomic_interaction_cutoff() );
	Real const atom_pair_cutoff2( atom_pair_cutoff*atom_pair_cutoff );
	Real dcut = ( atom_pair_cutoff + rsd2.nbr_radius() );
	Real d2cut = dcut * dcut;
	core::Vector r2_nb_xyz( rsd2.nbr_atom_xyz());
	for ( chemical::AtomIndices::const_iterator
			hnum  = rsd1.Hpos_polar().begin(),
			hnume = rsd1.Hpos_polar().end();
			hnum != hnume; ++hnum ) {
		Size const don_h_atom( *hnum );
		Size const don_base_atom( rsd1.atom_base( don_h_atom ) );
		if ( rsd1.xyz( don_h_atom ).distance_squared( r2_nb_xyz ) > d2cut ) continue; 
		for ( Size jj = 1; jj <= rsd2.nheavyatoms(); ++jj ) {
			if ( rsd1.xyz( don_h_atom ).distance_squared( rsd2.xyz(jj) ) > atom_pair_cutoff2 ) continue;
			get_atom_atom_occ_solvation( don_h_atom, don_base_atom, rsd1, jj, rsd2, energy, true, weights[ occ_sol_fitted ],
				r1_atom_derivs[ don_base_atom ].f1(), r1_atom_derivs[ don_base_atom ].f2(),
				r1_atom_derivs[ don_h_atom ].f1(), r1_atom_derivs[ don_h_atom ].f2(),
				r2_atom_derivs[ jj ].f1(), r2_atom_derivs[ jj ].f2() );
		}
	}

	for ( chemical::AtomIndices::const_iterator
			anum = rsd1.accpt_pos().begin(),
			anume = rsd1.accpt_pos().end();
			anum != anume; ++anum ) {
		Size const acc_atom( *anum );
		Size const base_atom ( rsd1.atom_base( acc_atom ) );
		if ( rsd1.xyz( acc_atom ).distance_squared( r2_nb_xyz ) > d2cut ) continue; 
		for ( Size jj = 1; jj <= rsd2.nheavyatoms(); ++jj ) {
			if ( rsd1.xyz( acc_atom ).distance_squared( rsd2.xyz(jj) ) > atom_pair_cutoff2 ) continue;
			get_atom_atom_occ_solvation( acc_atom, base_atom, rsd1, jj, rsd2, energy, true, weights[ occ_sol_fitted ],
				r1_atom_derivs[ base_atom ].f1(), r1_atom_derivs[ base_atom ].f2(),
				r1_atom_derivs[ acc_atom ].f1(), r1_atom_derivs[ acc_atom ].f2(),
				r2_atom_derivs[ jj ].f1(), r2_atom_derivs[ jj ].f2() );
		}
	}

}


void
OccludedHbondSolEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const &,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	eval_residue_pair_derivatives_one_way( rsd, rsd, weights, atom_derivs, atom_derivs );
}



Real
OccludedHbondSolEnergy::res_res_occ_sol_one_way(
	conformation::Residue const & polar_rsd,
	conformation::Residue const & occ_rsd
) const
{

	// Rhiju importantly notes: for GeometricSolvation he originally had the code in
	// the following functions written out inside these loop -- and packing was faster.
	// Perhaps something to do with inlining or compiler optimization.
	// I've left it this way for now, because it helps prevent copying too
	// much of the code shared between residue pair scoring and for the derivatives.
	// However, if speed becomes important, here's a place to start.

	// jk note: moved the loop over occluding atoms into the next fxn, this could be the speed diff...

	Real geo_solE(0.), energy(0.);

	Real dcut = ( occ_hbond_sol_database_.atomic_interaction_cutoff() + occ_rsd.nbr_radius() );
	Real d2cut = dcut * dcut;
	core::Vector occ_nb_xyz( occ_rsd.nbr_atom_xyz());

	// cycle through donors in polar_rsd
	for ( chemical::AtomIndices::const_iterator hnum = polar_rsd.Hpos_polar().begin(), hnume = polar_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atom( *hnum );
		Size const don_base_atom( polar_rsd.atom_base( don_h_atom ) );
		if ( polar_rsd.xyz( don_h_atom ).distance_squared( occ_nb_xyz ) > d2cut ) continue; 
		for ( Size occ_atom = 1; occ_atom <= occ_rsd.natoms(); occ_atom++ ) {
			get_atom_atom_occ_solvation( don_h_atom, don_base_atom, polar_rsd, occ_atom, occ_rsd, energy );
			geo_solE += energy;
		}
	}

	// cycle through acceptors in polar_rsd
	for ( chemical::AtomIndices::const_iterator anum = polar_rsd.accpt_pos().begin(), anume = polar_rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const acc_atom( *anum );
		Size const base_atom ( polar_rsd.atom_base( acc_atom ) );
		if ( polar_rsd.xyz( acc_atom ).distance_squared( occ_nb_xyz ) > d2cut ) continue; 
		for ( Size occ_atom = 1; occ_atom <= occ_rsd.natoms(); occ_atom++ ) {
			get_atom_atom_occ_solvation( acc_atom, base_atom, polar_rsd, occ_atom, occ_rsd, energy );
			geo_solE += energy;
		}
	}

	return geo_solE;
}


void
OccludedHbondSolEnergy::get_atom_atom_occ_solvation(
	Size const polar_atom,
	Size const base_atom,
	conformation::Residue const & polar_rsd,
	Size const occ_atom,
	conformation::Residue const & occ_rsd,
	Real & energy,
	bool const update_deriv,  // = false
	Real const occ_sol_fitted_weight, // = 0.0
	//bool const update_deriv_base,  // = false
	//bool const update_deriv_occ,  // = false
	Vector & f1_base, // = dummy vector
	Vector & f2_base, // = dummy vector
	Vector & f1_polar, // = dummy vector
	Vector & f2_polar, // = dummy vector
	Vector & f1_occ, // = dummy vector
	Vector & f2_occ // = dummy vector
) const
{

	// In case of early return, initialize. Note that energy does NOT accumulate, but f1/f2 do.
	// Also note that f1 and f2 are returned unweighted.
	energy = 0.;

	// note: after testing, hydrogens need not occlude
	if ( occ_rsd.atom_is_hydrogen(occ_atom) ) return;
	//	if ( occ_atom > occ_rsd.nheavyatoms() ) return;

	// note: the lines above don't exclude Proline NV...
	// catch proline NV here (and other virtual atoms, etc.)
	if ( occ_rsd.atom_type(occ_atom).lj_radius() < 0.1 ) return;

	// can be occluded by atoms directly bonded to this group, but not by self
	if ( polar_rsd.seqpos() == occ_rsd.seqpos() ) {
		if ( polar_atom == occ_atom ) return;
		if ( base_atom == occ_atom ) return;
	}

	bool polar_atom_donates = false;
	if ( polar_rsd.atom_is_hydrogen(polar_atom) ) polar_atom_donates = true;

	if ( polar_atom_donates ) {
		// polar donor cannot be occluded by an acceptor (analogous to exact_occ_skip_Hbonders in exact model, but not quite the same)
		//for ( chemical::AtomIndices::const_iterator anum = occ_rsd.accpt_pos().begin(), anume = occ_rsd.accpt_pos().end(); anum != anume; ++anum ) {
		//	if ( occ_atom == *anum ) {
		//		return;
		//	}
		//}
		if ( occ_rsd.heavyatom_is_an_acceptor( occ_atom )) return;
	} else {
		// polar acceptor cannot be occluded by an donor base (analogous to exact_occ_skip_Hbonders in exact model, but not quite the same)
		//for ( chemical::AtomIndices::const_iterator hnum = occ_rsd.Hpos_polar().begin(), hnume = occ_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		//	Size const don_h_atom( *hnum );
		//	if ( occ_atom == occ_rsd.atom_base( don_h_atom ) ) {
		//		return;
		//	}
		//}
		if ( occ_rsd.heavyatom_has_polar_hydrogens( occ_atom )) return;
	}

	//bool update_deriv = false;
	//if ( update_deriv_polar || update_deriv_base || update_deriv_occ) update_deriv = true;
	// only one at a time
	//assert( ! ( update_deriv_polar && update_deriv_base ) );
	//assert( ! ( update_deriv_polar && update_deriv_occ ) );
	//assert( ! ( update_deriv_base && update_deriv_occ ) );

	assert( ( polar_atom_donates && atom_is_donor_h( polar_rsd, polar_atom ) ) ||
		( ( ! polar_atom_donates ) && atom_is_acceptor( polar_rsd, polar_atom ) ) );

	// If acceptor, do lookup on polar atom. If donor (ie. polar atom is a hydrogen), use the base atom instead
	Size polar_atom_type_lookup_index = polar_rsd.atom_type_index( polar_atom );
	if ( polar_atom_donates ) polar_atom_type_lookup_index = polar_rsd.atom_type_index( base_atom );

	Vector const & polar_atom_xyz( polar_rsd.atom( polar_atom ).xyz() );
	Vector const & base_atom_xyz( polar_rsd.atom( base_atom ).xyz() );

	Size const occ_atom_type_index = occ_rsd.atom_type_index( occ_atom );
	Vector const & occ_atom_xyz( occ_rsd.atom( occ_atom ).xyz() );

	// jumpout with no calculations if easy tests are violated, ie. no contribution to solvation energy
	Real const dist_sq = ( occ_atom_xyz - polar_atom_xyz).length_squared();
	if ( dist_sq > occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_max_sq_dist ) ) return;
	Real const curr_cos_angle = get_cos_angle( base_atom_xyz, polar_atom_xyz, occ_atom_xyz );
	if ( curr_cos_angle < occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_min_cos_angle ) ) return;

	// geometric filters are met, compute energy (and derivatives, if desired)
	// get the appropriate parameters
	Real const amp = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_amp );
	Real const dist_mu = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_dist_mu );
	Real const twice_dist_sigma_sq = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_twice_dist_sigma_sq );
	Real const cos_angle_mu = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_cos_angle_mu );
	Real const twice_cos_angle_sigma_sq = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_twice_cos_angle_sigma_sq );

	// Note: differences are in different order. Doesn't matter for scores, does for derivatives (and these make derivatives work).
	// Briefly, we're in the regime where dist energy contribution gets small as we get big values,
	// while cos_angle contribution gets small as we get smaller values
	Real const dist_diff = sqrt(dist_sq) - dist_mu;
	Real const cos_angle_diff = cos_angle_mu - curr_cos_angle;

	Real dist_exp, cos_angle_exp;
	if ( update_deriv ) {
		dist_exp = exp( - ( dist_diff * dist_diff / twice_dist_sigma_sq ) );
		cos_angle_exp = exp( - ( cos_angle_diff * cos_angle_diff / twice_cos_angle_sigma_sq ) );
		energy = amp * dist_exp * cos_angle_exp;
	} else {
		// do the calculation with a single exp
		energy = amp *
			exp( - ( ( dist_diff * dist_diff / twice_dist_sigma_sq ) + ( cos_angle_diff * cos_angle_diff / twice_cos_angle_sigma_sq ) ) );
	}

	if ( verbose_ && ( energy > 0. ) ) {
		tr <<"jk res "<< polar_rsd.name1() << I(3,polar_rsd.seqpos()) <<
			" atom " << polar_rsd.atom_name( polar_atom ) << " is occluded by occ_res " <<
			occ_rsd.name1() << I(3, occ_rsd.seqpos()) <<
			" atom " << occ_rsd.atom_name( occ_atom ) <<
			" with energy " << F(8,3,energy) << std::endl;
	}

	if ( ! update_deriv ) return;

	// compute angle f1/f2
	// note: energy is what was computed above, since it does NOT accumulate
	Real const curr_angle = numeric::arccos(curr_cos_angle); // note: radians
	Real const cos_angle_dfunc = -2. * cos_angle_diff * energy / twice_cos_angle_sigma_sq;
	Real const angle_dfunc = -std::sin(curr_angle) * cos_angle_dfunc * occ_sol_fitted_weight;
	Real theta(0.);
	Vector angle_f1_p1(0.), angle_f2_p1(0.);
	Vector angle_f1_p2(0.), angle_f2_p2(0.);
	Vector angle_f1_p3(0.), angle_f2_p3(0.);

	numeric::deriv::angle_p1_p2_p3_deriv( base_atom_xyz, polar_atom_xyz, occ_atom_xyz, theta, 
		angle_f1_p1, angle_f2_p1, angle_f1_p2, angle_f2_p2, angle_f1_p3, angle_f2_p3 );
	//numeric::deriv::angle_p2_deriv( base_atom_xyz, polar_atom_xyz, occ_atom_xyz, theta, angle_f1, angle_f2 );
	f1_polar += angle_dfunc * angle_f1_p2;
	f2_polar += angle_dfunc * angle_f2_p2;
	//numeric::deriv::angle_p1_deriv( base_atom_xyz, polar_atom_xyz, occ_atom_xyz, theta, angle_f1, angle_f2 );
	f1_base += angle_dfunc * angle_f1_p1;
	f2_base += angle_dfunc * angle_f2_p1;
	//numeric::deriv::angle_p1_deriv( occ_atom_xyz, polar_atom_xyz, base_atom_xyz, theta, angle_f1, angle_f2 );
	f1_occ += angle_dfunc * angle_f1_p3;
	f2_occ += angle_dfunc * angle_f2_p3;

	// compute distance f1/f2
	// note: energy is what was computed above, since it does NOT accumulate
	Real const dist_dfunc = -2. * dist_diff * energy * occ_sol_fitted_weight / twice_dist_sigma_sq;
	Real dist(0.);
	Vector dist_f1(0.), dist_f2(0.);

	numeric::deriv::distance_f1_f2_deriv( polar_atom_xyz, occ_atom_xyz, dist, dist_f1, dist_f2 );
	f1_polar += dist_dfunc * dist_f1;
	f2_polar += dist_dfunc * dist_f2;

	f1_occ -= dist_dfunc * dist_f1;
	f2_occ -= dist_dfunc * dist_f2;

	//} else {
	//	numeric::deriv::distance_f1_f2_deriv( occ_atom_xyz, polar_atom_xyz, dist, dist_f1, dist_f2 );
	//}
	//f1 += dist_dfunc * dist_f1;
	//f2 += dist_dfunc * dist_f2;

	return;

}

/*void
OccludedHbondSolEnergy::deprecated_eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const &,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
 	) const
{

	Size const curr_resnum( atom_id.rsd() );
	Size const curr_atomno( atom_id.atomno() );
	conformation::Residue const & curr_rsd( pose.residue( curr_resnum ) );
	int const curr_res_map( domain_map( curr_resnum ) );
	bool const curr_res_fixed( curr_res_map != 0 );

	// Loop over all atoms of neighboring residues, INCLUDING SELF
	utility::vector1 <core::Size> neighborlist;
	neighborlist.push_back( curr_resnum );
	EnergyGraph const & energy_graph( pose.energies().energy_graph() );
	for ( graph::Graph::EdgeListConstIter
			iru  = energy_graph.get_node( curr_resnum )->const_edge_list_begin(),
			irue = energy_graph.get_node( curr_resnum )->const_edge_list_end();
			iru != irue; ++iru ) {
		Size const other_resnum( (*iru)->get_other_ind( curr_resnum ) );
		if ( curr_res_fixed && curr_res_map == domain_map( other_resnum ) ) continue; // fixed wrt one another
		neighborlist.push_back( other_resnum );
	}

	Vector f1(0.), f2(0.);
	core::Real energy(0.); // dummy variable

	// If this atom is a polar atom, consider its occlusion by all other (neighboring) residues
	if ( atom_is_donor_h( curr_rsd, curr_atomno ) || atom_is_acceptor( curr_rsd, curr_atomno ) ) {
		Size const curr_base_atom ( curr_rsd.atom_base( curr_atomno ) );
		for ( Size other_res_inx = 1; other_res_inx <= neighborlist.size(); ++other_res_inx ) {
			Size const other_resnum( neighborlist[other_res_inx] );
			conformation::Residue const & other_rsd( pose.residue( other_resnum ) );
			for ( Size occ_atom = 1; occ_atom <= other_rsd.natoms(); occ_atom++ ) {
				get_atom_atom_occ_solvation( curr_atomno, curr_base_atom, curr_rsd, occ_atom, other_rsd, energy, true, false, false, f1, f2 );
			}
		}
	}

	// If this atom is a base atom, consider occlusion of its polar atom by all other (neighboring) residues
	if ( atom_is_valid_base( curr_rsd, curr_atomno ) ) {
		// because a base atom can have multiple polar atoms, we need to loop over them here
		for ( chemical::AtomIndices::const_iterator hnum = curr_rsd.Hpos_polar().begin(), hnume = curr_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			Size const don_h_atom( *hnum );
			Size const base_atom ( curr_rsd.atom_base( don_h_atom ) );
			if ( base_atom != curr_atomno ) continue;
			for ( Size other_res_inx = 1; other_res_inx <= neighborlist.size(); ++other_res_inx ) {
				Size const other_resnum( neighborlist[other_res_inx] );
				conformation::Residue const & other_rsd( pose.residue( other_resnum ) );
				for ( Size occ_atom = 1; occ_atom <= other_rsd.natoms(); occ_atom++ ) {
					get_atom_atom_occ_solvation( don_h_atom, curr_atomno, curr_rsd, occ_atom, other_rsd, energy, false, true, false, f1, f2 );
				}
			}
		}
		for ( chemical::AtomIndices::const_iterator anum = curr_rsd.accpt_pos().begin(), anume = curr_rsd.accpt_pos().end(); anum != anume; ++anum ) {
			Size const acc_atom( *anum );
			Size const base_atom ( curr_rsd.atom_base( acc_atom ) );
			if ( base_atom != curr_atomno ) continue;
			for ( Size other_res_inx = 1; other_res_inx <= neighborlist.size(); ++other_res_inx ) {
				Size const other_resnum( neighborlist[other_res_inx] );
				conformation::Residue const & other_rsd( pose.residue( other_resnum ) );
				for ( Size occ_atom = 1; occ_atom <= other_rsd.natoms(); occ_atom++ ) {
					get_atom_atom_occ_solvation( acc_atom, curr_atomno, curr_rsd, occ_atom, other_rsd, energy, false, true, false, f1, f2 );
				}
			}
		}
	}

	// Now consider occlusion of polar groups on neighboring residues by this atom
	for ( Size other_res_inx = 1; other_res_inx <= neighborlist.size(); ++other_res_inx ) {
		Size const other_resnum( neighborlist[other_res_inx] );
		conformation::Residue const & other_rsd( pose.residue( other_resnum ) );
		for ( chemical::AtomIndices::const_iterator hnum = other_rsd.Hpos_polar().begin(), hnume = other_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			Size const don_h_atom( *hnum );
			Size const don_base_atom( other_rsd.atom_base( don_h_atom ) );
			get_atom_atom_occ_solvation( don_h_atom, don_base_atom, other_rsd, curr_atomno, curr_rsd, energy, false, false, true, f1, f2 );
		}
		for ( chemical::AtomIndices::const_iterator anum = other_rsd.accpt_pos().begin(), anume = other_rsd.accpt_pos().end(); anum != anume; ++anum ) {
			Size const acc_atom( *anum );
			Size const base_atom ( other_rsd.atom_base( acc_atom ) );
			get_atom_atom_occ_solvation( acc_atom, base_atom, other_rsd, curr_atomno, curr_rsd, energy, false, false, true, f1, f2 );
		}
	}

	// F1/F2 accumulate
	F1 += weights[ occ_sol_fitted ] * f1;
	F2 += weights[ occ_sol_fitted ] * f2;

}*/


Distance
OccludedHbondSolEnergy::atomic_interaction_cutoff() const
{
	//	tr << "atomic_interaction_cutoff is:  " << occ_hbond_sol_database_.atomic_interaction_cutoff() << std::endl;
	// jk max interaction distance is computed using the hydrogen for donors - is this okay? or should we add one to get a heavyatom distance?
	// probably is doesn't matter, since at worst we'll just end up using an acceptor-based distance, which is fine...
	return occ_hbond_sol_database_.atomic_interaction_cutoff();
}

Real
OccludedHbondSolEnergy::get_cos_angle(
	Vector const & base_atom_xyz,
	Vector const & polar_atom_xyz,
	Vector const & occluding_atom_xyz ) const
{
	return dot( (polar_atom_xyz - base_atom_xyz).normalize(), (occluding_atom_xyz - polar_atom_xyz).normalize() );
}



// Helper function that should live inside conformation::Residue (Rhiju's comment)
bool OccludedHbondSolEnergy::atom_is_donor_h( conformation::Residue const & rsd, Size const atom ) const {
	for ( chemical::AtomIndices::const_iterator hnum = rsd.Hpos_polar().begin(), hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atom( *hnum );
		if ( don_h_atom == atom ) return true;
	}
	return false;
}

// Helper function that should live inside conformation::Residue (Rhiju's comment)
bool OccludedHbondSolEnergy::atom_is_acceptor( conformation::Residue const & rsd, Size const atom ) const {
	for ( chemical::AtomIndices::const_iterator anum = rsd.accpt_pos().begin(), anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const acc_atom( *anum );
		if ( acc_atom == atom ) return true;
	}
	return false;
}

// Helper function that should live inside conformation::Residue (Rhiju's comment)
bool OccludedHbondSolEnergy::atom_is_valid_base( conformation::Residue const & rsd, Size const atom ) const {
	for ( chemical::AtomIndices::const_iterator hnum = rsd.Hpos_polar().begin(), hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atom( *hnum );
		Size const base_atom ( rsd.atom_base( don_h_atom ) );
		if ( base_atom == atom ) return true;
	}
	for ( chemical::AtomIndices::const_iterator anum = rsd.accpt_pos().begin(), anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const acc_atom( *anum );
		Size const base_atom ( rsd.atom_base( acc_atom ) );
		if ( base_atom == atom ) return true;
	}
	return false;
}
core::Size
OccludedHbondSolEnergy::version() const
{
	return 1; // Initial versioning
}



} // geometric_solvation
} // scoring
} // core

