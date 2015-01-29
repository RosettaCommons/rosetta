// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/geometric_solvation/OccludedHbondSolEnergy_onebody.cc
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
#include <core/scoring/geometric_solvation/OccludedHbondSolEnergy_onebody.hh>
#include <core/scoring/geometric_solvation/OccludedHbondSolEnergy_onebodyCreator.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/geometric_solvation/DatabaseOccSolEne.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/prof.hh>

// Package headers

// Project headers
#include <numeric/trig.functions.hh>
// AUTO-REMOVED #include <numeric/deriv/distance_deriv.hh>
// AUTO-REMOVED #include <numeric/deriv/angle_deriv.hh>

// Utility headers
// AUTO-REMOVED #include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


static thread_local basic::Tracer tr( "core.scoring.geometric_solvation.OccludedHbondSolEnergy_onebody" );

namespace core {
namespace scoring {
namespace geometric_solvation {

using namespace ObjexxFCL::format;

/// @details This must return a fresh instance of the OccludedHbondSolEnergy_onebody class,
/// never an instance already in use
methods::EnergyMethodOP
OccludedHbondSolEnergy_onebodyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return methods::EnergyMethodOP( new geometric_solvation::OccludedHbondSolEnergy_onebody( options ) );
}

ScoreTypes
OccludedHbondSolEnergy_onebodyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( occ_sol_fitted_onebody );
	return sts;
}


// jumpouts will apply if this is the best possible energy.
// this value corresponds to the discontinuity we'll deem acceptable.
// deriv_check starts to give bad results with a value of 0.05 (dist=6.6), but is mostly acceptable with 0.01 (dist=7.5)
core::Real const MIN_OCC_ENERGY = { 0.01 };


OccludedHbondSolEnergy_onebody::OccludedHbondSolEnergy_onebody(
	methods::EnergyMethodOptions const & options,
	bool const verbose )
:
	parent( methods::EnergyMethodCreatorOP( new OccludedHbondSolEnergy_onebodyCreator ) ),
	occ_hbond_sol_database_( ScoringManager::get_instance()->get_DatabaseOccSolEne( options.etable_type(), MIN_OCC_ENERGY ) ),
	verbose_( verbose )
{
	if ( verbose_ ) tr <<"OccludedHbondSolEnergy_onebody constructor" << std::endl;
}

OccludedHbondSolEnergy_onebody::OccludedHbondSolEnergy_onebody( OccludedHbondSolEnergy_onebody const & src ):
	parent( src ),
	occ_hbond_sol_database_( src.occ_hbond_sol_database_ ),
	verbose_( src.verbose_ )
{
	if ( verbose_ ) tr <<"OccludedHbondSolEnergy_onebody constructor" << std::endl;
}

methods::EnergyMethodOP
OccludedHbondSolEnergy_onebody::clone() const
{
	return methods::EnergyMethodOP( new OccludedHbondSolEnergy_onebody( *this ) );
}

void
OccludedHbondSolEnergy_onebody::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}

void
OccludedHbondSolEnergy_onebody::setup_for_packing(
	pose::Pose & pose,
	utility::vector1< bool > const &,
	utility::vector1< bool > const &
) const
{
	pose.update_residue_neighbors();
}

void
OccludedHbondSolEnergy_onebody::setup_for_derivatives( pose::Pose & , ScoreFunction const & ) const
{
	tr << "Error - no derivatives yet for OccludedHbondSolEnergy_onebody (occ_sol_fitted_onebody)" << std::endl;
debug_assert(false);
	exit(1);
}

void
OccludedHbondSolEnergy_onebody::setup_for_minimizing( pose::Pose & , ScoreFunction const & , kinematics::MinimizerMapBase const & ) const
{
	tr << "Error - no derivatives yet for OccludedHbondSolEnergy_onebody (occ_sol_fitted_onebody)" << std::endl;
debug_assert(false);
	exit(1);
}

Distance
OccludedHbondSolEnergy_onebody::atomic_interaction_cutoff() const
{
	tr << "atomic_interaction_cutoff is:  " << occ_hbond_sol_database_.atomic_interaction_cutoff() << std::endl;
	// jk max interaction distance is computed using the hydrogen for donors - is this okay? or should we add one to get a heavyatom distance?
	// probably is doesn't matter, since at worst we'll just end up using an acceptor-based distance, which is fine...
	return occ_hbond_sol_database_.atomic_interaction_cutoff();
}


void OccludedHbondSolEnergy_onebody::residue_energy(
	conformation::Residue const & polar_rsd,
	pose::Pose const & pose,
	EnergyMap & emap
) const {

	core::Size polar_resnum = (core::Size) polar_rsd.seqpos();
	core::Real residue_geosol(0.), energy(0.);

	// loop over all atoms of neighboring residues, INCLUDING SELF
	core::scoring::TenANeighborGraph const & graph = pose.energies().tenA_neighbor_graph();
	utility::vector1 <core::Size> neighborlist;
	neighborlist.push_back( polar_resnum);
	for ( core::graph::Graph::EdgeListConstIter
					neighbor_iter = graph.get_node( polar_resnum )->const_edge_list_begin(),
					neighbor_iter_end = graph.get_node( polar_resnum )->const_edge_list_end();
				neighbor_iter != neighbor_iter_end; ++neighbor_iter ) {
		neighborlist.push_back( (*neighbor_iter)->get_other_ind( polar_resnum ) );
	}

	// cycle through donors in polar_rsd
	for ( chemical::AtomIndices::const_iterator hnum = polar_rsd.Hpos_polar().begin(), hnume = polar_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atom( *hnum );
		Size const base_atom( polar_rsd.atom_base( don_h_atom ) );
		core::Real polar_group_energy = 0.;
		for ( Size occ_inx = 1; occ_inx <= neighborlist.size(); ++occ_inx ) {
			core::Size const occ_resnum( neighborlist[occ_inx] );
			conformation::Residue const occ_rsd = pose.residue(occ_resnum);
			for ( Size occ_atom = 1; occ_atom <= occ_rsd.natoms(); occ_atom++ ) {
				get_atom_atom_occ_solvation( don_h_atom, base_atom, polar_rsd, occ_atom, occ_rsd, energy );
				polar_group_energy += energy;
			}
		}
		residue_geosol += polar_group_energy;
		std::string const base_atom_name = polar_rsd.atom_name( base_atom );
		//		std::cout << "jk FITTED_ONEBODY Donor " << base_atom_name << "  " << pose.residue(polar_resnum).aa() << " " << polar_resnum << "  " << polar_group_energy << std::endl;
	}

	// cycle through acceptors in polar_rsd
	for ( chemical::AtomIndices::const_iterator anum = polar_rsd.accpt_pos().begin(), anume = polar_rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const acc_atom( *anum );
		Size const base_atom ( polar_rsd.atom_base( acc_atom ) );
		core::Real polar_group_energy = 0.;
		for ( Size occ_inx = 1; occ_inx <= neighborlist.size(); ++occ_inx ) {
			core::Size const occ_resnum( neighborlist[occ_inx] );
			conformation::Residue const occ_rsd = pose.residue(occ_resnum);
			for ( Size occ_atom = 1; occ_atom <= occ_rsd.natoms(); occ_atom++ ) {
				get_atom_atom_occ_solvation( acc_atom, base_atom, polar_rsd, occ_atom, occ_rsd, energy );
				polar_group_energy += energy;
			}
		}
		residue_geosol += polar_group_energy;
		std::string const base_atom_name = polar_rsd.atom_name( base_atom );
		//		std::cout << "jk FITTED_ONEBODY Acceptor " << base_atom_name << "  " << pose.residue(polar_resnum).aa() << " " << polar_resnum << "  " << polar_group_energy << std::endl;
	}

	emap[ occ_sol_fitted_onebody ] += residue_geosol;

}



Real
OccludedHbondSolEnergy_onebody::res_res_occ_sol_one_way(
	conformation::Residue const & polar_rsd,
	conformation::Residue const & occ_rsd ) const
{

	// Rhiju importantly notes: for GeometricSolvation he originally had the code in
	// the following functions written out inside these loop -- and packing was faster.
	// Perhaps something to do with inlining or compiler optimization.
	// I've left it this way for now, because it helps prevent copying too
	// much of the code shared between residue pair scoring and for the derivatives.
	// However, if speed becomes important, here's a place to start.

	// jk note: moved the loop over occluding atoms into the next fxn, this could be the speed diff...

	Real geo_solE(0.), energy(0.);

	// cycle through donors in polar_rsd
	for ( chemical::AtomIndices::const_iterator hnum = polar_rsd.Hpos_polar().begin(), hnume = polar_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atom( *hnum );
		Size const don_base_atom( polar_rsd.atom_base( don_h_atom ) );
		for ( Size occ_atom = 1; occ_atom <= occ_rsd.natoms(); occ_atom++ ) {
			get_atom_atom_occ_solvation( don_h_atom, don_base_atom, polar_rsd, occ_atom, occ_rsd, energy );
			geo_solE += energy;
		}
	}

	// cycle through acceptors in polar_rsd
	for ( chemical::AtomIndices::const_iterator anum = polar_rsd.accpt_pos().begin(), anume = polar_rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const acc_atom( *anum );
		Size const base_atom ( polar_rsd.atom_base( acc_atom ) );
		for ( Size occ_atom = 1; occ_atom <= occ_rsd.natoms(); occ_atom++ ) {
			get_atom_atom_occ_solvation( acc_atom, base_atom, polar_rsd, occ_atom, occ_rsd, energy );
			geo_solE += energy;
		}
	}

	return geo_solE;
}


void
OccludedHbondSolEnergy_onebody::get_atom_atom_occ_solvation(
	Size const polar_atom,
	Size const base_atom,
	conformation::Residue const & polar_rsd,
	Size const occ_atom,
	conformation::Residue const & occ_rsd,
	Real & energy
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
		for ( chemical::AtomIndices::const_iterator anum = occ_rsd.accpt_pos().begin(), anume = occ_rsd.accpt_pos().end(); anum != anume; ++anum ) {
			if ( occ_atom == *anum ) {
				return;
			}
		}
	} else {
		// polar acceptor cannot be occluded by an donor base (analogous to exact_occ_skip_Hbonders in exact model, but not quite the same)
		for ( chemical::AtomIndices::const_iterator hnum = occ_rsd.Hpos_polar().begin(), hnume = occ_rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
			Size const don_h_atom( *hnum );
			if ( occ_atom == occ_rsd.atom_base( don_h_atom ) ) {
				return;
			}
		}
	}

debug_assert( ( polar_atom_donates && atom_is_donor_h( polar_rsd, polar_atom ) ) ||
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

	// geometric filters are met, compute energy (no derivatives!)
	// get the appropriate parameters
	Real const amp = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_amp );
	Real const dist_mu = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_dist_mu );
	Real const twice_dist_sigma_sq = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_twice_dist_sigma_sq );
	Real const cos_angle_mu = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_cos_angle_mu );
	Real const twice_cos_angle_sigma_sq = occ_hbond_sol_database_( polar_atom_donates, polar_atom_type_lookup_index, occ_atom_type_index, OccFitParam_twice_cos_angle_sigma_sq );

	// Note: differences are in different order. Doesn't matter for scores, does for derivatives
	// Briefly, we're in the regime where dist energy contribution gets small as we get big values,
	// while cos_angle contribution gets small as we get smaller values
	Real const dist_diff = sqrt(dist_sq) - dist_mu;
	Real const cos_angle_diff = cos_angle_mu - curr_cos_angle;
	energy = amp *
		exp( - ( ( dist_diff * dist_diff / twice_dist_sigma_sq ) + ( cos_angle_diff * cos_angle_diff / twice_cos_angle_sigma_sq ) ) );

	return;

}


Real
OccludedHbondSolEnergy_onebody::get_cos_angle( Vector const & base_atom_xyz,
																	 Vector const & polar_atom_xyz,
																	 Vector const & occluding_atom_xyz ) const
{
	return dot( (polar_atom_xyz - base_atom_xyz).normalize(), (occluding_atom_xyz - polar_atom_xyz).normalize() );
}



// Helper function that should live inside conformation::Residue (Rhiju's comment)
bool OccludedHbondSolEnergy_onebody::atom_is_donor_h( conformation::Residue const & rsd, Size const atom ) const {
	for ( chemical::AtomIndices::const_iterator hnum = rsd.Hpos_polar().begin(), hnume = rsd.Hpos_polar().end(); hnum != hnume; ++hnum ) {
		Size const don_h_atom( *hnum );
		if ( don_h_atom == atom ) return true;
	}
	return false;
}

// Helper function that should live inside conformation::Residue (Rhiju's comment)
bool OccludedHbondSolEnergy_onebody::atom_is_acceptor( conformation::Residue const & rsd, Size const atom ) const {
	for ( chemical::AtomIndices::const_iterator anum = rsd.accpt_pos().begin(), anume = rsd.accpt_pos().end(); anum != anume; ++anum ) {
		Size const acc_atom( *anum );
		if ( acc_atom == atom ) return true;
	}
	return false;
}

// Helper function that should live inside conformation::Residue (Rhiju's comment)
bool OccludedHbondSolEnergy_onebody::atom_is_valid_base( conformation::Residue const & rsd, Size const atom ) const {
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
OccludedHbondSolEnergy_onebody::version() const
{
	return 1; // Initial versioning
}



} // geometric_solvation
} // scoring
} // core

