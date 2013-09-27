// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/GeometricSolEnergy.fwd.hh
/// @brief  Geometric solvation energy.
/// @author Rhiju Das

// Unit Headers
#include <core/scoring/geometric_solvation/GeometricSolEnergy.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyCreator.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh>

// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/types.hh>


// Project headers

// Utility headers
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/prof.hh>

#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/StringVectorOption.hh>
#include <ObjexxFCL/FArray3D.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>

static basic::Tracer tr( "core.scoring.geometric_solvation.GeometricSolEnergy" );

//////////////////////////////////////////////////////////////////
// All the good stuff is now in GeometricSolEnergyEvaluator, which is
// shared with (faster) ContextIndependentGeometricSolEnergy.
//////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace geometric_solvation {

/// @details This must return a fresh instance of the GeometricSolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
GeometricSolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new GeometricSolEnergy( options );
}

ScoreTypes
GeometricSolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( geom_sol );
	sts.push_back( geom_sol_intra_RNA );
	return sts;
}


///@brief copy c-tor
GeometricSolEnergy::GeometricSolEnergy( methods::EnergyMethodOptions const & opts
) :
	parent( new GeometricSolEnergyCreator ),
	options_( new methods::EnergyMethodOptions( opts ) ),
	evaluator_( new GeometricSolEnergyEvaluator( opts ) ),
	using_extended_method_( false )
{
}

/// copy ctor
GeometricSolEnergy::GeometricSolEnergy( GeometricSolEnergy const & src ):
	ContextDependentTwoBodyEnergy( src ),
	options_( new methods::EnergyMethodOptions( *src.options_ ) ),
	evaluator_( src.evaluator_ ),
	using_extended_method_( src.using_extended_method_ )
{}

/// clone
methods::EnergyMethodOP
GeometricSolEnergy::clone() const
{
	return new GeometricSolEnergy( *this );
}

///
void
GeometricSolEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn ) const
{
	using core::scoring::EnergiesCacheableDataType::HBOND_SET;

	// We need the H-bond set -- well, at least the backbone/backbone h-bonds
	// when computing geometric solvation scores.
	// Since this is probably being computed elsewhere, might make sense
	// to have a "calculated" flag.
	// But, anyway, the geometric sol calcs take way longer than this.
	pose.update_residue_neighbors();

	hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( options_->hbond_options() ) );
	hbond_set->setup_for_residue_pair_energies( pose );
	pose.energies().data().set( HBOND_SET, hbond_set );
	
	if ( pose.energies().use_nblist() ) {
		NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::GEOM_SOLV_NBLIST ) );
		nblist.prepare_for_scoring( pose, scfxn, *this );
	}
}

/////////////////////////////////////////////////////////////////////////////
void
GeometricSolEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const
{
	using core::scoring::EnergiesCacheableDataType::HBOND_SET;

// same setup as for HBondEnergy.cc. Note that this is probably repeating some work
// that has already occurred in hbonds calculation. Probably could have a "calculated"
// flag to save the computation ... but the geometric sol. calculation takes so
// much longer that this initial hbond loop isn't too bad.
 	pose.update_residue_neighbors();

	hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( options_->hbond_options() ) );
	hbonds::fill_hbond_set( pose, true /*calc derivs*/, *hbond_set );
	hbond_set->resize_bb_donor_acceptor_arrays( pose.total_residue() );
	pose.energies().data().set( HBOND_SET, hbond_set );
}


//void
//GeometricSolEnergy::setup_for_derivatives( pose::Pose & /*pose*/, ScoreFunction const & ) const
//{
//}

///////////////////////////////////////////////////////////////////////////////
void
GeometricSolEnergy::setup_for_minimizing(
    pose::Pose & pose,
    ScoreFunction const & sfxn,
    kinematics::MinimizerMapBase const & min_map
) const
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    
    //set_nres_mono(pose);
    
    if ( pose.energies().use_nblist() ) {
	//if ( true ) {
        // stash our nblist inside the pose's energies object
        Energies & energies( pose.energies() );
        
        // setup the atom-atom nblist
        NeighborListOP nblist;
        Real const tolerated_motion = pose.energies().use_nblist_auto_update() ? option[ run::nblist_autoupdate_narrow ] : 1.5;
        Real const XX = 5.2 + 2 * tolerated_motion;
        nblist = new NeighborList( min_map.domain_map(), XX*XX, XX*XX, XX*XX);
        if ( pose.energies().use_nblist_auto_update() ) {
            nblist->set_auto_update( tolerated_motion );
        }
        // this partially becomes the EtableEnergy classes's responsibility
        nblist->setup( pose, sfxn, *this);
        energies.set_nblist( EnergiesCacheableDataType::GEOM_SOLV_NBLIST, nblist );
    }
}

///////////////////////////////////////////////////////////////////////////////
bool
GeometricSolEnergy::defines_score_for_residue_pair(
    conformation::Residue const & rsd1,
    conformation::Residue const & rsd2,
    bool res_moving_wrt_eachother
) const
{
    if ( rsd1.seqpos() == rsd2.seqpos() ) {
        return false;
    }
    return res_moving_wrt_eachother;
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
GeometricSolEnergy::get_count_pair_function(
    Size const res1,
    Size const res2,
    pose::Pose const & pose,
    ScoreFunction const &
) const
{
    using namespace etable::count_pair;
    if ( res1 == res2 ) {
        return new CountPairNone;
    }
    
    conformation::Residue const & rsd1( pose.residue( res1 ) );
    conformation::Residue const & rsd2( pose.residue( res2 ) );
    return get_count_pair_function( rsd1, rsd2 );
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
GeometricSolEnergy::get_count_pair_function(
    conformation::Residue const & rsd1,
    conformation::Residue const & rsd2
) const
{
    using namespace etable::count_pair;
    
    if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return new CountPairNone;
    
    if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
        return CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
    }
    return new CountPairAll;
    
}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
GeometricSolEnergy::get_intrares_countpair(
    conformation::Residue const & res,
    pose::Pose const &,
    ScoreFunction const &
) const
{
	using namespace etable::count_pair;
	return CountPairFactory::create_intrares_count_pair_function( res, CP_CROSSOVER_3 );
}

/////////////////////////////
bool
GeometricSolEnergy::minimize_in_whole_structure_context( pose::Pose const & pose ) const
{
	return pose.energies().use_nblist_auto_update();
}


///////////////////////////////////////////////////////////////////////////////
bool
GeometricSolEnergy::use_extended_residue_pair_energy_interface() const
{
    return true;
}

////////////////////////////////////////////////////////////////////////////////
void
GeometricSolEnergy::setup_for_minimizing_for_residue_pair(
    conformation::Residue const & rsd1,
    conformation::Residue const & rsd2,
    pose::Pose const & pose,
    ScoreFunction const &,
    kinematics::MinimizerMapBase const &,
    ResSingleMinimizationData const &,
    ResSingleMinimizationData const &,
    ResPairMinimizationData & pair_data
) const
{
    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    if ( pose.energies().use_nblist_auto_update() ) return;
    
    etable::count_pair::CountPairFunctionCOP count_pair = get_count_pair_function( rsd1, rsd2 );
    //assert( rsd1.seqpos() < rsd2.seqpos() );
    
    // update the existing nblist if it's already present in the min_data object
    ResiduePairNeighborListOP nblist( static_cast< ResiduePairNeighborList * > (pair_data.get_data( geom_solv_pair_nblist )() ));
    if ( ! nblist ) nblist = new ResiduePairNeighborList;
    
    /// STOLEN CODE!
    Real const tolerated_narrow_nblist_motion = 0.75; //option[ run::nblist_autoupdate_narrow ];
    Real const XX2 = std::pow( 5.2 + 2*tolerated_narrow_nblist_motion, 2 );
    
    nblist->initialize_from_residues( XX2, XX2, XX2, rsd1, rsd2, count_pair );
    
    pair_data.set_data( geom_solv_pair_nblist, nblist );
}

////////////////////////////////////////////////////////////////////////////////
void
GeometricSolEnergy::residue_pair_energy_ext(
    conformation::Residue const & rsd1,
    conformation::Residue const & rsd2,
    ResPairMinimizationData const & min_data,
    pose::Pose const & pose,
    ScoreFunction const &,
    EnergyMap & emap
) const
{
    using_extended_method_ = true;
	//return;
	if ( pose.energies().use_nblist_auto_update() ) return;
    Real score( 0.0 );
    Real energy( 0.0 );
    
    //assert( rsd1.seqpos() < rsd2.seqpos() );
    
    //assert( dynamic_cast< ResiduePairNeighborList const * > (min_data.get_data( geom_solv_pair_nblist )() ));
    ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( geom_solv_pair_nblist ) ) );
    utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );
    Size m = 0;
    Size n = 0;
        
    for ( Size ii = 1, iiend = neighbs.size(); ii <= iiend; ++ii ) {
        m = neighbs[ ii ].atomno1();
        n = neighbs[ ii ].atomno2();
        if ( evaluator_->atom_is_heavy( rsd2, n ) ) {
            if ( evaluator_->atom_is_donor_h( rsd1, m ) ) {
                evaluator_->get_atom_atom_geometric_solvation_for_donor( m, rsd1, n, rsd2, pose, energy );
                score += energy;
            } else if ( evaluator_->atom_is_acceptor( rsd1, m ) ) {
                evaluator_->get_atom_atom_geometric_solvation_for_acceptor( m, rsd1, n, rsd2, pose, energy );
                score += energy;
            }
        }
        if ( evaluator_->atom_is_heavy ( rsd1, m ) ) {
            if ( evaluator_->atom_is_donor_h( rsd2, n ) ) {
                evaluator_->get_atom_atom_geometric_solvation_for_donor( n, rsd2, m, rsd1, pose, energy );
                score += energy;
            } else if ( evaluator_->atom_is_acceptor( rsd2, n ) ) {
                evaluator_->get_atom_atom_geometric_solvation_for_acceptor( n, rsd2, m, rsd1, pose, energy );
                score += energy;
            }
        }
    }
    emap[ geom_sol ] += score;
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

/// Everything in here.
void
GeometricSolEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap ) const
{
	using_extended_method_ = false;
	if ( pose.energies().use_nblist() ) return;
    evaluator_->residue_pair_energy( rsd1, rsd2, pose, scorefxn, emap );
}
    
void
GeometricSolEnergy::finalize_total_energy(
    pose::Pose & pose,
    ScoreFunction const &,
    EnergyMap & totals ) const
{
    //if ( !using_extended_method_ ) return;
    if ( ! pose.energies().use_nblist() || ! pose.energies().use_nblist_auto_update() ) return;
    NeighborList const & nblist
        ( pose.energies().nblist( EnergiesCacheableDataType::GEOM_SOLV_NBLIST ) );
    nblist.check_domain_map( pose.energies().domain_map() );
    utility::vector1< conformation::Residue const * > resvect;
    resvect.reserve( pose.total_residue() );
    for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
        resvect.push_back( & pose.residue( ii ) );
    }
    Real score( 0.0 );
    Real energy( 0.0 );
    
    for ( Size i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
        conformation::Residue const & ires( *resvect[i] );
        for ( Size ii=1, ii_end=ires.natoms(); ii<= ii_end; ++ii ) {
            AtomNeighbors const & nbrs( nblist.upper_atom_neighbors(i,ii) );
            for ( AtomNeighbors::const_iterator nbr_iter=nbrs.begin(),
                    nbr_end=nbrs.end(); nbr_iter!= nbr_end; ++nbr_iter ) {
                AtomNeighbor const & nbr( *nbr_iter );
                
                Size const  j( nbr.rsd() );
                if ( i==j ) continue;
                Size const jj( nbr.atomno() );
                // could reorder the nbr lists so that we dont need this check:
                //if ( ( j < i ) || ( j == i && jj <= ii ) ) continue;
                
                conformation::Residue const & jres( *resvect[j] );
                
                if ( evaluator_->atom_is_heavy( jres, jj ) ) {
                    if ( evaluator_->atom_is_donor_h( ires, ii ) ) {
                        evaluator_->get_atom_atom_geometric_solvation_for_donor( ii, ires, jj, jres, pose, energy );
                        score += energy;
                    } else if ( evaluator_->atom_is_acceptor( ires, ii ) ) {
                        evaluator_->get_atom_atom_geometric_solvation_for_acceptor( ii, ires, jj, jres, pose, energy );
                        score += energy;
                    }
                }
                if ( evaluator_->atom_is_heavy ( ires, ii ) ) {
                    if ( evaluator_->atom_is_donor_h( jres, jj ) ) {
                        evaluator_->get_atom_atom_geometric_solvation_for_donor( jj, jres, ii, ires, pose, energy );
                        score += energy;
                    } else if ( evaluator_->atom_is_acceptor( jres, jj ) ) {
                        evaluator_->get_atom_atom_geometric_solvation_for_acceptor( jj, jres, ii, ires, pose, energy );
                        score += energy;
                    }
                }
            }
        }
    }
    
    totals[ geom_sol ] += score;
}

void
GeometricSolEnergy::eval_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	ScoreFunction const & scorefxn,
	EnergyMap const & weights,
	Vector & F1,
	Vector & F2
) const
{
	if ( !use_extended_residue_pair_energy_interface() ) {
        evaluator_->eval_atom_derivative( atom_id, pose, domain_map, scorefxn, weights, F1, F2 );
        return;
    }
	
	if ( ! pose.energies().use_nblist_auto_update() ) return;
    
    Real energy( 0.0 );
	hbonds::HBondDerivs deriv;
    
    if ( defines_intrares_energy( weights ) ) {
        evaluator_->eval_atom_derivative_intra_RNA(atom_id, pose, weights, F1, F2);
    }
    
    Size const i( atom_id.rsd() );
	conformation::Residue const & ires( pose.residue( i ) );
	Size const ii( atom_id.atomno() );
    
	//	Size const nres = pose.total_residue();
	static bool const update_deriv( true );
    //assert( pose.energies().use_nblist() );
	NeighborList const & nblist
        ( pose.energies().nblist( EnergiesCacheableDataType::GEOM_SOLV_NBLIST ) );
	AtomNeighbors const & nbrs( nblist.atom_neighbors(i,ii) );
    
    for ( scoring::AtomNeighbors::const_iterator it2=nbrs.begin(),
         it2e=nbrs.end(); it2 != it2e; ++it2 ) {
		scoring::AtomNeighbor const & nbr( *it2 );
		Size const j( nbr.rsd() );
        if (i == j) continue;
		Size const jj( nbr.atomno() );
		conformation::Residue const & jres( pose.residue( j ) );
        
        if ( evaluator_->atom_is_heavy( jres, jj ) ) {
            if ( evaluator_->atom_is_donor_h( ires, ii ) ) {
                evaluator_->get_atom_atom_geometric_solvation_for_donor( ii, ires, jj, jres, pose, energy, update_deriv, deriv );
                F1 += weights[ geom_sol ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
                F2 += weights[ geom_sol ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
            } else if ( evaluator_->atom_is_acceptor( ires, ii ) ) {
                evaluator_->get_atom_atom_geometric_solvation_for_acceptor( ii, ires, jj, jres, pose, energy, update_deriv, deriv );
                F1 += weights[ geom_sol ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
                F2 += weights[ geom_sol ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
            }
        }
        if ( evaluator_->atom_is_heavy ( ires, ii ) ) {
            if ( evaluator_->atom_is_donor_h( jres, jj ) ) {
                evaluator_->get_atom_atom_geometric_solvation_for_donor( jj, jres, ii, ires, pose, energy, update_deriv, deriv );
                F1 -= weights[ geom_sol ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
                F2 -= weights[ geom_sol ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
            } else if ( evaluator_->atom_is_acceptor( jres, jj ) ) {
                evaluator_->get_atom_atom_geometric_solvation_for_acceptor( jj, jres, ii, ires, pose, energy, update_deriv, deriv );
                F1 -= weights[ geom_sol ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
                F2 -= weights[ geom_sol ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
            }
        }
    }

}

Distance
GeometricSolEnergy::atomic_interaction_cutoff() const
{
	return evaluator_->atomic_interaction_cutoff();
}

///////////////////////////////////////////////////////////////////////////////////////
bool
GeometricSolEnergy::defines_intrares_energy( EnergyMap const & weights ) const
{
	//Change to this on Feb 06, 2012. Ensure that the function returns false if weights[geom_sol_intra_RNA] == 0.0
	bool condition_1 = (weights[geom_sol_intra_RNA] > 0.0001) ? true : false;
	return condition_1;
}

///////////////////////////////////////////////////////////////////////////////////////
void
GeometricSolEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap ) const
{
	evaluator_->eval_intrares_energy( rsd, pose, scorefxn, emap );
}

//////////////////////////////
	////////////////////////////////////////////////////
void
GeometricSolEnergy::eval_residue_pair_derivatives(
	conformation::Residue const & ires,
	conformation::Residue const & jres,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData const & min_data,
	pose::Pose const & pose, // provides context
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & r1_atom_derivs,
	utility::vector1< DerivVectorPair > & r2_atom_derivs
) const
{
	if ( pose.energies().use_nblist_auto_update() ) return;
	
	//assert( rsd1.seqpos() < rsd2.seqpos() );
	//assert( dynamic_cast< ResiduePairNeighborList const * > (min_data.get_data( elec_pair_nblist )() ));
	
	ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > ( min_data.get_data_ref( geom_solv_pair_nblist ) ) );
	utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );
	
	Real energy( 0.0 );
	Vector F1;
	Vector F2;
	hbonds::HBondDerivs deriv;
	Size ii;
	Size jj;
	static bool const update_deriv( true );
	
	for ( Size k = 1, kend = neighbs.size(); k <= kend; ++k ) {
		ii = neighbs[ k ].atomno1();
		jj = neighbs[ k ].atomno2();
		F1 = 0.0;
		F2 = 0.0;
		
		if ( evaluator_->atom_is_heavy( jres, jj ) ) {
			if ( evaluator_->atom_is_donor_h( ires, ii ) ) {
				evaluator_->get_atom_atom_geometric_solvation_for_donor( ii, ires, jj, jres, pose, energy, update_deriv, deriv );
				F1 += weights[ geom_sol_fast ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
				F2 += weights[ geom_sol_fast ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
			} else if ( evaluator_->atom_is_acceptor( ires, ii ) ) {
				evaluator_->get_atom_atom_geometric_solvation_for_acceptor( ii, ires, jj, jres, pose, energy, update_deriv, deriv );
				F1 += weights[ geom_sol_fast ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
				F2 += weights[ geom_sol_fast ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
			}
		}
		if ( evaluator_->atom_is_heavy ( ires, ii ) ) {
			if ( evaluator_->atom_is_donor_h( jres, jj ) ) {
				evaluator_->get_atom_atom_geometric_solvation_for_donor( jj, jres, ii, ires, pose, energy, update_deriv, deriv );
				F1 -= weights[ geom_sol_fast ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
				F2 -= weights[ geom_sol_fast ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
			} else if ( evaluator_->atom_is_acceptor( jres, jj ) ) {
				evaluator_->get_atom_atom_geometric_solvation_for_acceptor( jj, jres, ii, ires, pose, energy, update_deriv, deriv );
				F1 -= weights[ geom_sol_fast ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
				F2 -= weights[ geom_sol_fast ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
			}
		}
		
		r1_atom_derivs[ ii ].f1() += F1;
		r1_atom_derivs[ ii ].f2() += F2;
		r2_atom_derivs[ jj ].f1() -= F1;
		r2_atom_derivs[ jj ].f2() -= F2;
	}
}

///////////////////////////////////////////////////////////////////////////////////////
///@brief GeometricSolEnergy is context sensitive
void
GeometricSolEnergy::indicate_required_context_graphs(
	utility::vector1< bool > & context_graphs_required
) const
{
	context_graphs_required[ ten_A_neighbor_graph ] = true;
}

core::Size
GeometricSolEnergy::version() const
{
	return 2; // Initial versioning
}


} // hbonds
} // scoring
} // core

