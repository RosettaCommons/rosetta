// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/ContextIndependentGeometricSolEnergy.cc
/// @author Parin Sripakdeevong (sripakpa@stanford.edu), Rhiju Das (rhiju@stanford.edu)
/// @brief  Similar to the standard version of GeometricSolEnergy.cc BUT without the CONTEXT_DEPENDENT stuff. ALOT OF CODE DUPLICATION!
/// @brief  Significantly speed up when used in src/protocol/swa/rna/StepwiseRNA_Sampler.cc!

// Unit Headers
#include <core/scoring/geometric_solvation/ContextIndependentGeometricSolEnergy.hh>
#include <core/scoring/geometric_solvation/ContextIndependentGeometricSolEnergyCreator.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.hh>


// Package headers
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/kinematics/MinimizerMapBase.hh>

#include <core/id/AtomID.hh>

// Project headers
#include <core/conformation/Residue.hh>

// Utility headers
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>


static basic::Tracer tr("core.scoring.geometric_solvation.ContextIndependentGeometricSolEnergy" );

//////////////////////////////////////////////////
//////////////////////////////////////////////////
// Added on July. 22, 2011, Parin Sripakdeevong (sripakpa@stanford.edu).
//
// This copies huge amounts of code from GeometricSolEnergy.cc.
// Should instead make a GeometricSolPotential.cc, which holds *all* the core functions,
// and then GeometricSolEnergy and ContextIndependentGeometricSolEnergy can both call those core functions.
///////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace geometric_solvation {


/// @details This must return a fresh instance of the GeometricSolEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
ContextIndependentGeometricSolEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & options
) const {
	return new ContextIndependentGeometricSolEnergy( options );
}

ScoreTypes
ContextIndependentGeometricSolEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( geom_sol_fast );
	sts.push_back( geom_sol_fast_intra_RNA );


	return sts;
}

///@brief copy c-tor
ContextIndependentGeometricSolEnergy::ContextIndependentGeometricSolEnergy( methods::EnergyMethodOptions const & opts) :
	parent( new ContextIndependentGeometricSolEnergyCreator ),
	options_( new methods::EnergyMethodOptions( opts ) ),
	evaluator_( new GeometricSolEnergyEvaluator( opts ) ),
  precalculated_bb_bb_energy_(0.0f),
	using_extended_method_(false)
{
	if ( options_->hbond_options().use_hb_env_dep() ) {
		utility_exit_with_message( "Environment dependent hydrogen bonds are not compatible with geom_sol_fast. You can turn off hydrogen-bond environment dependence (e.g., with NO_HB_ENV_DEP in score file), and protein & RNA structure prediction won't get worse! (Or if you insist on using environment dependence, use geom_sol instead of geom_sol_fast.)" );
	}
}

/// copy ctor
ContextIndependentGeometricSolEnergy::ContextIndependentGeometricSolEnergy( ContextIndependentGeometricSolEnergy const & src ):
	ContextIndependentTwoBodyEnergy( src ),
	options_( new methods::EnergyMethodOptions( *src.options_ ) ),
	evaluator_( src.evaluator_ ),
    precalculated_bb_bb_energy_(src.precalculated_bb_bb_energy_),
	using_extended_method_(src.using_extended_method_)
{
}

/// clone
methods::EnergyMethodOP
ContextIndependentGeometricSolEnergy::clone() const
{
	return new ContextIndependentGeometricSolEnergy( *this );
}


void
ContextIndependentGeometricSolEnergy::precalculate_bb_bb_energy_for_design(
 pose::Pose const & pose
) const {

  Size const total_residue = pose.total_residue();

  EnergyGraph const & energy_graph( pose.energies().energy_graph() );

  for (Size i = 1; i <= total_residue; i++ ){

    conformation::Residue const & res_i( pose.residue( i ) );

    for( graph::Graph::EdgeListConstIter
        iter = energy_graph.get_node( i )->const_edge_list_begin();
        iter != energy_graph.get_node( i )->const_edge_list_end();
        ++iter ){

      Size j( (*iter)->get_other_ind( i ) );


      conformation::Residue const & res_j( pose.residue( j ) );

      //only need to do it one way since will sample the reverse when res_i is at j
      precalculated_bb_bb_energy_ += evaluator_->geometric_sol_one_way_bb_bb(res_i, res_j, pose);


    }

  }

}

void
ContextIndependentGeometricSolEnergy::setup_for_packing(
  pose::Pose  & pose,
  utility::vector1< bool > const &,
  utility::vector1< bool > const & designing_residues
) const
{

  bool might_be_designing = false;

	for ( Size ii = 1; ii <= designing_residues.size(); ++ii ) {
		if ( designing_residues[ ii ] ) {
			might_be_designing = true;
			break;
		}
	}

  precalculated_bb_bb_energy_ = 0.0f;

  //if nothing can be designed no reason to precalculate backbone/backbone geom_solv
  if(!might_be_designing) return;

  precalculate_bb_bb_energy_for_design(pose);



}

//void
//ContextIndependentGeometricSolEnergy::setup_for_minimizing(
//										 pose::Pose & pose,
//										 ScoreFunction const & sfxn,
//										 kinematics::MinimizerMapBase const & min_map
//										 ) const
//{
//	using namespace basic::options;
//	using namespace basic::options::OptionKeys;
//
//	//set_nres_mono(pose);
//
//	if ( pose.energies().use_nblist() ) {
//	//if ( true ) {
//		// stash our nblist inside the pose's energies object
//		Energies & energies( pose.energies() );
//
//		// setup the atom-atom nblist
//		NeighborListOP nblist;
//		Real const tolerated_motion = pose.energies().use_nblist_auto_update() ? option[ run::nblist_autoupdate_narrow ] : 1.5;
//		Real const XX = 5.2 + 2 * tolerated_motion;
//		nblist = new NeighborList( min_map.domain_map(), XX*XX, XX*XX, XX*XX);
//		if ( pose.energies().use_nblist_auto_update() ) {
//		//if ( true ) {
//			nblist->set_auto_update( tolerated_motion );
//		}
//		// this partially becomes the EtableEnergy classes's responsibility
//		nblist->setup( pose, sfxn, *this);
//		energies.set_nblist( EnergiesCacheableDataType::GEOM_SOLV_NBLIST, nblist );
//	}
//}

/////////////////////////////////////////////////////////////////////////////////
//void
//ContextIndependentGeometricSolEnergy::setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sfxn ) const
//{
//	pose.update_residue_neighbors();
//}

///
void
ContextIndependentGeometricSolEnergy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn ) const
{
	pose.update_residue_neighbors();
	if ( false ) {
	//if ( pose.energies().use_nblist() ) {
		NeighborList const & nblist( pose.energies().nblist( EnergiesCacheableDataType::GEOM_SOLV_NBLIST ) );
		nblist.prepare_for_scoring( pose, scfxn, *this );
	}
}

///////////////////////////////////////////////////////////////////////////////
bool
ContextIndependentGeometricSolEnergy::defines_score_for_residue_pair(
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
ContextIndependentGeometricSolEnergy::get_count_pair_function(
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
ContextIndependentGeometricSolEnergy::get_count_pair_function(
											conformation::Residue const & rsd1,
											conformation::Residue const & rsd2
											) const
{
	using namespace etable::count_pair;

	if ( ! defines_score_for_residue_pair(rsd1, rsd2, true) ) return new CountPairNone;

	// PUT THIS IN PROPERLY LATER.
	if ( rsd1.is_bonded( rsd2 ) || rsd1.is_pseudo_bonded( rsd2 ) ) {
		//		return CountPairFactory::create_count_pair_function( rsd1, rsd2, CP_CROSSOVER_4 );
	}
	return new CountPairAll;

}

///////////////////////////////////////////////////////////////////////////////
etable::count_pair::CountPairFunctionCOP
ContextIndependentGeometricSolEnergy::get_intrares_countpair(
										   conformation::Residue const & res,
										   pose::Pose const &,
										   ScoreFunction const &
										   ) const
{
	using namespace etable::count_pair;
	return CountPairFactory::create_intrares_count_pair_function( res, CP_CROSSOVER_3 );
}

///////////////////////////////////////////////////////////////////////////////
bool
ContextIndependentGeometricSolEnergy::use_extended_residue_pair_energy_interface() const
{
	return true;
}

///////////////////////////
void
ContextIndependentGeometricSolEnergy::setup_for_minimizing_for_residue(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData &
) const
{}


////////////////////////////////////////////////////////////////////////////////
void
ContextIndependentGeometricSolEnergy::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const &,
	ScoreFunction const &,
	kinematics::MinimizerMapBase const &,
	ResSingleMinimizationData const &,
	ResSingleMinimizationData const &,
	ResPairMinimizationData & pair_data
) const
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//if ( pose.energies().use_nblist_auto_update() ) return;

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

/////////////
bool
ContextIndependentGeometricSolEnergy::requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose::Pose const & ) const
{
	return false;
}

////////////////////
void
ContextIndependentGeometricSolEnergy::eval_intrares_derivatives(
	conformation::Residue const & rsd,
	ResSingleMinimizationData const &,
	pose::Pose const & pose,
	EnergyMap const & weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) const
{
	if ( !defines_intrares_energy( weights ) ) {
		return;
	}

	for ( Size ii=1, ii_end=rsd.natoms(); ii<= ii_end; ++ii ) {
		Vector F1( 0.0 );
		Vector F2( 0.0 );
		id::AtomID atom_id( ii, rsd.seqpos() );
		evaluator_->eval_atom_derivative_intra_RNA(atom_id, pose, weights, F1, F2);
		atom_derivs[ ii ].f1() += F1;
		atom_derivs[ ii ].f2() += F2;
    }
}


////////////////////////////////////////////////////
void
ContextIndependentGeometricSolEnergy::eval_residue_pair_derivatives(
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
	//if ( pose.energies().use_nblist_auto_update() ) return;

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

////////////////////////////////////////////////////////////////////////////////
void
ContextIndependentGeometricSolEnergy::residue_pair_energy_ext(
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
	//if ( pose.energies().use_nblist_auto_update() ) return;
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
	emap[ geom_sol_fast ] += score;
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
/// Everything in here.
void
ContextIndependentGeometricSolEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const
{
  using_extended_method_ = false;
  //if ( pose.energies().use_nblist() ) return;
  //if the backbone/backbone energy has already been calculated in setup_for_packing
  //this is done only if we are doing fix backbone design and the backbone/backbone
  //energy cannot change (Joseph Yesselman 9/11/13)

  if(precalculated_bb_bb_energy_ > 0.0f) {

    emap[ geom_sol_fast ] += evaluator_->geometric_sol_one_way_sc(rsd1, rsd2, pose) +
                             evaluator_->geometric_sol_one_way_sc(rsd2, rsd1, pose);

  }


  else {

    EnergyMap emap_local;
    evaluator_->residue_pair_energy( rsd1, rsd2, pose, scorefxn, emap_local );
    emap[ geom_sol_fast ] += emap_local[ geom_sol ];

  }

}

bool
ContextIndependentGeometricSolEnergy::minimize_in_whole_structure_context( pose::Pose const & ) const
{
	//return pose.energies().use_nblist_auto_update();
	return false;
}

void
ContextIndependentGeometricSolEnergy::finalize_total_energy(
  pose::Pose & pose,
  ScoreFunction const &,
  EnergyMap & totals
) const
{
	//if ( ! pose.energies().use_nblist() || ! pose.energies().use_nblist_auto_update() ) {
	if ( !using_extended_method_ ) {
  		totals[ geom_sol_fast ] += precalculated_bb_bb_energy_;
		return;
  	}

	return;
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

	totals[ geom_sol_fast ] += score;

}

//void
//ContextIndependentGeometricSolEnergy::eval_atom_derivative(
//		 id::AtomID const & atom_id,
//		 pose::Pose const & pose,
//		 kinematics::DomainMap const & domain_map,
//		 ScoreFunction const & scorefxn,
//		 EnergyMap const & weights,
//		 Vector & F1,
//		 Vector & F2
//) const
//{
//	if ( ! pose.energies().use_nblist_auto_update() ) return;
////	//	A hack. Hopefully not a big slow down.
////	EnergyMap weights_local;
////	weights_local[ geom_sol ]           = weights[ geom_sol_fast ];
////	weights_local[ geom_sol_intra_RNA ] = weights[ geom_sol_fast_intra_RNA ];
////
////
////	evaluator_->eval_atom_derivative( atom_id, pose, domain_map, scorefxn, weights_local, F1, F2 );
//    Real energy( 0.0 );
//	hbonds::HBondDerivs deriv;
//
//    if ( defines_intrares_energy( weights ) ) {
//        evaluator_->eval_atom_derivative_intra_RNA(atom_id, pose, weights, F1, F2);
//    }
//
//    Size const i( atom_id.rsd() );
//	conformation::Residue const & ires( pose.residue( i ) );
//	Size const ii( atom_id.atomno() );
//
//	//	Size const nres = pose.total_residue();
//	static bool const update_deriv( true );
//    //assert( pose.energies().use_nblist() );
//	NeighborList const & nblist
//		( pose.energies().nblist( EnergiesCacheableDataType::GEOM_SOLV_NBLIST ) );
//	AtomNeighbors const & nbrs( nblist.atom_neighbors(i,ii) );
//
//    for ( scoring::AtomNeighbors::const_iterator it2=nbrs.begin(),
//         it2e=nbrs.end(); it2 != it2e; ++it2 ) {
//		scoring::AtomNeighbor const & nbr( *it2 );
//		Size const j( nbr.rsd() );
//        if (i == j) continue;
//		Size const jj( nbr.atomno() );
//		conformation::Residue const & jres( pose.residue( j ) );
//
//        if ( evaluator_->atom_is_heavy( jres, jj ) ) {
//            if ( evaluator_->atom_is_donor_h( ires, ii ) ) {
//                evaluator_->get_atom_atom_geometric_solvation_for_donor( ii, ires, jj, jres, pose, energy, update_deriv, deriv );
//                F1 += weights[ geom_sol_fast ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
//                F2 += weights[ geom_sol_fast ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
//            } else if ( evaluator_->atom_is_acceptor( ires, ii ) ) {
//                evaluator_->get_atom_atom_geometric_solvation_for_acceptor( ii, ires, jj, jres, pose, energy, update_deriv, deriv );
//                F1 += weights[ geom_sol_fast ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
//                F2 += weights[ geom_sol_fast ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
//            }
//        }
//        if ( evaluator_->atom_is_heavy ( ires, ii ) ) {
//            if ( evaluator_->atom_is_donor_h( jres, jj ) ) {
//                evaluator_->get_atom_atom_geometric_solvation_for_donor( jj, jres, ii, ires, pose, energy, update_deriv, deriv );
//                F1 -= weights[ geom_sol_fast ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
//                F2 -= weights[ geom_sol_fast ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
//            } else if ( evaluator_->atom_is_acceptor( jres, jj ) ) {
//                evaluator_->get_atom_atom_geometric_solvation_for_acceptor( jj, jres, ii, ires, pose, energy, update_deriv, deriv );
//                F1 -= weights[ geom_sol_fast ] * ( deriv.h_deriv.f1() +  deriv.don_deriv.f1() );
//                F2 -= weights[ geom_sol_fast ] * ( deriv.h_deriv.f2() +  deriv.don_deriv.f2() );
//            }
//        }
//    }
//}

Distance
ContextIndependentGeometricSolEnergy::atomic_interaction_cutoff() const
{
	return evaluator_->atomic_interaction_cutoff();
}

///////////////////////////////////////////////////////////////////////////////////////
bool
ContextIndependentGeometricSolEnergy::defines_intrares_energy( EnergyMap const & weights ) const
{
	//Change to this on Feb 06, 2012. Ensure that the function returns false if weights[geom_sol_intra_RNA] == 0.0
	bool condition_1 = (weights[geom_sol_fast_intra_RNA] > 0.0001) ? true : false;
	return condition_1;
}


///////////////////////////////////////////////////////////////////////////////////////
void
ContextIndependentGeometricSolEnergy::eval_intrares_energy(
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & emap
) const{

	EnergyMap emap_local;
	evaluator_->eval_intrares_energy( rsd, pose, scorefxn, emap_local );
	emap[ geom_sol_fast_intra_RNA ] += emap_local[ geom_sol_intra_RNA ];

}

///@brief ContextIndependentGeometricSolEnergy is not context sensitive, of course.
void
ContextIndependentGeometricSolEnergy::indicate_required_context_graphs(utility::vector1< bool > &) const
{}


core::Size
ContextIndependentGeometricSolEnergy::version() const
{
	return 2; // Initial versioning
}


} // geometric_solvation
} // scoring
} // core

