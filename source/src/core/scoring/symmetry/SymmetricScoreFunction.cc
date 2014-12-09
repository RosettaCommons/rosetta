// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.cc
/// @brief  ScoreFunction class definition.
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Modified by Sergey Lyskov


// Unit headers
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>

// Package headers
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/symmetry/SymmetricEnergies.hh>

#include <core/scoring/methods/EnergyMethod.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/scoring/methods/TwoBodyEnergy.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/LREnergyContainer.hh>

// // Project headers
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/id/TorsionID.hh>
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <basic/prof.hh>

// Symmetry extras
// AUTO-REMOVED #include <core/conformation/RotamerSetBase.hh> // TMP HACK
#include <core/scoring/EnvPairPotential.hh>
#include <basic/datacache/CacheableData.hh> //TMP HACK
#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>
//#include <core/scoring/symmetry/NBListCache.hh>
//#include <core/scoring/symmetry/NBListCache.fwd.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>

namespace core {
namespace scoring {
namespace symmetry {

using namespace conformation::symmetry;

static thread_local basic::Tracer TR( "core.scoring.symmetry.SymmetricScoreFunction" );
///////////////////////////////////////////////////////////////////////////////
SymmetricScoreFunction::SymmetricScoreFunction():
	ScoreFunction()
	{}

///////////////////////////////////////////////////////////////////////////////
ScoreFunctionOP
SymmetricScoreFunction::clone() const
{
	SymmetricScoreFunctionOP newscorefxn( new SymmetricScoreFunction );
	newscorefxn->assign( *this );
	return newscorefxn;
}

///////////////////////////////////////////////////////////////////////////////

///@detail INTERNAL USE ONLY
void
SymmetricScoreFunction::assign( SymmetricScoreFunction const & src )
{
	ScoreFunction::assign( src );
	// Add symmetry specific copy items here.
}

///@detail INTERNAL USE ONLY
void
SymmetricScoreFunction::assign( ScoreFunction const & src )
{
	ScoreFunction::assign( src );
}

///////////////////////////////////////////////////////////////////////////////

// to start out, just thinking fullatom energies
//
// NOTE: no freakin rotamer trials inside scoring!
Real
SymmetricScoreFunction::operator()( pose::Pose & pose ) const
{
	// This is tricky! If we are non-symmetric we want to use the regular score function
	// instead.
	if ( !core::pose::symmetry::is_symmetric( pose ) ) {
		TR << "Warning!!! Using a symmetric score function on a non-symmetric pose" << std::endl;
		ScoreFunctionOP asym_score( this->clone_as_base_class() );
		return ( *asym_score )(pose);
	}

	// completely unnecessary temporary hack to force refold if nec. for profiling
	pose.residue( pose.total_residue() );

	PROF_START( basic::SCORE );
	//std::cout << "ScoreFunction::operator()\n";

	// notify the pose that we are starting a score evaluation.
	// also tells the cached energies object about the scoring
	// parameters, for the purposes of invalidating cached data
	// if necessary
	//
	// at this point the dof/xyz-moved information will be converted
	// to a domain map. Energy/neighbor links between pair-moved residues
	// will be deleted, and 1d energies of res-moved residues will be
	// cleared.
	//
	// further structure modification will be prevented until scoring is
	// completed
	//
  PROF_START( basic::SCORE_BEGIN_NOTIFY );
	pose.scoring_begin( *this );
	//std::cout << "ScoreFunction::operator() 1\n";

	if ( pose.energies().total_energy() != 0.0 ) {
		TR.Error << "STARTING SCORE NON-ZERO!" << std::endl;
	}
	PROF_STOP( basic::SCORE_BEGIN_NOTIFY );

	// ensure that the total_energies are zeroed out -- this happens in Energies.scoring_begin()
	// unneccessary pose.energies().total_energies().clear();
	//std::cout << "ScoreFunction::operator() 2\n";

  PROF_START( basic::SCORE_SETUP );
	// do any setup necessary
	setup_for_scoring( pose );
	//std::cout << "ScoreFunction::operator() 3\n";

	// Make some arrays symmetrical before scoring
	correct_arrays_for_symmetry( pose );

	PROF_STOP( basic::SCORE_SETUP );

	// evaluate the residue-residue energies that only exist between
	// neighboring residues
	PROF_START( basic::SCORE_NEIGHBOR_ENERGIES );

	eval_twobody_neighbor_energies( pose );

	PROF_STOP ( basic::SCORE_NEIGHBOR_ENERGIES );

	// Ingemar, I put this in.... -Will
	PROF_START( basic::SCORE_LONG_RANGE_ENERGIES );
	eval_long_range_twobody_energies( pose );
	PROF_STOP ( basic::SCORE_LONG_RANGE_ENERGIES );

  PROF_START( basic::SCORE_ONEBODY_ENERGIES );
	// evaluate the onebody energies -- rama, dunbrack, ...
	eval_onebody_energies( pose );
	PROF_STOP( basic::SCORE_ONEBODY_ENERGIES );

  PROF_START( basic::SCORE_FINALIZE );

	// give energyfunctions a chance update/finalize energies
	// etable nblist calculation is performed here
	for ( AllMethodsIterator iter=all_energies_begin(),
        iter_end= all_energies_end(); iter != iter_end; ++iter ) {
		(*iter)->finalize_total_energy( pose, *this, pose.energies().finalized_energies() );
	}

	correct_finalize_score( pose );

	pose.energies().total_energies() += pose.energies().finalized_energies();

	if ( pose.energies().use_nblist() ) {
		pose.energies().total_energies() += pose.energies().minimization_graph()->fixed_energies();
	}

  PROF_STOP( basic::SCORE_FINALIZE );

	//std::cout << "SymSfxn: " <<  pose.energies().use_nblist() << " ";
	//pose.energies().total_energies().show_if_nonzero_weight( std::cout, weights() );
	//std::cout << std::endl;

	PROF_START( basic::SCORE_DOT );

	// dot the weights with the scores
	pose.energies().total_energy() = pose.energies().total_energies().dot( weights() );
	pose.energies().total_energies()[ total_score ] = pose.energies().total_energy();

	PROF_STOP( basic::SCORE_DOT );

	PROF_START( basic::SCORE_END_NOTIFY );

	// notify that scoring is over
	pose.scoring_end( *this );

  PROF_STOP( basic::SCORE_END_NOTIFY );


	PROF_STOP ( basic::SCORE );

	return pose.energies().total_energy();
}

void
SymmetricScoreFunction::setup_for_minimizing(
	pose::Pose & pose,
	kinematics::MinimizerMapBase const & min_map
) const
{
	bool const new_sym_min( basic::options::option[ basic::options::OptionKeys::optimization::new_sym_min ] );
	/// 1. Initialize the nodes of the minimization graph
	/// 2. Initialize the edges with the short-ranged two-body energies
	/// 3. Initialize the edges with the long-ranged two-body energies
	/// 4. Run setup-for-minimization on the edges of the mingraph; dropping edges with no active 2b enmeths.
	/// 5. Let whole-structure energies initialize themselves

	SymmetricConformation & symm_conf (
		dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfo const & symm_info( * symm_conf.Symmetry_Info() );

	SymmetricEnergies & symm_energies(
		dynamic_cast< SymmetricEnergies & > ( pose.energies()) );


	MinimizationGraphOP g( new MinimizationGraph( pose.total_residue() ) );
	MinimizationGraphOP dg( new MinimizationGraph( pose.total_residue() ) ); // derivative graph

	std::list< methods::EnergyMethodCOP > eval_derivs_with_pose_enmeths;
	for ( AllMethods::const_iterator iter=all_methods().begin(),
			iter_end= all_methods().end(); iter != iter_end; ++iter ) {
		if ((*iter)->defines_high_order_terms( pose ) || (*iter)->minimize_in_whole_structure_context( pose ) )
			eval_derivs_with_pose_enmeths.push_back( *iter );
	}

	EnergyMap fixed_energies; // portions of the score function that will not change over the course of minimization.

	/// Accumulate the portions of the scoring function for those residues with non-zero domain-map values
	/// but also let all energy methods register with each node, since some may require a setup-for-scoring opportunity
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {

		bool accumulate_fixed_energies( symm_info.fa_is_independent(ii) &&
			ii > symm_info.num_total_residues_without_pseudo() );

		setup_for_minimizing_for_node( * g->get_minimization_node( ii ), pose.residue( ii ),
			min_map, pose, accumulate_fixed_energies, fixed_energies );
		g->get_minimization_node( ii )->weight( symm_info.score_multiply_factor() );
		setup_for_minimizing_for_node( * dg->get_minimization_node( ii ), pose.residue( ii ),
			min_map, pose, false, fixed_energies ); // only accumulate once
	}
	g->copy_connectivity(  pose.energies().energy_graph() );
	dg->copy_connectivity( pose.energies().energy_graph() );

	kinematics::DomainMap const & domain_map( min_map.domain_map() );
	for ( core::graph::Graph::EdgeListIter
			edge_iter = g->edge_list_begin(),
			edge_iter_end = g->edge_list_end(),
			dedge_iter = dg->edge_list_begin(),
			dedge_iter_end = dg->edge_list_end(),
			ee_edge_iter = pose.energies().energy_graph().edge_list_begin();
			edge_iter != edge_iter_end; ++edge_iter, ++dedge_iter, ++ee_edge_iter ) {
		Size const node1 = (*edge_iter)->get_first_node_ind();
		Size const node2 = (*edge_iter)->get_second_node_ind();
		assert( node1 == (*ee_edge_iter)->get_first_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering
		assert( node2 == (*ee_edge_iter)->get_second_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering
		assert( node1 == (*dedge_iter)->get_first_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering
		assert( node2 == (*dedge_iter)->get_second_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering
		assert( symm_info.bb_follows( node1 ) == 0 || symm_info.bb_follows( node2 ) == 0 );

		// domain map check here?
		bool const res_moving_wrt_eachother(
			domain_map( node1 ) == 0 ||
			domain_map( node2 ) == 0 ||
			domain_map( node1 ) != domain_map( node2 ) );

		Real edge_weight = symm_info.score_multiply( node1, node2 );
		Real edge_dweight = symm_info.deriv_multiply( node1, node2 );
		if ( new_sym_min ) { // new way
			edge_dweight = edge_weight; // NOTE
			if ( edge_weight != 0.0 ) {
				MinimizationEdge & minedge( static_cast< MinimizationEdge & > (**edge_iter) );
				setup_for_minimizing_sr2b_enmeths_for_minedge(
					pose.residue( node1 ), pose.residue( node2 ),
					minedge, min_map, pose, res_moving_wrt_eachother, true,
					static_cast< EnergyEdge const * > (*ee_edge_iter), fixed_energies, edge_weight );
				minedge.weight( edge_weight );
				minedge.dweight( edge_dweight );
			}
		} else { // classic way
			if ( edge_weight != 0.0 ) {
				MinimizationEdge & minedge( static_cast< MinimizationEdge & > (**edge_iter) );
				setup_for_minimizing_sr2b_enmeths_for_minedge(
					pose.residue( node1 ), pose.residue( node2 ),
					minedge, min_map, pose, res_moving_wrt_eachother, true,
					static_cast< EnergyEdge const * > (*ee_edge_iter), fixed_energies, edge_weight );
				minedge.weight( edge_weight );
				minedge.dweight( edge_dweight );
			} else {
				MinimizationEdge & minedge( static_cast< MinimizationEdge & > (**dedge_iter) );
				setup_for_minimizing_sr2b_enmeths_for_minedge(
					pose.residue( node1 ), pose.residue( node2 ),
					minedge, min_map, pose, res_moving_wrt_eachother, false,
					static_cast< EnergyEdge const * > (*ee_edge_iter), fixed_energies ); // edge weight of 1
				minedge.dweight( edge_dweight );
			}
		}
	}

	/// 3. Long range energies need time to get included into the graph, which may require the addition of new edges
	/// 3a: CILR2B energies
	///    i.   Iterate across all ci2b long range energy methods
	///    ii.  Iterate across all residue pairs indicated in the long-range energy containers
	///    iii. If two residues have the same non-zero domain-map coloring, then accumulate their interaction energy and move on
	///    iv.  otherwise, find the corresponding minimization-graph edge for this residue pair
	///    v.   adding a new edge if necessary,
	///    vi.  and prepare the minimization data for this edge

	for ( LR_2B_MethodIterator
			iter = long_range_energies_begin(),
			iter_end = long_range_energies_end();
			iter != iter_end; ++iter ) {
		// NO! add these terms to the minimization graph for scoring, even if not for derivative evaluation
		///if ( (*iter)->minimize_in_whole_structure_context( pose ) ) continue;

		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue;

		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( ii ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				Size const jj = rni->upper_neighbor_id();
				bool const res_moving_wrt_eachother(
					domain_map( ii ) == 0 ||
					domain_map( jj ) == 0 ||
					domain_map( ii ) != domain_map( jj ) );

				Real edge_weight = symm_info.score_multiply( ii, jj );
				//Real edge_dweight = symm_info.deriv_multiply( ii, jj );
				if ( new_sym_min ) {
					if ( edge_weight != 0.0 ) {
						// adjust/add the edge to the scoring graph
						// use the edge_weight as the deriv weight, too
						setup_for_lr2benmeth_minimization_for_respair(
							pose.residue( ii ), pose.residue( jj ), *iter, *g, min_map, pose,
							res_moving_wrt_eachother, true, rni, fixed_energies, edge_weight, edge_weight );
					}
				} else {
					if ( edge_weight != 0.0 ) {
						// adjust/add the edge to the scoring graph
						setup_for_lr2benmeth_minimization_for_respair(
							pose.residue( ii ), pose.residue( jj ), *iter, *g, min_map, pose,
							res_moving_wrt_eachother, true, rni, fixed_energies, edge_weight );
					} else {
						/// adjust/add this edge to the derivative graph
						setup_for_lr2benmeth_minimization_for_respair(
							pose.residue( ii ), pose.residue( jj ), *iter, *dg, min_map, pose,
							res_moving_wrt_eachother, false, rni, fixed_energies ); // no edge weight
					}
				}
			}
		}
	}

	/// 4a. drop unused edges from the scoring graph; call setup for minimizing on those remaining
	for ( core::graph::Graph::EdgeListIter	edge_iter = g->edge_list_begin(),
			edge_iter_end = g->edge_list_end(); edge_iter != edge_iter_end; /* no increment */ ) {
		Size const node1 = (*edge_iter)->get_first_node_ind();
		Size const node2 = (*edge_iter)->get_second_node_ind();

		MinimizationEdge & minedge( static_cast< MinimizationEdge & > (**edge_iter) );

		core::graph::Graph::EdgeListIter	edge_iter_next( edge_iter );
		++edge_iter_next;

		if ( minedge.any_active_enmeths() ) {
			//std::cout << " active scoring graph edge: " << node1 << " " << node2 << std::endl;
			minedge.setup_for_minimizing( pose.residue(node1), pose.residue(node2), pose, *this, min_map );
		} else {
			/// The edge will not contribute anything to scoring during minimization,
			/// so delete it from the graph, so we don't have to pay the expense of traversing
			/// through it.
			g->delete_edge( *edge_iter );
		}
		edge_iter = edge_iter_next;
	}

	/// 4b. drop unused edges from the derivatives graph; call setup for minimizing on those remaining
	for ( core::graph::Graph::EdgeListIter	edge_iter = dg->edge_list_begin(),
			edge_iter_end = dg->edge_list_end(); edge_iter != edge_iter_end; /* no increment */ ) {
		Size const node1 = (*edge_iter)->get_first_node_ind();
		Size const node2 = (*edge_iter)->get_second_node_ind();

		MinimizationEdge & minedge( static_cast< MinimizationEdge & > (**edge_iter) );

		core::graph::Graph::EdgeListIter	edge_iter_next( edge_iter );
		++edge_iter_next;

		if ( minedge.any_active_enmeths() ) {
			//std::cout << " active deriv graph edge: " << node1 << " " << node2 << std::endl;
			minedge.setup_for_minimizing( pose.residue(node1), pose.residue(node2), pose, *this, min_map );
		} else {
			/// The edge will not contribute anything to derivative evaluation during minimization;
			/// it either represents an interaction that is not changing as the result of minimization,
			/// or the interaction is handled by the scoring graph. Delete it from the graph so
			/// we don't have to pay the expense of traversing through it.
			g->delete_edge( *edge_iter );
		}
		edge_iter = edge_iter_next;
	}

	/// 5.  Whole structure energies and energies that are opting out of the MinimizationGraph
	/// routines get a chance to setup for minimizing (using the entire pose as context) and
	for ( std::list< methods::EnergyMethodCOP >::const_iterator
			iter     = eval_derivs_with_pose_enmeths.begin(),
			iter_end = eval_derivs_with_pose_enmeths.end();
			iter != iter_end; ++iter ) {
		(*iter)->setup_for_minimizing( pose, *this, min_map );
		g->add_whole_pose_context_enmeth( *iter );
	}


	//std::cout << "Fixed energies: ";
	//fixed_energies.show_if_nonzero_weight( std::cout, weights() );
	//std::cout << std::endl;

	g->set_fixed_energies( fixed_energies );
	symm_energies.set_minimization_graph( g );
	symm_energies.set_derivative_graph( dg );

}


///////////////////////////////////////////////////////////////////////////////

void
SymmetricScoreFunction::eval_twobody_neighbor_energies(
	pose::Pose & pose
) const
{
	// cached energies object
	Energies & energies( pose.energies() );

	EnergyMap & total_energies( const_cast< EnergyMap & > ( energies.total_energies() ) );

	SymmetricConformation & SymmConf (
		dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	// the neighbor/energy links
	EnergyGraph & energy_graph( energies.energy_graph() );

	// are we using the atom-atom nblist?
	// renaming for true purpose: are we minimizing -- if so,
	// zero the energies stored on edges in the Energy graph, but do
	// not mark the edges as having had their energies computed
	bool const minimizing( energies.use_nblist() );

	if ( minimizing ) {
		/// When minimizing, do not touch the EnergyGraph -- leave it fixed
		MinimizationGraphCOP g = energies.minimization_graph();
		EnergyMap scratch_emap;
		for ( core::graph::Graph::EdgeListConstIter
				edge_iter = g->const_edge_list_begin(),
				edge_iter_end = g->const_edge_list_end();
				edge_iter != edge_iter_end; ++edge_iter ) {
			Size const node1 = (*edge_iter)->get_first_node_ind();
			Size const node2 = (*edge_iter)->get_second_node_ind();
			conformation::Residue const & rsd1( pose.residue( node1 ));
			conformation::Residue const & rsd2( pose.residue( node2 ));
			MinimizationEdge const & minedge( static_cast< MinimizationEdge const & > (**edge_iter) );

			eval_weighted_res_pair_energy_for_minedge( minedge, rsd1, rsd2, pose, *this, total_energies, scratch_emap );
		}

	} else {
		EnergyMap tbemap;

		for ( Size i=1, i_end = pose.total_residue(); i<= i_end; ++i ) {
			conformation::Residue const & resl( pose.residue( i ) );
			for ( graph::Graph::EdgeListIter
					iru  = energy_graph.get_node(i)->upper_edge_list_begin(),
					irue = energy_graph.get_node(i)->upper_edge_list_end();
					iru != irue; ++iru ) {
				EnergyEdge & edge( static_cast< EnergyEdge & > (**iru) );

				Size const j( edge.get_second_node_ind() );
				conformation::Residue const & resu( pose.residue( j ) );
				// the pair energies cached in the link
				///EnergyMap & emap( energy_edge->energy_map() );
				tbemap.zero( cd_2b_types() );
				tbemap.zero( ci_2b_types() );

				// the context-dependent guys can't be cached, so they are always reevaluated
				eval_cd_2b( resl, resu, pose, tbemap );

				for ( Size ii = 1; ii <= cd_2b_types().size(); ++ii ) {
					tbemap[ cd_2b_types()[ ii ]] *= symm_info->score_multiply(i,j);
				}

				if ( edge.energies_not_yet_computed() ) {
					// energies not yet computed? <-> (during minimization w/ nblist) moving rsd pair
					/// TEMP zero portions;

					if ( minimizing ) {
						// confirm that this rsd-rsd interaction will be included
						// in the atom pairs on the nblist:
						//assert( ( pose.energies().domain_map(i) !=
						//					pose.energies().domain_map(j) ) ||
						//				( pose.energies().res_moved(i) ) );
						// ensure that these are zeroed, since we will hit them at the
						// end inside the nblist calculation
						// they almost certainly should be, since the energies have not
						// yet been computed...
						eval_ci_2b( resl, resu, pose, tbemap );
						for ( Size ii = 1; ii <= ci_2b_types().size(); ++ii ) {
							tbemap[ ci_2b_types()[ ii ]] *= symm_info->score_multiply(i,j);
						}

						edge.store_active_energies( tbemap );
						// do not mark energies as computed!!!!!!!!!!!!!
					} else {
						eval_ci_2b( resl, resu, pose, tbemap );
						for ( Size ii = 1; ii <= ci_2b_types().size(); ++ii ) {
							tbemap[ ci_2b_types()[ ii ]] *= symm_info->score_multiply(i,j);
						}

						edge.store_active_energies( tbemap );
						edge.mark_energies_computed();
					}
				} else {
					/// Read the CI energies from the edge, as they are still valid;
					for ( Size ii = 1; ii <= ci_2b_types().size(); ++ii ) {
						tbemap[ ci_2b_types()[ ii ]] = edge[ ci_2b_types()[ ii ] ];
					}

					/// Save the freshly computed CD energies on the edge
					edge.store_active_energies( tbemap, cd_2b_types() );
				}

				total_energies.accumulate( tbemap, ci_2b_types() );
				total_energies.accumulate( tbemap, cd_2b_types() );
			}
		} // nbrs of i
	} // i=1,nres

}

void
SymmetricScoreFunction::eval_long_range_twobody_energies( pose::Pose & pose ) const
{
	bool const minimizing( pose.energies().use_nblist() );
	if ( minimizing ) return; // long range energies are handled as part of the 2-body energies in the minimization graph

	// find SymmInfo
	SymmetricConformation & SymmConf (
		dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	EnergyMap & total_energies(const_cast< EnergyMap & > (pose.energies().total_energies()));

	for ( CI_LR_2B_MethodIterator iter=ci_lr_2b_methods_begin(),
			iter_end = ci_lr_2b_methods_end(); iter != iter_end; ++iter ) {

		LREnergyContainerOP lrec = pose.energies().nonconst_long_range_container( (*iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.

		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			for ( ResidueNeighborIteratorOP
					rni = lrec->upper_neighbor_iterator_begin( ii ),
					rniend = lrec->upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				EnergyMap emap;
				if ( ! rni->energy_computed() ) {
					Size jj = rni->upper_neighbor_id();
					(*iter)->residue_pair_energy( pose.residue(ii), pose.residue( jj ), pose, *this, emap );
					emap *= symm_info->score_multiply( ii, jj );
					rni->save_energy( emap );

					/// DANGER DANGER DANGER.  use_nblist() now means "In the process of a minimization". There is
					/// no such thing as a non-neighborlist minimization.  This will confuse people at some point.
					/// We ought to have some "are we minimizing currently" flag who's meaning can be decoupled
					/// from the neighborlist idea.
					if ( ! pose.energies().use_nblist() ) {
						// Do we need to do this symmetrically?
						rni->mark_energy_computed();
					}
				} else {
					rni->retrieve_energy( emap ); // pbmod
				}
				total_energies += emap;
			}
		}
	}

	//fpd CD LR methods should always be computed
	for ( CD_LR_2B_MethodIterator iter = cd_lr_2b_methods_begin(),
			iter_end = cd_lr_2b_methods_end(); iter != iter_end; ++iter ) {

		LREnergyContainerOP lrec
			= pose.energies().nonconst_long_range_container( (*iter)->long_range_type() );

		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			for ( ResidueNeighborIteratorOP
					rni = lrec->upper_neighbor_iterator_begin( ii ),
					rniend = lrec->upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {

				EnergyMap emap;
				Size jj = rni->upper_neighbor_id();
				(*iter)->residue_pair_energy( pose.residue(ii), pose.residue(jj), pose, *this, emap );
				emap *= symm_info->score_multiply( ii, jj );

				rni->save_energy( emap );
				//rni->mark_energy_computed();

				total_energies += emap;

			}
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
SymmetricScoreFunction::eval_onebody_energies( pose::Pose & pose ) const
{
	// context independent onebody energies
	Energies & energies( pose.energies() );
	EnergyMap & totals( energies.total_energies() );

	 // find SymmInfo
	SymmetricConformation & SymmConf (
		dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	bool const minimizing( energies.use_nblist() );

	if ( minimizing ) {
		EnergyMap scratch_emap;
		MinimizationGraphCOP mingraph = energies.minimization_graph();
		assert( mingraph );
		for ( Size ii = 1; ii <= symm_info->num_total_residues_without_pseudo(); ++ii ) {
			if ( !symm_info->fa_is_independent(ii) ) continue;

			conformation::Residue const & rsd = pose.residue( ii );

			MinimizationNode const & iiminnode =  * mingraph->get_minimization_node( ii );
			eval_weighted_res_onebody_energies_for_minnode( iiminnode, rsd, pose, *this, totals, scratch_emap );

		}
	} else {
		for ( Size i=1; i<= pose.total_residue(); ++i ) {
			if ( !symm_info->fa_is_independent(i) ||
					i > symm_info->num_total_residues_without_pseudo() ) continue;

			EnergyMap & emap( energies.onebody_energies( i ) );

			// 1body intxns ///////////////////////////
			if ( energies.res_moved( i ) ) {
				// have to recalculate
				emap.clear(); // should already have been done when the domain map was internalized
				eval_ci_1b( pose.residue(i), pose, emap );
				emap *= symm_info->score_multiply(i,i);
			}

			emap.zero( cd_1b_types() ); // cant be cached
			EnergyMap cd_1b_emap;
			eval_cd_1b( pose.residue(i), pose, cd_1b_emap );
			cd_1b_emap *= symm_info->score_multiply_factor();
			emap += cd_1b_emap;

			// 2body energy methods are allowed to define 1body intxns ///////////////////////////
			if ( any_intrares_energies() ) {
				// context independent:
				if ( energies.res_moved( i ) ) {
					EnergyMap ci_intrares_emap;
					eval_ci_intrares_energy( pose.residue(i), pose, ci_intrares_emap );
					ci_intrares_emap *= symm_info->score_multiply_factor();
					emap += ci_intrares_emap;
				}

				// context dependent
				EnergyMap cd_intrares_emap;
				eval_cd_intrares_energy( pose.residue(i), pose, cd_intrares_emap );
				cd_intrares_emap *= symm_info->score_multiply_factor();
				emap += cd_intrares_emap;
			}

			totals += emap;

			energies.reset_res_moved( i ); // mark one body energies as having been calculated
			for ( std::vector< Size>::const_iterator
				clone     = symm_info->bb_clones( i ).begin(),
			clone_end = symm_info->bb_clones( i ).end();
			clone != clone_end; ++clone ){
					energies.reset_res_moved( *clone );
			}
		}
		//std::cout << "totals: "<<  i  << totals;
	}
}


void
SymmetricScoreFunction::setup_for_derivatives( pose::Pose & pose ) const
{
	SymmetricConformation & symm_conf (
		dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	SymmetricEnergies & symm_energies( dynamic_cast< SymmetricEnergies & > ( pose.energies()) );

	/// 1. Call setup_for_derivatives on the (active) nodes of the scoring graph; a node is active
	/// if it's part of the asymmetric unit, or if it has any edge
	MinimizationGraphOP g  = symm_energies.minimization_graph();
	for ( Size ii = 1; ii <= symm_info->num_total_residues_without_pseudo(); ++ii ) {
		if ( symm_info->bb_is_independent( ii ) || g->get_node( ii )->num_edges() != 0 ) {
			g->get_minimization_node( ii )->setup_for_derivatives( pose.residue(ii), pose, *this );
		}
	}
	/// 2. Call setup_for_derivatives on the edges of the scoring graph
	for ( graph::Graph::EdgeListIter
			edgeit = g->edge_list_begin(), edgeit_end = g->edge_list_end();
			edgeit != edgeit_end; ++edgeit ) {
		MinimizationEdge & minedge = static_cast< MinimizationEdge & > ( (**edgeit) );
		minedge.setup_for_derivatives(
			pose.residue(minedge.get_first_node_ind()),
			pose.residue(minedge.get_second_node_ind()),
			pose, *this );
	}

	MinimizationGraphOP dg = symm_energies.derivative_graph();
	/// 3. Call setup_for_derivatives on the (active) nodes of the derivatives graph; here, if a node in the asymmetric unit
	/// has no edges, then all of its derivaties will be handled by the scoring graph; we don't need to call
	/// setup for derivatives on that node.
	for ( Size ii = 1; ii <= symm_info->num_total_residues_without_pseudo(); ++ii ) {
		if ( dg->get_node( ii )->num_edges() != 0 ) {
			dg->get_minimization_node( ii )->setup_for_derivatives( pose.residue(ii), pose, *this );
		}
	}
	/// 4. Call setup_for_derivatives on the edges of the derivatives graph
	for ( graph::Graph::EdgeListIter
			edgeit = dg->edge_list_begin(), edgeit_end = dg->edge_list_end();
			edgeit != edgeit_end; ++edgeit ) {
		MinimizationEdge & minedge = static_cast< MinimizationEdge & > ( (**edgeit) );
		minedge.setup_for_derivatives(
			pose.residue(minedge.get_first_node_ind()),
			pose.residue(minedge.get_second_node_ind()),
			pose, *this );
 	}

	/// 5. Whole-pose-context energies should be allowed to setup for derivatives now.
	for ( MinimizationGraph::Energies::const_iterator
			iter     = g->whole_pose_context_enmeths_begin(),
			iter_end = g->whole_pose_context_enmeths_end();
				iter != iter_end; ++iter ) {
		(*iter)->setup_for_derivatives( pose, *this );
	}
}

///////////////////////////////////////////////////////////////////////////////
void
SymmetricScoreFunction::eval_npd_atom_derivative(
	id::AtomID const & atom_id,
	pose::Pose const & pose,
	kinematics::DomainMap const & domain_map,
	Vector & F1,
	Vector & F2
) const
{
	//std::cout << "SymmetricScoreFunction::eval_atom_derivative " << atom_id.rsd() << " " << atom_id.atomno() << " " << F1.x() << " " << F2.x() << std::endl;

	SymmetricEnergies const & symm_energies( dynamic_cast< SymmetricEnergies const & > (pose.energies()) );

	assert( pose.energies().minimization_graph() );
	MinimizationGraphCOP mingraph  = symm_energies.minimization_graph();
	//MinimizationGraphCOP dmingraph = symm_energies.derivative_graph();

	/*Size const rsdno = atom_id.rsd();
	Size const atomno = atom_id.atomno();

	conformation::Residue const & rsd = pose.residue( rsdno );

	MinimizationNode const & minnode =  * mingraph->get_minimization_node( rsdno );
	/// 1. eval intra-residue derivatives
	eval_atom_derivative_for_minnode( minnode, atomno, rsd, pose, domain_map, *this, weights(), F1, F2 );

	ResSingleMinimizationData const & ressingle_min_data( minnode.res_min_data() );
	/// 2a. eval inter-residue derivatives using edges from the scoring graph
	for ( graph::Node::EdgeListConstIter
			edgeit = minnode.const_edge_list_begin(), edgeit_end = minnode.const_edge_list_end();
			edgeit != edgeit_end; ++edgeit ) {
		Vector f1(0.0), f2(0.0);
		MinimizationEdge const & minedge = static_cast< MinimizationEdge const & > ( (**edgeit) );
		Size const other_rsdno = minedge.get_other_ind( rsdno );
		conformation::Residue const & other_rsd( pose.residue( other_rsdno ));
		ResSingleMinimizationData const & other_ressingle_min_data( mingraph->get_minimization_node( other_rsdno )->res_min_data() );

		eval_atom_deriv_for_minedge( minedge, atomno, rsd, other_rsd,
			ressingle_min_data, other_ressingle_min_data,
			pose, domain_map, *this, weights(), f1, f2 );
		//std::cout << "  scoring minedge with " << other_rsdno << " " << f1.x() << " " << f2.x() << std::endl;
		F1 += f1;
		F2 += f2;
	}

	/// 2b. eval inter-residue derivatives using edges from the scoring graph
	ResSingleMinimizationData const & dressingle_min_data( dmingraph->get_minimization_node( rsdno )->res_min_data() );
	for ( graph::Node::EdgeListConstIter
			edgeit = dmingraph->get_minimization_node( rsdno )->const_edge_list_begin(),
			edgeit_end = dmingraph->get_minimization_node( rsdno )->const_edge_list_end();
			edgeit != edgeit_end; ++edgeit ) {
		Vector f1(0.0), f2(0.0);
		MinimizationEdge const & minedge = static_cast< MinimizationEdge const & > ( (**edgeit) );
		Size const other_rsdno = minedge.get_other_ind( rsdno );
		conformation::Residue const & other_rsd( pose.residue( other_rsdno ));
		ResSingleMinimizationData const & other_dressingle_min_data( dmingraph->get_minimization_node( other_rsdno )->res_min_data() );

		eval_atom_deriv_for_minedge( minedge, atomno, rsd, other_rsd,
			dressingle_min_data, other_dressingle_min_data,
			pose, domain_map, *this, weights(), f1, f2 );
		//std::cout << "  deriv minedge with " << other_rsdno << " " << f1.x() << " " << f2.x() << std::endl;
		F1 += f1;
		F2 += f2;
	}*/

	/// 3. Whole-pose-context energies should have their contribution calculated here.
	for ( MinimizationGraph::Energies::const_iterator
			iter     = mingraph->whole_pose_context_enmeths_begin(),
			iter_end = mingraph->whole_pose_context_enmeths_end();
			iter != iter_end; ++iter ) {
		(*iter)->eval_atom_derivative( atom_id, pose, domain_map, *this, weights(), F1, F2 );
	}

}

///////////////////////////////////////////////////////////////////////////////
Real
SymmetricScoreFunction::eval_dof_derivative(
	id::DOF_ID const & dof_id,
	id::TorsionID const & torsion_id,
	pose::Pose const & pose
) const
{

	assert( pose.energies().minimization_graph() );
	MinimizationGraphCOP mingraph = pose.energies().minimization_graph();

	Size const rsdno = torsion_id.valid() ? torsion_id.rsd() : dof_id.atom_id().rsd();
	conformation::Residue const & rsd = pose.residue( rsdno );

	MinimizationNode const & minnode = * mingraph->get_minimization_node( rsdno );

	return eval_weighted_dof_deriv_for_minnode( minnode, rsd, pose, dof_id, torsion_id, *this, weights() );
}

void
SymmetricScoreFunction::intersubunit_hbond_energy(
	pose::Pose & pose,
	EnergyMap & intersubunit_energy
) const
{

	using EnergiesCacheableDataType::HBOND_SET;

	hbonds::HBondSet const & hbond_set(
		static_cast< hbonds::HBondSet const & > ( pose.energies().data().get( HBOND_SET )));

	SymmetricConformation & SymmConf (
		dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );


	for ( Size i = 1; i <= hbond_set.nhbonds(); ++i ) {
		hbonds::HBond const & hbond( hbond_set.hbond(i) );

		Real sr_bb_hbenergy = 0.0;
		Real lr_bb_hbenergy = 0.0;

		switch ( get_hbond_weight_type( hbond.eval_type() ) ) {
			case hbonds::hbw_SR_BB:
				sr_bb_hbenergy = hbond.energy() * hbond.weight();
				break;
			case hbonds::hbw_LR_BB:
				lr_bb_hbenergy = hbond.energy() * hbond.weight();
				break;
			case hbonds::hbw_SR_BB_SC:
			case hbonds::hbw_LR_BB_SC:
			case hbonds::hbw_SC:
				// dont care about these
				break;
			default:
				TR.Warning << "Warning: energy from unexpected HB type ignored " << hbond.eval_type() << std::endl;
				break;
		}
		Real factor (0.0);

		// get the score factor for this edge
		factor = symm_info->score_multiply( hbond.don_res() , hbond.acc_res() );

		// adjust for self-hbonds
		//fpd  is this necessary?
		if ( symm_info->bb_follows( hbond.acc_res() ) == hbond.don_res() ||
				symm_info->bb_follows( hbond.don_res() ) == hbond.acc_res() ) {
			factor = symm_info->score_multiply_factor() - 1;
		}
		intersubunit_energy[ hbond_lr_bb ] += factor*lr_bb_hbenergy;
		intersubunit_energy[ hbond_sr_bb ] += factor*sr_bb_hbenergy;
	}
}

void
SymmetricScoreFunction::symmetrical_allow_hbonds( pose::Pose & pose ) const
{
	using EnergiesCacheableDataType::HBOND_SET;

	hbonds::HBondSet & hbond_set
		( static_cast< hbonds::HBondSet & > ( pose.energies().data().get( HBOND_SET )));

	SymmetricConformation & SymmConf (
  	dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	for ( Size i = 1; i <= hbond_set.nhbonds(); ++i ) {
		hbonds::HBond const & hbond( hbond_set.hbond(i) );
		Size acc( hbond.acc_res() );
		Size don( hbond.don_res() );
		if ( symm_info->fa_is_independent( acc ) ) {
			for ( std::vector<Size>::const_iterator clone=symm_info->bb_clones(acc).begin(), clone_end=symm_info->bb_clones(acc).end();
			      clone != clone_end; ++clone ) {
				hbond_set.set_backbone_backbone_acceptor( *clone, hbond_set.acc_bbg_in_bb_bb_hbond( acc ) );
			}
		}
		if ( symm_info->fa_is_independent( don ) ) {
			for ( std::vector<Size>::const_iterator clone=symm_info->bb_clones(don).begin(), clone_end=symm_info->bb_clones(don).end();
			      clone != clone_end; ++clone ) {
				hbond_set.set_backbone_backbone_acceptor( *clone, hbond_set.acc_bbg_in_bb_bb_hbond( don ) );
			}
		}
	}
	pose.energies().data().set( HBOND_SET, basic::datacache::CacheableDataOP( hbond_set.get_self_ptr() ) );
}

void
SymmetricScoreFunction::set_symmetric_residue_neighbors_hbonds( pose::Pose & pose ) const
{

	using EnergiesCacheableDataType::HBOND_SET;

	hbonds::HBondSet & hbond_set
		( static_cast< hbonds::HBondSet & >
     ( pose.energies().data().get( HBOND_SET )));

  SymmetricConformation & SymmConf (
        dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
  SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	  for ( uint res = 1; res <= pose.total_residue(); ++res ) {
    if ( symm_info->get_use_symmetry() ) {
      if ( !symm_info->fa_is_independent( res ) ) {
        int symm_res ( symm_info->bb_follows( res ) );
				int neighbors_symm ( hbond_set.nbrs( symm_res ) );
				hbond_set.set_nbrs( res, neighbors_symm );
			}
		}
  }
	pose.energies().data().set( HBOND_SET, basic::datacache::CacheableDataOP( hbond_set.get_self_ptr() ) );
}

void
SymmetricScoreFunction::set_symmetric_cenlist( pose::Pose & pose ) const
{
//	EnvPairPotential pairpot;
//	pairpot.compute_centroid_environment_symmetric( pose );

	//using core::pose::datacache::CacheableDataType::CEN_LIST_INFO;

	CenListInfoOP cenlist(
		utility::pointer::static_pointer_cast< CenListInfo >( 
			pose.data().get_ptr( core::pose::datacache::CacheableDataType::CEN_LIST_INFO )
		)
	);

  SymmetricConformation & SymmConf (
        dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	  for ( uint res = 1; res <= pose.total_residue(); ++res ) {
    if ( symm_info->get_use_symmetry() ) {
      if ( !symm_info->fa_is_independent( res ) ) {
        int symm_res ( symm_info->bb_follows( res ) );
				double fcen6_symm ( cenlist->fcen6( symm_res) );
				double fcen10_symm ( cenlist->fcen10( symm_res ) );
				double fcen12_symm ( cenlist->fcen12( symm_res ) );
				cenlist->set_fcen6( res, fcen6_symm );
				cenlist->set_fcen10( res, fcen10_symm );
				cenlist->set_fcen12( res, fcen12_symm );
			}
		}
  }
	pose.data().set( core::pose::datacache::CacheableDataType::CEN_LIST_INFO, cenlist );
}

void
SymmetricScoreFunction::correct_arrays_for_symmetry( pose::Pose & pose ) const
{
	//using core::pose::datacache::CacheableDataType::CEN_LIST_INFO;

	if ( pose.data().has( core::pose::datacache::CacheableDataType::CEN_LIST_INFO ) ) {
		set_symmetric_cenlist( pose );
	}

	if ( has_nonzero_weight( hbond_lr_bb ) || has_nonzero_weight( hbond_sr_bb )  ||
			 has_nonzero_weight( hbond_bb_sc ) || has_nonzero_weight( hbond_sc ) ) {
		symmetrical_allow_hbonds( pose );
		set_symmetric_residue_neighbors_hbonds( pose );
	}

}


// energy methods that compute scores in 'finalize_total_energy' (incl all whole structure energies)
//    should not be scaled by # subunits
// this method deals with these methods
void
SymmetricScoreFunction::correct_finalize_score( pose::Pose & pose ) const
{
	EnergyMap delta_energy( pose.energies().total_energies() );
	EnergyMap interface_energies;//, etable_energy;

	SymmetricConformation & SymmConf (
		dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
	Real const factor ( symm_info->score_multiply_factor() - 1 );

	if ( ! pose.energies().use_nblist() ) {
		/// apl mod -- hbonds now calculate bb/bb hbonds during residue_pair_energy_ext, so this is no longer necessary during minimization
		interface_energies[ hbond_lr_bb ] = pose.energies().finalized_energies()[hbond_lr_bb];
		interface_energies[ hbond_sr_bb ] = pose.energies().finalized_energies()[hbond_sr_bb];
	}
	interface_energies[ interchain_contact ] = pose.energies().finalized_energies()[interchain_contact];

	if ( pose.energies().use_nblist()  ) {
		// apl mod -- old style neighborlist has been depricated
		//etable_energy[fa_atr] = pose.energies().finalized_energies()[fa_atr];
		//etable_energy[fa_rep] = pose.energies().finalized_energies()[fa_rep];
		//etable_energy[fa_sol] = pose.energies().finalized_energies()[fa_sol];
		//etable_energy[fa_intra_rep] = pose.energies().finalized_energies()[fa_intra_rep];
		//etable_energy *= factor;
	}

	delta_energy = pose.energies().finalized_energies();
	delta_energy -= interface_energies;
	delta_energy *= factor;
	//delta_energy -= etable_energy;

	if ( (has_nonzero_weight( hbond_lr_bb ) || has_nonzero_weight( hbond_sr_bb ) ) && ! pose.energies().use_nblist() ) {
		EnergyMap new_interface_energy;
		intersubunit_hbond_energy(pose,new_interface_energy);
		delta_energy += new_interface_energy;
	}

	if ( has_nonzero_weight( interchain_contact ) ) {
		EnergyMap interchain_contact_energy;
		interchain_contact_energy[interchain_contact] = pose.energies().finalized_energies()[interchain_contact];
		delta_energy += interchain_contact_energy;
	}

	// whole structure energy correction
	for ( WS_MethodIterator iter=ws_methods_begin(),
			iter_end = ws_methods_end(); iter != iter_end; ++iter ) {
		ScoreTypes type_i = (*iter)->score_types();
		for ( Size j=1; j<=type_i.size(); ++j) {
			EnergyMap deldel_energy;
			deldel_energy[ type_i[j] ] = pose.energies().finalized_energies()[ type_i[j] ]*factor;
			delta_energy -= deldel_energy;
		}
	}

	pose.energies().finalized_energies() -= interface_energies;
	pose.energies().finalized_energies() += delta_energy;

}

} // symmetry
} // namespace scoring
} // namespace core
