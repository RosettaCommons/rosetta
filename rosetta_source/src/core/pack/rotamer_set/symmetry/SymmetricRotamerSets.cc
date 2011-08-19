// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/pack/RotamerSet/RotamerSets.cc
/// @brief  RotamerSets class implementation, for symmetric packing
/// @author Ingemar Andre

// Unit Headers
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
// AUTO-REMOVED #include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>


// C++
// AUTO-REMOVED #include <fstream>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

using namespace ObjexxFCL;


namespace core {
namespace pack {
namespace rotamer_set {
namespace symmetry {

SymmetricRotamerSets::SymmetricRotamerSets() {}
SymmetricRotamerSets::~SymmetricRotamerSets() {}

// @details For now, this function knows that all RPEs are computed before annealing begins
// In not very long, this function will understand that there are different ways to compute
// energies and to store them, but for now, it does not.
// symmetrical version of RotamerSets::compute_energies.
void
SymmetricRotamerSets::compute_energies(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scfxn,
	graph::GraphCOP packer_neighbor_graph,
	interaction_graph::InteractionGraphBaseOP ig
)
{
	using namespace interaction_graph;
	using namespace scoring;

	//basic::Tracer tt("core.pack.rotamer_set", basic::t_trace );

	ig->initialize( *this );
	compute_one_body_energies( pose, scfxn, packer_neighbor_graph, ig );

	PrecomputedPairEnergiesInteractionGraphOP pig =
		dynamic_cast< PrecomputedPairEnergiesInteractionGraph * > ( ig.get() );
	if ( pig ) {
		precompute_two_body_energies( pose, scfxn, packer_neighbor_graph, pig );
	} else {
		/// is this an on the fly graph?
		OnTheFlyInteractionGraphOP otfig = dynamic_cast< OnTheFlyInteractionGraph * > ( ig.get() );
		if ( otfig ) {
			utility_exit_with_message("Cannot use symmetry with on-the-fly interaction graph yet!");
		} else {
			utility_exit_with_message("Unknown interaction graph type encountered in RotamerSets::compute_energies()");
		}
	}

}

// @details compute all rotamer one-body interactions for a symmetrical rotamer_sets
void
SymmetricRotamerSets::compute_one_body_energies(
  pose::Pose const & pose,
  scoring::ScoreFunction const & scfxn,
  graph::GraphCOP packer_neighbor_graph,
  interaction_graph::InteractionGraphBaseOP ig
)
{
  // One body energies -- shared between pigs and otfigs
  for ( uint ii = 1; ii <= nmoltenres(); ++ii )
  {
    utility::vector1< core::PackerEnergy > one_body_energies( rotamer_set_for_moltenresidue( ii )->num_rotamers() );
    rotamer_set_for_moltenresidue( ii )->compute_one_body_energies(
      pose, scfxn, *task(), packer_neighbor_graph, one_body_energies );
    ig->add_to_nodes_one_body_energy( ii, one_body_energies );
  }
}

// @details calculate all rotamer two body interactions and place them in a ig. Adds in
// energies for symmetrical clone residues into the energy of the master residues. The
// idea is the following:
void
SymmetricRotamerSets::precompute_two_body_energies(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scfxn,
	graph::GraphCOP packer_neighbor_graph,
	interaction_graph::PrecomputedPairEnergiesInteractionGraphOP pig,
	bool const finalize_edges
)
{
	using namespace interaction_graph;
	using namespace scoring;

	// find SymmInfo
  SymmetricConformation const & SymmConf (
    dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
  SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	// Two body energies
	for ( uint ii = 1; ii <= nmoltenres(); ++ ii )
	{
		//tt << "pairenergies for ii: " << ii << '\n';
		uint ii_resid = moltenres_2_resid( ii );
		// observe that we will loop over only one subunit here
		for ( graph::Graph::EdgeListConstIter
			uli  = packer_neighbor_graph->get_node( ii_resid )->const_upper_edge_list_begin(),
			ulie = packer_neighbor_graph->get_node( ii_resid )->const_upper_edge_list_end();
			uli != ulie; ++uli )
		{
			uint jj_resid = (*uli)->get_second_node_ind();
			uint jj = resid_2_moltenres( jj_resid ); //pretend we're iterating over jj >= ii
//			if ( jj == 0 ) continue; // Andrew, remove this magic number!
			// if jj_resid is not repackable and jj_resid is not in a different subunit continue
			if ( jj == 0 && symm_info->chi_follows( jj_resid ) == 0  ) continue;
			// if jj_resid is in a different subunit we need to know its master
			uint jj_resid_master = jj_resid;
			if ( symm_info->chi_follows( jj_resid ) != 0 ) jj_resid_master = symm_info->chi_follows( jj_resid );
			// find out the moltenres id for the master
			uint jj_master = resid_2_moltenres( jj_resid_master );
			// if the master is not repackable continue
			if ( jj_master == 0 ) continue;
			// if we the ii_resid and jj_resid_master have an edge this intersubunit
			// interaction is already being calculated
//			if ( packer_neighbor_graph->get_edge_exists(ii_resid,jj_resid_master ) && jj == 0 ) continue;

			uint ii_master = ii;
			uint ii_resid_master = ii_resid;

			// the self interaction energy is already calculated in the one-body ter,
			if ( ii_master == jj_master ) continue;
			// swap the order of the nodes that form the edge if jj_resid_master < ii_resid.
			// Edges are always stored as pairs (a, b) where a > b. This only happens if we
			// are calculation interactions from the controling subunit to another subunit
			bool swap( false );
			if ( jj_resid_master < ii_resid ) {
				swap = true;
				uint temp_ii = ii_master;
				uint temp_ii_resid = ii_resid_master;
				ii_master = jj_master;
				ii_resid_master = jj_resid_master;
				jj_master = temp_ii;
				jj_resid_master = temp_ii_resid;
			}

			// make a pair energy table to store use jj_master and ii_master
			// to size the array ( if we have have a intersubunit interaction)
			FArray2D< core::PackerEnergy > pair_energy_table(
			nrotamers_for_moltenres( jj_master ),
			nrotamers_for_moltenres( ii_master ), 0.0 );

			RotamerSetOP ii_rotset = rotamer_set_for_moltenresidue( ii_master );
			RotamerSetOP jj_rotset = rotamer_set_for_moltenresidue( jj_master );

			// we have a intersubunit interaction then calculate the interaction
			// here instead of the intrasubunit interaction. If jj == 0 then we will calculate
			// intrasubunit interaction. If we swapped the order we have to orient ii instead of jj
			if ( symm_info->chi_follows( jj_resid ) != 0 && jj == 0 ) {
				if ( swap ) {
					RotamerSetOP rotated_ii_rotset(
						orient_rotamer_set_to_symmetric_partner(pose,ii_resid_master,jj_resid) );
					ii_rotset = rotated_ii_rotset;
					scfxn.prepare_rotamers_for_packing( pose, *ii_rotset );
				} else {
					RotamerSetOP rotated_jj_rotset(
						orient_rotamer_set_to_symmetric_partner(pose,jj_resid_master,jj_resid) );
					jj_rotset = rotated_jj_rotset;
					scfxn.prepare_rotamers_for_packing( pose, *jj_rotset );
				}
			}
			scfxn.evaluate_rotamer_pair_energies(
				*ii_rotset, *jj_rotset, pose, pair_energy_table );

			// Apply the multiplication factors
			pair_energy_table *= symm_info->score_multiply(ii_resid,jj_resid);
			if ( !pig->get_edge_exists( ii_master, jj_master ) ) {
				pig->add_edge( ii_master, jj_master );
			}
			pig->add_to_two_body_energies_for_edge( ii_master, jj_master, pair_energy_table );

			// finalize the edge
			if ( finalize_edges && ! scfxn.any_lr_residue_pair_energy(pose, ii_master, jj_master )
						&& final_visit_to_edge( pose, packer_neighbor_graph, ii_resid, jj_resid ) ){
				pig->declare_edge_energies_final( ii_master, jj_master );
			}
		}
	}

	// Iterate across the long range energy functions and use the iterators generated
	// by the LRnergy container object
	// The logic here is exactly the same as for the non-lr energy calculation in terms of symmetry
	for ( ScoreFunction::LR_2B_MethodIterator
			lr_iter = scfxn.long_range_energies_begin(),
			lr_end  = scfxn.long_range_energies_end();
			lr_iter != lr_end; ++lr_iter ) {
		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*lr_iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.
		// Potentially O(N^2) operation...

		for ( uint ii = 1; ii <= nmoltenres(); ++ ii ) {
			uint const ii_resid = moltenres_2_resid( ii );

			for ( ResidueNeighborConstIteratorOP
						rni = lrec->const_upper_neighbor_iterator_begin( ii_resid ),
						rniend = lrec->const_upper_neighbor_iterator_end( ii_resid );
						(*rni) != (*rniend); ++(*rni) ) {
				Size jj_resid = rni->upper_neighbor_id();

				uint jj = resid_2_moltenres( jj_resid ); //pretend we're iterating over jj >= ii
//				if ( jj == 0 ) continue; // Andrew, remove this magic number! (it's the signal that jj_resid is not "molten")
			// if jj_resid is not repackable and jj_resid is not in a different subunit continue
			if ( jj == 0 && symm_info->chi_follows( jj_resid ) == 0  ) continue;
			// if jj_resid is in a different subunit we need to know its master
			uint jj_resid_master = jj_resid;
			if ( symm_info->chi_follows( jj_resid ) != 0 ) jj_resid_master = symm_info->chi_follows( jj_resid );
			// find out the moltenres id for the master
			uint jj_master = resid_2_moltenres( jj_resid_master );
			// if the master is not repackable continue
			if ( jj_master == 0 ) continue;
			// if we the ii_resid and jj_resid_master have an edge this intersubunit
			// interaction is already being calculated
//			if ( packer_neighbor_graph->get_edge_exists(ii_resid,jj_resid_master ) && jj == 0 ) continue;

			uint ii_master = ii;
			uint ii_resid_master = ii_resid;

			// the self interaction energy is already calculated in the one-body ter,
	     if ( ii_master == jj_master ) continue;
			// swap the order of the nodes that form the edge if jj_resid_master < ii_resid.
			// Edges are always stored as pairs (a, b) where a > b. This only happens if we
			// are calculation interactions from the controling subunit to another subunit
			bool swap( false );
			if ( jj_resid_master < ii_resid ) {
				swap = true;
				uint temp_ii = ii_master;
				uint temp_ii_resid = ii_resid_master;
				ii_master = jj_master;
				ii_resid_master = jj_resid_master;
				jj_master = temp_ii;
				jj_resid_master = temp_ii_resid;
			}

			// make a pair energy table to store use jj_master and ii_master
			// to size the array ( if we have have a intersubunit interaction)
			FArray2D< core::PackerEnergy > pair_energy_table(
			nrotamers_for_moltenres( jj_master ),
			nrotamers_for_moltenres( ii_master ), 0.0 );

			RotamerSetOP ii_rotset = rotamer_set_for_moltenresidue( ii_master );
			RotamerSetOP jj_rotset = rotamer_set_for_moltenresidue( jj_master );

			// we have a intersubunit interaction then calculate the interaction
			// here instead of the intrasubunit interaction. If jj == 0 then we will calculate
			// intrasubunit interaction. If we swapped the order we have to orient ii instead of jj
			if ( symm_info->chi_follows( jj_resid ) != 0 && jj == 0 ) {
				if ( swap ) {
					RotamerSetOP rotated_ii_rotset(
						orient_rotamer_set_to_symmetric_partner(pose,ii_resid_master,jj_resid) );
					ii_rotset = rotated_ii_rotset;
					scfxn.prepare_rotamers_for_packing( pose, *ii_rotset );
				} else {
					RotamerSetOP rotated_jj_rotset(
						orient_rotamer_set_to_symmetric_partner(pose,jj_resid_master,jj_resid) );
					jj_rotset = rotated_jj_rotset;
					scfxn.prepare_rotamers_for_packing( pose, *jj_rotset );
				}
			}
				(*lr_iter)->evaluate_rotamer_pair_energies(
					*ii_rotset, *jj_rotset, pose, scfxn, scfxn.weights(), pair_energy_table );

			pair_energy_table *= symm_info->score_multiply(ii_resid,jj_resid);

				if ( ! pig->get_edge_exists( ii_master, jj_master ) ) { pig->add_edge( ii_master, jj_master ); }
				pig->add_to_two_body_energies_for_edge( ii_master, jj_master, pair_energy_table );
				if ( finalize_edges && final_visit_to_edge( pose, packer_neighbor_graph, ii_resid, jj_resid ) )
					pig->declare_edge_energies_final( ii_master, jj_master );
			}
		}
	}
}

// @details orients all rotamers in a rotamer_set to a different (symmetrical) position
RotamerSetOP
SymmetricRotamerSets::orient_rotamer_set_to_symmetric_partner(
  pose::Pose const & pose,
  uint const & seqpos,
  uint const & sympos
)
{

  RotamerSetCOP rotset_in = rotamer_set_for_residue( seqpos );
  SymmetricRotamerSetFactory rsf;
  RotamerSetOP sym_rotamer_set = rsf.create_rotamer_set( pose.residue( seqpos ) );
	sym_rotamer_set->set_resid(sympos);
  for (Rotamers::const_iterator
    rot     = rotset_in->begin(),
    rot_end = rotset_in->end();
    rot != rot_end; ++rot) {
      conformation::ResidueOP target_rsd( (*rot)->clone() /*pose.residue( seqpos )*/ );
      target_rsd->orient_onto_residue( pose.residue( sympos ) );
      sym_rotamer_set->add_rotamer( *target_rsd );
  }
  return sym_rotamer_set;
}

// @details is this the final visit to an edge?
bool
SymmetricRotamerSets::final_visit_to_edge(
	pose::Pose const & pose,
	graph::GraphCOP packer_neighbor_graph,
	uint ii_resid,
	uint jj_resid
	)
{
	// find SymmInfo
  SymmetricConformation const & SymmConf (
    dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
  SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	// find the highest jj clone sequence number we could have
	uint jj_highest(jj_resid);
	for ( std::vector< Size>::const_iterator
	  clone     = symm_info->chi_clones( jj_resid ).begin(),
		clone_end = symm_info->chi_clones( jj_resid ).end();
		clone != clone_end; ++clone ){
		if (  *clone > jj_highest &&
					packer_neighbor_graph->get_edge_exists(ii_resid, *clone) )
			jj_highest = *clone;
    }
		// find the highest jj clone sequence number we could have
  uint ii_highest(ii_resid);
  for ( std::vector< Size>::const_iterator
    clone     = symm_info->chi_clones( ii_resid ).begin(),
    clone_end = symm_info->chi_clones( ii_resid ).end();
    clone != clone_end; ++clone ){
    if (  *clone > ii_highest &&
          packer_neighbor_graph->get_edge_exists(*clone, jj_resid) )
      ii_highest = *clone;
    }
		if ( ii_resid == ii_highest && jj_resid == jj_highest ) return true;
		else return false;
}

} // namespace symmetry
} // namespace rotamer_set
} // namespace pack
} // namespace core

