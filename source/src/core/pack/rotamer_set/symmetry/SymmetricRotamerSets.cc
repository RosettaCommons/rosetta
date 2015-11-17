// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/pack/interaction_graph/SymmOnTheFlyInteractionGraph.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>


// C++

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>


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
		utility::pointer::dynamic_pointer_cast< core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraph > ( ig );
	if ( pig ) {
		precompute_two_body_energies( pose, scfxn, packer_neighbor_graph, pig );
	} else {
		SymmOnTheFlyInteractionGraphOP symotfig =
			utility::pointer::dynamic_pointer_cast< core::pack::interaction_graph::SymmOnTheFlyInteractionGraph > ( ig );
		if ( symotfig ) {
			prepare_symm_otf_interaction_graph( pose, scfxn, packer_neighbor_graph, symotfig );
		} else {
			utility_exit_with_message( "Encountered incompatible interaction graph type in SymmetricRotamerSets::compute_energies" );
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
	for ( uint ii = 1; ii <= nmoltenres(); ++ii ) {
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
	for ( uint ii = 1; ii <= nmoltenres(); ++ii ) {
		//tt << "pairenergies for ii: " << ii << '\n';
		uint ii_resid = moltenres_2_resid( ii );
		// observe that we will loop over only one subunit here
		for ( graph::Graph::EdgeListConstIter
				uli  = packer_neighbor_graph->get_node( ii_resid )->const_edge_list_begin(),
				ulie = packer_neighbor_graph->get_node( ii_resid )->const_edge_list_end();
				uli != ulie; ++uli ) {
			uint jj_resid = (*uli)->get_other_ind( ii_resid );
			uint jj = resid_2_moltenres( jj_resid ); //pretend we're iterating over jj >= ii
			// if jj_resid is not repackable and jj_resid is not in a different subunit continue
			uint jj_resid_master(0);
			if ( symm_info->chi_follows( jj_resid ) == 0  ) {
				// jj_resid is independent like ii_resid
				if ( jj == 0 ) continue;
				if ( jj_resid <= ii_resid ) continue; // we will hit this in the other order
				jj_resid_master = jj_resid;
			} else {
				// if jj_resid is in a different subunit we need to know its master
				jj_resid_master = symm_info->chi_follows( jj_resid );
			}
			// find out the moltenres id for the master
			uint jj_master = resid_2_moltenres( jj_resid_master );
			// if the master is not repackable continue
			if ( jj_master == 0 ) continue;
			// if we the ii_resid and jj_resid_master have an edge this intersubunit
			// interaction is already being calculated
			//   if ( packer_neighbor_graph->get_edge_exists(ii_resid,jj_resid_master ) && jj == 0 ) continue;

			if ( fabs( symm_info->score_multiply(ii_resid,jj_resid) )<1e-3 ) continue;
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
			//if ( symm_info->chi_follows( jj_resid ) != 0 ) { runtime_assert( jj == 0 ); } // otherwise clone pos packable
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
					rni = lrec->const_neighbor_iterator_begin( ii_resid ),
					rniend = lrec->const_neighbor_iterator_end( ii_resid );
					(*rni) != (*rniend); ++(*rni) ) {
				Size jj_resid = ( rni->upper_neighbor_id() == ii_resid ? rni->lower_neighbor_id() : rni->upper_neighbor_id() );

				uint jj = resid_2_moltenres( jj_resid ); //pretend we're iterating over jj >= ii
				//    if ( jj == 0 ) continue; // Andrew, remove this magic number! (it's the signal that jj_resid is not "molten")
				// if jj_resid is not repackable and jj_resid is not in a different subunit continue
				uint jj_resid_master(0);
				if ( symm_info->chi_follows( jj_resid ) == 0  ) {
					// jj_resid is independent like ii_resid
					if ( jj == 0 ) continue;
					if ( jj_resid <= ii_resid ) continue; // we will hit this in the other order
					jj_resid_master = jj_resid;
				} else {
					// if jj_resid is in a different subunit we need to know its master
					jj_resid_master = symm_info->chi_follows( jj_resid );
				}
				// find out the moltenres id for the master
				uint jj_master = resid_2_moltenres( jj_resid_master );
				// if the master is not repackable continue
				if ( jj_master == 0 ) continue;
				// if we the ii_resid and jj_resid_master have an edge this intersubunit
				// interaction is already being calculated
				//   if ( packer_neighbor_graph->get_edge_exists(ii_resid,jj_resid_master ) && jj == 0 ) continue;

				if ( fabs( symm_info->score_multiply(ii_resid,jj_resid)) <1e-3 ) continue;

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

			}
		}
	}


	/// yes, this is O(N^2) but the logic for doing it correctly in the loops above is really hairy
	/// for the time being, this works.
	if ( !finalize_edges ) return;

	for ( Size ii=1; ii<= nmoltenres(); ++ii ) {
		for ( Size jj=ii+1; jj<= nmoltenres(); ++jj ) {
			if ( pig->get_edge_exists( ii, jj ) ) {
				pig->declare_edge_energies_final( ii, jj );
			}
		}
	}
}

/// @details Add edges between all adjacent nodes in the
/// interaction graph, and note which of the subunit pairs are interacting.
void
SymmetricRotamerSets::prepare_symm_otf_interaction_graph(
	pose::Pose const & pose,
	scoring::ScoreFunction const & sfxn,
	graph::GraphCOP packer_neighbor_graph,
	interaction_graph::SymmOnTheFlyInteractionGraphOP ig
)
{
	// find SymmInfo
	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	/// Mark all-amino-acid nodes as worthy of distinguishing between bb & sc
	for ( Size ii = 1; ii <= nmoltenres(); ++ii ) {
		RotamerSetCOP ii_rotset = rotamer_set_for_moltenresidue( ii );
		bool all_canonical_aas( true );
		for ( Size jj = 1, jje = ii_rotset->get_n_residue_types(); jj <= jje; ++jj ) {
			conformation::ResidueCOP jj_rotamer = ii_rotset->rotamer( ii_rotset->get_residue_type_begin( jj ) );
			if ( jj_rotamer->aa() > chemical::num_canonical_aas ) {
				all_canonical_aas = false;
			}
		}
		if ( all_canonical_aas ) {
			ig->distinguish_backbone_and_sidechain_for_node( ii, true );
		}
	}


	for ( Size ii = 1; ii <= nmoltenres(); ++ii ) {
		Size ii_resid = moltenres_2_resid( ii );
		//std::cout << "ii_resid: " << ii_resid << std::endl;
		for ( graph::Graph::EdgeListConstIter
				li = packer_neighbor_graph->get_node( ii_resid )->const_edge_list_begin(),
				lie = packer_neighbor_graph->get_node( ii_resid )->const_edge_list_end(); li != lie; ++li ) {
			Size jj_resid = (*li)->get_other_ind( ii_resid );
			Size jj = resid_2_moltenres( jj_resid );
			Size jj_resid_master = symm_info->chi_follows( jj_resid ) == 0 ? jj_resid : symm_info->chi_follows( jj_resid );
			Size jj_master = jj == 0 ? resid_2_moltenres( jj_resid_master ) : jj;

			// skip jj if it is a background residue.  This is true
			// if jj == 0 (i.e. not a moltenresidue) and if jj_resid
			// is in the asymmetric unit (i.e. chi_follows(jj_resid) == 0 )
			// If jj is not in the asymmetric unit, and its master is
			// also not a molten residue
			if ( jj == 0 && symm_info->chi_follows( jj_resid ) == 0 ) { /*std::cout << "jj==0 && symm_info->chi_follows( jj_resid ) == 0" << std::endl;*/ continue; }
			// ok, jj_resid's master is also not a molten residue
			if ( jj_master == 0 ) { /*std::cout << "jj_master == 0" << std::endl;*/ continue; }
			// or if jj_master == ii, then we're looking at a residue interacting with its symmetric clone, which
			// in the context of packing, is qualified as a one-body interaction,
			if ( jj_master == ii ) { /*std::cout << "jj_master == ii" << std::endl;*/ continue; }

			// If jj is in the asymmetric unit, then only consider it if it's an upper neighbor of ii
			if ( symm_info->chi_follows( jj_resid ) == 0 && jj_resid < ii_resid ) continue;

			// OK: ii_resid interacts with jj_resid and their interaction
			// needs to be counted as part of the interaction graph.
			Size ii_master( ii ), ii_resid_master( ii_resid );
			bool swap( false );
			if ( jj_resid_master < ii_resid_master ) {
				Size temp;
				temp = ii_master; ii_master = jj_master; jj_master = temp;
				temp = ii_resid_master; ii_resid_master =jj_resid_master; jj_resid_master = temp;
				swap = true;
			}
			//std::cout << "OK we have a winner: " << ii_master << " " << jj_master << " from " << ii_resid << " " << jj_resid << std::endl;

			// next, check if the edge already exists, and, if it does,
			// then continue
			if ( ! ig->get_edge_exists( ii_master, jj_master ) ) {
				ig->add_edge( ii_master, jj_master );
				ig->note_short_range_interactions_exist_for_edge( ii_master, jj_master );
			}
			// now let's inform the graph that, for this ii_resid / jj_resid pair that
			// an interaction exists.

			ig->set_residues_adjacent_for_subunit_pair_for_edge( ii_master, jj_master,
				swap ? 2 : 1,
				symm_info->subunit_index( jj_resid ));
		}
	}

	// Iterate across the long range energy functions and use the iterators generated
	// by the LRnergy container object
	for ( scoring::ScoreFunction::LR_2B_MethodIterator
			lr_iter = sfxn.long_range_energies_begin(),
			lr_end  = sfxn.long_range_energies_end();
			lr_iter != lr_end; ++lr_iter ) {
		scoring::LREnergyContainerCOP lrec = pose.energies().long_range_container( (*lr_iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue; // only score non-empty energies.
		// Potentially O(N^2) operation...

		for ( uint ii = 1; ii <= nmoltenres(); ++ ii ) {
			Size const ii_resid = moltenres_2_resid( ii );

			for ( scoring::ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( ii_resid ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii_resid );
					(*rni) != (*rniend); ++(*rni) ) {

				Size jj_resid = rni->upper_neighbor_id();
				Size jj = resid_2_moltenres( jj_resid );

				// skip jj if this is a background residue.  This is true
				// if jj == 0 (i.e. not a moltenresidue) and if jj_resid
				// is in the asymmetric unit (i.e. chi_follows(jj_resid) == 0 )
				// or if jj is not in the asymmetric unit, and its master is
				// also not a molten residue
				if ( jj == 0 && symm_info->chi_follows( jj_resid ) == 0 ) { /*std::cout << "jj==0 && symm_info->chi_follows( jj_resid ) == 0" << std::endl;*/ continue; }
				Size jj_resid_master = symm_info->chi_follows( jj_resid ) == 0 ? jj_resid : symm_info->chi_follows( jj_resid );
				Size jj_master = jj == 0 ? resid_2_moltenres( jj_resid_master ) : jj;
				// ok, jj_resid's master is also not a molten residue
				if ( jj_master == 0 ) { /*std::cout << "jj_master == 0" << std::endl;*/ continue; }
				// or if jj_master == ii, then we're looking at a residue interacting with its symmetric clone, which
				// in the context of packing, is qualified as a one-body interaction,
				if ( jj_master == ii ) { /*std::cout << "jj_master == ii" << std::endl;*/ continue; }

				// OK: ii_resid interacts with jj_resid and their interaction
				// needs to be counted as part of the interaction graph.
				Size ii_master( ii ), ii_resid_master( ii_resid );
				if ( jj_resid_master < ii_resid_master ) {
					std::swap( ii_master, jj_master );
				}
				//std::cout << "OK we have a winner: " << ii_master << " " << jj_master << " from " << ii_resid << " " << jj_resid << std::endl;

				/// figure out how to notify the lrec that it should only consider these two residues neighbors
				/// for the sake of packing...
				if ( ! ig->get_edge_exists( ii_master, jj_master ) ) {
					ig->add_edge( ii_master, jj_master );
				}
				ig->note_long_range_interactions_exist_for_edge( ii_master, jj_master );
			}
		}
	}

	compute_proline_correction_energies_for_otf_graph( pose, symm_info, sfxn, packer_neighbor_graph, ig );
}

void
SymmetricRotamerSets::compute_proline_correction_energies_for_otf_graph(
	pose::Pose const & pose,
	conformation::symmetry::SymmetryInfoCOP symm_info,
	scoring::ScoreFunction const & scfxn,
	graph::GraphCOP packer_neighbor_graph,
	interaction_graph::SymmOnTheFlyInteractionGraphOP otfig
)
{
	using namespace conformation;
	using namespace conformation::symmetry;
	using namespace interaction_graph;

	Size const nsubunits = symm_info->subunits();

	utility::vector1< Size > example_gly_rotamer_inds( nmoltenres(), 0 );
	utility::vector1< Size > example_pro_rotamer_inds( nmoltenres(), 0 );
	utility::vector1< utility::vector1< ResidueOP > > example_gly_rotamers( nsubunits );
	utility::vector1< utility::vector1< ResidueOP > > example_pro_rotamers( nsubunits );
	for ( Size ii = 1; ii <= nsubunits; ++ii ) {
		example_gly_rotamers[ ii ].resize( nmoltenres() );
		example_pro_rotamers[ ii ].resize( nmoltenres() );
	}

	/// 0. Find example proline and glycine residues where possible.
	for ( Size ii = 1; ii <= nmoltenres(); ++ii ) {
		RotamerSetCOP ii_rotset = rotamer_set_for_moltenresidue( ii );

		Size potential_gly_replacement( 0 );
		for ( Size jj = 1, jje = ii_rotset->get_n_residue_types(); jj <= jje; ++jj ) {
			ResidueCOP jj_rotamer = ii_rotset->rotamer( ii_rotset->get_residue_type_begin( jj ) );
			if ( jj_rotamer->aa() == chemical::aa_gly ) {
				example_gly_rotamer_inds[ ii ] = ii_rotset->get_residue_type_begin( jj );
			} else if ( jj_rotamer->aa() == chemical::aa_pro ) {
				example_pro_rotamer_inds[ ii ] = ii_rotset->get_residue_type_begin( jj );
			} else if ( jj_rotamer->is_protein() ) {
				potential_gly_replacement = ii_rotset->get_residue_type_begin( jj );
			}
		}
		/// failed to find a glycine backbone for this residue, any other
		/// backbone will do, so long as its a protein backbone.
		if ( example_gly_rotamer_inds[ ii ] == 0 ) {
			example_gly_rotamer_inds[ ii ] = potential_gly_replacement;
		}
		/// clone a representative rotamer from each of the subunits from the symotofig
		for ( Size jj = 1; jj <= nsubunits; ++jj ) {
			if ( example_pro_rotamer_inds[ ii ] != 0 ) {
				example_pro_rotamers[ jj ][ ii ] = otfig->get_on_the_fly_node( ii )->get_rotamer( example_pro_rotamer_inds[ ii ], jj ).clone();
			}
			if ( example_gly_rotamer_inds[ ii ] != 0 ) {
				example_gly_rotamers[ jj ][ ii ] = otfig->get_on_the_fly_node( ii )->get_rotamer( example_gly_rotamer_inds[ ii ], jj ).clone();
			}
		}
	}

	/// 1. Iterate across all edges in the graph
	for ( Size ii = 1; ii <= nmoltenres(); ++ ii ) {
		Size const ii_resid = moltenres_2_resid( ii );
		if ( ! otfig->distinguish_backbone_and_sidechain_for_node( ii ) ) continue;
		for ( graph::Graph::EdgeListConstIter
				li  = packer_neighbor_graph->get_node( ii_resid )->const_edge_list_begin(),
				lie = packer_neighbor_graph->get_node( ii_resid )->const_edge_list_end();
				li != lie; ++li ) {
			Size jj_resid = (*li)->get_other_ind( ii_resid );
			Size const iijj_scale = Size(symm_info->score_multiply( ii_resid, jj_resid ));
			//std::cout << "pro_corr: ii_resid " << ii_resid << " jj_resid " << jj_resid << " scale: " << iijj_scale << std::endl;
			if ( iijj_scale == 0 ) continue;

			Size jj = resid_2_moltenres( jj_resid ); //pretend we're iterating over jj >= ii

			// if jj_resid isn't a molten residue, but it is part of the asymmetric unit, pass;
			// jj_resid represents a background residue.
			if ( jj == 0 && symm_info->chi_follows( jj_resid ) == 0 ) continue;

			// for residue pairs in the ASU, only consider pairs where ii < jj
			if ( symm_info->chi_follows( jj_resid ) == 0 && jj < ii ) continue;

			Size jj_resid_master = symm_info->chi_follows( jj_resid ) == 0 ? jj_resid : symm_info->chi_follows( jj_resid );
			Size jj_master = jj == 0 ? resid_2_moltenres( jj_resid_master ) : jj;
			if ( jj_master == 0 ) continue;
			if ( jj_master == ii ) continue;
			if ( ! otfig->distinguish_backbone_and_sidechain_for_node( jj_master ) ) continue;

			Size ii_master( ii ), ii_resid_master( ii_resid );
			Size ii_subunit = symm_info->subunit_index( ii_resid ); Size jj_subunit = symm_info->subunit_index( jj_resid );
			//std::cout << "proline correction; " << ii << " (" << ii_subunit << "), " << jj << " (" << jj_subunit << ") from residues " << ii_resid << " " << jj_resid << " with scale " << iijj_scale << std::endl;
			//bool swap( false );
			if ( jj_resid_master < ii_resid_master ) {
				Size temp;
				temp = ii_master; ii_master = jj_master; jj_master = temp;
				temp = ii_resid_master; ii_resid_master = jj_resid_master; jj_resid_master = temp;
				temp = ii_subunit; ii_subunit = jj_subunit; jj_subunit = temp;
				//swap = true;  // set but never used ~Labonte
			}

			/// 2. Calculate proline-correction terms between neighbors

			RotamerSetCOP ii_rotset = rotamer_set_for_moltenresidue( ii_master );
			RotamerSetCOP jj_rotset = rotamer_set_for_moltenresidue( jj_master );

			//std::cout << "SymmetricRotamerSets::compute_proline_corrections " << (*uli)->get_first_node_ind() << " " << (*uli)->get_second_node_ind() << " " << ii_master << " " << jj_master << " " << ii_subunit << " " << jj_subunit << std::endl;

			SymmOnTheFlyNode * iinode = otfig->get_on_the_fly_node( ii_master );
			SymmOnTheFlyNode * jjnode = otfig->get_on_the_fly_node( jj_master );

			for ( Size kk = 1; kk <= ii_rotset->num_rotamers(); ++kk ) {
				//bool target_interaction = ( ii_master == 1 && ( kk == 1 || kk == 9 ) );

				core::PackerEnergy bb_bbnonproE( iijj_scale ), bb_bbproE( iijj_scale );
				core::PackerEnergy sc_npbb_energy( iijj_scale ), sc_probb_energy( iijj_scale );
				//calc sc_npbb_energy;
				if ( example_gly_rotamer_inds[ jj_master ] != 0 ) {
					bb_bbnonproE   *= get_bb_bbE( pose, scfxn, iinode->get_rotamer( kk, ii_subunit ), *example_gly_rotamers[ jj_subunit ][ jj_master ] );
					sc_npbb_energy *= get_sc_bbE( pose, scfxn, iinode->get_rotamer( kk, ii_subunit ), *example_gly_rotamers[ jj_subunit ][ jj_master ] );
				} else {
					bb_bbnonproE = sc_npbb_energy = 0;
				}
				//calc sc_probb_energy
				if ( example_pro_rotamer_inds[ jj_master ] != 0 ) {
					bb_bbproE       *= get_bb_bbE( pose, scfxn, iinode->get_rotamer( kk, ii_subunit ), *example_pro_rotamers[ jj_subunit ][ jj_master ]  );
					sc_probb_energy *= get_sc_bbE( pose, scfxn, iinode->get_rotamer( kk, ii_subunit ), *example_pro_rotamers[ jj_subunit ][ jj_master ] );
				} else {
					bb_bbproE = sc_probb_energy = 0;
				}
				//if ( target_interaction ) { std::cout << "  procorr: " << ii_master << " " << jj_master << " 1 kk: " << kk << " " << bb_bbnonproE << " " << sc_npbb_energy << " " << bb_bbproE << " " << sc_probb_energy << " add1b: " << sc_npbb_energy +  0.5 * bb_bbnonproE << std::endl; }
				otfig->add_to_one_body_energy_for_node_state( ii_master, kk, sc_npbb_energy +  0.5 * bb_bbnonproE  );
				otfig->add_ProCorrection_values_for_edge( ii_master, jj_master, ii_master, kk,
					bb_bbnonproE, bb_bbproE, sc_npbb_energy, sc_probb_energy );
			}

			for ( Size kk = 1; kk <= jj_rotset->num_rotamers(); ++kk ) {
				// bool target_interaction = ( ii_master == 1 );
				core::PackerEnergy bb_bbnonproE( iijj_scale ), bb_bbproE( iijj_scale );
				core::PackerEnergy sc_npbb_energy( iijj_scale ), sc_probb_energy( iijj_scale );
				//calc sc_npbb_energy;
				if ( example_gly_rotamer_inds[ ii_master ] != 0 ) {
					bb_bbnonproE   *= get_bb_bbE( pose, scfxn, jjnode->get_rotamer( kk, jj_subunit ), *example_gly_rotamers[ ii_subunit ][ ii_master ] );
					sc_npbb_energy *= get_sc_bbE( pose, scfxn, jjnode->get_rotamer( kk, jj_subunit ), *example_gly_rotamers[ ii_subunit ][ ii_master ] );
				} else {
					bb_bbnonproE = sc_npbb_energy = 0;
				}
				//calc sc_probb_energy
				if ( example_pro_rotamer_inds[ ii_master ] != 0 ) {
					bb_bbproE       *= get_bb_bbE( pose, scfxn, jjnode->get_rotamer( kk, jj_subunit ), *example_pro_rotamers[ ii_subunit ][ ii_master ] );
					sc_probb_energy *= get_sc_bbE( pose, scfxn, jjnode->get_rotamer( kk, jj_subunit ), *example_pro_rotamers[ ii_subunit ][ ii_master ] );
				} else {
					bb_bbproE = sc_probb_energy = 0;
				}
				//if ( target_interaction ) { std::cout << "  procorr: " << ii_master << " " << jj_master << " 2 kk: " << kk << " " << bb_bbnonproE << " " << sc_npbb_energy << " " << bb_bbproE << " " << sc_probb_energy << " add1b: " << sc_npbb_energy +  0.5 * bb_bbnonproE << std::endl; }
				otfig->add_to_one_body_energy_for_node_state( jj_master, kk, sc_npbb_energy + 0.5 * bb_bbnonproE );

				otfig->add_ProCorrection_values_for_edge( ii_master, jj_master, jj_master, kk,
					bb_bbnonproE, bb_bbproE, sc_npbb_energy, sc_probb_energy );
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
	for ( Rotamers::const_iterator
			rot     = rotset_in->begin(),
			rot_end = rotset_in->end();
			rot != rot_end; ++rot ) {
		conformation::ResidueOP target_rsd( (*rot)->clone() /*pose.residue( seqpos )*/ );

		// peptoids have a different orientation function due to the placement of the first side chain atom
		if ( target_rsd->type().is_peptoid() ) {
			target_rsd->orient_onto_residue_peptoid( pose.residue( sympos ), pose.conformation() );
		} else {
			target_rsd->orient_onto_residue( pose.residue( sympos ) );
		}

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
			clone != clone_end; ++clone ) {
		if (  *clone > jj_highest &&
				packer_neighbor_graph->get_edge_exists(ii_resid, *clone) ) {
			jj_highest = *clone;
		}
	}
	// find the highest jj clone sequence number we could have
	uint ii_highest(ii_resid);
	for ( std::vector< Size>::const_iterator
			clone     = symm_info->chi_clones( ii_resid ).begin(),
			clone_end = symm_info->chi_clones( ii_resid ).end();
			clone != clone_end; ++clone ) {
		if (  *clone > ii_highest &&
				packer_neighbor_graph->get_edge_exists(*clone, jj_resid) ) {
			ii_highest = *clone;
		}
	}
	if ( ii_resid == ii_highest && jj_resid == jj_highest ) return true;
	else return false;
}

} // namespace symmetry
} // namespace rotamer_set
} // namespace pack
} // namespace core

