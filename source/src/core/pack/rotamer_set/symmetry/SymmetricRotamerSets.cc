// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/RotamerSet/RotamerSets.cc
/// @brief  RotamerSets class implementation, for symmetric packing
/// @author Ingemar Andre

// Unit Headers
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>
#include <core/pack/interaction_graph/SymmOnTheFlyInteractionGraph.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/MirrorSymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <utility/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>

#include <basic/thread_manager/RosettaThreadManager.hh>
#include <basic/thread_manager/RosettaThreadAssignmentInfo.hh>


// C++

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>

static basic::Tracer TR( "core.pack.rotamer_set.symmetry.SymmetricRotamerSets", basic::t_info );

using namespace ObjexxFCL;

// The multiplier cutoff below which we ignore interactions.
#define SYMM_INFO_SCORE_MULTIPLY_CUTOFF 1e-3

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace rotamer_set {
namespace symmetry {

SymmetricRotamerSets::SymmetricRotamerSets() = default;
SymmetricRotamerSets::~SymmetricRotamerSets() = default;

// @details For now, this function knows that all RPEs are computed before annealing begins
// In not very long, this function will understand that there are different ways to compute
// energies and to store them, but for now, it does not.
// symmetrical version of RotamerSets::compute_energies.
void
SymmetricRotamerSets::compute_energies(
	pose::Pose const & pose,
	scoring::ScoreFunction const & scfxn,
	utility::graph::GraphCOP packer_neighbor_graph,
	interaction_graph::InteractionGraphBaseOP ig,
	core::Size const threads_to_request
) {
	using namespace interaction_graph;
	using namespace scoring;

	ig->initialize( *this );
	//We create a vector of computations to do.
	basic::thread_manager::RosettaThreadAssignmentInfo thread_assignment_info( basic::thread_manager::RosettaThreadRequestOriginatingLevel::CORE_PACK );
	utility::vector1< basic::thread_manager::RosettaThreadFunction > work_vector;
	work_vector = construct_one_body_energy_work_vector( pose, scfxn, packer_neighbor_graph, ig, thread_assignment_info );

	PrecomputedPairEnergiesInteractionGraphOP pig =
		utility::pointer::dynamic_pointer_cast< core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraph > ( ig );
	if ( pig ) {
		//We append to the vector of computations to do.
		append_two_body_energy_computations_to_work_vector( pose, scfxn, packer_neighbor_graph, pig, work_vector, thread_assignment_info );
		basic::thread_manager::RosettaThreadManager::get_instance()->do_work_vector_in_threads( work_vector, threads_to_request, thread_assignment_info );

		//In a single thread, finalize the edges (not threadsafe):
		pig->declare_all_edge_energies_final();
	} else {
		SymmOnTheFlyInteractionGraphOP symotfig =
			utility::pointer::dynamic_pointer_cast< core::pack::interaction_graph::SymmOnTheFlyInteractionGraph > ( ig );
		if ( symotfig ) {
			if ( (threads_to_request == 0 || threads_to_request > 1) && basic::thread_manager::RosettaThreadManager::total_threads() > 1 ) {
				TR.Warning << "Cannot pre-compute an on-the-fly interaction graph in threads.  Only precomputing single-body energies in threads." << std::endl;
			}
			basic::thread_manager::RosettaThreadManager::get_instance()->do_work_vector_in_threads( work_vector, threads_to_request, thread_assignment_info );
			prepare_symm_otf_interaction_graph( pose, scfxn, packer_neighbor_graph, symotfig );
		} else {
			utility_exit_with_message( "Encountered incompatible interaction graph type in SymmetricRotamerSets::compute_energies" );
		}
	}

#ifdef MULTI_THREADED
	TR << "Completed symmetric interaction graph pre-calculation in " << thread_assignment_info.get_assigned_total_thread_count() << " available threads (" << (thread_assignment_info.get_requested_thread_count() == 0 ? std::string( "all available" ) : std::to_string( thread_assignment_info.get_requested_thread_count() ) ) << " had been requested)." << std::endl;
#endif
}

/// @brief Append to a vector of calculations that already contains one-body energy calculations additonal work units
/// for two-body energy calculations, for subsequent multi-threaded evaluation.
/// @details Each individual calculation is for the interaction of all possible rotamer pairs at two positions.  Again,
/// not as granular as conceivably possible, but easier to set up.
/// @note The work_vector vector is extended by this operation.  This version is for the symmetric case.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitue.org).
void
SymmetricRotamerSets::append_two_body_energy_computations_to_work_vector(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & scfxn,
	utility::graph::GraphCOP packer_neighbor_graph,
	interaction_graph::PrecomputedPairEnergiesInteractionGraphOP pig,
	utility::vector1< basic::thread_manager::RosettaThreadFunction > & work_vector,
	basic::thread_manager::RosettaThreadAssignmentInfo const & thread_assignment_info
) const {
	// find SymmInfo
	auto const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );
	runtime_assert( symm_info != nullptr );

	// Two body energies
	for ( core::Size ii(1); ii <= nmoltenres(); ++ii ) {
		core::Size const ii_resid( moltenres_2_resid( ii ) );
		// We will loop over only one subunit here.
		for ( utility::graph::Graph::EdgeListConstIter
				uli  = packer_neighbor_graph->get_node( ii_resid )->const_edge_list_begin(),
				ulie = packer_neighbor_graph->get_node( ii_resid )->const_edge_list_end();
				uli != ulie; ++uli ) {
			core::Size jj_resid = (*uli)->get_other_ind( ii_resid );
			core::Size jj = resid_2_moltenres( jj_resid ); //pretend we're iterating over jj >= ii
			// if jj_resid is not repackable and jj_resid is not in a different subunit continue
			core::Size jj_resid_master(0);
			if ( symm_info->chi_follows( jj_resid ) == 0  ) {
				// jj_resid is independent like ii_resid
				if ( jj == 0 ) continue;
				if ( jj_resid <= ii_resid ) continue; // we will hit this in the other order
				jj_resid_master = jj_resid;
			} else {
				// if jj_resid is in a different subunit we need to know its master
				jj_resid_master = symm_info->chi_follows( jj_resid );
			}
			// Find out the moltenres id for the master:
			core::Size jj_master = resid_2_moltenres( jj_resid_master );
			// If the master is not repackable continue:
			if ( jj_master == 0 ) continue;
			// If we the ii_resid and jj_resid_master have an edge this intersubunit.

			if ( fabs( symm_info->score_multiply(ii_resid,jj_resid) ) <  SYMM_INFO_SCORE_MULTIPLY_CUTOFF ) continue;
			core::Size ii_master( ii );

			// The self interaction energy is already calculated in the one-body term,
			if ( ii_master == jj_master ) continue;
			// Swap the order of the nodes that form the edge if jj_resid_master < ii_resid.
			// Edges are always stored as pairs (a, b) where a > b. This only happens if we
			// are calculating interactions from the controling subunit to another subunit.
			bool swap( false );
			if ( jj_resid_master < ii_resid ) {
				swap = true;
				core::Size temp_ii = ii_master;
				ii_master = jj_master;
				jj_master = temp_ii;
			}

			RotamerSetCOP ii_rotset( rotamer_set_for_moltenresidue( ii_master ) );
			RotamerSetCOP jj_rotset( rotamer_set_for_moltenresidue( jj_master ) );

			// Add the edge to the interaction graph:
			if ( !pig->get_edge_exists( ii_master, jj_master ) ) {
				pig->add_edge( ii_master, jj_master );
			}

			//NOTE: The behaviour of std::bind is to pass EVERYTHING by copy UNLESS std::cref or std::ref is used.  So in
			//the following, all of the std::pairs are passed by copy.
			work_vector.push_back(
				std::bind(
				&interaction_graph::PrecomputedPairEnergiesInteractionGraph::set_twobody_energies_multithreaded, pig.get(),
				std::make_pair( ii_master, jj_master ), std::make_pair( ii_rotset, jj_rotset ), std::cref(pose), std::cref(scfxn),
				std::cref(thread_assignment_info), false, symm_info, std::make_pair( ii_resid, jj_resid ),
				symm_info->chi_follows( jj_resid ) != 0 && jj == 0, swap
				)
			);
		}
	}

	// Iterate across the long range energy functions and use the iterators generated
	// by the LRnergy container object
	// The logic here is exactly the same as for the non-lr energy calculation in terms of symmetry
	for ( auto
			lr_iter( scfxn.long_range_energies_begin() ),
			lr_end( scfxn.long_range_energies_end() );
			lr_iter != lr_end; ++lr_iter ) {
		core::scoring::LREnergyContainerCOP lrec( pose.energies().long_range_container( (*lr_iter)->long_range_type() ) );
		if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.

		for ( core::Size ii = 1; ii <= nmoltenres(); ++ ii ) {
			core::Size const ii_resid = moltenres_2_resid( ii );

			for ( core::scoring::ResidueNeighborConstIteratorOP
					rni = lrec->const_neighbor_iterator_begin( ii_resid ),
					rniend = lrec->const_neighbor_iterator_end( ii_resid );
					(*rni) != (*rniend); ++(*rni) ) {
				core::Size jj_resid = ( rni->upper_neighbor_id() == ii_resid ? rni->lower_neighbor_id() : rni->upper_neighbor_id() );

				core::Size jj = resid_2_moltenres( jj_resid ); //pretend we're iterating over jj >= ii

				// if jj_resid is not repackable and jj_resid is not in a different subunit continue
				core::Size jj_resid_master(0);
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
				core::Size jj_master = resid_2_moltenres( jj_resid_master );
				// if the master is not repackable continue
				if ( jj_master == 0 ) continue;
				// if the ii_resid and jj_resid_master have an edge this intersubunit.

				if ( fabs( symm_info->score_multiply(ii_resid,jj_resid)) < SYMM_INFO_SCORE_MULTIPLY_CUTOFF ) continue;

				core::Size ii_master = ii;

				// the self interaction energy is already calculated in the one-body ter,
				if ( ii_master == jj_master ) continue;
				// swap the order of the nodes that form the edge if jj_resid_master < ii_resid.
				// Edges are always stored as pairs (a, b) where a > b. This only happens if we
				// are calculation interactions from the controling subunit to another subunit
				bool swap( false );
				if ( jj_resid_master < ii_resid ) {
					swap = true;
					core::Size temp_ii = ii_master;
					ii_master = jj_master;
					jj_master = temp_ii;
				}

				RotamerSetCOP ii_rotset = rotamer_set_for_moltenresidue( ii_master );
				RotamerSetCOP jj_rotset = rotamer_set_for_moltenresidue( jj_master );

				if ( ! pig->get_edge_exists( ii_master, jj_master ) ) {
					pig->add_edge( ii_master, jj_master );
				}

				//NOTE: The behaviour of std::bind is to pass EVERYTHING by copy UNLESS std::cref or std::ref is used.  So in
				//the following, all of the std::pairs are passed by copy.
				work_vector.push_back(
					std::bind(
					&interaction_graph::PrecomputedPairEnergiesInteractionGraph::add_longrange_twobody_energies_multithreaded, pig.get(),
					std::make_pair( ii_master, jj_master ), std::make_pair( ii_rotset, jj_rotset ), std::cref(pose), *lr_iter, std::cref(scfxn),
					std::cref(thread_assignment_info), false, symm_info, std::make_pair( ii_resid, jj_resid ),
					symm_info->chi_follows( jj_resid ) != 0 && jj == 0, swap
					)
				);
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
	utility::graph::GraphCOP packer_neighbor_graph,
	interaction_graph::SymmOnTheFlyInteractionGraphOP ig
)
{
	// find SymmInfo
	auto const & SymmConf (
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
			//std::cout << "temp" << all_canonical_aas << std::endl;
		}
	}


	for ( Size ii = 1; ii <= nmoltenres(); ++ii ) {
		Size ii_resid = moltenres_2_resid( ii );
		//std::cout << "ii_resid: " << ii_resid << std::endl;
		for ( utility::graph::Graph::EdgeListConstIter
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
				//temp = ii_resid_master; ii_resid_master =jj_resid_master; jj_resid_master = temp; // swapped values are unused
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
	for ( auto
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
	utility::graph::GraphCOP packer_neighbor_graph,
	interaction_graph::SymmOnTheFlyInteractionGraphOP otfig
)
{
	using namespace conformation;
	using namespace conformation::symmetry;
	using namespace interaction_graph;

	Size const nsubunits = symm_info->subunits();

	utility::vector1< Size > example_reg_rotamer_inds( nmoltenres(), 0 );
	utility::vector1< Size > example_gly_rotamer_inds( nmoltenres(), 0 );
	utility::vector1< Size > example_pro_rotamer_inds( nmoltenres(), 0 );
	utility::vector1< utility::vector1< ResidueOP > > example_reg_rotamers( nsubunits );
	utility::vector1< utility::vector1< ResidueOP > > example_gly_rotamers( nsubunits );
	utility::vector1< utility::vector1< ResidueOP > > example_pro_rotamers( nsubunits );
	for ( Size ii = 1; ii <= nsubunits; ++ii ) {
		example_reg_rotamers[ ii ].resize( nmoltenres() );
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
				example_reg_rotamer_inds[ ii ] = ii_rotset->get_residue_type_begin( jj );
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
			if ( example_reg_rotamer_inds[ ii ] != 0 ) {
				example_reg_rotamers[ jj ][ ii ] = otfig->get_on_the_fly_node( ii )->get_rotamer( example_reg_rotamer_inds[ ii ], jj ).clone();
			}
		}
	}

	/// 1. Iterate across all edges in the graph
	for ( Size ii = 1; ii <= nmoltenres(); ++ ii ) {
		Size const ii_resid = moltenres_2_resid( ii );
		if ( ! otfig->distinguish_backbone_and_sidechain_for_node( ii ) ) continue;
		for ( utility::graph::Graph::EdgeListConstIter
				li  = packer_neighbor_graph->get_node( ii_resid )->const_edge_list_begin(),
				lie = packer_neighbor_graph->get_node( ii_resid )->const_edge_list_end();
				li != lie; ++li ) {
			Size jj_resid = (*li)->get_other_ind( ii_resid );
			auto const iijj_scale = Size(symm_info->score_multiply( ii_resid, jj_resid ));
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
				//temp = ii_resid_master; ii_resid_master = jj_resid_master; jj_resid_master = temp; // swapped values are never used
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

				core::PackerEnergy bb_bbregE( iijj_scale ), bb_bbproE( iijj_scale ), bb_bbglyE( iijj_scale );
				core::PackerEnergy sc_regbb_energy( iijj_scale ), sc_probb_energy( iijj_scale ), sc_glybb_energy( iijj_scale );
				//calc sc_regbb_energy;
				if ( example_reg_rotamer_inds[ jj_master ] != 0 ) {
					bb_bbregE       *= get_bb_bbE( pose, scfxn, iinode->get_rotamer( kk, ii_subunit ), *example_reg_rotamers[ jj_subunit ][ jj_master ] );
					sc_regbb_energy *= get_sc_bbE( pose, scfxn, iinode->get_rotamer( kk, ii_subunit ), *example_reg_rotamers[ jj_subunit ][ jj_master ] );
				} else {
					bb_bbregE = sc_regbb_energy = 0;
				}

				if ( example_gly_rotamer_inds[ jj_master ] != 0 ) {
					bb_bbglyE       *= get_bb_bbE( pose, scfxn, iinode->get_rotamer( kk, ii_subunit ), *example_gly_rotamers[ jj_subunit ][ jj_master ] );
					sc_glybb_energy *= get_sc_bbE( pose, scfxn, iinode->get_rotamer( kk, ii_subunit ), *example_gly_rotamers[ jj_subunit ][ jj_master ] );
				} else {
					bb_bbglyE = sc_glybb_energy = 0;
				}
				//calc sc_probb_energy
				if ( example_pro_rotamer_inds[ jj_master ] != 0 ) {
					bb_bbproE       *= get_bb_bbE( pose, scfxn, iinode->get_rotamer( kk, ii_subunit ), *example_pro_rotamers[ jj_subunit ][ jj_master ]  );
					sc_probb_energy *= get_sc_bbE( pose, scfxn, iinode->get_rotamer( kk, ii_subunit ), *example_pro_rotamers[ jj_subunit ][ jj_master ] );
				} else {
					bb_bbproE = sc_probb_energy = 0;
				}
				//if ( target_interaction ) { std::cout << "  procorr: " << ii_master << " " << jj_master << " 1 kk: " << kk << " " << bb_bbnonproE << " " << sc_npbb_energy << " " << bb_bbproE << " " << sc_probb_energy << " add1b: " << sc_npbb_energy +  0.5 * bb_bbnonproE << std::endl; }

				otfig->add_to_one_body_energy_for_node_state( ii_master, kk, sc_regbb_energy +  0.5 * bb_bbregE  );
				otfig->add_ProCorrection_values_for_edge( ii_master, jj_master, ii_master, kk,
					bb_bbregE, bb_bbproE, sc_regbb_energy, sc_probb_energy );
				otfig->add_GlyCorrection_values_for_edge( ii_master, jj_master, ii_master, kk,
					bb_bbregE, bb_bbglyE, sc_regbb_energy, sc_glybb_energy );
			}

			for ( Size kk = 1; kk <= jj_rotset->num_rotamers(); ++kk ) {
				// bool target_interaction = ( ii_master == 1 );
				core::PackerEnergy bb_bbregE( iijj_scale ), bb_bbproE( iijj_scale ), bb_bbglyE( iijj_scale );
				core::PackerEnergy sc_regbb_energy( iijj_scale ), sc_probb_energy( iijj_scale ), sc_glybb_energy( iijj_scale );
				//calc sc_regbb_energy;
				if ( example_reg_rotamer_inds[ ii_master ] != 0 ) {
					bb_bbregE       *= get_bb_bbE( pose, scfxn, jjnode->get_rotamer( kk, jj_subunit ), *example_reg_rotamers[ ii_subunit ][ ii_master ] );
					sc_regbb_energy *= get_sc_bbE( pose, scfxn, jjnode->get_rotamer( kk, jj_subunit ), *example_reg_rotamers[ ii_subunit ][ ii_master ] );
				} else {
					bb_bbregE = sc_regbb_energy = 0;
				}

				//calc sc_glybb_energy;
				if ( example_gly_rotamer_inds[ ii_master ] != 0 ) {
					bb_bbglyE       *= get_bb_bbE( pose, scfxn, jjnode->get_rotamer( kk, jj_subunit ), *example_gly_rotamers[ ii_subunit ][ ii_master ] );
					sc_glybb_energy *= get_sc_bbE( pose, scfxn, jjnode->get_rotamer( kk, jj_subunit ), *example_gly_rotamers[ ii_subunit ][ ii_master ] );
				} else {
					bb_bbglyE = sc_glybb_energy = 0;
				}

				//calc sc_probb_energy
				if ( example_pro_rotamer_inds[ ii_master ] != 0 ) {
					bb_bbproE       *= get_bb_bbE( pose, scfxn, jjnode->get_rotamer( kk, jj_subunit ), *example_pro_rotamers[ ii_subunit ][ ii_master ] );
					sc_probb_energy *= get_sc_bbE( pose, scfxn, jjnode->get_rotamer( kk, jj_subunit ), *example_pro_rotamers[ ii_subunit ][ ii_master ] );
				} else {
					bb_bbproE = sc_probb_energy = 0;
				}
				//if ( target_interaction ) { std::cout << "  procorr: " << ii_master << " " << jj_master << " 2 kk: " << kk << " " << bb_bbnonproE << " " << sc_npbb_energy << " " << bb_bbproE << " " << sc_probb_energy << " add1b: " << sc_npbb_energy +  0.5 * bb_bbnonproE << std::endl; }
				otfig->add_to_one_body_energy_for_node_state( jj_master, kk, sc_regbb_energy + 0.5 * bb_bbregE );

				otfig->add_ProCorrection_values_for_edge( ii_master, jj_master, jj_master, kk,
					bb_bbregE, bb_bbproE, sc_regbb_energy, sc_probb_energy );
				otfig->add_GlyCorrection_values_for_edge( ii_master, jj_master, jj_master, kk,
					bb_bbregE, bb_bbglyE, sc_regbb_energy, sc_glybb_energy );
			}
		}
	}

}

/// @brief Generate a new rotamer set oriented onto a reference rotamer set by cloning the reference set and reorienting
/// using symmetry information.  Return an owning pointer to the new, reoriented rotamer set.
/// @details Orients all rotamers in a rotamer_set to a different (symmetrical) position
/// @note If set_up_mirror_symmetry_if_has_mirror_symmetry_ is true, then residues in mirrored subunits have their
/// ResidueTypes switched to the types of the opposite chirality.  If false (the default), then they keep the same
/// types, and only have their geometry mirrored.  The default is false, which is computationally less expensive
/// at the expense of having the incorrect types in mirrored subunits.
RotamerSetOP
SymmetricRotamerSets::orient_rotamer_set_to_symmetric_partner(
	pose::Pose const & pose,
	RotamerSetCOP rotset_for_seqpos,
	uint const & sympos,
	bool const set_up_mirror_types_if_has_mirror_symmetry/*=false*/
) {

	RotamerSetOP sym_rotamer_set = RotamerSetFactory::create_rotamer_set( pose );
	sym_rotamer_set->set_resid(sympos);

	auto const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation() ) );

	core::conformation::symmetry::MirrorSymmetricConformationCOP mirrorconf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::MirrorSymmetricConformation const>( pose.conformation_ptr() ) );

	for ( const auto & rot : *rotset_for_seqpos ) {
		bool mirrored(false);
		if ( mirrorconf ) { //If this is a mirror symmetric conformation, figure out whether this subunit is mirrored:
			mirrored = mirrorconf->res_is_mirrored( sympos );
		}

		// fpd  the commented line turns on chirality flipping in the packer
		// fpd  right now, no scorefunctions need it however and it is _very_ slow (triples time of design)
		// VKM -- I've left it off, but added a function interface to turn it back on.  In the special case the the interaction graph for non-
		// pairwise score terms, it needs to be on, but for the default interaction graph, it's probably a waste of clock cycles.  So the
		// ResidueArrayAnnealingEvaluator calls this with set_up_mirror_types-if_has_mirror_symmetry=true, and everything else calls it with
		// set_up_mirror_types_if_has_mirror_symmetry=false.
		conformation::ResidueOP target_rsd;
		if ( set_up_mirror_types_if_has_mirror_symmetry && mirrored ) {
			target_rsd = rot->clone_flipping_chirality( *pose.residue_type_set_for_pose( rot->type().mode() ) );
		} else {
			target_rsd = rot->clone();
		}

		// peptoids have a different orientation function due to the placement of the first side chain atom
		//if ( target_rsd->type().is_peptoid() ) {
		// target_rsd->orient_onto_residue_peptoid( pose.residue( sympos ), pose.conformation() );
		//} else {
		// target_rsd->orient_onto_residue( pose.residue( sympos ) );
		//}

		//fpd use the stored transforms instead!  Will work for anything!
		for ( int i=1; i<=(int)target_rsd->natoms(); ++i ) {
			target_rsd->set_xyz(i ,
				SymmConf.apply_transformation_norecompute( target_rsd->xyz(i), target_rsd->seqpos(), sympos ) );
		}

		sym_rotamer_set->add_rotamer( *target_rsd );
	}
	return sym_rotamer_set;
}

// @details is this the final visit to an edge?
bool
SymmetricRotamerSets::final_visit_to_edge(
	pose::Pose const & pose,
	utility::graph::GraphCOP packer_neighbor_graph,
	uint ii_resid,
	uint jj_resid
)
{
	// find SymmInfo
	auto const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	// find the highest jj clone sequence number we could have
	uint jj_highest(jj_resid);
	for ( unsigned long clone : symm_info->chi_clones( jj_resid ) ) {
		if (  clone > jj_highest &&
				packer_neighbor_graph->get_edge_exists(ii_resid, clone) ) {
			jj_highest = clone;
		}
	}
	// find the highest jj clone sequence number we could have
	uint ii_highest(ii_resid);
	for ( unsigned long clone : symm_info->chi_clones( ii_resid ) ) {
		if (  clone > ii_highest &&
				packer_neighbor_graph->get_edge_exists(clone, jj_resid) ) {
			ii_highest = clone;
		}
	}
	if ( ii_resid == ii_highest && jj_resid == jj_highest ) return true;
	else return false;
}

//fpd function to set some pose data needed SymmetricRotamerSets
void
SymmetricRotamerSets::initialize_pose_for_rotsets_creation(
	pose::Pose & pose
) const {
	auto & SymmConf (
		dynamic_cast<SymmetricConformation &> ( pose.conformation() ) );
	SymmConf.recalculate_transforms();
}


} // namespace symmetry
} // namespace rotamer_set
} // namespace pack
} // namespace core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::rotamer_set::symmetry::SymmetricRotamerSets::save( Archive & arc ) const {
	arc( cereal::base_class< core::pack::rotamer_set::RotamerSets >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::rotamer_set::symmetry::SymmetricRotamerSets::load( Archive & arc ) {
	arc( cereal::base_class< core::pack::rotamer_set::RotamerSets >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::rotamer_set::symmetry::SymmetricRotamerSets );
CEREAL_REGISTER_TYPE( core::pack::rotamer_set::symmetry::SymmetricRotamerSets )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_rotamer_set_symmetry_SymmetricRotamerSets )
#endif // SERIALIZATION
