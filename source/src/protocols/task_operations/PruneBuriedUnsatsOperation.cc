// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/PruneBuriedUnsatsOperation.cc
/// @brief  Eliminate rotamers that contain a sidechain buried unsat and make no other sidechain h-bonds
/// @detail This is intended to speed up packing by eliminating useless polar rotamers. This can
///         also help to reduce the number of buried unsats in designs because Rosetta can't pack them.
/// @author Longxing Cao -- original idea
/// @author Brian Coventry ( bcov@uw.edu ) -- code implementation


// Unit Headers
#include <protocols/task_operations/PruneBuriedUnsatsOperation.hh>
#include <protocols/task_operations/PruneBuriedUnsatsOperationCreator.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pack/hbonds/HBondGraph_util.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/scoring/atomic_depth/AtomicDepth.hh>
#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

static basic::Tracer TR( "protocols.task_operations.PruneBuriedUnsatsOperation" );

namespace protocols {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

PruneBuriedUnsats_RotamerSetsOperation::PruneBuriedUnsats_RotamerSetsOperation(
	bool allow_even_trades,
	core::Real atomic_depth_probe_radius,
	core::Real atomic_depth_resolution,
	core::Real atomic_depth_cutoff,
	core::Real minimum_hbond_energy,
	core::scoring::ScoreFunctionOP scorefxn_sc,
	core::scoring::ScoreFunctionOP scorefxn_bb
) :
	core::pack::rotamer_set::RotamerSetsOperation(),
	allow_even_trades_( allow_even_trades ),
	atomic_depth_probe_radius_( atomic_depth_probe_radius ),
	atomic_depth_resolution_( atomic_depth_resolution ),
	atomic_depth_cutoff_( atomic_depth_cutoff ),
	minimum_hbond_energy_( minimum_hbond_energy ),
	scorefxn_sc_( scorefxn_sc ),
	scorefxn_bb_( scorefxn_bb )
{}

PruneBuriedUnsats_RotamerSetsOperation::~PruneBuriedUnsats_RotamerSetsOperation() = default;

core::pack::rotamer_set::RotamerSetsOperationOP
PruneBuriedUnsats_RotamerSetsOperation::clone() const
{
	return core::pack::rotamer_set::RotamerSetsOperationOP( new PruneBuriedUnsats_RotamerSetsOperation( *this ) );
}


// Steps
//
// 1. Enumerate all atoms and hbonds
// 2. Identify buried atoms
// 2.5. Prune buried atoms satisfied by backbone atoms (or unpackable rotamers)
// 3. Identify all atoms that hbond to a buried atom
// 4. Identify unsatisfiable atoms
// 5. Sum unsatisfiable atoms and atoms that hbond per rotamer
// 6. Prune anything negative





void
PruneBuriedUnsats_RotamerSetsOperation::alter_rotamer_sets(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const &,
	core::pack::task::PackerTask const &,
	utility::graph::GraphCOP,
	core::pack::rotamer_set::RotamerSets & rotamer_sets_in
)
{

	using namespace core;

	Size starting_rotamers = rotamer_sets_in.nrotamers();

	// 1. Enumerate all atoms and hbonds

	utility::vector1<bool> position_had_rotset;
	pack::rotamer_set::RotamerSetsOP complete_rotsets = nullptr;

	TR << "Building hbond graph " << std::endl;

	scoring::hbonds::graph::AtomLevelHBondGraphOP hb_graph;
	hb_graph = core::pack::hbonds::hbond_graph_from_partial_rotsets( pose, rotamer_sets_in, scorefxn_sc_, scorefxn_bb_,
		complete_rotsets, position_had_rotset, minimum_hbond_energy_ );

	TR << "Hbond graph has: " << hb_graph->num_edges() << " edges requiring: " << hb_graph->getTotalMemoryUsage() << " bytes" << std::endl;


	// 2. Identify buried atoms

	// Not quite atom ids as the residue is actually the rotamer
	core::id::AtomID_Map<bool> is_buried( complete_rotsets->nrotamers(), false );
	core::id::AtomID_Map<bool> hbonds_to_buried( complete_rotsets->nrotamers(), false );
	core::id::AtomID_Map<bool> is_satisfied( complete_rotsets->nrotamers(), false );

	scoring::atomic_depth::AtomicDepth atomic_depth( pose, atomic_depth_probe_radius_, true, atomic_depth_resolution_ );
	chemical::AtomTypeSet const & atom_type_set = pose.residue(1).type().atom_type_set();


	for ( Size ihbnode = 1; ihbnode <= hb_graph->num_nodes(); ihbnode++ ) {

		scoring::hbonds::graph::AtomLevelHBondNode * hbnode = hb_graph->get_node( ihbnode );
		runtime_assert( hbnode );
		conformation::ResidueCOP rotamer = complete_rotsets->rotamer( ihbnode );

		is_buried.resize( ihbnode, rotamer->nheavyatoms() );
		hbonds_to_buried.resize( ihbnode, rotamer->nheavyatoms() ); // setting this guy up for later
		is_satisfied.resize( ihbnode, rotamer->nheavyatoms() ); // setting this guy up for later

		utility::vector1< scoring::hbonds::graph::AtomInfo > const & hbnode_atoms =
			hbnode->polar_sc_atoms_not_satisfied_by_background();

		for ( scoring::hbonds::graph::AtomInfo const & atom_info : hbnode_atoms ) {
			if ( atom_info.is_hydrogen() ) continue;

			Real depth = atomic_depth.calcdepth( rotamer->atom( atom_info.local_atom_id() ), atom_type_set );
			if ( depth < atomic_depth_cutoff_ ) continue;

			is_buried( ihbnode, atom_info.local_atom_id() ) = true;
			// std::cout << "BUR: " << ihbnode << " at: " << atom_info.local_atom_id() << std::endl;
		}
	}


	// 2.5. Prune buried atoms satisfied by backbone atoms (or unpackable rotamers)

	for ( Size ihbnode = 1; ihbnode <= hb_graph->num_nodes(); ihbnode++ ) {

		scoring::hbonds::graph::AtomLevelHBondNode * hbnode = hb_graph->get_node( ihbnode );
		runtime_assert( hbnode );
		conformation::ResidueCOP rotamer = complete_rotsets->rotamer( ihbnode );
		Size node_resnum = complete_rotsets->res_for_rotamer( ihbnode );

		// Do this in a smart way. Only look at nodes that might have backbone or unpackable atoms
		//  i.e. The first rotamer at every position -- backbone atoms and potentially unpackable
		//         Then ! position_had_rotset

		// BB is only on first moltenres. Unpackable positions only have 1 residue
		if ( complete_rotsets->rotid_on_moltenresidue( ihbnode ) != 1 ) continue;

		bool all_atoms_unpackable = ! position_had_rotset[ node_resnum ];


		for ( utility::graph::LowMemEdgeListIter it = hbnode->edge_list_begin( *hb_graph );
				it != hbnode->edge_list_end( *hb_graph );
				++it ) {

			scoring::hbonds::graph::AtomLevelHBondEdge * hb_edge =
				static_cast< scoring::hbonds::graph::AtomLevelHBondEdge * >( *it );
			Size first_node = hb_edge->get_first_node_ind();
			Size second_node = hb_edge->get_second_node_ind();

			bool we_are_the_first_node = first_node == ihbnode;
			Size other_node = we_are_the_first_node ? second_node : first_node;

			for ( scoring::hbonds::graph::HBondInfo const & hb_info : hb_edge->hbonds() ) {
				//          first is donor   0        1
				//  we are first   0        don      acc
				//                 1        acc      don
				//                                           xor
				bool we_are_the_donor = ! ( we_are_the_first_node ^ hb_info.first_node_is_donor() );

				Size our_heavy_atom = we_are_the_donor ? hb_info.local_atom_id_D() : hb_info.local_atom_id_A();

				if ( ! all_atoms_unpackable ) {
					if ( ! rotamer->atom_is_backbone( our_heavy_atom ) ) continue;
				}

				// We are now looking at an hbond where our side is unpackable

				Size other_atom = we_are_the_donor ? hb_info.local_atom_id_A() : hb_info.local_atom_id_D();
				is_buried( other_node, other_atom ) = false;


				// std::cout << "UNBUR: " << other_node << " at: " << other_atom << std::endl;
			}
		}
	}

	// 3. Identify all atoms that hbond to a buried atom
	// 4. Identify unsatisfiable atoms

	for ( Size ihbnode = 1; ihbnode <= hb_graph->num_nodes(); ihbnode++ ) {

		scoring::hbonds::graph::AtomLevelHBondNode * hbnode = hb_graph->get_node( ihbnode );
		runtime_assert( hbnode );
		conformation::ResidueCOP rotamer = complete_rotsets->rotamer( ihbnode );

		for ( utility::graph::LowMemEdgeListIter it = hbnode->edge_list_begin( *hb_graph );
				it != hbnode->edge_list_end( *hb_graph );
				++it ) {

			scoring::hbonds::graph::AtomLevelHBondEdge * hb_edge =
				static_cast< scoring::hbonds::graph::AtomLevelHBondEdge * >( *it );
			Size first_node = hb_edge->get_first_node_ind();
			Size second_node = hb_edge->get_second_node_ind();

			bool we_are_the_first_node = first_node == ihbnode;
			if ( ! we_are_the_first_node ) continue;   // Don't look at every hbond twice!

			Size other_node = we_are_the_first_node ? second_node : first_node;

			for ( scoring::hbonds::graph::HBondInfo const & hb_info : hb_edge->hbonds() ) {
				//          first is donor   0        1
				//  we are first   0        don      acc
				//                 1        acc      don
				//                                           xor
				bool we_are_the_donor = ! ( we_are_the_first_node ^ hb_info.first_node_is_donor() );

				Size our_heavy_atom = we_are_the_donor ? hb_info.local_atom_id_D() : hb_info.local_atom_id_A();

				Size other_atom = we_are_the_donor ? hb_info.local_atom_id_A() : hb_info.local_atom_id_D();

				bool we_are_buried = is_buried( ihbnode, our_heavy_atom );
				bool other_is_buried = is_buried( other_node, other_atom );

				// std::cout << "HB: " << ihbnode << " at: " << our_heavy_atom << " " << other_node << " at: " << other_atom << std::endl;

				if ( we_are_buried ) {
					is_satisfied( ihbnode, our_heavy_atom ) = true;
					hbonds_to_buried( other_node, other_atom ) = true;
					// std::cout << "SAT: is_sat: " << ihbnode << " at: " << our_heavy_atom << " by: " << other_node << " at: " << other_atom << std::endl;
				}
				if ( other_is_buried ) {
					is_satisfied( other_node, other_atom ) = true;
					hbonds_to_buried( ihbnode, our_heavy_atom ) = true;
					// std::cout << "SAT: is_sat: " << other_node << " at: " << other_atom << " by: " << ihbnode << " at: " << our_heavy_atom << std::endl;
				}
			}
		}
	}


	// 5. Sum unsatisfiable atoms and atoms that hbond per rotamer (excluding backbone)

	utility::vector1< utility::vector1<int> > resnum_rotamer_score( pose.size() );

	int even_trade_penalty = allow_even_trades_ ? 0 : -1;

	for ( Size ihbnode = 1; ihbnode <= hb_graph->num_nodes(); ihbnode++ ) {

		scoring::hbonds::graph::AtomLevelHBondNode * hbnode = hb_graph->get_node( ihbnode );
		runtime_assert( hbnode );
		conformation::ResidueCOP rotamer = complete_rotsets->rotamer( ihbnode );
		Size node_resnum = complete_rotsets->res_for_rotamer( ihbnode );
		Size moltenres = complete_rotsets->moltenres_for_rotamer( ihbnode );
		Size rotid_on_moltenres = complete_rotsets->rotid_on_moltenresidue( ihbnode );

		utility::vector1<int> & rotamer_scores = resnum_rotamer_score[ node_resnum ];
		if ( rotamer_scores.size() == 0 ) rotamer_scores.resize( complete_rotsets->nrotamers_for_moltenres( moltenres ) );

		bool penalty_applied = false;
		int satisfying_minus_unsatisfied = 0;
		for ( Size iatom = 1; iatom <= hbonds_to_buried.n_atom( ihbnode ); iatom++ ) {
			if ( rotamer->atom_is_backbone( iatom ) ) continue;
			if ( hbonds_to_buried( ihbnode, iatom ) ) {
				satisfying_minus_unsatisfied++;
			}
			if ( is_buried( ihbnode, iatom ) && ! is_satisfied( ihbnode, iatom ) ) {
				satisfying_minus_unsatisfied--;
				if ( ! penalty_applied ) {
					satisfying_minus_unsatisfied += even_trade_penalty;
					penalty_applied = true;
				}
			}
		}
		// std::cout << "SCORE: " << ihbnode << " " << satisfying_minus_unsatisfied << std::endl;
		rotamer_scores[ rotid_on_moltenres ] = satisfying_minus_unsatisfied;
	}

	// 6. Prune anything negative (unless everything is negative in which case take the least negative)
	for ( Size resnum = 1; resnum <= resnum_rotamer_score.size(); resnum++ ) {
		if ( ! position_had_rotset[ resnum ] ) continue;

		utility::vector1<int> const & rotamer_scores = resnum_rotamer_score[ resnum ];

		auto min_max_element = std::minmax_element( rotamer_scores.begin(), rotamer_scores.end() );
		int min_score = *(min_max_element.first);
		int max_score = *(min_max_element.second);

		// std::cout << "SCORES: resnum: " << resnum << " low: " << min_score << " max: " << max_score << std::endl;

		if ( min_score >= 0 ) continue;

		// If max score < 0, we keep everything that has max_score
		int prune_if_worse_than = std::min<int>( 0, max_score );

		// Shortcut the situation where nothing will be pruned
		if ( min_score == prune_if_worse_than ) continue;

		utility::vector1<bool> to_prune( rotamer_scores.size(), false );

		for ( Size i = 1; i <= to_prune.size(); i++ ) {
			to_prune[i] = rotamer_scores[i] < prune_if_worse_than;
		}

		rotamer_sets_in.rotamer_set_for_residue( resnum )->drop_rotamers( to_prune );
	}

	rotamer_sets_in.update_offset_data();

	Size ending_rotamers = rotamer_sets_in.nrotamers();

	TR << "Pruned " << starting_rotamers - ending_rotamers << " / " << starting_rotamers << " rotamers." << std::endl;


}

///////////////////////////////////////////////////////////////////////////////////////////
core::pack::task::operation::TaskOperationOP
PruneBuriedUnsatsOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new PruneBuriedUnsatsOperation );
}

void PruneBuriedUnsatsOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PruneBuriedUnsatsOperation::provide_xml_schema( xsd );
}

std::string PruneBuriedUnsatsOperationCreator::keyname() const
{
	return PruneBuriedUnsatsOperation::keyname();
}

PruneBuriedUnsatsOperation::PruneBuriedUnsatsOperation() :
	core::pack::task::operation::TaskOperation(),
	allow_even_trades_( false ),
	atomic_depth_probe_radius_( 2.3f ),
	atomic_depth_resolution_( 0.5f ),
	atomic_depth_cutoff_( 4.5f ),
	minimum_hbond_energy_( -0.2f )
{
	{
		scorefxn_sc_ = core::scoring::ScoreFunctionOP ( new core::scoring::ScoreFunction() );
		core::scoring::methods::EnergyMethodOptions opts = scorefxn_sc_->energy_method_options();
		core::scoring::hbonds::HBondOptions hbopts = opts.hbond_options();
		hbopts.use_hb_env_dep(false);
		opts.hbond_options( hbopts );
		scorefxn_sc_->set_energy_method_options( opts );
		scorefxn_sc_->set_weight( core::scoring::hbond_bb_sc, 1.0 );
		scorefxn_sc_->set_weight( core::scoring::hbond_sc, 1.0 );
	}
	{
		scorefxn_bb_ = core::scoring::ScoreFunctionOP ( new core::scoring::ScoreFunction() );
		core::scoring::methods::EnergyMethodOptions opts = scorefxn_bb_->energy_method_options();
		core::scoring::hbonds::HBondOptions hbopts = opts.hbond_options();
		hbopts.use_hb_env_dep(false);
		opts.hbond_options( hbopts );
		scorefxn_bb_->set_energy_method_options( opts );
		scorefxn_bb_->set_weight( core::scoring::hbond_sr_bb, 1.0 );
		scorefxn_bb_->set_weight( core::scoring::hbond_lr_bb, 1.0 );
	}
}

/// @brief default constructor
PruneBuriedUnsatsOperation::PruneBuriedUnsatsOperation( PruneBuriedUnsatsOperation const & other ) :
	TaskOperation( other ),
	allow_even_trades_( other.allow_even_trades_ ),
	atomic_depth_probe_radius_( other.atomic_depth_probe_radius_ ),
	atomic_depth_resolution_( other.atomic_depth_resolution_ ),
	atomic_depth_cutoff_( other.atomic_depth_cutoff_ ),
	minimum_hbond_energy_( other.minimum_hbond_energy_),
	scorefxn_sc_( other.scorefxn_sc_->clone() ),
	scorefxn_bb_( other.scorefxn_bb_->clone() )
{}

/// @brief destructor
PruneBuriedUnsatsOperation::~PruneBuriedUnsatsOperation()= default;

PruneBuriedUnsatsOperation &
PruneBuriedUnsatsOperation::operator=( PruneBuriedUnsatsOperation const & other ) {
	TaskOperation::operator=(other);
	allow_even_trades_ = other.allow_even_trades_;
	atomic_depth_probe_radius_ = other.atomic_depth_probe_radius_;
	atomic_depth_resolution_ = other.atomic_depth_resolution_;
	atomic_depth_cutoff_ = other.atomic_depth_cutoff_;
	minimum_hbond_energy_ = other.minimum_hbond_energy_;
	scorefxn_sc_ = other.scorefxn_sc_->clone();
	scorefxn_bb_ = other.scorefxn_bb_->clone();
	return *this;
}

/// @brief clone
core::pack::task::operation::TaskOperationOP
PruneBuriedUnsatsOperation::clone() const {
	return core::pack::task::operation::TaskOperationOP( new PruneBuriedUnsatsOperation( *this ) );
}

/// @brief
void
PruneBuriedUnsatsOperation::apply( Pose const & /*pose*/, core::pack::task::PackerTask & task ) const
{
	PruneBuriedUnsats_RotamerSetsOperationOP rso( new PruneBuriedUnsats_RotamerSetsOperation(
		allow_even_trades_, atomic_depth_probe_radius_, atomic_depth_resolution_, atomic_depth_cutoff_,
		minimum_hbond_energy_, scorefxn_sc_, scorefxn_bb_ ) );
	task.append_rotamersets_operation( rso );
}

void
PruneBuriedUnsatsOperation::parse_tag( TagCOP tag , DataMap & )
{
	allow_even_trades( tag->getOption< bool >( "allow_even_trades", false ) );
	atomic_depth_probe_radius( tag->getOption< core::Real >( "atomic_depth_probe_radius", 2.3f ) );
	atomic_depth_resolution( tag->getOption< core::Real >( "atomic_depth_resolution", 0.5f ) );
	atomic_depth_cutoff( tag->getOption< core::Real >( "atomic_depth_cutoff", 4.5f ) );
	minimum_hbond_energy( tag->getOption< core::Real >( "minimum_hbond_energy", -0.2f ) );
}

void PruneBuriedUnsatsOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::attribute_w_default(
		"allow_even_trades", xsct_rosetta_bool,
		"Allow residues that satisfy an unsat and create a new unsatisfiable one.",
		"false"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"atomic_depth_probe_radius", xsct_real,
		"Probe radius for atomic depth calculation to determine burial.",
		"2.3"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"atomic_depth_resolution", xsct_real,
		"Voxel resolution with which to calculate atomic depth.",
		"0.5"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"atomic_depth_cutoff", xsct_real,
		"Atomic depth at which atoms are considered buried.",
		"4.5"  )
		+ XMLSchemaAttribute::attribute_w_default(
		"minimum_hbond_energy", xsct_real,
		"Minimum energy (out of the typical rosetta -2.0) for a hbond to be considered to satisfy a polar.",
		"-0.2"  );

	task_op_schema_w_attributes(
		xsd, keyname(), attributes,
		"Eliminate rotamers that contain a sidechain buried unsat and make no other sidechain h-bonds.");
}


} // TaskOperations
} // protocols

