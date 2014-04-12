// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinPacker
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/sampling/protein/StepWiseProteinPacker.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinUtil.hh>
#include <protocols/stepwise/StepWiseUtil.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <basic/Tracer.hh>

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>

#include <utility/exit.hh>

#include <string>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

using namespace core;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of proteins (and probably other
// biopolymers soon). Take a starting pose and a list of residues to sample,
//  and comprehensively sample all backbone torsion angles by recursion.
//
// Probably should be folded into a sample-and-screen framework!
// See StepWiseRNA_Modeler for example.
//   -- Rhiju, 2014.
//
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.stepwise.protein.StepWiseProteinPacker" ) ;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseProteinPacker::StepWiseProteinPacker( utility::vector1< Size > const & moving_residues ):
		moving_residues_( moving_residues ),
		green_packer_( new protocols::simple_moves::GreenPacker ),
		use_green_packer_( false ),
		use_packer_instead_of_rotamer_trials_( false ),
		pack_at_neighbors_only_( true ),
		allow_virtual_side_chains_( false ),
		rescore_only_( false )
  {
  }

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseProteinPacker::~StepWiseProteinPacker()
  {}

	/////////////////////
	std::string
	StepWiseProteinPacker::get_name() const {
		return "StepWiseProteinPacker";
	}

	////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::initialize( pose::Pose & pose ) {
		if ( pose.is_fullatom() ){
			if ( use_green_packer_ ) {
				initialize_green_packer( pose.total_residue() );
			} else {
				initialize_for_regular_packer( pose );
			}
		}
	}

	////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::apply( pose::Pose & pose ) {
		if ( pose.is_fullatom() && !rescore_only_ ){
			if ( use_green_packer_ ) {
				green_packer_->apply( pose );
			} else {
				apply_regular_packer( pose, pack_at_neighbors_only_ );
			}
		}
		(*scorefxn_)( pose );
	}

	////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::initialize_green_packer( Size const & nres )
	{
		using namespace protocols::moves;
		using namespace core::pack;
		using namespace core::pack::task;
		using namespace core::pack::task::operation;

		protocols::simple_moves::UserDefinedGroupDiscriminatorOP user_defined_group_discriminator( new protocols::simple_moves::UserDefinedGroupDiscriminator);
		utility::vector1< Size > group_ids;

		Size current_group = 0;
		Size spectator_group = 0;
		for (Size i = 1; i <= nres; i++ ) {
			bool found_it( false );
			for (Size k = 1; k <= moving_residues_.size(); k++ ) {
				if ( i == moving_residues_[k] ) {
					found_it = true;
					break;
				}
			}
			if (found_it ) {
				current_group = 0;
				TR << "GREENPACKER SAMPLER " << i << std::endl;
			} else {
				if ( current_group == 0 ) spectator_group++;
				current_group = spectator_group;
				TR << "GREENPACKER SPECTATOR   " << i <<  " --> group " << spectator_group << std::endl;
			}
			group_ids.push_back( current_group );
		}

		user_defined_group_discriminator->set_group_ids( group_ids );
		green_packer_->set_scorefunction( *scorefxn_ );
		green_packer_->set_group_discriminator( user_defined_group_discriminator );

		TaskFactoryOP initial_task_factory( new TaskFactory );
		initial_task_factory->push_back( new InitializeFromCommandline );
		initial_task_factory->push_back( new RestrictToRepacking );
		green_packer_->set_reference_round_task_factory( initial_task_factory );

		TaskFactoryOP general_task_factory( new TaskFactory );
		general_task_factory->push_back( new InitializeFromCommandline );
		general_task_factory->push_back( new RestrictToRepacking );
		green_packer_->set_task_factory( general_task_factory );

		//green_packer_->reset();
	}


	////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::initialize_for_regular_packer( core::pose::Pose const & pose ){
		using namespace core::id;

		// When setting up rotamers, currently RotamerSet_ assumes that concrete residues do *not*
		// have virtual side chains.
		// This hack is kind of silly, and should be fixed by TaskFactory handling virtual side-chains.
		Pose pose_without_virtual_side_chains = pose;
		for ( Size i = 1; i <= pose.total_residue(); i++ ) remove_variant_type_from_pose_residue( pose_without_virtual_side_chains, "VIRTUAL_SIDE_CHAIN", i );

		pack_task_ = pack::task::TaskFactory::create_packer_task( pose_without_virtual_side_chains );
		pack_task_->restrict_to_repacking();
		for (Size i = 1; i <= pose.total_residue(); i++) {

			pack_task_->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
			pack_task_->nonconst_residue_task(i).or_include_current( true );

			if ( pose.residue(i).is_protein() ) {
				pack_task_->nonconst_residue_task(i).or_ex1( true );
				pack_task_->nonconst_residue_task(i).or_ex2( true );
				//			pack_task_->nonconst_residue_task(i).or_ex3( true );
				//			pack_task_->nonconst_residue_task(i).or_ex4( true );
				pack_task_->nonconst_residue_task(i).or_include_virtual_side_chain( allow_virtual_side_chains_ );
			} else if ( pose.residue(i).is_RNA() ) {
				pack_task_->nonconst_residue_task(i).or_ex4( true );
			}

			if ( pose.residue(i).has_variant_type( "VIRTUAL_RESIDUE" ) ) pack_task_->nonconst_residue_task(i).prevent_repacking();
		}

		// save this pose for later.
		pose_init_ = pose.clone();

	}


  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::figure_out_neighbors( core::pose::Pose & pose,
																											 utility::vector1< bool > & residues_allowed_to_be_packed ){

		using namespace core::scoring;

		(*scorefxn_)( pose );
		EnergyGraph const & energy_graph( pose.energies().energy_graph() );

		for ( Size n = 1; n <= moving_residues_.size(); n++ ) {

			Size const i = moving_residues_[ n ];
			residues_allowed_to_be_packed[ i ] = true;

			for( graph::Graph::EdgeListConstIter
						 iter = energy_graph.get_node( i )->const_edge_list_begin();
					 iter != energy_graph.get_node( i )->const_edge_list_end();
					 ++iter ){
				Size j( (*iter)->get_other_ind( i ) );
				residues_allowed_to_be_packed[ j ] = true;
			}
		}

	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::reinstate_side_chain_angles( pose::Pose & pose, pose::Pose const & src_pose ){

		for ( Size i = 1; i <= pose.total_residue(); i++ ) {
			make_variants_match( pose, src_pose, i, "VIRTUAL_SIDE_CHAIN" );
			for ( Size n = 1; n <= pose.residue_type( i ).nchi(); n++ ) {
				pose.set_chi( n, i, src_pose.chi( n, i ) );
			}
		}
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::apply_regular_packer( core::pose::Pose & pose, bool const pack_at_neighbors_only ){

		if ( pack_at_neighbors_only ) {

			// need to reinstate side-chain angles at all positions to "pre-packed values"
			reinstate_side_chain_angles( pose, *pose_init_ );

			// figure out neighbors
			utility::vector1< bool > residues_allowed_to_be_packed( pose.total_residue(), false );
			figure_out_neighbors( pose, residues_allowed_to_be_packed );

			//set up new task
			pack_task_->restrict_to_repacking();

			for (Size i = 1; i <= pose.total_residue(); i++) {
				if (  residues_allowed_to_be_packed[ i ] )  {
					pack_task_->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
					pack_task_->nonconst_residue_task(i).or_include_current( true );
					if ( pose.residue(i).is_protein() ) {
						pack_task_->nonconst_residue_task(i).or_ex1( true );
						pack_task_->nonconst_residue_task(i).or_ex2( true );
					} else if ( pose.residue(i).is_RNA() ){
						pack_task_->nonconst_residue_task(i).or_ex4( true );
					}
				} else {
					pack_task_->nonconst_residue_task(i).prevent_repacking();
				}
			}
		}


		// OK, pack!!
		if ( use_packer_instead_of_rotamer_trials_ ) {
			pack::pack_rotamers(  pose, *scorefxn_, pack_task_ );
		} else {
			pack::rotamer_trials( pose, *scorefxn_, pack_task_ );
		}

	}

  //////////////////////////////////////////////////////////////////////////
	// splits into far-apart partitions and packs. Packing can include
	// virtual side chains (if -allow_virtual_side_chains) is on.
	void
	StepWiseProteinPacker::do_prepack( pose::Pose & pose ){

		pose::Pose pose_to_split = pose;
		pose_to_split.remove_constraints(); // floating point errors if coordinate constraints are in there.
		split_pose( pose_to_split, moving_res_list_ );

		initialize_for_regular_packer( pose_to_split );
		apply_regular_packer( pose_to_split, false /*pack_at_neighbors_only*/ );

		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			make_variants_match( pose, pose_to_split, n, "VIRTUAL_SIDE_CHAIN" );
			std::map< Size, Size > res_map;
			res_map[ n ] = n;
			copy_dofs_match_atom_names( pose, pose_to_split, res_map, false /*backbone_only*/, false /*ignore_virtual*/ );
		}

	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
		scorefxn_ = scorefxn;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::set_use_green_packer( bool const & setting ){
		use_green_packer_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::set_use_packer_instead_of_rotamer_trials( bool const & setting ){
		use_packer_instead_of_rotamer_trials_ = setting;
	}


} //protein
} //sampling
} //stepwise
} //protocols
