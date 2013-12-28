// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinResidueSampler
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/enumerate/protein/StepWiseProteinResidueSampler.hh>
#include <protocols/stepwise/enumerate/protein/StepWiseProteinUtil.hh>

//////////////////////////////////
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
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
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <basic/Tracer.hh>

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>


#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/exit.hh>

#include <string>

//Auto Headers
#include <core/id/TorsionID.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
//Auto using namespaces end


using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of proteins (and probably other
// biopolymers soon). Take a starting pose and a list of residues to sample,
//  and comprehensively sample all backbone torsion angles by recursion.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.stepwise.protein.StepWiseProteinResidueSampler" ) ;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseProteinResidueSampler::StepWiseProteinResidueSampler(
																								 utility::vector1< Size > const & moving_residues,
																								 utility::vector1< core::id::TorsionID > const & which_torsions,
																								 utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_set_lists ):
		moving_residues_( moving_residues ),
		which_torsions_( which_torsions ),
		main_chain_torsion_set_lists_( main_chain_torsion_set_lists ),
		scorefxn_( core::scoring::getScoreFunction() ),
		green_packer_( new protocols::simple_moves::GreenPacker ),
		use_green_packer_( false ),
		use_packer_instead_of_rotamer_trials_( false ),
		pack_at_neighbors_only_( true ),
		do_prepack_( false ),
		silent_file_( "" ),
		sfd_( new core::io::silent::SilentFileData),
		which_jump_( 0 )
  {
  }

  StepWiseProteinResidueSampler::StepWiseProteinResidueSampler(
																								 utility::vector1< Size > const & moving_residues,
																								 utility::vector1< core::id::TorsionID > const & which_torsions,
																								 utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_set_lists,
																								 Size const which_jump,
																								 utility::vector1< core::kinematics::Jump > const & jumps ):
		moving_residues_( moving_residues ),
		which_torsions_( which_torsions ),
		main_chain_torsion_set_lists_( main_chain_torsion_set_lists ),
		scorefxn_( core::scoring::getScoreFunction() ),
		green_packer_( new protocols::simple_moves::GreenPacker ),
		use_green_packer_( false ),
		use_packer_instead_of_rotamer_trials_( false ),
		pack_at_neighbors_only_( true ),
		do_prepack_( false ),
		silent_file_( "" ),
		sfd_( new core::io::silent::SilentFileData),
		which_jump_( which_jump ),
		jumps_( jumps )
  {
		if ( which_torsions.size() == 0 ){
			for ( Size n = 1; n <= jumps_.size(); n++ ) {
				utility::vector1< core::Real > dummy_main_chain_torsion_set;
				main_chain_torsion_set_lists_.push_back( dummy_main_chain_torsion_set );
			}
		}
  }

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseProteinResidueSampler::~StepWiseProteinResidueSampler()
  {}
/////////////////////
std::string
StepWiseProteinResidueSampler::get_name() const {
return "StepWiseProteinResidueSampler";
}


  //////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinResidueSampler::apply( core::pose::Pose & pose )
	{

		//Size which_res( 1 );
		//		Size count( 1 );

		clock_t const time_start( clock() );

		if ( use_green_packer_ ) {
			initialize_green_packer( pose.total_residue() );
		} else {
			initialize_for_regular_packer( pose );
		}

		sample_residues( pose );

		std::cout << "Total time in StepWiseProteinResidueSampler: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

	}


	////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::sample_residues( core::pose::Pose & pose )
	{

		 using namespace core::chemical;
		 using namespace core::scoring;
		 using namespace core::pose;
		 using namespace core::kinematics;

		 for ( Size k = 1; k <= main_chain_torsion_set_lists_.size(); k++ ) {

			 utility::vector1< Real > const & main_chain_torsion_set_list( main_chain_torsion_set_lists_[ k ] );

			 assert( main_chain_torsion_set_list.size() == which_torsions_.size() );

			 for ( Size i = 1; i <= which_torsions_.size(); i++ ) {
				 pose.set_torsion( which_torsions_[ i ], main_chain_torsion_set_list[ i ] );
			 }

			 if ( which_jump_ > 0 )	pose.set_jump( which_jump_, jumps_[ k ] );

			 std::string const tag = "S_"+ lead_zero_string_of( k-1, 5 );
			 std::cout << " Decoy " << tag << " : " << k << " out of " << main_chain_torsion_set_lists_.size() << std::endl;

			 if ( use_green_packer_ ) {
				 green_packer_->apply( pose );
			 } else {
				 apply_regular_packer( pose );
			 }

			 (*scorefxn_)( pose );

			 output_silent_struct( pose, get_native_pose(), silent_file_, tag, sfd_,
														 calc_rms_res_ );

		 }

		 if (sfd_->size() == 0 ){ //perhaps just did a prepack, no actual sampling. OK, save it.
			 assert( do_prepack_ );
			 output_silent_struct( pose, get_native_pose(), silent_file_, "S_PREPACK", sfd_,
														 calc_rms_res_ );
		 }

	}


	////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::initialize_green_packer( Size const & nres )
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
	StepWiseProteinResidueSampler::initialize_for_regular_packer( core::pose::Pose & pose ){
		using namespace core::id;

		/////////////////////////////////////////////////
		// Need to do a full "prepack". Want all side-chains outside the loop
		/// in their favorite states.
		/////////////////////////////////////////////////

		/////////////////////////////////////////////////
		// First space out the moving residues (will return them later)
		/////////////////////////////////////////////////
		DOF_ID dof_id;
		utility::vector1< DOF_ID > dofs_to_reset;
		utility::vector1< Real > dof_values;

		if ( do_prepack_ ){
			for ( Size i = 1; i <= moving_residues_.size(); i++ ) {
				dof_id = DOF_ID( named_atom_id_to_atom_id( NamedAtomID( " N  ", moving_residues_[ i ]), pose ), D );
				dofs_to_reset.push_back( dof_id );
				dof_values.push_back( pose.dof( dof_id ) );
				pose.set_dof( dof_id, 50.0 ); //dummy, extended value.

				dof_id = DOF_ID( named_atom_id_to_atom_id( NamedAtomID( " C  ", moving_residues_[ i ]), pose ), D );
				dofs_to_reset.push_back( dof_id );
				dof_values.push_back( pose.dof( dof_id ) );
				pose.set_dof( dof_id, 50.0 ); //dummy, extended value.
			}
			//pose.dump_pdb( "extended_before_pack.pdb" );
		}

		/////////////////////////////////////////////////
		// Then do the pack
		/////////////////////////////////////////////////
		pack_task_ = pack::task::TaskFactory::create_packer_task( pose );
		pack_task_->restrict_to_repacking();
		for (Size i = 1; i <= pose.total_residue(); i++) {
			if ( !pose.residue(i).is_protein() ) continue;
			pack_task_->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
			pack_task_->nonconst_residue_task(i).or_ex1( true );
			pack_task_->nonconst_residue_task(i).or_ex2( true );
			//			pack_task_->nonconst_residue_task(i).or_ex3( true );
			//			pack_task_->nonconst_residue_task(i).or_ex4( true );
			pack_task_->nonconst_residue_task(i).or_include_current( true );
		}

		if ( do_prepack_ ){
			Size const ntrials_prepack( 1 ); //is this necessary?
			// each repack doesn't take very long -- best to do multiple trials
			// and pick the lowest energy one.
			pack::pack_rotamers_loop( pose, *scorefxn_, pack_task_, ntrials_prepack );
			//pack::rotamer_trials( pose, *scorefxn_, pack_task_ );
		}

		/////////////////////////////////////////////////
		// Should we minimize the side chains too?
		/////////////////////////////////////////////////

		/////////////////////////////////////////////////
		// return spaced out moving residues.
		/////////////////////////////////////////////////
		if ( do_prepack_ ){
			for ( Size n = 1; n <= dofs_to_reset.size(); n++ ) {
				pose.set_dof( dofs_to_reset[ n ], dof_values[ n ] );
			}
		}

		// save this pose for later.
		pose_prepack_ = new Pose;
		*pose_prepack_ = pose;

	}


  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::figure_out_neighbors( core::pose::Pose & pose,
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
	reinstate_side_chain_angles( pose::Pose & pose, pose::Pose const & src_pose ){
		for ( Size i = 1; i <= pose.total_residue(); i++ ) {
			for ( Size n = 1; n <= pose.residue_type( i ).nchi(); n++ ) {
				pose.set_chi( n, i, src_pose.chi( n, i ) );
			}
		}
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::reinstate_side_chain_angles( pose::Pose & pose, pose::Pose const & src_pose ){
		for ( Size i = 1; i <= pose.total_residue(); i++ ) {
			for ( Size n = 1; n <= pose.residue_type( i ).nchi(); n++ ) {
				pose.set_chi( n, i, src_pose.chi( n, i ) );
			}
		}
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::apply_regular_packer( core::pose::Pose & pose ){

		if ( pack_at_neighbors_only_ ) {

			// need to reinstate side-chain angles at all positions to "pre-packed values"
			reinstate_side_chain_angles( pose, *pose_prepack_ );

			// figure out neighbors
			utility::vector1< bool > residues_allowed_to_be_packed( pose.total_residue(), false );
			figure_out_neighbors( pose, residues_allowed_to_be_packed );

			//set up new task
			pack_task_->restrict_to_repacking();

			//			std::cout << "PACKING ==> ";
			for (Size i = 1; i <= pose.total_residue(); i++) {
				if ( pose.residue(i).is_protein() && residues_allowed_to_be_packed[ i ] )  {
					//					std::cout << i << " " ;
					pack_task_->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
					pack_task_->nonconst_residue_task(i).or_ex1( true );
					pack_task_->nonconst_residue_task(i).or_ex2( true );
					pack_task_->nonconst_residue_task(i).or_include_current( true );
				} else {
					pack_task_->nonconst_residue_task(i).prevent_repacking();
				}
			}
			//			std::cout << std::endl;
		}


		// OK, pack!!
		if ( use_packer_instead_of_rotamer_trials_ ) {
			pack::pack_rotamers(  pose, *scorefxn_, pack_task_ );
		} else {
			pack::rotamer_trials( pose, *scorefxn_, pack_task_ );
		}

	}

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinResidueSampler::set_silent_file( std::string const & silent_file ){
    silent_file_ = silent_file;
  }

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
		scorefxn_ = scorefxn;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::set_use_green_packer( bool const & setting ){
		use_green_packer_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::set_use_packer_instead_of_rotamer_trials( bool const & setting ){
		use_packer_instead_of_rotamer_trials_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::set_do_prepack( bool const & setting ){
		do_prepack_ = setting;
	}

  //////////////////////////////////////////////////////////////////////////
	core::io::silent::SilentFileDataOP &
	StepWiseProteinResidueSampler::silent_file_data(){
		return sfd_;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinResidueSampler::set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res ){
		calc_rms_res_ = calc_rms_res;
	}


} //protein
} //enumerate
} //stepwise
} //protocols
