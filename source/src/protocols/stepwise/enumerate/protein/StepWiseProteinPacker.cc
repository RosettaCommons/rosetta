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
#include <protocols/stepwise/enumerate/protein/StepWiseProteinPacker.hh>
#include <protocols/stepwise/enumerate/protein/StepWiseProteinUtil.hh>
#include <protocols/stepwise/enumerate/protein/sample_generators/StepWisePoseSampleGenerator.hh>
#include <protocols/stepwise/enumerate/protein/PoseFilter.hh>

//////////////////////////////////
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.hh>
#include <core/types.hh>
#include <core/graph/Graph.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pose/Pose.hh>

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


#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/exit.hh>

#include <string>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto using namespaces
using namespace ObjexxFCL; // AUTO USING NS


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

static basic::Tracer TR( "protocols.stepwise.protein.StepWiseProteinPacker" ) ;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseProteinPacker::StepWiseProteinPacker(
													 utility::vector1< Size > const & moving_residues,
													 protocols::stepwise::enumerate::protein::sample_generators::StepWisePoseSampleGeneratorOP sample_generator ):
		moving_residues_( moving_residues ),
		scorefxn_( core::scoring::getScoreFunction() ),
		green_packer_( new protocols::simple_moves::GreenPacker ),
		use_green_packer_( false ),
		use_packer_instead_of_rotamer_trials_( false ),
		pack_at_neighbors_only_( true ),
		rescore_only_( false ),
		silent_file_( "" ),
		sfd_( new core::io::silent::SilentFileData),
		sample_generator_( sample_generator )
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


  //////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinPacker::apply( core::pose::Pose & pose )
	{

		//Size which_res( 1 );
		//		Size count( 1 );

		clock_t const time_start( clock() );

		if ( pose.is_fullatom() ){
			if ( use_green_packer_ ) {
				initialize_green_packer( pose.total_residue() );
			} else {
				initialize_for_regular_packer( pose );
			}
		}

		sample_residues( pose );

		std::cout << "Total time in StepWiseProteinPacker: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

	}


	////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::sample_residues( core::pose::Pose & pose )
	{

		 using namespace core::chemical;
		 using namespace core::scoring;
		 using namespace core::pose;
		 using namespace core::kinematics;
		 using namespace protocols::stepwise;

		 sample_generator_->reset();

		 Size k( 0 );

		 while( sample_generator_->has_another_sample() ){

			 sample_generator_->get_next_sample( pose );
			 k++;

			 if ( pose_filter_ && !pose_filter_->passes_filter( pose ) ) continue;

			 std::string const tag = "S_"+ lead_zero_string_of( k-1, 5 );

			 if ( pose.is_fullatom() || ( k % 100 == 0 ) ) print_tag( tag, k );

			 if ( pose.is_fullatom() && !rescore_only_ ){
				 if ( use_green_packer_ ) {
					 green_packer_->apply( pose );
				 } else {
					 apply_regular_packer( pose );
				 }
			 }

			 (*scorefxn_)( pose );

			 output_silent_struct( pose, get_native_pose(), silent_file_, tag, sfd_,
														 calc_rms_res_ );

		 }

		 //Nothing found? At least produce one pose...
		 if ( sfd_->size() == 0 ) {
			 TR << "Warning -- nothing passed filter -- outputting a placeholder pose." << std::endl;;
			 std::string const tag = "S_"+ lead_zero_string_of( k-1, 5 );
			 output_silent_struct( pose, get_native_pose(), silent_file_, tag, sfd_,
														 calc_rms_res_ );
		 }


	}

	////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::print_tag( std::string const & tag, Size const k ) {
		std::cout << " Decoy " << tag << " : ";
		std::string const pack_or_trials = use_packer_instead_of_rotamer_trials_ ? "PACK" : "ROT_TRIALS" ;
		std::cout << pack_or_trials;
		std::cout << "   Number " << k;
		if ( sample_generator_->size() > 0 ) std::cout  << " out of " << sample_generator_->size();
		std::cout << std::endl;
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
	StepWiseProteinPacker::initialize_for_regular_packer( core::pose::Pose & pose ){
		using namespace core::id;

		pack_task_ = pack::task::TaskFactory::create_packer_task( pose );
		pack_task_->restrict_to_repacking();
		for (Size i = 1; i <= pose.total_residue(); i++) {

			pack_task_->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
			pack_task_->nonconst_residue_task(i).or_include_current( true );

			if ( pose.residue(i).is_protein() ) {
				pack_task_->nonconst_residue_task(i).or_ex1( true );
				pack_task_->nonconst_residue_task(i).or_ex2( true );
				//			pack_task_->nonconst_residue_task(i).or_ex3( true );
				//			pack_task_->nonconst_residue_task(i).or_ex4( true );
			} else if ( pose.residue(i).is_RNA() ) {
			// Following could be useful...
				pack_task_->nonconst_residue_task(i).or_ex4( true ); //extra rotamers?? Parin S. Jan 28, 2010
				// new -- in case RNA backbone only is input?
				//pack_task_->nonconst_residue_task(i).sample_rna_chi( true );
			}

			if ( pose.residue(i).has_variant_type( "VIRTUAL_RESIDUE" ) ) pack_task_->nonconst_residue_task(i).prevent_repacking();
		}

		// save this pose for later.
		pose_init_ = new Pose;
		*pose_init_ = pose;

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
			for ( Size n = 1; n <= pose.residue_type( i ).nchi(); n++ ) {
				pose.set_chi( n, i, src_pose.chi( n, i ) );
			}
		}
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::apply_regular_packer( core::pose::Pose & pose ){

		if ( pack_at_neighbors_only_ ) {

			// need to reinstate side-chain angles at all positions to "pre-packed values"
			reinstate_side_chain_angles( pose, *pose_init_ );

			// figure out neighbors
			utility::vector1< bool > residues_allowed_to_be_packed( pose.total_residue(), false );
			figure_out_neighbors( pose, residues_allowed_to_be_packed );

			//set up new task
			pack_task_->restrict_to_repacking();

			//			std::cout << "PACKING ==> ";
			for (Size i = 1; i <= pose.total_residue(); i++) {
				if (  residues_allowed_to_be_packed[ i ] )  {
					pack_task_->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
					pack_task_->nonconst_residue_task(i).or_include_current( true );
					if ( pose.residue(i).is_protein() ) {
						pack_task_->nonconst_residue_task(i).or_ex1( true );
						pack_task_->nonconst_residue_task(i).or_ex2( true );
					} else if ( pose.residue(i).is_RNA() ){
						pack_task_->nonconst_residue_task(i).or_ex4( true ); //extra rotamers?? Parin S. Jan 28, 2010
						// new -- in case RNA backbone only is input?
						//pack_task_->nonconst_residue_task(i).sample_rna_chi( true );
					}
				} else {
					pack_task_->nonconst_residue_task(i).prevent_repacking();
				}
			}

			// NEW -- more rotamers at moving residues...
			//			for ( Size n = 1; n <= moving_residues_.size(); n++ ) {
			//				Size const i = moving_residues_[ n ];
			//std::cout << "MORE ROTAMERS TO " << i << std::endl;
			//				pack_task_->nonconst_residue_task(i).or_ex1_sample_level( core::pack::task::EX_SIX_QUARTER_STEP_STDDEVS );
			//				pack_task_->nonconst_residue_task(i).or_ex2_sample_level( core::pack::task::EX_SIX_QUARTER_STEP_STDDEVS );
			//			}

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
  StepWiseProteinPacker::set_silent_file( std::string const & silent_file ){
    silent_file_ = silent_file;
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

  //////////////////////////////////////////////////////////////////////////
	core::io::silent::SilentFileDataOP &
	StepWiseProteinPacker::silent_file_data(){
		return sfd_;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res ){
		calc_rms_res_ = calc_rms_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPacker::set_pose_filter( protocols::stepwise::enumerate::protein::PoseFilterOP pose_filter ){
		pose_filter_ = pose_filter;
	}

} //protein
} //enumerate
} //stepwise
} //protocols
