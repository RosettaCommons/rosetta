// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinPoseMinimizer
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/enumerate/protein/StepWiseProteinPoseMinimizer.hh>
#include <protocols/stepwise/enumerate/protein/StepWiseProteinUtil.hh>
#include <protocols/stepwise/StepWiseUtil.hh>

//////////////////////////////////
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/types.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>
#include <ObjexxFCL/format.hh>


#include <string>

//Auto Headers
#include <utility/vector1.hh>
using namespace core;
using core::Real;
using ObjexxFCL::format::F;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseProteinPoseMinimizer::StepWiseProteinPoseMinimizer( core::io::silent::SilentFileDataOP const sfd, utility::vector1< Size > const & moving_residues ):
		Mover(),
    moving_residues_( moving_residues )
  {
		input_silent_file_data_ = sfd;
		initialize_parameters();
	}

  //////////////////////////////////////////////////////////////////////////
  StepWiseProteinPoseMinimizer::StepWiseProteinPoseMinimizer( PoseList & pose_list, utility::vector1< Size > const & moving_residues ):
		Mover(),
    moving_residues_( moving_residues )
  {
		initialize_protein_input_silent_file_data_from_pose_list( pose_list );
		initialize_parameters();
  }

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseProteinPoseMinimizer::~StepWiseProteinPoseMinimizer()
  {}

	/////////////////////
	std::string
	StepWiseProteinPoseMinimizer::get_name() const {
		return "StepWiseProteinPoseMinimizer";
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseMinimizer::initialize_parameters(){
		Mover::type( "StepWiseProteinPoseMinimizer" );
		move_takeoff_torsions_ = true;
		rescore_only_ = false;
		move_jumps_between_chains_ = false;
		silent_file_ = "";
    fa_scorefxn_ = core::scoring::getScoreFunction();
		min_type_ = "dfpmin_armijo_nonmonotone"; // used to be dfpmin
		cartesian_ = true;
		min_tolerance_ = 0.000025 ; // used to be 0.00000025
	}


  //////////////////////////////////////////////////////////////////////////
	void StepWiseProteinPoseMinimizer::initialize_protein_input_silent_file_data_from_pose_list( PoseList & /*pose_list*/ ){
    using namespace core::pose;
    using namespace core::io::silent;

		input_silent_file_data_->clear();
    for ( PoseList::iterator iter = pose_list_.begin(); iter != pose_list_.end(); iter++ ) {
			PoseOP & pose_op( iter->second );
			BinaryProteinSilentStruct s( *pose_op, iter->first /*tag*/ );
			input_silent_file_data_->add_structure( s );
		}
	}

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinPoseMinimizer::apply( core::pose::Pose & pose )
  {
    using namespace core::optimization;
    using namespace core::scoring;
    using namespace core::scoring::constraints;
    using namespace core::pose;
    using namespace core::io::silent;

		clock_t const time_start( clock() );

		ConstraintSetOP cst_set = pose.constraint_set()->clone();

		utility::vector1< std::pair<core::Size,core::Size> > disulfides;
		core::conformation::disulfide_bonds(pose.conformation(), disulfides);

		//if ( cartesian_ ){ // this is messy -- move to its own function if it works.
		//			for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		//				if ( pose.residue_type(i).has_variant_type( "CUTPOINT_UPPER" ) ) remove_variant_type_from_pose_residue( pose, "CUTPOINT_UPPER", i );
		//				if ( pose.residue_type(i).has_variant_type( "CUTPOINT_LOWER" ) ) remove_variant_type_from_pose_residue( pose, "CUTPOINT_LOWER", i );
		//			}
		//		}

		CartesianMinimizer cart_minimizer;
		AtomTreeMinimizer minimizer;
    bool const use_nblist( true );
    MinimizerOptions options( min_type_, min_tolerance_, use_nblist, false, false );
    options.nblist_auto_update( true );

    kinematics::MoveMap mm_start, mm;
		std::cout << "MOVE TAKEOFF TORSIONS: " << move_takeoff_torsions_ << std::endl;
		protocols::stepwise::figure_out_moving_residues( mm_start, pose, fixed_res_, move_takeoff_torsions_, move_jumps_between_chains_ );
		mm = mm_start;

		//		using namespace core::id;
		//		for ( Size i = 1; i <= pose.total_residue(); i++ ){
		//			std::cout << " MM: " << i << ' ' << mm.get( TorsionID( i, BB, 1 ) ) << ' ' << mm.get( TorsionID( i, BB, 2 ) ) << ' ' << mm.get( TorsionID( i, BB, 3 ) ) << std::endl;
		//		}

    Size count( 1 );
		Real const original_coordinate_cst_weight = fa_scorefxn_->get_weight( coordinate_constraint );

		sfd_ = new SilentFileData;
		core::import_pose::pose_stream::SilentFilePoseInputStream input;
		input.set_silent_file_data( input_silent_file_data_ );

		while ( input.has_another_pose() ) {

      std::cout << "Minimizing decoy " << count++ << " out of " << input_silent_file_data_->size() << std::endl;

			input.fill_pose( pose );

			// Following are necessary because poses from clustering went thorugh silent struct and lost their constraints & disulfide information.
			pose.constraint_set( cst_set );
			pose.conformation().fix_disulfides( disulfides );

			Real const score_original = (*fa_scorefxn_)( pose );

			// The movemap has all dofs for "non-fixed residues" free to move.
			// We can also let sidechains minimize in fixed-residues -- for
			// speed only look at neighbors of moving residues.
			mm = mm_start;
			let_neighboring_chis_minimize( mm, pose );

			if ( !rescore_only_ ){

				// One minimize with loose coordinate tethers to make sure the pose doesn't blow up.
				core::scoring::constraints::add_coordinate_constraints( pose );
				if ( fa_scorefxn_->has_zero_weight( coordinate_constraint) ) fa_scorefxn_->set_weight( coordinate_constraint, 1.0 );
				minimizer.run( pose, mm, *fa_scorefxn_, options );

				// Now a regular minimize.
				pose.constraint_set( cst_set ); // return original constraints (no added coordinate constraints)
				fa_scorefxn_->set_weight( coordinate_constraint, original_coordinate_cst_weight );

				// for poses with chainbreaks, do an initial minimization with a weak linear_chainbreak term. (anneal it in.)
				if ( pose_has_chainbreak( pose ) ){

					Real const linear_chainbreak_weight_original = fa_scorefxn_->get_weight( linear_chainbreak );
					if ( linear_chainbreak_weight_original < 20.0 ) std::cout << "WARNING!! Your linear_chainbreak weight is " << F(8,3,linear_chainbreak_weight_original ) << ", which is less than recommended (20.0) " << std::endl;

					fa_scorefxn_->set_weight( linear_chainbreak, linear_chainbreak_weight_original * 0.25 );
					if ( !cartesian_) {
						minimizer.run( pose, mm, *fa_scorefxn_, options );
					}
					fa_scorefxn_->set_weight( linear_chainbreak, linear_chainbreak_weight_original );
				}


				if (cartesian_) {
					cart_minimizer.run( pose, mm, *fa_scorefxn_, options );
				} else {
					minimizer.run( pose, mm, *fa_scorefxn_, options );
				}

			}

			setPoseExtraScores( pose, "score_orig", score_original );
      std::string const & tag( tag_from_pose( pose ) );
      protocols::stepwise::enumerate::protein::output_silent_struct( pose, get_native_pose(), silent_file_, tag, sfd_, calc_rms_res_ );

			std::cout << "Score minimized from " <<F(8,3, score_original) << " to " << F(8,3,(*fa_scorefxn_)( pose )) << std::endl;

			// Running into file locking issues
			//			utility::sys_sleep( 0.5 );
			//exit( 0 );

      // Might was well replace pose in original list.
      //*pose_op = pose;

    }

		std::cout << "Total time in StepWiseProteinPoseMinimizer: " <<
			static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;

  }

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseMinimizer::let_neighboring_chis_minimize(
																															core::kinematics::MoveMap & mm,
																															core::pose::Pose & pose ){

		using namespace core::scoring;

		(*fa_scorefxn_)( pose );
		EnergyGraph const & energy_graph( pose.energies().energy_graph() );

		for ( Size n = 1; n <= moving_residues_.size(); n++ ) {

			Size const i = moving_residues_[ n ];

			for( graph::Graph::EdgeListConstIter
						 iter = energy_graph.get_node( i )->const_edge_list_begin();
					 iter != energy_graph.get_node( i )->const_edge_list_end();
					 ++iter ){

				Size j( (*iter)->get_other_ind( i ) );
				if ( pose.residue(j).has_variant_type( "VIRTUAL_RESIDUE" ) ) continue;

				if ( pose.residue(j).is_protein() ){
					mm.set_chi( j, true );
				} else if ( pose.residue(j).is_RNA() ){
					mm.set( id::TorsionID( j, id::CHI, 4), true ); // 2'-OH.
				}

			}
		}

	}

  //////////////////////////////////////////////////////////////////////////
	bool
	StepWiseProteinPoseMinimizer::pose_has_chainbreak( pose::Pose const & pose ){
		// this is pretty conservative -- actually  there might be
		// cases where the pose has a chainbreak but the minimized dofs would
		// not affect the relative positions of the chainbreak residues.
		for ( Size i = 1; i <= pose.total_residue(); i++ ){
			if ( pose.residue_type(i).has_variant_type( "CUTPOINT_UPPER" ) ) return true;
			if ( pose.residue_type(i).has_variant_type( "CUTPOINT_LOWER" ) ) return true;
		}
		return false;
	}

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinPoseMinimizer::set_silent_file( std::string const & silent_file ){
    silent_file_ = silent_file;
  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinPoseMinimizer::set_min_type( std::string const & min_type ){
    min_type_ = min_type;
  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinPoseMinimizer::set_min_tolerance( Real const & min_tolerance ){
    min_tolerance_ = min_tolerance;
  }

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseMinimizer::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
		fa_scorefxn_ = scorefxn;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseMinimizer::set_fixed_res( utility::vector1< core::Size > const & fixed_res ){
		fixed_res_ = fixed_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseMinimizer::set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res ){
		calc_rms_res_ = calc_rms_res;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinPoseMinimizer::set_cartesian( bool const setting ){
		cartesian_ = setting;
		if (cartesian_) min_type_ = "lbfgs_armijo_nonmonotone";
	}

  //////////////////////////////////////////////////////////////////////////
	core::io::silent::SilentFileDataOP &
	StepWiseProteinPoseMinimizer::silent_file_data(){
		return sfd_;
	}

	//	void
	//	StepWiseProteinPoseMinimizer::set_constraint_set( core::scoring::constraints::ConstraintSetOP const & cst_set ){
	//		cst_set_ = cst_set;
	//	}


} //protein
} //enumerate
} //stepwise
} //protocols
