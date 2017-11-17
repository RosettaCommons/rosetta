// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/sevya/msd.cc
/// @brief  Experimental multi state design implementation. Written to encourage convergence in protein designs occurring
/// simultaneously, rather than enforce that they have identical sequences.
/// @author Alex Sevy


//core library
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <protocols/toolbox/task_operations/RestrictToInterfaceVectorOperation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/constraints/ResidueTypeLinkingConstraint.hh>
#include <core/chemical/ResidueType.hh>
#include <protocols/simple_moves/MinPackMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/NegativePackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <devel/init.hh>
#include <protocols/protein_interface_design/design_utils.hh>
// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <utility/io/izstream.hh>
#include <utility/json_spirit/json_spirit_reader.h>
#include <protocols/jobdist/standard_mains.hh>
#include <basic/Tracer.hh>
#include <sstream>
#include <utility/string_util.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>


namespace basic { namespace options { namespace OptionKeys { namespace msd
{
basic::options::StringVectorOptionKey positive_states( "msd:positive_states" );
basic::options::StringVectorOptionKey negative_states( "msd:negative_states" );
basic::options::StringOptionKey upper_chains( "msd:upper_chains" );
basic::options::StringOptionKey lower_chains( "msd:lower_chains" );
basic::options::IntegerOptionKey n_iterations( "msd:n_iterations" );
basic::options::BooleanOptionKey find_consensus( "msd:find_consensus" );
basic::options::BooleanOptionKey debug( "msd:debug" );
basic::options::BooleanOptionKey restrict_design_to_interface( "msd:restrict_design_to_interface" );
basic::options::BooleanOptionKey restrict_min_to_interface( "msd:restrict_min_to_interface" );
basic::options::IntegerOptionKey design_rounds( "msd:design_rounds" );
basic::options::StringOptionKey ramp( "msd:weight_ramp");
basic::options::RealOptionKey favor_native( "msd:favor_native" );
basic::options::BooleanOptionKey soft_design( "msd:soft_design" );
basic::options::BooleanOptionKey fixed_bb( "msd:fixed_bb" );
basic::options::FileVectorOptionKey multiple_resfiles(" msd:multiple_resfiles");
basic::options::FileVectorOptionKey multiple_constraints(" msd:multiple_constraints" );
}}}}


///////////////////////////////////////////////////////////////////////////////

static basic::Tracer TR("msd_pilot.main");

std::string
get_out_number(core::Size in) {
	std::stringstream num ("");
	if ( in < 10 ) {
		num << "000" << in;
	} else if ( in < 100 ) {
		num << "00" << in;
	} else if ( in < 1000 ) {
		num << "0" << in;
	} else {
		num << in;
	}
	return num.str();
}

void
add_fnr_constraints( core::pose::PoseOP pose,
	core::pose::PoseOP reference_pose,
	utility::vector1< core::Size > corr_pairs,
	core::Real penalty, bool positive ) {
	for ( core::Size k = 1; k <= corr_pairs.size(); ++k ) {
		pose->add_constraint( new core::scoring::constraints::ResidueTypeConstraint(
			*pose,
			corr_pairs[ k ],
			reference_pose->residue_type( corr_pairs[ k ] ).name3(),
			(positive ? penalty : -penalty ) ) );
	}
}

void
apply_linked_constraints(utility::vector1< std::pair< core::pose::PoseOP, bool > > poses,
	utility::vector1< utility::vector1< core::Size > > corr_pairs,
	core::Real positive_penalty, core::Real negative_penalty,
	core::Size current_pose )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	/* Iterate through poses in the outer loop */
	for ( core::Size i = 1; i <= poses.size(); ++i ) {
		/* Iterate through residues to be linked */
		for ( core::Size k = 1; k <= corr_pairs[ current_pose ].size(); ++k ) {
			if ( i == current_pose ) continue;

			poses[ current_pose ].first->add_constraint(new core::scoring::constraints::ResidueTypeConstraint(
				*poses[ current_pose ].first, //pose
				corr_pairs[ current_pose ][ k ], //seqpos
				poses[ i ].first->residue_type( corr_pairs[ i ][ k ] ).name3(), //AAname
				(poses[ current_pose ].second ? positive_penalty : negative_penalty ) //favor native bonus
				) );
			if ( option[ msd::debug ] ) {
				TR << "linking position " << corr_pairs[ current_pose ][ k ] << " on pose " << current_pose << " to position " << corr_pairs[ i ][ k ] << " on pose " << i << std::endl;
				TR << "linking residue " << poses[ current_pose ].first->residue_type( corr_pairs[ current_pose ][ k ]).name3() << " to residue " << poses[ i ].first->residue_type( corr_pairs[ i ][ k ] ).name3() << std::endl;
			}
		}
	}
}


utility::vector1< core::Size >
parse_resfile ( core::pack::task::PackerTaskCOP design_task, core::pose::Pose & pose )
{
	utility::vector1< core::Size > vector;
	utility::vector1<bool> designing = design_task->designing_residues();
	core::pose::PDBInfoCOP pdb_info = pose.pdb_info();
	for ( core::Size i = 1; i <= designing.size(); ++i ) {
		if ( designing[ i ] ) {
			vector.push_back( i );
		}
	}
	return vector;
}

utility::vector1< utility::vector1< core::Size > >
parse_resfiles ( utility::vector1< core::pack::task::TaskFactoryOP > task_factories, utility::vector1< std::pair < core::pose::PoseOP, bool > > poses )
{
	utility::vector1< utility::vector1< core::Size > > return_vector;
	core::Size designable_residues = -1;
	for ( core::Size i = 1; i <= task_factories.size(); ++i ) {
		core::pack::task::PackerTaskOP design_task = task_factories[i]->create_task_and_apply_taskoperations( *(poses[i].first) );
		return_vector.push_back( parse_resfile( design_task, *( poses[1].first ) ) );
		if ( i == 1 )  designable_residues = return_vector[ 1 ].size();
		else {
			if ( return_vector[ i ].size() != designable_residues ) {
				utility_exit_with_message( "All resfiles must have the same number of designable residues");
			}
		}
	}
	return return_vector;
}

void
update_constraints( utility::vector1< std::pair< core::pose::PoseOP, bool > > poses,
	utility::vector1< std::pair< core::pose::PoseOP, bool > > reference_poses,
	utility::vector1< utility::vector1< core::Size > > res_links,
	core::Real positive_linked_penalty,
	core::Real negative_linked_penalty,
	core::Real fnr_penalty,
	core::Size current_pose) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;

	// Need to remove res type linking constraints between each round, since residue ids have changed thru design

	ConstraintSetOP cst_set = poses[ current_pose ].first->constraint_set()->clone();
	utility::vector1< ConstraintCOP > constraints = cst_set->get_all_constraints();
	for ( core::Size i = 1; i <= constraints.size(); ++i ) {
		if ( constraints[ i ]->score_type() == core::scoring::res_type_constraint ) {
			cst_set->remove_constraint( constraints[i], 0 );
		}
	}
	poses[ current_pose ].first->constraint_set( cst_set );


	add_fnr_constraints( poses[ current_pose ].first, reference_poses[ current_pose ].first, res_links[ current_pose ], fnr_penalty, poses[ current_pose ].second );
	apply_linked_constraints( poses, res_links, positive_linked_penalty, negative_linked_penalty, current_pose );

}

utility::vector1< core::Size >
parse_chains ( std::string chains ) {
	utility::vector1< core::Size > ret;
	for ( core::Size i = 1; i <= chains.size(); ++i ) {
		ret.push_back( static_cast<core::Size>( atoi( chains.substr( i-1, 1 ).c_str() ) ) );
	}
	return ret;
}

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		//Add local options
		option.add( msd::positive_states, "Positive states for MSD").def("");
		option.add( msd::negative_states, "Negative states for MSD").def("");
		option.add( msd::upper_chains, "Upper chains of the interface (i.e. 12). Only needed if restricting design or minimization to interface").def("34");
		option.add( msd::lower_chains, "Lower chains of the interface (i.e. 34). Only needed if restricting design or minimization to interface").def("12");
		option.add( msd::n_iterations, "Number of iterations for negative design").def(5);
		option.add( msd::find_consensus, "Find consensus mutations at the end of the protocol" ).def(1);
		option.add( msd::debug, "Show debug messages").def(0);
		option.add( msd::design_rounds, "Rounds of design to perform" ).def(2);
		option.add( msd::ramp, "Weights to enforce sequence convergence for rounds of design" ).def("0.5, 1.0, 2.0");
		option.add( msd::favor_native, "Bonus for favor native residue" ).def(0.5);
		option.add( msd::restrict_design_to_interface, "Restrict design to interface vector" ).def(0);
		option.add( msd::restrict_min_to_interface, "Restrict minimization to interface vector" ).def(0);
		option.add( msd::soft_design, "Perform initial round of design with soft repulsive force" ).def(1);
		option.add( msd::fixed_bb, "Fixed backbone design (no backbone minimization)").def(0);
		option.add( msd::multiple_resfiles, "Option to input a unique resfile for each input structure (must have same number of resfiles as states)");
		option.add( msd::multiple_constraints, "Option to input a unique constraint file for each input structure (must have same number of constraint files as states)");
		devel::init( argc, argv );

		bool debug = option[ msd::debug ];
		utility::vector1< core::Size > upper_chains = parse_chains ( option[ msd::upper_chains ] );
		utility::vector1< core::Size > lower_chains = parse_chains ( option[ msd::lower_chains ] );

		utility::vector1< std::string > pdb_positive, pdb_negative;
		if ( option[ msd::positive_states ].user() ) {
			pdb_positive = option[ msd::positive_states ];
		}
		if ( option[ msd::negative_states ].user() ) {
			pdb_negative = option[ msd::negative_states ];
		}
		if ( pdb_positive.size() + pdb_negative.size() < 2 ) {
			utility_exit_with_message( "Expected at least two pdbs to be specified from the -msd:positive_states and -msd:negative_states flags" );
		}


		//Poses plus bool representing whether it's a positive state
		utility::vector1<std::pair < core::pose::PoseOP, bool > > poses;

		for ( core::Size i = 1; i <= pdb_positive.size(); ++i ) {
			poses.push_back(std::make_pair( core::import_pose::pose_from_file( pdb_positive[ i ], false), 1 ) , core::import_pose::PDB_file);
		}
		for ( core::Size i = 1; i <= pdb_negative.size(); ++i ) {
			poses.push_back(std::make_pair( core::import_pose::pose_from_file( pdb_negative[ i ], false), 0 ) , core::import_pose::PDB_file);
		}

		// Get multiple resfiles if these are provided, or single resfile
		utility::vector1< std::string > resfiles;
		if ( option[ msd::multiple_resfiles ].user() ) {
			resfiles = option[ msd::multiple_resfiles ]();
			// Sanity check for resfiles
			if ( resfiles.size() != poses.size() ) {
				utility_exit_with_message( "Must have same number of resfiles as input structures" );
			}
		} else {
			resfiles = option[ packing::resfile ]();
		}

		// Initialize user specified constraints

		// Add constraints to pose if unique files are specified - careful as these will supercede normal command line constraints
		if ( option[ msd::multiple_constraints ].user() ) {
			utility::vector1< std::string > constraint_files = option[ msd::multiple_constraints ]();
			// Sanity check for constraint files
			if ( constraint_files.size() != poses.size() ) {
				utility_exit_with_message( "Must have same number of constraint files as input structures" );
			}
			for ( core::Size current_pose = 1; current_pose <= constraint_files.size(); ++current_pose ) {
				core::scoring::constraints::ConstraintSetOP cstset_ = core::scoring::constraints::ConstraintIO::get_instance()->read_constraints(
					constraint_files[ current_pose ], new core::scoring::constraints::ConstraintSet, *poses[ current_pose ].first );
				poses[ current_pose ].first->constraint_set( cstset_ );
			}

		} else {
			for ( core::Size current_pose = 1; current_pose <= poses.size(); ++current_pose ) {
				// Since there's an internal check if the constraints flag is passed I don't need to check here
				core::scoring::constraints::add_constraints_from_cmdline_to_pose( *poses[ current_pose ].first );
				core::scoring::constraints::add_fa_constraints_from_cmdline_to_pose( *poses[ current_pose ].first );
			}
		}

		core::Size nstruct = 1;
		if ( option[ out::nstruct ].user() ) {
			nstruct = option[ out::nstruct ];
		}

		std::string min_type = "lbfgs_armijo_nonmonotone";
		if ( option[ run::min_type ].user() ) {
			min_type = option[ run::min_type];
		}

		// Get weights for ramping
		utility::vector1<core::Size> ramp = utility::string_split( option[ msd::ramp ], ',', core::Real() );
		if ( ramp.size() != ( option[ msd::soft_design ] ? (core::Size) option[ msd::design_rounds ]+1 : (core::Size) option[ msd::design_rounds ] ) ) {
			utility_exit_with_message( "There must be the same number of ramping weights as design rounds" );
		}


		/* Set up task operations and assign to the right task factory */
		core::pack::task::operation::InitializeFromCommandlineOP ifcl = new core::pack::task::operation::InitializeFromCommandline;
		core::pack::task::operation::RestrictToRepackingOP rtr = new core::pack::task::operation::RestrictToRepacking;

		utility::vector1< core::pack::task::operation::ReadResfileOP > rrf_vector;
		for ( core::Size i = 1; i <= poses.size(); ++i ) {
			if ( resfiles.size() > 1 ) {
				rrf_vector.push_back( new core::pack::task::operation::ReadResfile( resfiles[ i ] ) );
			} else {
				rrf_vector.push_back( new core::pack::task::operation::ReadResfile );
			}
		}


		protocols::toolbox::task_operations::RestrictToInterfaceVectorOperationOP rtiv = new protocols::toolbox::task_operations::RestrictToInterfaceVectorOperation;
		rtiv->CB_dist_cutoff(10);
		rtiv->nearby_atom_cutoff(5.5);
		rtiv->vector_angle_cutoff(75);
		rtiv->vector_dist_cutoff(9.0);
		rtiv->upper_chain(upper_chains);
		rtiv->lower_chain(lower_chains);

		utility::vector1< core::pack::task::TaskFactoryOP > design_task_factories;
		for ( core::Size i = 1; i <= rrf_vector.size(); ++i ) {
			core::pack::task::TaskFactoryOP design_task_factory = new core::pack::task::TaskFactory;
			design_task_factory->push_back( ifcl );
			if ( option[ msd::restrict_design_to_interface ] ) design_task_factory->push_back( rtiv );
			design_task_factory->push_back( rrf_vector[i] );
			design_task_factories.push_back( design_task_factory );
		}

		core::pack::task::TaskFactoryOP packing_min_task_factory = new core::pack::task::TaskFactory;
		packing_min_task_factory->push_back( ifcl );
		if ( option[ msd::restrict_min_to_interface ] ) packing_min_task_factory->push_back( rtiv );
		packing_min_task_factory->push_back( rtr );

		core::pack::task::TaskFactoryOP min_task_factory = new core::pack::task::TaskFactory;
		if ( option[ msd::restrict_min_to_interface ] ) min_task_factory->push_back( rtiv );

		/* Create score functions */
		core::scoring::ScoreFunctionOP sfxn_clean = core::scoring::getScoreFunction();
		core::scoring::ScoreFunctionOP sfxn_design = sfxn_clean->clone();
		core::scoring::ScoreFunctionOP soft_rep_clean = core::scoring::ScoreFunctionFactory::create_score_function("soft_rep");
		core::scoring::ScoreFunctionOP soft_rep_design = soft_rep_clean->clone();

		// Set up constraint weights for scoring functions
		soft_rep_design->set_weight( core::scoring::res_type_constraint, 1.0 );
		sfxn_design->set_weight( core::scoring::res_type_constraint, 1.0 );
		if ( option[ constraints::cst_file ].user() ) {
			core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *soft_rep_clean );
			core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *soft_rep_design );
			core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *sfxn_clean );
			core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *sfxn_design );
		} else if ( option[ constraints::cst_fa_file ].user() ) {
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *soft_rep_clean );
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *soft_rep_design );
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *sfxn_clean );
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *sfxn_design );
		}

		/* Create movers to do the design */
		utility::vector1< protocols::simple_moves::PackRotamersMoverOP > soft_des_vector;
		for ( core::Size i = 1; i <= rrf_vector.size(); ++i ) {
			protocols::simple_moves::PackRotamersMoverOP soft_des = new protocols::simple_moves::PackRotamersMover;
			soft_des->score_function( soft_rep_design );
			soft_des->task_factory( design_task_factories[i] );
			soft_des_vector.push_back( soft_des );
		}

		utility::vector1< protocols::simple_moves::NegativePackRotamersMoverOP > neg_design_vector;
		for ( core::Size i = 1; i <= rrf_vector.size(); ++i ) {
			protocols::simple_moves::NegativePackRotamersMoverOP neg_design = new protocols::simple_moves::NegativePackRotamersMover;
			neg_design->score_function ( sfxn_design );
			neg_design->task_factory ( design_task_factories[i] );
			neg_design->n_iterations( option[ msd::n_iterations ] );
			neg_design_vector.push_back( neg_design );
		}

		utility::vector1< protocols::simple_moves::PackRotamersMoverOP > pos_design_vector;
		for ( core::Size i = 1; i <= rrf_vector.size(); ++i ) {
			protocols::simple_moves::PackRotamersMoverOP pos_design = new protocols::simple_moves::PackRotamersMover;
			pos_design->score_function( sfxn_design );
			pos_design->task_factory( design_task_factories[ i ] );
			pos_design_vector.push_back( pos_design );
		}


		/* Create pack rotamer movers to use during minimization */
		protocols::simple_moves::PackRotamersMoverOP soft_rp = new protocols::simple_moves::PackRotamersMover;
		soft_rp->score_function( soft_rep_clean );
		soft_rp->task_factory( packing_min_task_factory );

		protocols::simple_moves::PackRotamersMoverOP rp = new protocols::simple_moves::PackRotamersMover;
		rp->score_function( sfxn_clean );
		rp->task_factory( packing_min_task_factory );

		//create task aware minimization movers to do minimization
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;

		protocols::simple_moves::MinMoverOP soft_min_mover  = new protocols::simple_moves::MinMover;
		soft_min_mover->min_type( min_type );
		soft_min_mover->tolerance( 0.001 );
		soft_min_mover->score_function( soft_rep_clean );

		protocols::simple_moves::TaskAwareMinMoverOP soft_min = new protocols::simple_moves::TaskAwareMinMover(
			soft_min_mover, min_task_factory );
		soft_min->chi(1);
		soft_min->bb(1);
		soft_min->jump(1);

		protocols::simple_moves::MinMoverOP soft_min_sc_mover = new protocols::simple_moves::MinMover;
		soft_min_sc_mover->min_type( min_type );
		soft_min_mover->tolerance( 0.001 );
		soft_min_mover->score_function( soft_rep_clean );

		protocols::simple_moves::TaskAwareMinMoverOP soft_min_sc = new protocols::simple_moves::TaskAwareMinMover(
			soft_min_sc_mover, min_task_factory );
		soft_min_sc->chi(1);
		soft_min_sc->bb(0);
		soft_min_sc->jump(0);

		protocols::simple_moves::MinMoverOP min_sc_mover = new protocols::simple_moves::MinMover;
		min_sc_mover->min_type( min_type );
		soft_min_mover->tolerance( 0.001 );
		soft_min_mover->score_function( sfxn_clean );

		protocols::simple_moves::TaskAwareMinMoverOP min_sc = new protocols::simple_moves::TaskAwareMinMover(
			min_sc_mover, min_task_factory );
		min_sc->chi(1);
		min_sc->bb(0);
		min_sc->jump(0);

		protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover;
		min_mover->min_type( min_type );
		min_mover->tolerance( 0.001 );
		min_mover->score_function( sfxn_clean );

		protocols::simple_moves::TaskAwareMinMoverOP min = new protocols::simple_moves::TaskAwareMinMover(
			min_mover, min_task_factory );
		min->chi(1);
		min->bb(1);
		min->jump(1);

		utility::vector1< utility::vector1< core::Size > > residue_links = parse_resfiles( design_task_factories, poses );

		// Where the magic happens

		utility::vector1< std::pair < core::pose::PoseOP, bool > > modifiable_poses;
		for ( core::Size n = 1; n <= nstruct; ++n ) {
			modifiable_poses.clear();

			/* Create copy of each pose so original list stays unchanged */
			for ( core::Size i = 1; i <= poses.size(); ++i ) {
				modifiable_poses.push_back( std::make_pair( poses[i].first->clone(), poses[i].second ) );
			}

			/* Perform round of soft design if requested */
			if ( option [ msd::soft_design ] ) {
				for ( core::Size i = 1; i <= modifiable_poses.size(); ++i ) {
					if ( debug ) TR << "applying constraints round 1 for pose " << i << std::endl;
					update_constraints( modifiable_poses, poses, residue_links, ramp[ 1 ], -1.0, option[ msd::favor_native ], i);

					if ( modifiable_poses[ i ].second ) {
						if ( debug ) TR << "design round 1 - pose " << i << std::endl;
						soft_des_vector[ i ]->apply( *modifiable_poses[ i ].first );
						//add minimize parsed protocol mover
						if ( !option[ msd::fixed_bb] ) {
							if ( debug ) TR << "minimize round 1 - pose " << i << std::endl;
							if ( debug ) TR << "minimize softrp round 1 - pose " << i << std::endl;
							soft_rp->apply( *modifiable_poses[ i ].first );
							if ( debug ) TR << "minimize softminsc round 1 - pose " << i << std::endl;
							soft_min_sc->apply( *modifiable_poses[ i ].first );
							if ( debug ) TR << "minimize softrp round 1 - pose " << i << std::endl;
							soft_rp->apply( *modifiable_poses[ i ].first );
							if ( debug ) TR << "minimize softmin round 1 - pose " << i << std::endl;
							soft_min->apply( *modifiable_poses[ i ].first );
							if ( debug ) TR << "minimize softrp round 1 - pose " << i << std::endl;
							soft_rp->apply( *modifiable_poses[ i ].first );
							if ( debug ) TR << "minimize minsc round 1 - pose " << i << std::endl;
							min_sc->apply( *modifiable_poses[ i ].first );
							if ( debug ) TR << "minimize softrp round 1 - pose " << i << std::endl;
							soft_rp->apply( *modifiable_poses[ i ].first );
							if ( debug ) TR << "minimize min round 1 - pose " << i << std::endl;
							min->apply( *modifiable_poses[ i ].first );
							if ( debug ) TR << "minimize rp round 1 - pose " << i << std::endl;
							rp->apply( *modifiable_poses[ i ].first );
							if ( debug ) TR << "minimize minsc round 1 - pose " << i << std::endl;
							min_sc->apply( *modifiable_poses[ i ].first );
							if ( debug ) TR << "minimize rp round 1 - pose " << i << std::endl;
							rp->apply( *modifiable_poses[ i ].first );
							if ( debug ) TR << "minimize min round 1 - pose " << i << std::endl;
							min->apply( *modifiable_poses[ i ].first );
						}
					} else {
						neg_design_vector[ i ]->apply( *modifiable_poses[ i ].first );
					}
				}
			}

			for ( core::Size j = 1; j <= (core::Size) option[ msd::design_rounds ]; ++j ) {
				/* Perform subsequent rounds of design */
				for ( core::Size i = 1; i <= modifiable_poses.size(); ++i ) {
					core::Real linking_wt = ( option[ msd::soft_design ] ? ramp[ j+1 ] : ramp[ j ] );
					if ( debug ) TR << "applying constraints round "<< j << " for pose " << i << std::endl;
					update_constraints( modifiable_poses, poses, residue_links, linking_wt, -2, option[ msd::favor_native ], i);

					if ( modifiable_poses [ i ].second ) {
						if ( debug ) TR << "design round " << j << " - pose " << i << std::endl;
						pos_design_vector[ i ]->apply( *modifiable_poses[ i ].first );
						//add minimize parsed protocol mover
						if ( !option[ msd::fixed_bb ] ) {
							if ( debug ) TR << "minimize round " << j << " - pose " << i << std::endl;
							soft_rp->apply( *modifiable_poses[ i ].first );
							soft_min_sc->apply( *modifiable_poses[ i ].first );
							soft_rp->apply( *modifiable_poses[ i ].first );
							soft_min->apply( *modifiable_poses[ i ].first );
							soft_rp->apply( *modifiable_poses[ i ].first );
							min_sc->apply( *modifiable_poses[ i ].first );
							soft_rp->apply( *modifiable_poses[ i ].first );
							min->apply( *modifiable_poses[ i ].first );
							rp->apply( *modifiable_poses[ i ].first );
							min_sc->apply( *modifiable_poses[ i ].first );
							rp->apply( *modifiable_poses[ i ].first );
							min->apply( *modifiable_poses[ i ].first );
						}
					} else {
						neg_design_vector[ i ]->apply( *modifiable_poses[ i ].first );
					}


				}
			}


			if ( option[ msd::find_consensus ] ) {
				TR << "finding consensus mutations" << std::endl;
				/* Revert mutations to find the best aa at each spot */
				core::Real scaling_factor = ( (core::Real) pdb_positive.size() )/pdb_negative.size();
				for ( core::Size i = 1; i <= residue_links[ 1 ].size(); ++i ) {
					core::Size seqpos1 = residue_links[ 1 ][ i ];
					//check to see if they are all the same
					bool diff = false;
					for ( core::Size j = 2; j <= modifiable_poses.size(); ++j ) {
						core::Size seqpos2 = residue_links[ j ][ i ];
						if ( modifiable_poses[ j ].first->residue( seqpos2 ).aa() != modifiable_poses[ 1 ].first->residue( seqpos1 ).aa() ) {
							diff = true;
							break;
						}
					}
					if ( diff ) {
						core::Size min_index = 0;
						core::Real min_score = 0;
						utility::vector1< core::Real > scores ( modifiable_poses.size(), 0 );
						for ( core::Size ref = 1; ref <= modifiable_poses.size(); ++ref ) {
							for ( core::Size comp = 1; comp <= modifiable_poses.size(); ++comp ) {
								//if residues at this position are same
								if ( modifiable_poses[ ref ].first->residue( residue_links[ ref ][ i ] ).aa() ==
										modifiable_poses[ comp ].first->residue( residue_links[ comp ][ i ] ).aa() ) {
									scores[ ref ] += (modifiable_poses[ comp ].second ?
										modifiable_poses[ comp ].first->energies().total_energy() :
										std::min( scaling_factor*-modifiable_poses[ comp ].first->energies().total_energy(), 0.0) );
									continue;
								} else { //if they're different

									core::pose::PoseOP mutpose = modifiable_poses[ comp ].first->clone();
									mutpose->replace_residue( residue_links[ comp ][ i ], modifiable_poses[ ref ].first->residue( residue_links[ ref ][ i ] ), true );
									rp->apply( *mutpose );
									scores[ ref ] += (modifiable_poses[ comp ].second ?
										mutpose->energies().total_energy() :
										std::min( scaling_factor*-mutpose->energies().total_energy(), 0.0 ) );
								}

							}
							if ( scores[ ref ] < min_score || ref == 1 ) {
								min_index = ref;
								min_score = scores[ ref ];
							}
						}
						for ( core::Size pose = 1; pose <= modifiable_poses.size(); ++pose ) {
							if ( modifiable_poses[ pose ].first->residue( residue_links[ pose ][ i ] ).aa() !=
									modifiable_poses[ min_index ].first->residue( residue_links[ min_index ][ i ] ).aa() ) {
								modifiable_poses[ pose ].first->replace_residue( residue_links[ pose ][ i ],
									modifiable_poses[ min_index ].first->residue( residue_links[ min_index ][ i ] ),
									true );
								rp->apply( *modifiable_poses[ pose ].first );
							}
						}
					}
				}
			}


			/* Finished - dump out the scored pdbs to the output file */
			core::Size j = 1;
			for ( core::Size i = 1; i <= pdb_positive.size(); ++i ) {
				std::string file_name = pdb_positive[ i ];
				file_name = file_name.substr(0, file_name.length()-4);
				std::string out_name = (option[ out::suffix ]()=="") ?
					file_name + "_" + get_out_number(n) + ".pdb" :
					file_name + "_" + option[ out::suffix ]() + "_" + get_out_number(n) + ".pdb";
				modifiable_poses[ j ].first->dump_scored_pdb( out_name, *sfxn_clean );
				j++;
			}
			for ( core::Size i = 1; i <= pdb_negative.size(); ++i ) {
				std::string file_name = pdb_negative[i];
				file_name = file_name.substr(0, file_name.length()-4);
				std::string out_name = (option[ out::suffix ]()=="") ?
					file_name + "_" + get_out_number(n) + ".pdb" :
					file_name + "_" + option[ out::suffix ]() + "_" + get_out_number(n) + ".pdb";
				modifiable_poses[ j ].first->dump_scored_pdb( out_name, *sfxn_clean );
				j++;
			}

		}
	} catch (utility::excn::Exception const & e ) {
		utility_exit_with_message("caught exception " + e.msg());
	}
}
