// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/selection.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/VariantType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>

// Mover headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/oop/OopRandomPuckMover.hh>
#include <protocols/simple_moves/oop/OopPuckMover.hh>
#include <protocols/simple_moves/oop/OopMover.hh>
#include <protocols/simple_moves/oop/OopRandomSmallMover.hh>
#include <protocols/simple_moves/oop/OopPatcher.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/ncbb/oop/OopCreatorMover.hh>
#include <protocols/ncbb/oop/OopCreatorMoverCreator.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>

#include <numeric/conversions.hh>

//Basic headers
#include <basic/resource_manager/ResourceManager.hh>


// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/tag/Tag.hh>

// C++ headers
#include <string>
#include <sstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::rigid;
using namespace protocols::toolbox;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::chiral;
using namespace protocols::simple_moves::oop;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;


// tracer - used to replace cout
static THREAD_LOCAL basic::Tracer TR( "OOP_Creator" );

namespace protocols {
namespace ncbb {
namespace oop {

OopCreatorMover::OopCreatorMover():
	oop_plus_positions_ (utility::tools::make_vector1(0)),
	oop_minus_positions_ (utility::tools::make_vector1(0)),
	oop_d_plus_positions_ (utility::tools::make_vector1(0)),
	oop_d_minus_positions_ (utility::tools::make_vector1(0)),
	oop_low_e_puck_positions_ (utility::tools::make_vector1(0)),
	prepend_n_residues_ (0),
	append_n_residues_ (0),
	final_repack_ (false),
	final_minimize_ (false),
	final_mc_ (false),
	final_correct_oop_post_ (false)
{
	Mover::type("OopCreatorMover");
}

OopCreatorMover::OopCreatorMover(
	utility::vector1<core::Size> oop_plus_positions,
	utility::vector1<core::Size> oop_minus_positions,
	utility::vector1<core::Size> oop_d_plus_positions,
	utility::vector1<core::Size> oop_d_minus_positions,
	utility::vector1<core::Size> oop_low_e_puck_positions,
	core::Size prepend_n_residues,
	core::Size append_n_residues,
	bool final_repack,
	bool final_minimize,
	bool final_mc,
	bool final_correct_oop_post
): Mover("OopCreatorMover"),
	oop_plus_positions_ (oop_plus_positions),
	oop_minus_positions_ (oop_minus_positions),
	oop_d_plus_positions_ (oop_d_plus_positions),
	oop_d_minus_positions_ (oop_d_minus_positions),
	oop_low_e_puck_positions_ (oop_low_e_puck_positions),
	prepend_n_residues_ (prepend_n_residues),
	append_n_residues_ (append_n_residues),
	final_repack_ (final_repack),
	final_minimize_ (final_minimize),
	final_mc_ (final_mc),
	final_correct_oop_post_ (final_correct_oop_post)
{
}

void
OopCreatorMover::apply(
	core::pose::Pose & pose
)
{

	// create score function
	//kdrew: old standard scoring function, using MM scoring function now because of NCAAs
	//scoring::ScoreFunctionOP score_fxn( get_score_function() );
	scoring::ScoreFunctionOP score_fxn( ScoreFunctionFactory::create_score_function( scoring::MM_STD_WTS) );
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

	scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

	utility::vector1< Size > all_positions;

	utility::vector1< Size > const plus_positions = oop_plus_positions_;//option[ oop_creator::oop_plus_positions].value();
	all_positions.insert( all_positions.end(), plus_positions.begin(), plus_positions.end() );
	for ( Size i = 1; i <= plus_positions.size(); i++ ) {
		OopPatcherOP oop_patcher( new OopPatcher( plus_positions[i] ) );
		oop_patcher->apply( pose );

		OopPuckPlusMoverOP opm_plus( new OopPuckPlusMover( plus_positions[i] ) );
		opm_plus->apply( pose );
	}

	utility::vector1< Size > const minus_positions = oop_minus_positions_;//option[ oop_creator::oop_minus_positions].value();
	all_positions.insert( all_positions.end(), minus_positions.begin(), minus_positions.end() );
	for ( Size i = 1; i <= minus_positions.size(); i++ ) {
		OopPatcherOP oop_patcher( new OopPatcher( minus_positions[i] ) );
		oop_patcher->apply( pose );

		OopPuckMinusMoverOP opm_minus( new OopPuckMinusMover( minus_positions[i] ) );
		opm_minus->apply( pose );
	}

	utility::vector1< Size > const d_plus_positions = oop_d_plus_positions_;//option[ oop_creator::oop_d_plus_positions].value();
	all_positions.insert( all_positions.end(), d_plus_positions.begin(), d_plus_positions.end() );
	for ( Size i = 1; i <= d_plus_positions.size(); i++ ) {
		OopPatcherOP oop_patcher( new OopPatcher( d_plus_positions[i] ) );
		oop_patcher->apply( pose );

		OopDPuckPlusMoverOP opm_plus( new OopDPuckPlusMover( d_plus_positions[i] ) );
		opm_plus->apply( pose );
	}

	utility::vector1< Size > const d_minus_positions = oop_d_minus_positions_;//option[ oop_creator::oop_d_minus_positions].value();
	all_positions.insert( all_positions.end(), d_minus_positions.begin(), d_minus_positions.end() );
	for ( Size i = 1; i <= d_minus_positions.size(); i++ ) {
		OopPatcherOP oop_patcher( new OopPatcher( d_minus_positions[i] ) );
		oop_patcher->apply( pose );

		OopDPuckMinusMoverOP opm_minus( new OopDPuckMinusMover( d_minus_positions[i] ) );
		opm_minus->apply( pose );
	}

	utility::vector1< Size > const low_e_puck_positions = oop_low_e_puck_positions_;//option[ oop_creator::oop_low_e_puck_positions ].value();
	all_positions.insert( all_positions.end(), low_e_puck_positions.begin(), low_e_puck_positions.end() );
	for ( Size i = 1; i <= low_e_puck_positions.size(); ++i ) {
		OopPatcherOP oop_patcher( new OopPatcher( low_e_puck_positions[i] ) );
		oop_patcher->apply( pose );
		//kdrew: poor man's version for now, if L use PuckPlus, if D use DPuckPlus
		if ( pose.residue_type( low_e_puck_positions[i] ).is_l_aa() ) {
			//kdrew: use PuckPlus
			OopPuckPlusMoverOP opm_plus( new OopPuckPlusMover( low_e_puck_positions[i] ) );
			opm_plus->apply( pose );
		} else if ( pose.residue_type( low_e_puck_positions[i] ).is_d_aa() ) {
			//kdrew: use DPuckPlus
			OopDPuckPlusMoverOP opm_plus( new OopDPuckPlusMover( low_e_puck_positions[i] ) );
			opm_plus->apply( pose );
		} else {
			TR << " residue: " << pose.residue_type( low_e_puck_positions[i] ).name() << " not found in chiral map" <<  std::endl;
			TR << " possibly achiral (ex GLY) or not listed in map" <<  std::endl;
			TR << " not changing" << std::endl;
		}
	}

	//kdrew: sets oop_post phi/psi near low energy well
	if ( final_correct_oop_post_ ) { //option[ oop_creator::correct_oop_post ].value() )
		for ( Size i = 1; i <= all_positions.size(); ++i ) {
			//kdrew: the +1 is to get the oop_post position
			if ( pose.residue_type( all_positions[i] +1 ).is_d_aa() ) {
				pose.set_phi( all_positions[i] +1, 135.0 ) ;
				pose.set_psi( all_positions[i] +1, -70.0 ) ;
			} else {
				pose.set_phi( all_positions[i] +1, -135.0 ) ;
				pose.set_psi( all_positions[i] +1, 70.0 ) ;
			}
		}

	}


	//kdrew: create glycine residue
	ResidueOP gly( ResidueFactory::create_residue( *core::pose::get_restype_for_pose( pose, "GLY" ) ) );

	Size pep_begin( pose.conformation().chain_begin( 1 ) );
	Size pep_end( pose.conformation().chain_end( 1 ) );

	//kdrew: since we probably added new connection types (i.e. oop CYP and CZP atoms) above, need to reset connections
	pose.conformation().detect_bonds();
	pose.conformation().detect_pseudobonds();
	for ( Size i=1; i<=pose.size(); ++i ) {
		pose.conformation().update_polymeric_connection(i);
	}

	//kdrew: grabbed code from chrisk pep_spec
	//kdrew: append residues , hard coded to glycine
	for ( Size i = 1; i <= append_n_residues_ /*Size( option[ oop_creator::append_n_residues ].value() ) */; ++i ) {
		TR << "in append: " << pep_end << std::endl;
		pose.conformation().safely_append_polymer_residue_after_seqpos( *gly, pep_end, true );
		pep_end = pep_end + 1;
		pose.set_omega( pep_end - 1, 180.0 );
		pose.conformation().update_polymeric_connection( pep_end );
		pose.conformation().update_polymeric_connection( pep_end - 1 );
	}
	//kdrew: prepend residues , hard coded to glycine
	for ( Size i = 1; i <= append_n_residues_ /*Size( option[ oop_creator::prepend_n_residues ].value() )*/ ; ++i ) {
		TR << "in prepend: " << pep_begin << std::endl;
		pose.conformation().safely_prepend_polymer_residue_before_seqpos( *gly, pep_begin, true );
		pep_end = pep_end + 1;
		pep_begin =  pose.conformation().chain_begin( 1 ) ; //reset pep beginning
		pose.set_omega( pep_begin, 180.0 );
		pose.conformation().update_polymeric_connection( pep_begin );
		pose.conformation().update_polymeric_connection( pep_begin + 1 );
	}


	//kdrew: monte carlo phi/psi of oop to find low energy
	if ( final_mc_ ) { //option[ oop_creator::final_mc ].value() )
		moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
		moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, 0.2 ) );

		kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
		simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, 0.2, 1 ) );
		pert_pep_small->angle_max( 'H', 2.0 );
		pert_pep_small->angle_max( 'L', 2.0 );
		pert_pep_small->angle_max( 'E', 2.0 );

		utility::vector1< Size > oop_pre_positions;
		//kdrew: load all oop_pre positions into vector and make all non-oop_pre positions movable by small mover
		for ( Size i = 1; i <= pose.size(); ++i ) {
			TR << "resid: " << i << " is OOP_PRE: " << pose.residue(i).has_variant_type(chemical::OOP_PRE) << std::endl;

			if ( pose.residue(i).has_variant_type(chemical::OOP_PRE) != 1 ) {
				if ( pose.residue_type( i ).is_l_aa() ) {
					TR << "setting small movable resid: "<< i<<std::endl;
				} else {
					oop_pre_positions.push_back(i);
				}
			}

		}

		pert_sequence->add_mover( pert_pep_small );

		//kdrew: add all oop_pre positions to random small mover
		if ( oop_pre_positions.size() > 0 ) {
			OopRandomSmallMoverOP opm( new OopRandomSmallMover ( oop_pre_positions, 2.0 ) );
			moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( opm, oop_pre_positions.size() * 1000 ) );
			pert_sequence->add_mover( pert_pep_repeat );
		}
		moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

		pert_trial->apply( pose );
		pert_mc->recover_low( pose );
	}

	if ( final_repack_ ) {
		// create a task factory and task operations
		using core::pack::task::operation::TaskOperationCOP;
		TaskFactoryOP tf( new TaskFactory() );
		tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );

		using namespace basic::resource_manager;
		if ( ResourceManager::get_instance()->has_option( packing::resfile ) ||  option[ packing::resfile ].user() ) {
			operation::ReadResfileOP rrop( new operation::ReadResfile() );
			rrop->default_filename();
			tf->push_back( rrop );
		} else {
			//kdrew: do not do design, makes NATAA if res file is not specified
			operation::RestrictToRepackingOP rtrp( new operation::RestrictToRepacking() );
			tf->push_back( rtrp );
		}

		// create a pack rotamers mover
		simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover() );
		packer->task_factory( tf );
		packer->score_function( score_fxn );
		packer->apply(pose);
	}

	//kdrew: monte carlo phi/psi of oop to find low energy
	if ( final_mc_ ) {
		moves::SequenceMoverOP pert_sequence( new moves::SequenceMover() );
		moves::MonteCarloOP pert_mc( new moves::MonteCarlo( pose, *score_fxn, 0.2 ) );

		kinematics::MoveMapOP pert_pep_mm( new kinematics::MoveMap() );
		simple_moves::SmallMoverOP pert_pep_small( new simple_moves::SmallMover( pert_pep_mm, 0.2, 1 ) );
		pert_pep_small->angle_max( 'H', 2.0 );
		pert_pep_small->angle_max( 'L', 2.0 );
		pert_pep_small->angle_max( 'E', 2.0 );

		utility::vector1< Size > oop_pre_positions;
		//kdrew: load all oop_pre positions into vector and make all non-oop_pre positions movable by small mover
		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( pose.residue(i).has_variant_type(chemical::OOP_PRE) != 1 ) {
				pert_pep_mm->set_bb( i );
			} else {
				oop_pre_positions.push_back(i);
			}
		}
		pert_sequence->add_mover( pert_pep_small );

		//kdrew: add all oop_pre positions to random small mover
		if ( oop_pre_positions.size() > 0 ) {
			OopRandomSmallMoverOP opm( new OopRandomSmallMover ( oop_pre_positions, 2.0 ) );
			moves::RepeatMoverOP pert_pep_repeat( new moves::RepeatMover( opm, oop_pre_positions.size() * 1000 ) );
			pert_sequence->add_mover( pert_pep_repeat );
		}
		moves::TrialMoverOP pert_trial( new moves::TrialMover( pert_sequence, pert_mc ) );

		pert_trial->apply( pose );
		pert_mc->recover_low( pose );
	}

	if ( final_repack_ ) {
		// create a task factory and task operations
		using core::pack::task::operation::TaskOperationCOP;
		TaskFactoryOP tf( new TaskFactory() );
		tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );

		using namespace basic::resource_manager;
		if ( ResourceManager::get_instance()->has_option( packing::resfile ) ||  option[ packing::resfile ].user() ) {
			operation::ReadResfileOP rrop( new operation::ReadResfile() );
			rrop->default_filename();
			tf->push_back( rrop );
		} else {
			//kdrew: do not do design, makes NATAA if res file is not specified
			operation::RestrictToRepackingOP rtrp( new operation::RestrictToRepacking() );
			tf->push_back( rtrp );
		}

		// create a pack rotamers mover
		simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover() );
		packer->task_factory( tf );
		packer->score_function( score_fxn );

		packer->apply(pose);
	}


	if ( final_minimize_ ) {
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::func;
		using namespace core::scoring::constraints;
		using namespace numeric::conversions;

		//kdrew: add constraints to omega angle, (this problem might have been fixed and these constraints are unnecessary)
		for ( Size i = 1; i < pose.conformation().chain_end( 1 ); ++i ) {
			AtomID id1,id2,id3,id4;
			TorsionID torsion_id = TorsionID( i, BB, 3 ); //kdrew: 3 is omega angle
			//kdrew: put constraint on omega angle
			pose.conformation().get_torsion_angle_atom_ids( torsion_id, id1, id2, id3, id4 );
			Real torsion_value( pose.torsion( torsion_id ) );
			CircularHarmonicFuncOP circularharm_func( new CircularHarmonicFunc( radians( torsion_value ), radians( 10.0 ) ) );
			ConstraintCOP dihedral1( ConstraintOP( new DihedralConstraint( id1, id2, id3, id4, circularharm_func ) ) );
			pose.add_constraint( dihedral1 );
		}
		//kdrew: if constraint weight is not set on commandline or elsewhere, set to 1.0
		if ( score_fxn->has_zero_weight( dihedral_constraint ) ) {
			score_fxn->set_weight( dihedral_constraint, 1.0 );
		}
		if ( score_fxn->has_zero_weight( atom_pair_constraint ) ) {
			score_fxn->set_weight( atom_pair_constraint, 1.0 );
		}

		// create move map for minimization
		kinematics::MoveMapOP mm( new kinematics::MoveMap() );
		mm->set_bb( true );
		mm->set_chi( true );
		mm->set_jump( 1, true );

		// create minimization mover
		simple_moves::MinMoverOP minM( new protocols::simple_moves::MinMover( mm, score_fxn, option[ OptionKeys::run::min_type ].value(), 0.01, true ) );
		minM->apply( pose );
	}
}

protocols::moves::MoverOP
OopCreatorMover::clone() const
{
	return protocols::moves::MoverOP( new OopCreatorMover () );
}

void
OopCreatorMover::parse_my_tag
(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose
) {
	using namespace utility;

	if ( tag->hasOption( "oop_plus_positions") ) {
		oop_plus_positions_ = core::pose::get_resnum_list(tag, "oop_plus_positions", pose);
	} else {
		oop_plus_positions_ = *(new vector1<core::Size>(0));
	}

	if ( tag->hasOption( "oop_minus_positions") ) {
		oop_minus_positions_ = core::pose::get_resnum_list(tag, "oop_minus_positions", pose);
	} else {
		oop_minus_positions_ = *(new vector1<core::Size>(0));
	}

	if ( tag->hasOption( "oop_d_plus_positions") ) {
		oop_d_plus_positions_ = core::pose::get_resnum_list(tag, "oop_d_plus_positions", pose);
	} else {
		oop_d_plus_positions_ = *(new vector1<core::Size>(0));
	}

	if ( tag->hasOption( "oop_d_minus_positions") ) {
		oop_d_minus_positions_ = core::pose::get_resnum_list(tag, "oop_d_minus_positions", pose);
	} else {
		oop_d_minus_positions_ = *(new vector1<core::Size>(0));
	}

	if ( tag->hasOption( "oop_low_e_puck_positions") ) {
		oop_low_e_puck_positions_ = core::pose::get_resnum_list(tag, "oop_low_e_puck_positions", pose);
	} else {
		oop_low_e_puck_positions_ = *(new vector1<core::Size>(0));
	}

	if ( tag->hasOption( "prepend_n_residues") ) {
		prepend_n_residues_ = tag->getOption<core::Size>("prepend_n_residues", prepend_n_residues_);
	} else {
		prepend_n_residues_ = 0;
	}

	if ( tag->hasOption( "append_n_residues") ) {
		append_n_residues_ = tag->getOption<core::Size>("append_n_residues", append_n_residues_);
	} else {
		append_n_residues_ = 0;
	}

	if ( tag->hasOption( "final_repack") ) {
		final_repack_ = tag->getOption<bool>("final_repack", final_repack_);
	} else {
		final_repack_ = false;
	}

	if ( tag->hasOption( "final_minimize") ) {
		final_minimize_ = tag->getOption<bool>("final_minimize", final_minimize_);
	} else {
		final_minimize_ = false;
	}

	if ( tag->hasOption( "final_mc") ) {
		final_mc_ = tag->getOption<bool>("final_mc", final_mc_);
	} else {
		final_mc_ = false;
	}

	if ( tag->hasOption( "final_correct_oop_post") ) {
		final_correct_oop_post_ = tag->getOption<bool>("final_correct_oop_post", final_correct_oop_post_);
	} else {
		final_correct_oop_post_ = false;
	}
}

// MoverCreator
// XRW TEMP moves::MoverOP
// XRW TEMP OopCreatorMoverCreator::create_mover() const {
// XRW TEMP  return moves::MoverOP( new OopCreatorMover() );
// XRW TEMP }

// XRW TEMP std::string OopCreatorMoverCreator::keyname() const {
// XRW TEMP  return OopCreatorMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string OopCreatorMover::mover_name(){
// XRW TEMP  return "OopCreatorMover";
// XRW TEMP }

std::string OopCreatorMover::get_name() const {
	return mover_name();
}

std::string OopCreatorMover::mover_name() {
	return "OopCreatorMover";
}

void OopCreatorMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "oop_plus_positions", xsct_int_cslist, "Positions where the oop should have a positive pucker" )
		+ XMLSchemaAttribute( "oop_minus_positions", xsct_int_cslist, "Positions where the oop should have a negative pucker" )
		+ XMLSchemaAttribute( "oop_d_plus_positions", xsct_int_cslist, "Positions where the oop has a D amino acid and should have a positive pucker" )
		+ XMLSchemaAttribute( "oop_d_minus_positions", xsct_int_cslist, "Positions where the oop has a D amino acid and should have a negative pucker" )
		+ XMLSchemaAttribute( "oop_low_e_puck_positions", xsct_int_cslist, "Positions where the oop should have the lower energy pucker option" )
		+ XMLSchemaAttribute::attribute_w_default( "prepend_n_residues", xsct_non_negative_integer, "Prepend this many residues to the input pose", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "append_n_residues", xsct_non_negative_integer, "Append this many residues to the input pose", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "final_repack", xsct_rosetta_bool, "Repack the pose at the very end", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "final_minimize", xsct_rosetta_bool, "Minimize the pose at the very end", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "final_mc", xsct_rosetta_bool, "Run a short MC simulation on the pose at the very end", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "final_correct_oop_post", xsct_rosetta_bool, "Do special corrections to the oop_post positions", "false" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string OopCreatorMoverCreator::keyname() const {
	return OopCreatorMover::mover_name();
}

protocols::moves::MoverOP
OopCreatorMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new OopCreatorMover );
}

void OopCreatorMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	OopCreatorMover::provide_xml_schema( xsd );
}




}
}
}

