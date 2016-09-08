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
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/RotamerSampleOptions.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/MMAtomType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/ncbb/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/AtomCoordinateCstMover.hh>
#include <protocols/relax/cst_util.hh>

#include <protocols/residue_optimization/MetapatchEnumeration.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <utility/vector1.functions.hh>

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
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <sstream>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::rigid;
using namespace protocols::toolbox;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static THREAD_LOCAL basic::Tracer TR("protocols.residue_optimization.MetapatchEnumeration" );

namespace protocols {
namespace residue_optimization {

/// OPPORTUNITIES FOR IMPROVEMENT
// You should be able to pass this guy a mover or entire protocol with which to do your pre/post sampling
// MPI?

bool
MetapatchEnumeration::tabooed(
	utility::vector1< std::string > const & patch_names
) {
	// Loop through tabooed_ and see if patch_names contains all of any element thereof.
	for ( Size ii = 1; ii <= tabooed_.size(); ++ii ) {
		utility::vector1< std::string > taboo_set = tabooed_[ ii ];

		if ( patch_names.size() < taboo_set.size() ) continue;
		bool matched_whole_set = true;
		for ( Size jj = 1; jj <= taboo_set.size(); ++jj ) {
			// Is the jj-th element of the tabooed set present?
			bool found = false;
			for ( Size kk = 1; kk <= patch_names.size(); ++kk ) {
				if ( taboo_set[ jj ] == patch_names[ kk ] ) {
					found = true;
				}
			}
			if ( !found ) {
				matched_whole_set = false;
				break;
			}
		}
		if ( matched_whole_set ) {
			TR << "Matched taboo set. Skipping." << std::endl;
			return true;
		}
	}
	return false;
}

void
MetapatchEnumeration::generate_derived_types(
	Size resi,
	Size tp,
	utility::vector1< std::string > & types_considered
) {
	ResidueTypeSetCAP rts_cap = chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );
	ResidueTypeSetCOP rts     = ResidueTypeSetCOP( rts_cap );
	ResidueType const current_type = rts->name_map( types_considered[ tp ] );

	TR << "Considering derivatives of " << current_type.name() << "." << std::endl;
	// If patch names contain a tabooed mutation set, move on.
	if ( tabooed( get_patch_names( current_type ) ) ) return;

	for ( Size mp = 1; mp <= metapatch_names_.size(); ++mp ) {
		//// Temporary storage for types created in the process of metapatch application.
		utility::vector1< std::string > newly_viable_types;

		utility::vector1< std::string > good_atoms = rts->metapatch( metapatch_names_[ mp ] )->atoms( current_type );
		for ( utility::vector1< std::string >::const_iterator at = good_atoms.begin(); at != good_atoms.end(); ++at ) {

			std::string const trimmed_atom = *at;
			std::string const base_name    = residue_type_base_name( current_type );
			std::string const patch_name   = base_name + "-" + trimmed_atom + "-" + metapatch_names_[ mp ];

			utility::vector1< std::string > patch_names = get_patch_names( current_type );

			//// If, by some perversion, you can find current atom name in patch names, move on.
			// If trimmed atom is less than any already present atom, move on
			bool move_on = false;
			for ( Size pn = 1; pn <= patch_names.size(); ++pn ) {
				utility::vector1< std::string > components = utility::string_split( patch_names[pn], '-' );
				if ( components[1] != base_name )     continue;
				if ( components[2] >= trimmed_atom )  move_on = true;
			}
			if ( move_on ) continue;

			patch_names.push_back( patch_name );

			// Sort metapatch derived names if needed.
			//utility::vector1< std::string >::iterator it = patch_names.begin();
			//for ( ; it != patch_names.end(); ++it ) {
			// if ( it->find( base_name ) ) break;
			//}
			//std::sort( it, patch_names.end() );

			// Reassemble patch names plus base name into a name.
			std::string full_name = base_name;
			for ( Size pn = 1; pn <= patch_names.size(); ++pn ) {
				full_name += ":" + patch_names[ pn ];
			}
			ResidueType const & rt = rts->name_map( full_name );

			if ( rt.name() != full_name ) {
				TR << "ERROR: name mismatch " << rt.name() << " " << full_name << std::endl;
				continue;
			}
			newly_viable_types.push_back( full_name );

			Pose mut_pose = pose_;
			ResidueOP new_res = ResidueFactory::create_residue( rt, mut_pose.residue( resi ), mut_pose.conformation(), true );
			if ( new_res == nullptr ) continue;

			TR << "\tEvaluating mutation: " << full_name << "." << std::endl;

			core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( mut_pose.residue( resi ), *new_res, mut_pose.conformation() );
			mut_pose.conformation().replace_residue( resi, *new_res, false );

			final_sampling( mut_pose, resi );

			Real mut_score = ( *evaluation_score_fxn_ )( mut_pose );
			mut_score -= mut_pose.energies().onebody_energies( resi )[ mm_twist ] * evaluation_score_fxn_->get_weight( mm_twist );

			Real delta = ( mode_ == 1 ? mut_score : binding( mut_pose ) ) - score_;

			TR << "\tMutation changed energy by: " << delta << "." << std::endl;
			if ( delta < threshold_ ) {
				std::stringstream fn;
				fn << "pose_" << resi << "_" << rt.name() << ".pdb";
				final_scores_.push_back( delta );
				summary_lines_.push_back( fn.str() );
				mut_pose.dump_scored_pdb( fn.str(), *evaluation_score_fxn_ );
			}

			if ( delta > taboo_ ) {
				TR << "\tTabooing that combination of modifications." << std::endl;
				tabooed_.push_back( patch_names );
			}
		}

		for ( Size nvt = 1; nvt <= newly_viable_types.size(); ++nvt ) {
			types_considered.push_back( newly_viable_types[nvt] );
		}
	}
}

void
MetapatchEnumeration::generate_metapatched_variants( Size resi ) {
	std::string starting_name = pose_.residue( resi ).type().name();
	ResidueType const starting_type = pose_.residue( resi ).type();

	// Reset mm bb except for area near resi
	for ( Size ii = 1; ii <= pose_.size(); ++ii ) {
		if ( ii >= resi - 2 || ii <= resi + 2  ) {
			mm_->set_bb( ii, false );
		} else  {
			mm_->set_bb( ii, true );
		}
	}

	// Instead of doing this, just eliminate mm_twist scores from the variable
	// residue
	//evaluation_score_fxn_->set_weight( mm_twist, 0 );
	// Binding mode can use sampling score function to score
	initial_sampling( pose_ );
	score_ = ( mode_ == 1 ) ? ( *evaluation_score_fxn_ )( pose_ ) : binding( pose_ );

	pose_.dump_scored_pdb( "initial_sampled.pdb", *evaluation_score_fxn_ );

	utility::vector1< std::string > types_considered;
	types_considered.push_back( starting_type.name() );

	for ( Size tp = 1; tp <= types_considered.size(); ++tp ) {
		TR << "Generating derived types at " << resi << ", number " << tp << " of " << types_considered.size() << "." << std::endl;
		generate_derived_types( resi, tp, types_considered );
	}
}

Real
MetapatchEnumeration::binding( Pose pose ) {
	Real v1 = ( *evaluation_score_fxn_ )( pose );

	protocols::rigid::RigidBodyTransMover trans_mover( pose, 1 );
	trans_mover.step_size( 1000 );
	trans_mover.apply( pose );
	Real v2 = ( *evaluation_score_fxn_ )( pose );

	return ( v1 - v2 );
}

void
MetapatchEnumeration::initial_sampling(
	Pose & pose
) {
	simple_moves::MinMoverOP min( new simple_moves::MinMover( mm_, sampling_score_fxn_, "lbfgs_armijo_nonmonotone", 0.01, true )  );

	if ( pack_ ) {
		PackerTaskOP setup_pack = TaskFactory::create_packer_task( pose );
		setup_pack->initialize_from_command_line().or_include_current( true );
		for ( core::Size resj = 1; resj <= pose.size(); ++resj ) {
			setup_pack->nonconst_residue_task( resj ).or_ex1_sample_level( pack::task::EX_SIX_QUARTER_STEP_STDDEVS );
			setup_pack->nonconst_residue_task( resj ).or_ex2_sample_level( pack::task::EX_SIX_QUARTER_STEP_STDDEVS );
			setup_pack->nonconst_residue_task( resj ).restrict_to_repacking();
		}

		protocols::simple_moves::PackRotamersMoverOP pack( new protocols::simple_moves::PackRotamersMover( sampling_score_fxn_, setup_pack ) );

		Real score = 100000;
		for ( Size pc = 1; pc <= 10; ++pc ) {
			Pose packpose = pose;
			pack->apply( packpose );
			min->apply( packpose );
			Real newscore = ( *evaluation_score_fxn_ )( packpose );
			if ( newscore < score ) {
				score = newscore;
				pose = packpose;
			}
		}
	} else {

		// Add virtual root if one doesn't exist
		if ( pose.residue( pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			core::pose::addVirtualResAsRoot( pose );
		}

		protocols::relax::AtomCoordinateCstMover coord_cst;
		//pose::PoseOP ref_pose( new Pose( pose ) );
		//coord_cst.set_refstruct( ref_pose );
		coord_cst.apply( pose );


		protocols::relax::FastRelax fast_relax( sampling_score_fxn_, 3 );
		fast_relax.set_movemap( mm_ );
		fast_relax.apply( pose );

		protocols::relax::delete_virtual_residues( pose );
	}

	min->apply( pose );
}

void MetapatchEnumeration::final_sampling(
	Pose & mut_pose,
	Size resi
) {
	kinematics::MoveMapOP mut_mm = mm_;
	//for ( Size ii = 1; ii <= mut_pose.size(); ++ii ) {
	// mut_mm->set_bb(  ii, true );
	// mut_mm->set_chi( ii, true );
	//}

	simple_moves::MinMoverOP mut_min( new simple_moves::MinMover( mut_mm, sampling_score_fxn_, "lbfgs_armijo_nonmonotone", 0.01, true )  );

	if ( pack_ ) {
		TaskFactoryOP mutant_tf( new TaskFactory );
		protocols::toolbox::task_operations::DesignAroundOperationOP dao( new protocols::toolbox::task_operations::DesignAroundOperation );///restrict repacking to 8.0A around target res to save time
		dao->include_residue( resi );
		mutant_tf->push_back( dao );
		PackerTaskOP mutant_pt = mutant_tf->create_task_and_apply_taskoperations( mut_pose );
		mutant_pt->initialize_from_command_line().or_include_current( true );
		for ( core::Size resj = 1; resj <= mut_pose.size(); ++resj ) {
			mutant_pt->nonconst_residue_task( resj ).or_ex1_sample_level( pack::task::EX_SIX_QUARTER_STEP_STDDEVS );
			mutant_pt->nonconst_residue_task( resj ).or_ex2_sample_level( pack::task::EX_SIX_QUARTER_STEP_STDDEVS );
			mutant_pt->nonconst_residue_task( resj ).restrict_to_repacking();
		}
		protocols::simple_moves::PackRotamersMoverOP pack_mut( new protocols::simple_moves::PackRotamersMover( sampling_score_fxn_, mutant_pt ) );

		Real score = 100000;
		for ( Size pc = 1; pc <= 10; ++pc ) {
			Pose packpose = mut_pose;
			pack_mut->apply( packpose );
			mut_min->apply( packpose );
			Real newscore = ( *evaluation_score_fxn_ )( packpose );
			if ( newscore < score ) {
				score = newscore;
				mut_pose = packpose;
			}
		}
	} else {
		// Add virtual root if one doesn't exist
		if ( mut_pose.residue( mut_pose.fold_tree().root() ).aa() != core::chemical::aa_vrt ) {
			core::pose::addVirtualResAsRoot( mut_pose );
		}

		protocols::relax::AtomCoordinateCstMover coord_cst;
		//pose::PoseOP ref_pose( new Pose( mut_pose ) );
		//coord_cst.set_refstruct( ref_pose );
		coord_cst.apply( mut_pose );

		protocols::relax::FastRelax fast_relax( sampling_score_fxn_, 5 );
		fast_relax.set_movemap( mut_mm );
		fast_relax.apply( mut_pose );
		protocols::relax::delete_virtual_residues( mut_pose );

		/*Real score = 100000;
		for ( Size pc = 1; pc <= 3; ++pc ) {
		Pose relaxpose = mut_pose;
		fast_relax.apply( mut_pose );
		Real newscore = ( *evaluation_score_fxn_ )( relaxpose );
		if ( newscore < score ) {
		score = newscore;
		mut_pose = relaxpose;
		}
		}*/

		mut_pose.remove_constraints();
	}

	mut_min->apply( mut_pose );
}

}
}

