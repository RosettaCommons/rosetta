// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LoopCreationMover.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <devel/loop_creation/LoopCreationMover.hh>
#include <devel/loop_creation/LoopCreationMoverCreator.hh>

//Basic
#include <basic/MetricValue.hh>
#include <basic/datacache/DiagnosticData.hh>

//Core
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rms_util.hh>

//Protocols
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/analysis/LoopAnalyzerMover.hh>
#include <protocols/loops/loops_main.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>


#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <protocols/filters/BasicFilters.hh>

#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/rosetta_scripts/util.hh>

//Devel
#include <devel/loop_creation/LoopInserter.hh>
#include <devel/loop_creation/LoopCloser.hh>
#include <devel/loop_creation/KICLoopCloser.hh>

//Utility
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/sort_predicates.hh>

//Basic
#include <basic/Tracer.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


namespace devel {
namespace loop_creation {

static THREAD_LOCAL basic::Tracer TR( "devel.loop_creation.LoopCreationMover" );

using namespace core;

//****CREATOR METHODS****//
std::string
LoopCreationMoverCreator::keyname() const
{
	return LoopCreationMoverCreator::mover_name();
}

protocols::moves::MoverOP
LoopCreationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopCreationMover );
}

std::string
LoopCreationMoverCreator::mover_name()
{
	return "LoopCreationMover";
}
//****END CREATOR METHODS****//

/// @brief default constructor
LoopCreationMover::LoopCreationMover():
	loop_inserter_(/* NULL */),
	loop_closer_(/* NULL */),
	loop_filter_(/* NULL */),
	attempts_per_anchor_(1),
	refine_(true),
	design_loops_(false),
	include_neighbors_(false),
	minimize_loops_(true),
	filter_by_lam_(false),
	lam_score_cutoff_(0),
	asym_size_(0),
	dump_pdbs_(false)
{
	init();
}

/// @brief explicit constructor
LoopCreationMover::LoopCreationMover(
	LoopInserterOP loop_inserter,
	LoopCloserOP loop_closer,
	core::Size attempts_per_anchor,
	bool refine,
	bool design_loops,
	bool include_neighbors,
	bool minimize_loops
):
	Mover("LoopCreationMover"),
	loop_inserter_(loop_inserter),
	loop_closer_(loop_closer),
	attempts_per_anchor_(attempts_per_anchor),
	refine_(refine),
	design_loops_(design_loops),
	include_neighbors_(include_neighbors),
	minimize_loops_(minimize_loops),
	filter_by_lam_(false),
	dump_pdbs_(false)
{
	init();
}

void
LoopCreationMover::init()
{
	loop_filter_ = protocols::filters::FilterOP( new protocols::filters::TrueFilter );
}

protocols::moves::MoverOP
LoopCreationMover::clone() const {
	return( protocols::moves::MoverOP( new LoopCreationMover( *this ) ) );
}
protocols::moves::MoverOP
LoopCreationMover::fresh_instance() const {
	return protocols::moves::MoverOP( new LoopCreationMover );
}

std::string
LoopCreationMover::get_name() const {
	return "LoopCreationMover";
}

protocols::loops::Loop
LoopCreationMover::get_last_created_loop() const
{
	return last_created_loop_;
}

//Attempt to build residues to connect residues anchor point to anchor point+1.
//Do this num_iterations time and score each one that passes the loop filter.
//Set the best scoring to the current pose.
void
LoopCreationMover::apply(
	pose::Pose & pose
){
	using namespace std;
	using namespace core;
	using utility::vector1;

	//Sanity checks
	if ( loop_closer_ == 0 ) {
		utility_exit_with_message("No LoopCloser given to LoopCreationMover!");
	}
	if ( loop_inserter_ == 0 ) {
		utility_exit_with_message("No LoopInserter given to LoopCreationMover!");
	}
	if ( loop_filter_ == 0 ) {
		loop_filter_ = protocols::filters::FilterOP( new protocols::filters::TrueFilter );//default to true filter
	}
	TR.Debug << "Fold tree prior to LoopCreationMover: " << pose.fold_tree() << endl;

	protocols::jd2::JobOP const job_me ( protocols::jd2::JobDistributor::get_instance()->current_job() );
	std::string const job_name ( protocols::jd2::JobDistributor::get_instance()->job_outputter()->output_name(job_me) );

	//If not loop anchors were given, get the single anchor from the loop inserter
	if ( loop_anchors_.size()==0 ) {
		//loop_inserter_->loop_anchor(loop_anchor_);
		loop_anchors_.push_back(loop_inserter_->loop_anchor());
	}

	protocols::loops::Loop created_loop;

	pose::Pose edit_pose;
	for ( core::Size i=1; i<=loop_anchors_.size(); ++i ) {
		bool loop_closed = false;
		bool loop_passed = false;
		pose::remove_upper_terminus_type_from_pose_residue(pose, loop_anchors_[i]);
		pose::remove_lower_terminus_type_from_pose_residue(pose, loop_anchors_[i]+1);

		loop_inserter_->loop_anchor(loop_anchors_[i]);
		for ( Size j=1; j<=attempts_per_anchor_; ++j ) {
			//already succeeded for this anchor
			if ( loop_passed ) break;

			//copy pose to get a working copy
			edit_pose = pose;

			//Insert new loop residues
			TR << "Beginning loop creation attempt " << j << " of " << attempts_per_anchor_
				<< ", for loop anchor " << i << " of " << loop_anchors_.size() << std::endl;

			clock_t start_time = clock();
			loop_inserter_->apply(edit_pose);
			clock_t insert_time = clock() - start_time;
			TR << "Clocks - loop inserter run " << j << " of " << attempts_per_anchor_
				<< " - total time: " << insert_time << std::endl;

			created_loop = loop_inserter_->get_created_loop();
			TR.Debug << "New loop: " << created_loop << endl;

			//DEBUG
			if ( dump_pdbs_ ) {
				std::stringstream filename;
				filename << job_name << "_insertion_" << j << ".pdb";
				edit_pose.dump_pdb(filename.str());
			}
			//END DEBUG

			start_time = clock();
			loop_closer_->loop(created_loop);
			loop_closer_->apply(edit_pose);
			clock_t closure_time = clock() - start_time;
			TR << "Clocks - loop closer run " << j << " of " << attempts_per_anchor_
				<< " - total time: " << closure_time << std::endl;

			if ( loop_closer_->success() ) {
				loop_closed=true;
				TR.Debug << "Loop insertion/closure " << j << " succeeded" << std::endl;

				//ensure that a chemical bond exists now that we've closed the loop
				edit_pose.conformation().declare_chemical_bond(created_loop.cut(), "C", created_loop.cut()+1,"N");

				//For some reason the creation/closing of loops really f***s the hydrogen placement.
				//This is an attempt to fix that
				for ( core::Size i=created_loop.start()-1; i<=created_loop.stop()+1; ++i ) {
					core::conformation::ResidueOP new_res = edit_pose.residue(i).clone();
					core::conformation::idealize_hydrogens(*new_res, edit_pose.conformation());
					edit_pose.conformation().replace_residue(i, *new_res, false);
				}

				//DEBUG
				if ( dump_pdbs_ ) {
					std::stringstream filename;
					filename << job_name << "_closure_" << j << ".pdb";
					edit_pose.dump_pdb(filename.str());
				}
				//END DEBUG

				if ( refine_ ) {
					start_time = clock();
					refine_loop(edit_pose, created_loop);
					clock_t refine_time = clock()-start_time;
					TR.Debug << "Clocks - loop refinement completed in total time: " << refine_time << std::endl;

					//DEBUG
					if ( dump_pdbs_ ) {
						std::stringstream filename;
						filename << job_name << "_creation_" << j << "_refined.pdb";
						edit_pose.dump_pdb(filename.str());
					}
					//END DEBUG
				}

				if ( loop_filter_->apply(edit_pose) ) {

					TR << "Running LoopAnalyzerMover for job: " << job_name
						<< ", loop: " << created_loop << std::endl;

					protocols::loops::Loops loops;
					//loops.add_loop(last_created_loop_);
					loops.add_loop(protocols::loops::Loop(loop_inserter_->modified_range().first, loop_inserter_->modified_range().second));
					protocols::analysis::LoopAnalyzerMoverOP lam( new protocols::analysis::LoopAnalyzerMover(loops, true) );
					lam->apply(edit_pose);

					core::Real total_loop_score = lam->get_total_score();
					core::Real max_rama = lam->get_max_rama();
					core::Real max_omega = lam->get_max_omega();
					core::Real max_pbond = lam->get_max_pbond();
					core::Real max_chainbreak = lam->get_max_chainbreak();

					TR << "Total loop analyzer score: " << total_loop_score << std::endl;
					TR << "LAM score cutoff: " << lam_score_cutoff_ << std::endl;
					TR << "Max rama: " << max_rama << std::endl;
					TR << "Max omega: " << max_omega << std::endl;
					TR << "Max peptide bond: " << max_pbond << std::endl;
					TR << "Max chainbreak: " << max_chainbreak << std::endl;

					//Use the LoopAnalyzerMover to filter loops
					if ( filter_by_lam_ ) {
						if ( total_loop_score <= lam_score_cutoff_ ) {
							loop_passed=true;
							pose=edit_pose;
							last_created_loop_=created_loop;
							job_me->add_string_real_pair("total_loop_score", total_loop_score);
							job_me->add_string_real_pair("max_loop_rama", max_rama);
							job_me->add_string_real_pair("max_loop_omega", max_omega);
							job_me->add_string_real_pair("max_loop_pbond", max_pbond);
							job_me->add_string_real_pair("max_loop_chainbreak", max_chainbreak);
						}
					} else {
						//done with this anchor!
						loop_passed=true;
						pose=edit_pose;
						last_created_loop_=created_loop;
					}
				}
			}
		}
		if ( !loop_closed ) {
			TR << "Unable to create a closed loop between residues "
				<< loop_inserter_->loop_anchor() << " and " << loop_inserter_->loop_anchor()+1
				<< ", consider longer loop lengths or a different LoopInserter and/or LoopCloser." << endl;
			set_last_move_status(protocols::moves::FAIL_RETRY);
			return;
			//utility_exit_with_message(err.str());
		}
		if ( !loop_passed ) {
			TR << "No closed loops between residues "
				<< loop_inserter_->loop_anchor() << " and " << loop_inserter_->loop_anchor()+1
				<< " passed filters. Try relaxing filters." << endl;
			//utility_exit_with_message(err.str());
			set_last_move_status(protocols::moves::FAIL_RETRY);
			return;
		}

		update_anchors(loop_anchors_, last_created_loop_, i);

		core::Size temp_asym_size = asym_size_;
		if ( asym_size_ != 0 ) {
			core::Size cur_anchor=loop_anchors_[i];
			temp_asym_size+=last_created_loop_.size();
			core::Size dup_anchor=cur_anchor+temp_asym_size;

			TR << "first dup anchor is: " << dup_anchor << std::endl;
			TR << "asym size is: " << temp_asym_size << std::endl;
			TR << "loop residue size is: " << last_created_loop_.size() << std::endl;

			core::Size last_protein_resnum(0);
			for ( core::Size i(1); i<=pose.total_residue(); ++i ) {
				if ( pose.residue(i).is_protein() ) {
					last_protein_resnum = i;
				}
			}

			utility::vector1<core::Size> chain_endings;
			chain_endings.push_back( cur_anchor + last_created_loop_.size() );
			while ( dup_anchor < last_protein_resnum ) //strictly less than so we don't try to build a loop that connects to nothing
					{
				TR << "Duplication loop build on anchor " << cur_anchor << " to new anchor: " << dup_anchor << std::endl;
				copy_last_loop_to_new_anchor(pose, dup_anchor);
				chain_endings.push_back( dup_anchor + last_created_loop_.size() );
				dup_anchor+=temp_asym_size;

			}

			pose.conformation().chain_endings( chain_endings );


		}
	}//loop anchors
}

/// @brief update the loops to reflect
///the position changes after an insertion
void
LoopCreationMover::update_anchors(
	utility::vector1<core::Size> & loop_anchors,
	protocols::loops::Loop const & new_loop,
	core::Size index_of_new_loop
){
	for ( core::Size i=index_of_new_loop+1; i<=loop_anchors.size(); ++i ) {
		loop_anchors_[i]+=new_loop.size();
	}
}

void
LoopCreationMover::copy_last_loop_to_new_anchor(
	core::pose::Pose & pose,
	core::Size new_anchor
){
	using namespace core;
	using namespace core::chemical;

	TR.Debug << "Loop to copy: " << last_created_loop_ << std::endl;
	pose::remove_upper_terminus_type_from_pose_residue(pose, new_anchor);
	pose::remove_lower_terminus_type_from_pose_residue(pose, new_anchor+1);

	core::Size old_anchor = loop_inserter_->loop_anchor();

	core::Size modifications_begin = loop_inserter_->modified_range().first;
	core::Size mod_begin_offset = old_anchor-modifications_begin;
	core::Size new_modifications_begin=new_anchor-mod_begin_offset;

	core::Size modifications_end = loop_inserter_->modified_range().second;
	core::Size mod_end_offset = old_anchor-modifications_end;
	core::Size new_modifications_end=new_anchor-mod_end_offset;
	new_modifications_end-=last_created_loop_.size();

	TR.Debug << "Old modifications range: " << modifications_begin << " " << modifications_end << std::endl;
	TR.Debug << "New modification range " << new_modifications_begin << " " << new_modifications_end << std::endl;

	core::Size modified_size = modifications_end-modifications_begin+1;
	core::Size loop_size = last_created_loop_.size();
	core::Size flanking_size = (modified_size-loop_size)/2;
	TR.Debug << "Size of modified region " << modified_size << std::endl;
	TR.Debug << "Size of loop region " << loop_size << std::endl;
	TR.Debug << "Size of flanking region " << flanking_size << std::endl;
	runtime_assert((flanking_size*2)+loop_size == modified_size);

	core::pose::Pose loop_pose(pose, modifications_begin-1, modifications_end);
	for ( core::Size i=1; i<=loop_pose.total_residue(); ++i ) {
		pose::remove_upper_terminus_type_from_pose_residue(loop_pose, i);
		pose::remove_lower_terminus_type_from_pose_residue(loop_pose, i);
	}

	core::id::AtomID_Map< core::id::AtomID > atom_map;
	atom_map.clear();
	core::pose::initialize_atomid_map( atom_map, loop_pose, core::id::BOGUS_ATOM_ID );

	core::id::AtomID const id1( pose.residue(new_modifications_begin-1).atom_index("CA"), new_modifications_begin-1);
	core::id::AtomID const id2( loop_pose.residue(1).atom_index("CA"), 1);
	atom_map[ id2 ] = id1;

	core::id::AtomID const id3( pose.residue(new_modifications_begin-1).atom_index("C"), new_modifications_begin-1);
	core::id::AtomID const id4( loop_pose.residue(1).atom_index("C"), 1);
	atom_map[ id4 ] = id3;

	core::id::AtomID const id5( pose.residue(new_modifications_begin-1).atom_index("N"), new_modifications_begin-1);
	core::id::AtomID const id6( loop_pose.residue(1).atom_index("N"), 1);
	atom_map[ id6 ] = id5;

	core::id::AtomID const id7( pose.residue(new_modifications_begin-1).atom_index("O"), new_modifications_begin-1);
	core::id::AtomID const id8( loop_pose.residue(1).atom_index("O"), 1);
	atom_map[ id8 ] = id7;

	core::scoring::superimpose_pose(loop_pose, pose/*const*/, atom_map);
	//loop_pose.dump_pdb("loop_pose_aligned.pdb");

	for ( core::Size i=0; i<modified_size; ++i ) {
		if ( i >= flanking_size && i < flanking_size+loop_size ) {
			TR.Debug << "appending loop residue " << i+2 << " to residue " << new_modifications_begin+i-1 << std::endl;
			pose.append_polymer_residue_after_seqpos(loop_pose.residue(i+2), new_modifications_begin+i-1, false);
		} else {
			TR.Debug << "replacing residue " << new_modifications_begin+i << " with loop_pose residue " << i+2 << std::endl;
			pose.replace_residue(new_modifications_begin+i, loop_pose.residue(i+2), false);
		}
	}
	//pose.dump_pdb("cloned_loop.pdb");
}

core::Real
LoopCreationMover::refine_loop(
	core::pose::Pose & pose,
	protocols::loops::Loop loop
){
	using namespace std;
	using namespace core;
	using utility::vector1;


	TR << "REFINE LOOP" << std::endl;

	core::scoring::ScoreFunctionOP scorefxn_min;
	if ( scorefxn_ == 0 ) {
		scorefxn_ = scoring::get_score_function();
		scorefxn_min = scorefxn_->clone();
	} else {
		scorefxn_min = scorefxn_->clone();
		TR << "found input scorefxn" << std::endl;
	}

	scorefxn_min->set_weight( scoring::chainbreak, 100.0 );//loop should already be closed, so this will just prevent re-breaking by minimization

	//Initialize a mover if we haven't already
	protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover );
	pack_mover->score_function( scorefxn_ );

	//Setup task factory for minimization and packing
	pack::task::TaskFactoryOP task_factory( new pack::task::TaskFactory );

	set<Size> loop_residues;
	for ( Size i=loop.start(); i<=loop.stop(); ++i ) {
		loop_residues.insert(i);
	}

	//If we're including loop neighbors in packing/redesign then calculate them
	set<Size> residues_to_pack;
	if ( include_neighbors_ ) {
		string const nb_calc("neighbor_calculator");
		pose::metrics::CalculatorFactory::Instance().register_calculator( nb_calc,
			core::pose::metrics::PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( loop_residues ) ) );

		basic::MetricValue< set< Size > > neighbor_mv;
		pose.metric( nb_calc, "neighbors", neighbor_mv);
		residues_to_pack = neighbor_mv.value();
		pose::metrics::CalculatorFactory::Instance().remove_calculator(nb_calc);
	} else {
		residues_to_pack=loop_residues;
	}

	//Restrict all non-packable residues and allow movements at all packable residues
	pack::task::operation::PreventRepackingOP prevent_packing( new pack::task::operation::PreventRepacking );
	pack::task::operation::RestrictResidueToRepackingOP repack_res( new pack::task::operation::RestrictResidueToRepacking() );
	for ( Size i=1; i<=pose.total_residue(); ++i ) {
		if ( residues_to_pack.find(i)==residues_to_pack.end() ) {
			prevent_packing->include_residue(i);
		} else {
			//If we're not designing then restrict packable residues to repack-only
			if ( !design_loops_ ) {
				repack_res->include_residue(i);
			}
		}
	}
	task_factory->push_back(prevent_packing);
	task_factory->push_back(repack_res);

	pack_mover->task_factory(task_factory);
	pack_mover->apply(pose);

	if ( minimize_loops_ ) {
		core::Size modifications_begin = loop_inserter_->modified_range().first;
		core::Size modifications_end = loop_inserter_->modified_range().second;
		protocols::loops::Loop extended_loop(modifications_begin, modifications_end, loop.cut());
		protocols::loops::set_single_loop_fold_tree(pose, extended_loop);

		kinematics::MoveMapOP movemap( new kinematics::MoveMap );
		for ( core::Size i=extended_loop.start(); i<=extended_loop.stop(); ++i ) {
			movemap->set_bb(i, true);
			movemap->set_chi(i, true);
		}
		TR.Debug << "Movemap for minimization: " << *movemap << std::endl;

		protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover(movemap, scorefxn_min, "lbfgs_armijo_nonmonotone", 0.01, false ) );

		TR << "Score prior to minimization: " << scorefxn_min->score(pose) << std::endl;
		//  pose.dump_pdb("pre_minimization.pdb");

		min_mover->apply(pose);
	}

	// utility::vector1< bool > residues_to_score(pose.total_residue(), false);
	// for(core::Size i=loop.start(); i<=loop.stop(); ++i)
	// {
	//  residues_to_score[i]=true;
	// }
	//
	// core::Real loop_score = scorefxn->get_sub_score(pose, residues_to_score);
	core::Real loop_score = scorefxn_min->score(pose);
	return (loop_score/loop.size());
}

LoopInserterOP
LoopCreationMover::loop_inserter() const
{
	return loop_inserter_;
}

LoopCloserOP
LoopCreationMover::loop_closer() const
{
	return loop_closer_;
}

void
LoopCreationMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & movers,
	Pose const & /*pose*/
){
	using namespace std;

	TR << "Parsing tag for LoopCreationMover" << std::endl;
	if ( tag->hasOption("scorefxn") ) {
		std::string const scorefxn_name = tag->getOption<std::string>("scorefxn");
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data, scorefxn_name );
		TR << "Found input scorefxn: " << scorefxn_name << std::endl;
	}


	//****REQUIRED TAGS****//
	if ( tag->hasOption("loop_inserter") ) {
		string const loop_inserter_name( tag->getOption< string >( "loop_inserter" ) );
		protocols::moves::Movers_map::const_iterator find_mover( movers.find( loop_inserter_name ) );
		bool const mover_found( find_mover != movers.end() );
		if ( mover_found ) {
			loop_inserter_ = utility::pointer::dynamic_pointer_cast< devel::loop_creation::LoopInserter > ( find_mover->second );
		} else {
			utility_exit_with_message( "Mover " + loop_inserter_name + " not found" );
		}
	} else {
		utility_exit_with_message("loop_inserter tag is required");
	}

	if ( tag->hasOption("loop_closer") ) {
		string const loop_closer_name( tag->getOption< string >( "loop_closer" ) );
		protocols::moves::Movers_map::const_iterator find_mover( movers.find( loop_closer_name ) );
		bool const mover_found( find_mover != movers.end() );
		if ( mover_found ) {
			loop_closer_ = utility::pointer::dynamic_pointer_cast< devel::loop_creation::LoopCloser > ( find_mover->second );
		} else {
			utility_exit_with_message( "Mover " + loop_closer_name + " not found" );
		}
	} else {
		utility_exit_with_message("loop_closer tag is required");
	}

	//****OPTIONAL TAGS****//

	if ( tag->hasOption("dump_pdbs") ) {
		dump_pdbs_ = tag->getOption<bool>("dump_pdbs");
	}

	if ( tag->hasOption("refine") ) {
		refine_ = tag->getOption<bool>("refine");
	}

	if ( tag->hasOption("attempts_per_anchor") ) {
		attempts_per_anchor_ = tag->getOption< core::Size >("attempts_per_anchor");
	}

	if ( tag->hasOption("minimize_loops") ) {
		minimize_loops_ = tag->getOption<bool>("minimize_loops");
	}

	if ( tag->hasOption("design_loops") ) {
		design_loops_ = tag->getOption<bool>("design_loops");
	}

	if ( tag->hasOption("include_neighbors") ) {
		include_neighbors_ = tag->getOption<bool>("include_neighbors");
	}

	if ( tag->hasOption("loop_anchors") ) {
		string const loop_anchors_string = tag->getOption<string>("loop_anchors");
		utility::vector1<string> loop_anchor_strings=utility::string_split(loop_anchors_string, ',');
		for ( core::Size i=1; i<=loop_anchor_strings.size(); ++i ) {
			loop_anchors_.push_back(utility::string2int(loop_anchor_strings[i]));
		}
	}

	if ( tag->hasOption("asym_size") ) {
		asym_size_ = tag->getOption<core::Size>("asym_size");
	}

	if ( tag->hasOption("filter_by_lam") ) {
		filter_by_lam_ = tag->getOption<bool>("filter_by_lam");
	}
	if ( tag->hasOption("lam_score_cutoff") ) {
		lam_score_cutoff_ = tag->getOption<core::Real>("lam_score_cutoff");
	}
}

} //loop creation
} //devel
