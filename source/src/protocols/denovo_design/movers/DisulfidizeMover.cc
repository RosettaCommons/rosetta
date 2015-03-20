// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/movers/DisulfidizeMover.cc
/// @brief The DisulfidizeMover
/// @detailed
/// @author Tom Linsky (tlinsky@uw.edu) -- Adapting code from remodelmover into a mover
/// @author Gabe Rocklin (grocklin@uw.edu) -- Disulfide code


//Unit Headers
#include <protocols/denovo_design/movers/DisulfidizeMover.hh>
#include <protocols/denovo_design/movers/DisulfidizeMoverCreator.hh>

//Project Headers
#include <protocols/denovo_design/util.hh>

//Protocol Headers
#include <protocols/forge/remodel/RemodelDesignMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>

//Core Headers
#include <core/conformation/util.hh>
#include <core/pack/task/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/util/disulfide_util.hh>

//Basic Headers
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

//C++ Headers

static thread_local basic::Tracer TR( "protocols.denovo_design.DisulfidizeMover" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
////////////////////////////////////////////////////////////////////////////////////////////////////

std::string
DisulfidizeMoverCreator::keyname() const
{
  return DisulfidizeMoverCreator::mover_name();
}

protocols::moves::MoverOP
DisulfidizeMoverCreator::create_mover() const
{
  return protocols::moves::MoverOP( new DisulfidizeMover() );
}

std::string
DisulfidizeMoverCreator::mover_name()
{
  return "DisulfidizeMover";
}

///  ---------------------------------------------------------------------------------
///  DisulfidizeMover main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
DisulfidizeMover::DisulfidizeMover() :
	protocols::moves::Mover(),
	match_rt_limit_( 1.0 ),
	max_disulf_score_( -0.25 ),
	min_loop_( 8 ),
	min_disulfides_( 1 ),
	max_disulfides_( 3 ),
	include_current_ds_( false ),
	keep_current_ds_( false ),
	last_pose_()
{
	accumulator_.clear();
}

/// @brief destructor - this class has no dynamic allocation, so
DisulfidizeMover::~DisulfidizeMover()
{
}

std::string
DisulfidizeMover::get_name() const
{
	return DisulfidizeMoverCreator::mover_name();
}

/// Return a copy of ourselves
protocols::moves::MoverOP
DisulfidizeMover::clone() const
{
	return protocols::moves::MoverOP( new DisulfidizeMover(*this) );
}

void
DisulfidizeMover::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const & ,
		protocols::moves::Movers_map const & ,
		core::pose::Pose const & )
{
	if ( tag->hasOption( "match_rt_limit" ) )
		match_rt_limit_ = tag->getOption< core::Real >( "match_rt_limit" );
	if ( tag->hasOption( "min_disulfides" ) )
		min_disulfides_ = tag->getOption< core::Size >( "min_disulfides" );
	if ( tag->hasOption( "max_disulfides" ) )
		max_disulfides_ = tag->getOption< core::Size >( "max_disulfides" );
	if ( tag->hasOption( "keep_current_disulfides" ) )
		keep_current_ds_ = tag->getOption< bool >( "keep_current_disulfides" );
	if ( tag->hasOption( "include_current_disulfides" ) )
		include_current_ds_ = tag->getOption< bool >( "include_current_disulfides" );
	if ( tag->hasOption( "min_loop" ) )
		min_loop_ = tag->getOption< core::Size >( "min_loop" );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
DisulfidizeMover::apply( core::pose::Pose & pose )
{
	bool const newpose = ( !last_pose_ || !same_pose( pose, *last_pose_ ) );

	if ( newpose ) {
		if ( accumulator_.size() ) {
			// warn the user if the accumulator still has poses and we are starting with a new input
			TR << "Warning: the accumulator still has poses, but the DisulfidizeMover has been passed a new input. Poses in the accumulator will be lost." << std::endl;
			accumulator_.clear();
		}
		generate_results( pose );
		last_pose_ = core::pose::PoseCOP( new core::pose::Pose(pose) );
	}

	// if this is the same pose we saw last time, pop a result off the accumulator
	if ( accumulator_.size() ) {
		pose = *(pop_result());
		TR << "Returning a result from the accumulator -- there are " << accumulator_.size() << " poses left." << std::endl;
		set_last_move_status( protocols::moves::MS_SUCCESS );
	} else {
		TR << "There are no poses left in the accumulator for this input pose. Gracefully failing." << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
	}
}

/// @brief populate the internally cached list of results for a given pose
void
DisulfidizeMover::generate_results( core::pose::Pose const & pose )
{
	assert( accumulator_.size() == 0 );

	// create initial list of possible disulfides between residue subset 1 and subset 2
	core::pack::task::residue_selector::ResidueSubset subset1( pose.total_residue(), true );
	core::pack::task::residue_selector::ResidueSubset subset2( pose.total_residue(), true );
	DisulfideList disulf_partners = find_disulfides_in_the_neighborhood( pose, subset1, subset2 );
	if ( disulf_partners.size() == 0 ) {
		if ( min_disulfides_ == 0 ) {
			push_result( pose.clone() );
			return;
		}
		TR << "Failed to build any disulfides." << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
		return;
	}

	//Use the recursive multiple disulfide former
	DisulfideList empty_disulfide_list;
	utility::vector1< DisulfideList > disulfide_configurations =
		recursive_multiple_disulfide_former( empty_disulfide_list, disulf_partners );

	// iterate over disulfide configurations
	for ( utility::vector1< DisulfideList >::const_iterator ds_config = disulfide_configurations.begin();
			ds_config != disulfide_configurations.end();
			++ds_config) {
		if ( (*ds_config).size() >= min_disulfides_ && (*ds_config).size() <= max_disulfides_ ) {

			//form all the disulfides in the disulfide configuration
			TR << "Building disulfide configuration ";
			for ( DisulfideList::const_iterator my_ds = (*ds_config).begin(); my_ds != (*ds_config).end(); ++my_ds) {
				TR << (*my_ds).first << "-" << (*my_ds).second << " ";
			}
			TR << std::endl;

			core::pose::PoseOP disulf_copy_pose( new core::pose::Pose(pose) );
			make_disulfides( *disulf_copy_pose, *ds_config, false );
			push_result( disulf_copy_pose );
		}
	}
}

/// @brief pushes a result onto the accumulator
void
DisulfidizeMover::push_result( core::pose::PoseCOP pose )
{
	accumulator_.push_back( pose );
}

/// @brief pops a result off of the front of the accumulator (FIFO), decreasing its size by one
core::pose::PoseCOP
DisulfidizeMover::pop_result()
{
	assert( accumulator_.size() );
	core::pose::PoseCOP retval = accumulator_.front();
	accumulator_.pop_front();
	assert( retval );
	return retval;
}

/// @brief Function for recursively creating multiple disulfides
utility::vector1< DisulfidizeMover::DisulfideList >
DisulfidizeMover::recursive_multiple_disulfide_former(
		DisulfideList const & disulfides_formed,
		DisulfideList const & disulfides_possible ) const
{
	utility::vector1< DisulfideList > final_configurations;

	if ( disulfides_formed.size() < max_disulfides_ ) {

		//select one primary new disulfide to be added
		for ( DisulfideList::const_iterator new_disulfide = disulfides_possible.begin();
				new_disulfide != disulfides_possible.end();
				++new_disulfide) {

			//add the configuration with the new disulfide
			DisulfideList new_disulfides_formed = disulfides_formed;
			new_disulfides_formed.push_back(*new_disulfide);
			final_configurations.push_back(new_disulfides_formed);

			//storage for possible disulfides which do not clash with the new one we have chosen
			DisulfideList disulfides_to_be_added;

			//identify new secondary disulfides which do not clash with the primary
			DisulfideList::const_iterator potential_further_disulfide = new_disulfide;
			for ( potential_further_disulfide++; potential_further_disulfide != disulfides_possible.end();
					++potential_further_disulfide ) {

				if ((*potential_further_disulfide).first != (*new_disulfide).first &&
						(*potential_further_disulfide).second != (*new_disulfide).first &&
						(*potential_further_disulfide).first != (*new_disulfide).second &&
						(*potential_further_disulfide).second != (*new_disulfide).second) {

					disulfides_to_be_added.push_back(*potential_further_disulfide);
				}

			} //end identification of possible secondary disulfides

			if ( disulfides_to_be_added.size() > 0 ) { //Add new disulfides if new ones can be added

				utility::vector1< DisulfideList > new_disulfide_configurations =
					recursive_multiple_disulfide_former(new_disulfides_formed, disulfides_to_be_added);

				for ( utility::vector1< DisulfideList >::const_iterator new_configuration = new_disulfide_configurations.begin();
						new_configuration != new_disulfide_configurations.end();
						++new_configuration) {
					final_configurations.push_back(*new_configuration);
				}
			} //end adding new disulfides recursively AFTER the first selected new one
		} //end addition of primary disulfide + all secondary possibilities
	}
	return final_configurations;
}

void
add_to_list( DisulfidizeMover::DisulfideList & disulf_partners, core::Size const r1, core::Size const r2 )
{
	std::pair< core::Size, core::Size > const temp_pair = std::make_pair( r1, r2 );
	std::pair< core::Size, core::Size > const alt_pair = std::make_pair( r2, r1 );
	if (std::find(disulf_partners.begin(), disulf_partners.end(), temp_pair) == disulf_partners.end() &&
			std::find(disulf_partners.begin(), disulf_partners.end(), alt_pair) == disulf_partners.end()) {
		disulf_partners.push_back( temp_pair );
	}
}

/// @brief find disulfides in the given neighborhood between residues in set 1 and residues in set 2 
DisulfidizeMover::DisulfideList
DisulfidizeMover::find_disulfides_in_the_neighborhood(
		core::pose::Pose const & pose,
		core::pack::task::residue_selector::ResidueSubset const & residueset1,
		core::pack::task::residue_selector::ResidueSubset const & residueset2 ) const
{
	TR << "FINDING DISULF" << std::endl;

	// figure out which positions are "central" positions - I presume these are positions from which DS bonds can emanate.
	// then figure out which positions are not "central" but still "modeled". I assume these are the disulfide landing range
	// positions.
	std::set< core::Size > resid_set1;
	std::set< core::Size > resid_set2;
	TR << "subset1: [ ";
	for ( core::Size ii=1, endii=residueset1.size(); ii<=endii; ++ii ) {
		if ( residueset1[ ii ] ) {
			TR << ii << " ";
			resid_set1.insert( ii );
		}
	}
	TR << "]" << std::endl;

	TR << "subset2: [ ";
	for ( core::Size ii=1, endii=residueset2.size(); ii<=endii; ++ii ) {
		if ( residueset2[ ii ] ) {
			TR << ii << " ";
			resid_set2.insert( ii );
		}
	}
	TR << "]" << std::endl;
	return find_disulfides_in_the_neighborhood( pose, resid_set1, resid_set2 );
}

/// @brief find disulfides in the given neighborhood between residues in set 1 and residues in set 2 
DisulfidizeMover::DisulfideList
DisulfidizeMover::find_disulfides_in_the_neighborhood(
		core::pose::Pose const & pose,
		std::set< core::Size > const & set1,
		std::set< core::Size > const & set2 ) const
{
	DisulfideList disulf_partners;

	core::scoring::disulfides::DisulfideMatchingPotential disulfPot;

	core::pose::Pose pose_copy = pose;
	construct_poly_ala_pose( pose_copy, false );

	core::scoring::ScoreFunctionOP sfxn_disulfide_only = core::scoring::ScoreFunctionOP(new core::scoring::ScoreFunction());
	sfxn_disulfide_only->set_weight(core::scoring::dslf_fa13, 1.0);

	for ( std::set< core::Size >::const_iterator itr=set1.begin(), end=set1.end(); itr!=end; itr++ ) {
		//gly/non-protein check
		if ( !check_residue_type( pose, *itr ) ) {
			continue;
		}
		for ( std::set< core::Size >::const_iterator itr2=set2.begin(), end2=set2.end(); itr2!=end2 ; itr2++ ) {
			if ( *itr == *itr2 ) {
				continue;
			}

			TR.Debug << "DISULF trying disulfide between " << *itr << " and " << *itr2 << std::endl;

			// gly/non-protein check
			if ( !check_residue_type( pose, *itr2 ) ) {
				continue;
			}

			// existing disulfides check
			if ( keep_current_ds_ && ( pose.residue(*itr).name() == "CYD" || pose.residue(*itr2).name() == "CYD" ) ) {
				TR <<"DISULF \tkeep_current_ds set to True, keeping residues that are already in disulfides." << std::endl;
				continue;
			}

			// seqpos distance check
			if ( !check_disulfide_seqpos(*itr, *itr2) ) {
				continue;
			}

			// cartesian CB-CB distance check
			if ( !check_disulfide_cb_distance( pose, *itr, *itr2 ) ) {
				continue;
			}

			// disulfide score check
			if ( !check_disulfide_score( pose_copy, *itr, *itr2, sfxn_disulfide_only ) ) {
				if ( include_current_ds_ && pose.residue(*itr).is_bonded( pose.residue(*itr2) ) ) {
					TR << "DISULF \tIncluding pre-existing disulfide despite failed disulf_fa_max check." << std::endl;
					add_to_list( disulf_partners, *itr, *itr2 );
				}
				continue;
			}

			// disulfide potential scoring
			if ( !check_disulfide_match_rt( pose, *itr, *itr2, disulfPot ) ) {
				if ( include_current_ds_ && pose.residue(*itr).is_bonded(pose.residue(*itr2)) ) {
					TR << "DISULF \tIncluding pre-existing disulfide despite failed match_rt_limit check." << std::endl;
					add_to_list( disulf_partners, *itr, *itr2 );
				}
				continue;
			}

			TR << "DISULF " <<  *itr << "x" << *itr2 << std::endl;
			add_to_list( disulf_partners, *itr, *itr2 );
		}
	}
	return disulf_partners;
}

/// @brief forms a disulfide between res1 and res2, optionally allowing backbone movement
void
DisulfidizeMover::make_disulfide(
		core::pose::Pose & pose,
		core::Size const res1,
		core::Size const res2,
		bool const relax_bb ) const
{
	TR << "build_disulf between " << res1 << " and " << res2 << std::endl;
	// create movemap which allows only chi to move
	core::kinematics::MoveMapOP mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap());
	mm->set_bb( res1, relax_bb );
	mm->set_bb( res2, relax_bb );
	mm->set_chi( res1, true );
	mm->set_chi( res2, true );

	core::conformation::form_disulfide( pose.conformation(), res1, res2 );
	core::util::rebuild_disulfide( pose, res1, res2,
			NULL, //task
			NULL, //sfxn
			mm, // movemap
			NULL ); // min sfxn
}

/// @brief creates disulfides given the list of pairs given
void
DisulfidizeMover::make_disulfides(
		core::pose::Pose & pose,
		DisulfidizeMover::DisulfideList const & disulf,
		bool const relax_bb ) const
{
	for ( DisulfideList::const_iterator itr=disulf.begin(); itr != disulf.end(); ++itr) {
		make_disulfide( pose, (*itr).first, (*itr).second, relax_bb );
	}
}

/// @brief temporarily tries building a disulfide between the given positions, scores, and restores the pose
core::Real
DisulfidizeMover::build_and_score_disulfide(
		core::pose::Pose & blank_pose,
		core::scoring::ScoreFunctionOP sfxn,
		const bool relax_bb,
		core::Size const res1,
		core::Size const res2 ) const
{
	assert( sfxn );
	assert( res1 );
	assert( res2 );
	assert( res1 <= blank_pose.total_residue() );
	assert( res2 <= blank_pose.total_residue() );
	assert( res1 != res2 );

	TR << "building and scoring " << res1 << " to " << res2 << std::endl;
	// save existing residues
	core::conformation::Residue old_res1 = blank_pose.residue(res1);
	core::conformation::Residue old_res2 = blank_pose.residue(res2);

	make_disulfide( blank_pose, res1, res2, relax_bb );

	core::Real const score = (*sfxn)(blank_pose) * 0.50;

	// restore saved residues
	blank_pose.replace_residue(res1, old_res1, false);
	blank_pose.replace_residue(res2, old_res2, false);

	return score;
}

/// @brief checks residue type to ensure that we can mutate this to CYD
bool
DisulfidizeMover::check_residue_type( core::pose::Pose const & pose, core::Size const res ) const
{
	bool const retval = ( pose.residue(res).is_protein() &&
		( pose.residue(res).aa() != core::chemical::aa_gly ) );
	if ( !retval ) {
			TR.Debug << "DISULF \tSkipping residue " << res << ". Residue of this type (" << pose.residue(res).name() << " ) cannot be mutated to CYD." << std::endl;
	}
	return retval;
}

/// @brief checks seqpos to ensure that min_loop is satisfied
bool
DisulfidizeMover::check_disulfide_seqpos( core::Size const res1, core::Size const res2 ) const
{
	core::SSize const seqpos_diff = std::abs( static_cast< core::SSize >(res2) - static_cast< core::SSize >(res1) );
	bool const retval = ( ( seqpos_diff + 1 ) >= std::abs( static_cast< core::SSize >(min_loop_) ) );
	if ( !retval ) {
		TR.Debug << "Residues " << res1 << " and " << res2 << " are too close in sequence (min_loop=" << min_loop_ << "). Skipping." << std::endl;
	}
	return retval;
}

/// @brief checks disulfide CB-CB distance
bool
DisulfidizeMover::check_disulfide_cb_distance(
		core::pose::Pose const & pose,
		core::Size const res1,
		core::Size const res2 ) const
{
	core::Real const dist_squared = pose.residue(res1).nbr_atom_xyz().distance_squared(pose.residue(res2).nbr_atom_xyz());
	bool const retval = ( dist_squared <= 25 );
	if ( !retval ) {
		TR.Debug << "DISULF \tTOO FAR. CB-CB distance squared: " << dist_squared << std::endl;
	}
	return retval;
}

/// @brief checks disulfide rosetta score
bool
DisulfidizeMover::check_disulfide_score(
		core::pose::Pose & pose,
		core::Size const res1,
		core::Size const res2,
		core::scoring::ScoreFunctionOP sfxn ) const
{
	core::Real const disulfide_fa_score =
		build_and_score_disulfide( pose, sfxn,
				false, // relax bb
				res1, res2 );
	TR << "DISULF FA SCORE RES " << res1 << " " << res2 << " " << disulfide_fa_score << std::endl;
	bool const retval = ( disulfide_fa_score <= max_disulf_score_ );
	if ( !retval ) {
		TR << "DISULF \tFailed disulf_fa_max check." << std::endl;
	}
	return retval;
}

/// @brief checks disulfide match rt
bool 
DisulfidizeMover::check_disulfide_match_rt(
			core::pose::Pose const & pose,
			core::Size const res1,
			core::Size const res2,
			core::scoring::disulfides::DisulfideMatchingPotential const & disulfPot ) const
{
	core::Energy match_t = 0.0;
	core::Energy match_r = 0.0;
	core::Energy match_rt = 0.0;
	disulfPot.score_disulfide( pose.residue(res1), pose.residue(res2), match_t, match_r, match_rt );
	TR << "DISULF \tmatch_t: " << match_t << ", match_r: " << match_r << ", match_rt: " << match_rt << std::endl;
	bool const retval = ( match_rt <= match_rt_limit_  );
	if ( !retval ) {
		TR << "DISULF \tFailed match_rt_limit check." << std::endl;
	}
	return retval;
}

} // namespace denovo_design
} // namespace protocols
