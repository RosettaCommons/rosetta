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
/// @details
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
#include <core/chemical/ChemicalManager.hh>
#include <core/pack/task/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/util/disulfide_util.hh>

//Basic Headers
#include <basic/datacache/DataMap.hh>
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
	return "Disulfidize";
}

///  ---------------------------------------------------------------------------------
///  DisulfidizeMover main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
DisulfidizeMover::DisulfidizeMover() :
	protocols::rosetta_scripts::MultiplePoseMover(),
	match_rt_limit_( 2.0 ),
	max_disulf_score_( -0.25 ),
	min_loop_( 8 ),
	min_disulfides_( 1 ),
	max_disulfides_( 3 ),
	include_current_ds_( false ),
	keep_current_ds_( false ),
	score_or_matchrt_( true ),
	set1_selector_(),
	set2_selector_()
{
	set_rosetta_scripts_tag( utility::tag::TagOP( new utility::tag::Tag() ) );
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

/// @brief sets the selector for set 1 -- disulfides will connect residues in set 1 to residues in set 2
void
DisulfidizeMover::set_set1_selector( core::pack::task::residue_selector::ResidueSelectorCOP selector )
{
	set1_selector_ = selector;
}

/// @brief sets the selector for set 2 -- disulfides will connect residues in set 1 to residues in set 2
void
DisulfidizeMover::set_set2_selector( core::pack::task::residue_selector::ResidueSelectorCOP selector )
{
	set2_selector_ = selector;
}

/// @brief sets the min_loop value (number of residues between disulfide-joined residues) (default=8)
void
DisulfidizeMover::set_min_loop( core::Size const minloopval )
{
	min_loop_ = minloopval;
}

/// @brief sets the maximum allowed per-disulfide dslf_fa13 score (default=0.0)
void
DisulfidizeMover::set_max_disulf_score( core::Real const maxscoreval )
{
	max_disulf_score_ = maxscoreval;
}

/// @brief sets the maximum allowed "match-rt-limit" (default=2.0)
void
DisulfidizeMover::set_match_rt_limit( core::Real const matchrtval )
{
	match_rt_limit_ = matchrtval;
}

void
DisulfidizeMover::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & ,
		protocols::moves::Movers_map const & ,
		core::pose::Pose const & )
{
	if ( tag->hasOption( "match_rt_limit" ) )
		set_match_rt_limit( tag->getOption< core::Real >( "match_rt_limit" ) );
	if ( tag->hasOption( "min_disulfides" ) )
		min_disulfides_ = tag->getOption< core::Size >( "min_disulfides" );
	if ( tag->hasOption( "max_disulfides" ) )
		max_disulfides_ = tag->getOption< core::Size >( "max_disulfides" );
	if ( tag->hasOption( "keep_current_disulfides" ) )
		keep_current_ds_ = tag->getOption< bool >( "keep_current_disulfides" );
	if ( tag->hasOption( "include_current_disulfides" ) )
		include_current_ds_ = tag->getOption< bool >( "include_current_disulfides" );
	if ( tag->hasOption( "min_loop" ) )
		set_min_loop( tag->getOption< core::Size >( "min_loop" ) );
	if ( tag->hasOption( "max_disulf_score" ) )
		set_max_disulf_score( tag->getOption< core::Real >( "max_disulf_score" ) );
	// by default, a disulfide is valid if it passes score OR matchrt
	// if this option is false, disulfides must pass score AND matchrt
	if ( tag->hasOption( "score_or_matchrt" ) )
		score_or_matchrt_ = tag->getOption< bool >( "score_or_matchrt" );

	if ( tag->hasOption( "set1" ) ) {
		set_set1_selector( get_residue_selector( data, tag->getOption< std::string >( "set1" ) ) );
	}
	if ( tag->hasOption( "set2" ) ) {
		set_set2_selector( get_residue_selector( data, tag->getOption< std::string >( "set2" ) ) );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief does the work -- scans pose for disulfides, stores first result in pose, adds others to additional_poses
bool
DisulfidizeMover::process_pose(
		core::pose::Pose & pose,
		utility::vector1 < core::pose::PoseOP > & additional_poses )
{
	DisulfideList current_ds = find_current_disulfides( pose );
	TR << "Current disulfides are: " << current_ds << std::endl;

	if ( !keep_current_ds_ ) {
		mutate_disulfides_to_ala( pose, current_ds );
	}

	// get two sets of residues which will be connected by disulfides
	core::pack::task::residue_selector::ResidueSubset subset1( pose.total_residue(), true );
	core::pack::task::residue_selector::ResidueSubset subset2( pose.total_residue(), true );
	if ( set1_selector_ ) {
		subset1 = set1_selector_->apply( pose );
	}
	if ( set2_selector_ ) {
		subset2 = set2_selector_->apply( pose );
	}

	// create initial list of possible disulfides between residue subset 1 and subset 2
	DisulfideList disulf_partners = find_possible_disulfides( pose, subset1, subset2 );
	if ( include_current_ds_ ) {
		for ( DisulfideList::const_iterator ds=current_ds.begin(), endds=current_ds.end(); ds!=endds; ++ds ) {
			disulf_partners.push_back( *ds );
		}
	}

	if ( disulf_partners.size() == 0 ) {
		if ( min_disulfides_ == 0 ) {
			return true;
		}
		TR << "Failed to build any disulfides." << std::endl;
		return false;
	}

	//Use the recursive multiple disulfide former
	DisulfideList empty_disulfide_list;
	utility::vector1< DisulfideList > disulfide_configurations =
		recursive_multiple_disulfide_former( empty_disulfide_list, disulf_partners );

	TR << "disulfide_configurations=" << disulfide_configurations << std::endl;
	PoseList results;
	if ( min_disulfides_ == 0 ) {
		results.push_back( pose.clone() );
	}

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

			core::pose::PoseOP disulf_copy_pose( pose.clone() );
			make_disulfides( *disulf_copy_pose, *ds_config, false );
			tag_disulfides( *disulf_copy_pose, *ds_config );
			results.push_back( disulf_copy_pose );
		}
	}

	TR << "Found " << results.size() << " total results." << std::endl;
	if ( results.size() == 0 ) {
		return false;
	}

	core::Size count = 1;
	for ( PoseList::const_iterator p=results.begin(), endp=results.end(); p!=endp; ++p ) {
		if ( count == 1 ) {
			pose = **p;
		} else {
			additional_poses.push_back( *p );
		}
		++count;
	}
	return true;
}

/// @brief finds disulfides within a pose
DisulfidizeMover::DisulfideList
DisulfidizeMover::find_current_disulfides( core::pose::Pose const & pose ) const
{
	DisulfideList retval;
	std::set< core::Size > cyds;
	for ( core::Size i=1, endi=pose.total_residue(); i<=endi; ++i ) {
		if ( pose.residue(i).type().is_disulfide_bonded() )
			cyds.insert( i );
	}
	for ( std::set< core::Size >::const_iterator cyd1=cyds.begin(), endcyds=cyds.end(); cyd1!=endcyds; ++cyd1 ) {
		for ( std::set< core::Size >::const_iterator cyd2=cyd1; cyd2!=endcyds; ++cyd2 ) {
			if ( pose.residue(*cyd1).is_bonded( pose.residue(*cyd2) ) ) {
				retval.push_back( std::make_pair( *cyd1, *cyd2 ) );
			}
		}
	}
	return retval;
}

/// @brief mutates the given disulfides to ALA
void
DisulfidizeMover::mutate_disulfides_to_ala(
		core::pose::Pose & pose,
		DisulfideList const & current_ds ) const
{
	TR << "Mutating current disulfides to ALA" << std::endl;
	// mutate current disulfides to alanine if we aren't keeping or including them
	for ( DisulfideList::const_iterator ds=current_ds.begin(), endds=current_ds.end(); ds!=endds; ++ds ) {
		protocols::simple_moves::MutateResidue mut( ds->first, "ALA" );
		mut.apply( pose );
		protocols::simple_moves::MutateResidue mut2( ds->second, "ALA" );
		mut2.apply( pose );
	}
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
			DisulfideList::const_iterator potential_further_disulfide = new_disulfide, end = disulfides_possible.end();
			for ( ++potential_further_disulfide; potential_further_disulfide != end;
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
DisulfidizeMover::find_possible_disulfides(
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
	return find_possible_disulfides( pose, resid_set1, resid_set2 );
}

/// @brief find disulfides in the given neighborhood between residues in set 1 and residues in set 2
DisulfidizeMover::DisulfideList
DisulfidizeMover::find_possible_disulfides(
		core::pose::Pose const & pose,
		std::set< core::Size > const & set1,
		std::set< core::Size > const & set2 ) const
{
	DisulfideList disulf_partners;

	// for "match-rt" scoring
	core::scoring::disulfides::DisulfideMatchingPotential disulfPot;

	// work on a poly-ala copy of the input pose
	core::pose::Pose pose_copy = pose;
	construct_poly_ala_pose( pose_copy, false, set1, set2 );

	core::scoring::ScoreFunctionOP sfxn_disulfide_only = core::scoring::ScoreFunctionOP(new core::scoring::ScoreFunction());
	sfxn_disulfide_only->set_weight(core::scoring::dslf_fa13, 1.0);

	for ( std::set< core::Size >::const_iterator itr=set1.begin(), end=set1.end(); itr!=end; ++itr ) {
		//gly/non-protein check
		if ( !check_residue_type( pose, *itr ) ) {
			continue;
		}
		for ( std::set< core::Size >::const_iterator itr2=set2.begin(), end2=set2.end(); itr2!=end2; ++itr2 ) {
			if ( *itr == *itr2 ) {
				continue;
			}

			TR.Debug << "DISULF trying disulfide between " << *itr << " and " << *itr2 << std::endl;

			// gly/non-protein check
			if ( !check_residue_type( pose, *itr2 ) ) {
				continue;
			}

			// existing disulfides check
			if ( keep_current_ds_ && ( pose.residue(*itr).type().is_disulfide_bonded() || pose.residue(*itr2).type().is_disulfide_bonded() ) ) {
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

			// disulfide rosetta score
			bool good_score = check_disulfide_score( pose_copy, *itr, *itr2, sfxn_disulfide_only );

			// stop if we need good score AND matchrt
			if ( !score_or_matchrt_ && !good_score ) {
				continue;
			}

			// disulfide potential scoring
			bool good_match = check_disulfide_match_rt( pose, *itr, *itr2, disulfPot );
			if ( score_or_matchrt_ ) {
				if ( !good_score && !good_match ) {
					continue;
				}
			} else {
				if ( !good_score || !good_match ) {
					continue;
				}
			}

			TR << "DISULF " <<  *itr << "x" << *itr2 << std::endl;
			add_to_list( disulf_partners, *itr, *itr2 );
		}
	}
	return disulf_partners;
}

/// @brief creates a residue tags on disulfides to inform users that this disulfide was created by disulfidize
void
DisulfidizeMover::tag_disulfide(
		core::pose::Pose & pose,
		core::Size const res1,
		core::Size const res2 ) const
{
	if ( !pose.pdb_info() ) {
		pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( pose ) ) );
	}
	debug_assert( pose.pdb_info() );
	pose.pdb_info()->add_reslabel( res1, "DISULFIDIZE" );
	pose.pdb_info()->add_reslabel( res2, "DISULFIDIZE" );
}

/// @brief creates a residue tags on disulfides to inform users that this disulfide was created by disulfidize
void
DisulfidizeMover::tag_disulfides(
		core::pose::Pose & pose,
		DisulfidizeMover::DisulfideList const & disulf ) const
{
	for ( DisulfideList::const_iterator itr=disulf.begin(); itr != disulf.end(); ++itr) {
		tag_disulfide( pose, (*itr).first, (*itr).second );
	}
}

/// @brief forms a disulfide between res1 and res2, optionally allowing backbone movement
void
DisulfidizeMover::make_disulfide(
		core::pose::Pose & pose,
		core::Size const res1,
		core::Size const res2,
		bool const relax_bb ) const
{
	TR.Debug << "build_disulf between " << res1 << " and " << res2 << std::endl;
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
