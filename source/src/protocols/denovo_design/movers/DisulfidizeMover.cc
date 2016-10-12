// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/movers/DisulfidizeMover.cc
/// @brief The DisulfidizeMover
/// @details
/// @author Tom Linsky (tlinsky@uw.edu) -- Adapting code from remodelmover into a mover
/// @author Gabe Rocklin (grocklin@uw.edu) -- Disulfide code
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Generalizing for D-cysteine, plus various bugfixes.

//Unit Headers
#include <protocols/denovo_design/movers/DisulfidizeMover.hh>
#include <protocols/denovo_design/movers/DisulfidizeMoverCreator.hh>

//Project Headers
#include <protocols/denovo_design/util.hh>

//Protocol Headers
#include <protocols/forge/remodel/RemodelDesignMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>

//Core Headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/util/disulfide_util.hh>
#include <protocols/rosetta_scripts/util.hh>

//Basic Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/stream_util.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

//C++ Headers

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.DisulfidizeMover" );

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

/// @brief Default constructor
///
DisulfidizeMover::DisulfidizeMover() :
	protocols::rosetta_scripts::MultiplePoseMover(),
	match_rt_limit_( 2.0 ),
	max_disulf_score_( 1.5 ),
	max_dist_sq_( 25 ),
	min_loop_( 8 ),
	min_disulfides_( 1 ),
	max_disulfides_( 3 ),
	include_current_ds_( false ),
	keep_current_ds_( false ),
	score_or_matchrt_( false ),
	allow_l_cys_( true ),
	allow_d_cys_( false ),
	mutate_gly_( false ),
	mutate_pro_( false ),
	set1_selector_(),
	set2_selector_(),
	sfxn_()
{
	set_rosetta_scripts_tag( utility::tag::TagOP( new utility::tag::Tag() ) );
}

/// @brief Copy constructor
///
DisulfidizeMover::DisulfidizeMover( DisulfidizeMover const &src ) :
	protocols::rosetta_scripts::MultiplePoseMover( src ),
	match_rt_limit_( src.match_rt_limit_ ),
	max_disulf_score_( src.max_disulf_score_ ),
	max_dist_sq_( src.max_dist_sq_ ),
	min_loop_( src.min_loop_ ),
	min_disulfides_( src.min_disulfides_ ),
	max_disulfides_( src.max_disulfides_ ),
	include_current_ds_( src.include_current_ds_ ),
	keep_current_ds_( src.keep_current_ds_ ),
	score_or_matchrt_( src.score_or_matchrt_ ),
	allow_l_cys_( src.allow_l_cys_ ),
	allow_d_cys_( src.allow_d_cys_ ),
	mutate_gly_( src.mutate_gly_ ),
	mutate_pro_( src.mutate_pro_ ),
	set1_selector_( src.set1_selector_ ), //Copies the const owning pointer -- points to same object!
	set2_selector_( src.set2_selector_ ), //Copies the const owning pointer -- points to same object!
	sfxn_() //Cloned below, if it exists.  NULL pointer is possible, too.
{
	if ( src.sfxn_ ) {
		sfxn_=src.sfxn_->clone(); //Copy the source scorefunction, if it exists.
	}
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
DisulfidizeMover::set_set1_selector( core::select::residue_selector::ResidueSelectorCOP selector )
{
	set1_selector_ = selector;
}

/// @brief sets the selector for set 2 -- disulfides will connect residues in set 1 to residues in set 2
void
DisulfidizeMover::set_set2_selector( core::select::residue_selector::ResidueSelectorCOP selector )
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

/// @brief Set the types of cysteines that we design with:
/// @details By default, we use only L-cysteine (not D-cysteine).
void
DisulfidizeMover::set_cys_types( bool const lcys, bool const dcys )
{
	allow_l_cys_ = lcys;
	allow_d_cys_ = dcys;
	runtime_assert_string_msg( allow_l_cys_ || allow_d_cys_, "Error in protocols::denovo_design::movers::DisulfidizeMover::set_cys_types():  The Disulfidize mover must use at least one of L-cysteine or D-cysteine.  The user has specified that both are prohibited.  Check your input, please." );
	return;
}

/// @brief Set the scorefunction to use for scoring disulfides, minimizing, and repacking.
/// @details Clones the input scorefunction; does not copy it.
void
DisulfidizeMover::set_scorefxn( core::scoring::ScoreFunctionCOP sfxn_in )
{
	sfxn_ = sfxn_in->clone();
	return;
}

/// @brief If true, GLY --> CYS mutations will be allowed. Default=false
void
DisulfidizeMover::set_mutate_gly( bool const mutate_gly )
{
	mutate_gly_ = mutate_gly;
}

/// @brief If true, PRO --> CYS mutations will be allowed. Default=false
void
DisulfidizeMover::set_mutate_pro( bool const mutate_pro )
{
	mutate_pro_ = mutate_pro;
}

/// @brief Given a list of disulfides and a symmetric pose, prune the list to remove symmetry
/// duplicates.
/// @details Does nothing if the pose is not symmetric.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu).
void
DisulfidizeMover::prune_symmetric_disulfides (
	core::pose::Pose const & pose,
	DisulfideList & disulf_list ) const
{
	//Number of disulfides in the input list:
	core::Size const ndisulf(disulf_list.size());
	if ( ndisulf < 2 ) return; //Nothing to prune.

	//Output list (starts empty):
	DisulfideList disulf_list_out;

	//Get an owning pointer to the symmetric conformation:
	core::conformation::symmetry::SymmetricConformationCOP symm_conf( utility::pointer::dynamic_pointer_cast< core::conformation::symmetry::SymmetricConformation const >( pose.conformation_ptr() ) );
	if ( !symm_conf ) return; //Do nothing if the conformation is not symmetric.

	for ( core::Size i=1; i<ndisulf; ++i ) { //Loop through the input list, starting with the second entry.
		bool n1_is_indep = symm_conf->Symmetry_Info()->bb_is_independent(  disulf_list[i].first );
		bool n2_is_indep = symm_conf->Symmetry_Info()->bb_is_independent(  disulf_list[i].second );

		if ( n1_is_indep || n2_is_indep ) {
			disulf_list_out.push_back( disulf_list[i] );
		}
	}

	disulf_list = disulf_list_out; //Replace the old disulf list.

	return;
}

void
DisulfidizeMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	if ( tag->hasOption( "scorefxn" ) ) {
		set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	}
	set_match_rt_limit( tag->getOption< core::Real >( "match_rt_limit", match_rt_limit_ ) );
	min_disulfides_ = tag->getOption< core::Size >( "min_disulfides", min_disulfides_ );
	max_disulfides_ = tag->getOption< core::Size >( "max_disulfides", max_disulfides_ );
	keep_current_ds_ = tag->getOption< bool >( "keep_current_disulfides", keep_current_ds_ );
	include_current_ds_ = tag->getOption< bool >( "include_current_disulfides", include_current_ds_ );
	set_min_loop( tag->getOption< core::Size >( "min_loop", min_loop_  ) );
	set_max_disulf_score( tag->getOption< core::Real >( "max_disulf_score", max_disulf_score_ ) );
	set_mutate_gly( tag->getOption< bool >( "mutate_gly", mutate_gly_ ) );
	set_mutate_pro( tag->getOption< bool >( "mutate_pro", mutate_pro_ ) );

	if ( tag->hasOption( "max_cb_dist" ) ) {
		core::Real const max_dist = tag->getOption< core::Real >( "max_cb_dist" );
		max_dist_sq_ = max_dist*max_dist;
	}

	// by default, a disulfide is valid if it passes score OR matchrt
	// if this option is false, disulfides must pass score AND matchrt
	score_or_matchrt_ = tag->getOption< bool >( "score_or_matchrt", score_or_matchrt_ );
	if ( tag->hasOption( "set1" ) ) {
		set_set1_selector( get_residue_selector( data, tag->getOption< std::string >( "set1" ) ) );
	}
	if ( tag->hasOption( "set2" ) ) {
		set_set2_selector( get_residue_selector( data, tag->getOption< std::string >( "set2" ) ) );
	}

	//By default, we only design with L-cystine:
	set_cys_types(
		tag->getOption< bool >( "use_l_cys", allow_l_cys_ ),
		tag->getOption< bool >( "use_d_cys", allow_d_cys_ ) );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief does the work -- scans pose for disulfides, stores first result in pose, adds others to additional_poses
bool
DisulfidizeMover::process_pose(
	core::pose::Pose & pose,
	utility::vector1 < core::pose::PoseOP > & additional_poses )
{
	bool const is_symmetric( core::pose::symmetry::is_symmetric(pose) ); //Is the pose symmetric?

	// remove any jump atoms from the fold tree which might cause problems when converting pose to poly-ala
	pose.fold_tree( remove_all_jump_atoms( pose.fold_tree() ) );

	core::scoring::ScoreFunctionOP sfxn( sfxn_ ? sfxn_->clone() : core::scoring::get_score_function() );
	if ( is_symmetric ) {
		runtime_assert_string_msg( core::pose::symmetry::is_symmetric(*sfxn), "Error in protocols::denovo_design::movers::DisulfidizeMover::process_pose().  A symmetric scorefunction must be provided for symmetric poses." );
	}

	// get two sets of residues which will be connected by disulfides
	core::select::residue_selector::ResidueSubset subset1;
	core::select::residue_selector::ResidueSubset subset2;
	if ( set1_selector_ ) {
		subset1 = set1_selector_->apply( pose );
	} else {
		subset1.assign( pose.size(), true );
	}
	if ( set2_selector_ ) {
		subset2 = set2_selector_->apply( pose );
	} else {
		subset2.assign( pose.size(), true );
	}

	DisulfideList current_ds = find_current_disulfides( pose, subset1, subset2 );
	if ( current_ds.size() > 0 ) { TR << "Current disulfides are: " << current_ds << std::endl; }
	else { TR << "No disulfides were already present in the pose." << std::endl; }

	if ( !keep_current_ds_ ) {
		mutate_disulfides_to_ala( pose, current_ds ); //Updated for D-cys, VKM 17 Aug 2015.
	}

	// create initial list of possible disulfides between residue subset 1 and subset 2
	DisulfideList disulf_partners = find_possible_disulfides( pose, subset1, subset2, sfxn ); //Updated for D-residues
	if ( include_current_ds_ ) {
		for ( DisulfideList::const_iterator ds=current_ds.begin(); ds!=current_ds.end(); ++ds ) {
			disulf_partners.push_back( *ds );
		}
	}
	if ( disulf_partners.size() > 0 ) TR << "Potential disulfides are: " << disulf_partners << std::endl;
	else TR << "No potential disulfides found." << std::endl;

	//If the pose is symmetric, prune the disulfide list to eliminate redundant disulfide pairs:
	if ( disulf_partners.size() > 0  && is_symmetric ) {
		prune_symmetric_disulfides( pose, disulf_partners );
		TR << "After pruning disulfide symmetry copies, potential disulfides are: " << disulf_partners << std::endl;
	}

	if ( disulf_partners.size() == 0 ) {
		if ( min_disulfides_ == 0 ) {
			return true;
		}
		if ( TR.visible() ) TR << "Failed to build any disulfides." << std::endl;
		return false;
	}

	//Use the recursive multiple disulfide former
	DisulfideList empty_disulfide_list;
	utility::vector1< DisulfideList > disulfide_configurations =
		recursive_multiple_disulfide_former( empty_disulfide_list, disulf_partners );

	if ( TR.visible() ) TR << "disulfide_configurations=" << disulfide_configurations << std::endl;
	PoseList results;
	if ( min_disulfides_ == 0 ) {
		results.push_back( pose.clone() );
	}

	// iterate over disulfide configurations
	for ( utility::vector1< DisulfideList >::const_iterator ds_config = disulfide_configurations.begin();
			ds_config != disulfide_configurations.end();
			++ds_config ) {
		if ( (*ds_config).size() >= min_disulfides_ && (*ds_config).size() <= max_disulfides_ ) {

			//form all the disulfides in the disulfide configuration
			if ( TR.visible() ) {
				TR << "Building disulfide configuration ";
				for ( DisulfideList::const_iterator my_ds = (*ds_config).begin(); my_ds != (*ds_config).end(); ++my_ds ) {
					TR << (*my_ds).first << "-" << (*my_ds).second << " ";
				}
				TR << std::endl;
			}

			core::pose::PoseOP disulf_copy_pose( pose.clone() );
			make_disulfides( *disulf_copy_pose, *ds_config, false, sfxn ); //Should be D-residue compatible.
			tag_disulfides( *disulf_copy_pose, *ds_config );
			results.push_back( disulf_copy_pose );
		}
	}

	if ( TR.visible() ) TR << "Found " << results.size() << " total results." << std::endl;
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

/// @brief finds disulfides within a pose subset
DisulfidizeMover::DisulfideList
DisulfidizeMover::find_current_disulfides(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & subset1,
	core::select::residue_selector::ResidueSubset const & subset2 ) const
{
	debug_assert( pose.size() == subset1.size() );
	debug_assert( pose.size() == subset2.size() );
	DisulfideList retval;
	std::set< core::Size > cyds;
	for ( core::Size resi=1; resi<=pose.size(); ++resi ) {
		if ( ( !subset1[ resi ] ) && ( !subset2[ resi ] ) ) continue;
		if ( pose.residue(resi).type().is_disulfide_bonded() ) {
			cyds.insert( resi );
		}
	}

	for ( std::set< core::Size >::const_iterator cyd1=cyds.begin(); cyd1!=cyds.end(); ++cyd1 ) {
		for ( std::set< core::Size >::const_iterator cyd2=cyd1; cyd2!=cyds.end(); ++cyd2 ) {
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
	if ( TR.visible() ) TR << "Mutating current disulfides to ALA" << std::endl;
	// mutate current disulfides to alanine if we aren't keeping or including them
	for ( DisulfideList::const_iterator ds=current_ds.begin(), endds=current_ds.end(); ds!=endds; ++ds ) {

		core::conformation::break_disulfide( pose.conformation(), ds->first, ds->second ); //First, break the existing disulfide

		if ( !pose.residue(ds->first).type().is_d_aa() ) {
			protocols::simple_moves::MutateResidue mut( ds->first, "ALA" );
			mut.apply( pose );
		} else {
			protocols::simple_moves::MutateResidue mut( ds->first, "DALA" );
			mut.apply( pose );
		}

		if ( !pose.residue(ds->second).type().is_d_aa() ) {
			protocols::simple_moves::MutateResidue mut2( ds->second, "ALA" );
			mut2.apply( pose );
		} else {
			protocols::simple_moves::MutateResidue mut2( ds->second, "DALA" );
			mut2.apply( pose );
		}

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
				++new_disulfide ) {

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

				if ( (*potential_further_disulfide).first != (*new_disulfide).first &&
						(*potential_further_disulfide).second != (*new_disulfide).first &&
						(*potential_further_disulfide).first != (*new_disulfide).second &&
						(*potential_further_disulfide).second != (*new_disulfide).second ) {

					disulfides_to_be_added.push_back(*potential_further_disulfide);
				}

			} //end identification of possible secondary disulfides

			if ( disulfides_to_be_added.size() > 0 ) { //Add new disulfides if new ones can be added

				utility::vector1< DisulfideList > new_disulfide_configurations =
					recursive_multiple_disulfide_former(new_disulfides_formed, disulfides_to_be_added);

				for ( utility::vector1< DisulfideList >::const_iterator new_configuration = new_disulfide_configurations.begin();
						new_configuration != new_disulfide_configurations.end();
						++new_configuration ) {
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
	if ( std::find(disulf_partners.begin(), disulf_partners.end(), temp_pair) == disulf_partners.end() &&
			std::find(disulf_partners.begin(), disulf_partners.end(), alt_pair) == disulf_partners.end() ) {
		disulf_partners.push_back( temp_pair );
	}
}

/// @brief find disulfides in the given neighborhood between residues in set 1 and residues in set 2
DisulfidizeMover::DisulfideList
DisulfidizeMover::find_possible_disulfides(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & residueset1,
	core::select::residue_selector::ResidueSubset const & residueset2,
	core::scoring::ScoreFunctionOP sfxn
) const {
	if ( TR.visible() ) TR << "FINDING DISULF" << std::endl;

	// figure out which positions are "central" positions - I presume these are positions from which DS bonds can emanate.
	// then figure out which positions are not "central" but still "modeled". I assume these are the disulfide landing range
	// positions.
	std::set< core::Size > resid_set1;
	std::set< core::Size > resid_set2;
	if ( TR.visible() ) TR << "subset1: [ ";
	for ( core::Size ii=1, endii=residueset1.size(); ii<=endii; ++ii ) {
		if ( residueset1[ ii ] ) {
			if ( TR.visible() ) TR << ii << " ";
			resid_set1.insert( ii );
		}
	}
	if ( TR.visible() ) TR << "]" << std::endl;

	if ( TR.visible() ) TR << "subset2: [ ";
	for ( core::Size ii=1, endii=residueset2.size(); ii<=endii; ++ii ) {
		if ( residueset2[ ii ] ) {
			if ( TR.visible() ) TR << ii << " ";
			resid_set2.insert( ii );
		}
	}
	if ( TR.visible() ) TR << "]" << std::endl;
	return find_possible_disulfides( pose, resid_set1, resid_set2, sfxn );
}

/// @brief find disulfides in the given neighborhood between residues in set 1 and residues in set 2
DisulfidizeMover::DisulfideList
DisulfidizeMover::find_possible_disulfides(
	core::pose::Pose const & pose,
	std::set< core::Size > const & set1,
	std::set< core::Size > const & set2,
	core::scoring::ScoreFunctionOP sfxn
) const {
	DisulfideList disulf_partners;

	// for "match-rt" scoring
	core::scoring::disulfides::DisulfideMatchingPotential disulfPot;

	// work on a poly-ala copy of the input pose
	core::pose::Pose pose_copy = pose;

	// combine sets
	std::set< core::Size > all_residues = set1;
	all_residues.insert( set2.begin(), set2.end() );
	construct_poly_ala_pose( pose_copy, keep_current_ds_, all_residues ); //Updated to work with D-residues.

	core::scoring::ScoreFunctionOP sfxn_disulfide_only (
		core::pose::symmetry::is_symmetric(pose) ?
		core::scoring::ScoreFunctionOP( new core::scoring::symmetry::SymmetricScoreFunction() ) :
		core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction() )
	);
	sfxn_disulfide_only->set_weight(core::scoring::dslf_fa13, 1.0);

	for ( std::set< core::Size >::const_iterator itr=set1.begin(), end=set1.end(); itr!=end; ++itr ) {
		//gly/pro/non-protein check
		if ( !check_residue_type( pose, *itr ) ) {
			continue;
		}
		for ( std::set< core::Size >::const_iterator itr2=set2.begin(), end2=set2.end(); itr2!=end2; ++itr2 ) {
			if ( *itr == *itr2 ) {
				continue;
			}

			if ( TR.Debug.visible() ) TR.Debug << "DISULF trying disulfide between " << *itr << " and " << *itr2 << std::endl;

			// gly/pro/non-protein check
			if ( !check_residue_type( pose, *itr2 ) ) {
				continue;
			}

			// existing disulfides check
			if ( keep_current_ds_ && ( pose.residue(*itr).type().is_disulfide_bonded() || pose.residue(*itr2).type().is_disulfide_bonded() ) ) {
				if ( TR.visible() ) TR <<"DISULF \tkeep_current_ds set to True, keeping residues that are already in disulfides." << std::endl;
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
			bool good_score = check_disulfide_score( pose_copy, *itr, *itr2, sfxn_disulfide_only, sfxn ); //Updated for D-residues

			// stop if we need good score AND matchrt
			if ( !score_or_matchrt_ && !good_score ) {
				continue;
			}

			// disulfide potential scoring -- check_disulfide_match_rt can now score mirror-image pairs, too
			bool good_match(false);
			bool const mixed( mixed_disulfide( pose, *itr, *itr2) );
			if ( !mixed ) {
				bool good_match1(false), good_match_inv(false);
				if ( allow_l_cys_ ) good_match1=check_disulfide_match_rt( pose, *itr, *itr2, disulfPot, false );
				if ( allow_d_cys_ ) good_match_inv=check_disulfide_match_rt( pose, *itr, *itr2, disulfPot, true );
				good_match = good_match1 || good_match_inv;
			} else {
				good_match=true; //Don't apply the match_rt test if this is a mixed disulfide
			}
			if ( ( !mixed ) && score_or_matchrt_ ) {
				if ( !good_score && !good_match ) {
					continue;
				}
			} else {
				if ( !good_score || !good_match ) {
					continue;
				}
			}

			if ( TR.visible() ) TR << "DISULF " <<  *itr << "x" << *itr2 << std::endl;
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
	std::stringstream disulf_comment;
	for ( DisulfideList::const_iterator ds=disulf.begin(); ds!=disulf.end(); ++ds ) {
		if ( !disulf_comment.str().empty() ) disulf_comment << ";";
		disulf_comment << ds->first << "," << ds->second;
	}
	core::pose::add_comment( pose, "DISULFIDIZE", disulf_comment.str() );
	for ( DisulfideList::const_iterator itr=disulf.begin(); itr != disulf.end(); ++itr ) {
		tag_disulfide( pose, (*itr).first, (*itr).second );
	}
}

/// @brief forms a disulfide between res1 and res2, optionally allowing backbone movement
void
DisulfidizeMover::make_disulfide(
	core::pose::Pose & pose,
	core::Size const res1,
	core::Size const res2,
	bool const relax_bb,
	core::scoring::ScoreFunctionOP sfxn ) const
{
	if ( TR.Debug.visible() ) TR.Debug << "build_disulf between " << res1 << " and " << res2 << std::endl;
	// create movemap which allows only chi to move
	core::kinematics::MoveMapOP mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap());
	mm->set_bb( res1, relax_bb );
	mm->set_bb( res2, relax_bb );
	mm->set_chi( res1, true );
	mm->set_chi( res2, true );

	core::conformation::form_disulfide( pose.conformation(), res1, res2, allow_d_cys_, !allow_l_cys_ ); //Updated for D-residues
	core::util::rebuild_disulfide( pose, res1, res2,
		NULL, //task
		sfxn, //sfxn
		mm, // movemap
		sfxn // min sfxn
	); //Seems already to work with D-residues
}

/// @brief creates disulfides given the list of pairs given
void
DisulfidizeMover::make_disulfides(
	core::pose::Pose & pose,
	DisulfidizeMover::DisulfideList const & disulf,
	bool const relax_bb,
	core::scoring::ScoreFunctionOP sfxn ) const
{
	for ( DisulfideList::const_iterator itr=disulf.begin(); itr != disulf.end(); ++itr ) {
		make_disulfide( pose, (*itr).first, (*itr).second, relax_bb, sfxn );
	}
}

/// @brief temporarily tries building a disulfide between the given positions, scores, and restores the pose
core::Real
DisulfidizeMover::build_and_score_disulfide(
	core::pose::Pose & blank_pose,
	core::scoring::ScoreFunctionOP sfxn_disulfonly,
	core::scoring::ScoreFunctionOP sfxn_full,
	const bool relax_bb,
	core::Size const res1,
	core::Size const res2 ) const
{
	assert( sfxn_disulfonly );
	assert( sfxn_full );
	assert( res1 );
	assert( res2 );
	assert( res1 <= blank_pose.size() );
	assert( res2 <= blank_pose.size() );
	assert( res1 != res2 );

	if ( TR.visible() ) TR << "building and scoring " << res1 << " to " << res2 << std::endl;
	// save existing residues
	core::conformation::Residue old_res1 = blank_pose.residue(res1);
	core::conformation::Residue old_res2 = blank_pose.residue(res2);

	make_disulfide( blank_pose, res1, res2, relax_bb, sfxn_full );

	core::Real const score = (*sfxn_disulfonly)(blank_pose) * 0.50;

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
		( mutate_pro_ ||
		( ( pose.residue(res).aa() != core::chemical::aa_pro ) && ( pose.residue(res).aa() != core::chemical::aa_dpr ) ) ) &&
		( mutate_gly_ || ( pose.residue(res).aa() != core::chemical::aa_gly ) ) );
	if ( !retval && TR.Debug.visible() ) {
		TR.Debug << "DISULF \tSkipping residue " << res << ". Residue of this type (" << pose.residue(res).name() << ") cannot be mutated to CYD." << std::endl;
	}
	return retval;
}

/// @brief checks seqpos to ensure that min_loop is satisfied
bool
DisulfidizeMover::check_disulfide_seqpos( core::Size const res1, core::Size const res2 ) const
{
	core::SSize const seqpos_diff = std::abs( static_cast< core::SSize >(res2) - static_cast< core::SSize >(res1) );
	bool const retval = ( ( seqpos_diff + 1 ) >= std::abs( static_cast< core::SSize >(min_loop_) ) );
	if ( !retval && TR.Debug.visible() ) {
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
	bool const retval = ( dist_squared <= max_dist_sq_ );
	if ( !retval && TR.Debug.visible() ) {
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
	core::scoring::ScoreFunctionOP sfxn_disulfonly,
	core::scoring::ScoreFunctionOP sfxn_full ) const
{
	core::Real const disulfide_fa_score ( build_and_score_disulfide( pose, sfxn_disulfonly, sfxn_full, false /*relax_bb*/, res1, res2 ) );
	if ( TR.visible() ) TR << "DISULF FA SCORE RES " << res1 << " " << res2 << " " << disulfide_fa_score << std::endl;
	bool const retval = ( disulfide_fa_score <= max_disulf_score_ );
	if ( TR.visible() && !retval ) {
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
	core::scoring::disulfides::DisulfideMatchingPotential const & disulfPot,
	bool const mirror ) const
{
	core::Energy match_t = 0.0;
	core::Energy match_r = 0.0;
	core::Energy match_rt = 0.0;
	disulfPot.score_disulfide( pose.residue(res1), pose.residue(res2), match_t, match_r, match_rt, mirror );
	if ( TR.visible() ) TR << "DISULF \tmatch_t: " << match_t << ", match_r: " << match_r << ", match_rt: " << match_rt << std::endl;
	bool const retval = ( match_rt <= match_rt_limit_  );
	if ( !retval && TR.visible() ) {
		TR << "DISULF \tFailed match_rt_limit check." << std::endl;
	}
	return retval;
}

/// @brief Returns true if this is a mixed D/L disulfide, false otherwise.
///
bool DisulfidizeMover::mixed_disulfide (
	core::pose::Pose const & pose,
	core::Size const res1,
	core::Size const res2 ) const
{
	if ( allow_l_cys_ && !allow_d_cys_ ) return false;
	if ( allow_d_cys_ && !allow_l_cys_ ) return false;

	if ( pose.residue(res1).type().is_d_aa() && !pose.residue(res2).type().is_d_aa() ) return true;
	if ( !pose.residue(res1).type().is_d_aa() && pose.residue(res2).type().is_d_aa() ) return true;

	return false;
}

} // namespace denovo_design
} // namespace protocols
