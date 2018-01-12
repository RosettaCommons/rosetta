// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/cyclic_peptide/TryDisulfPermutations.cc
/// @brief  Tries all permutations of disulfide bonds between disulfide-forming residues in a pose.
/// @details  Performs a repack/minimize of
/// the disulfide-forming residues with a simplified score function.  Returns the permutation with
/// the lowest disulfide energy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/cyclic_peptide/TryDisulfPermutations.hh>
#include <protocols/cyclic_peptide/TryDisulfPermutationsCreator.hh>
#include <utility/tag/Tag.hh>

#include <numeric/constants.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <basic/datacache/DataMap.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using basic::Error;
using basic::Warning;

namespace protocols {
namespace cyclic_peptide {

static basic::Tracer TR("protocols.cyclic_peptide.TryDisulfPermutations");

// XRW TEMP std::string
// XRW TEMP TryDisulfPermutationsCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return TryDisulfPermutations::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP TryDisulfPermutationsCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new TryDisulfPermutations );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP TryDisulfPermutations::mover_name()
// XRW TEMP {
// XRW TEMP  return "TryDisulfPermutations";
// XRW TEMP }

/// @brief Creator for TryDisulfPermutations mover.
TryDisulfPermutations::TryDisulfPermutations():
	Mover("TryDisulfPermutations"),
	consider_already_bonded_(true),
	mintype_("dfpmin"),
	mintolerance_(0.00001),
	selector_()
	//TODO -- initialize data here.
{}


/// @brief Copy constructor for TryDisulfPermutations mover.
///
/// @brief Creator for TryDisulfPermutations mover.
TryDisulfPermutations::TryDisulfPermutations( TryDisulfPermutations const & src ):
	protocols::moves::Mover(src),
	consider_already_bonded_( src.consider_already_bonded_ ),
	mintype_( src.mintype_ ),
	mintolerance_( src.mintolerance_ ),
	selector_() //Cloned below
	//TODO -- copy data here.
{
	if ( src.selector_ ) selector_ = src.selector_->clone();
}


/// @brief Destructor for TryDisulfPermutations mover.
TryDisulfPermutations::~TryDisulfPermutations() = default;


/// @brief Clone operator to create a pointer to a fresh TryDisulfPermutations object that copies this one.
protocols::moves::MoverOP TryDisulfPermutations::clone() const {
	return protocols::moves::MoverOP( new TryDisulfPermutations( *this ) );
}


/// @brief Fresh_instance operator to create a pointer to a fresh TryDisulfPermutations object that does NOT copy this one.
protocols::moves::MoverOP TryDisulfPermutations::fresh_instance() const {
	return protocols::moves::MoverOP( new TryDisulfPermutations );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Actually apply the mover to the pose.
///
void TryDisulfPermutations::apply (
	core::pose::Pose & pose
) {
	core::Size const nres( pose.size() ); //Number of residues in the pose
	utility::vector1 < core::Size > disulf_res; //List of disulfide-forming residues

	core::select::residue_selector::ResidueSubset mask; //Used only if there's a ResidueSelector
	if ( selector_ ) { //If there's a ResidueSelector, apply it
		mask = selector_->apply(pose);
	}

	for ( core::Size ir=1; ir<=nres; ++ir ) {
		if ( selector_ && !mask[ir] ) continue; //Skip unselected residues, if there's a selector.
		if ( pose.residue(ir).type().get_disulfide_atom_name() == "NONE" ) continue;
		if ( !consider_already_bonded() && pose.residue(ir).type().is_disulfide_bonded() ) continue; //Skip residues that are already disulfide-bonded if this option is set.
		disulf_res.push_back( ir ) ; //This IS a disulfide former, so store it.
	}

	// Exit if there are no disulfide-forming residues:
	if ( disulf_res.size() == 0 ) {
		if ( TR.Warning.visible() ) {
			TR.Warning << "The TryDisulfPermutations mover was applied to a pose that has no disulfide-forming residues.  The input pose is returned." << std::endl;
			TR.Warning.flush();
		}
		return;
	}

	utility::vector1 < utility::vector1 < std::pair <core::Size, core::Size> > > disulf_permutations; // The disulfide permutations: the outer vector is the permutations, and the inner vector is a list of pairs of residues representing disulfide pairs.

	if ( disulf_res.size() % 2 == 0 ) { //If we have an even number of disulfide residues, we can generate permutations easily:
		generate_disulf_permutations( disulf_res, disulf_permutations ); //Generate the list of disulfide permutations from the list of disulfide-forming residues.
	} else { //If the number of disulfide formers is odd, we need to generate permuations for all subsets in which we leave one out.
		for ( core::Size i=1, imax=disulf_res.size(); i<=imax; ++i ) {
			if ( TR.Debug.visible() ) {
				TR.Debug << "Generating permutations -- omitting residue " << disulf_res[i] << "." << std::endl;
				TR.Debug.flush();
			}
			utility::vector1 < core::Size> disulf_res_omit_one;
			disulf_res_omit_one.reserve( imax - 1 );
			for ( core::Size j=1; j<=imax; ++j ) {
				if ( j!=i ) disulf_res_omit_one.push_back(disulf_res[j]); //Make a list that leaves out one disulfide-forming residue.
			}
			generate_disulf_permutations( disulf_res_omit_one, disulf_permutations ); //APPENDS permutations to the list already generated
		}
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Disulfide permutations:" << std::endl;
		for ( core::Size i=1, imax=disulf_permutations.size(); i<=imax; ++i ) {
			for ( core::Size j=1, jmax=disulf_permutations[i].size(); j<=jmax; ++j ) {
				TR.Debug << "[" << disulf_permutations[i][j].first << "," << disulf_permutations[i][j].second << "]";
				if ( j<jmax ) TR.Debug << " ";
			}
			TR.Debug << std::endl;
		}
		TR.Debug.flush();
	}


	core::pose::PoseOP lowest_energy_pose; //Pointer to the lowest-energy pose encountered so far.
	core::Real lowest_energy(0.0); //The lowest disulfide energy encountered so far.

	for ( core::Size i=1, imax=disulf_permutations.size(); i<=imax; ++i ) { //Loop through all permutations of disulfides.
		core::pose::PoseOP pose_copy( pose.clone() );

		generate_disulfides( pose_copy, disulf_permutations[i] ); //Set up the disulfides.

		core::Real const cur_energy( repack_minimize_disulfides( pose_copy, disulf_res ) ); //Repack and minimize the disulfide-forming residues, and return the energy.

		if ( i == 1 || cur_energy < lowest_energy ) { //If this is the first permutation or the lowest energy encountered so far, store it.
			lowest_energy = cur_energy;
			lowest_energy_pose = pose_copy;
		}
	}

	runtime_assert(lowest_energy_pose); //Ensure that the pointer points at something.  Should always be true.
	pose = (*lowest_energy_pose); //Return the lowest energy pose.

	return;
}

////////////////////////////////////////////////////////////////////////////////


/// @brief Returns the name of this mover ("TryDisulfPermutations").
// XRW TEMP std::string TryDisulfPermutations::get_name() const{
// XRW TEMP  return "TryDisulfPermutations";
// XRW TEMP }

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////

/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
void
TryDisulfPermutations::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &datamap,
	protocols::filters::Filters_map const &/*filters*/,
	protocols::moves::Movers_map const &/*movers*/,
	core::pose::Pose const & /*pose*/
) {

	if ( tag->getName() != "TryDisulfPermutations" ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible -- the tag name does not match the mover name.");
	}

	if ( TR.visible() ) TR << "Parsing options for TryDisulfPermutations (\"" << tag->getOption<std::string>("name" ,"") << "\") mover." << std::endl;

	if ( tag->hasOption("consider_already_bonded") ) {
		set_consider_already_bonded( tag->getOption<bool>( "consider_already_bonded" ) );
		if ( TR.visible() ) TR << "Mover will consider disulfides involving " << (consider_already_bonded() ? "all disulfide-forming residues, including those already in disulfides." : "only those disulfide-forming residues not already in disulfides." ) << std::endl ;
	}

	if ( tag->hasOption( "min_type" ) ) {
		set_mintype( tag->getOption<std::string>("min_type") );
		if ( TR.visible() ) TR << "Set minimization type to " << mintype() << "." << std::endl;
	}

	if ( tag->hasOption( "min_tolerance" ) ) {
		set_mintolerance( tag->getOption<core::Real>("min_tolerance") );
		if ( TR.visible() ) TR << "Set minimization tolerance to " << mintolerance() << std::endl;
	}

	if ( tag->hasOption( "selector" ) ) {
		std::string const selector_name ( tag->getOption< std::string >( "selector" ) );
		try {
			set_selector( datamap.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selector_name ) );
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the Datamap from ReadResfile::parse_tag()\n" + e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
		}
		if ( TR.visible() ) TR << "Added ResidueSelector \"" << selector_name << "\"." << std::endl;
	}


	if ( TR.visible() ) TR.flush();
	return;
} //parse_my_tag

/// @brief Set the residue selector.
/// @details CLONES the input.
void TryDisulfPermutations::set_selector( core::select::residue_selector::ResidueSelectorCOP selector_in )
{
	runtime_assert(selector_in); //Must point to something.
	selector_ = selector_in->clone();
	return;
}


////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

/// @brief Given a list of disulfide-forming residues, generate all possible permtutations.
/// @param[in] residues The list of disulfide-forming positions.
/// @param[out] permutations The list of lists of disulfide pairs representing all possible permutations.
/// @details Note that this function does NOT clear the permutations vector; it only appends to it.
void TryDisulfPermutations::generate_disulf_permutations(
	utility::vector1 < core::Size > const &residues,
	utility::vector1 < utility::vector1 < std::pair< core::Size, core::Size > > > &permutations
) const {
	core::Size const nres( residues.size() );
	runtime_assert_string_msg(nres % 2 == 0, "Error in protocols::cyclic_peptide::TryDisulfPermutations::generate_disulf_permutations(): The number of disulfide-forming residues must be even."); //A restriction that I'll remove eventually
	core::Size const npairs( nres / 2 );

	//Count number of permutations:
	core::Size npermutations( 1 );
	core::Size countvar( 1 );
	for ( core::Size i=1, imax=npairs; i<=imax; ++i ) {
		npermutations *= countvar;
		countvar += 2;
	}

	utility::vector1 < core::Size > permutation_index(npairs, 0);
	for ( core::Size i=1; i<=npermutations; ++i ) {
		utility::vector1 < bool > placed( residues.size(), false );
		utility::vector1 < std::pair < core::Size, core::Size > > cur_permutation;
		recursively_generate_permutation( cur_permutation, residues, placed, permutation_index, 1 );

		permutations.push_back(cur_permutation);

		increment_index( permutation_index, npairs, npairs );
		if ( TR.Debug.visible() ) {
			TR.Debug << "Permutation index: ";
			for ( core::Size i=1, imax=permutation_index.size(); i<=imax; ++i ) {
				TR.Debug << permutation_index[i];
				if ( i<imax ) TR.Debug << ",";
			}
			TR.Debug << std::endl;
			TR.Debug.flush();
		}
	}

	return;
}


/// @brief Given a permutation index, generate a permutation.
/// @details Assumes that cur_permutation is an empty vector.  Recursively calls itself.
void TryDisulfPermutations::recursively_generate_permutation(
	utility::vector1 < std::pair < core::Size, core::Size > > &cur_permutation,
	utility::vector1 < core::Size > const &residues,
	utility::vector1 < bool > &placed,
	utility::vector1 < core::Size > const &index,
	core::Size const &level
) const {
	runtime_assert( residues.size() == placed.size() ); //Should always be true.
	core::Size first(0), second(0), counter(0);
	for ( core::Size i=1, imax=residues.size(); i<=imax; ++i ) {
		if ( !placed[i] ) {
			if ( first==0 ) {
				placed[i] = true;
				first = residues[i];
			} else {
				if ( counter == index[level] ) {
					second=residues[i];
					placed[i]=true;
					break;
				}
				++counter;
			}
		}
	}

	cur_permutation.push_back( std::pair<core::Size, core::Size>( first, second ) ); //Add the pair to the current permutation

	//Check whether all have been placed, and recursively call this function again if not:
	bool all_placed(true);
	for ( core::Size i=1, imax=residues.size(); i<=imax; ++i ) {
		if ( !placed[i] ) {
			all_placed = false;
			break;
		}
	}
	if ( !all_placed ) {
		recursively_generate_permutation( cur_permutation, residues, placed, index, level+1 );
	}

	return;
}

/// @brief Increment the perturbation index.
/// @details Recursively calls itself.
void TryDisulfPermutations::increment_index(
	utility::vector1 < core::Size > &index,
	core::Size const index_index,
	core::Size const max_levels
) const {
	++index[index_index];
	if ( index[index_index] >= (max_levels - index_index + 1) * 2 - 1 ) {
		index[index_index] = 0;
		if ( index_index > 1 ) {
			increment_index( index, index_index - 1, max_levels);
		} else {
			for ( core::Size i=1; i<=max_levels; ++i ) index[i] = 0;
		}
	}
	return;
}

/// @brief Given a set of disulfides to form, make these disulfide bonds.
///
void TryDisulfPermutations::generate_disulfides(
	core::pose::PoseOP pose,
	utility::vector1 < std::pair <core::Size, core::Size> > disulf_pairs
) const {
	for ( core::Size i=1, imax=disulf_pairs.size(); i<=imax; ++i ) {
		core::conformation::form_disulfide( pose->conformation(), disulf_pairs[i].first, disulf_pairs[i].second, true, false );
	}
	return;
}

/// @brief Repack and minimize a set of disulfides, using only the fa_dslf and fa_dun score terms.
/// @details Returns the fa_dslf energy after repacking and minimization.
core::Real TryDisulfPermutations::repack_minimize_disulfides(
	core::pose::PoseOP pose,
	utility::vector1 < core::Size > const &disulf_res
) const {
	//A disulfide plus Dunbrack-only scorefunction:
	core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction() );
	sfxn->set_weight( core::scoring::dslf_fa13, 1.25 );
	sfxn->set_weight( core::scoring::fa_dun, 0.7 );

	//PackRotamers mover:
	protocols::minimization_packing::PackRotamersMoverOP packrot( new protocols::minimization_packing::PackRotamersMover( sfxn ) );
	//Set up TaskOperations:
	core::pack::task::TaskFactoryOP taskfact( new core::pack::task::TaskFactory() );
	taskfact->push_back(core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking() )); //Prevent design
	core::pack::task::operation::PreventRepackingOP turn_off_packing( new core::pack::task::operation::PreventRepacking() );
	for ( core::Size i = 1; i <= pose->size(); ++i ) {
		if ( !is_in_list(i, disulf_res) ) {
			turn_off_packing->include_residue(i); //Prevent residues that aren't disulfide-bondable from repacking.
		}
	}
	taskfact->push_back(turn_off_packing);
	packrot->task_factory(taskfact);

	//Apply the PackRotamers mover:
	packrot->apply(*pose);

	//MinMover:
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb(false);
	movemap->set_chi(true);
	movemap->set_jump(false);
	for ( core::Size i=1, imax=pose->size(); i<=imax; ++i ) {
		movemap->set_chi(i, is_in_list(i, disulf_res) ); //Set disulfide-forming residues as minimizable, all else as not.
	}
	protocols::minimization_packing::MinMoverOP minmover( new protocols::minimization_packing::MinMover( movemap, sfxn, mintype(), mintolerance(), true ) );
	minmover->apply(*pose); //Minimize the pose.

	//Score the pose with the disulfide term only:
	sfxn->set_weight( core::scoring::fa_dun, 0.0);
	sfxn->set_weight( core::scoring::dslf_fa13, 1.0);
	(*sfxn)(*pose);
	return pose->energies().total_energy();
}

std::string TryDisulfPermutations::get_name() const {
	return mover_name();
}

std::string TryDisulfPermutations::mover_name() {
	return "TryDisulfPermutations";
}

void TryDisulfPermutations::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "consider_already_bonded", xsct_rosetta_bool, "Also sample permutations of disulfides that are already bonded" )
		+ XMLSchemaAttribute( "min_type", xsct_minimizer_type, "Type of minimization to perform" )
		+ XMLSchemaAttribute( "min_tolerance", xsct_real, "tolerance for the minimizer" )
		+ XMLSchemaAttribute( "selector", xs_string, "Residue selector specifying which residues to consider" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Mover that tries different sets of possible disulfide bonds between existing residues capable of forming disulfides", attlist );
}

std::string TryDisulfPermutationsCreator::keyname() const {
	return TryDisulfPermutations::mover_name();
}

protocols::moves::MoverOP
TryDisulfPermutationsCreator::create_mover() const {
	return protocols::moves::MoverOP( new TryDisulfPermutations );
}

void TryDisulfPermutationsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TryDisulfPermutations::provide_xml_schema( xsd );
}


} //namespace cyclic_peptide
} //namespace protocols
