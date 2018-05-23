// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/CrosslinkerMover.cc
/// @brief This mover links two or more residues with a (possibly symmetric) cross-linker.  It adds the crosslinker, sets up constraints,
/// optionally packs and energy-mimizes it into place (packing/minimizing only the crosslinker and the side-chains to which it connects),
/// and then optionally relaxes the whole structure.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// Unit headers
#include <protocols/cyclic_peptide/CrosslinkerMover.hh>
#include <protocols/cyclic_peptide/CrosslinkerMoverCreator.hh>
#include <protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.hh>
#include <protocols/cyclic_peptide/crosslinker/TBMB_Helper.hh>
#include <protocols/cyclic_peptide/crosslinker/TMA_Helper.hh>
#include <protocols/cyclic_peptide/crosslinker/TetrahedralMetal_Helper.hh>
#include <protocols/cyclic_peptide/crosslinker/OctahedralMetal_Helper.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

// Protocols headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/mover_schemas.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "protocols.cyclic_peptide.CrosslinkerMover" );

namespace protocols {
namespace cyclic_peptide {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
CrosslinkerMover::CrosslinkerMover():
	protocols::moves::Mover( CrosslinkerMover::class_name() ),
	residue_selector_(), //Defaults to NULL pointer.
	linker_(no_crosslinker), //Defaults to no_crosslinker
	add_linker_(true),
	constrain_linker_(true),
	pack_and_minimize_linker_and_sidechains_(true),
	do_final_fastrelax_(false),
	sfxn_(),
	sidechain_frlx_rounds_(3),
	final_frlx_rounds_(3),
	filter_by_sidechain_distance_(true),
	filter_by_constraints_energy_(true),
	filter_by_total_score_(false),
	filter_by_total_score_cutoff_energy_(0.0),
	sidechain_distance_filter_multiplier_(1.0),
	constraints_energy_filter_multiplier_(1.0),
	symm_type_('A'),
	symm_count_(1),
	metal_type_("Zn")
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
////////////////////////////////////////////////////////////////////////////////
CrosslinkerMover::CrosslinkerMover( CrosslinkerMover const & /*src*/ ) = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer
/// members)
////////////////////////////////////////////////////////////////////////////////
CrosslinkerMover::~CrosslinkerMover()= default;

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
CrosslinkerMover::apply( core::pose::Pose& pose){
	//Check for a residue selector, then apply it to the pose:
	runtime_assert_string_msg( residue_selector(), "Error in protocols::cyclic_peptide::CrosslinkerMover::apply():  A residue selector must be specified before calling this function." );
	runtime_assert_string_msg( scorefxn(), "Error in protocols::cyclic_peptide::CrosslinkerMover::apply():  A scorefunction must be specified before calling this function." );

	core::select::residue_selector::ResidueSubset const selection( residue_selector()->apply(pose) );

	//Create the helper, which has the functions that set up specific types of crosslinkers:
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperOP helper;
	switch( linker_ ) {
	case TBMB :
		helper = protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperOP( protocols::cyclic_peptide::crosslinker::TBMB_HelperOP( new protocols::cyclic_peptide::crosslinker::TBMB_Helper ) );
		break;
	case TMA :
		helper = protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperOP( protocols::cyclic_peptide::crosslinker::TMA_HelperOP( new protocols::cyclic_peptide::crosslinker::TMA_Helper ) );
		break;
	case tetrahedral_metal :
		helper = protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperOP( protocols::cyclic_peptide::crosslinker::TetrahedralMetal_HelperOP( new protocols::cyclic_peptide::crosslinker::TetrahedralMetal_Helper( metal_type() ) ) );
		break;
	case octahedral_metal :
		helper = protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperOP( protocols::cyclic_peptide::crosslinker::OctahedralMetal_HelperOP( new protocols::cyclic_peptide::crosslinker::OctahedralMetal_Helper( metal_type() ) ) );
		break;

	default :
		utility_exit_with_message( "Error in protocols::cyclic_peptide::CrosslinkerMover::apply(): Invalid crosslinker specified." );
	}

	helper->set_symmetry( symm_type(), symm_count() );

	// Decide whether to call symmetric or asymmetric functions from here on:
	if ( core::conformation::symmetry::is_symmetric( pose.conformation() ) ) {
		symmetric_apply( pose, selection, helper );
	} else {
		asymmetric_apply( pose, selection, helper );
	}
}

/// @brief Given a CrossLinker enum, get its name.
///
std::string
CrosslinkerMover::get_crosslinker_name(
	CrossLinker const crosslinker
) {
	switch( crosslinker) {
	case no_crosslinker :
		return "no_crosslinker";
	case TBMB :
		return "TBMB";
	case TMA :
		return "TMA";
	case tetrahedral_metal :
		return "tetrahedral_metal";
	case octahedral_metal :
		return "octahedral_metal";
	default :
		break;
	}
	return "unknown_crosslinker";
}

/// @brief Given a CrossLinker name, get its enum.
///
CrossLinker
CrosslinkerMover::get_crosslinker_enum(
	std::string const &name
) {
	for ( core::Size i(2); i < static_cast<core::Size>(end_of_crosslinker_list); ++i ) {
		if ( !name.compare( get_crosslinker_name( static_cast<CrossLinker>(i) ) ) ) {
			return static_cast<CrossLinker>(i);
		}
	}
	return unknown_crosslinker;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
////////////////////////////////////////////////////////////////////////////////
void
CrosslinkerMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

//////////////////////////////
/// RosettaScripts Support ///
//////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
CrosslinkerMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	runtime_assert_string_msg( tag->hasOption("residue_selector"), "Error in protocols::cyclic_peptide::CrosslinkerMover::parse_my_tag(): A residue selector MUST be supplied with the \"residue_selector\" option." );
	set_residue_selector( protocols::rosetta_scripts::parse_residue_selector( tag, datamap ) );
	runtime_assert_string_msg( residue_selector(), "Error in protocols::cyclic_peptide::CrosslinkerMover::parse_my_tag(): Setting residue selector failed." );

	runtime_assert_string_msg( tag->hasOption("linker_name"), "Error in protocols::cyclic_peptide::CrosslinkerMover::parse_my_tag(): The name of the linker residue MUST be supplied with the \"linker_name\" option." );
	set_linker_name( tag->getOption<std::string>("linker_name") );

	set_symmetry( tag->getOption<std::string>("symmetry", "A1") );

	set_metal_type( tag->getOption<std::string>("metal_type", "Zn") );


	set_behaviour(
		tag->getOption<bool>( "add_linker", add_linker() ),
		tag->getOption<bool>( "constrain_linker", constrain_linker() ),
		tag->getOption<bool>( "pack_and_minimize_linker_and_sidechains", pack_and_minimize_linker_and_sidechains() ),
		tag->getOption<bool>( "do_final_fastrelax", do_final_fastrelax() )
	);

	set_filter_behaviour(
		tag->getOption<bool>( "filter_by_sidechain_distance", filter_by_sidechain_distance() ),
		tag->getOption<bool>( "filter_by_constraints_energy", filter_by_constraints_energy() ),
		tag->getOption<bool>( "filter_by_final_energy", filter_by_total_score() ),
		tag->getOption<core::Real>( "final_energy_cutoff", filter_by_total_score_cutoff_energy() ),
		tag->getOption<core::Real>( "sidechain_distance_filter_multiplier", sidechain_distance_filter_multiplier() ),
		tag->getOption<core::Real>( "constraints_energy_filter_multiplier", constraints_energy_filter_multiplier() )
	);

	runtime_assert_string_msg( tag->hasOption("scorefxn"), "Error in protocols::cyclic_peptide::CrosslinkerMover::parse_my_tag(): A scorefunction must be supplied with the \"scorefxn\" option." );
	set_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );

	set_sidechain_frlx_rounds( tag->getOption<core::Size>("sidechain_fastrelax_rounds", sidechain_frlx_rounds()) );
	set_final_frlx_rounds( tag->getOption<core::Size>("final_fastrelax_rounds", final_frlx_rounds()) );
}

/// @brief Provide information on what options are available in XML tag.
///
void
CrosslinkerMover::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	using namespace utility::tag;

	XMLSchemaRestriction linker_names_allowed;
	std::string const linker_possibles("TBMB|TMA|tetrahedral_metal|octahedral_metal");
	linker_names_allowed.name("linker_names_allowed");
	linker_names_allowed.base_type( xs_string );
	linker_names_allowed.add_restriction( xsr_pattern, linker_possibles + "(," + linker_possibles + ")+" );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "name", xs_string, "A unique name for this instance of the CrosslinkerMover." )
		+ XMLSchemaAttribute::required_attribute( "linker_name", xs_string, "The name of the type of linker to use.  For example, use TBMB for 1,3,5-tris(bromomethyl)benzene. (Allowed options are " + linker_possibles + ".)" )
		+ XMLSchemaAttribute( "symmetry", xs_string, "The symmetry of the input pose.  For example, \"C3\" for cyclic, threefold symmetry.  The symmetry must be a character followed by an integer.  Allowed characters are 'A' (asymmetric), 'C' (cyclic), 'D' (dihedral), and 'S' (mirror cyclic)." )
		+ XMLSchemaAttribute( "metal_type", xs_string, "For crosslinks mediated by metals, which metal is mediating the crosslink?  Defaults to \"Zn\".  (Should be written as an element name -- e.g. \"Cu\", \"Ca\", \"Fe\", etc.)" )
		+ XMLSchemaAttribute( "add_linker", xsct_rosetta_bool, "Should the linker geometry be added to the pose?  Default true." )
		+ XMLSchemaAttribute( "constrain_linker", xsct_rosetta_bool, "Should constraints for the linker be added to the pose?  Default true." )
		+ XMLSchemaAttribute( "pack_and_minimize_linker_and_sidechains", xsct_rosetta_bool, "Should the linker and the connecting sidechains be repacked, and should the jump to the linker, and the linker and connnecting side-chain degrees of torsional freedom, be energy-minimized?  Default true." )
		+ XMLSchemaAttribute( "sidechain_fastrelax_rounds", xs_integer, "The number of rounds of FastRelax to apply when packing and minimizing side-chains and the liker.  Default 3." )
		+ XMLSchemaAttribute( "do_final_fastrelax", xsct_rosetta_bool, "Should the whole pose be subjected to a FastRelax?  Default false." )
		+ XMLSchemaAttribute( "final_fastrelax_rounds", xs_integer, "The number of rounds of FastRelax to apply when relaxing the whole pose.  Default 3." )
		+ XMLSchemaAttribute( "filter_by_sidechain_distance", xsct_rosetta_bool, "Prior to adding the linker geometry, should this mover abort with failure status if the selected side-chains are too far apart to connect to the linker?  Default true." )
		+ XMLSchemaAttribute( "sidechain_distance_filter_multiplier", xsct_real, "This is a multiplier for the sidechain distance cutoff filter.  Higher values make the filter less stringent.  Default 1.0." )
		+ XMLSchemaAttribute( "filter_by_constraints_energy", xsct_rosetta_bool, "After adding the linker geometry, adding constraints, and repacking and minimizing the linker and the connecting side-chains, should ths mover abort with failure status if the constraints energy is too high (i.e. the energy-minimized linker geometry is bad)?  Default true." )
		+ XMLSchemaAttribute( "constraints_energy_filter_multiplier", xsct_real, "This is a multiplier for the constraints energy cutoff filter.  Higher values make the filter less stringent.  Default 1.0." )
		+ XMLSchemaAttribute( "filter_by_final_energy", xsct_rosetta_bool, "At the end of this protocol, should this mover exit with error status if the final energy is above a user-defined cutoff?  Default false." )
		+ XMLSchemaAttribute( "final_energy_cutoff", xsct_real, "If we are exiting with error status if the final energy is too high, this is the energy cutoff.  Default 0.0." )
		;
	//get attributes for parse residue selector
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "residue_selector", "A previously-defined residue selector that has been set up to select the residues that will be cross-linked." );
	protocols::rosetta_scripts::attributes_for_parse_score_function_w_description_when_required(attlist, "scorefxn", "A scorefunction to use for packing, energy-minimization, and filtering.  If constraints are turned off in this score function, they will be turned on automatically at apply time." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Adds a crosslinker linking two or more user-specified side-chains.", attlist );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
////////////////////////////////////////////////////////////////////////////////
moves::MoverOP
CrosslinkerMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new CrosslinkerMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
CrosslinkerMover::clone() const
{
	return protocols::moves::MoverOP( new CrosslinkerMover( *this ) );
}

/// @brief Get the name of the Mover
std::string
CrosslinkerMover::get_name() const
{
	return CrosslinkerMover::class_name();
}

std::string
CrosslinkerMover::class_name()
{
	return "CrosslinkerMover";
}

/// @brief Returns the name of this Mover.
std::string
CrosslinkerMover::mover_name() {
	return "CrosslinkerMover";
}


/// @brief Set the residue selector to use.
///
void
CrosslinkerMover::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	runtime_assert_string_msg( selector_in, "Error in protocols::cyclic_peptide::CrosslinkerMover::set_residue_selector(): A null pointer was passed to this function." );
	residue_selector_ = selector_in;
}

/// @brief Set the linker name.
///
void
CrosslinkerMover::set_linker_name(
	std::string const &name_in
) {
	runtime_assert_string_msg( !name_in.empty(), "Error in protocols::cyclic_peptide::CrosslinkerMover::set_linker_name(): An empty string was passed to this function." );
	CrossLinker linker_in( get_crosslinker_enum( name_in ) );
	runtime_assert_string_msg( linker_in != unknown_crosslinker, "Error in protocols::cyclic_peptide::CrosslinkerMover::set_linker_name(): Could not interpret the linker name \"" + name_in + "\"." );
	linker_ = linker_in;
}

/// @brief Get the linker name.
///
std::string
CrosslinkerMover::linker_name() const {
	return get_crosslinker_name( linker_ );
}

/// @brief Set the behaviour of this mover.
///
void
CrosslinkerMover::set_behaviour(
	bool const add_linker,
	bool const constrain_linker,
	bool const pack_and_minimize_linker_and_sidechains,
	bool const do_final_fastrelax
) {
	add_linker_ = add_linker;
	constrain_linker_ = constrain_linker;
	pack_and_minimize_linker_and_sidechains_ = pack_and_minimize_linker_and_sidechains;
	do_final_fastrelax_ = do_final_fastrelax;
}

/// @brief Set the filtering behaviour of this mover.
///
void
CrosslinkerMover::set_filter_behaviour(
	bool const filter_by_sidechain_distance,
	bool const filter_by_constraints_energy,
	bool const filter_by_total_score,
	core::Real const &filter_by_total_score_cutoff_energy,
	core::Real const &sidechain_distance_filter_multiplier,
	core::Real const &constraints_energy_filter_multiplier
) {
	filter_by_sidechain_distance_ = filter_by_sidechain_distance;
	filter_by_constraints_energy_ = filter_by_constraints_energy;
	filter_by_total_score_ = filter_by_total_score;
	filter_by_total_score_cutoff_energy_ = filter_by_total_score_cutoff_energy;
	sidechain_distance_filter_multiplier_ = sidechain_distance_filter_multiplier;
	constraints_energy_filter_multiplier_ = constraints_energy_filter_multiplier;
}

/// @brief Set the scorefunction to use for packing and minimization.
/// @details Cloned at apply time.  (That is, the scorefunction is shared until apply time).
void
CrosslinkerMover::set_scorefxn(
	core::scoring::ScoreFunctionCOP sfxn_in
) {
	runtime_assert_string_msg( sfxn_in, "Error in Error in protocols::cyclic_peptide::CrosslinkerMover::set_scorefxn(): A null pointer was passed to this function." );
	sfxn_ = sfxn_in;
}

/// @brief Get the scorefunction to use for packing and minimization.
///
core::scoring::ScoreFunctionCOP
CrosslinkerMover::scorefxn() const {
	return sfxn_;
}

/// @brief Set the number of rounds of FastRelax to apply when minimizing the linker and the
/// side-chains that connect to it.
void
CrosslinkerMover::set_sidechain_frlx_rounds(
	core::Size const rounds_in
) {
	runtime_assert_string_msg(rounds_in > 0, "Error in protocols::cyclic_peptide::CrosslinkerMover::set_sidechain_frlx_rounds(): The number of rounds must be greater than zero.");
	sidechain_frlx_rounds_ = rounds_in;
}

/// @brief Set the number of rounds of FastRelax to apply at the end.
///
void
CrosslinkerMover::set_final_frlx_rounds(
	core::Size const rounds_in
) {
	runtime_assert_string_msg(rounds_in > 0, "Error in protocols::cyclic_peptide::CrosslinkerMover::set_final_frlx_rounds(): The number of rounds must be greater than zero.");
	final_frlx_rounds_ = rounds_in;
}

/// @brief Parse a string with a symmetry type (e.g. "C3") and set the symmetry accordingly.
///
void
CrosslinkerMover::set_symmetry(
	std::string const &symmetry_in
) {
	runtime_assert( symmetry_in.length() > 1 );
	set_symm_type( symmetry_in.c_str()[0] );
	set_symm_count( static_cast< int >( std::atoi( &(symmetry_in.c_str()[1]) ) ) );
}

/// @brief Set the symmety type.
/// @details 'C' for cylic, 'S' for mirror cyclic, 'D' for dihedral, 'A' for asymmetric.
/// @note 'A' (asymmetric) by default.
void
CrosslinkerMover::set_symm_type(
	char const type_in
) {
	runtime_assert_string_msg( type_in == 'A' || type_in == 'C' || type_in == 'S' || type_in == 'D',
		"Error in protocols::cyclic_peptide::CrosslinkerMover::set_symm_type(): The type must be one of 'A' (aymmetric), 'C' (cyclic), 'D' (dihedral), or 'S' (mirror cyclic).");
	symm_type_ = type_in;
}

/// @brief Set the symmetry copy count.
/// @details For example, symm_type_='C' and symm_count_=3 would
/// specify C3 symmetry.  A value of 1 means asymmetry.  1 by default.
/// @note Deliberately an int.
void
CrosslinkerMover::set_symm_count(
	signed int const count_in
) {
	if ( symm_type() == 'A' ) {
		runtime_assert_string_msg( count_in == 1, " Error in protocols::cyclic_peptide::CrosslinkerMover::set_symm_count(): The symmetry count must be 1 for the asymmetric ('A'-type) case." );
	} else {
		runtime_assert_string_msg( count_in > 1,
			"Error in protocols::cyclic_peptide::CrosslinkerMover::set_symm_count(): The symmetry count must be greater than 1 for non-asymmetric types ('C', 'D', or 'S')." );
		if ( symm_type() == 'S' ) { runtime_assert_string_msg( count_in % 2 == 0, "Error in protocols::cyclic_peptide::CrosslinkerMover::set_symm_count(): With 'S'-type symmetry, the symmetry count must be divisible by two." ); }
	}
	symm_count_ = static_cast< core::Size >( count_in );
}

/// @brief For metal-mediated crosslinkers, set what metal mediates the crosslink.
void
CrosslinkerMover::set_metal_type(
	std::string const &metal_in
) {
	metal_type_ = metal_in;
}

std::ostream &
operator<<( std::ostream & os, CrosslinkerMover const & mover )
{
	mover.show(os);
	return os;
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
CrosslinkerMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new CrosslinkerMover );
}

std::string
CrosslinkerMoverCreator::keyname() const
{
	return CrosslinkerMover::class_name();
}

void
CrosslinkerMoverCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	CrosslinkerMover::provide_xml_schema( xsd );
}


///////////////////////
/// private methods ///
///////////////////////

/// @brief Apply the mover to a symmetric pose.
/// @details Requires symmetry in the pose matching the expected symmetry.
void
CrosslinkerMover::symmetric_apply(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperCOP helper
) {
	runtime_assert_string_msg( symm_type() != 'A' && symm_count() > 1, "Error in protocols::cyclic_peptide::CrosslinkerMover::symmetric_apply(): The function was called, but the mover has not been set for symmetry.");
	runtime_assert_string_msg( helper->selection_is_symmetric( selection, pose, helper->symm_subunits_expected() ), "Error in protocols::cyclic_peptide::CrosslinkerMover::symmetric_apply(): Either the pose lacks symmetry matching the expected symmetry, or the selector does not select equivalent, symmetric residues." );

	core::pose::PoseOP pose_copy( pose.clone() );

	bool failed(false);
	if ( filter_by_sidechain_distance() ) {
		failed=filter_by_sidechain_distance_symmetric( *pose_copy, selection, helper );
		TR << "Symmetric sidechain distance filter " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
	}
	if ( !failed ) {
		if ( add_linker() ) {
			add_linker_symmetric( *pose_copy, selection, helper );
		}
		if ( constrain_linker() ) {
			add_linker_constraints_symmetric( *pose_copy, selection, helper, add_linker() );
		}
		if ( pack_and_minimize_linker_and_sidechains() ) {
			pack_and_minimize_linker_and_sidechains( *pose_copy, selection, helper, false, true);
		}
		if ( filter_by_constraints_energy() ) {
			failed = filter_by_constraints_energy_symmetric( *pose_copy, selection, helper, add_linker() );
			TR << "Symmetric linker constraints filter " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
		}
		if ( !failed ) {
			if ( do_final_fastrelax() ) {
				pack_and_minimize_linker_and_sidechains( *pose_copy, selection, helper, true, true);
				if ( filter_by_constraints_energy() ) {
					failed = filter_by_constraints_energy_symmetric( *pose_copy, selection, helper, add_linker() );
					TR << "Symmetric linker constraints filter (after full relaxation) " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
				}
			}
			if ( !failed ) {
				if ( filter_by_total_score() ) {
					failed = filter_by_total_score( *pose_copy );
					TR << "Symmetric total score + constraints filter (after full relaxation) " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
				}
			}
		}
	}

	if ( !failed ) {
		TR << "Symmetric CrosslinkerMover reports SUCCESS.  Updating pose." << std::endl;
		pose = *pose_copy;
		set_last_move_status( protocols::moves::MS_SUCCESS );
	} else {
		TR << "Symmetric CrosslinkerMover reports FAILURE.  Returning input pose." << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
	}

}

/// @brief Determine whether the residues to be crosslinked are too far apart.  This version is for symmetric poses.
/// @details Returns TRUE for failure (too far apart), FALSE for success.
bool
CrosslinkerMover::filter_by_sidechain_distance_symmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperCOP helper
) const {
	return helper->filter_by_sidechain_distance_symmetric( pose, selection, sidechain_distance_filter_multiplier() );
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
bool
CrosslinkerMover::filter_by_constraints_energy_symmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const &selection,
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperCOP helper,
	bool const linker_was_added
) const {
	return helper->filter_by_constraints_energy_symmetric( pose, selection, linker_was_added, constraints_energy_filter_multiplier() );
}

/// @brief Given a selection of residues, add a crosslinker, align it crudely to the
/// selected residues, and set up covalent bonds.  This version is for symmetric poses.
void
CrosslinkerMover::add_linker_symmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperCOP helper
) const {
	helper->add_linker_symmetric(pose, selection);
}

/// @brief Given a selection of residues that have already been connected to a crosslinker,
/// add constraints for the crosslinker.  This version is for symmetric poses.
void
CrosslinkerMover::add_linker_constraints_symmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperCOP helper,
	bool const linker_was_added
) const {
	helper->add_linker_constraints_symmetric( pose, selection, linker_was_added );
}

/// @brief Apply the mover to an asymmetric pose.
/// @details Requires and asymmetric pose, and no symmetry.
void
CrosslinkerMover::asymmetric_apply(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperCOP helper
) {
	runtime_assert_string_msg( symm_type() == 'A' && symm_count() == 1,
		"Error in protocols::cyclic_peptide::CrosslinkerMover::asymmetric_apply(): The function was called, but the mover is set up for symmetry.");

	core::pose::PoseOP pose_copy( pose.clone() );

	bool failed(false);

	if ( filter_by_sidechain_distance() ) {
		failed = filter_by_sidechain_distance_asymmetric( *pose_copy, selection, helper );
		TR << "Sidechain distance filter " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
	}
	if ( !failed ) {
		if ( add_linker() ) {
			add_linker_asymmetric( *pose_copy, selection, helper );
		}
		if ( constrain_linker() ) {
			add_linker_constraints_asymmetric( *pose_copy, selection, helper );
		}
		if ( pack_and_minimize_linker_and_sidechains() ) {
			pack_and_minimize_linker_and_sidechains( *pose_copy, selection, helper, false, false);
		}
		if ( filter_by_constraints_energy() ) {
			failed = filter_by_constraints_energy_asymmetric( *pose_copy, selection, helper );
			TR << "Linker constraints filter " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
		}
		if ( !failed ) {
			if ( do_final_fastrelax() ) {
				pack_and_minimize_linker_and_sidechains( *pose_copy, selection, helper, true, false);
				if ( filter_by_constraints_energy() ) {
					failed = filter_by_constraints_energy_asymmetric( *pose_copy, selection, helper );
					TR << "Linker constraints filter (after full relaxation) " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
				}
			}
			if ( !failed ) {
				if ( filter_by_total_score() ) {
					failed = filter_by_total_score( *pose_copy );
					TR << "Total score + constraints filter (after full relaxation) " << (failed ? "FAILED" : "PASSED" ) << "." << std::endl;
				}
			}
		}
	}

	if ( !failed ) {
		TR << "CrosslinkerMover reports SUCCESS.  Updating pose." << std::endl;
		pose = *pose_copy;
		set_last_move_status( protocols::moves::MS_SUCCESS );
	} else {
		TR << "CrosslinkerMover reports FAILURE.  Returning input pose." << std::endl;
		set_last_move_status( protocols::moves::FAIL_RETRY );
	}
}

/// @brief Determine whether the residues to be crosslinked are too far apart.
/// @details Returns TRUE for failure (too far apart), FALSE for success.
bool
CrosslinkerMover::filter_by_sidechain_distance_asymmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperCOP helper
) const {
	return helper->filter_by_sidechain_distance_asymmetric( pose, selection, sidechain_distance_filter_multiplier() );
}

/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
bool
CrosslinkerMover::filter_by_constraints_energy_asymmetric(
	core::pose::Pose const &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperCOP helper
) const {
	return helper->filter_by_constraints_energy_asymmetric( pose, selection, constraints_energy_filter_multiplier() );
}

/// @brief Determine whether the overall system has too high an overall score (including constraints) at the end of the protocol.
/// @details Returns TRUE for failure (too high an overall score) and FALSE for success.
bool
CrosslinkerMover::filter_by_total_score(
	core::pose::Pose const &pose
) const {
	core::pose::Pose pose_copy(pose);

	core::scoring::ScoreFunctionOP sfxn( scorefxn()->clone() ); //Local copy of scorefunction.
	if ( sfxn->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0); }
	if ( sfxn->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::angle_constraint, 1.0); }
	if ( sfxn->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::dihedral_constraint, 1.0); }

	(*sfxn)(pose_copy);
	core::Real const total_eng( pose_copy.energies().total_energy() );
	bool const failed( total_eng > filter_by_total_score_cutoff_energy() );

	if ( TR.Debug.visible() ) {
		TR.Debug << "Total energy (including constraints) at end of CrosslinkerMover protocol: " << total_eng << std::endl;
		TR.Debug << "Energy cutoff for CrosslinkerMover protocol: " << filter_by_total_score_cutoff_energy() << std::endl;
		if ( failed ) { TR.Debug << "Energy filter FAILED." << std::endl; }
		else { TR.Debug << "Energy filter PASSED." << std::endl; }
	}

	return failed;
}

/// @brief Given a selection of residues, add a crosslinker, align it crudely to the
/// selected residues, and set up covalent bonds.
void
CrosslinkerMover::add_linker_asymmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperCOP helper
) const {
	helper->add_linker_asymmetric(pose, selection);
}

/// @brief Given a selection of residues that have already been connected to a crosslinker,
/// add constraints for the crosslinker.
void
CrosslinkerMover::add_linker_constraints_asymmetric(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperCOP helper
) const {
	helper->add_linker_constraints_asymmetric( pose, selection );
}

/// @brief Repack and minimize the sidechains.
/// @details Also repacks and minimzes the linker, letting all jumps vary.
void
CrosslinkerMover::pack_and_minimize_linker_and_sidechains(
	core::pose::Pose &pose,
	core::select::residue_selector::ResidueSubset const & selection,
	protocols::cyclic_peptide::crosslinker::CrosslinkerMoverHelperCOP helper,
	bool const whole_structure,
	bool const symmetric
) const {

	core::scoring::ScoreFunctionOP sfxn( scorefxn()->clone() ); //Local copy of scorefunction.
	if ( sfxn->get_weight( core::scoring::atom_pair_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::atom_pair_constraint, 1.0); }
	if ( sfxn->get_weight( core::scoring::angle_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::angle_constraint, 1.0); }
	if ( sfxn->get_weight( core::scoring::dihedral_constraint ) == 0.0 ) { sfxn->set_weight( core::scoring::dihedral_constraint, 1.0); }

	core::Size const rounds( whole_structure ? final_frlx_rounds() : sidechain_frlx_rounds() );

	for ( core::Size i=1; i<=rounds; ++i ) {
		helper->pre_relax_round_update_steps(pose, selection, whole_structure, symmetric, add_linker());
		protocols::relax::FastRelax frlx( sfxn, 1 );

		if ( !whole_structure ) {
			utility::vector1< core::Size > res_indices, linker_indices;

			helper->get_sidechain_indices( selection, res_indices );
			if ( symmetric ) {
				if ( add_linker() && helper->helper_adds_linker_residue() ) {
					for ( core::Size i(2), imax(res_indices.size()); i<=imax; ++i ) {
						res_indices[i] += (i-1); // This may need to be made a virtual function, overridden by derived classes, at some point.
					}
				}
				helper->get_linker_indices_symmetric(pose, res_indices, linker_indices);
			} else {
				if ( helper->adds_crosslinker_residue() ) {
					linker_indices.resize(1);
					linker_indices[1] = helper->get_linker_index_asymmetric( pose, res_indices );
				}
			}

			core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
			movemap->set_bb(false);
			movemap->set_chi(false);
			movemap->set_jump(false);
			for ( core::Size i(1), imax(res_indices.size()); i<=imax; ++i ) {
				movemap->set_chi( res_indices[i], true );
			}

			for ( core::Size i(1), imax(linker_indices.size()); i<=imax; ++i ) {
				movemap->set_chi(linker_indices[i], true);
			}
			if ( symmetric ) {
				if ( helper->adds_crosslinker_residue() ) {
					utility::vector1< core::Size > jump_indices;
					get_jump_indices_for_symmetric_crosslinker( pose, linker_indices, jump_indices );
					for ( core::Size i(1), imax(jump_indices.size()); i<=imax; ++i ) {
						movemap->set_jump(jump_indices[i], true);
					}
				}
			} else {
				if ( helper->adds_crosslinker_residue() ) {
					movemap->set_jump( get_jump_index_for_crosslinker( pose, linker_indices[1] ), true );
				}
			}

			frlx.set_movemap(movemap);
		}

		frlx.apply(pose);
		helper->post_relax_round_update_steps(pose, selection, whole_structure, symmetric, add_linker());
	}
}

/// @brief Given a pose and the index of a crosslinker, figure out the jump in the foldtree that moves the crosslinker.
/// @details Throws an error if the foldtree isn't set up so that a unique jump moves the crosslinker.
core::Size
CrosslinkerMover::get_jump_index_for_crosslinker(
	core::pose::Pose const &pose,
	core::Size const linker_index
) const {
	runtime_assert_string_msg( linker_index > 0  && linker_index <= pose.total_residue(),
		"Error in protocols::cyclic_peptide::CrosslinkerMover::get_jump_index_for_crosslinker(): The linker index is out of range." );

	core::kinematics::FoldTree const &foldtree( pose.fold_tree() ); //Get a reference to the fold tree.
	runtime_assert_string_msg( foldtree.is_jump_point( linker_index ),
		"Error in protocols::cyclic_peptide::CrosslinkerMover::get_jump_index_for_crosslinker(): The linker index is not a jump point." );
	return foldtree.get_jump_that_builds_residue( linker_index );
}

/// @brief Given a pose and the indices of the pieces of a symmetric crosslinker, figure out the jumps in the foldtree that move the crosslinker.
/// @details Throws an error if the foldtree isn't set up so that a unique jump moves the crosslinker.
void
CrosslinkerMover::get_jump_indices_for_symmetric_crosslinker(
	core::pose::Pose const &pose,
	utility::vector1< core::Size > const & linker_indices_in,
	utility::vector1< core::Size > & jump_indices_out
) const {
	core::Size const linker_count( linker_indices_in.size() );
	if ( linker_count == 0 ) { //There's no linker, so return no jump indices.
		jump_indices_out.clear();
		return;
	}
	jump_indices_out.resize(linker_count);

	core::kinematics::FoldTree const &foldtree( pose.fold_tree() ); //Get a reference to the fold tree.
	for ( core::Size i(1); i<=linker_count; ++i ) {
		runtime_assert_string_msg( foldtree.is_jump_point( linker_indices_in[i] ),
			"Error in protocols::cyclic_peptide::CrosslinkerMover::get_jump_indices_for_symmetric_crosslinker(): A provided linker index is not a jump point." );
		jump_indices_out[i] = foldtree.get_jump_that_builds_residue( linker_indices_in[i] );
	}

}

} //protocols
} //cyclic_peptide
