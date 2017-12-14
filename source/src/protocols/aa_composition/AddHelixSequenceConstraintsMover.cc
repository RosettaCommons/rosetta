// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/AddHelixSequenceConstraintsMover.cc
/// @brief This mover adds sequence constraints to the ends of each helix, requiring at least one positively-charged residue in the three C-terminal residues, and at least one negatively-charged resiude in the three N-terminal residues.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu -- wrote the mover)
/// @author Gabriel Rocklin (came up with the rules that this mover should apply, based on high-throughput experiments)

// Unit headers
#include <protocols/aa_composition/AddHelixSequenceConstraintsMover.hh>
#include <protocols/aa_composition/AddHelixSequenceConstraintsMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/aa_composition_energy/AACompositionConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>

// Protocols headers
#include <protocols/aa_composition/ClearCompositionConstraintsMover.hh>
#include <protocols/aa_composition/util.hh>
#include <protocols/rosetta_scripts/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.aa_composition.AddHelixSequenceConstraintsMover" );

namespace protocols {
namespace aa_composition {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
AddHelixSequenceConstraintsMover::AddHelixSequenceConstraintsMover():
	protocols::moves::Mover( AddHelixSequenceConstraintsMover::mover_name() ),
	reset_mode_(false),
	min_helix_length_(8),
	add_n_terminal_constraints_(true),
	add_c_terminal_constraints_(true),
	add_overall_constraints_(true),
	add_alanine_constraints_(true),
	add_hydrophobic_constraints_(true),
	min_n_terminal_charges_(2),
	n_terminus_size_(3),
	n_terminal_constraint_strength_(15.0),
	min_c_terminal_charges_(2),
	c_terminus_size_(3),
	c_terminal_constraint_strength_(15.0),
	types_to_avoid_(),
	overall_max_count_(0),
	overall_constraints_strength_(5.0),
	desired_ala_fraction_(0.1),
	ala_constraint_under_strength_(0.2),
	ala_constraint_over_strength_(0.2),
	desired_min_hydrophobic_fraction_(0.25),
	hydrophobic_constraint_strength_(0.2),
	residue_selector_()
{
	initialize_types_to_avoid();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
AddHelixSequenceConstraintsMover::AddHelixSequenceConstraintsMover( AddHelixSequenceConstraintsMover const & /*src*/ ) = default;

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
AddHelixSequenceConstraintsMover::~AddHelixSequenceConstraintsMover()= default;

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover.
/// @details Calls find_helices(), do_terminal_constraints(), do_overall_or_alanine_constraints().
void
AddHelixSequenceConstraintsMover::apply( core::pose::Pose &pose ){
	check_lengths_sensible();

	if ( reset_mode() ) delete_old_sequence_constraints( pose );

	utility::vector1 < std::pair < core::Size, core::Size > > helices;
	find_helices( pose, helices );

	std::string n_terminal_constraint_setup, c_terminal_constraint_setup, overall_constraint_setup, alanine_constraint_setup, hydrophobic_constraint_setup;
	if ( add_n_terminal_constraints() ) set_up_terminal_constraints( n_terminal_constraint_setup, true );
	if ( add_c_terminal_constraints() ) set_up_terminal_constraints( c_terminal_constraint_setup, false );
	if ( add_overall_constraints() ) set_up_overall_constraints( overall_constraint_setup );
	if ( add_alanine_constraints() ) set_up_alanine_constraints( alanine_constraint_setup );
	if ( add_hydrophobic_constraints() ) set_up_hydrophobic_constraints( hydrophobic_constraint_setup );

	if ( has_residue_selector() ) {
		remove_unselected_helices( pose, helices );
	}

	for ( core::Size i(1), imax(helices.size()); i<=imax; ++i ) {
		TR << "Applying sequence constraints for helix " << i << ", running from residue " << helices[i].first << " through residue " << helices[i].second << "." << std::endl;
		if ( add_n_terminal_constraints() ) {
			TR << "\tAdding N-terminal sequence constraints to helix." << std::endl;
			do_terminal_constraints( pose, true, helices[i].first, helices[i].second, n_terminal_constraint_setup );
		}
		if ( add_c_terminal_constraints() ) {
			TR << "\tAdding C-terminal sequence constraints to helix." << std::endl;
			do_terminal_constraints( pose, false, helices[i].first, helices[i].second, c_terminal_constraint_setup );
		}
		if ( add_overall_constraints() ) {
			TR << "\tAdding sequence constraints to helix penalizing helix-disfavouring residue types." << std::endl;
			do_overall_or_alanine_constraints( pose, helices[i].first, helices[i].second, overall_constraint_setup );
		}
		if ( add_alanine_constraints() ) {
			TR << "\tAdding sequence constraints to helix to control fractional alanine content." << std::endl;
			do_overall_or_alanine_constraints( pose, helices[i].first, helices[i].second, alanine_constraint_setup );
		}
		if ( add_hydrophobic_constraints() ) {
			TR << "\tAdding sequence constraints to helix to enforce a minimum hydrophobic content." << std::endl;
			do_overall_or_alanine_constraints( pose, helices[i].first, helices[i].second, hydrophobic_constraint_setup );
		}
		TR.flush();
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
AddHelixSequenceConstraintsMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
AddHelixSequenceConstraintsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	core::select::residue_selector::ResidueSelectorCOP res_sel( protocols::rosetta_scripts::parse_residue_selector( tag, data ) );
	if ( res_sel && res_sel != nullptr ) {
		set_residue_selector( res_sel );
	}

	set_reset_mode( tag->getOption<bool>( "reset", reset_mode() ) );
	set_min_helix_length( tag->getOption<core::Size>( "min_helix_length", min_helix_length() ) );

	set_add_n_terminal_constraints(
		tag->getOption<bool>( "add_n_terminal_constraints", add_n_terminal_constraints() ),
		tag->getOption<core::Size>( "min_n_terminal_charges", min_n_terminal_charges() ),
		tag->getOption<core::Size>( "n_terminal_residues", n_terminus_size() ),
		tag->getOption<core::Real>( "n_terminal_constraint_strength", n_terminal_constraint_strength() )
	);

	set_add_c_terminal_constraints(
		tag->getOption<bool>( "add_c_terminal_constraints", add_c_terminal_constraints() ),
		tag->getOption<core::Size>( "min_c_terminal_charges", min_c_terminal_charges() ),
		tag->getOption<core::Size>( "c_terminal_residues", c_terminus_size() ),
		tag->getOption<core::Real>( "c_terminal_constraint_strength", c_terminal_constraint_strength() )
	);

	set_add_overall_constraints(
		tag->getOption<bool>( "add_overall_constraints", add_overall_constraints() ),
		tag->getOption<std::string>( "types_to_avoid", "ASN ASP SER GLY THR VAL" ),
		tag->getOption<core::Size>( "overall_max_count", overall_max_count() ),
		tag->getOption<core::Real>( "overall_constraints_strength", overall_constraints_strength() )
	);

	set_add_alanine_constraints(
		tag->getOption<bool>( "add_alanine_constraints", add_alanine_constraints() ),
		tag->getOption<core::Real>( "desired_alanine_fraction", desired_ala_fraction() ),
		tag->getOption<core::Real>( "ala_constraint_under_strength", ala_constraint_under_strength() ),
		tag->getOption<core::Real>( "ala_constraint_over_strength", ala_constraint_over_strength() )
	);

	set_add_hydrophobic_constraints(
		tag->getOption<bool>( "add_hydrophobic_constraints", add_hydrophobic_constraints() ),
		tag->getOption<core::Real>( "desired_min_hydrophobic_fraction", desired_min_hydrophobic_fraction() ),
		tag->getOption<core::Real>( "hydrophobic_constraint_strength", hydrophobic_constraint_strength() )
	);

	check_lengths_sensible();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
AddHelixSequenceConstraintsMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new AddHelixSequenceConstraintsMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
AddHelixSequenceConstraintsMover::clone() const
{
	return protocols::moves::MoverOP( new AddHelixSequenceConstraintsMover( *this ) );
}

std::string AddHelixSequenceConstraintsMover::get_name() const {
	return mover_name();
}

std::string AddHelixSequenceConstraintsMover::mover_name() {
	return "AddHelixSequenceConstraints";
}

void AddHelixSequenceConstraintsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	core::select::residue_selector::attributes_for_parse_residue_selector(
		attlist, "residue_selector",
		"An optional, previously-defined ResidueSelector.  If provided, only helices that contain at least one residue that is selected by the residue selector will have constraints applied.  If not used, constraints are applied to all helices in the pose." );

	attlist + XMLSchemaAttribute::attribute_w_default(
		"reset", xsct_rosetta_bool,
		"If true, all sequence constraints in the pose will be cleared (deleted) before this mover is applied.  If false, the mover will append to existing sequence constraints.  False by default.",
		"false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"min_helix_length", xsct_non_negative_integer,
		"The minimum number of residues that a helix must have for this mover to act on it.  By default, helices smaller than 8 residues are ignored since they have negligible helix macrodipoles.",
		"8");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"add_n_terminal_constraints", xsct_rosetta_bool,
		"If true, this mover will add sequence constraints requiring a user-specified minimum number of negatively-charged residues at the N-terminus of each helix.  True by default.",
		"true");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"min_n_terminal_charges", xsct_non_negative_integer,
		"The minimum number of negatively-charged residues required at the N-terminus of helices.  Defaults to 2 residues.",
		"2");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"n_terminal_residues", xsct_non_negative_integer,
		"The length of the stretch of residues that must contain negative charges at the N-terminus of a helix.  Defaults to 3 residues.",
		"3");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"n_terminal_constraint_strength", xsct_real,
		"The strength of the sequence constraint requiring negative charges at the N-termini of helices.  If set to be too weak, Rosetta's packer may sometimes put in too few negative charges.  15.0 by default.\n(For advanced users, this is the energetic penalty applied when there is one fewer than the desired number of negatively-charged residues.  The penalty ramps quadratically for two fewer, three fewer, etc.)",
		"15.0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"add_c_terminal_constraints", xsct_rosetta_bool,
		"If true, this mover will add sequence constraints requiring a user-specified minimum number of positively-charged residues at the C-terminus of each helix.  True by default.",
		"true");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"min_c_terminal_charges", xsct_non_negative_integer,
		"The minimum number of positively-charged residues required at the C-terminus of helices.  Defaults to 2 residues.",
		"2");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"c_terminal_residues", xsct_non_negative_integer,
		"The length of the stretch of residues that must contain positive charges at the C-terminus of a helix.  Defaults to 3 residues.",
		"3");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"c_terminal_constraint_strength", xsct_real,
		"The strength of the sequence constraint requiring positive charges at the N-termini of helices.  If set to be too weak, Rosetta's packer may sometimes put in too few positive charges.  15.0 by default.\n(For advanced users, this is the energetic penalty applied when there is one fewer than the desired number of positively-charged residues.  The penalty ramps quadratically for two fewer, three fewer, etc.)",
		"15.0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"add_overall_constraints", xsct_rosetta_bool,
		"If true, this mover will add sequence constraints penalizing more than a user-specified maximum number of helix-disfavouring residues in each helix.  True by default.",
		"true");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"types_to_avoid", xs_string,
		"The list of helix-disfavouring residue types that should be penalized at all helix positions by the \"overall\" constraints.  This must be a whitespace-separated list of three-letter residue type codes, with no commas.",
		"ASN ASP SER GLY THR VAL");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"overall_max_count", xsct_non_negative_integer,
		"The maximum allowed number of helix-disfavouring residue types.  Default zero (though the penalty for having one is small by default).",
		"0");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"overall_constraints_strength", xsct_real,
		"The strength of the sequence constraint penalizing helix-disfavouring residue types.  If set to be too weak, Rosetta's packer may sometimes put in too few helix-disfavouring residues.  5.0 by default.\n(For advanced users, this is the energetic penalty applied when there is one more than the maximum allowed number of helix-disfavouring residues.  The penalty ramps quadratically for two more, three more, etc.)",
		"5.0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"add_alanine_constraints", xsct_rosetta_bool,
		"If true, this mover will add sequence constraints penalizing too many or too few alanine residues in each helix.  The user can set a user-defined desired fractional alanine content.  Note that this constraint is usually set to be weak, so that some deviation from the desired fractional alanine content is tolerated.  True by default.",
		"true");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"desired_alanine_fraction", xsct_real,
		"The desired fractional alanine content in each helix.  Defaults to 0.1 (10 percent alanine).  Note that deviation from this fraction is possible, if the alanine constraints are weak (and they are weak by default).",
		"0.1");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"ala_constraint_under_strength", xsct_real,
		"The alanine constraint penalty that is imposed if the fractional alanine content is 1 percent less than the desired fractional content.  The penalty ramps quadratically as the alanine content falls below desired.  Note that this penalty is weak, by default, to allow some deviation from the desired fractional alanine content on a case-by-case basis.",
		"0.2");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"ala_constraint_over_strength", xsct_real,
		"The alanine constraint penalty that is imposed if the fractional alanine content is 1 percent more than the desired fractional content.  The penalty ramps quadratically as the alanine content rises above desired.  Note that this penalty is weak, by default, to allow some deviation from the desired fractional alanine content on a case-by-case basis.",
		"0.2");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"add_hydrophobic_constraints", xsct_rosetta_bool,
		"If true, this mover will add sequence constraints penalizing too few hydrophobic residues in each helix.  The user can set a user-defined minimum desired fractional hydrophobic content.  Note that this constraint is usually set to be weak, so that some designs with fewer than the minimum hydrophobic count will be returned.  Note that alanine is NOT considered hydrophobic.  True by default.",
		"true");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"desired_min_hydrophobic_fraction", xsct_real,
		"The desired minimum fractional hydrophobic content in each helix.  Defaults to 0.25 (25 percent hydrophobic).  There is no penalty for more hydrophobic residues; only for fewer than desired.  Note that helices slightly below this fraction are possible, if the hydrophobic constraints are weak (and they are weak by default).",
		"0.25");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"hydrophobic_constraint_strength", xsct_real,
		"The hydrophobic constraint penalty that is imposed if the fractional hydrophobic content is 1 percent less than the desired fractional content.  The penalty ramps quadratically as the hydrophobic content falls below desired.  Note that this penalty is weak, by default, to allow some deviation from the desired fractional hydrophobic content on a case-by-case basis.",
		"0.2");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"The AddHelixSequenceConstraints mover sets up sequence constraints for each helix in a pose or in a selection.  "
		"It can require negative and positive charges at the N- and C-termini, respectively, can limit the number of "
		"helix-disfavouring residues in each helix, can require that the helix be a user-specified fraction alanine, "
		"and can require a minimum fractional hydrophobic content in each helix.  "
		"Note that these constraints remain attached to the pose, and are intended to be used during design with the aa_composition "
		"score term.  Helices are detected using DSSP when this mover is applied, so if the secondary structure changes between "
		"application of this mover and design, the constraints will applied to out-of-date residue indices.  (In such a case, the "
		"sequence constraints can be re-applied with this mover after first clearing the old constraints with either the "
		"ClearCompositionConstraintsMover, or by setting \"reset=true\" in this mover's options.)\n\n"
		"Note that this mover's defaults have been set so that it can be applied without manually setting anything, and still "
		"produce reasonable behaviour.  For advanced users, all settings can be tweaked manually, but this shouldn't be necessary "
		"in many cases."
		, attlist
	);

}

/***********************************************
Public setters
***********************************************/

/// @brief Set a residue selector if specified by the user
///
void AddHelixSequenceConstraintsMover::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in) { residue_selector_ = selector_in; }


/// @brief Set whether old composition constraints should be deleted before applying this mover.
///
void AddHelixSequenceConstraintsMover::set_reset_mode( bool const setting ) { reset_mode_ = setting; }

/// @brief Set the minimum number of residues in a helix for the helix to be considered.
/// @details Very short helices (e.g. one turn, five residues) have negligible helix macrodipoles,
/// and can be disregarded.  Defaults to 8.
void AddHelixSequenceConstraintsMover::set_min_helix_length( core::Size const setting ) { min_helix_length_ = setting; }

/// @brief Set whether this mover should add constraints to the N-termini of helices requiring negative charges there?
///
void
AddHelixSequenceConstraintsMover::set_add_n_terminal_constraints(
	bool const setting,
	core::Size const min_charges,
	core::Size const nterm_size,
	core::Real const &strength
) {
	add_n_terminal_constraints_ = setting;
	runtime_assert_string_msg( min_charges <= nterm_size,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::set_add_n_terminal_constraints(): The minimum number of charges cannot be greater than the number of N-terminal residues being considered." );
	min_n_terminal_charges_ = min_charges;
	n_terminus_size_ = nterm_size;
	runtime_assert_string_msg( strength >= 0,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::set_add_n_terminal_constraints(): The strength of the N-terminal charge constraint must be greater than or equal to zero." );
	n_terminal_constraint_strength_ = strength;
}

/// @brief Set whether this mover should add constraints to the C-termini of helices requiring negative charges there?
/// @details Default true.
void
AddHelixSequenceConstraintsMover::set_add_c_terminal_constraints(
	bool const setting,
	core::Size const min_charges,
	core::Size const cterm_size,
	core::Real const &strength
) {
	add_c_terminal_constraints_ = setting;
	runtime_assert_string_msg( min_charges <= cterm_size,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::set_add_c_terminal_constraints(): The minimum number of charges cannot be greater than the number of C-terminal residues being considered." );
	min_c_terminal_charges_ = min_charges;
	c_terminus_size_ = cterm_size;
	runtime_assert_string_msg( strength >= 0,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::set_add_c_terminal_constraints(): The strength of the C-terminal charge constraint must be greater than or equal to zero." );
	c_terminal_constraint_strength_ = strength;
}

/// @brief Set whether this mover should add constraints to the helix penalizing helix-disfavouring residue types?
/// @details Default true.
void AddHelixSequenceConstraintsMover::set_add_overall_constraints( bool const setting ) { add_overall_constraints_ = setting; }

/// @brief Set whether this mover should add constraints to the helix penalizing helix-disfavouring residue types.
/// @details Default true.
void AddHelixSequenceConstraintsMover::set_add_overall_constraints(
	bool const setting,
	std::string const &types_to_avoid,
	core::Size const overall_max_count,
	core::Real const &overall_constraints_strength
) {
	add_overall_constraints_ = setting;

	types_to_avoid_.clear();
	types_to_avoid_ = utility::split_whitespace( types_to_avoid );
	runtime_assert_string_msg(
		types_to_avoid.size() > 0,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::set_add_overall_constraints():  At least one residue type must be specified."
	);

	overall_max_count_ = overall_max_count;
	runtime_assert_string_msg(
		overall_constraints_strength >= 0.0,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::set_add_overall_constraints():  The overall constraint strength must be greater than or equal to zero."
	);
}


/// @brief Set whether this mover should add constraints to the helix requiring a certain alanine fraction?
/// @details Default true.
void AddHelixSequenceConstraintsMover::set_add_alanine_constraints( bool const setting ) { add_alanine_constraints_ = setting; }

/// @brief Set whether this mover should add constraints to the helix requiring a certain alanine fraction.
/// @details Default true.
void
AddHelixSequenceConstraintsMover::set_add_alanine_constraints(
	bool const setting,
	core::Real const &desired_ala_fraction,
	core::Real const &ala_constraint_under_strength,
	core::Real const &ala_constraint_over_strength
) {
	runtime_assert_string_msg(
		ala_constraint_under_strength >= 0.0,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::set_add_alanine_constraints():  The alanine-under constraint strength must be greater than or equal to zero."
	);
	runtime_assert_string_msg(
		ala_constraint_over_strength >= 0.0,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::set_add_alanine_constraints():  The alanine-over constraint strength must be greater than or equal to zero."
	);
	runtime_assert_string_msg(
		desired_ala_fraction >= 0.0 && desired_ala_fraction <= 1.0,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::set_add_alanine_constraints():  The desired per-helix alanine fraction must be between zero and one."
	);

	add_alanine_constraints_ = setting;
	desired_ala_fraction_ = desired_ala_fraction;
	ala_constraint_under_strength_ = ala_constraint_under_strength;
	ala_constraint_over_strength_ = ala_constraint_over_strength;
}

/// @brief Set whether this mover should add constraints to the helix requiring at minimum fraction of hydrophobic residues.
/// @details Default true.
void
AddHelixSequenceConstraintsMover::set_add_hydrophobic_constraints(
	bool const setting
) {
	add_hydrophobic_constraints_ = setting;
}

/// @brief Set whether this mover should add constraints to the helix requiring at minimum fraction of hydrophobic residues.
/// @details Default true.
void
AddHelixSequenceConstraintsMover::set_add_hydrophobic_constraints(
	bool const setting,
	core::Real const &desired_min_hydrophobic_fraction,
	core::Real const &hydrophobic_constraint_strength
) {
	add_hydrophobic_constraints_ = setting;

	runtime_assert_string_msg(
		desired_min_hydrophobic_fraction >= 0.0 && desired_min_hydrophobic_fraction <= 1.0,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::set_add_hydrophobic_constraints():  The desired per-helix hydrophobic fraction must be between zero and one."
	);
	desired_min_hydrophobic_fraction_ = desired_min_hydrophobic_fraction;

	runtime_assert_string_msg(
		hydrophobic_constraint_strength >= 0.0,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::set_add_hydrophobic_constraints():  The hydrophobic constraint penalty must be greater than or equal to zero."
	);
	hydrophobic_constraint_strength_ = hydrophobic_constraint_strength;
}


/***********************************************
Private methods
***********************************************/

/// @brief Once the mover has been set up, check that the lengths of termini is sensible compared to the min_helix_length setting.
/// @details Throws an error if add_n_terminal_constraints or add_c_terminal_constraints is true, and the min_helix_length is less
/// than the greater of the active termini lenghts.
void
AddHelixSequenceConstraintsMover::check_lengths_sensible() const {
	core::Size const minlen( std::max( add_n_terminal_constraints() ? n_terminus_size() : 0, add_c_terminal_constraints() ? c_terminus_size() : 0 ) );
	runtime_assert_string_msg( min_helix_length() >= minlen,
		"Error in protocols::aa_composition::AddHelixSequenceConstraintsMover::check_lengths_sensible(): If N- and/or C-terminal sequence constraints are added, then the terminus size for the active option(s) must be less than or equal to the minimum helix length."
	);
}


/// @brief Initialize the types to avoid to ASN, ASP, SER, GLY, THR, VAL.
///
void
AddHelixSequenceConstraintsMover::initialize_types_to_avoid() {
	types_to_avoid_.resize(6);
	types_to_avoid_[1] = "ASN";
	types_to_avoid_[2] = "ASP";
	types_to_avoid_[3] = "SER";
	types_to_avoid_[4] = "GLY";
	types_to_avoid_[5] = "THR";
	types_to_avoid_[6] = "VAL";
}

/// @brief Given a pose, remove all sequence constraints from it.
///
void
AddHelixSequenceConstraintsMover::delete_old_sequence_constraints(
	core::pose::Pose &pose
) const {
	ClearCompositionConstraintsMover clear_comp_csts;
	clear_comp_csts.apply( pose );
}

/// @brief Set up strings defining terminal constraints.
///
void
AddHelixSequenceConstraintsMover::set_up_terminal_constraints(
	std::string &setup_string_out,
	bool const n_terminal
) const {
	std::stringstream outstr;

	outstr << "PENALTY_DEFINITION\n";
	outstr << "PROPERTIES " << (n_terminal ? "NEGATIVE_CHARGE\n" : "POSITIVE_CHARGE\n");
	outstr << "DELTA_START -1\n";
	outstr << "DELTA_END 1\n";
	outstr << "PENALTIES " << (n_terminal ? n_terminal_constraint_strength() : c_terminal_constraint_strength() ) << " 0 0\n";
	outstr << "ABSOLUTE " << (n_terminal ? min_n_terminal_charges() : min_c_terminal_charges() ) << "\n";
	outstr << "BEFORE_FUNCTION QUADRATIC\n";
	outstr << "AFTER_FUNCTION CONSTANT\n";
	outstr << "END_PENALTY_DEFINITION\n";

	setup_string_out = outstr.str();
}

/// @brief Set up strings defining constraints penalizing helix-disfavouring residue types.
///
void
AddHelixSequenceConstraintsMover::set_up_overall_constraints(
	std::string &setup_string_out
) const {
	std::stringstream outstr;

	outstr << "PENALTY_DEFINITION\n";
	outstr << "TYPE ";
	for ( core::Size i(1), imax( types_to_avoid().size() ); i<=imax; ++i ) {
		outstr << types_to_avoid(i);
		if ( i<imax ) { outstr << " "; }
		else { outstr << "\n"; }
	}
	outstr << "DELTA_START -1\n";
	outstr << "DELTA_END 1\n";
	outstr << "PENALTIES 0 0 " << overall_constraints_strength() << "\n";
	outstr << "ABSOLUTE " << overall_max_count() << "\n";
	outstr << "BEFORE_FUNCTION CONSTANT\n";
	outstr << "AFTER_FUNCTION QUADRATIC\n";
	outstr << "END_PENALTY_DEFINITION\n";

	setup_string_out = outstr.str();
}

/// @brief Set up strings defining constraints penalizing too few or too many alanines.
///
void
AddHelixSequenceConstraintsMover::set_up_alanine_constraints(
	std::string &setup_string_out
) const {
	std::stringstream outstr;

	outstr << "PENALTY_DEFINITION\n";
	outstr << "TYPE ALA\n";
	outstr << "FRACT_DELTA_START -0.01\n";
	outstr << "FRACT_DELTA_END 0.01\n";
	outstr << "PENALTIES " << ala_constraint_under_strength() << " 0 " << ala_constraint_over_strength() << "\n";
	outstr << "FRACTION " << desired_ala_fraction() << "\n";
	outstr << "BEFORE_FUNCTION " << ( ala_constraint_under_strength() == 0.0 ? "CONSTANT\n" : "QUADRATIC\n" );
	outstr << "AFTER_FUNCTION " << ( ala_constraint_over_strength() == 0.0 ? "CONSTANT\n" : "QUADRATIC\n" );
	outstr << "END_PENALTY_DEFINITION\n";

	setup_string_out = outstr.str();
}

/// @brief Set up strings defining constraints penalizing too few hydrophobics.
///
void
AddHelixSequenceConstraintsMover::set_up_hydrophobic_constraints(
	std::string &setup_string_out
) const {
	std::stringstream outstr;

	outstr << "PENALTY_DEFINITION\n";
	outstr << "PROPERTIES HYDROPHOBIC\n";
	outstr << "FRACT_DELTA_START -0.01\n";
	outstr << "FRACT_DELTA_END 0.01\n";
	outstr << "PENALTIES " << hydrophobic_constraint_strength() << " 0 0\n";
	outstr << "FRACTION " << desired_min_hydrophobic_fraction() << "\n";
	outstr << "BEFORE_FUNCTION QUADRATIC\n";
	outstr << "AFTER_FUNCTION CONSTANT\n";
	outstr << "END_PENALTY_DEFINITION\n";

	setup_string_out = outstr.str();
}


/// @brief Given a pose and an empty vector (that will be cleared by this operation and repopulated), run DSSP and
/// identify all helices longer than a specified length.
/// @details Calls protocols::aa_composition::find_helices_over_length(), in util.cc/hh.
void
AddHelixSequenceConstraintsMover::find_helices(
	core::pose::Pose const &pose,
	utility::vector1 < std::pair < core::Size, core::Size > > &helices
) const {
	helices.clear();
	find_helices_over_length( pose, helices, min_helix_length() ); //In util.cc/hh.
}

/// @brief Given a pose and an already-populated list of helices, remove helices that do not have at least one residue selected by the
/// residue selector associated with this mover.
/// @details Does nothing if there's no residue selector.
void
AddHelixSequenceConstraintsMover::remove_unselected_helices(
	core::pose::Pose const &pose,
	utility::vector1 < std::pair < core::Size, core::Size > > &helices
) const {
	using namespace core::select::residue_selector;

	if ( !has_residue_selector() ) return; //Do nothing of there's no residue selector.

	utility::vector1 < std::pair < core::Size, core::Size > > new_helices_list;
	ResidueSubset const selected( residue_selector()->apply(pose) );

	for ( core::Size i=1, imax=helices.size(); i<=imax; ++i ) {
		bool helix_selected(false);
		for ( core::Size j=helices[i].first, jmax=helices[i].second; j<=jmax; ++j ) {
			if ( selected[j] ) {
				helix_selected = true;
				break;
			}
		}
		if ( helix_selected ) { new_helices_list.push_back( helices[i] ); }
	}

	helices = new_helices_list;
}

/// @brief Given a pose, add constraints for the termini requiring charges in the first or last N residues.
/// @details If n_terminus is true, this adds the requirement that the N-terminus has negative charges; if it
/// is false, this adds the requirement that the C-terminus has positive charges.
void
AddHelixSequenceConstraintsMover::do_terminal_constraints(
	core::pose::Pose &pose,
	bool const n_terminus,
	core::Size const start_res,
	core::Size const end_res,
	std::string const &constraint_setup
) const {
	core::select::residue_selector::ResidueIndexSelectorOP selector( new core::select::residue_selector::ResidueIndexSelector );
	if ( n_terminus ) {
		for ( core::Size i(start_res), imax(start_res + n_terminus_size() - 1); i<=imax; ++i ) {
			selector->append_index( i );
		}
	} else {
		for ( core::Size i(end_res - c_terminus_size() + 1), imax(end_res); i<=imax; ++i ) {
			selector->append_index( i );
		}
	}

	core::scoring::aa_composition_energy::AACompositionConstraintOP constraint( new core::scoring::aa_composition_energy::AACompositionConstraint );
	constraint->initialize_from_file_contents( constraint_setup );
	constraint->set_selector( selector );

	pose.add_constraint( utility::pointer::dynamic_pointer_cast < core::scoring::constraints::Constraint const >( constraint ) );
}

/// @brief Given a pose and helix start and end points, add constraints for the helix limiting helix-disfavouring residues.
///
void
AddHelixSequenceConstraintsMover::do_overall_or_alanine_constraints(
	core::pose::Pose &pose,
	core::Size const start_res,
	core::Size const end_res,
	std::string const &constraint_setup
) const {
	core::select::residue_selector::ResidueIndexSelectorOP selector( new core::select::residue_selector::ResidueIndexSelector );
	for ( core::Size i(start_res); i<=end_res; ++i ) {
		selector->append_index( i );
	}

	core::scoring::aa_composition_energy::AACompositionConstraintOP constraint( new core::scoring::aa_composition_energy::AACompositionConstraint );
	constraint->initialize_from_file_contents( constraint_setup );
	constraint->set_selector( selector );

	pose.add_constraint( utility::pointer::dynamic_pointer_cast < core::scoring::constraints::Constraint const >( constraint ) );
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
AddHelixSequenceConstraintsMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new AddHelixSequenceConstraintsMover );
}

std::string
AddHelixSequenceConstraintsMoverCreator::keyname() const
{
	return AddHelixSequenceConstraintsMover::mover_name();
}

void AddHelixSequenceConstraintsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddHelixSequenceConstraintsMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, AddHelixSequenceConstraintsMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //aa_composition
