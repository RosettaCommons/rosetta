// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/aa_composition/AddHelixSequenceConstraintsMover.hh
/// @brief This mover adds sequence constraints to the ends of each helix, requiring at least one positively-charged residue in the three C-terminal residues, and at least one negatively-charged resiude in the three N-terminal residues.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_aa_composition_AddHelixSequenceConstraintsMover_HH
#define INCLUDED_protocols_aa_composition_AddHelixSequenceConstraintsMover_HH

// Unit headers
#include <protocols/aa_composition/AddHelixSequenceConstraintsMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace aa_composition {

///@brief This mover adds sequence constraints to the ends of each helix, requiring at least one positively-charged residue in the three C-terminal residues, and at least one negatively-charged resiude in the three N-terminal residues.
class AddHelixSequenceConstraintsMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	AddHelixSequenceConstraintsMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	AddHelixSequenceConstraintsMover( AddHelixSequenceConstraintsMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~AddHelixSequenceConstraintsMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	//AddHelixSequenceConstraintsMover & operator=( AddHelixSequenceConstraintsMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

public: // getters

	/// @brief An optional residue selector.  If provided, only helices that have at least one residue in the selection have
	/// constraints applied.  If unused, all heliecs receive constraints.
	inline core::select::residue_selector::ResidueSelectorCOP residue_selector() const { return residue_selector_; }

	/// @brief Has a residue selector been specified for this mover?
	///
	inline bool has_residue_selector() const { return (residue_selector_ != nullptr); }

	/// @brief Get whether old composition constraints should be deleted before applying this mover.
	/// @details Default false (i.e. new composition constraints are appended to the old).
	inline bool reset_mode() const { return reset_mode_; }

	/// @brief Get the minimum number of residues in a helix for the helix to be considered.
	/// @details Very short helices (e.g. one turn, five residues) have negligible helix macrodipoles,
	/// and can be disregarded.  Defaults to 8.
	inline core::Size min_helix_length() const { return min_helix_length_; }

	/// @brief Should this mover add constraints to the N-termini of helices requiring negative charges there?
	/// @details Default true.
	inline bool add_n_terminal_constraints() const { return add_n_terminal_constraints_; }

	/// @brief Should this mover add constraints to the C-termini of helices requiring positive charges there?
	/// @details Default true.
	inline bool add_c_terminal_constraints() const { return add_c_terminal_constraints_; }

	/// @brief Should this mover add constraints to the helix penalizing helix-disfavouring residue types?
	/// @details Default true.
	inline bool add_overall_constraints() const { return add_overall_constraints_; }

	/// @brief Should this mover add constraints to the helix requiring a certain alanine fraction?
	/// @details Default true.
	inline bool add_alanine_constraints() const { return add_alanine_constraints_; }

	/// @brief Should this mover add constraints to the helix requiring at minimum fraction of hydrophobic residues?
	/// @details Default true.
	inline bool add_hydrophobic_constraints() const { return add_hydrophobic_constraints_; }

	/// @brief The number of negative charges required at the N-terminus of each helix.
	/// @details Default 2.
	inline core::Size min_n_terminal_charges() const { return min_n_terminal_charges_; }

	/// @brief The number of N-terminal residues in each helix, in which negative charges must be found.
	/// @brief Default 3.
	inline core::Size n_terminus_size() const { return n_terminus_size_; }

	/// @brief The weight given to the N-terminal constraints.
	/// @details The mover applies a quadratic penalty for having fewer than the desired number of negative
	/// charges at the N-terminus of each helix.  This is the penalty for having one fewer than the desired number.
	/// There is no penalty or bonus for having more than the desired number.
	inline core::Real const & n_terminal_constraint_strength() const { return n_terminal_constraint_strength_; }

	/// @brief The number of positive charges required at the C-terminus of each helix.
	/// @details Default 2.
	inline core::Size min_c_terminal_charges() const { return min_c_terminal_charges_; }

	/// @brief The number of C-terminal residues in each helix, in which positive charges must be found.
	/// @brief Default 3.
	inline core::Size c_terminus_size() const { return c_terminus_size_; }

	/// @brief The weight given to the C-terminal constraints.
	/// @details The mover applies a quadratic penalty for having fewer than the desired number of positive
	/// charges at the C-terminus of each helix.  This is the penalty for having one fewer than the desired number.
	/// There is no penalty or bonus for having more than the desired number.
	inline core::Real const & c_terminal_constraint_strength() const { return c_terminal_constraint_strength_; }

	/// @brief Access the vector of types to avoid.
	///
	inline utility::vector1< std::string > const & types_to_avoid() const { return types_to_avoid_; }

	/// @brief Access an element in the vector of types to avoid.
	/// @details Note: no bounds checking!
	inline std::string const & types_to_avoid( core::Size const index ) const { return types_to_avoid_[index]; }

	/// @brief Max number of residues from the types to avoid list per helix.
	///
	inline core::Size overall_max_count() const { return overall_max_count_; }

	/// @brief Get the strength of the helix constraints that penalize helix-disfavouring residue types.
	/// @details This is the penalty for having one more than the maximum count of the disfavoured residues.
	/// Past this, the penalty ramps quadratically.
	inline core::Real const & overall_constraints_strength() const { return overall_constraints_strength_; }

	/// @brief Get the desired fractional content of alanine in each helix.
	/// @details Penalties ramp quadratically.
	inline core::Real const & desired_ala_fraction() const { return desired_ala_fraction_; }

	/// @brief Get the penalty for having 1% less than the desired fractional content of alanine.
	/// @details Penalties ramp quadratically.
	inline core::Real const & ala_constraint_under_strength() const { return ala_constraint_under_strength_; }

	/// @brief Get the penalty for having 1% more than the desired fractional content of alanine.
	/// @details Penalties ramp quadratically.
	inline core::Real const & ala_constraint_over_strength() const { return ala_constraint_over_strength_; }

	/// @brief Get the desired minimum hydrophobic fraction.
	///
	inline core::Real const & desired_min_hydrophobic_fraction() const { return desired_min_hydrophobic_fraction_; }

	/// @brief Get the strength of the hydrophobic constraint.  (The penalty for having more than 1% less than the
	/// desired minimum hydrophobic content.
	inline core::Real const & hydrophobic_constraint_strength() const { return hydrophobic_constraint_strength_; }

public: // setters

	/// @brief Set a residue selector if specified by the user
	///
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in);

	/// @brief Set whether old composition constraints should be deleted before applying this mover.
	///
	void set_reset_mode( bool const setting );

	/// @brief Set the minimum number of residues in a helix for the helix to be considered.
	/// @details Very short helices (e.g. one turn, five residues) have negligible helix macrodipoles,
	/// and can be disregarded.  Defaults to 8.
	void set_min_helix_length( core::Size const setting );

	/// @brief Set whether this mover should add constraints to the N-termini of helices requiring negative charges there.
	///
	void set_add_n_terminal_constraints( bool const setting, core::Size const min_charges, core::Size const nterm_size, core::Real const &strength );

	/// @brief Set whether this mover should add constraints to the C-termini of helices requiring positive charges there.
	/// @details Default true.
	void set_add_c_terminal_constraints( bool const setting, core::Size const min_charges, core::Size const cterm_size, core::Real const &strength );

	/// @brief Set whether this mover should add constraints to the helix penalizing helix-disfavouring residue types.
	/// @details Default true.
	void set_add_overall_constraints( bool const setting );

	/// @brief Set whether this mover should add constraints to the helix penalizing helix-disfavouring residue types.
	/// @details Default true.
	void set_add_overall_constraints( bool const setting, std::string const &types_to_avoid, core::Size const overall_max_count, core::Real const &overall_constraints_strength );

	/// @brief Set whether this mover should add constraints to the helix requiring a certain alanine fraction.
	/// @details Default true.
	void set_add_alanine_constraints( bool const setting );

	/// @brief Set whether this mover should add constraints to the helix requiring a certain alanine fraction.
	/// @details Default true.
	void set_add_alanine_constraints( bool const setting, core::Real const &desired_ala_fraction, core::Real const &ala_constraint_under_strength, core::Real const &ala_constraint_over_strength );

	/// @brief Set whether this mover should add constraints to the helix requiring at minimum fraction of hydrophobic residues.
	/// @details Default true.
	void set_add_hydrophobic_constraints( bool const setting );

	/// @brief Set whether this mover should add constraints to the helix requiring at minimum fraction of hydrophobic residues.
	/// @details Default true.
	void set_add_hydrophobic_constraints( bool const setting, core::Real const &desired_min_hydrophobic_fraction, core::Real const &hydrophobic_constraint_strength );

private: // methods

	/// @brief Once the mover has been set up, check that the lengths of termini is sensible compared to the min_helix_length setting.
	/// @details Throws an error if add_n_terminal_constraints or add_c_terminal_constraints is true, and the min_helix_length is less
	/// than the lesser of the active termini lenghts.
	void check_lengths_sensible() const;

	/// @brief Initialize the types to avoid to ASN, ASP, SER, GLU, THR, VAL.
	///
	void initialize_types_to_avoid();

	/// @brief Given a pose, remove all sequence constraints from it.
	///
	void delete_old_sequence_constraints( core::pose::Pose &pose ) const;

	/// @brief Set up strings defining terminal constraints.
	///
	void set_up_terminal_constraints( std::string &setup_string_out, bool const n_terminal ) const;

	/// @brief Set up strings defining constraints penalizing helix-disfavouring residue types.
	///
	void set_up_overall_constraints( std::string &setup_string_out ) const;

	/// @brief Set up strings defining constraints penalizing too few or too many alanines.
	///
	void set_up_alanine_constraints( std::string &setup_string_out ) const;

	/// @brief Set up strings defining constraints penalizing too few hydrophobics.
	///
	void set_up_hydrophobic_constraints( std::string &setup_string_out ) const;

	/// @brief Given a pose and an empty vector (that will be cleared by this operation and repopulated), run DSSP and
	/// identify all helices longer than a specified length.
	/// @details Calls protocols::aa_composition::find_helices_over_length(), in util.cc/hh.
	void find_helices( core::pose::Pose const &pose, utility::vector1 < std::pair < core::Size, core::Size > > &helices ) const;

	/// @brief Given a pose and an already-populated list of helices, remove helices that do not have at least one residue selected by the
	/// residue selector associated with this mover.
	/// @details Does nothing if there's no residue selector.
	void remove_unselected_helices( core::pose::Pose const &pose, utility::vector1 < std::pair < core::Size, core::Size > > &helices ) const;

	/// @brief Given a pose and helix start and end points, add constraints for the helix termini requiring charges in the first or last N residues.
	/// @details If n_terminus is true, this adds the requirement that the N-terminus has negative charges; if it
	/// is false, this adds the requirement that the C-terminus has positive charges.
	void do_terminal_constraints( core::pose::Pose &pose, bool const n_terminus, core::Size const start_res, core::Size const end_res, std::string const &constraint_setup ) const;

	/// @brief Given a pose and helix start and end points, add constraints for the helix limiting helix-disfavouring residues.
	///
	void do_overall_or_alanine_constraints( core::pose::Pose &pose, core::Size const start_res, core::Size const end_res, std::string const &constraint_setup ) const;

private: // data

	/// @brief Should old composition constraints be deleted before applying this mover?
	/// @details Default false (i.e. new composition constraints are appended to the old).
	bool reset_mode_;

	/// @brief The minimum number of residues in a helix for the helix to be considered.
	/// @details Very short helices (e.g. one turn, five residues) have negligible helix macrodipoles,
	/// and can be disregarded.  Defaults to 8.
	core::Size min_helix_length_;

	/// @brief Should this mover add constraints to the N-termini of helices requiring negative charges there?
	/// @details Default true.
	bool add_n_terminal_constraints_;

	/// @brief Should this mover add constraints to the C-termini of helices requiring positive charges there?
	/// @details Default true.
	bool add_c_terminal_constraints_;

	/// @brief Should this mover add constraints to the helix penalizing helix-disfavouring residue types?
	/// @details Default true.
	bool add_overall_constraints_;

	/// @brief Should this mover add constraints to the helix requiring a certain alanine fraction?
	/// @details Default true.
	bool add_alanine_constraints_;

	/// @brief Should this mover add constraints to the helix requiring at minimum fraction of hydrophobic residues?
	/// @details Default true.
	bool add_hydrophobic_constraints_;

	/// @brief The number of negative charges required at the N-terminus of each helix.
	/// @details Default 2.
	core::Size min_n_terminal_charges_;

	/// @brief The number of N-terminal residues in each helix, in which negative charges must be found.
	/// @brief Default 3.
	core::Size n_terminus_size_;

	/// @brief The weight given to the N-terminal constraints.
	/// @details The mover applies a quadratic penalty for having fewer than the desired number of negative
	/// charges at the N-terminus of each helix.  This is the penalty for having one fewer than the desired number.
	/// There is no penalty or bonus for having more than the desired number.
	core::Real n_terminal_constraint_strength_;

	/// @brief The number of positive charges required at the C-terminus of each helix.
	/// @details Default 2.
	core::Size min_c_terminal_charges_;

	/// @brief The number of C-terminal residues in each helix, in which positive charges must be found.
	/// @brief Default 3.
	core::Size c_terminus_size_;

	/// @brief The weight given to the C-terminal constraints.
	/// @details The mover applies a quadratic penalty for having fewer than the desired number of positive
	/// charges at the C-terminus of each helix.  This is the penalty for having one fewer than the desired number.
	/// There is no penalty or bonus for having more than the desired number.
	core::Real c_terminal_constraint_strength_;

	/// @brief The types, listed as three-letter codes, that should be avoided in helices.
	/// @details The overall constraints cap the number of residues of these types that are allowed.
	/// Defaults to ASN, ASP, SER, GLY, THR, VAL.
	utility::vector1 < std::string > types_to_avoid_;

	/// @brief The maximum number of residues from the types_to_avoid_ list that will be tolerated per helix.
	/// @details Default 0.
	core::Size overall_max_count_;

	/// @brief The penalty strength of the overall penalty, which penalizes types from the types_to_avoid_ list.
	/// @details  The mover applies a quadratic penalty for having more than the maximum number of residues of these types.
	/// There is no penalty or bonus for having fewer.  This is the penalty for having one more than the max tolerated.
	/// Default 5.0.
	core::Real overall_constraints_strength_;

	/// @brief The desired fractional content of alanine in each helix.
	///
	core::Real desired_ala_fraction_;

	/// @brief The penalty for having 1% less than the desired fractional content of alanine.
	/// @details Default 0.2.
	core::Real ala_constraint_under_strength_;

	/// @brief The penalty for having 1% more than the desired fractional content of alanine.
	/// @details Default 0.2.
	core::Real ala_constraint_over_strength_;

	/// @brief The desired minimum hydrophobic content of helices.
	/// @details Defaults to 0.25.
	core::Real desired_min_hydrophobic_fraction_;

	/// @brief The strength of the hydrophobic constraint.  (The penalty for having more than 1% less than the
	/// desired minimum hydrophobic content.
	core::Real hydrophobic_constraint_strength_;

	/// @brief An optional residue selector.  If provided, only helices that have at least one residue in the selection have
	/// constraints applied.  If unused, all heliecs receive constraints.
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;


};

std::ostream &
operator<<( std::ostream & os, AddHelixSequenceConstraintsMover const & mover );

} //protocols
} //aa_composition

#endif //protocols_aa_composition_AddHelixSequenceConstraintsMover_HH
