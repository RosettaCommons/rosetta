// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/ThreefoldLinkerMover.hh
/// @brief This mover links three cysteine residues with a three-way cross-linker.  It adds the crosslinker,
/// sets up constraints, optionally packs and energy-mimizes it into place (packing/minimizing only the crosslinker
/// and the side-chains to which it connects), andthen optionally relaxes the whole structure.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_ThreefoldLinkerMover_hh
#define INCLUDED_protocols_cyclic_peptide_ThreefoldLinkerMover_hh

// Unit headers
#include <protocols/cyclic_peptide/ThreefoldLinkerMover.fwd.hh>
#include <protocols/cyclic_peptide/threefold_linker/ThreefoldLinkerMoverHelper.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace cyclic_peptide {

/// @brief Enumeration of allowed crosslinkers.
enum CrossLinker {
	//When adding an effect to this enum, add its name to the get_crosslinker_name() function, and its helper to the switch statement
	//in the apply() function.
	no_crosslinker = 1, //Keep this first

	TBMB,

	unknown_crosslinker, //Keep this second-to-last.
	end_of_crosslinker_list = unknown_crosslinker //Keep this last.
};

///@brief This mover links three cysteine residues with a threefld-symmetric cross-linker.  It adds the crosslinker, sets up constraints, optionally packs and energy-mimizes it into place (packing/minimizing only the crosslinker and the side-chains to which it connects), andthen optionally relaxes the whole structure.
class ThreefoldLinkerMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ThreefoldLinkerMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	ThreefoldLinkerMover( ThreefoldLinkerMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~ThreefoldLinkerMover();

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	virtual void
	apply( core::pose::Pose & pose );

	/// @brief Show the contents of the Mover
	static std::string
	class_name();

	/// @brief Returns the name of this Mover.
	static
	std::string
	mover_name();

	virtual void
	show( std::ostream & output = std::cout ) const;

	/// @brief Get the name of the Mover
	virtual std::string
	get_name() const;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief Parse XML tag (to use this Mover in Rosetta Scripts).
	///
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	/// @brief Provide information on what options are available in XML tag.
	///
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	//ThreefoldLinkerMover & operator=( ThreefoldLinkerMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

	/// @brief Set the residue selector to use.
	///
	void set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	/// @brief Get the residue selector that this mover uses.
	inline core::select::residue_selector::ResidueSelectorCOP residue_selector() const { return residue_selector_; }

	/// @brief Set the linker name.
	///
	void set_linker_name( std::string const &name_in );

	/// @brief Get the linker name.
	///
	std::string linker_name() const;

	/// @brief Set the behaviour of this mover.
	///
	void set_behaviour( bool const add_linker, bool const constrain_linker, bool const pack_and_minimize_linker_and_sidechains, bool const do_final_fastrelax, bool const threefold_symmetric );

	/// @brief Set the filtering behaviour of this mover.
	///
	void set_filter_behaviour( bool const filter_by_sidechain_distance, bool const filter_by_constraints_energy, bool const filter_by_total_score, core::Real const &filter_by_total_score_cutoff_energy, core::Real const & sidechain_distance_filter_multiplier, core::Real const & constraints_energy_filter_multiplier );

	/// @brief Get whether we're adding the linker.
	///
	inline bool add_linker() const { return add_linker_; }

	/// @brief Get whether we're adding the linker constraints.
	///
	inline bool constrain_linker() const { return constrain_linker_; }

	/// @brief Get whether we're packing and minimizing the linker and the side-chains to which it is connected.
	///
	inline bool pack_and_minimize_linker_and_sidechains() const { return pack_and_minimize_linker_and_sidechains_; }

	/// @brief Get whether we're doing a final FastRelax on the whole structure.
	///
	inline bool do_final_fastrelax() const { return do_final_fastrelax_; }

	/// @brief Get whether this is being applied to a threefold-symmetric system.
	///
	inline bool threefold_symmetric() const { return threefold_symmetric_; }

	/// @brief Set the scorefunction to use for packing and minimization.
	/// @details Cloned at apply time.  (That is, the scorefunction is shared until apply time).
	void set_scorefxn( core::scoring::ScoreFunctionCOP sfxn_in );

	/// @brief Get the scorefunction to use for packing and minimization.
	///
	core::scoring::ScoreFunctionCOP scorefxn() const;

	/// @brief Set the number of rounds of FastRelax to apply when minimizing the linker and the
	/// side-chains that connect to it.
	void set_sidechain_frlx_rounds( core::Size const rounds_in );

	/// @brief Set the number of rounds of FastRelax to apply at the end.
	///
	void set_final_frlx_rounds( core::Size const rounds_in );

	/// @brief Get the number of rounds of FastRelax to apply when minimizing the linker and the
	/// side-chains that connect to it.
	inline core::Size sidechain_frlx_rounds() const { return sidechain_frlx_rounds_; }

	/// @brief Get the number of rounds of FastRelax to apply at the end.
	///
	inline core::Size final_frlx_rounds() const { return final_frlx_rounds_; }


private: // methods

	/// @brief Apply the mover to a symmetric pose.
	/// @details Requires threefold symmetry in the pose, and threefold_symmetric_ = true.
	void symmetric_apply( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper );

	/// @brief Determine whether the residues to be crosslinked are too far apart.  This version is for symmetric poses.
	/// @details Returns TRUE for failure (too far apart), FALSE for success.
	bool filter_by_sidechain_distance_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper ) const;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	bool filter_by_constraints_energy_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper, bool const linker_was_added ) const;

	/// @brief Given a selection of exactly three residues, add a crosslinker, align it crudely to the
	/// selected residues, and set up covalent bonds.  This version is for symmetric poses.
	void add_linker_symmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper ) const;

	/// @brief Given a selection of exactly three residues that have already been connected to a crosslinker,
	/// add constraints for the crosslinker.  This version is for symmetric poses.
	void add_linker_constraints_symmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper, bool const linker_was_added ) const;

	/// @brief Apply the mover to an asymmetric pose.
	/// @details Requires and asymmetric pose, and threefold_symmetric_ = false.
	void asymmetric_apply( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper );

	/// @brief Determine whether the residues to be crosslinked are too far apart.
	/// @details Returns TRUE for failure (too far apart), FALSE for success.
	bool filter_by_sidechain_distance_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper ) const;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	bool filter_by_constraints_energy_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper ) const;

	/// @brief Determine whether the overall system has too high an overall score (including constraints) at the end of the protocol.
	/// @details Returns TRUE for failure (too high an overall score) and FALSE for success.
	bool filter_by_total_score( core::pose::Pose const &pose ) const;

	/// @brief Given a selection of exactly three residues, add a crosslinker, align it crudely to the
	/// selected residues, and set up covalent bonds.
	void add_linker_asymmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper ) const;

	/// @brief Given a selection of exactly three residues that have already been connected to a crosslinker,
	/// add constraints for the crosslinker.
	void add_linker_constraints_asymmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper ) const;

	/// @brief Repack and minimize the sidechains.
	/// @details Also repacks and minimzes the linker, letting the jump vary.
	void pack_and_minimize_linker_and_sidechains( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, protocols::cyclic_peptide::threefold_linker::ThreefoldLinkerMoverHelperCOP helper, bool const whole_structure, bool const symmetric ) const;

	/// @brief Given a selection from a ResidueSelector, check that exactly three residues have been selected.
	/// @details Returns true if exactly three have been selected, false otherwise.
	bool exactly_three_selected( core::select::residue_selector::ResidueSubset const & selection ) const;

	/// @brief Given a ResidueSubset with exactly three residues selected, pull out the three indices.
	///
	void get_cys_indices( core::select::residue_selector::ResidueSubset const & selection, core::Size &cys1, core::Size &cys2, core::Size &cys3 ) const;

	/// @brief Given a CrossLinker enum, get its name.
	///
	std::string get_crosslinker_name( CrossLinker const crosslinker ) const;

	/// @brief Given a CrossLinker name, get its enum.
	///
	CrossLinker get_crosslinker_enum( std::string const &name ) const;

	/// @brief Are we filtering by sidechain distance before placing the linker?
	///
	inline bool filter_by_sidechain_distance() const { return filter_by_sidechain_distance_; }

	/// @brief Are we filtering by constraints energy after placing and minimizing the linker?
	///
	inline bool filter_by_constraints_energy() const { return filter_by_constraints_energy_; }

	/// @brief Are we filtering by total score plus constraints energy at the end of the protocol?
	///
	inline bool filter_by_total_score() const { return filter_by_total_score_; }

	/// @brief The cutoff energy to use, if we're filtering by the total score plus constraints energy at the end of the protocol.
	///
	inline core::Real const & filter_by_total_score_cutoff_energy() const { return filter_by_total_score_cutoff_energy_; }

	/// @brief The distance filter multiplier.
	///
	inline core::Real const & sidechain_distance_filter_multiplier() const { return sidechain_distance_filter_multiplier_; }

	/// @brief The constraints energy filter multiplier.
	///
	inline core::Real const & constraints_energy_filter_multiplier() const { return constraints_energy_filter_multiplier_; }

private: // data

	/// @brief A residue selector to select exactly three cysteine residues to cross-link.
	///
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;

	/// @brief The three-way crosslinker to use.
	/// @details Must be set.  See definition of CrossLinker enum in ThrefoldLinkerMover.hh.
	CrossLinker linker_;

	/// @brief Should we add the linker?
	/// @details Default true.
	bool add_linker_;

	/// @brief Should we set up constraints for the linker?
	/// @details Default true.
	bool constrain_linker_;

	/// @brief Should we pack and minimize the linker and the side-chains to which it connects?
	/// @details Default true.
	bool pack_and_minimize_linker_and_sidechains_;

	/// @brief Should we FastRelax the whole structure?
	/// @details Default false.
	bool do_final_fastrelax_;

	/// @brief Is this a threefold-symmetric system?
	/// @details Default false.
	bool threefold_symmetric_;

	/// @brief The scorefunction to use for packing and minimization.
	/// @details Cloned at apply time.
	core::scoring::ScoreFunctionCOP sfxn_;

	/// @brief The number of rounds of FastRelax to apply when packing and minimizing
	/// side-chains and the linker.
	/// @details Default 3.
	core::Size sidechain_frlx_rounds_;

	/// @brief The number of rounds of FastRelax to apply at the end of the protocol.
	/// @details Default 3.
	core::Size final_frlx_rounds_;

	/// @brief Should we filter based on initial distance between residues that are to be linked?
	/// @details Default true.  Each linker helper has its own evaluator to decide whether or not to discard poses.
	bool filter_by_sidechain_distance_;

	/// @brief Should we filter based on constraints energy after the initial minimization of the linker?
	/// @details Default true.  Each linker helper has its own evaluator for this.
	bool filter_by_constraints_energy_;

	/// @brief Should we filter by total score (including constraints energy) at the end?
	/// @details Default false.
	bool filter_by_total_score_;

	/// @brief If we're filtering by total score, this is the cutoff energy.
	/// @details Defaults to 0.
	core::Real filter_by_total_score_cutoff_energy_;

	/// @brief Multiplier to affect stringency of sidechain distance filter.  Higher values are more permissive.
	/// @details Default 1.0.
	core::Real sidechain_distance_filter_multiplier_;

	/// @brief Multiplier to affect stringency of constraints energy filter.  Higher values are more permissive.
	/// @details Default 1.0.
	core::Real constraints_energy_filter_multiplier_;


};

std::ostream &
operator<<( std::ostream & os, ThreefoldLinkerMover const & mover );

} //protocols
} //cyclic_peptide

#endif //protocols/cyclic_peptide_ThreefoldLinkerMover_hh
