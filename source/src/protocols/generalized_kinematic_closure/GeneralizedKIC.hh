// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/generalized_kinematic_closure/GeneralizedKIC.hh
/// @brief  Headers for GeneralizedKIC mover class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_generalized_kinematic_closure_GeneralizedKIC_hh
#define INCLUDED_protocols_generalized_kinematic_closure_GeneralizedKIC_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.fwd.hh>
#include <protocols/generalized_kinematic_closure/perturber/GeneralizedKICperturber.hh>
#include <protocols/generalized_kinematic_closure/filter/GeneralizedKICfilter.hh>
#include <protocols/generalized_kinematic_closure/selector/GeneralizedKICselector.hh>
#include <core/scoring/Ramachandran.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/filters/ContingentFilter.fwd.hh>
#include <protocols/filters/ContingentFilter.hh>

// Project Headers
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.fwd.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace generalized_kinematic_closure {

class GeneralizedKIC : public protocols::moves::Mover
{
public:
	GeneralizedKIC();
	GeneralizedKIC(GeneralizedKIC const &src);
	~GeneralizedKIC() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;


	/// @brief Actually apply the mover to the pose.
	void apply(core::pose::Pose & pose) override;

	// XRW TEMP  std::string get_name() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const &
	) override;


	/// @brief Add a residue (by index in the pose) to the list of residues making up the loop to be closed.
	void add_loop_residue( core::Size const residue_index );

	/// @brief Clears the list of loop residues
	void clear_loop_residues();

	/// @brief Clears the list of perturber residues
	void clear_perturber_residue_lists();

	/// @brief Clears the list of perturber residues for perturber perturber_idx
	void clear_perturber_residue_list( core::Size const perturber_idx );

	/// @brief Add a residue (by index in the pose) to the list of residues making up the tails attached to
	/// the loop to be closed.
	/// @details  "Tails" are residues that are not part of the loop to be closed, but which are attached to
	/// the loop and which "come along for the ride" as the loop moves.
	void add_tail_residue( core::Size const residue_index );

	/// @brief Check that loop residues haven't been specified multiple times.
	/// @details Also checks that the loop residues are all within the pose.
	void check_loop_residues_sensible( core::pose::Pose const &pose);

	/// @brief Check that tail residues haven't been specified multiple times,
	/// and that the tail residues don't overlap with loop residues.
	/// @details Also checks that the tail residues are all within the pose.
	void check_tail_residues_sensible( core::pose::Pose const &pose );

	// @brief Function to set the effect of this mover on parts of the pose that are covalently attached to
	// the loop to be closed, but which aren't part of it.  Settings are:
	// 0 -- Moves the loop only; can pull apart covalent bonds to anything outside of the loop that isn't
	//      an anchor point.
	// 1 -- Moves the loop and anything downstream of the loop in the foldtree.  Can still pull apart
	//    connections to non-child geometry.
	// CURRENTLY DEPRECATED
	//void set_mover_effect_on_bonded_geometry( core::Size const effect ); //TODO -- Make this do something.

	// @brief Function to get the effect of this mover on parts of the pose that are covalently attached to
	// the loop to be closed, but which aren't part of it.  Return values are:
	// 0 -- Moves the loop only; can pull apart covalent bonds to anything outside of the loop that isn't
	//      an anchor point.
	// 1 -- Moves the loop and anything downstream of the loop in the foldtree.  Can still pull apart
	//    connections to non-child geometry.
	// CURRENTLY DEPRECATED
	//core::Size get_mover_effect_on_bonded_geometry() { return effect_on_bonded_geometry_; }

	/// @brief Set whether or not this mover builds ideal geometry for the loop, or uses the existing geometry,
	/// imperfect bond lengths/angles and all.
	void set_build_ideal_geometry( bool const buildideal ) { build_ideal_geometry_ = buildideal; return; }


	/// @brief Function to set the pivot atoms for kinematic closure:
	void set_pivot_atoms(
		core::Size const rsd1,
		std::string const &at1,
		core::Size const rsd2,
		std::string const &at2,
		core::Size const rsd3,
		std::string const &at3
	);

	/// @brief Get whether or not this mover builds ideal geometry for the loop, or uses the existing geometry,
	/// imperfect bond lengths/angles and all.
	inline bool get_build_ideal_geometry() { return build_ideal_geometry_; }

	/// @brief Returns whether or not the last call to apply() sucessfully closed the loop.  Defaults to false if
	/// apply() hasn't yet been called.
	inline bool last_run_successful() { return last_run_successful_; }

	/// @brief Tells GeneralizedKIC to close a bond, setting bond length, bond angle, and bond torsion values.  This
	/// actually just adds appropriate set_dihedral, set_bondangle, and set_bondlength perturbers to the perturber
	/// list.  Note that subsequent perturbers OR the closure itself can overwrite the bond length, bond angle, or
	/// torsion angles set here.
	/// @details
	/// @param[in] rsd1 -- The index of the first atom's residue (indexed based on residue indices in the original pose).
	/// @param[in] at1 -- The name of the first atom defining the bond to be closed.
	/// @param[in] rsd2 -- The index of the second atom's residue (indexed based on residue indices in the original pose).
	/// @param[in] at1 -- The name of the second atom defining the bond to be closed.
	/// @param[in] bondlength -- The length of the bond between the two atoms.
	/// @param[in] bondangle1 -- The bond angle defined by (atom preceding at1 in the chain to be closed), (atm1), (atm2).
	/// @param[in] bondangle1 -- The bond angle defined by (atm1), (atm2), (atom following at2 in the chain to be closed).
	/// @param[in] torsion -- The torsion angle defined by (atom preceding at1 in the chain to be closed), (atm1), (atm2), (atom following at2 in the chain to be closed).
	void close_bond (
		core::Size const rsd1,
		std::string const &at1,
		core::Size const rsd2,
		std::string const &at2,
		core::Size const rsd1_before,
		std::string const &at1_before,
		core::Size const rsd2_after,
		std::string const &at2_after,
		core::Real const &bondlength,
		core::Real const &bondangle1,
		core::Real const &bondangle2,
		core::Real const &torsion,
		bool const randomize_this_torsion,
		bool const randomize_flanking_torsions
	);

	/// @brief Set the selector (the algorithm controlling how a solution will be chosen
	/// from among the solutions passing filters).
	void set_selector_type ( selector::selector_type const &stype);

	/// @brief Set the selector (the algorithm controlling how a solution will be chosen
	/// from among the solutions passing filters).  This sets the selector by name.
	void set_selector_type ( std::string const &stypename);

	/// @brief Set the selector's scorefunction.
	///
	void set_selector_scorefunction ( core::scoring::ScoreFunctionOP sfxn );

	/// @brief Set the selector's Boltzmann temperature value.
	///
	void set_selector_kbt ( core::Real const &kbt );

	/// @brief Get the selector (the algorithm controlling how a solution will be chosen
	/// from among the solutions passing filters).
	inline selector::selector_type get_selector_type () const { return selector_->get_selector_type(); }

	/// @brief Get the name of the selector (the algorithm controlling how a solution will be chosen
	/// from among the solutions passing filters).
	inline std::string get_selector_type_name () const { return selector_->get_selector_type_name( selector_->get_selector_type() ); }


	/// @brief Add a new perturber to the list of perturbers.
	void add_perturber ();


	/// @brief Add a new perturber to the list of perturbers, setting the effect.
	void add_perturber ( perturber::perturber_effect const &effect );


	/// @brief Add a new perturber to the list of perturbers, setting the effect by effect name string.
	void add_perturber ( std::string const &effectname );

	/// @brief Set a perturber's effect.
	/// @details
	///
	/// @param[in] perturber_index -- The index in the list of perturbers already added.
	/// @param[in] effect -- The perturber effect type, based on the perturber::perturber_effect enum (e.g. set_dihedral, randomize_backbone, etc.).
	void set_perturber_effect ( core::Size const perturber_index, perturber::perturber_effect const &effect );

	/// @brief Initialize a perturber's BinTransitionCalculator object, and load a bin_params file.
	///
	void load_perturber_bin_params( core::Size const perturber_index, std::string const &bin_params_file );

	/// @brief Initialize a perturber's BinTransitionCalculator object, and load a bin_params file.
	/// @details This acts on the last perturber in the perturber list.
	void load_perturber_bin_params( std::string const &bin_params_file );

	/// @brief Set the number of iterations for a perturber.
	///
	void set_perturber_iterations( core::Size const perturber_index, core::Size const val );

	/// @brief Set the number of iterations for a perturber.
	/// @details This acts on the last perturber in the perturber list.
	void set_perturber_iterations( core::Size const val );

	/// @brief Set whether the perturb_backbone_by_bins perturber requires residues to change their torsion
	/// bins every move, or whether they can stay within the same bin.
	void set_perturber_must_switch_bins( core::Size const perturber_index, bool const val );

	/// @brief Set whether the perturb_backbone_by_bins perturber requires residues to change their torsion
	/// bins every move, or whether they can stay within the same bin.
	/// @details This acts on the last perturber in the perturber list.
	void set_perturber_must_switch_bins( bool const val );

	/// @brief Set the bin for the set_backbone_bin perturber.
	///
	void set_perturber_bin( core::Size const perturber_index, std::string const &bin );

	/// @brief Set the bin for the set_backbone_bin perturber.
	/// @details This acts on the last perturber in the perturber list.
	void set_perturber_bin( std::string const &bin );

	/// @brief Set the custom Ramachandran table for the randomize_alpha_backbone_by_rama perturber.
	///
	void set_perturber_custom_rama_table( core::Size const perturber_index, std::string const &table_name );

	/// @brief Set the custom Ramachandran table for the randomize_alpha_backbone_by_rama perturber.
	/// @details This version works by Rama_Table_Type.
	void set_perturber_custom_rama_table( core::Size const perturber_index, core::scoring::Rama_Table_Type const table_type );

	/// @brief Set the curstom Ramachandran table for the randomize_alpha_backbone_by_rama perturber.
	/// @details This acts on the last perturber in the perturber list.
	void set_perturber_custom_rama_table( std::string const &table_name );

	/// @brief Set the curstom Ramachandran table for the randomize_alpha_backbone_by_rama perturber.
	/// @details This acts on the last perturber in the perturber list, and works by Rama_Table_Type.
	void set_perturber_custom_rama_table( core::scoring::Rama_Table_Type const table_type );

	/// @brief Set whether the perturber's generated poses should be used for BOINC graphics.
	/// @details Does nothing outside of the BOINC build.
	void set_perturber_attach_boinc_ghost_observer( core::Size const perturber_index, bool const setting );

	/// @brief Set whether the perturber's generated poses should be used for BOINC graphics.
	/// @details Does nothing outside of the BOINC build.  This version acts on the last perturber in the perturber list.
	void set_perturber_attach_boinc_ghost_observer( bool const setting );

	/// @brief Add a value to the list of values that a perturber takes.
	void add_value_to_perturber_value_list ( core::Size const perturber_index, core::Real const &val );


	/// @brief Add a value to the list of values that a perturber takes.  This operates on the last perturber in the perturber list.
	void add_value_to_perturber_value_list ( core::Real const &val );


	/// @brief Add a residue to the list of residues that a perturber takes.  Note that residue_index is based on indices of the ORIGINAL POSE,
	/// not the loop in isolation.
	void add_residue_to_perturber_residue_list ( core::Size const perturber_index, core::Size const residue_index );

	/// @brief Add a residue to the list of residues that a perturber takes.  Note that residue_index is based on indices of the ORIGINAL POSE,
	/// not the loop in isolation.  This version acts on the last perturber added.
	void add_residue_to_perturber_residue_list ( core::Size const residue_index );


	/// @brief Add a set of AtomIDs to the list of sets of AtomIDs that a perturber takes.
	void add_atomset_to_perturber_atomset_list ( core::Size const perturber_index, utility::vector1 < core::id::NamedAtomID > const &atomset );


	/// @brief Add a set of AtomIDs to the list of sets of AtomIDs that a perturber takes.  This operates on the last perturber in the perturber list.
	void add_atomset_to_perturber_atomset_list ( utility::vector1 < core::id::NamedAtomID > const &atomset );


	/// @brief Add a new filter to the list of filters.
	void add_filter ();


	/// @brief Add a new filter to the list of filters, setting the filter type.
	void add_filter ( filter::filter_type const &filtertype );


	/// @brief Add a new filter to the list of filters, setting the filter type by name.
	void add_filter ( std::string const &filtertypename );


	/// @brief Add a real-valued parameter to a filter's parameter list.
	void add_filter_parameter (core::Size const filter_index, std::string const &param_name, core::Real const &value);


	/// @brief Add an integer-valued parameter to a filter's parameter list.
	void add_filter_parameter (core::Size const filter_index, std::string const &param_name, core::Size const value);


	/// @brief Add a Boolean-valued parameter to a filter's parameter list.
	void add_filter_parameter (core::Size const filter_index, std::string const &param_name, bool const value);


	/// @brief Add a string-valued parameter to a filter's parameter list.
	void add_filter_parameter (core::Size const filter_index, std::string const &param_name, std::string const &value);


	/// @brief Add a real-valued parameter to the last filter's parameter list.
	void add_filter_parameter (std::string const &param_name, core::Real const &value);


	/// @brief Add an integer-valued parameter to the last filter's parameter list.
	void add_filter_parameter (std::string const &param_name, core::Size const value);


	/// @brief Add a Boolean-valued parameter to the last filter's parameter list.
	void add_filter_parameter (std::string const &param_name, bool const value);


	/// @brief Add a string-valued parameter to the last filter's parameter list.
	void add_filter_parameter (std::string const &param_name, std::string const &value);

	/// @brief Set the residue number that a backbone_bin filter is acting on.
	///
	void set_filter_resnum( core::Size const filter_index, core::Size const value );

	/// @brief Set the residue number that a backbone_bin filter is acting on.
	/// @details This version acts on the last filter in the filter list.
	void set_filter_resnum( core::Size const value );

	/// @brief Set the bin name that a backbone_bin filter is looking for.
	///
	void set_filter_bin( core::Size const filter_index, std::string const &name_in );

	/// @brief Set the bin name that a backbone_bin filter is looking for.
	/// @details This version acts on the last filter in the filter list.
	void set_filter_bin( std::string const &name_in );

	/// @brief Set the rama term cutoff energy for the alpha_aa_rama_check filter.
	///
	void set_filter_rama_cutoff_energy( core::Size const filter_index, core::Real const &cutoff_energy );

	/// @brief Set the rama term cutoff energy for the alpha_aa_rama_check filter.
	/// @details This version acts on the last filter in the filter list.
	void set_filter_rama_cutoff_energy( core::Real const &cutoff_energy );

	/// @brief Set whether the filter's generated poses should be used for BOINC graphics.
	/// @details Does nothing outside of the BOINC build.
	void set_filter_attach_boinc_ghost_observer( core::Size const filter_index, bool const setting );

	/// @brief Set whether the filter's generated poses should be used for BOINC graphics.
	/// @details Does nothing outside of the BOINC build.  This version acts on the last filter in the filter list.
	void set_filter_attach_boinc_ghost_observer( bool const setting );

	/// @brief Initialize a filter's BinTransitionCalculator object, and load a bin_params file.
	///
	void load_filter_bin_params( core::Size const filter_index, std::string const &bin_params_file );

	/// @brief Initialize a filter's BinTransitionCalculator object, and load a bin_params file.
	/// @details This acts on the last filter in the filter list.
	void load_filter_bin_params( std::string const &bin_params_file );


	/// @brief Set the number of closure attempts.
	/// @details Perturbation, closure, and filtering is carried out for every closure
	/// attempt.  Successful closures from ALL attempts are then selected from by
	/// selectors.
	void set_closure_attempts( core::Size const attempts);


	/// @brief Returns the number of closure attempts.
	/// @details Perturbation, closure, and filtering is carried out for every closure
	/// attempt.  Successful closures from ALL attempts are then selected from by
	/// selectors.
	inline core::Size get_closure_attempts() const { return n_closure_attempts_; }

	/// @brief Sets number of tries before giving up.
	/// @details If this is set to 0, then no such check is made.
	/// The algorithm tries n_closure_attempts_ times if and only if at least one solution is found in the first
	/// ntries_before_giving_up_ attempts.
	void set_ntries_before_giving_up ( core::Size const ntries );

	/// @brief Gets number of tries before giving up.
	/// @details The algorithm tries n_closure_attempts_ times if and only if at least one solution is found in the first
	/// ntries_before_giving_up_ attempts.
	inline core::Size get_ntries_before_giving_up () const { return ntries_before_giving_up_; }


	/// @brief Gets the number of solutions that must be found before stopping prematurely (i.e.
	/// before n_closure_attempts_ is reached).
	/// @details If nonzero, attempts are made until the specified number of solutions are found OR
	/// n_closure_attempts_ is reached (or we give up because ntries_before_giving_up_ is reached
	/// without finding any solutions).  If zero, n_closure_attempts_ attempts are made, and then
	/// a solution is chosen from among the successful solutions (or we give up because
	/// ntries_before_giving_up_ is reached without finding any solutions.
	inline core::Size min_solution_count() const { return min_solution_count_; }


	/// @brief Sets the number of solutions that must be found before stopping prematurely (i.e.
	/// before n_closure_attempts_ is reached).
	/// @details If nonzero, attempts are made until the specified number of solutions are found OR
	/// n_closure_attempts_ is reached (or we give up because ntries_before_giving_up_ is reached
	/// without finding any solutions).  If zero, n_closure_attempts_ attempts are made, and then
	/// a solution is chosen from among the successful solutions (or we give up because
	/// ntries_before_giving_up_ is reached without finding any solutions.
	void set_min_solution_count(core::Size const count_in) { min_solution_count_ = count_in; return; }

	/// @brief Clear the stored solution poses.
	///
	void clear_stored_solutions();

	/// @brief Sets the mover that will be applied to all solutions that pass filters prior to applying the selector.
	///
	void set_preselection_mover ( protocols::moves::MoverOP mover );

	/// @brief Sets whether we are attaching a BOINC "ghost" pose observer.
	/// @details Only does anything in the BOINC graphics build.
	void set_attach_boinc_ghost_observer( bool const setting ) { attach_boinc_ghost_observer_ = setting; return; }

	/// @brief Gets whether we are attaching a BOINC "ghost" pose observer.
	/// @details Only does anything in the BOINC graphics build.
	inline bool attach_boinc_ghost_observer() const { return attach_boinc_ghost_observer_; }

	/// @brief Sets whether we're using low-memory mode.
	/// @details If false (the default) then a vector of result poses is stored.  This can use a lot of memory, though.
	/// If true, then only loop DOFs are stored.  This loses the results of applying any preselection movers, though.
	inline void set_low_memory_mode( bool const setting ) { low_memory_mode_=setting; return; }

	/// @brief Gets whether we're using low-memory mode.
	/// @details If false (the default) then a vector of result poses is stored.  This can use a lot of memory, though.
	/// If true, then only loop DOFs are stored.  This loses the results of applying any preselection movers, though.
	inline bool low_memory_mode() const { return low_memory_mode_; }

	/// @brief Sets whether the mover sets its status to failure of no solution is found.
	/// @details True means it does NOT fail if no solution is found.  False is the default (fail if no solution found).
	inline void set_dont_fail_if_no_solution_found( bool const setting ) { dont_fail_if_no_solution_found_=setting; return; }

	/// @brief Gets whether the mover sets its status to failure of no solution is found.
	/// @details True means it does NOT fail if no solution is found.  False is the default (fail if no solution found).
	inline bool dont_fail_if_no_solution_found() const { return dont_fail_if_no_solution_found_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Sets whether this mover corrects positions of polymer-dependent atoms.
	///
	inline void set_correct_polymer_dependent_atoms( bool const setting ) { correct_polymer_dependent_atoms_ = setting; }

	/// @brief Gets whether this mover corrects positions of polymer-dependent atoms.
	///
	inline bool correct_polymer_dependent_atoms() const { return correct_polymer_dependent_atoms_; }

private:

	/// @brief The list of residues (as inidices of the original pose) making up the loop to be closed.
	utility::vector1 < core::Size > loopresidues_;


	/// @brief The list of tail residues (as indices of the original pose) that will "come along for the ride" as the
	/// loop moves.  These must be attached to the original loop, either directly, or through other tail residues.
	utility::vector1 < core::Size > tailresidues_;

	/// @brief The lower end of the loop to be closed is presumably connected to geometry that is not
	/// moved by this mover.  However, the connection could be through backbone or sidechain.  The
	/// mover needs to know the connection ID on the lower end of the loop to the anchor geometry.
	/// Note: This will be set to 0 if the only connection that this residue can form is to the loop.
	core::Size lower_anchor_connID_;

	/// @brief The upper end of the loop to be closed is presumably connected to geometry that is not
	/// moved by this mover.  However, the connection could be through backbone or sidechain.  The
	/// mover needs to know the connection ID on the upper end of the loop to the anchor geometry.
	/// Note: This will be set to 0 if the only connection that this residue can form is to the loop.
	core::Size upper_anchor_connID_;


	/// @brief Should the mover build ideal geometry for the loop to be closed?  Default true.
	bool build_ideal_geometry_;


	/// @brief Was the last apply() successful in generating a new kinematically-closed structure?
	bool last_run_successful_;

	/// @brief List of the atoms in the segment to be closed, and their x,y,z coordinates.  Note that
	/// the atomIDs refer to residue indices that correspond to the temporary pose created by the
	/// apply() function, not to the indices in the original pose.
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > atomlist_;


	/// @brief Pivot residue 1 (index based on original pose).
	core::Size pivot_1_rsd_;


	/// @brief Pivot atom 1 name
	std::string pivot_1_atmname_;


	/// @brief Pivot residue 2 (index based on original pose).
	core::Size pivot_2_rsd_;


	/// @brief Pivot atom 2 name
	std::string pivot_2_atmname_;


	/// @brief Pivot residue 3 (index based on original pose).
	core::Size pivot_3_rsd_;


	/// @brief Pivot atom 3 name
	std::string pivot_3_atmname_;


	/// @brief List of GeneralizedKICperturbers to apply prior to kinematic closure.
	perturber::GeneralizedKICperturberOPs perturberlist_;


	/// @brief List of GeneralizedKICfilters to apply to eliminate bad solutions.
	filter::GeneralizedKICfilterOPs filterlist_;

	/// @brief The selection algorithm to use to pick a solution from the set of solutions passing
	/// the filters.
	selector::GeneralizedKICselectorOP selector_;


	/// @brief The number of closure attempts to carry out.
	/// @details Perturbation, closure, and filtering is carried out for every closure
	/// attempt.  Successful closures from ALL attempts are then selected from by
	/// selectors.  Zero means no limit.
	core::Size n_closure_attempts_;


	/// @brief If greater than zero, GenKIC will seek solutinons until it has at least
	/// this number of solutions.  If it is zero, then GenKIC will seek solutions until
	/// n_closure_attempts_ is reached, or until ntries_before_giving_up_ is reached
	/// without finding a solution.
	/// @details Zero by default.
	core::Size min_solution_count_;


	/// @brief Owning pointer for the RosettaScripts ContingentFilter associated with this mover, if there is one.
	protocols::filters::ContingentFilterOP rosettascripts_filter_;


	/// @brief Bool determining whether there exists a RosettaScripts ContingentFilter associated with this mover.
	bool rosettascripts_filter_exists_;


	/// @brief Owning pointer for a pre-selection mover applied to all solutions that will be passed to the selector.
	protocols::moves::MoverOP pre_selection_mover_;


	/// @brief Bool determining whether there exists a pre-selection mover that wlil be applied.
	bool pre_selection_mover_exists_;


	/// @brief The number of tries to make before giving up if no solution has been found yet.
	/// @details If this is set to 0, then no such check is made.
	/// The algorithm tries n_closure_attempts_ times if and only if at least one solution is found in the first
	/// ntries_before_giving_up_ attempts.
	core::Size ntries_before_giving_up_;

	/// @brief Vector of owning pointers to solutions.  Used for temporary storage.
	///
	utility::vector1 < core::pose::PoseOP > solutions_;

	/// @brief Should we attach a BOINC "ghost" pose observer for BOINC graphics output?
	/// @details Default false.  Only does anything in the BOINC build.
	bool attach_boinc_ghost_observer_;

	/// @brief Use low memory mode?  Default false.
	/// @details If false (the default) then a vector of result poses is stored.  This can use a lot of memory, though.
	/// If true, then only loop DOFs are stored.  This loses the results of applying any preselection movers, though.
	bool low_memory_mode_;

	/// @brief By default, the mover sets its status to failure if no solution is found.  If dont_fail_if_no_solution_found_ is true (not
	/// the default) then this doesn't happen.
	bool dont_fail_if_no_solution_found_;

	/// @brief Vector of loop bond angles from final solutions.
	/// @details Only stored in low-memory mode.
	utility::vector1 < utility::vector1 < core::Real > > bondangle_solutions_;

	/// @brief Vector of loop bond lengths from final solutions.
	/// @details Only stored in low-memory mode.
	utility::vector1 < utility::vector1 < core::Real > > bondlength_solutions_;

	/// @brief Should we correct positions of polymer bond-dependent atoms?
	/// @details Default false.
	bool correct_polymer_dependent_atoms_;

private:

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Gets a path through atoms in the residue from one connection atom to another within a residue.
	/// @details  This is not necessarily the shortest path, geometrically, but hopefully it will do for our purposes.  (I
	/// didn't want to reimplement Dijkstra's algorithm, since that would tax my programming abilities and would result in
	/// unnecessarily slow performance.)  This algorithm works by stepping back (child->parent) from the higher-index atom
	/// until (a) it finds the lower-index atom, or (b) it reaches the root.  If it reaches the root, it then starts at the
	/// lower-index atom and traces back until it reaches an atom already in the path.  The chosen path is therefore the path
	/// from one atom, to the nearest common ancestor, to the other atom.
	///
	/// @param[in] first_atom -- The index of the first atom in the path.
	/// @param[in] second_atom -- The index of the first atom in the path.
	/// @param[in] rsd -- The residue object (const instance).
	/// @param[out] path_indices -- The list of indices of the atoms found in the path.
	///
	/// @author Vikram K. Mulligan
	void get_path (
		core::Size const first_atom,
		core::Size const second_atom,
		core::conformation::Residue const &rsd,
		utility::vector1 <core::Size> &path_indices
	);


	/// @brief Function to get the FIRST connection in res_with_connection that connects it to other_res:
	core::Size get_connection(
		core::Size const res_with_connection,
		core::Size const other_res,
		core::pose::Pose const &pose
	);

	/// @brief Function that returns true if two residues in a pose have a direct geometric connection,
	/// and false otherwise.
	bool has_geometric_connection (
		core::Size const residue1,
		core::Size const residue2,
		core::pose::Pose const &pose
	);

	/// @brief As the list of residues in the loop to be closed is updated, we need to figure out how
	/// that loop is connected to the geometry outside of the loop (i.e. what's considred the connection
	/// to stationary geometry).  This function loops through all connIDs on the terminal residues of
	/// the loop to be closed and picks the first one that links to geometry not in the loop as the
	/// anchor connnection.  TODO: Add a manual override to specifiy that a different connection is the
	/// anchor.
	void infer_anchor_connIDs(core::pose::Pose const &pose);


	/// @brief Find the residue that is the anchor of the lower end of the loop that we're about to close and add it to the loop pose by a jump.
	void addloweranchor(
		core::pose::Pose &perturbedloop_pose,
		core::pose::Pose const &pose
	);


	/// @brief Add the loop geometry from the starting pose to the temporary pose used for kinematic closure.  This will build ideal geometry if build_ideal==true (NOT YET TESTED).
	void addloopgeometry(
		core::pose::Pose &perturbedloop_pose,
		core::pose::Pose const &pose,
		bool const build_ideal,
		//core::Size const effect_on_bonded_geom,
		utility::vector1 < std::pair < core::Size, core::Size > > &residue_map
	);

	/// @brief Add the tail geometry from the starting pose to the temporary pose used for kinematic closure.  This will build ideal geometry if build_ideal==true (NOT YET TESTED).
	/// @details The tails are residues that are not part of the loop to be closed, but which are attached to them and which "come along for the ride" as the loop to be closed moves.
	/// This function MUST be called AFTER addloopgeometry().
	void addtailgeometry(
		core::pose::Pose &perturbedloop_pose,
		core::pose::Pose const &pose,
		bool const build_ideal,
		//core::Size const effect_on_bonded_geom,
		utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map,
		utility::vector1 < std::pair < core::Size, core::Size > > &tail_residue_map
	);


	/// @brief Find the residue that is the anchor of the upper end of the loop that we're about to close and add it to the loop pose by a bond.
	void addupperanchor(
		core::pose::Pose &perturbedloop_pose,
		core::pose::Pose const &pose
	);

	/// @brief Do the actual kinematic closure.
	/// @details Inputs are pose (the loop to be closed), original_pose (the reference pose, unchanged by operation), residue_map (the mapping of
	/// residues from pose to original_pose) and tail_residue_map (the mapping of tail residues from pose to original_pose).  Output is the index
	/// of the solution in the solutions_ vector.
	/// @notes The selected_torsions, selected_bondangles, and selected_bondlengths vectors are output only in low-memory mode.
	bool doKIC(
		core::pose::Pose const &pose,
		core::pose::Pose const &original_pose,
		utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map,
		utility::vector1 < std::pair < core::Size, core::Size > > const &tail_residue_map,
		core::Size &solution_index,
		utility::vector1 < core::Real > & selected_torsions,
		utility::vector1 < core::Real > & selected_bondangles,
		utility::vector1 < core::Real > & selected_bondlengths
	);


	/// @brief Generate the list of atomIDs for the chain that will be closed by kinematic closure.
	void generate_atomlist(
		core::pose::Pose const &pose,
		utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map
	);


	/// @brief Generate the numeric::kinematic_closure::bridgeObjects data from the atomlist_ object.
	void generate_bridgeobjects_data_from_atomlist(
		utility::vector1< utility::fixedsizearray1< core::Real,3 > > &atoms, //atom xyz
		utility::vector1< core::Real > &dt, //desired torsions for each atom
		utility::vector1< core::Real > &da, //desired bond angle for each atom
		utility::vector1< core::Real > &db, //desired bond length for each atom
		utility::vector1< core::Size > &order //use 1, 2, 3
	);

	/// @brief Given a residue_map vector of pairs, where each pair is < residue_index_in_perturbedloop_pose, residue_index_in_original_pose >,
	/// and a residue index in the original pose, return the corresponding residue index in the perturbed loop pose.
	core::Size get_perturbedloop_rsd ( core::Size const original_pose_rsd, utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map );


	/// @brief Pick the pivots for kinematic closure.
	void pick_pivots(
		core::pose::Pose const &original_pose,
		utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
		utility::vector1 < core::Size > &pivots
	);

	/// @brief Apply the list of perturbers (everything in perturberlist_) to alter the desired torsion, desired angle, and desired bond length lists prior to calling bridgeObjects.
	/// @details
	///
	/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
	/// @param[in] original_pose -- A pose consisting of the full, original structure.
	/// @param[in] residue_map -- The mapping of (residue in loop_pose, residue in original_pose).
	/// @param[in] tail_residue_map -- The mapping of (tail residue in loop_pose, tail residue in original_pose).
	/// @param[in,out] torsions -- The desired torsion angles, potentially altered or overwritten by the perturbers.
	/// @param[in,out] bondangles -- The desired bond angles, potentially altered or overwritten by the perturbers.
	/// @param[in,out] bondlenghts -- The desired bond lengths, potentially altered or overwritten by the perturbers.
	void apply_perturbations(
		core::pose::Pose const &loop_pose,
		core::pose::Pose const &original_pose,
		utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map,
		utility::vector1 < std::pair < core::Size, core::Size > > const &tail_residue_map,
		utility::vector1 < core::Real > &torsions,
		utility::vector1 < core::Real > &bondangles,
		utility::vector1 < core::Real > &bondlengths
	) const;

	/// @brief Apply filters to the list of solutions, and proceed to eliminate any and all that fail to pass filters.
	/// @details This removes entries from the torsions, bondangles, and bondlengths lists, and decrements nsol (the number of solutions) appropriately.
	/// @param[in] original_pose -- A pose consisting of the full, original structure.
	/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
	/// @param[in] residue_map -- The mapping of (residue in loop_pose, residue in original_pose).
	/// @param[in] tail_residue_map -- The mapping of (tail residue in loop_pose, tail residue in original_pose).
	/// @param[in] atomlist -- The list of atomIDs in the chain that was closed, plus their xyz coordinates from the original pose.
	/// @param[in,out] torsions -- The torsion angles returned by bridgeObjects, as a matrix of [solution #][torsion index].  Columns can be deleted by filters.
	/// @param[in,out] bondangles -- The bond angles returned by bridgeObjects, as a matrix of [solution #][bondangle index].  Columns can be deleted by filters.
	/// @param[in,out] bondlenghts -- The bond lengths returned by bridgeObjects, as a matrix of [solution #][bondlength index].  Columns can be deleted by filters.
	/// @param[in,out] nsol -- The number of solutions, which can be decremented as filters delete solutions.  Note that this must be an int, not a core::Size.
	void filter_solutions(
		core::pose::Pose const &original_pose,
		core::pose::Pose const &loop_pose,
		utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
		utility::vector1 < std::pair <core::Size, core::Size> > const &tail_residue_map,
		utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 <utility::vector1<core::Real> > &torsions,
		utility::vector1 <utility::vector1<core::Real> > &bondangles,
		utility::vector1 <utility::vector1<core::Real> > &bondlengths,
		int &nsol
	) const;

	/// @brief Applies the selector to choose a solution.
	/// @details  If the selector could not select a solution (e.g. if the preselection mover returned failed status for every solution), this function returns 0;
	/// otherwise, returns the index of the solution in the solutions_ vector.
	/// @param[in,out] pose -- The loop to be closed.
	/// @param[in] original_pose -- The original pose.  Can be used for reference by selectors.
	/// @param[in] residue_map -- Mapping of (loop residue, original pose residue).
	/// @param[in] tail_residue_map -- Mapping of (tail residue index in pose, tail residue index in original_pose).
	/// @param[in] atomlist -- The list of (AtomID, original XYZ coordinates of atoms) representing the chain that was closed.
	/// @param[in] torsions -- Matrix of [closure attempt #][solution #][torsion #] with torsion values for each torsion angle in the chain.  A selector will pick one solution.
	/// @param[in] bondangles -- Matrix of [closure attempt #][solution #][angle #] with bond angle values for each bond angle in the chain.  A selector will pick one solution.
	/// @param[in] bondlengths -- Matrix of [closure attempt #][solution #][bondlength #] with bond length for each bond in the chain.  A selector will pick one solution.
	/// @param[in] nsol_for_attempt -- List of the number of solutions for each attempt.
	/// @param[in] total_solutions -- Total number of solutions found.
	/// @param[in] energies_for_solution -- Vectors of the energies from each attempt.  Used only in low-memory mode.
	core::Size select_solution (
		core::pose::Pose const &pose,
		core::pose::Pose const &original_pose, //The original pose
		utility::vector1 <std::pair <core::Size, core::Size> > const &residue_map, //mapping of (loop residue, original pose residue)
		utility::vector1 <std::pair <core::Size, core::Size> > const &tail_residue_map, //mapping of (tail residue index in pose, tail residue index in original_pose)
		utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, //torsions for each atom
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondangles, //bond angle for each atom
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondlengths, //bond length for each atom
		utility::vector1 <core::Size> const &nsol_for_attempt,
		core::Size const total_solutions,
		utility::vector1 <core::Real> const &energies_for_solution
	) const;

	/// @brief Trims extra atoms from the start and end of the atom list, if the first and last pivots are not the fifth and fifth-last atoms, respectively.
	///
	void prune_extra_atoms( utility::vector1 <core::Size> &pivots );

	/// @brief Returns whether a preselection mover has been specified.
	///
	bool preselection_mover_exists() const { return pre_selection_mover_exists_; }

	/// @brief Add a solution to the solutions list.
	/// @details Stores an owning pointer to a pose.  Not used in low-memory mode (see function overload).
	void add_solution( core::pose::PoseOP pose_in );

	/// @brief Sets the pose to a particular solution from the solutions list.
	/// @details This is for regular mode, so it just grabs this from the stored solution pose vector.
	void set_final_solution( core::pose::Pose &output_pose, core::Size const solution_index ) const;

	/// @brief Sets the pose to a particular solution from the solutions list.
	/// @details This is for low-memory mode, so it constructs the pose from DoF vectors and reapplies preselection movers.
	protocols::moves::MoverStatus
	set_final_solution(
		core::pose::Pose &pose,
		core::pose::Pose &looppose,
		utility::vector1 < std::pair < core::Size, core::Size > > const & residue_map,
		utility::vector1 < std::pair < core::Size, core::Size > > const & tail_residue_map,
		utility::vector1< core::Real > const &torsions,
		utility::vector1< core::Real > const &bondangles,
		utility::vector1< core::Real > const &bondlengths
	) const;

	/// @brief Return the number of stored solutions
	///
	inline core::Size total_stored_solutions() const { return solutions_.size(); }

};

} //namespace generalized_kinematic_closure
} //namespace protocols

#endif //INCLUDED_protocols_generalized_kinematic_closure_GeneralizedKIC_hh
