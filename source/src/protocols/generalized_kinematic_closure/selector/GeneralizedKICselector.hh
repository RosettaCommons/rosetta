// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/generalized_kinematic_closure/selector/GeneralizedKICselector.hh
/// @brief  Headers for GeneralizedKICselector class (helper class for the GeneralizedKIC mover that defines how closure
/// solutions that pass filters are chosen).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_generalized_kinematic_closure_selector_GeneralizedKICselector_hh
#define INCLUDED_protocols_generalized_kinematic_closure_selector_GeneralizedKICselector_hh

// Unit Headers
//#include <protocols/moves/Mover.hh>
#include <protocols/generalized_kinematic_closure/selector/GeneralizedKICselector.fwd.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

//// Project Headers
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/grid/CartGrid.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>

// How to add new selector types:
// 1. Add a new entry in the selector_type enum list.
// 2. Add the selector name to the get_selector_type_name function.
// 3. Add the selector to the switch() statement in apply().
// 4. Create an apply_<selector_type_name>() function as a private method of the GeneralizedKICselector class.

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace generalized_kinematic_closure {
namespace selector {

enum selector_type {
	//When adding selector types, add the names to the get_selector_type_name() function.

	no_selector = 1,
	random_selector,
	lowest_energy_selector, //Uses whatever scorefunction is passed to the selector
	boltzmann_energy_selector, //Randomly picks weighted by exp(-E/kbt).
	lowest_rmsd_selector, //Picks the loop conformation closest to the original
	lowest_delta_torsion_selector, //Picks the loop conformation closest to the original (in torsion space), avoid flip at ends

	unknown_selector, //Keep this second-to-last.
	end_of_selector_list = unknown_selector //Keep this last.

};


///////////////////////////////////////////////////////////////////////
//   GENERALIZED KIC SELECTOR CLASS                                  //
///////////////////////////////////////////////////////////////////////

class GeneralizedKICselector : public utility::pointer::ReferenceCount
{
public:
	GeneralizedKICselector();
	~GeneralizedKICselector();

	///
	/// @brief Returns the name of this class.
	std::string get_name() const;

	///
	/// @brief Given a selector type, return its name.  Returns "unknown_selector" if not recognized.
	std::string get_selector_type_name( core::Size const selector_type ) const;
	
	///
	/// @brief Given the name of a selector type, return the selector type enum.  Returns unknown_selector if not recognized.
	selector_type get_selector_type_by_name( std::string const &selectorname ) const;

	///
	/// @brief Sets the selector type for this selector.
	void set_selector_type( selector_type const &stype);

	///
	/// @brief Sets the selector type for this selector by name.
	void set_selector_type( std::string const &stypename);

	///
	/// @brief Returns the selector type for this selector.
	selector_type get_selector_type () const { return selectortype_; }

	///
	/// @brief Set the scorefunction used by this selector.
	void set_scorefunction( core::scoring::ScoreFunctionOP sfxn ) { selector_sfxn_=sfxn; return; }

	///
	/// @brief Set the Boltzmann temperature used by this selector.
	void set_boltzmann_temp( core::Real const &temp) { boltzmann_kbt_=temp; return; }

	///
	/// @brief Returns the Boltzmann temperature used by this selector.
	core::Real get_boltzmann_temp() const { return boltzmann_kbt_; }

	/// @brief Applies a selector type to choose a solution and set a loop pose.
	/// @details
	/// @param[in,out] pose -- The loop to be closed.  This function puts it into its new, closed conformation.
	/// @param[in] original_pose -- The original pose.  Can be used for reference by selectors.
	/// @param[in] residue_map -- Mapping of (loop residue, original pose residue).
	/// @param[in] tail_residue_map -- Mapping of (tail residue index in pose, tail residue index in original_pose).
	/// @param[in] atomlist -- The list of (AtomID, original XYZ coordinates of atoms) representing the chain that was closed.
	/// @param[in] torsions -- Matrix of [closure attempt #][solution #][torsion #] with torsion values for each torsion angle in the chain.  A selector will pick one solution.
	/// @param[in] bondangles -- Matrix of [closure attempt #][solution #][angle #] with bond angle values for each bond angle in the chain.  A selector will pick one solution.
	/// @param[in] bondlengths -- Matrix of [closure attempt #][solution #][bondlength #] with bond length for each bond in the chain.  A selector will pick one solution.
	/// @param[in] nsol_for_attempt -- List of the number of solutions for each attempt.
	/// @param[in] total_solutions -- Total number of solutions found.
	/// @param[in] pre_selectoin_mover -- Pointer to a mover applied to each solution before applying the selector.
	/// @param[in] preselection_mover_exists -- Boolean that determines whether a mover has been specified.
	void apply (
		core::pose::Pose &pose,
		core::pose::Pose const &original_pose, //The original pose
		utility::vector1 <std::pair <core::Size, core::Size> > const &residue_map, //mapping of (loop residue, original pose residue)
		utility::vector1 <std::pair <core::Size, core::Size> > const &tail_residue_map, //mapping of (tail residue index in pose, tail residue index in original_pose)
		utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, //torsions for each atom 
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondangles, //bond angle for each atom
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondlengths, //bond length for each atom
		utility::vector1 <core::Size> const &nsol_for_attempt,
		core::Size const total_solutions,
		protocols::moves::MoverOP pre_selection_mover,
		bool const preselection_mover_exists
	) const;


private:
////////////////////////////////////////////////////////////////////////////////
//          PRIVATE VARIABLES                                                 //
////////////////////////////////////////////////////////////////////////////////

	///
	/// @brief The selector type for this selector (see the selector_type enum for all types).
	selector_type selectortype_;

	///
	/// @brief An owning pointer to a scoring function that a selector can use.
	/// @details This must be set explicitly; otherwise, it's set to NULL by default.
	core::scoring::ScoreFunctionOP selector_sfxn_;

	///
	/// @brief A Boltzmann temperature (kbt, in Rosetta energy units) that some selectors can use.
	core::Real boltzmann_kbt_;

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE APPLY FUNCTIONS FOR EACH FILTER                           //
////////////////////////////////////////////////////////////////////////////////

	/// @brief Applies a random_selector selector.
	/// @details This picks a solution randomly from the solutions that passed filters.
	void apply_random_selector(
		utility::vector1<core::Size> const &nsol_for_attempt,
		core::Size const total_solutions,
		core::Size &chosen_attempt_number,
		core::Size &chosen_solution
	) const;

	/// @brief Applies a lowest_energy_selector selector.
	/// @details This picks the lowest-energy solution, as scored with sfxn.  It's a good idea to use a modified
	/// scorefunction for this (something that just has the backbone conformation and H-bonding terms, for
	/// example, since side-chains will not be repacked by default prior to invoking this selector).
	void apply_lowest_energy_selector(
		utility::vector1<core::Size> const &nsol_for_attempt,
		core::Size const total_solutions,
		core::Size &chosen_attempt_number,
		core::Size &chosen_solution,
		core::scoring::ScoreFunctionOP sfxn,
		utility::vector1 <std::pair <core::Size, core::Size> > const &residue_map,
		utility::vector1 <std::pair <core::Size, core::Size> > const &tail_residue_map,
		utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, 
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondangles, 
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondlengths, 
		core::pose::Pose const &ref_loop_pose,
		core::pose::Pose const &ref_pose,
		core::Real const &boltzmann_kbt,
		protocols::moves::MoverOP pre_selection_mover,
		bool const preselection_mover_exists, 
		bool const use_boltzmann
	) const;

	/// @brief Applies a lowest_rmsd_selector selector.
	/// @details This picks the solution with the lowest RMSD from the starting pose.
	void apply_lowest_rmsd_selector( 
		utility::vector1<core::Size> const &nsol_for_attempt,
		core::Size const total_solutions,
		core::Size &chosen_attempt_number,
		core::Size &chosen_solution,
		utility::vector1 <std::pair <core::Size, core::Size> > const &residue_map,
		utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, 
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondangles, 
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondlengths, 
		core::pose::Pose const &ref_loop_pose,
		bool const preselection_mover_exists 
	) const;

	/// @brief Applies a lowest_delta_torsion_selector.
	/// @details This picks the solution with the lowest RMSD from the starting pose.
	void apply_lowest_delta_torsion_selector(
		utility::vector1<core::Size> const &nsol_for_attempt,
		core::Size const total_solutions,
		core::Size &chosen_attempt_number,
		core::Size &chosen_solution,
		utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, 
		core::pose::Pose const &pose,
		bool const preselection_mover_exists
	) const;

}; //GeneralizedKICselector class

} //namespace selector
} //namespace generalized_kinematic_closure
} //namespace protocols

#endif //INCLUDED_protocols_generalized_kinematic_closure_selector_GeneralizedKICselector_hh
