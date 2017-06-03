// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.hh
/// @brief A base class for helper objects that the CrosslinkerMover uses to set up specific types
/// of threefold linkers.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_crosslinker_CrosslinkerMoverHelper_hh
#define INCLUDED_protocols_cyclic_peptide_crosslinker_CrosslinkerMoverHelper_hh

// Unit headers
#include <protocols/cyclic_peptide/crosslinker/CrosslinkerMoverHelper.fwd.hh>

// Protocol headers

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

/// @brief A base class for helper objects that the CrosslinkerMover uses to set up specific types
/// of threefold linkers.
class CrosslinkerMoverHelper : public utility::pointer::ReferenceCount {

public: //Constructors

	/// @brief Default constructor
	CrosslinkerMoverHelper();

	/// @brief Copy constructor
	CrosslinkerMoverHelper( CrosslinkerMoverHelper const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~CrosslinkerMoverHelper();


public: // public pure virtual methods

	/// @brief Given a pose and a selection of exactly three residues, add the linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	/// @details Must be defined by derived classes.  Version for asymmetric poses.
	virtual void add_linker_asymmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const = 0;

	/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
	/// @details Can be called by add_linker_asymmetric().  Must be defined by derived classes (pure virtual).  Version for asymmetric poses.
	virtual void add_linker_bonds_asymmetric(core::pose::Pose &pose, core::Size const res1, core::Size const res2, core::Size const res3, core::Size const linker_index ) const = 0;

	/// @brief Given a pose and a selection of exactly three residues, add the linker,
	/// align it crudely to the selected residues, and set up covalent bonds.
	/// @details Must be defined by derived classes.  Version for symmetric poses.
	virtual void add_linker_symmetric(core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const = 0;

	/// @brief Given a pose and a linker, add bonds between the linker and the residues that coordinate the linker.
	/// @details Can be called by add_linker_symmetric().  Must be defined by derived classes (pure virtual).  Version for symmetric poses.
	virtual void add_linker_bonds_symmetric(core::pose::Pose &pose, core::Size const res1, core::Size const linker_index1, core::Size const linker_index2 ) const = 0;

	/// @brief Given a selection of exactly three residues that have already been connected to a crosslinker,
	/// add constraints for the crosslinker.
	/// @details Must be defined by derived classes.  Version for asymmetric poses.
	virtual void add_linker_constraints_asymmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection) const = 0;

	/// @brief Given a selection of exactly three residues that have already been connected to a crosslinker,
	/// add constraints for the crosslinker.
	/// @details Must be defined by derived classes.  Version for symmetric poses.
	virtual void add_linker_constraints_symmetric( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const & selection, bool const linker_was_added) const = 0;

	/// @brief Given indices of three residues that are already linked to a linker, get the index
	/// of the linker.
	/// @details Throws an error if the three residues are not all linked to the same linker residue.  Must be defined by derived classes.
	virtual core::Size get_linker_index_asymmetric( core::pose::Pose const &pose, core::Size const res1, core::Size const res2, core::Size const res3 ) const = 0;

	/// @brief Given indices of three cysteine residues that are already linked to pieces of a linker, get
	/// of the indices of the symmetric pieces of the linker.
	/// @details Throws an error if a residue is not linked to something.  Must be defined by derived classes.
	virtual void get_linker_indices_symmetric( core::pose::Pose const &pose, core::Size const res1, core::Size const res2, core::Size const res3, core::Size &linker_index1, core::Size &linker_index2, core::Size &linker_index3 ) const = 0;

	/// @brief Given a pose with residues selected to be linked by a linker, determine whether the residues are too far apart.
	/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	virtual bool filter_by_sidechain_distance_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const & filter_multiplier ) const = 0;

	/// @brief Given a pose with residues selected to be linked by a linker, determine whether the residues are too far apart.
	/// @details Returns TRUE for failure (residues too far apart) and FALSE for success.  This version is for symmetric poses.
	/// @note Higher values of the filter multiplier make it more permissive.
	virtual bool filter_by_sidechain_distance_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const & filter_multiplier ) const = 0;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	virtual bool filter_by_constraints_energy_asymmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, core::Real const & filter_multiplier ) const = 0;

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  This version is for symmetric poses.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	/// @note Higher values of the filter multiplier make it more permissive.
	virtual bool filter_by_constraints_energy_symmetric( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, bool const linker_was_added, core::Real const & filter_multiplier ) const = 0;

public: // public defined methods

	/// @brief Given a ResidueSubset with exactly three residues selected, pull out the three indices.
	/// @details Overwrites res1, res2, and res3.
	virtual void get_sidechain_indices( core::select::residue_selector::ResidueSubset const & selection, core::Size &res1, core::Size &res2, core::Size &res3 ) const;

	/// @brief Determine whether a selection is symmetric.
	/// @details Returns true if and only if (a) the pose is symmetric, (b) there are three symmetry copies, and (c) the selected residues are equivalent residues in different
	/// symmetry copies.  Note that, ideally, I'd like to test for c3 symmetry, but this is as close as was feasible.
	/// @note Can be overriden.
	virtual bool selection_is_symmetric( core::select::residue_selector::ResidueSubset const & selection, core::pose::Pose const &pose ) const;

	/// @brief Optional steps that the helper can apply before every relaxation round.
	/// @details Defaults to doing nothing; can be overriden.  (One example is the TMA helper, which uses this to update amide bond-dependent atom positions).
	virtual void pre_relax_round_update_steps( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const &selection, bool const whole_structure, bool const symmetric, bool const linker_was_added ) const;

	/// @brief Optional steps that the helper can apply after every relaxation round.
	/// @details Defaults to doing nothing; can be overriden.  (One example is the TMA helper, which uses this to update amide bond-dependent atom positions).
	virtual void post_relax_round_update_steps( core::pose::Pose &pose, core::select::residue_selector::ResidueSubset const &selection, bool const whole_structure, bool const symmetric, bool const linker_was_added ) const;

protected: // protected methods

	/// @brief Determine whether the sidechain-crosslinker system has too high a constraints score.  Both the symmetric and asymmetric versions call this code.
	/// @details Returns TRUE for failure (too high a constraints score) and FALSE for success.
	virtual bool filter_by_constraints_energy( core::pose::Pose const &pose, core::select::residue_selector::ResidueSubset const & selection, bool const symmetric, bool const linker_was_added, core::Real const &filter_multiplier ) const;


private: // private methods


private: // data



};

} //crosslinker
} //protocols
} //cyclic_peptide

#endif //protocols/cyclic_peptide_crosslinker_CrosslinkerMoverHelper_hh
