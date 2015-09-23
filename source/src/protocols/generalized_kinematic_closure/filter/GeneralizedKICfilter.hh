// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/generalized_kinematic_closure/filter/GeneralizedKICfilter.hh
/// @brief  Headers for GeneralizedKICfilter class (helper class for the GeneralizedKIC mover that defines filters to accept/reject closure solutions).
///         Filters must be pass/fail (i.e. there are no shades of grey, here.)
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_generalized_kinematic_closure_filter_GeneralizedKICfilter_hh
#define INCLUDED_protocols_generalized_kinematic_closure_filter_GeneralizedKICfilter_hh

// Unit Headers
//#include <protocols/moves/Mover.hh>
#include <protocols/generalized_kinematic_closure/filter/GeneralizedKICfilter.fwd.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

//// Project Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/bin_transitions/BinTransitionCalculator.fwd.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>

// How to add new filter types:
// 1. Add a new entry in the filter_type enum list.
// 2. Add the perturber effect name to the get_filter_type_name function.
// 3. Add the perturber effect to the switch() statement in apply().
// 4. Create an apply_<filter_type_name>() function as a private method of the GeneralizedKICfilter class.

///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace generalized_kinematic_closure {
namespace filter {

enum filter_type {
	//When adding filter types, add the names to the get_filter_type_name() function.

	no_filter = 1,
	loop_bump_check,
	atom_pair_distance,
	backbone_bin,
	alpha_aa_rama_check,

	unknown_filter, //Keep this second-to-last.
	end_of_filter_list = unknown_filter //Keep this last.

};


///////////////////////////////////////////////////////////////////////
//   GENERALIZED KIC FILTER CLASS                                    //
///////////////////////////////////////////////////////////////////////

class GeneralizedKICfilter : public utility::pointer::ReferenceCount
{
public:
	GeneralizedKICfilter();
	GeneralizedKICfilter(GeneralizedKICfilter const &src);
	~GeneralizedKICfilter();
	GeneralizedKICfilterOP clone() const;


	/// @brief Returns the name of this class.
	std::string get_name() const;


	/// @brief Given a filter type, return its name.  Returns "unknown_filter" if not recognized.
	std::string get_filter_type_name( core::Size const filter_type ) const;

	///
	/// @brief Given the name of a filter type, return the filter type enum.  Returns unknown_filter if not recognized.
	filter_type get_filter_type_by_name( std::string const &filtername ) const;


	/// @brief Sets the filter type for this filter.
	void set_filter_type( filter_type const &ftype);


	/// @brief Sets the filter type for this filter by name.
	void set_filter_type( std::string const &ftypename);


	/// @brief Gets the filter type name for THIS filter.
	std::string get_this_filter_type_name () const;


	/// @brief Add a real-valued filter parameter.
	void add_filter_param( std::string const &param_name, core::Real const &value );


	/// @brief Add a integer-valued filter parameter.
	void add_filter_param( std::string const &param_name, core::Size const value );


	/// @brief Add a Boolean-valued filter parameter.
	void add_filter_param( std::string const &param_name, bool const value );


	/// @brief Add a string-valued filter parameter.
	void add_filter_param( std::string const &param_name, std::string const &value );

	/// @brief Get a real-valued filter parameter.
	/// @details Returns false if the parameter couldn't be found.
	bool get_filter_param( std::string const &param_name, core::Real &outvalue ) const;

	/// @brief Get a integer-valued filter parameter.
	/// @details Returns false if the parameter couldn't be found.
	bool get_filter_param( std::string const &param_name, core::Size &outvalue ) const;

	/// @brief Get a Boolean-valued filter parameter.
	/// @details Returns false if the parameter couldn't be found.
	bool get_filter_param( std::string const &param_name, bool &outvalue ) const;

	/// @brief Get a string-valued filter parameter.
	/// @details Returns false if the parameter couldn't be found.
	bool get_filter_param( std::string const &param_name, std::string &outvalue ) const;

	/// @brief Set the residue number that this filter acts on.
	/// @details Only used by some filters.
	inline void set_resnum( core::Size const val ) { resnum_ = val; return; }

	/// @brief Get the residue number that this filter acts on.
	/// @details Only used by some filters.
	inline core::Size resnum( ) const { return resnum_; }

	/// @brief Set the bin name for this filter.
	/// @details Only used by some filters.
	inline void set_binname( std::string const &name_in ) {bin_ = name_in; return; }

	/// @brief Set the alpha_aa_rama_check filter's rama term cutoff energy.
	/// @details Only used by alpha_aa_rama_check filter.
	inline void set_rama_cutoff_energy( core::Real const &val ) { rama_threshhold_=val; return; }

	/// @brief Get the alpha_aa_rama_check filter's rama term cutoff energy.
	/// @details Only used by alpha_aa_rama_check filter.
	inline core::Real rama_cutoff_energy( ) const { return rama_threshhold_; }

	/// @brief Get the bin name for this filter.
	/// @details Only used by some filters.
	inline std::string binname( ) const { return bin_; }

	/// @brief Initializes the BinTransitionCalculator object and loads a bin_params file.
	///
	void load_bin_params( std::string const &bin_params_file );

	/// @brief Apply this filter to ONE of the kinematic closure solutions produced by the bridgeObjects function,
	/// and return pass or fail.
	/// @details
	/// @param[in] original_pose -- The full, initial pose.
	/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
	/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
	/// @param[in] tail_residue_map -- The mapping of (tail residue index in loop_pose, tail residue index in original_pose).
	/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
	/// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
	/// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
	/// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
	bool apply(
		core::pose::Pose const &original_pose,
		core::pose::Pose const &loop_pose,
		utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
		utility::vector1 < std::pair <core::Size, core::Size> > const &tail_residue_map,
		utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 < core::Real > const &torsions,
		utility::vector1 < core::Real > const &bondangles,
		utility::vector1 < core::Real > const &bondlengths
	) const;


private:
	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE VARIABLES                                                 //
	////////////////////////////////////////////////////////////////////////////////


	/// @brief The filter type for this filter (see the filter_type enum for all types).
	filter_type filtertype_;


	/// @brief Real-valued filter parameters
	utility::vector1 < std::pair<std::string, core::Real > > filter_params_real_;


	/// @brief Integer-valued filter parameters
	utility::vector1 < std::pair<std::string, core::Size > > filter_params_size_;


	/// @brief Boolean-valued filter parameters
	utility::vector1 < std::pair<std::string, bool > > filter_params_bool_;


	/// @brief String-valued filter parameters
	utility::vector1 < std::pair<std::string, std::string > > filter_params_string_;

	/// @brief A BinTransitionCalculatorOP.  This will be null by default,
	/// and will only point to a BinTransitionCalculator object in the case
	/// of those filters that use torsion bin transition probabilities.
	core::scoring::bin_transitions::BinTransitionCalculatorOP bin_transition_calculator_;

	/// @brief A parameter specifically for the backbone_bin filter.  The bin
	/// that the residue must lie within.
	std::string bin_;

	/// @brief A parameter specifically for the backbone_bin filter.  The residue
	/// that must lie within the mainchain torsion bin specified.
	core::Size resnum_;

	/// @brief A parameter specifically for the alpha_aa_rama_check filter.  Rama energy
	/// above which the solution is rejected.  Set to 0.3 by default.
	core::Real rama_threshhold_;

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Given an index in the original pose and a mapping from loop to pose,
	/// return the index in the loop.
	core::Size get_loop_index (
		core::Size const original_pose_index,
		utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map
	) const;

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE APPLY FUNCTIONS FOR EACH FILTER                           //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Applies the loop_bump_check filter, which checks for clashes between the atoms in the chain
	/// to be closed and the rest of the structure (or for clashes within these atoms).
	/// @details Returns "true" for pass and "false" for fail.  Does NOT check for clashes with tail residues.
	/// @param[in] original_pose -- The full, initial pose.
	/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
	/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
	/// @param[in] tail_residue_map -- The mapping of (tail residue index in loop_pose, tail residue index in original_pose).
	/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
	/// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
	/// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
	/// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
	bool apply_loop_bump_check(
		core::pose::Pose const &original_pose,
		core::pose::Pose const &loop_pose,
		utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
		utility::vector1 < std::pair <core::Size, core::Size> > const &tail_residue_map,
		utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 < core::Real > const &torsions,
		utility::vector1 < core::Real > const &bondangles,
		utility::vector1 < core::Real > const &bondlengths
	) const;

	/// @brief Applies the atom_pair_distance filter, checking that the distance between two atoms is less than
	/// a given threshhold (or greater than a given threshhold if the user so specifies with the "greater_than"
	/// option).
	/// @details Returns "true" for pass and "false" for fail.  The user can set the following options:
	/// "distance" (real-valued, mandatory)
	/// "atom1" (string-valued, mandatory)
	/// "atom2" (string-valued, mandatory)
	/// "res1" (integer-valued, mandatory, based on original pose numbering)
	/// "res2" (integer-valued, mandatory, based on original pose numbering)
	/// "greater_than" (boolean, optional, false by default)
	/// @param[in] original_pose -- The full, initial pose.
	/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
	/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
	/// @param[in] tail_residue_map -- The mapping of (tail residue index in loop_pose, tail residue index in original_pose).
	/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
	/// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
	/// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
	/// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
	bool apply_atom_pair_distance(
		core::pose::Pose const &original_pose,
		core::pose::Pose const &loop_pose,
		utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
		utility::vector1 < std::pair <core::Size, core::Size> > const &tail_residue_map,
		utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 < core::Real > const &torsions,
		utility::vector1 < core::Real > const &bondangles,
		utility::vector1 < core::Real > const &bondlengths
	) const;

	/// @brief Applies the backbone_bin filter, checking that a given residue lies within a defined
	/// mainchain torsion bin and failing if it does not.
	/// @details Returns "true" for pass and "false" for fail.  The user needs to have set a bin
	/// transition probabilities file, a bin, and a residue.
	/// @param[in] original_pose -- The full, initial pose.
	/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
	/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
	/// @param[in] tail_residue_map -- The mapping of (tail residue index in loop_pose, tail residue index in original_pose).
	/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
	/// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
	/// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
	/// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
	bool apply_backbone_bin(
		core::pose::Pose const &original_pose,
		core::pose::Pose const &loop_pose,
		utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
		utility::vector1 < std::pair <core::Size, core::Size> > const &tail_residue_map,
		utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 < core::Real > const &torsions,
		utility::vector1 < core::Real > const &bondangles,
		utility::vector1 < core::Real > const &bondlengths
	) const;

	/// @brief Calculates Ramachandran energy for an alpha-amino acid based on its phi/psi values.
	/// @details Returns "true" for pass (below threshhold) and "false" for fail.
	/// @param[in] original_pose -- The full, initial pose.
	/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
	/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
	/// @param[in] tail_residue_map -- The mapping of (tail residue index in loop_pose, tail residue index in original_pose).
	/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
	/// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
	/// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
	/// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
	bool apply_alpha_aa_rama_check(
		core::pose::Pose const &original_pose,
		core::pose::Pose const &loop_pose,
		utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
		utility::vector1 < std::pair <core::Size, core::Size> > const &tail_residue_map,
		utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 < core::Real > const &torsions,
		utility::vector1 < core::Real > const &bondangles,
		utility::vector1 < core::Real > const &bondlengths
	) const;

}; //GeneralizedKICfilter class

} //namespace filter
} //namespace generalized_kinematic_closure
} //namespace protocols

#endif //INCLUDED_protocols_generalized_kinematic_closure_filter_GeneralizedKICfilter_hh
