// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/protocols/generalized_kinematic_closure/filter/GeneralizedKICfilter.hh
/// @brief  Headers for GeneralizedKICfilter class (helper class for the GeneralizedKIC mover that defines filters to accept/reject closure solutions).
///	        Filters must be pass/fail (i.e. there are no shades of grey, here.)
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
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/grid/CartGrid.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

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
	~GeneralizedKICfilter();

	///
	/// @brief Returns the name of this class.
	std::string get_name() const;

	///
	/// @brief Given a filter type, return its name.  Returns "unknown_filter" if not recognized.
	std::string get_filter_type_name( core::Size const filter_type ) const;
	
	///
	/// @brief Given the name of a filter type, return the filter type enum.  Returns unknown_filter if not recognized.
	filter_type get_filter_type_by_name( std::string const &filtername ) const;

	///
	/// @brief Sets the filter type for this filter.
	void set_filter_type( filter_type const &ftype);

	///
	/// @brief Sets the filter type for this filter by name.
	void set_filter_type( std::string const &ftypename);

	///
	/// @brief Gets the filter type name for THIS filter.
	std::string get_this_filter_type_name () const;

	/// @brief Apply this filter to ONE of the kinematic closure solutions produced by the bridgeObjects function,
	/// and return pass or fail.
	/// @details
	/// @param[in] original_pose -- The full, initial pose.
  /// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
	/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
	/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
  /// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
  /// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
  /// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
	bool apply(
		core::pose::Pose const &original_pose,
		core::pose::Pose const &loop_pose,
		utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
		utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 < core::Real > const &torsions,
		utility::vector1 < core::Real > const &bondangles,
		utility::vector1 < core::Real > const &bondlengths
	) const;


private:
////////////////////////////////////////////////////////////////////////////////
//          PRIVATE VARIABLES                                                 //
////////////////////////////////////////////////////////////////////////////////

	///
	/// @brief The filter type for this filter (see the filter_type enum for all types).
	filter_type filtertype_;

	utility::vector1 < std::pair<std::string, core::Real > > filter_params_real_;

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE APPLY FUNCTIONS FOR EACH FILTER                           //
////////////////////////////////////////////////////////////////////////////////

	/// @brief Applies the loop_bump_check filter, which checks for clashes between the atoms in the chain
	/// to be closed and the rest of the structure (or for clashes within these atoms).
	/// @details Returns "true" for pass and "false" for fail.
	/// @param[in] original_pose -- The full, initial pose.
  /// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
	/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
	/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
  /// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
  /// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
  /// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
	bool apply_loop_bump_check(
		core::pose::Pose const &original_pose,
		core::pose::Pose const &loop_pose,
		utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
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
