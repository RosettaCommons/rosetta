// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/FitSimpleHelix.hh
/// @brief  Headers for FitSimpleHelix mover class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_helical_bundle_FitSimpleHelix_hh
#define INCLUDED_protocols_helical_bundle_FitSimpleHelix_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/helical_bundle/FitSimpleHelix.fwd.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/ContingentFilter.fwd.hh>
#include <protocols/filters/ContingentFilter.hh>

// Project Headers
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/grid/CartGrid.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>

#include <set>

#include <core/grid/CartGrid.fwd.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace helical_bundle {

class FitSimpleHelix : public protocols::moves::Mover
{
public:
	FitSimpleHelix();
	virtual ~FitSimpleHelix();

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	///
	/// @brief Actually apply the mover to the pose.
	virtual void apply(core::pose::Pose & pose);

	virtual std::string get_name() const;

	/*virtual void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const &
	);*/

	void set_initial_guesses (core::Real const &r1_initial, core::Real const &omega1_initial, core::Real const &dz1_initial) {
		r1_initial_ = r1_initial;
		omega1_initial_=omega1_initial;
		dz1_initial_=dz1_initial;
		return;
	}

	void set_range (core::Size const start, core::Size const end) {
		start_index_ = start;
		end_index_ = end;
		runtime_assert_string_msg( start + 1 < end, "In FitSimpleHelix::set_range(): the end residue must be at least two residues after the start residue." );
		return;
	}

	///@brief Set the minimization type (e.g. dfpmin, linmin, etc.)
	///@details Defaults to dfpmin if not set.
	void set_min_type ( std::string const &min_type)
	{
		min_type_=min_type;
		return;
	}

	///@brief Set the minimizer tolerance (defaults to 1E-7, the default for many other protocols).
	///
	void set_min_tolerance (core::Real const &min_tol_in)
	{
		min_tolerance_ = min_tol_in;
		return;
	}

	/// @brief Set the mainchain atom that will be fit first, and used as the reference for other mainchain atoms.
	/// @details If not set, this defaults to "CA".
	void set_reference_atom (std::string const &ref_atom_in)
	{
		reference_atom_ = ref_atom_in;
		return;
	}

	/// @brief Function to retrieve the final values of the helical parameters, post-fit.
	/// @details Call this function after the "apply" function, and pass it containers for
	/// the data to be retrieved.
	void get_crick_parameters (
		utility::vector1 < core::Real > &r1_out,
		core::Real &omega1_out,
		core::Real &z1_out,
		utility::vector1 < core::Real > &delta_omega1_out,
		utility::vector1 < core::Real > &delta_z1_out		
	) const;

private:
////////////////////////////////////////////////////////////////////////////////
//          PRIVATE DATA                                                      //
////////////////////////////////////////////////////////////////////////////////

	/// @brief Initial guess for r1.
	///
	core::Real r1_initial_;

	/// @brief Initial guess for omega1.
	///
	core::Real omega1_initial_;

	/// @brief Initial guess for dz1, the rise per residue.
	///
	core::Real dz1_initial_;

	/// @brief Index of first residue in helix.
	///
	core::Size start_index_;

	/// @brief Index of last residue in helix.
	///
	core::Size end_index_;

	/// @brief The type of minimization that will be used.
	///
	std::string min_type_;

	/// @brief The minimizer tolerance.
	///
	core::Real min_tolerance_;

	/// @brief The mainchain torsion atom to be used as the reference atom.
	/// @details The fitter will first fit r1, omega1, and dz1 (the radius, turn per
	/// residue, and rise per residue) of the reference atom.  For all other mainchain
	/// atoms, the fitter will then solve for r1' and the omega and z offsets, keeping
	/// omega and dz1 constant across all mainchain atoms.
	std::string reference_atom_;

	/// @brief Output radii values (one for each atom fitted by the fitter).
	///
	utility::vector1 < core::Real > r1_vals_output_;

	/// @brief Output omega1 value (turn per residue, shared by all atoms fitted by the fitter).
	///
	core::Real omega1_val_output_;

	/// @brief Output z1 value (rise per residue, shared by all atoms fitted by the fitter).
	///
	core::Real z1_val_output_;

	/// @brief Output delta_omega1 values (radial offset, one for each atom fitted by the fitter).
	///
	utility::vector1 < core::Real > delta_omega1_vals_output_;

	/// @brief Output delta_z1 values (z offset, one for each atom fitted by the fitter).
	///
	utility::vector1 < core::Real > delta_z1_vals_output_;


////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

};

} //namespace helical_bundle
} //namespace protocols

#endif //INCLUDED_protocols_helical_bundle_FitSimpleHelix_hh
