// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/SetTorsion.hh
/// @brief Sets the value of a desired torsion.  Header files for the mover.
/// @author Modified 4 June 2015 by Vikram K. Mulligan (vmullig@uw.edu), Baker Laboratory,
/// to add perturb torsion option (I didn't write this file, though).

#ifndef INCLUDED_protocols_simple_moves_SetTorsion_hh
#define INCLUDED_protocols_simple_moves_SetTorsion_hh

#include <protocols/simple_moves/SetTorsion.fwd.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/id/NamedAtomID.hh>

//parsing
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map

#include <utility/vector1.hh>


// Utility headers

// C++ headers

// Unit headers

namespace protocols {
namespace simple_moves {

enum TorsionPerturbType {
	perturbtorsion_uniform=1,
	perturbtorsion_gaussian,
	perturbtorsion_unknown //KEEP THIS LAST
};

/// @brief A mover to change one torsion angle
class SetTorsion : public protocols::moves::Mover
{
private:
	typedef protocols::moves::Mover parent;
public:
	/// @brief default ctor
	SetTorsion();
	virtual ~SetTorsion();

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const {
		return (protocols::moves::MoverOP( new protocols::simple_moves::SetTorsion( *this ) ) );
	}
	virtual protocols::moves::MoverOP fresh_instance() const {
		return protocols::moves::MoverOP( new SetTorsion );
	}

	core::Size n_torsion_sets() const {
		return residues_.size();
	}

	utility::vector1<core::Size> residue_list(core::Size iset, core::pose::Pose const & pose);

	/// @brief Actually get the value that the torsion will be set to.
	/// @details Depending on settings, this will look up a value, generate a random value, or perturb an input value.
	core::Real angle(
		core::Size const iset,
		core::Real const &old_angle
	) const;

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	std::string torsion_name(core::Size const iset) {
		return torsion_name_[iset];
	}

	/// @brief Sets a residue index that will serve as the root of the FoldTree for the SetTorsion operation.
	/// @details The FoldTree is reset afterwards (i.e. the mover does not permanently change the FoldTree).
	void set_fold_tree_root(core::Size const root) { fold_tree_root_ = root; return; }


	/// @brief Returns the residue index that will serve as the root of the FoldTree for the SetTorsion operation.
	///
	core::Size get_fold_tree_root() const { return fold_tree_root_; }

	/// @brief Add a perturbation type to the list of perturbation types.
	/// @details Currently only allows "uniform" or "gaussian".  Checks for proper input.  After this operation, the
	/// perturbation_type_ vector is one entry longer.
	void add_perturbation_type( std::string const &type_in ) {
		if ( type_in=="uniform" ) {
			perturbation_type_.push_back(perturbtorsion_uniform);
			return;
		} else if ( type_in=="gaussian" ) {
			perturbation_type_.push_back(perturbtorsion_gaussian);
			return;
		} else {
			utility_exit_with_message("User input error caught by protocols::simple_moves::SetTorsion::add_perturbation_type():  Type not recognized.  The perturbation type must be \"uniform\" or \"gaussian\".\n");
		}
		return;
	}

	/// @brief Get a perturbation type.
	///
	inline TorsionPerturbType perturbation_type( core::Size const index ) const {
		return perturbation_type_[index];
	}

	/// @brief Add a perturbation magnitude to the list of perturbation magnitudes.
	/// @details Checks for non-negative magnitude.  After this operation, the
	/// perturbation_magnitude_ vector is one entry longer.
	void add_perturbation_magnitude( core::Real const &mag_in ) {
		runtime_assert_string_msg( mag_in>=0, "User input error caught by protocols::simple_moves::SetTorsion::add_perturbation_magnitude():  The perturbation magnitude must be zero or greater; negative values are not allowed." );
		perturbation_magnitude_.push_back(mag_in);
		return;
	}

	/// @brief Get a perturbation magnitude.
	///
	inline core::Real perturbation_magnitude( core::Size const index ) const {
		return perturbation_magnitude_[index];
	}


private:
	bool random_set_;
	utility::vector1< std::string > angle_;
	utility::vector1< std::string > residues_;
	utility::vector1< std::string > torsion_name_; // phi/psi etc.

	/// @brief Used for custom rama sampling.  One entry per <Torsion> block; defaults to "".
	///
	utility::vector1< std::string > custom_rama_map_;

	utility::vector1< Size > extending_;
	utility::vector1< utility::vector1< core::id::NamedAtomID > > torsion_atoms_;

	/// @brief The type of perturbation, if the "perturb" option is used.
	/// @details Current options are "uniform" and "gaussian".  Defaults to "gaussian".
	utility::vector1< TorsionPerturbType  > perturbation_type_;

	/// @brief The perturbation magnitude, if the "perturb" option is used.
	/// @details Defaults to 1.0.
	utility::vector1< core::Real  > perturbation_magnitude_;

	/// @brief The root for the FoldTree during the operation.  If 0, the default FoldTree is used.  The FoldTree is reset to the input FoldTree after the operation.
	///
	core::Size fold_tree_root_;

};

} // moves
} // protocols

#endif //INCLUDED_protocols_simple_moves_SetTorsion_HH_
