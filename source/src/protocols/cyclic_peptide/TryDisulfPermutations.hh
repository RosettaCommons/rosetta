// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/cyclic_peptide/TryDisulfPermutations.hh
/// @brief  Headers for TryDisulfPermutations.cc.  Tries all permutations of disulfide bonds between
/// disulfide-forming residues in a pose.
/// @details  Performs a repack/minimize of the disulfide-forming residues with a simplified score function.
/// Returns the permutation with the lowest disulfide energy.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_TryDisulfPermutations_hh
#define INCLUDED_protocols_cyclic_peptide_TryDisulfPermutations_hh

// Unit Headers
#include <protocols/cyclic_peptide/TryDisulfPermutations.fwd.hh>
#include <protocols/moves/Mover.hh>

// Scripter Headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/constants.hh>


///////////////////////////////////////////////////////////////////////

namespace protocols {
namespace cyclic_peptide {

class TryDisulfPermutations : public protocols::moves::Mover
{
public: //Typedefs


public: //Constructors, destructors, clone operator, etc.

	TryDisulfPermutations();

	TryDisulfPermutations( TryDisulfPermutations const &src );

	~TryDisulfPermutations() override;

	protocols::moves::MoverOP clone() const override;

	protocols::moves::MoverOP fresh_instance() const override;


	/// @brief Actually apply the mover to the pose.
	void apply(core::pose::Pose & pose) override;

	// XRW TEMP  std::string get_name() const override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const &
	) override;

	/// @brief Will this mover consider alternative disulfides involving residues currently in a disulfide bond?
	///
	inline bool consider_already_bonded() const { return consider_already_bonded_; }

	/// @brief Set whether this mover will consider alternative disulfides involving residues currently in a disulfide bond.
	///
	inline void set_consider_already_bonded( bool const setting ) { consider_already_bonded_=setting; return; }

	/// @brief Get the minmization flavour.
	///
	inline std::string const &mintype() const { return mintype_; }

	/// @brief Set the minimization flavour.
	///
	inline void set_mintype( std::string const &type) { mintype_ = type; return; }

	/// @brief Get the minimizer tolerance.
	///
	inline core::Real const &mintolerance() const { return mintolerance_; }

	/// @brief Set the minimizer tolerance.
	///
	inline void set_mintolerance( core::Real const &val ) { mintolerance_ = val; return; };

	/// @brief Set the residue selector.
	/// @details CLONES the input.
	void set_selector( core::select::residue_selector::ResidueSelectorCOP selector_in );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE DATA                                                      //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Should we consider residues that are already disulfide-bonded (i.e. let them form other
	/// disulfides), or should we just consider unbonded disulfide-formers?
	/// @details Default true (considers everything).
	bool consider_already_bonded_;

	/// @brief The minimizer behaviour.
	///
	std::string mintype_;

	/// @brief The minimizer tolerance.
	///
	core::Real mintolerance_;

	/// @brief An optional ResidueSelector, for applying this mover to a subset of a pose.
	///
	core::select::residue_selector::ResidueSelectorOP selector_;

	////////////////////////////////////////////////////////////////////////////////
	//          PRIVATE FUNCTIONS                                                 //
	////////////////////////////////////////////////////////////////////////////////

	/// @brief Given a list of disulfide-forming residues, generate all possible permtutations.
	/// @param[in] residues The list of disulfide-forming positions.
	/// @param[out] permutations The list of lists of disulfide pairs representing all possible permutations.
	/// @details Note that this function does NOT clear the permutations vector; it only appends to it.
	void generate_disulf_permutations(
		utility::vector1 < core::Size > const &residues,
		utility::vector1 < utility::vector1 < std::pair< core::Size, core::Size > > > &permutations
	) const;


	/// @brief Given a permutation index, generate a permutation.
	/// @details Assumes that cur_permutation is an empty vector.  Recursively calls itself.
	void recursively_generate_permutation (
		utility::vector1 < std::pair < core::Size, core::Size > > &cur_permutation,
		utility::vector1 < core::Size > const &residues,
		utility::vector1 < bool > &placed,
		utility::vector1 < core::Size > const &index,
		core::Size const &level
	) const;

	/// @brief Increment the perturbation index.
	/// @details Recursively calls itself.
	void increment_index(
		utility::vector1 < core::Size > &index,
		core::Size const index_index,
		core::Size const max_levels
	) const;

	/// @brief Given a set of disulfides to form, make these disulfide bonds.
	///
	void generate_disulfides(
		core::pose::PoseOP pose,
		utility::vector1 < std::pair <core::Size, core::Size> > disulf_pairs
	) const;

	/// @brief Repack and minimize a set of disulfides, using only the fa_dslf and fa_dun score terms.
	/// @details Returns the fa_dslf energy after repacking and minimization.
	core::Real repack_minimize_disulfides(
		core::pose::PoseOP pose,
		utility::vector1 < core::Size > const &disulf_res
	) const;

	/// @brief Is a value in a list?
	///
	inline bool is_in_list (core::Size const val, utility::vector1 <core::Size> const &list) const {
		for ( core::Size i=1, imax=list.size(); i<=imax; ++i ) {
			if ( list[i] == val ) return true;
		}
		return false;
	}

};

} //namespace cyclic_peptide
} //namespace protocols

#endif //INCLUDED_protocols_cyclic_peptide_TryDisulfPermutations_hh
