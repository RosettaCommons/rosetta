// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/hbonds/hbtrie/HBondsTrieVsTrieCachedDataContainer.hh
/// @brief A class for passing data to the trie-vs-trie calculation for hydrogen bonds, without having to cache it in mutable data in the HBondEnergy method or whatnot.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_core_scoring_hbonds_hbtrie_HBondsTrieVsTrieCachedDataContainer_hh
#define INCLUDED_core_scoring_hbonds_hbtrie_HBondsTrieVsTrieCachedDataContainer_hh

#include <core/scoring/hbonds/hbtrie/HBondsTrieVsTrieCachedDataContainer.fwd.hh>
#include <core/scoring/trie/TrieVsTrieCachedDataContainerBase.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Core headers
#include <core/scoring/EnergyMap.hh>
#include <core/types.hh>

namespace core {
namespace scoring {
namespace hbonds {
namespace hbtrie {

/// @brief A class for passing data to the trie-vs-trie calculation for hydrogen bonds, without
/// having to cache it in mutable data in the HBondEnergy method or whatnot.
/// @details The EnergyMap passed to this object is stored by raw pointer.  THIS OBJECT DOES NOT OWN THE ENERGY MAP.
/// Calling code must guarantee that the EnergyMap continues to exist until this object is destroyed.  (This is
/// unfortunately needed for speed.)
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class HBondsTrieVsTrieCachedDataContainer : public core::scoring::trie::TrieVsTrieCachedDataContainerBase {

public:

	/// @brief Default constructor.
	HBondsTrieVsTrieCachedDataContainer() = delete;

	/// @brief Weights constructor.
	HBondsTrieVsTrieCachedDataContainer(
		core::scoring::EnergyMap const & setting
	) :
		core::scoring::trie::TrieVsTrieCachedDataContainerBase(),
		weights_( &setting )
	{}

	/// @brief Copy constructor.
	HBondsTrieVsTrieCachedDataContainer(HBondsTrieVsTrieCachedDataContainer const & src);

	/// @brief Destructor.
	~HBondsTrieVsTrieCachedDataContainer() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	core::scoring::trie::TrieVsTrieCachedDataContainerBaseOP clone() const override;

public: //Setters

	/// @brief Set the rotamer sequence separation.
	inline void set_rotamer_seq_sep( signed long const setting ) { rotamer_seq_sep_ = setting; }

	/// @brief Set residue 1.
	inline void set_res1( core::Size const setting ) { res1_ = setting; }

	/// @brief Set residue 2.
	inline void set_res2( core::Size const setting ) { res2_ = setting; }

	/// @brief Set residue 1_nb.
	inline void set_res1_nb( core::Size const setting ) { res1_nb_ = setting; }

	/// @brief Set residue 2_nb.
	inline void set_res2_nb( core::Size const setting ) { res2_nb_ = setting; }

public: //Getters

	/// @brief Get the weights vector.
	inline core::scoring::EnergyMap const & weights() const { return *weights_; }

	/// @brief Get the rotamer sequence separation.
	inline signed long rotamer_seq_sep() const { return rotamer_seq_sep_; }

	/// @brief Get residue 1.
	inline core::Size res1() const { return res1_; }

	/// @brief Get residue 2.
	inline core::Size res2() const { return res2_; }

	/// @brief Get residue 1_nb.
	inline core::Size res1_nb() const { return res1_nb_; }

	/// @brief Get residue 2_nb.
	inline core::Size res2_nb() const { return res2_nb_; }

private:

	core::scoring::EnergyMap const * weights_ = nullptr;
	signed long rotamer_seq_sep_ = 0;
	core::Size res1_ = 0;
	core::Size res2_ = 0;
	core::Size res1_nb_ = 0;
	core::Size res2_nb_ = 0;

};

} //hbtrie
} //hbonds
} //scoring
} //core

#endif //INCLUDED_core_scoring_hbonds_hbtrie_HBondsTrieVsTrieCachedDataContainer_hh
