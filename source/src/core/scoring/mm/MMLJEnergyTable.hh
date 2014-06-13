// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMLJEnergyTable.hh
/// @brief  Molecular mechanics lj score class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_scoring_mm_MMLJEnergyTable_hh
#define INCLUDED_core_scoring_mm_MMLJEnergyTable_hh

// Unit headers
#include <core/scoring/mm/MMLJEnergyTable.fwd.hh>
#include <core/scoring/mm/MMLJScore.fwd.hh>
#include <core/scoring/mm/MMLJScore.hh>
#include <core/scoring/mm/MMLJLibrary.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility header
#include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers

namespace core {
namespace scoring {
namespace mm {

/// @brief blah
class MMLJEnergyTable : public utility::pointer::ReferenceCount
{

public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~MMLJEnergyTable();

  /// @briefs typedefs
  typedef utility::vector1< Real > EnergyVector;
  typedef utility::vector1< utility::vector1< EnergyVector* > > MMLJScoreTable;

  /// @brief Default ctor
  MMLJEnergyTable();

  /// @brief blah
  MMLJScore const &
  mm_lj_score() const
  { return mm_lj_score_; }

	/// @brief blah
	Real
	max_dist() const
	{ return max_dist_; }

  /// @brief blah
  void
  score( Size atom1, Size atom2, Size & path_distance, Real & squared_distance, Real & rep, Real & atr ) const;

  /// @brief blah
  void
  deriv_score( Size atom1, Size atom2, Size & path_distance, Real & squared_distance, Real & rep, Real & atr  ) const;

private:

  /// @brief Local MMLJLibrary for looking up lj parameters
  mm::MMLJScore mm_lj_score_;

  MMLJScoreTable mm_atom_pair_rep_energy_table_;
  MMLJScoreTable mm_atom_pair_atr_energy_table_;
  MMLJScoreTable mm_atom_pair_rep_deriv_table_;
  MMLJScoreTable mm_atom_pair_atr_deriv_table_;

  MMLJScoreTable mm_atom_pair_rep_three_bond_energy_table_;
  MMLJScoreTable mm_atom_pair_atr_three_bond_energy_table_;
  MMLJScoreTable mm_atom_pair_rep_three_bond_deriv_table_;
  MMLJScoreTable mm_atom_pair_atr_three_bond_deriv_table_;

  Real max_dist_; // in angstrom squared
  Real bin_dist_; // in angstrom squared
  Size bins_per_angstrom_squared_;

	Real linear_switch_point; // percent of distance at minimum
};

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_mm_lj_energy_table_HH
