// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CrossPeakList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_LabelResonance_hh
#define INCLUDED_protocols_noesy_assign_LabelResonance_hh


// Unit Headers
#include <protocols/noesy_assign/LabelResonance.fwd.hh>
#include <protocols/noesy_assign/ProtonResonance.fwd.hh>
#include <protocols/noesy_assign/Resonance.hh>

// Package Headers
#include <protocols/noesy_assign/ResonanceList.fwd.hh>
#include <protocols/noesy_assign/PeakCalibrator.fwd.hh>


// Project Headers
#include <core/chemical/AA.hh>
#include <core/types.hh>
#include <core/id/NamedAtomID.hh>

// Utility headers
#include <utility/vector1.hh>

//// C++ headers
#include <set>

namespace protocols {
namespace noesy_assign {
/*! @detail
LabelResonance combines resonanceID (label), chemical shift (freq), tolerance (error), and the assigned atom (atom, name, resid)
(provided accessor methods of "LabelResonance": label, atom, resid, name, freq, error, tolerance, calibration_atom_type )
*/

class LabelResonance : public Resonance {
  //typedefs
public:
private:
  typedef Resonance Parent;
  //methods
public:
  LabelResonance();

  LabelResonance(
      core::Size label,
      core::Real freq,
      core::Real error,
      core::id::NamedAtomID const& id,
      core::chemical::AA,
			core::Real intensity
  );

  ~LabelResonance();

  virtual ResonanceOP clone() {
    return new LabelResonance( *this );
  }

  ///@brief match the proton and corresponding label atom at same time
  virtual bool match2D(
    core::Real proton_freq,
    core::Real proton_error,
    FoldResonance const& proton_folder,
    core::Real label_freq,
    core::Real label_error,
    FoldResonance const& label_folder,
		ResonancePairs& matches
  ) const;

};

}
}
#endif
