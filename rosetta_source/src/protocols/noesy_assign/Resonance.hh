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

#ifndef INCLUDED_protocols_noesy_assign_Resonance_hh
#define INCLUDED_protocols_noesy_assign_Resonance_hh


// Unit Headers
//#include <devel/NoesyAssign/ResonanceList.fwd.hh>
#include <core/types.hh>

// Package Headers
#include <core/id/NamedAtomID.hh>
#include <protocols/noesy_assign/PeakCalibrator.fwd.hh>

// Project Headers
// AUTO-REMOVED #include <core/chemical/AA.hh>

// Utility headers
// AUTO-REMOVED #include <utility/exit.hh>
// #include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
//#include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>
//#include <basic/Tracer.hh>
// #include <basic/options/option.hh>
// #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
// AUTO-REMOVED #include <string>
// AUTO-REMOVED #include <map>

namespace protocols {
namespace noesy_assign {
/*! @detail
Resonance combines resonanceID (label), chemical shift (freq), tolerance (error), and the assigned atom (atom, name, resid)
(provided accessor methods of "Resonance": label, atom, resid, name, freq, error, tolerance, calibration_atom_type )
*/

class Resonance {
public:
  Resonance();
  Resonance( core::Size label, core::Real freq, core::Real error, core::id::NamedAtomID id );
  ~Resonance();

	///@brief output
  void write_to_stream( std::ostream& ) const;

	///@brief ResonanceID
  core::Size label() const { return label_; }

	///@brief Atom
  core::id::NamedAtomID const& atom() const { return atom_; }
  core::Size resid() const { return atom_.rsd(); }
  std::string const& name() const { return atom_.atom(); }

	///@brief resonance frequency (chemical shift)
  core::Real freq() const { return freq_; }
	core::Real error() const { return error_; }
	core::Real tolerance() const { return error_; }

	///@brief classification for calibration... e.g., Backbone, sidechain, etc..
	CALIBRATION_ATOM_TYPE calibration_atom_type() const { return calibration_atom_type_; }

private:
  core::Size label_;
  core::Real freq_;
  core::Real error_;
  core::id::NamedAtomID atom_;
	CALIBRATION_ATOM_TYPE calibration_atom_type_;
};

}
}
#endif
