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

#ifndef INCLUDED_protocols_noesy_assign_FloatingResonance_hh
#define INCLUDED_protocols_noesy_assign_FloatingResonance_hh


// Unit Headers
#include <protocols/noesy_assign/FloatingResonance.fwd.hh>
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
FloatingResonance combines resonanceID (label), chemical shift (freq), tolerance (error), and the assigned atom (atom, name, resid)
(provided accessor methods of "FloatingResonance": label, atom, resid, name, freq, error, tolerance, calibration_atom_type )
*/

class FloatingResonance : public Resonance {

	//typedefs
public:
	typedef std::set< core::Size > FloatList;
private:
	typedef Resonance Parent;


	//methods
public:
  FloatingResonance();
  FloatingResonance( Resonance const& res, FloatList const&, ResonanceList* );
  ~FloatingResonance();

	virtual core::Real pmatch( core::Real peakfreq, core::Real error, FoldResonance const& folder ) const;
	virtual void write_to_stream( std::ostream& os ) const;
	virtual void write_to_stream( std::ostream&, core::chemical::AA aa ) const;

private:
	void _write_partner_ids( std::ostream& os ) const;
	FloatList partner_ids_;
	ResonanceList const* res_list_;
};

}
}
#endif
