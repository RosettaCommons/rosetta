// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
	~FloatingResonance() override;

	ResonanceOP clone() override {
		return ResonanceOP( new FloatingResonance( *this ) );
	}

	core::Real pmatch( core::Real peakfreq, core::Real error, FoldResonance const& folder ) const override;
	void write_to_stream( std::ostream& os ) const override;
	void write_to_stream( std::ostream&, core::chemical::AA aa ) const override;
	core::Size ambiguity() const override {
		return partner_ids_.size();
	}
	core::Size float_label( core::Size ifloat ) const override;

	/// @brief match the proton and corresponding label atom at same time
	bool match2D(
		core::Real proton_freq,
		core::Real proton_error,
		FoldResonance const& proton_folder,
		core::Real label_freq,
		core::Real label_error,
		FoldResonance const& label_folder,
		ResonancePairs& matches
	) const override;

private:
	void _write_partner_ids( std::ostream& os ) const;

	//only one resonance of the group should act as representative, we just take the one with smallest ID
	bool is_representative_resonance() const {
		return is_representative_resonance_;
	}

	FloatList partner_ids_;
	ResonanceList const* res_list_;
	bool is_representative_resonance_;
};

}
}
#endif
