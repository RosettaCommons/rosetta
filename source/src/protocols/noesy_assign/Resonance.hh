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
#include <protocols/noesy_assign/Resonance.fwd.hh>

// Package Headers
#include <protocols/noesy_assign/FoldResonance.fwd.hh>
#include <protocols/noesy_assign/PeakCalibrator.fwd.hh>

// Project Headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/chemical/AA.hh>
#include <core/types.hh>
#include <core/id/NamedAtomID.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/access_ptr.hh>
//// C++ headers
#include <deque>

namespace protocols {
namespace noesy_assign {
/*! @detail
Resonance combines resonanceID (label), chemical shift (freq), tolerance (error), and the assigned atom (atom, name, resid)
(provided accessor methods of "Resonance": label, atom, resid, name, freq, error, tolerance, calibration_atom_type )
*/

class Resonance : public utility::pointer::ReferenceCount {

public:
	typedef utility::vector1< core::Size > ResonanceIDs;
	typedef std::pair< core::Size, core::Size > ResonancePair;
	typedef utility::vector1< ResonancePair > ResonancePairs;
  typedef utility::vector1< ResonanceAP > ResonanceAPs;

  Resonance();
  Resonance( core::Size label, core::Real freq, core::Real error, core::id::NamedAtomID const& id, core::chemical::AA, core::Real intensity = 1.0 );
  ~Resonance();

	virtual ResonanceOP clone() = 0;

	/// @brief output
  virtual void write_to_stream( std::ostream& ) const;

	virtual void write_to_stream( std::ostream&, core::chemical::AA aa ) const;

	virtual core::Size ambiguity() const {
		return 1;
	}
	/// @brief ResonanceID
  core::Size label() const { return label_; }
	virtual core::Size float_label( core::Size /*ifloat*/ ) const { return label_; }
	/// @brief Atom
  core::id::NamedAtomID const& atom() const { return atom_; }
  core::Size resid() const { return atom_.rsd(); }
  std::string const& name() const { return atom_.atom(); }
	bool is_proton() const { return is_proton_; }
	/// @brief resonance frequency (chemical shift)
  core::Real freq() const { return freq_; }
	core::Real error() const { return error_; }
	core::Real tolerance() const { return error_; }

	/// @brief Resonance matches the given cross-peaks frequency
	bool match( core::Real freq, core::Real error, FoldResonance const& folder ) const {
		return pmatch( freq, error, folder ) <= 1.0;
	}

	/// @brief match the proton and corresponding label atom at same time
	virtual bool match2D(
	  core::Real proton_freq,
		core::Real proton_error,
		FoldResonance const& proton_folder,
		core::Real label_freq,
		core::Real label_error,
		FoldResonance const& label_folder,
		ResonancePairs& matches
	) const = 0;

	void add_connected_resonance( ResonanceAP ptr );
	void clear_connected_resonances();
	bool has_connected_resonances() const {
		return connected_resonance_ids_.size() > 0;
	}
	Resonance const& first_connected_resonance() const;
	ResonanceIDs const& connected_resonance_ids() const {
		return connected_resonance_ids_;
	}
	ResonanceAPs const& connected_resonances() const;

	virtual core::Real pmatch(  core::Real freq, core::Real error, FoldResonance const& folder ) const;

	void combine( std::deque< ResonanceOP >& last_resonances, bool drain );

	core::chemical::AA aa() const { return aa_; }
	/// @brief in ILV-labelled proteins, the both LV methyls are labelled randomly with 50% probability,
	/// whereas I delta methyls are labelled 100%
	core::Real intensity() const { return intensity_; }
	void set_intensity( core::Real setting ) {
		intensity_ = setting;
	}

	/// @brief classification for calibration... e.g., Backbone, sidechain, etc..
	CALIBRATION_ATOM_TYPE calibration_atom_type() const { return calibration_atom_type_; }

	core::Real _pmatch(  core::Real freq, core::Real error, FoldResonance const& folder ) const;

private:
	void _write_to_stream( std::ostream& ) const;
  core::Size label_;
  core::Real freq_;
  core::Real error_;
	bool is_proton_;
  core::id::NamedAtomID atom_;
	core::chemical::AA aa_;
	core::Real intensity_;
	CALIBRATION_ATOM_TYPE calibration_atom_type_;
	ResonanceIDs connected_resonance_ids_;
	ResonanceAPs connected_resonance_ptrs_;
};

}
}
#endif
