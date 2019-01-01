// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/ParaIon.hh
/// @brief   class to store properties of paramagnetic ions
///          such as charge and spin quantum numbers
/// @details last Modified: 05/11/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_io_nmr_ParaIon_HH
#define INCLUDED_core_io_nmr_ParaIon_HH

// Unit headers
#include <core/io/nmr/ParaIon.fwd.hh>

// Project headers
#include <core/types.hh>

// C++ headers
#include <string>
#include <iostream>

namespace core {
namespace io {
namespace nmr {

class ParaIon {

public: // Methods

	/// @brief default constructor
	ParaIon();

	/// @brief constructor with arguments
	ParaIon(
		std::string const& label,
		Real const charge,
		Real const s,
		Real const l,
		Real const j,
		Real const gj,
		Real const te
	);

	/// @brief copy constructor
	ParaIon( ParaIon const & src );

	/// @brief assignment operator
	ParaIon &
	operator=( ParaIon const & src );

	/// @brief destructor
	~ParaIon();

	/// Setters and Getters
	std::string get_ion_label() const { return ion_label_; }
	Real get_charge() const {return charge_; }
	Real get_S() const { return S_; }
	Real get_L() const { return L_; }
	Real get_J() const { return J_; }
	Real get_gJ() const { return gJ_; }
	Real get_tau_e() const { return tau_e_; }

	void set_ion_label(std::string const & label) { ion_label_ = label; }
	void set_charge(Real charge) {charge_ = charge; }
	void set_S(Real s) { S_ = s; }
	void set_L(Real l) { L_ = l; }
	void set_J(Real j) { J_ = j; }
	void set_gJ(Real gj) { gJ_ = gj; }
	void set_tau_e(Real te) { tau_e_ = te; }

	/// @brief calculate Lande g-factor
	void calc_gJ();

	/// @brief output ParaIon data
	void show(std::ostream & out) const;

private: // Data

	/// @brief label of metal ion
	std::string ion_label_;

	/// @brief charge of metal ion
	Real charge_;

	/// @brief spin quantum number
	Real S_;

	/// @brief orbital quantum number
	Real L_;

	/// @brief total quantum number
	Real J_;

	/// @brief Lande g-factor
	Real gJ_;

	/// @brief electron relaxation time (10-12 s)
	Real tau_e_;
};

} // namespace nmr
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_nmr_ParaIon_HH
