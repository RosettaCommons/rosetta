// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/PcsTensor.hh
///
/// @brief Hold chi-tensor information of PCS
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references JBNMR 2008  41:179-189 schmitz et all will explains the tensor convention used
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsTensor_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsTensor_hh

// Package headers

// Project headers
#include <core/types.hh>

// Utility headers

// Numeric headers

// Objexx headers

// C++ headers
#include <string>

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

class PcsTensor {
private:

	core::Real a_; //alpha
	core::Real b_; //beta
	core::Real g_; //gamma
	core::Real ax_; //axial component
	core::Real rh_; //rhombic component
	//ideally, alpha beta gamma ax and rh will follow the UTR convention
	//see numbat paper JBNMR 2008  41:179-189 schmitz et all

	core::Real chi_xx_;
	core::Real chi_yy_;
	core::Real chi_xy_;
	core::Real chi_xz_;
	core::Real chi_yz_;

	std::string label_; //most likely the filename of the pcs exepriment
	//Used to identify which data has been used in case of multiple lanthanides

private:

	/// @brief NOT READY. This should set alpha beta gamma axial rhombic
	void
	set_abgar();

public:
	PcsTensor(); //Construct

	~PcsTensor(); //destruct

	PcsTensor(PcsTensor const & other); //copy

	PcsTensor & //=
	operator=(PcsTensor const & other);

	PcsTensor(core::Real const chi_xx,
		core::Real const chi_yy,
		core::Real const chi_xy,
		core::Real const chi_xz,
		core::Real const chi_yz,
		std::string const label);

	/// @brief Give me delta chi_xx
	core::Real
	get_delta_X_xx() const;

	/// @brief Give me delta chi_yy
	core::Real
	get_delta_X_yy() const;

	/// @brief Give me delta chi_zz
	core::Real
	get_delta_X_zz() const;

	/// @brief Give me chi_zz
	core::Real
	get_chi_zz() const;

	/// @brief Give me chi_xx
	core::Real
	get_chi_xx() const;

	/// @brief Give me chi_yy
	core::Real
	get_chi_yy() const;

	/// @brief Give me chi_xy
	core::Real
	get_chi_xy() const;

	/// @brief Give me chi_xz
	core::Real
	get_chi_xz() const;

	/// @brief Give me chi_yz
	core::Real
	get_chi_yz() const;

	/// @brief Reset the tensor from the other tensor value
	void
	reset_from_ref(PcsTensor & other);

	/// @brief Give me the tensor label
	std::string const &
	get_label() const;

	/// @brief Reset the tensor from the given values.
	void
	reset_tensor(core::Real const chi_xx,
		core::Real const chi_yy,
		core::Real const chi_xy,
		core::Real const chi_xz,
		core::Real const chi_yz);

	/// @brief Print me
	friend std::ostream &
	operator<<(std::ostream& out, const PcsTensor &me);
};


}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
#endif
