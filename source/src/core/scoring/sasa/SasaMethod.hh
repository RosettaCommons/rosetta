// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/sasa/SasaMethod.hh
/// @brief Abstract class defining a 'SasaMethod'
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_scoring_sasa_SASAMETHOD_HH
#define INCLUDED_core_scoring_sasa_SASAMETHOD_HH

#include <core/pose/Pose.hh>
#include <core/scoring/sasa/SasaMethod.fwd.hh>
#include <core/id/AtomID_Map.hh>

namespace core {
namespace scoring {
namespace sasa {

	///@brief Type of Radii to use.
	///@details
	///       LJ:  Refers to Leonard Jones radii - Rosetta uses radii at the minimum of the potential (sigma2).
	///       Legacy:  Refers to radii optimized for a no longer in use term, but some protocols have been optimized to use it.
	///       naccess:  Refers to radii used in the program naccess.  Originally derived from Chothia.  Do not use for all-atom SASA as hydrogens are implicitly included.
	///                           'The Nature of the Accessible and Buried Surfaces in Proteins' J. Mol. Biol. (1976) 105, 1-14
	///       reduce:   Radii used by the program reduce.  Hydrogens are explicitly included in the radii.
	///
	enum SasaRadii {
		LJ = 1,
		legacy,
		naccess,
		reduce,

		chothia=naccess,
		SasaRadii_total = reduce
	};


///@brief Abstract base class for SasaMethods.  Feel free to edit as needed.
class SasaMethod : public utility::pointer::ReferenceCount {


public:

	SasaMethod(Real probe_radius, SasaRadii radii_set);
	virtual ~SasaMethod();

	///@brief Calculate Sasa.  Atoms not calculated have -1 sasa in AtomID_Map.  This is carried over for compatability purposes.
	virtual Real
	calculate(
			const pose::Pose & pose,
			const id::AtomID_Map<bool> & atom_subset,
			id::AtomID_Map< Real > & atom_sasa,
			utility::vector1< Real > & rsd_sasa) = 0;

	virtual std::string
	get_name() const = 0;


public:

	///@brief Include the probe radius in calc.  Typical for SASA.
	void
	set_include_probe_radius_in_calc(bool include_probe_radius);

	///@brief Set the probe radius.  Typical value is that of water at 1.4 A
	void
	set_probe_radius(Real probe_radius);

	///@brief Set the radii type.
	void
	set_radii_set(SasaRadii radii_set);



///////////Legacy Options ///////////
public:
	//void
	//set_expand_polar_radii(bool expand_polars, core::Size expansion_radius = 1.0);

	///@brief Legacy option to increase polar hydrogen radii to 1.08A.  Supported for now.
	void
	set_use_big_polar_hydrogen(bool big_polar_h);


protected:
	Real probe_radius_;
	SasaRadii radii_set_;

	bool include_probe_radius_;
	bool use_big_polar_H_;
	//vector1<std::string> radii_names_;
};



}
}
}

#endif	//#ifndef INCLUDED_protocols/antibody_design_SASAMETHOD_HH
