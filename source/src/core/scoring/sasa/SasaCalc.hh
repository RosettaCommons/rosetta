// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/sasa/SasaCalc.hh
/// @brief Class based interface to original sasa functions.  Should result in less incorrect uses of radii.
///
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_scoring_sasa_SASACALC_HH
#define INCLUDED_core_scoring_sasa_SASACALC_HH

#include <core/scoring/sasa/SasaMethodFactory.hh>
#include <core/scoring/sasa/SasaCalc.fwd.hh>

#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace sasa {

/// @brief Main interface to sasa calculations outside of pose metrics.
class SasaCalc : public utility::pointer::ReferenceCount {

public:

	SasaCalc();

	SasaCalc(SasaMethodEnum method);


	virtual ~SasaCalc();

	/// @brief Calculate Sasa.  Atoms not calculated have -1 sasa in AtomID_Map.  This is carried over for compatability purposes.
	Real
	calculate(const pose::Pose & pose);


	///// Legacy-style interfaces //////////
	Real
	calculate(const pose::Pose & pose,  id::AtomID_Map<Real>& atom_sasa);

	Real
	calculate(const pose::Pose & pose, utility::vector1< Real>&  rsd_sasa, utility::vector1< Real > & rsd_hsasa);

	Real
	calculate(const pose::Pose & pose, id::AtomID_Map<Real>& atom_sasa, utility::vector1< Real> & rsd_sasa, utility::vector1<Real> & rsd_hsasa);

	Real
	calculate(const pose::Pose & pose, id::AtomID_Map<Real> & atom_sasa, utility::vector1< Real > & rsd_sasa, utility::vector1<Real> & rsd_hsasa, utility::vector1< Real > & rsd_rel_hsasa);


	////////////////////////////////////////////////////////////////////////////////
	/// Data Access
	///
	///
public:

	///////Per Atom
	id::AtomID_Map< Real >
	get_atom_sasa() const {
		return atom_sasa_;
	};

	/// @brief Convenience function to fill most commonly used data.
	void
	fill_data(
		Real & total_hsasa,
		Real & total_rel_hsasa,
		id::AtomID_Map< Real > & atom_sasa,
		utility::vector1< Real > & rsd_sasa,
		utility::vector1< Real > & rsd_hsasa,
		utility::vector1< Real > & rel_hsasa);


	////////Per Residue
	utility::vector1<Real>
	get_residue_sasa() const {
		return rsd_sasa_;
	};

	utility::vector1< Real >
	get_residue_sasa_bb() const;

	utility::vector1< Real >
	get_residue_sasa_sc() const {
		return rsd_sasa_sc_;
	}


	utility::vector1< Real >
	get_residue_hsasa() const {
		return rsd_hsasa_;
	};

	utility::vector1< Real >
	get_residue_hsasa_bb() const;

	utility::vector1< Real >
	get_residue_hsasa_sc() const {
		return rsd_hsasa_sc_;
	}


	utility::vector1<Real >
	get_rel_hphobic_sasa_by_charge() const {
		return rel_hydrophobic_sasa_by_charge_;
	};


	/////////Totals
	Real
	get_total_sasa() const {
		return total_sasa_;
	}

	Real
	get_total_sasa_sc() const {
		return total_sasa_sc_;
	}

	Real
	get_total_sasa_bb() const {
		return total_sasa_ - total_sasa_sc_;
	}


	Real
	get_total_hsasa() const {
		return total_hsasa_;
	}

	Real
	get_total_hsasa_sc() const {
		return total_hsasa_sc_;
	}

	Real
	get_total_hsasa_bb() const {
		return total_hsasa_ - total_hsasa_sc_;
	}

	Real
	get_total_rel_hsasa() const {
		return total_rel_hsasa_;
	}


	////////////////////////////////////////////////////////////////////////////////
	/// Options
	///
	///
public:


	/////////// Common Options ///////////
	void
	set_defaults();

	void
	set_calculation_method(SasaMethodEnum method);

	/// @brief Include hydrogens explicitly
	void
	set_include_hydrogens_explicitly(bool include_hydrogens);

	/// @brief Probe radius of 1.4 (water) is typically used
	void
	set_probe_radius(Real probe_radius);

	/// @brief This is typically done.  Disabling it is more akin to obtaining the Surface Area than the SASA
	void
	set_include_probe_radius_in_atom_radii(bool include_probe_radius);


	/////////// Hydrophobic Calculation ///////////


	/// @brief Typically, only carbon or sulfers are included in the calculation.  If you are using ligands - this may not be good enough.
	void
	set_include_carbon_sulfer_only_in_hydrophobic_calc(bool include_c_s_only);

	/// @brief Polar carbons and other atoms should not be included in hydrophobic hSASA - though historically they were.
	/// .4 is a relative number.  This makes sure that carbonyl and carboxyl carbons are marked as polar, while others (protein-based) are non-polar
	void
	set_exclude_polar_atoms_by_charge(bool exclude_polar_all);


	void
	set_polar_charge_cutoff(Real cutoff);


	/////////// Radii Sets ///////////

	/// @brief Radii set to use when not including hydrogens (naccess/chothia, legacy)
	/// Do not use legacy unless you know what you are doing and why.
	void
	set_implicit_hydrogen_included_radii_set(SasaRadii radii_set);

	/// @brief Radii set to use when including hydrogens (LJ/reduce)
	void
	set_explicit_hydrogen_included_radii_set(SasaRadii radii_set);


	////////////////////////////////////////////////////////////////////////////////
	/// Legacy Options
	///
	///
public:


	//set_expand_polar_radii(bool expand_polars, core::Size expansion_radius = 1.0); (Not used anywhere)

	/// @brief Not for general use.  Used to calculate unsatisfied buried polars with legacy radii (which implicitly had included hydrogens.)
	void
	set_use_big_polar_hydrogen(bool big_polar_h);


private:

	void
	init(const pose::Pose & pose);

	/// @brief Creates and sets up the sasa method before calculate is called on the instance.
	void
	setup_sasa_method(SasaRadii radii_set);


	void
	calc_per_res_sasas(const pose::Pose & pose);


private:

	SasaMethodEnum method_type_;
	SasaMethodOP method_;

	SasaRadii implicit_radii_set_; //Radii to use when not including hydrogens implicitly
	SasaRadii explicit_radii_set_; //Radii to use when including hydrogens explicity

	bool include_hydrogens_;
	Real probe_radius_;

	bool include_c_s_only_in_hsasa_;
	bool exclude_polar_all_in_hsasa_;
	Real polar_charge_cutoff_;

	bool include_probe_radius_;
	//bool expand_polars_;
	//Real polar_expansion_radius_;

	bool big_polar_h_; //Legacy - no idea.

	id::AtomID_Map< bool > atom_subset_;
	id::AtomID_Map< Real > atom_sasa_;

	utility::vector1< Real > rsd_sasa_;
	utility::vector1< Real> rsd_sasa_sc_;

	utility::vector1< Real > rsd_hsasa_;
	utility::vector1< Real > rsd_hsasa_sc_;

	utility::vector1< Real > rel_hydrophobic_sasa_by_charge_;

	//Totals
	Real total_sasa_;
	Real total_sasa_sc_;

	Real total_hsasa_;
	Real total_hsasa_sc_;

	Real total_rel_hsasa_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //sasa
} //scoring
} //core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_sasa_SasaCalc )
#endif // SERIALIZATION


#endif //#ifndef INCLUDED_core_scoring_sasa_SASACALC_HH
