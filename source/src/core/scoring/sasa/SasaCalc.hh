// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core/scoring/sasa/SasaCalc.hh
/// @brief Class based interface to original sasa functions.  Should result in less incorrect uses of radii.
///  
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_scoring_sasa_SASACALC_HH
#define INCLUDED_core_scoring_sasa_SASACALC_HH

#include <core/scoring/sasa/SasaMethodFactory.hh>
#include <core/scoring/sasa/SasaCalc.fwd.hh>

namespace core {
namespace scoring {
namespace sasa {
	using utility::vector1;
	using namespace core;
	
///@brief Main 
class SasaCalc : public utility::pointer::ReferenceCount {
	
public:
	
	SasaCalc();
	
	SasaCalc(SasaMethodEnum method);
	
	
	virtual ~SasaCalc();
	
	///@brief Calculate Sasa.  Atoms not calculated have -1 sasa in AtomID_Map.  This is carried over for compatability purposes.  	
	Real
	calculate(const pose::Pose & pose);
	
	
	/////Legacy-style interface //////////
	//Real
	//calculate(const pose::Pose & pose,  id::AtomID_Map<Real>& atom_sasa);
	
	//Real
	//calculate(const pose::Pose & pose, vector1< Real>&  rsd_sasa, vector1< Real > & rsd_hsasa);
	
	//Real
	//calculate(const pose::Pose & pose, id::AtomID_Map<Real>& atom_sasa, vector1< Real> & rsd_sasa, vector1<Real> & rsd_hsasa);


	
//////Data Access //////////////////
public:

	///////Per Atom
	id::AtomID_Map< Real >
	get_atom_sasa() const {
		return atom_sasa_;
	};
	
	
	
	////////Per Residue
	vector1<Real>
	get_residue_sasa() const {
		return rsd_sasa_;
	};
	
	vector1< Real >
	get_residue_hsasa() const {
		return rsd_hsasa_;
	};
	
	vector1<Real >
	get_rel_hphobic_sasa_by_charge() const {
		return rel_hydrophobic_sasa_by_charge_;
	};
	
	///@brief Convenience function to fill all data. 
	void
	fill_all_data(
		Real & total_hsasa,
		id::AtomID_Map< Real > & atom_sasa,
		vector1< Real > & rsd_sasa,
		vector1< Real > & rsd_hsasa,
		vector1< Real > & rel_hsasa);
	
	/////////Totals
	Real
	get_total_sasa() const {
		return total_sasa_;
	}
	
	Real
	get_total_hsasa() const {
		return total_hsasa_;
	}
	
	Real
	get_total_rel_hsasa() const {
		return total_rel_hsasa_;
	}
	
	
	
	
/////Options /////////////
public:
	
	void
	set_defaults();
	
	void
	set_calculation_method(SasaMethodEnum method);
	
	///@brief Include hydrogens explicitly
	void
	set_include_hydrogens_explicitly(bool include_hydrogens);
	
	void
	set_probe_radius(Real probe_radius);
	
	
	void
	set_include_probe_radius_in_atom_radii(bool include_probe_radius); 
	
	
	///Hydrophobic Calc - Definition of 'polar' is by charge.  Better would be to use orbitals/chem if possible.
	
	
	
	///@brief Typically, only carbon or sulfers are included in the calculation.  If you are using ligands - this may not be good enough.
	void
	set_include_carbon_sulfer_only_in_hydrophobic_calc(bool include_c_s_only);
	
	///@brief Polar carbons and other atoms should not be included in hydrophobic hSASA - though historically they were.  
	/// .4 is a relative number.  This makes sure that carbonyl and carboxyl carbons are marked as polar, while others (protein-based) are non-polar
	void
	set_exclude_polar_atoms_by_charge(bool exclude_polar_all, Real charge_cutoff = .4);
	
	
	void
	set_polar_charge_cutoff(Real cutoff);
	
	
	

	
	
	///@brief Radii set to use when not including hydrogens (naccess/chothia, reduce, legacy)
	/// Do not use legacy unless you know what you are doing and why. Overrides default of naccess.
	void
	set_implicit_hydrogen_included_radii_set(SasaRadii radii_set);
	
	
	
	
///////////Legacy Options ///////////
public:
	
	
	//set_expand_polar_radii(bool expand_polars, core::Size expansion_radius = 1.0); (Not used anywhere)
	
	///@brief Not for general use.  Used to calculate unsatisfied buried polars with legacy radii (which implicitly had included hydrogens.)
	void
	set_use_big_polar_hydrogen(bool big_polar_h);
	
	///@brief Sets options to reproduce legacy default radii behavior.  Note this is wrong, and should only be used where needed.
	///@details
	///  Calculates all atom Sasa using Rosetta legacy radii - which were optimized for a no-longer used score-term and scorefunction.
	///  Note that the radii implicitly include hydrogens while the SASA will be calculated for all atoms including hydrogens.
	/// 
	void
	set_use_legacy_default_radii_with_all_atom_calc( bool legacy_defaults);
	
private:
	
	void
	init(const pose::Pose & pose);
	
	///@brief Creates and sets up the sasa method before calculate is called on the instance. 
	void
	setup_sasa_method(SasaRadii radii_set);
	
	
	void
	calc_per_res_hphobic_sasa(const pose::Pose & pose);
	
	
	
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
	vector1< Real > rsd_sasa_;
	vector1< Real > rsd_hsasa_;
	vector1< Real > rel_hydrophobic_sasa_by_charge_;
	
	//Totals
	Real total_sasa_;
	Real total_hsasa_;
	Real total_rel_hsasa_;
	
	bool legacy_defaults_;
	
};



} //sasa
} //scoring
} //core


#endif	//#ifndef INCLUDED_core_scoring_sasa_SASACALC_HH

