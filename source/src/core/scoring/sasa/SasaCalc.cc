// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/core/scoring/sasa/SasaCalc.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/scoring/sasa/SasaCalc.hh>
#include <core/scoring/sasa/SasaMethodFactory.hh>
#include <core/scoring/sasa/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/sasa.OptionKeys.gen.hh>

#include <core/pose/util.tmpl.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR("core.scoring.sasa.SasaCalc");

namespace core {
namespace scoring {
namespace sasa {
	using namespace core;
	using utility::vector1;
	
SasaCalc::SasaCalc():
	utility::pointer::ReferenceCount()
{
	method_type_ = LeGrand;
	set_defaults();
}

SasaCalc::SasaCalc(SasaMethodEnum method):
	utility::pointer::ReferenceCount()
{
	method_type_ = method;
	set_defaults();
}


SasaCalc::~SasaCalc(){}




void
SasaCalc::set_defaults() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	

	set_include_hydrogens_explicitly(option[OptionKeys::sasa::include_hydrogens_explicitly]()); 
	
	set_include_probe_radius_in_atom_radii(option[OptionKeys::sasa::include_probe_radius_in_atom_radii]());
	set_include_carbon_sulfer_only_in_hydrophobic_calc(option[OptionKeys::sasa::include_only_C_S_in_hsasa]());
	set_exclude_polar_atoms_by_charge(option[OptionKeys::sasa::exclude_polar_atoms_by_charge_in_hsasa]());
	set_implicit_hydrogen_included_radii_set(get_sasa_radii_set_from_string(option[OptionKeys::sasa::implicit_hydrogen_radii_set]()));
	explicit_radii_set_ = get_sasa_radii_set_from_string(option[OptionKeys::sasa::explicit_hydrogen_radii_set]());
	
	//set_expand_polar_radii(false);
	set_use_big_polar_hydrogen(false);
	
	set_use_legacy_default_radii_with_all_atom_calc(option[OptionKeys::sasa::use_legacy_behavior]());
	
}

void
SasaCalc::set_use_legacy_default_radii_with_all_atom_calc( bool legacy_defaults){
	
	legacy_defaults_ = legacy_defaults;

}

void
SasaCalc::set_include_hydrogens_explicitly(bool include_hydrogens) {
	include_hydrogens_ = include_hydrogens;
}

void
SasaCalc::set_probe_radius(core::Real probe_radius ) {
	probe_radius_ = probe_radius;
}


void
SasaCalc::set_polar_charge_cutoff(core::Real cutoff) {
	polar_charge_cutoff_= cutoff;
}

void
SasaCalc::set_include_probe_radius_in_atom_radii(bool include_probe_radius) {
	include_probe_radius_ = include_probe_radius;
}

void
SasaCalc::setup_sasa_method(SasaRadii radii_set) {
	method_ = create_sasa_method(method_type_, probe_radius_, radii_set);
	method_->set_include_probe_radius_in_calc(include_probe_radius_);
	method_->set_use_big_polar_hydrogen(big_polar_h_);
	
}

void
SasaCalc::set_implicit_hydrogen_included_radii_set(SasaRadii radii_set){
	implicit_radii_set_ = radii_set;
}

void
SasaCalc::set_exclude_polar_atoms_by_charge(bool exclude_polar_all, core::Real charge_cutoff) {
	exclude_polar_all_in_hsasa_ = exclude_polar_all;
	polar_charge_cutoff_ = charge_cutoff;
}

void
SasaCalc::set_include_carbon_sulfer_only_in_hydrophobic_calc(bool include_c_s_only) {
	include_c_s_only_in_hsasa_ = include_c_s_only;
}

//SasaCalc::set_expand_polar_radii(bool expand_polars, core::Size expansion_radius /*1.0*/) {
//	expand_polars_ = expand_polars;
//	polar_expansion_radius_ = expansion_radius;
//}

void
SasaCalc::set_use_big_polar_hydrogen(bool big_polar_h) {
	big_polar_h_ = big_polar_h;
}

void
SasaCalc::fill_all_data(Real& total_hsasa, id::AtomID_Map<Real>& atom_sasa, vector1<Real>& rsd_sasa, vector1<Real>& rsd_hsasa, vector1<Real>& rel_hsasa){
	total_hsasa = total_hsasa_;
	atom_sasa = atom_sasa_;
	rsd_sasa = rsd_sasa_;
	rsd_hsasa = rsd_hsasa_;
	rel_hsasa = rel_hydrophobic_sasa_by_charge_;
}

Real
SasaCalc::calculate(const pose::Pose& pose) {
	
	
	
	//Figure out which radii to use.
	SasaRadii radii_set;
	if (legacy_defaults_){
		TR << "Using legacy default method, radii, and atom subset.  Note that hSASA will be wrong. " << std::endl;
		include_hydrogens_ = true;
		radii_set = legacy;
	}
	else if (include_hydrogens_){
		radii_set = explicit_radii_set_;
	}
	else{
		radii_set = implicit_radii_set_;
	}
	
	
	
	init(pose);
	setup_sasa_method(radii_set);
	
	total_sasa_ = method_->calculate(pose, atom_subset_, atom_sasa_, rsd_sasa_);
	
	calc_per_res_hphobic_sasa(pose);
	
	return total_sasa_;
}


void
SasaCalc::init(const pose::Pose& pose) {
	
	atom_sasa_.clear();
	atom_subset_.clear();
	
	rsd_sasa_.clear();
	rsd_hsasa_.clear();
	rel_hydrophobic_sasa_by_charge_.clear();
	if (include_hydrogens_){
		core::pose::initialize_atomid_map(atom_sasa_, pose, 0.0);
		core::pose::initialize_atomid_map(atom_subset_, pose, true);
	}
	else {
		core::pose::initialize_atomid_map(atom_sasa_, pose, 0.0);
		core::pose::initialize_atomid_map(atom_subset_, pose, false);
		for (core::Size i = 1; i <= pose.total_residue(); ++i){
			for (core::Size x = 1; x <= pose.residue(i).nheavyatoms(); ++x){
				core::id::AtomID atom_id(x, i);
				atom_subset_[atom_id] = true;
			}
		}
	}
	
	rsd_sasa_.resize(pose.total_residue(), 0.0);
	rsd_hsasa_.resize(pose.total_residue(), 0.0);
	rel_hydrophobic_sasa_by_charge_.resize(pose.total_residue(), 0.0);
	
	total_sasa_ = 0;
	total_hsasa_ = 0;
	total_rel_hsasa_ = 0;
}

void
SasaCalc::calc_per_res_hphobic_sasa(const pose::Pose & pose) {
	
	for (Size i = 1; i <= atom_sasa_.n_residue(); ++i){
		if (rsd_sasa_[i] == 0.0) {
			rsd_hsasa_[i] = 0.0;
			rel_hydrophobic_sasa_by_charge_[i] = 0.0;
		}
		
		core::conformation::Residue const & res = pose.residue(i);

		core::Real sasa = 0.0;
		core::Real hsasa = 0.0;
		core::Real hsasa_rel = 0.0;
		
		for (Size x = 1; x <= atom_sasa_.n_atom(i); ++x) {
			
			core::id::AtomID atomid(x, i);
			
			if (! atom_subset_[atomid]) continue;

			//Options do not make sense together and will give wrong results.
			if (!include_c_s_only_in_hsasa_ && !exclude_polar_all_in_hsasa_){
				utility_exit_with_message("You cannot include atoms other than C and S in calculation without excluding them by charge.  All atoms would therefore be included in hSASA");
			}
			
			//Only include C or S in calculation if option is given
			if (include_c_s_only_in_hsasa_ && !(res.atom_type( x ).element() == "C" || res.atom_type( x ).element() == "S")) {
				continue;
			}
			
			
			//If we are excluding polars, skip the polar from inclusion into hsasa if charge is too great.
			if (exclude_polar_all_in_hsasa_   && std::abs(res.type().atom(x).charge()) > polar_charge_cutoff_){

				continue;
			}
			

			hsasa = hsasa + atom_sasa_[atomid];
			
			//If we are only including C and S; only look at hydrogens bound to them to include in hSASA.  Otherwise skip this or we will double count hydrogen contributions.
			//Skip this if our atom_subset does not already have hydrogens.
			if (include_hydrogens_ && res.atom_type(x).element() != "H" && include_c_s_only_in_hsasa_){
				utility::vector1< core::Size> bonded_indices = res.type().bonded_neighbor(x);
				for (core::Size nb = 1; nb <= bonded_indices.size(); ++nb){
					core::Size x_bonded = bonded_indices[nb];


					if (res.atom_type(x_bonded).element() == "H"){
						core::id::AtomID bonded_id(x_bonded, i);
						hsasa = hsasa + atom_sasa_[bonded_id];
					}
				} // for bonded atom neighbors
			} 
				
			
			core::Real charge = res.type().atom(x).charge();
			hsasa_rel = hsasa_rel + atom_sasa_[atomid] * (1- std::abs(charge));
			
		} //for atom
		
		rsd_sasa_[i] = sasa;
		rsd_hsasa_[i] = hsasa;
		rel_hydrophobic_sasa_by_charge_[i] = hsasa_rel;
		total_hsasa_ = total_hsasa_ + hsasa;
		
	} //for residue
}




//SasaCalc::calculate(const pose::Pose& pose, id::AtomID_Map<Real>& atom_sasa) {
	
//}

//SasaCalc::calculate(const pose::Pose& pose, vector1<Real>& rsd_sasa, vector1<Real>& rsd_hsasa) {
	
//}

//SasaCalc::calculate(const pose::Pose& pose, id::AtomID_Map<Real>& atom_sasa, vector1<Real>& rsd_sasa, vector1<Real>& rsd_hsasa){
	
//}


}
}
}
