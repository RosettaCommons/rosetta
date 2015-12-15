// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/core/scoring/sasa/SasaCalc.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/chemical/AtomType.hh>
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/scoring/sasa/SasaMethodFactory.hh>
#include <core/scoring/sasa/util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/sasa.OptionKeys.gen.hh>

#include <core/pose/util.tmpl.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.sasa.SasaCalc" );

#ifdef    SERIALIZATION
// Project serialization headers
#include <core/id/AtomID_Map.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

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
	set_defaults();
	method_type_ = method;
}


SasaCalc::~SasaCalc(){}


void
SasaCalc::set_defaults() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	set_calculation_method(get_sasa_method_from_string(option[OptionKeys::sasa::method]()));

	set_probe_radius(option[OptionKeys::sasa::probe_radius]());
	set_include_hydrogens_explicitly(option[OptionKeys::sasa::include_hydrogens_explicitly]());
	set_include_probe_radius_in_atom_radii(option[OptionKeys::sasa::include_probe_radius_in_atom_radii]());
	set_include_carbon_sulfer_only_in_hydrophobic_calc(option[OptionKeys::sasa::include_only_C_S_in_hsasa]());
	set_exclude_polar_atoms_by_charge(option[OptionKeys::sasa::exclude_polar_atoms_by_charge_in_hsasa]());
	set_polar_charge_cutoff(option[OptionKeys::sasa::polar_charge_cutoff]());
	set_implicit_hydrogen_included_radii_set(get_sasa_radii_set_from_string(option[OptionKeys::sasa::implicit_hydrogen_radii_set]()));
	set_explicit_hydrogen_included_radii_set(get_sasa_radii_set_from_string(option[OptionKeys::sasa::explicit_hydrogen_radii_set]()));

	//set_expand_polar_radii(false);
	set_use_big_polar_hydrogen(false);


}


void
SasaCalc::set_include_hydrogens_explicitly(bool include_hydrogens) {
	include_hydrogens_ = include_hydrogens;
}

void
SasaCalc::set_calculation_method(SasaMethodEnum method) {
	method_type_ = method;
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
SasaCalc::set_explicit_hydrogen_included_radii_set(SasaRadii radii_set) {
	explicit_radii_set_ = radii_set;
}
void
SasaCalc::set_exclude_polar_atoms_by_charge(bool exclude_polar_all) {
	exclude_polar_all_in_hsasa_ = exclude_polar_all;
}

void
SasaCalc::set_include_carbon_sulfer_only_in_hydrophobic_calc(bool include_c_s_only) {
	include_c_s_only_in_hsasa_ = include_c_s_only;
}

//SasaCalc::set_expand_polar_radii(bool expand_polars, core::Size expansion_radius /*1.0*/) {
// expand_polars_ = expand_polars;
// polar_expansion_radius_ = expansion_radius;
//}

void
SasaCalc::set_use_big_polar_hydrogen(bool big_polar_h) {
	big_polar_h_ = big_polar_h;
}

vector1<Real>
SasaCalc::get_residue_sasa_bb() const {

	vector1<Real> sasa_bb;
	for ( Size i = 1; i <= rsd_sasa_sc_.size(); ++i ) {
		sasa_bb[i] = rsd_sasa_[i] - rsd_sasa_sc_[i];
	}
	return sasa_bb;
}

vector1<Real>
SasaCalc::get_residue_hsasa_bb() const {

	vector1< Real > sasa_bb;
	for ( Size i = 1; i <= rsd_hsasa_sc_.size(); ++i ) {
		sasa_bb[i] = rsd_hsasa_[i] - rsd_hsasa_sc_[i];
	}
	return sasa_bb;
}
void
SasaCalc::fill_data(Real& total_hsasa,Real & total_rel_hsasa,  id::AtomID_Map<Real>& atom_sasa, vector1<Real>& rsd_sasa, vector1<Real>& rsd_hsasa, vector1<Real>& rel_hsasa){
	total_hsasa = total_hsasa_;
	total_rel_hsasa = total_rel_hsasa_;
	atom_sasa = atom_sasa_;
	rsd_sasa = rsd_sasa_;
	rsd_hsasa = rsd_hsasa_;
	rel_hsasa = rel_hydrophobic_sasa_by_charge_;
}

///// Legacy-style interfaces //////////
Real
SasaCalc::calculate(const pose::Pose& pose, id::AtomID_Map<Real>& atom_sasa, vector1<Real>& rsd_sasa, vector1<Real>& rsd_hsasa, vector1<Real>& rsd_rel_hsasa) {
	Real total_sasa = calculate(pose, atom_sasa, rsd_sasa, rsd_hsasa);
	rsd_rel_hsasa = rel_hydrophobic_sasa_by_charge_;
	return total_sasa;
}

Real
SasaCalc::calculate(const pose::Pose& pose, id::AtomID_Map<Real>& atom_sasa, vector1<Real>& rsd_sasa, vector1<Real>& rsd_hsasa) {
	Real total_sasa = calculate(pose, rsd_sasa, rsd_hsasa);
	atom_sasa = atom_sasa_;
	return total_sasa;

}

Real
SasaCalc::calculate(const pose::Pose& pose, vector1<Real>& rsd_sasa, vector1<Real>& rsd_hsasa) {
	Real total_sasa = calculate(pose);
	rsd_sasa = rsd_sasa_;
	rsd_hsasa = rsd_hsasa_;
	return total_sasa;
}

Real
SasaCalc::calculate(const pose::Pose& pose, id::AtomID_Map<Real>& atom_sasa){
	Real total_sasa = calculate(pose);
	atom_sasa = atom_sasa_;
	return total_sasa;
}

///// //////////

Real
SasaCalc::calculate(const pose::Pose& pose) {



	//Figure out which radii to use.
	SasaRadii radii_set;
	if ( include_hydrogens_ ) {
		radii_set = explicit_radii_set_;
	} else {
		radii_set = implicit_radii_set_;
	}



	init(pose);
	setup_sasa_method(radii_set);

	total_sasa_ = method_->calculate(pose, atom_subset_, atom_sasa_, rsd_sasa_);

	calc_per_res_sasas(pose);

	return total_sasa_;
}


void
SasaCalc::init(const pose::Pose& pose) {

	atom_sasa_.clear();
	atom_subset_.clear();

	rsd_sasa_.clear();
	rsd_sasa_sc_.clear();

	rsd_hsasa_.clear();
	rsd_hsasa_sc_.clear();

	rel_hydrophobic_sasa_by_charge_.clear();

	if ( include_hydrogens_ ) {
		core::pose::initialize_atomid_map(atom_sasa_, pose, 0.0);
		core::pose::initialize_atomid_map(atom_subset_, pose, true);
	} else {
		core::pose::initialize_atomid_map(atom_sasa_, pose, 0.0);
		core::pose::initialize_atomid_map(atom_subset_, pose, false);
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			for ( core::Size x = 1; x <= pose.residue(i).nheavyatoms(); ++x ) {
				core::id::AtomID atom_id(x, i);
				atom_subset_[atom_id] = true;
			}
		}
	}

	rsd_sasa_.resize(pose.total_residue(), 0.0);
	rsd_sasa_sc_.resize(pose.total_residue(), 0.0);

	rsd_hsasa_.resize(pose.total_residue(), 0.0);
	rsd_hsasa_sc_.resize(pose.total_residue(), 0.0);

	rel_hydrophobic_sasa_by_charge_.resize(pose.total_residue(), 0.0);

	total_sasa_ = 0;
	total_sasa_sc_ = 0;

	total_hsasa_ = 0;
	total_hsasa_sc_ = 0;

	total_rel_hsasa_ = 0;
}

void
SasaCalc::calc_per_res_sasas(const pose::Pose & pose) {

	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( rsd_sasa_[i] == 0.0 ) {
			continue;
		}

		core::conformation::Residue const & res = pose.residue(i);

		core::Real sasa_sc = 0.0;

		core::Real hsasa = 0.0;
		core::Real hsasa_sc = 0.0;

		core::Real hsasa_rel = 0.0;

		for ( Size x = 1; x <= pose.residue(i).natoms(); ++x ) {

			core::id::AtomID atomid(x, i);
			bool is_sc = ! pose.residue(i).atom_is_backbone(x);

			if ( ! atom_subset_[atomid] || atom_sasa_[atomid] <= 0.0 ) continue;

			if ( is_sc ) sasa_sc += atom_sasa_[atomid];

			core::Real charge = res.type().atom(x).charge();
			hsasa_rel = hsasa_rel + atom_sasa_[atomid] * (1- std::abs(charge));

			//Options do not make sense together and will give wrong results.
			if ( !include_c_s_only_in_hsasa_ && !exclude_polar_all_in_hsasa_ ) {
				utility_exit_with_message("You cannot include atoms other than C and S in calculation without excluding them by charge.  All atoms would therefore be included in hSASA");
			}

			//Only include C or S in calculation if option is given
			if ( include_c_s_only_in_hsasa_ && !(res.atom_type( x ).element() == "C" || res.atom_type( x ).element() == "S") ) {
				continue;
			}


			//If we are excluding polars, skip the polar from inclusion into hsasa if charge is too great.
			if ( exclude_polar_all_in_hsasa_   && std::abs(res.type().atom(x).charge()) > polar_charge_cutoff_ ) {

				continue;
			}


			hsasa = hsasa + atom_sasa_[atomid];
			if ( is_sc ) hsasa_sc = hsasa_sc + atom_sasa_[atomid];

			//If we are only including C and S; only look at hydrogens bound to them to include in hSASA.  Otherwise skip this or we will double count hydrogen contributions.
			//Skip this if our atom_subset does not already have hydrogens.
			if ( include_hydrogens_ && res.atom_type(x).element() != "H" && include_c_s_only_in_hsasa_ ) {
				utility::vector1< core::Size> bonded_indices = res.type().bonded_neighbor(x);
				for ( core::Size nb = 1; nb <= bonded_indices.size(); ++nb ) {
					core::Size x_bonded = bonded_indices[nb];


					if ( res.atom_type(x_bonded).element() == "H" ) {
						core::id::AtomID bonded_id(x_bonded, i);
						hsasa = hsasa + atom_sasa_[bonded_id];
					}
				} // for bonded atom neighbors
			}

		} //for atom

		//TR << sasa<<" " << hsasa << " " << hsasa_rel << std::endl;

		rsd_sasa_sc_[i] = sasa_sc;

		rsd_hsasa_[i] = hsasa;
		rsd_hsasa_sc_[i] = hsasa_sc;

		rel_hydrophobic_sasa_by_charge_[i] = hsasa_rel;

		total_sasa_sc_ = total_sasa_sc_ + sasa_sc;
		total_rel_hsasa_= total_rel_hsasa_ + hsasa_rel;
		total_hsasa_  = total_hsasa_ + hsasa;
		total_hsasa_sc_ = total_hsasa_sc_ + hsasa_sc;

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

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::sasa::SasaCalc::save( Archive & arc ) const {
	arc( CEREAL_NVP( method_type_ ) ); // enum core::scoring::sasa::SasaMethodEnum
	arc( CEREAL_NVP( method_ ) ); // SasaMethodOP
	arc( CEREAL_NVP( implicit_radii_set_ ) ); // enum core::scoring::sasa::SasaRadii
	arc( CEREAL_NVP( explicit_radii_set_ ) ); // enum core::scoring::sasa::SasaRadii
	arc( CEREAL_NVP( include_hydrogens_ ) ); // _Bool
	arc( CEREAL_NVP( probe_radius_ ) ); // Real
	arc( CEREAL_NVP( include_c_s_only_in_hsasa_ ) ); // _Bool
	arc( CEREAL_NVP( exclude_polar_all_in_hsasa_ ) ); // _Bool
	arc( CEREAL_NVP( polar_charge_cutoff_ ) ); // Real
	arc( CEREAL_NVP( include_probe_radius_ ) ); // _Bool
	arc( CEREAL_NVP( big_polar_h_ ) ); // _Bool
	arc( CEREAL_NVP( atom_subset_ ) ); // id::AtomID_Map<_Bool>
	arc( CEREAL_NVP( atom_sasa_ ) ); // id::AtomID_Map<Real>
	arc( CEREAL_NVP( rsd_sasa_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( rsd_sasa_sc_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( rsd_hsasa_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( rsd_hsasa_sc_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( rel_hydrophobic_sasa_by_charge_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( total_sasa_ ) ); // Real
	arc( CEREAL_NVP( total_sasa_sc_ ) ); // Real
	arc( CEREAL_NVP( total_hsasa_ ) ); // Real
	arc( CEREAL_NVP( total_hsasa_sc_ ) ); // Real
	arc( CEREAL_NVP( total_rel_hsasa_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::sasa::SasaCalc::load( Archive & arc ) {
	arc( method_type_ ); // enum core::scoring::sasa::SasaMethodEnum
	arc( method_ ); // SasaMethodOP
	arc( implicit_radii_set_ ); // enum core::scoring::sasa::SasaRadii
	arc( explicit_radii_set_ ); // enum core::scoring::sasa::SasaRadii
	arc( include_hydrogens_ ); // _Bool
	arc( probe_radius_ ); // Real
	arc( include_c_s_only_in_hsasa_ ); // _Bool
	arc( exclude_polar_all_in_hsasa_ ); // _Bool
	arc( polar_charge_cutoff_ ); // Real
	arc( include_probe_radius_ ); // _Bool
	arc( big_polar_h_ ); // _Bool
	arc( atom_subset_ ); // id::AtomID_Map<_Bool>
	arc( atom_sasa_ ); // id::AtomID_Map<Real>
	arc( rsd_sasa_ ); // utility::vector1<Real>
	arc( rsd_sasa_sc_ ); // utility::vector1<Real>
	arc( rsd_hsasa_ ); // utility::vector1<Real>
	arc( rsd_hsasa_sc_ ); // utility::vector1<Real>
	arc( rel_hydrophobic_sasa_by_charge_ ); // utility::vector1<Real>
	arc( total_sasa_ ); // Real
	arc( total_sasa_sc_ ); // Real
	arc( total_hsasa_ ); // Real
	arc( total_hsasa_sc_ ); // Real
	arc( total_rel_hsasa_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::sasa::SasaCalc );
CEREAL_REGISTER_TYPE( core::scoring::sasa::SasaCalc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_sasa_SasaCalc )
#endif // SERIALIZATION
