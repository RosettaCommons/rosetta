// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/TorsionSampler.hh
///
/// @brief Sample torsions based on experimental distributions
/// @author Guangfeng Zhou and Frank DiMaio

#ifndef INCLUDED_protocols_ligand_docking_GALigandDock_TorsionSampler_hh
#define INCLUDED_protocols_ligand_docking_GALigandDock_TorsionSampler_hh

#include <unordered_map>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <core/chemical/Bond.fwd.hh>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

class TorsionDistrParams{
public:
	TorsionDistrParams();
	TorsionDistrParams(std::string const& ttype_in, core::Size mult_in);
	void add_mode( core::Real A, core::Real mu, core::Real sigma);
	void set_torsion_type(std::string const& ttype_in) { torsion_type_ = ttype_in; }
	std::string torsion_type() const { return torsion_type_; }
	core::Size n_mode() const { return n_mode_; }
	void set_multiplicity( core::Real mult_in ) { mult_ = mult_in; }
	core::Size multiplicity() const { return mult_; }
	void get_mode_i(core::Size i, core::Real& A_out,
		core::Real& mu_out, core::Real& sigma_out ) const;
	void normalize();
	void sample_mode (core::Real& A_out,
		core::Real& mu_out, core::Real& sigma_out) const;
private:
	core::Size n_mode_;
	core::Size mult_;
	bool normalized_;
	std::string torsion_type_;
	utility::vector1< core::Real > A_vec_, mu_vec_, sigma_vec_, A_cumulative_vec_;

};

class TorsionSampler{
public:
	TorsionSampler();

	core::Real
	sample(core::chemical::BondName bn, core::chemical::BondRingness br,
		core::Size type1, core::Size type2, core::Size type3, core::Size type4) const;

	TorsionDistrParams const &
	lookup_tors_distr_params( core::chemical::BondName bn, core::chemical::BondRingness br,
		core::Size type1, core::Size type2, core::Size type3, core::Size type4
	) const;


private:
	void read_database( std::string const & filename);

private:
	std::unordered_map< uint64_t, core::Size > tors_lookup_;
	std::unordered_map< std::string, utility::vector1<core::Size> > name_index_map;
	utility::vector1< bool > defined_atom_types_;
	utility::vector1< TorsionDistrParams > tors_distr_;
	static const TorsionDistrParams null_tors_distr;
};




}
}
}

#endif
