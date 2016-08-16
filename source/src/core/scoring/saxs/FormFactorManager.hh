// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/saxs/FormFactorManager.hh
/// @brief Loads and manages atomic form factors
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_core_scoring_saxs_FormFactorManager_hh
#define INCLUDED_core_scoring_saxs_FormFactorManager_hh

#include <core/scoring/saxs/FormFactorManager.fwd.hh>

#include <core/scoring/saxs/FormFactor.hh>

// utility headers
#include <core/types.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <map>
#include <string>

namespace core {
namespace scoring {
namespace saxs {

/// @brief selects a given number of fragments using a quota scheme
class FormFactorManager {
public:

	/// @brief return singleton of the manager
	static FormFactorManager* get_manager();

	void register_ff(std::string atom_name,FormFactorOP new_ff);

	void load_ff(std::string config_file);

	void load_ff_from_db(std::string file_name);

	/// @brief returns true if the manager has form factor function for a given atom
	bool is_known_atom(std::string atom_name);

	/// @brief returns form factor function for a given atom index
	FormFactorOP get_ff(Size atom_id) { return ff_vector_[atom_id]; }

	/// @brief returns the number of form factors registered in this manager
	Size count_ff() { return ff_vector_.size(); }

	/// @brief returns form factor function for a given atom
	FormFactorOP get_ff(std::string atom_name);

	/// @brief returns a vector of know atom names
	utility::vector1<std::string> get_known_atoms() {

		return known_atoms_;
	}

	/// @brief returns an index of an atom type or 0 if teh atom is not registered
	Size get_atom_index(std::string atom_name) {

		if ( names_to_indexes_.find(atom_name) != names_to_indexes_.end() ) {
			return names_to_indexes_.find(atom_name)->second;
		} else {
			return 0;
		}
	}

	/// @brief asks all the registered form factors to tabulate their values for the new vector of q-arguments
	void tabulate(const utility::vector1<Real> & q) {

		for ( Size i=1; i<=ff_vector_.size(); ++i )  ff_vector_[i]->tabulate(q);
	}

private:
	/// @brief all known FF are stored here
	std::map<std::string,FormFactorOP> ff_map_;
	/// @brief vector of all known form factor, stored as OPs
	utility::vector1<FormFactorOP> ff_vector_;
	/// @brief all known atom names
	utility::vector1<std::string> known_atoms_;
	/// @brief maps all known atom names to integer indexes in ff_vector
	std::map<std::string,Size> names_to_indexes_;

	static FormFactorManager* singleton_;

	/// @brief  Constructor reads basic form factors from library
	FormFactorManager() {}
};


} // core
} // scoring
} // saxs

#endif /* INCLUDED_core_scoring_saxs_FormFactorManager_HH */
