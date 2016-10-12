// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/PhiPsiTalosIO.hh
/// @brief Class that reads and writes Phi-Psi in TALOS format
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_PhiPsiTalosIO_hh
#define INCLUDED_protocols_frag_picker_PhiPsiTalosIO_hh

// utility headers
#include <core/types.hh>

#include <string>
#include <map>


// boost headers
#include <boost/tuple/tuple.hpp>

#include <utility/vector1_bool.hh>


namespace protocols {
namespace frag_picker {

class PhiPsiTalosIO {
public:

	PhiPsiTalosIO() {
		first_residue_index_ = 1;
	}

	PhiPsiTalosIO(std::string file_name) {
		first_residue_index_ = 1;
		read(file_name);
	}

	core::Size get_first_residue_id() const {
		return first_residue_index_;
	}

	core::Size get_last_residue_id() const {
		return first_residue_index_ + sequence_.length() - 1;
	}

	inline
	std::string const &
	get_sequence() const {
		return sequence_;
	}

	inline
	bool
	has_entry(core::Size residue_id) const {
		if ( entries_.find(residue_id) != entries_.end() ) {
			return true;
		}
		return false;
	}

	inline const boost::tuple<core::Size, char, core::Real, core::Real, core::Real, core::Real, core::Real, core::Real,
	core::Size, std::string> get_entry(const core::Size res_id) {
		return entries_.find(res_id)->second;
	}
	inline core::Real phi(core::Size res_id) {
		return entries_.find(res_id)->second.get<2> ();
	}

	inline core::Real psi(core::Size res_id) {
		return entries_.find(res_id)->second.get<3> ();
	}

	inline core::Real d_phi(core::Size res_id) {
		return entries_.find(res_id)->second.get<4> ();
	}

	inline core::Real d_psi(core::Size res_id) {
		return entries_.find(res_id)->second.get<5> ();
	}

	inline core::Real dist(core::Size res_id) {
		return entries_.find(res_id)->second.get<6> ();
	}

	inline core::Real s2(core::Size res_id) {
		return entries_.find(res_id)->second.get<7> ();
	}

	inline std::string quality(core::Size res_id) {
		return entries_.find(res_id)->second.get<9> ();
	}

	void write(std::ostream&);
	void read(std::string const&);

private:
	std::string data_format_;
	std::string sequence_;
	utility::vector1<std::string> column_names_;
	core::Size first_residue_index_;
	core::Size last_residue_index_;
	std::map<core::Size, boost::tuple<core::Size, char, core::Real, core::Real, core::Real, core::Real, core::Real, core::Real,
		core::Size, std::string> > entries_;
	void print_sequence(std::string const &, std::ostream &) const;
};

} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_PhiPsiTalosIO_HH */

