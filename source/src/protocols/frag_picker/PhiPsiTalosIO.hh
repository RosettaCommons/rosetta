// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

using namespace core;

class PhiPsiTalosIO {
public:

	PhiPsiTalosIO() {
		first_residue_index_ = 1;
	}

	PhiPsiTalosIO(std::string file_name) {
		first_residue_index_ = 1;
		read(file_name);
	}

	Size get_first_residue_id() const {
		return first_residue_index_;
	}

	Size get_last_residue_id() const {
		return first_residue_index_ + sequence_.length() - 1;
	}

	inline
	std::string const &
	get_sequence() const {
		return sequence_;
	}

	inline
	bool
	has_entry(Size residue_id) const {
		if ( entries_.find(residue_id) != entries_.end() ) {
			return true;
		}
		return false;
	}

	inline const boost::tuple<Size, char, Real, Real, Real, Real, Real, Real,
	Size, std::string> get_entry(const Size res_id) {
		return entries_.find(res_id)->second;
	}
	inline Real phi(Size res_id) {
		return entries_.find(res_id)->second.get<2> ();
	}

	inline Real psi(Size res_id) {
		return entries_.find(res_id)->second.get<3> ();
	}

	inline Real d_phi(Size res_id) {
		return entries_.find(res_id)->second.get<4> ();
	}

	inline Real d_psi(Size res_id) {
		return entries_.find(res_id)->second.get<5> ();
	}

	inline Real dist(Size res_id) {
		return entries_.find(res_id)->second.get<6> ();
	}

	inline Real s2(Size res_id) {
		return entries_.find(res_id)->second.get<7> ();
	}

	inline std::string quality(Size res_id) {
		return entries_.find(res_id)->second.get<9> ();
	}

	void write(std::ostream&);
	void read(std::string const&);

private:
	std::string data_format_;
	std::string sequence_;
	utility::vector1<std::string> column_names_;
	Size first_residue_index_;
	Size last_residue_index_;
	std::map<Size, boost::tuple<Size, char, Real, Real, Real, Real, Real, Real,
		Size, std::string> > entries_;
	void print_sequence(std::string const &, std::ostream &) const;
};

} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_PhiPsiTalosIO_HH */

