// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/frag_picker/CSTalosIO.hh
/// @brief Class that reads and writes chemical shifts in TALOS format
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_CSTalosIO_hh
#define INCLUDED_protocols_frag_picker_CSTalosIO_hh

// utility headers
#include <core/types.hh>

#include <string>
#include <map>

// boost headers
#include <boost/tuple/tuple.hpp>
// AUTO-REMOVED #include <utility/vector1.hh>

#include <utility/vector1_bool.hh>


namespace protocols {
namespace frag_picker {

using namespace core;

class CSTalosIO {
public:

	CSTalosIO() {
		first_residue_index_ = 1;
		set_up_atom_order();
	}

	CSTalosIO(std::string file_name) {
		first_residue_index_ = 1;
		set_up_atom_order();
		read(file_name);
	}

	utility::vector1<utility::vector1<std::pair<Size, Real> > >
			repack_to_matrix();

	Size get_first_residue_id() const {
		return first_residue_index_;
	}

	Size get_last_residue_id() const {
		return first_residue_index_ + sequence_.length() - 1;
	}

	std::string get_sequence() const {
		return sequence_;
	}

	bool has_entry(Size residue_id) {
		if (resids_to_entries_map_.find(residue_id)
				!= resids_to_entries_map_.end())
			return true;
		return false;
	}

	void get_tuples(Size, utility::vector1<boost::tuple<Size, char,
			std::string, Real> >);

	void get_tuples(Size, utility::vector1<boost::tuple<Size, char,
			std::string, Real> > &) const;

	utility::vector1<boost::tuple<Size, char, std::string, Real> > get_entries() {
		return entries_;
	}

	void write(std::ostream&);
	void read(std::string const&);
	Real get_shift(Size, std::string const &) const;
	bool has_atom(Size, std::string const &) const;

private:
	std::string data_format_;
	std::string sequence_;
	utility::vector1<std::string> column_names_;
	Size first_residue_index_;
	utility::vector1<boost::tuple<Size, char, std::string, Real> > entries_;
	std::multimap<Size, Size> resids_to_entries_map_;
	std::map<std::string, Size> order_of_atoms_;
	void print_sequence(std::string const &, std::ostream &) const;
	utility::vector1<boost::tuple<Size, char, std::string, Real> >
			used_for_searching_;
	void set_up_atom_order();
};

} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_CSTalosIO_HH */

