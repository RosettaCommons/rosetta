// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/frag_picker/PhiPsiTalosIO.cc
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

// unit headers
#include <protocols/frag_picker/PhiPsiTalosIO.hh>

// project headers

#include <core/types.hh>
#include <basic/Tracer.hh>

// utility headers
#include <utility/exit.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>

// C++ headers
#include <string>
#include <map>

// boost headers
#include <boost/tuple/tuple.hpp>
// AUTO-REMOVED #include <boost/algorithm/string.hpp>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {

static basic::Tracer
		trPhiPsiTalosIO("protocols.frag_picker.PhiPsiTalosIO");

void PhiPsiTalosIO::read(std::string const & file_name) {

	utility::io::izstream data(file_name.c_str());
	trPhiPsiTalosIO.Info << "read talos data from " << file_name << std::endl;
	if (!data)
		utility_exit_with_message("[ERROR] Unable to open talos file: "
				+ file_name);

	std::string line;
	last_residue_index_ = 0;
	first_residue_index_ = 0;
	bool first_not_found = true;
	while (!data.eof()) {
		getline(data, line);
		std::istringstream line_stream(line);
		utility::vector1<std::string> strs;
		while (!line_stream.eof()) {
			std::string token;
			line_stream >> token;
			strs.push_back(token);
		}
		if (strs[1] == "DATA") {
			if (strs[2] == "SEQUENCE") {
				for (Size i = 3; i <= strs.size(); ++i) {
					sequence_ += strs[i];
				}
			} else if (strs[2] == "FIRST_RESID") {
				first_not_found = false;
			} else {
				trPhiPsiTalosIO.Warning << "Unrecognized DATA entry:" << line
						<< std::endl;
			}
		}
		if (strs[1] == "VARS") {
			for (Size i = 2; i <= strs.size(); ++i) {
				sequence_ += strs[i];
			}

		}
		if (strs[1] == "FORMAT") {
			data_format_ = line.substr(7);
		}

		if ((strs.size() == 10) && (strs[1] != "REMARK")) {
			char aa;
			std::istringstream line_stream(line);
			Real phi, psi, d_phi, d_psi, dist, s2;
			Size res_id, count;
			std::string cls;
			line_stream >> res_id >> aa >> phi >> psi >> d_phi >> d_psi >> dist
					>> s2 >> count >> cls;
			sequence_ += aa;

			boost::tuple<Size, char, Real, Real, Real, Real, Real, Real, Size,
					std::string> t(res_id, aa, phi, psi, d_phi, d_psi, dist,
					s2, count, cls);
			entries_.insert(std::pair<Size, boost::tuple<Size, char, Real,
					Real, Real, Real, Real, Real, Size, std::string> >(res_id,
					t));
			if (last_residue_index_ < res_id)
				last_residue_index_ = res_id;
		}
	}

	if (first_not_found)
		trPhiPsiTalosIO.Warning
				<< "FIRST_RESID keyword didn't show up in a file header\n\tAssuming the first residue id is 1"
				<< std::endl;
	if (sequence_.length() == 0) {
		for (Size i = first_residue_index_; i <= last_residue_index_; ++i) {
			if (has_entry(i))
				sequence_ += entries_.find(i)->second.get<1> ();
			else
				sequence_ += 'X';
		}
		trPhiPsiTalosIO.Warning
				<< "Could not find a SEQUENCE data in the input file\nSequence based on entires is:\n"
				<< sequence_ << std::endl;
	}
}

void PhiPsiTalosIO::write(std::ostream& out) {

	out << "DATA FIRST_RESID " << first_residue_index_ << "\n";
	out << "DATA SEQUENCE " << sequence_ << "\n";
	out << "VARS";
	for (Size i = 1; i <= column_names_.size(); i++)
		out << " " << column_names_[i];
	out << "\n" << "FORMAT " << data_format_ << std::endl << "\n";
	char buffer[100];
	char c1[2];
	c1[1] = 0;
	for (Size i = 1; i <= entries_.size(); i++) {
		c1[0] = entries_[i].get<1> ();
		sprintf(buffer, data_format_.c_str(), entries_[i].get<0> (), c1,
				entries_[i].get<2> (), entries_[i].get<3> (),
				entries_[i].get<4> (), entries_[i].get<5> (),
				entries_[i].get<6> (), entries_[i].get<7> (),
				entries_[i].get<8> (), entries_[i].get<9> ().c_str());
		out << buffer << "\n";
	}
	out << std::endl;
}

} // frag_picker
} // protocols
