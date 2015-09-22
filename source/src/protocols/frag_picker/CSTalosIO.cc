// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/CSTalosIO.cc
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

// unit headers
#include <protocols/frag_picker/CSTalosIO.hh>

// project headers
#include <core/types.hh>

#include <basic/Tracer.hh>

// utility headers
#include <utility/io/izstream.hh>

#include <string>
#include <map>
#include <utility/exit.hh>

// boost headers
#include <boost/tuple/tuple.hpp>
#include <utility/vector1.hh>

namespace protocols {
namespace frag_picker {

using namespace core;

static THREAD_LOCAL basic::Tracer tr( "protocols.frag_picker.TalosReader" );

utility::vector1<utility::vector1<std::pair<Size, Real> > > CSTalosIO::repack_to_matrix() {

	Size len = get_sequence().length();

	utility::vector1<utility::vector1<std::pair<Size, Real> > > data(len);
	for ( Size i = 1; i <= entries_.size(); ++i ) {
		// Size res_id = entries_[i].get<0> ();
		std::string& at_name = entries_[i].get<2> ();
		Real shift = entries_[i].get<3> ();
		Size atom_id = order_of_atoms_.find(at_name)->second;

		data[i].push_back(std::pair<Size, Real>(atom_id, shift));
	}

	return data;
}

void CSTalosIO::read(std::string const & file_name) {

	utility::io::izstream data(file_name.c_str());
	tr.Info << "read talos data from " << file_name << std::endl;
	if ( !data ) {
		utility_exit_with_message("[ERROR] Unable to open talos file: "
			+ file_name);
	}
	std::string line;
	std::string keyword;
	std::string subkeyword;
	std::string entry;
	sequence_ = "";
	bool header_done = false;
	bool first_not_found = false;
	while ( !header_done ) {
		getline(data, line);
		std::istringstream line_stream(line);
		line_stream >> keyword;
		if ( keyword == "DATA" ) {
			line_stream >> subkeyword;
			if ( subkeyword == "SEQUENCE" ) {
				while ( line_stream >> entry ) {
					sequence_ += entry;
				}
			} else if ( subkeyword == "FIRST_RESID" ) {
				first_not_found = true;
				line_stream >> first_residue_index_;
				tr.Info << "first residue is set to " << first_residue_index_ << " by FIRST_RESID entry in chemical shift file." << std::endl;
			} else {
				tr.Warning << "Unrecognized DATA entry:" << line
					<< std::endl;
			}
		}
		if ( keyword == "VARS" ) {
			while ( line_stream >> entry ) {
				column_names_.push_back(entry);
			}
		}
		if ( keyword == "FORMAT" ) {
			line_stream >> entry;
			data_format_ = line.substr(7);
			header_done = true;
		}
	}
	if ( !first_not_found ) {
		//this is quite normal...
		tr.Debug << "FIRST_RESID keyword didn't show up in a file header\n\tAssuming the first residue id is 1"
			<< std::endl;
	}

	while ( getline(data, line) ) {
		if ( line.length() > 7 ) {
			std::istringstream line_stream(line);
			Size ires;
			char aa;
			std::string atom_name;
			Real shift;
			line_stream >> ires >> aa >> atom_name >> shift;

			if ( atom_name == "H" ) {
				atom_name = "HN";
			}

			//std::cout << "READ_SHIFTS " << ires << " " << aa << " " << atom_name << " " << shift << std::endl;

			boost::tuple<Size, char, std::string, Real> t(ires, aa, atom_name,
				shift);
			entries_.push_back(t);
		}
	}

	for ( Size i = 1; i <= entries_.size(); i++ ) {
		//std::cout << "ENTRIES: " << entries_[i].get<0>() << " " << entries_[i].get<1>() << " " << entries_[i].get<2>() << std::endl;
		resids_to_entries_map_.insert(std::make_pair(entries_[i].get<0> (), i));
	}
}

void CSTalosIO::write(std::ostream& out) {

	out << "DATA FIRST_RESID " << first_residue_index_ << std::endl;
	out << "DATA SEQUENCE " << sequence_ << std::endl;
	out << "VARS";
	for ( Size i = 1; i <= column_names_.size(); i++ ) {
		out << " " << column_names_[i];
	}
	out << std::endl << "FORMAT " << data_format_ << std::endl << std::endl;
	char buffer[50];
	char c1[2];
	c1[1] = 0;
	for ( Size i = 1; i <= entries_.size(); i++ ) {
		c1[0] = entries_[i].get<1> ();
		sprintf(buffer, data_format_.c_str(), entries_[i].get<0> (), c1,
			entries_[i].get<2> ().c_str(), entries_[i].get<3> ());
		out << buffer << std::endl;
	}
}

void CSTalosIO::get_tuples(Size residue_id, utility::vector1<boost::tuple<Size,
	char, std::string, Real> > results) {

	std::multimap<Size, Size>::iterator iter = resids_to_entries_map_.find(
		residue_id);

	//while (iter != resids_to_entries_map_.end()) {
	while ( iter->first == residue_id ) {
		results.push_back(entries_[iter->second]);
		++iter;
	}
}

void CSTalosIO::get_tuples(Size residue_id, utility::vector1<boost::tuple<Size,
	char, std::string, Real> > & results) const {

	std::multimap<Size, Size>::const_iterator iter =
		resids_to_entries_map_.find(residue_id);

	//while (iter != resids_to_entries_map_.end()) {
	while ( iter->first == residue_id ) {
		//std::cout << "GET_TUPLES" << residue_id << " " << iter->first << " " << iter->second << std::endl;
		results.push_back(entries_[iter->second]);
		++iter;
	}
}

bool CSTalosIO::has_atom(Size residue_id, std::string const & atom_name) const {

	std::multimap<Size, Size>::const_iterator iterat =
		resids_to_entries_map_.find(residue_id);
	if ( iterat == resids_to_entries_map_.end() ) {
		return false;
	}

	utility::vector1<boost::tuple<Size, char, std::string, Real> >
		used_for_searching_;

	get_tuples(residue_id, used_for_searching_);
	for ( Size i = 1; i <= used_for_searching_.size(); i++ ) {
		//std::cout << "FROMTUPLE " << residue_id << " " << used_for_searching_[i].get<0>() << " " << used_for_searching_[i].get<2>() << " " << used_for_searching_[i].get<1>() << std::endl;

		//Seriously, get_tuples(residue_id,...) also fetches tuples for OTHER residue IDs.
		//I have no idea why.
		if ( (used_for_searching_[i].get<0>() == residue_id) && (used_for_searching_[i].get<2> () == atom_name) ) {
			used_for_searching_.clear();
			return true;
		}
	}
	used_for_searching_.clear();

	return false;
}

Real CSTalosIO::get_shift(Size residue_id, std::string const & atom_name) const {

	std::multimap<Size, Size>::const_iterator iterat =
		resids_to_entries_map_.find(residue_id);
	if ( iterat == resids_to_entries_map_.end() ) {
		return false;
	}

	utility::vector1<boost::tuple<Size, char, std::string, Real> >
		used_for_searching_;

	get_tuples(residue_id, used_for_searching_);
	for ( Size i = 1; i <= used_for_searching_.size(); i++ ) {
		//std::cout << "GET_SHIFT " << residue_id << " " << atom_name << " " << used_for_searching_[i].get<2>() << std::endl;

		//Seriously, get_tuples(residue_id,...) also fetches tuples for OTHER residue IDs.
		//I have no idea why.
		if ( (used_for_searching_[i].get<0>() == residue_id) && (used_for_searching_[i].get<2> () == atom_name) ) {
			Real ret = used_for_searching_[i].get<3> ();
			used_for_searching_.clear();
			return ret;
		}
	}

	utility_exit_with_message(
		"[ERROR] Unable locate chemical shift for an atom " + atom_name
		+ " within a residue");
	return 0.0;
}

void CSTalosIO::set_up_atom_order() {

	order_of_atoms_.insert(std::pair<std::string, Size>("N", 1));
	order_of_atoms_.insert(std::pair<std::string, Size>("HA", 2));
	order_of_atoms_.insert(std::pair<std::string, Size>("HA2", 5));
	order_of_atoms_.insert(std::pair<std::string, Size>("HA3", 2));
	order_of_atoms_.insert(std::pair<std::string, Size>("C", 3));
	order_of_atoms_.insert(std::pair<std::string, Size>("CA", 4));
	order_of_atoms_.insert(std::pair<std::string, Size>("CB", 5));
	order_of_atoms_.insert(std::pair<std::string, Size>("HN", 6));
	order_of_atoms_.insert(std::pair<std::string, Size>("H", 6));
}

} // frag_picker
} // protocols
