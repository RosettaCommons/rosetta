// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <utility/io/izstream.hh>

// C++ headers
#include <string>
#include <map>

// boost headers
#include <boost/tuple/tuple.hpp>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {

static basic::Tracer
tr("protocols.frag_picker.PhiPsiTalosIO");

void PhiPsiTalosIO::read(std::string const & file_name) {

	utility::io::izstream data(file_name.c_str());
	tr.Info << "read talos data from " << file_name << std::endl;
	if ( !data ) {
		utility_exit_with_message("[ERROR] Unable to open talos file: "
			+ file_name);
	}

	std::string line;
	last_residue_index_ = 0;
	first_residue_index_ = 1;
	bool first_not_found = true;
	utility::vector1< std::string > vars;
	while ( getline(data, line) ) {
		std::istringstream line_stream(line);
		utility::vector1<std::string> strs;
		std::string token;
		while ( line_stream >> token ) {
			strs.push_back(token);
		}
		if ( strs.size()==0 ) continue;
		tr.Trace << "token: " << strs.size() << " in line " << line << std::endl;
		if ( strs[1] == "DATA" ) {
			if ( strs[2] == "SEQUENCE" ) {
				for ( core::Size i = 3; i <= strs.size(); ++i ) {
					sequence_ += strs[i];
				}
			} else if ( strs[2] == "FIRST_RESID" ) {
				first_not_found = false;
				line_stream >> first_residue_index_;
				tr.Info << "FIRST_RESID entry in TALOS file. Setting first-residue to " << first_residue_index_ << std::endl;
			} else {
				tr.Warning << "Unrecognized DATA entry:" << line
					<< std::endl;
			}
		}
		if ( strs[1] == "VARS" ) {
			for ( core::Size i = 2; i <= strs.size(); ++i ) {
				vars.push_back(strs[i]);
			}
			if ( ( vars.size()!=10 && vars.size()!=11 ) || vars[1]!="RESID" || vars[2]!="RESNAME" || vars[9]!="COUNT" || vars.back()!="CLASS" ) {
				tr.Warning << "incompatible format in TALOS+ file "+file_name
					+".\n Expected VARS  RESID RESNAME PHI PSI DPHI DPSI DIST S2 COUNT CS_COUNT CLASS\n "
					+" or      VARS  RESID RESNAME PHI PSI DPHI DPSI DIST S2 COUNT CLASS\n "
					+" found instead: ";
				tr.Warning << vars.size() << " VARS: ";
				for ( core::Size i=1; i<=vars.size(); ++i ) {
					tr.Warning << vars[i] << " ";
				}
				tr.Warning << " LAST VAR: ->" << vars.back() << ":";
				tr.Warning << std::endl;
			}
		}
		if ( strs[1] == "FORMAT" ) {
			data_format_ = line.substr(7);
		}

		if ( (strs.size() == vars.size()) && (strs[1] != "REMARK") ) {
			char aa;
			std::istringstream line_stream(line);
			core::Real phi, psi, d_phi, d_psi, dist, s2;
			core::Size res_id, count;
			std::string cls;
			if ( vars.size()==10 ) {
				line_stream >> res_id >> aa >> phi >> psi >> d_phi >> d_psi >> dist
					>> s2 >> count >> cls;
			} else {
				core::Size cs_count;
				line_stream >> res_id >> aa >> phi >> psi >> d_phi >> d_psi >> dist
					>> s2 >> count >> cs_count >> cls;
			}
			data_format_ = " %4d %s %8.3f %8.3f %8.3f %8.3f %8.3f %5.3f %2d %s";
			//this adds sequence twice  sequence_ += aa;
			typedef boost::tuple<core::Size, char, core::Real,core::Real, core::Real, core::Real, core::Real, core::Real, core::Size, std::string> PhiPsiTalosLineEntry;
			PhiPsiTalosLineEntry  t(res_id, aa, phi, psi, d_phi, d_psi, dist,
				s2, count, cls);
			entries_.insert(std::pair<core::Size, PhiPsiTalosLineEntry > ( res_id, t) );
			if ( last_residue_index_ < res_id ) {
				last_residue_index_ = res_id;
			}
		}
	}

	if ( first_not_found ) {
		tr.Warning
			<< "FIRST_RESID keyword didn't show up in a file header\n\tAssuming the first residue id is 1"
			<< std::endl;
	}
	if ( sequence_.length() == 0 ) {
		for ( core::Size i = first_residue_index_; i <= last_residue_index_; ++i ) {
			if ( has_entry(i) ) {
				sequence_ += entries_.find(i)->second.get<1> ();
			} else {
				sequence_ += 'X';
			}
		}
		tr.Warning << "Could not find DATA SEQUENCE in file "<<file_name << std::endl
			<< "Sequence based on entires is:" << std::endl
			<< sequence_ << std::endl;
	}
}


void PhiPsiTalosIO::write(std::ostream& out) {

	out << "DATA FIRST_RESID " << first_residue_index_ << "\n";
	out << "DATA SEQUENCE " << sequence_ << "\n";
	out << "VARS";
	for ( core::Size i = 1; i <= column_names_.size(); i++ ) {
		out << " " << column_names_[i];
	}
	out << "\n" << "FORMAT " << data_format_ << std::endl << "\n";
	char buffer[100];
	char c1[2];
	c1[1] = 0;
	for ( core::Size i = 1; i <= entries_.size(); i++ ) {
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
