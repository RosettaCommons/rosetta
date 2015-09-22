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
/// @author Nikolas Sgourakis sgourn@u.w.edu

// unit headers
#include <protocols/frag_picker/JCouplingIO.hh>

// project headers
#include <core/types.hh>

#include <basic/Tracer.hh>

// utility headers
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

#include <string>
#include <map>
#include <utility/exit.hh>

// boost headers
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/option.hh>

namespace protocols {
namespace frag_picker {

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;

static THREAD_LOCAL basic::Tracer trJCouplingIO( "protocols.frag_picker.TalosReader" );


void JCouplingIO::read(std::string const & file_name) {

	utility::io::izstream data(file_name.c_str());
	trJCouplingIO.Info << "read Jcoupling data from " << file_name << std::endl;
	if ( !data ) {
		utility_exit_with_message("[ERROR] Unable to open Jcoupling file: "
			+ file_name);
	}

	std::string line;
	Size len;
	Real a, b, c, theta;

	getline(data, line);
	std::istringstream line_stream(line);
	line_stream >> len >> a >> b >> c >> theta;

	utility::vector1< Real > tmpv(len);
	utility::vector1< utility::vector1< Real > > temp_data( 1, tmpv );

	while ( !data.eof() ) {
		getline(data, line);
		std::istringstream line_stream(line);
		Real val, dev;
		Size res;
		line_stream >> res >> val >> dev;

		if ( !(option[frags::filter_JC].user() && (val >4) && (val<6)  ) ) {
			temp_data[res].push_back(val);
			temp_data[res].push_back(dev);
			//   std::cout << "DATA " << res << " " << val << " " << dev << std::endl;
		}
		//std::cout << "DATA " << res << " " << val << " " << dev << std::endl;

	}

	//Space for error checking and asserts

	data_ = temp_data;
	A_ = a;
	B_ = b;
	C_ = c;
	THETA_ = theta;
	sequence_length_ = len;
}

std::pair< Real, Real > JCouplingIO::get_data( Size const res_num, bool & has_data ){

	has_data = false;

	std::pair< Real, Real > temp;

	//std::cout << "GETTING " << res_num << " " << data_[res_num].size() << std::endl;

	if ( data_[res_num].size() == 2 ) {
		temp.first = data_[res_num][1];
		temp.second = data_[res_num][2];
		has_data = true;
	}
	return temp;
}

utility::vector1< Real > JCouplingIO::get_parameters() {

	utility::vector1<Real> params;

	params.push_back(A_);
	params.push_back(B_);
	params.push_back(C_);
	params.push_back(THETA_);

	return params;
}


} // frag_picker
} // protocols
