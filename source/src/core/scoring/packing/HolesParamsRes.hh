// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packing/HolesParamsRes.cc
/// @brief  Packing Score Params
/// @author Will Sheffler

#ifndef INCLUDED_core_scoring_packing_HolesParamsRes_hh
#define INCLUDED_core_scoring_packing_HolesParamsRes_hh

#include <core/types.hh>

#include <basic/database/open.hh>
#include <basic/options/option.hh>

#include <utility/io/izstream.hh>
#include <utility/exit.hh>

// option key includes

#include <basic/options/keys/holes.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace packing {


class HolesParamsRes {
public:

	HolesParamsRes( std::string fname = std::string("") ) {
		read_data_file( fname );
	}

	void read_data_file( std::string fname ) {
		using namespace std;
		using namespace utility;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace utility;

		res_size_["ALA"] =  5;
		res_size_["ARG"] = 11;
		res_size_["ASN"] =  8;
		res_size_["ASP"] =  8;
		res_size_["CYS"] =  6;
		res_size_["GLN"] =  9;
		res_size_["GLU"] =  9;
		res_size_["GLY"] =  4;
		res_size_["HIS"] = 10;
		res_size_["ILE"] =  8;
		res_size_["LEU"] =  8;
		res_size_["LYS"] =  9;
		res_size_["MET"] =  8;
		res_size_["PHE"] = 11;
		res_size_["PRO"] =  8;
		res_size_["SER"] =  6;
		res_size_["THR"] =  7;
		res_size_["TRP"] = 14;
		res_size_["TYR"] = 12;
		res_size_["VAL"] =  7;

		utility::io::izstream in;
		if ( "" != fname ) {
			in.close();
			in.clear();
			in.open( fname.c_str() );
		} else {
			std::string paramfile = (std::string)( basic::options::option[ basic::options::OptionKeys::holes::params ]() );
			if ( paramfile[0] == '/' || paramfile[0] == '.' || paramfile[0] == '~' ) {
				in.close();
				in.clear();
				in.open( paramfile.c_str() );
			} else {
				basic::database::open( in, paramfile );
			}
		}
		string s;
		core::Real r;
		in >> s;
		while ( s != "END" ) {
			if ( s == "INTERCEPT" ) {
				in >> r;
				rho_ = r;
			} else {
				if ( res_size_.find(s) == res_size_.end() ) {
					std::cerr << "HOLES: can't find res type '" << s << "'" << endl;
					utility_exit_with_message( "HOLES: can't find res type" );
				}
				residues_.push_back(s);
				int size = res_size_[s];
				params_[s].reserve(24*size+1);
				for ( int i = 1; i <= 24*size+1; i++ ) {
					in >> r;
					params_[s].push_back(r);
				}
			}
			in >> s;
		}
	}

	utility::vector1<core::Real> const & param(std::string res) const {
		return params_.find(res)->second;
	}

	core::Real rho() const { return rho_; }

	core::Real rho( std::string res ) const {
		if ( have_params(res) ) {
			return params_.find(res)->second[params_.find(res)->second.size()];
		} else return 0;
	}

	std::string residues( core::Size i ) const {
		return residues_[i];
	}

	bool have_params( std::string res ) const {
		return params_.find(res) != params_.end();
	}

private:
	utility::vector1<std::string> residues_;
	std::map<std::string,core::Size> res_size_;
	std::map<std::string,utility::vector1<core::Real> > params_;
	core::Real rho_;
};

} // end packing
} // end scoring
} // end core

#endif // INCLUDED_core_scoring_packing_HolesParamsRes_HH
