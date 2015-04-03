// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/packing/HolesParams.cc
/// @brief  Packing Score Params
/// @author Will Sheffler

#ifndef INCLUDED_core_scoring_packing_HolesParams_hh
#define INCLUDED_core_scoring_packing_HolesParams_hh

#include <core/types.hh>


#include <utility/io/izstream.hh>
#include <utility/exit.hh>

// option key includes


#include <utility/vector1.hh>
#include <map>


namespace core {
namespace scoring {
namespace packing {


class HolesParams {

	typedef std::pair< core::Size, char > Key;

public:

	HolesParams()
	 : max_atom_types_(21), sep_ss_(21,true) {
	}

	HolesParams( std::string fname )
	 : max_atom_types_(21), sep_ss_(21,true) {
		read_data_file( fname );
	}

	void read_data_file( std::string fname ) {
		using namespace std;
		using namespace core;
		using namespace utility;

		utility::io::izstream in;
		in.close();
		in.clear();
		in.open( fname.c_str() );

      if( !in ) { utility_exit_with_message("HolesParams::read_data_file can't find file"); }

		string tmp;
		char ss;
		Size at;
		for( Size i = 1; i <=24; ++i ) {
			in >> tmp;
			// std::cerr << "in: " << tmp << std::endl;
		}
		while( true ) {
			// std::cerr << "LOOP START" << std::endl;
			in >> at;
			if( 999 == at ) break;
			in >> ss;
			// std::cerr << "atype / ss " << "'" << at << "' '" << ss << "'" << std::endl;
			vector1<char> SS;
			if( '*'==ss ) { SS.push_back('E'); SS.push_back('H'); SS.push_back('L'); sep_ss_[at] = false; }
			else SS.push_back(ss);

			Key k(at,SS[1]);
			in >> nb_weights_[k];
			sa_weights_[k].resize(20,0.0);
			for( Size i = 1; i <= 20; ++i ) in >> sa_weights_[k][i];
			in >> intercepts_[k];

			// std::cerr << "nb:   " << nb_weights_[k] << std::endl;
			// for( Size i = 1; i <= 20; ++i ) std::cerr << sa_weights_[k][i] << " ";
			// std::cerr << std::endl;
			// std::cerr << "icpt: " << intercepts_[k] << std::endl;

			for( Size i = 2; i <= SS.size(); ++i ) {
				nb_weights_[Key(at,SS[i])] = nb_weights_[k];
				sa_weights_[Key(at,SS[i])] = sa_weights_[k];
				intercepts_[Key(at,SS[i])] = intercepts_[k];
			}

		}
		in >> intercept_;
	}

	core::Real nb_weight( core::Size const at, char const ss                ) const { return nb_weights_.find(Key(at,ss))->second; }
	core::Real sa_weight( core::Size const at, char const ss, core::Size sa ) const { return sa_weights_.find(Key(at,ss))->second[sa]; }
	core::Real intercept( core::Size const at, char const ss                ) const { return intercepts_.find(Key(at,ss))->second; }
	core::Real intercept(                                                   ) const { return intercept_; }

	core::Size max_atom_types() const { return max_atom_types_; }
	bool sep_ss( core::Size at ) const { return sep_ss_[at]; }

	bool have_params( core::Size const at, char const ss ) const {
		return sa_weights_.find(Key(at,ss)) != sa_weights_.end();
	}

	friend
	std::ostream &
	operator<< (std::ostream & out, HolesParams const & hp) {
		using namespace core;
		using namespace std;
		for( Size at = 1; at <= hp.max_atom_types(); ++at ) {
			if( hp.sep_ss(at) ) {
				out<<at<<" E "<<hp.nb_weight(at,'E'); for(Size i=1;i<21;++i) out<<" "<<hp.sa_weight(at,'E',i); out<<" "<<hp.intercept(at,'E')<<endl;
				out<<at<<" H "<<hp.nb_weight(at,'H'); for(Size i=1;i<21;++i) out<<" "<<hp.sa_weight(at,'H',i); out<<" "<<hp.intercept(at,'H')<<endl;
				out<<at<<" L "<<hp.nb_weight(at,'L'); for(Size i=1;i<21;++i) out<<" "<<hp.sa_weight(at,'L',i); out<<" "<<hp.intercept(at,'L')<<endl;
			} else {
				out<<at<<" * "<<hp.nb_weight(at,'E'); for(Size i=1;i<21;++i) out<<" "<<hp.sa_weight(at,'E',i); out<<" "<<hp.intercept(at,'E')<<endl;
			}
		}
		out << "999 " << hp.intercept() << endl;
		return out;
	}

private:

	core::Size const max_atom_types_;
	std::map<Key,utility::vector1<core::Real> > sa_weights_;
	std::map<Key,core::Real> nb_weights_;
	std::map<Key,core::Real> intercepts_;
	utility::vector1<bool> sep_ss_;
	core::Real intercept_;


};

} // end packing
} // end scoring
} // end core

#endif // INCLUDED_core_scoring_packing_HolesParams_HH
