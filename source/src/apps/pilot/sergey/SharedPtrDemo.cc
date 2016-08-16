// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
///
/// @brief  Demo for boost shared_ptr
/// @author Sergey Lyskov


#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/excn/Exceptions.hh>


#include <boost/shared_ptr.hpp>

#include <memory>

typedef boost::shared_ptr<int> intSP__;


typedef boost::shared_ptr<int> intSP;
typedef boost::shared_ptr<int const> intCSP;


class Base : public utility::pointer::ReferenceCount { };
typedef utility::pointer::owning_ptr< Base > BaseOP;
typedef utility::pointer::owning_ptr< Base const > BaseCOP;
typedef boost::shared_ptr<Base> BaseSP;
typedef boost::shared_ptr<Base const> BaseCSP;

class Son : public Base { };
typedef utility::pointer::owning_ptr< Son > SonOP;
typedef utility::pointer::owning_ptr< Son const > SonCOP;
typedef boost::shared_ptr<Son> SonSP;
typedef boost::shared_ptr<Son const> SonCSP;


void foo_baseOP(BaseOP) {}
void foo_baseCOP(BaseCOP) {}

void foo_baseSP(BaseSP) {}
void foo_baseCSP(BaseCSP) {}

void foo_const_intSP(intCSP) {}

int main( int argc, char * argv [] )
{

	try {

	intSP a;
	foo_const_intSP(a);

	BaseCOP base_cop;
	//foo_baseOP(base_cop);  <-- Does not work...
	foo_baseCOP(base_cop);

	SonCOP son_cop;
	//foo_baseOP(son_cop);  <-- Does not work...
	foo_baseCOP(son_cop);

	BaseCSP base_csp;
	//foo_baseSP(base_csp);  <-- Does not work...
	foo_baseCSP(base_csp);

	SonCSP son_csp;
	SonSP son_sp;
	//foo_baseSP(son_csp);  <-- Does not work...
	foo_baseSP( boost::dynamic_pointer_cast<Base>(son_sp) );
	foo_baseCSP(son_csp);

	BaseSP baseSP = boost::dynamic_pointer_cast<Base>(son_sp);
	SonSP  _2 = boost::dynamic_pointer_cast<Son>(baseSP);

    /*
	boost::weak_ptr<Son const> SonCWP(son_csp);
	SonCSP son_csp1 = SonCWP.lock();
	if(son_csp1) {
		// Object was not released!
	} */

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
