// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/delucasl/params_tester.cc
/// @author Sam DeLuca

#include <devel/init.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <iostream>
#include <utility/excn/Exceptions.hh>

int main(int argc, char* argv[])
{
	try {
		devel::init(argc,argv);
		core::chemical::ResidueTypeSetCAP residue_types =
			core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");

	}catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
