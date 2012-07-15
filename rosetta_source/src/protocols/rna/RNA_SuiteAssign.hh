// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rna/RNA_SuiteAssign.hh
/// @brief RNA suite assignment ported from suitename program
/// @author Fang-Chieh Chou


#ifndef INCLUDED_protocols_rna_RNA_SuiteAssign_HH
#define INCLUDED_protocols_rna_RNA_SuiteAssign_HH

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/io/ozstream.fwd.hh>

// Utility headers

// ObjexxFCL headers

//// C++ headers
#include <string>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>


using namespace core;

namespace protocols {
namespace rna {

class suite_info {
public:
	std::string name;
	Size classifier;
	utility::vector1 <Real> torsion;

	suite_info ( std::string const name_in, Size const classifier_in, utility::vector1 <Real>  const & torsion_in ) :
		name( name_in ),
		classifier( classifier_in ),
		torsion( torsion_in )
	{}

	suite_info () : name( "" ), classifier( 0 ) {}
};

class RNA_suite_list {
public:

	RNA_suite_list();
	~RNA_suite_list();

	suite_info name2suite( std::string const name );
	utility::vector1 <suite_info> full_list () const { return all_suites; };


private:
	void init_all_suite();
	utility::vector1 <Real> create_torsions( Real const delta1, Real const epsilon, Real const zeta, 
	Real const alpha, Real const beta, Real const gamma, Real const delta2);

	utility::vector1 <suite_info> all_suites;	
};

std::pair <std::string, std::pair <Size, Real> > suite_assign(pose::Pose const & pose, Size const res);
}
}

#endif
