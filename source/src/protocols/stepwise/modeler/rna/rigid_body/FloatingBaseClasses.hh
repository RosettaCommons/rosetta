// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/rigid_body/FloatingBaseClasses.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_FloatingBaseClasses_HH
#define INCLUDED_protocols_stepwise_rna_FloatingBaseClasses_HH

#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <utility>

#include <map>

typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace rigid_body {


struct BaseBin{

	int centroid_x;
	int centroid_y;
	int centroid_z;
	int euler_alpha;
	int euler_z;
	int euler_gamma;

};


struct
compare_base_bin{

	//The expression comp(a,b), where comp is an object of this comparison class and a and b are key values, shall return true if a is to be placed at an earlier position than b in a strict weak ordering operation

  bool
	operator() ( BaseBin const & first, BaseBin const & second ) const {

		if ( first.centroid_x != second.centroid_x ) return ( first.centroid_x < second.centroid_x ); //x
		if ( first.centroid_y != second.centroid_y ) return ( first.centroid_y < second.centroid_y ); //y
		if ( first.centroid_z != second.centroid_z ) return ( first.centroid_z < second.centroid_z ) ; //z
		if ( first.euler_alpha != second.euler_alpha ) return ( first.euler_alpha < second.euler_alpha );
		if ( first.euler_gamma != second.euler_gamma ) return ( first.euler_gamma < second.euler_gamma );
		if ( first.euler_z != second.euler_z ) return ( first.euler_z < second.euler_z );

		return false; //Equality case.
	}

};

	typedef	std::map< BaseBin, int, compare_base_bin > BaseBinMap;

struct
compare_int_pair{

	//The expression comp(a,b), where comp is an object of this comparison class and a and b are key values, shall return true if a is to be placed at an earlier position than b in a strict weak ordering operation


  bool
	operator() ( std::pair < int, int > const & pair_one, std::pair < int, int > const & pair_two ) const {

		if ( pair_one.first != pair_two.first ) return ( pair_one.first < pair_two.first );
		if ( pair_one.second != pair_two.second ) return ( pair_one.second < pair_two.second );

		return false; //Equality case.
	}

};


} //rigid_body
} //rna
} //modeler
} //stepwise
} //protocols

#endif
