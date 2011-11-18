// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Util.hh
/// @brief
/// @detailed
///
///  @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_swa_SWA_Dinucleotide_Sampler_Util_HH
#define INCLUDED_protocols_swa_SWA_Dinucleotide_Sampler_Util_HH


#include <protocols/swa/rna/StepWiseRNA_Dinucleotide_Sampler_Util.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <string>
#include <map>
#include <core/chemical/AA.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <numeric/angle.functions.hh> // Need this to prevent the compiling error: 'principal_angle_degrees' is not a member of 'numeric' Oct 14, 2009
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.hh>
#include <set>

typedef  numeric::xyzMatrix< core::Real > Matrix;

namespace protocols {
namespace swa {
namespace rna {


struct Anchor_ribose_stub{
	numeric::xyzVector<core::Real> origin;
  Matrix coordinate_matrix;
	Matrix invert_coordinate_matrix;
};

//Should make sure that alpha and gamma lies in the [-Pi:Pi] range.
struct Euler_angles{ 
	core::Real z; 
	core::Real alpha; //phi
	core::Real beta; //theta
	core::Real gamma; //psi
};

struct Base_bin{

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
	operator() (Base_bin const & first , Base_bin const & second) const {
		
		if(first.centroid_x != second.centroid_x) return (first.centroid_x<second.centroid_x); //x 
		if(first.centroid_y != second.centroid_y) return (first.centroid_y<second.centroid_y); //y
		if(first.centroid_z != second.centroid_z) return (first.centroid_z<second.centroid_z) ; //z
		if(first.euler_alpha != second.euler_alpha) return (first.euler_alpha < second.euler_alpha);
		if(first.euler_gamma != second.euler_gamma) return (first.euler_gamma < second.euler_gamma);
		if(first.euler_z != second.euler_z) return (first.euler_z < second.euler_z);

		return false; //Equality case.
	}

};

struct 
compare_int_pair{

	//The expression comp(a,b), where comp is an object of this comparison class and a and b are key values, shall return true if a is to be placed at an earlier position than b in a strict weak ordering operation


  bool 
	operator() (std::pair<int, int> const & pair_one , std::pair<int, int> const & pair_two) const {
		
		if(pair_one.first != pair_two.first) return (pair_one.first<pair_two.first); 
		if(pair_one.second != pair_two.second) return (pair_one.second<pair_two.second); 

		return false; //Equality case.
	}

};

struct 
compare_test{
  bool 
	operator() (core::Real const & first, core::Real const & second) const {
		return first< second;
	}
};
////////////////////////////////////////////////////////////////////////////

Base_bin
Get_base_bin(numeric::xyzVector<core::Real> const & centriod, Euler_angles const &  euler_angles);

Anchor_ribose_stub
Get_anchor_ribose_stub(core::conformation::Residue const & rsd, bool verbose=true);

numeric::xyzVector<core::Real> //This is equivalent to Vector
Convert_base_centroid_to_anchor_ribose_coord_system( numeric::xyzVector<core::Real> const & centroid, Anchor_ribose_stub const & anchor_ribose_stub);

numeric::xyzMatrix< core::Real >
Convert_base_coordinate_matrix_to_anchor_ribose_coord_system(numeric::xyzMatrix< core::Real > const & base_coordinate_matrix, Anchor_ribose_stub const & anchor_ribose_stub);
	
Euler_angles
Get_euler_angles( numeric::xyzMatrix< core::Real > const & base_coordinate_matrix_in_anchor_ribose_coord_system);

int
DOF_bin_value(std::map<Base_bin , int , compare_base_bin>::const_iterator const & base_bin_it, std::string const & DOF);

core::Real
DOF_bin_size(std::string const & DOF);

void
Analyze_base_bin_map(std::map<Base_bin , int , compare_base_bin> const & base_bin_map, std::string const foldername);

void
Analyze_base_bin_map(std::map<Base_bin , int , compare_base_bin> const & base_bin_map, std::string const & DOF_one, std::string const & DOF_two);

void
Analyze_base_bin_map_old(std::map<Base_bin , int , compare_base_bin> const & base_bin_map, bool const Is_dinucleotide);

void
Sample_base_test();

}
}
}

#endif


