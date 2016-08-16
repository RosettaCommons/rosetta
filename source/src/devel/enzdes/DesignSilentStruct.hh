// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/enzdes/DesignSilentStruct.hh
///
/// @brief protein silent-file structures for designs, also contains functionality to query
/// @brief a given pose and extract some more data to print
/// @author Florian Richter

#ifndef INCLUDED_devel_enzdes_DesignSilentStruct_hh
#define INCLUDED_devel_enzdes_DesignSilentStruct_hh


// mini headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>


#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>


#include <utility/pointer/ReferenceCount.hh>


// C++ Headers
#include <cstdlib>
#include <iostream>
//#include <utility/assert.hh>
//#include <vector>
#include <string>

#include <utility/vector1.hh>


//#include <map>
//#include <algorithm>

//eventually this will be moved to core/io/silent
namespace devel {
namespace enzdes {

class DesignSilentStruct : public core::io::silent::ProteinSilentStruct {

public:
	typedef core::io::silent::ProteinSilentStruct parent;

public:

	DesignSilentStruct(
		core::pose::Pose const & pose,
		std::string tag,
		bool const add_in,
		bool const onlyadd_in
	);

	DesignSilentStruct(
		core::pose::Pose const & pose,
		std::string tag, // = "empty_tag",
		utility::vector1<core::Size> const & spec_res_in,
		utility::vector1< std::string > const & rel_score_in,
		bool const add_in,
		bool const onlyadd_in );

	DesignSilentStruct(
		core::pose::Pose const & pose,
		std::string tag, // = "empty_tag",
		utility::vector1<core::Size> const & spec_res_in,
		utility::vector1< std::string > const & rel_score_in,
		std::map< core::Size, utility::vector1< std::pair< std::string, std::string > > > const & calculators,
		bool const add_in,
		bool const onlyadd_in );


	//this function will calculate the additional information

	void calculate_additional_info(
		core::pose::Pose const & pose,
		utility::vector1<core::Size> const & special_res,
		utility::vector1< std::string > const & score_terms,
		std::map< core::Size, utility::vector1< std::pair< std::string, std::string > > > const & calculators );

	void add_to_additional_silent_energies(
		utility::vector1< core::io::silent::SilentEnergy > const & silent_Es);

	//this function will go over the special_residues_ member object, and for each of these residues
	//write out the values for the relevant score terms as well as the evaluate the metric calculators
	void print_additional_info(std::ostream& out) const;

	void print_header( std::ostream& out ) const;

	void print_scores( std::ostream& out ) const;

	core::Real sum_constraint_terms( core::pose::Pose const & pose, int which_res);

private:

	utility::vector1< core::io::silent::SilentEnergy > additional_info_silent_energy_;

	bool print_additional_;
	bool print_only_additional_;

}; //class DesignSilentStruct

} //namespace enzdes
} //namespace devel


#endif
