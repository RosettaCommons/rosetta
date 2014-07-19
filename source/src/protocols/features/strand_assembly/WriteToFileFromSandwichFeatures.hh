// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SandwichFragment.hh
///
/// @brief Small helper class that stores the start and end of a strand secondary structure

/// @author Doo Nam Kim (started from Tim jacobs' code)

#ifndef INCLUDED_protocols_features_strand_assembly_WriteToFileFromSandwichFeatures_HH
#define INCLUDED_protocols_features_strand_assembly_WriteToFileFromSandwichFeatures_HH

//Devel
#include <protocols/features/strand_assembly/StrandAssemblyCommon.hh>
//#include <protocols/features/strand_assembly/SandwichFeatures.hh>

using namespace std;

namespace protocols {
namespace features {
namespace strand_assembly {

	utility::vector1<Size>
	get_vector_loop_AA_distribution(
		StructureID struct_id,
		utility::sql_database::sessionOP	db_session,
		string loop_kind);

	core::Size
	write_AA_distribution_without_direction_to_a_file(
		string	tag,
		StructureID	struct_id,
		utility::sql_database::sessionOP	db_session);

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif /* WriteToFileFromSandwichFeatures_HH_ */
