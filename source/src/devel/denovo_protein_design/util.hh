// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @details
/// @author Grant Murphy


#ifndef INCLUDE_src_devel_denovo_protein_design_util_hh
#define INCLUDE_src_devel_denovo_protein_design_util_hh

// Package Headers

// Project Headers

// Utility Headers
#include <core/types.hh>


// C++ Headers

#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>

//#include <sstream>
//#include <utility/vector1.hh>

namespace devel {
namespace denovo_protein_design {

void design_setup( core::pose::Pose & pose, core::pack::task::TaskFactoryOP designtaskfactory );
//
core::Size numberhelices( core::pose::Pose & pose);
//
void create_nucleated_sequence_from_template_pdb( std::string & nucleated_sequence , utility::vector1< bool > & KeepNativeTorsions );
// reads in a secondary structure file from the command line and places info in a vector
void read_ss_file_into_vector( std::string const & filename,  utility::vector1< char > & secstructs  );
//
void read_hp_file_into_vector( std::string const & filename,  utility::vector1< char > & hydrophobic_polar_sequence_v );
// creates a hydrophobic polar pattern from a secondary structure
void create_hydrophobic_polar_pattern_from_secstruct_vector( utility::vector1< char > & hydrophobic_polar_sequence_v , utility::vector1< char > const & secstructs );
// creates an amino acid sequence from a hydrophic polar pattern
void create_sequence_from_hydrophobic_polar_pattern( std::string & sequence, utility::vector1< char > const & hydrophobic_polar_sequence_v);
// returns at random  B P
char random_hydrophobic_polar();
// returns a random polar aa
char random_polar_aa();
// returns a random apolar aa
char random_apolar_aa();
// returns a random heptad helix position
char random_helix_position();

}
}

#endif
