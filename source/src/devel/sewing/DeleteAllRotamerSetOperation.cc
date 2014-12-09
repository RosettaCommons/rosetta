// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /rosetta/rosetta_source/src/devel/sewing/DeleteAllRotamersSetOperation.cc
/// @brief 
/// @author Tim Jacobs

//Unit
#include<devel/sewing/DeleteAllRotamerSetOperation.hh>

namespace devel {
namespace sewing {

core::pack::rotamer_set::RotamerSetOperationOP
DeleteAllRotamerSetOperation::clone() const
{
	return core::pack::rotamer_set::RotamerSetOperationOP( new DeleteAllRotamerSetOperation( *this ) );
}

void DeleteAllRotamerSetOperation::alter_rotamer_set(
			core::pose::Pose const &,
			core::scoring::ScoreFunction const &,
			core::pack::task::PackerTask const &,
			core::graph::GraphCOP,
			core::pack::rotamer_set::RotamerSet & rotamer_set)
{
	utility::vector1<bool> rotamer_vector(rotamer_set.num_rotamers(), true);
	rotamer_set.drop_rotamers(rotamer_vector);
}

} //sewing namespace
} //devel namespace
