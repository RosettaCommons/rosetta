// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

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
	return new DeleteAllRotamerSetOperation( *this );
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
