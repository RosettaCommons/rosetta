// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/file_data.hh
///
/// @brief
/// @author Sergey Lyskov

#ifndef INCLUDED_core_pose_Remarks_hh
#define INCLUDED_core_pose_Remarks_hh

#include <core/pose/Remarks.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <iostream>
#include <string>
#include <vector>

namespace core {
namespace pose {

typedef std::string String;

class RemarkInfo
{
public:
/// @brief default constructor to initialize all values
RemarkInfo() : num( 0 ), value() {}

/// For now, all member names have the same names as fields in PDB standard.
int num;
String value;

/// @brief Debug printing, serialazing to Tracer like object.
friend std::ostream& operator <<(std::ostream &os, RemarkInfo const & ri) {
	os << "<RemarkInfo>{" << "num=" << ri.num << " value=" << ri.value << "}";
		return os;
}

};

class Remarks : public	std::vector< RemarkInfo >, public utility::pointer::ReferenceCount
{
};

} // namespace pose
} // namespace core

#endif
