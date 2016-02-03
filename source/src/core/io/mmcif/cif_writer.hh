// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file core/io/mmcif/util.hh
/// @brief Functions for MMCIF writing.
/// @author Andy Watkins (andy.watkins2@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_core_io_mmcif_cif_writer_hh
#define INCLUDED_core_io_mmcif_cif_writer_hh

#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>

#include <core/pose/Pose.hh>

namespace core {
namespace io {
namespace mmcif {

void
dump_cif( core::pose::Pose const & pose, std::string const & cif_file);


void
dump_cif( std::string const & cif_file, StructFileRepOP sfr, StructFileReaderOptions const & options );



} //core
} //io
} //mmcif


#endif //core/io/mmcif_util_hh

