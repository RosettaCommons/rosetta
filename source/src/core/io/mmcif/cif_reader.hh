// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/mmcif/StructFileReader.hh
/// @author Andy Watkins (andy.watkins2@gmail.com)


#ifndef INCLUDED_core_io_mmcif_cif_reader_hh
#define INCLUDED_core_io_mmcif_cif_reader_hh


// Unit headers
#include <core/io/StructFileReaderOptions.fwd.hh>
#include <core/io/StructFileRep.fwd.hh>

// Project headers

// Utility headers
#include <utility/gemmi_util.fwd.hh>

// C++ headers
#ifdef WIN32
#include <string>
#endif



namespace core {
namespace io {
namespace mmcif {

typedef std::string String;


/// @details May raise a std::runtime_error from Gemmi parsing
StructFileRepOP create_sfr_from_cif_file( gemmi::cif::Document & cifdoc, StructFileReaderOptions const & options );


} // namespace mmcif
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_mmcif_mmCIFReader_hh
