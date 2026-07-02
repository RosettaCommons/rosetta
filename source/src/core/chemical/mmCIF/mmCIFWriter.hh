// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////////////////////////////
///
/// @brief
/// Class to convert a ResidueType to mmCIF.
///
/// @file   src/core/chemical/sdf/mmCIFWriter.hh
///
/// @details
///
/// @author Rocco Moretti (rmorettiase@gmail.com)
///
///
/////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_core_chemical_mmCIF_mmCIFWriter_hh
#define INCLUDED_core_chemical_mmCIF_mmCIFWriter_hh

#include <core/chemical/mmCIF/mmCIFWriter.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

#include <utility/vector1.hh>
#include <utility/VirtualBase.hh>
#include <map>
#include <core/types.hh>

#include <utility/gemmi_util.fwd.hh>

namespace core {
namespace chemical {
namespace mmCIF {

class mmCIFWriter : public utility::VirtualBase
{

public:

	mmCIFWriter() = default;
	~mmCIFWriter() override = default;

	void
	write_file( std::string const & filename, core::chemical::ResidueType const & restype );

	void
	write_stream( std::ostream & output_stream, core::chemical::ResidueType const & restype );

	void
	add_data_to_block( gemmi::cif::Block & block, core::chemical::ResidueType const & restype );

private:

};



}
}
}
#endif
