// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/tarparser.hh
/// @brief  utitlities for reading tar files
/// @author Guangfeng Zhou, Frank DiMaio


#ifndef INCLUDED_utility_io_tarparser_hh
#define INCLUDED_utility_io_tarparser_hh

#include <iostream>
#include <fstream>
#include <string>

#include <cstdint>

namespace utility {
namespace io {

constexpr int TAR_BLOCK_SIZE = 512;

/* tarfile reader
* derived from https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/util/compress/api/tar.cpp
*  ===========================================================================
*                            PUBLIC DOMAIN NOTICE
*               National Center for Biotechnology Information
*  This software/database is a "United States Government Work" under the
*  terms of the United States Copyright Act.  It was written as part of
*  the author's official duties as a United States Government employee and
*  thus cannot be copyrighted.  This software/database is freely available
*  to the public for use. The National Library of Medicine and the U.S.
*  Government have not placed any restriction on its use or reproduction.
*  ===========================================================================
*  Authors:  Vladimir Ivanov, Anton Lavrentiev
*/
struct TarFileHeader {        // byte offset
	char name[100];           //   0
	char mode[8];             // 100
	char uid[8];              // 108
	char gid[8];              // 116
	char size[12];            // 124
	char mtime[12];           // 136
	char checksum[8];         // 148
	char typeflag[1];         // 156
	char linkname[100];       // 157
	char magic[6];            // 257
	char version[2];          // 263
	char uname[32];           // 265
	char gname[32];           // 297
	char devmajor[8];         // 329
	char devminor[8];         // 337
	char prefix[155];         // 345
	char padding[12];         // 500
	// 512
};

class TarParser {
public:
	TarParser();

	void read(std::istream &datastream, std::string &filename,
		std::string &currentData);

	static uint64_t parseFileSize (TarFileHeader const &tarHeader);

private:
	bool validateHeader(TarFileHeader const &tarHeader);
	std::string const TAR_MAGIC = "ustar ";

};

}  // namespace io
}  // namespace utility

#endif  // INCLUDED_utility_io_tarparser_hh
