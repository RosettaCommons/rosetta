// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/io/tarparser.cc
/// @brief  utitlities for reading tar files
/// @author Guangfeng Zhou, Frank DiMaio

#include<utility/io/tarparser.hh>
#include <cstring>
#include <climits>
#include <cassert>
#include <vector>
#include <stdexcept>

namespace utility {
namespace io {

// Tarfile reader
// very simple tar file reader.

TarParser::TarParser() {}

void TarParser::read(std::istream &datastream, std::string &filename, std::string &currentData)
{
	char nullBlock[TAR_BLOCK_SIZE];
	std::memset(nullBlock, 0, TAR_BLOCK_SIZE);

	TarFileHeader header;
	datastream.read((char*) &header, TAR_BLOCK_SIZE);

	if ( std::memcmp(&header, nullBlock, TAR_BLOCK_SIZE)==0 ) {
		throw std::runtime_error("End of file reached");
	}
	bool isValid = validateHeader(header);
	if ( !isValid ) {
		throw std::runtime_error("Invalid tar header");
	}
	if ( header.typeflag[0] != '0' && header.typeflag[0] != '\0' ) {
		throw std::runtime_error("File is not regular");
	}
	if ( datastream ) {
		filename = header.name;
		uint64_t filesize = parseFileSize(header);
		std::vector<char> data(filesize+1);
		datastream.read(data.data(), filesize);
		data[filesize] = '\0';
		currentData = std::string(data.begin(), data.end());
		int padding = (TAR_BLOCK_SIZE - (filesize % TAR_BLOCK_SIZE)) % TAR_BLOCK_SIZE;
		datastream.ignore(padding);
	} else {
		throw std::runtime_error("Error reading file");
	}
}

bool TarParser::validateHeader(TarFileHeader const &tarHeader)
{
	if ( std::memcmp(tarHeader.magic, TAR_MAGIC.c_str(), 6) == 0 ) {
		return true;
	}
	return false;
}

uint64_t TarParser::parseFileSize (TarFileHeader const &tarHeader)
{
	uint64_t filesize = 0;
	// when the file size is > 8 GB, the first byte (index 0) is used as a flag,
	// so the index starts from 1
	// the large file is indicated by the 7th bit of the first byte set to 1
	int start = (tarHeader.size[0] & (0x01<<7)) ? 1 : 0;
	int base = start ? 256 : 8;
	for ( int i=start; i<12; i++ ) {
		if ( !start && (tarHeader.size[i]==0 || tarHeader.size[i]==' ') ) continue;
		filesize *=base;
		filesize += start ? tarHeader.size[i] : (tarHeader.size[i]-'0');
	}
	return filesize;
}

} // namespace io
} // namespace utility
