// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/BinaryRNASilentStruct.cc
///
/// @brief
/// @author Frank DiMaio
/// @author Mike Tyka
/// @author Rhiju Das

// mini headers
#include <utility/Binary_Util.hh>

// C++ Headers
#include <string>
#include <algorithm> //For std::min

namespace utility {

void swap4_aligned(void *v, long ndata) {
	uint32_t *data = (uint32_t *) v; // Unsigned as there's issues with bitshifting signed types
	for ( long i=0; i<ndata; i++ ) {
		uint32_t *N;
		N = data + i;
		*N=(((*N>>24)&0xff) | ((*N&0xff)<<24) | ((*N>>8)&0xff00) | ((*N&0xff00)<<8));
	}
}

/// @brief Given a block of memory (with memory pointing to the first byte) and a byte
/// count (length), convert every three bytes of memory into four bytes of ASCII characters,
/// and populate the string (jar) with those characters.
/// @details Note that decode6bit requires the length of the memory block into which
/// the contents of jar will be extracted, so the length of the memory location here
/// should be stored somehow if one hopes to accurately reproduce the bytes stored.
void encode6bit(const unsigned char* memory, unsigned int length, std::string &jar){
	jar = "";
	unsigned int i;
	unsigned int fourcount=0;
	unsigned int linewidth=15;
	for ( i=0; i<length; ) {
		unsigned char buffer[3] = {0,0,0};
		unsigned char outbuffer[4] = {0,0,0,0};
		int ibuf=0;
		for ( ; ((i<length)&&(ibuf<3)); i++ ) {
			buffer[ibuf] = memory[i];
			ibuf++;
		}
		encode_24_to_32(buffer[0],buffer[1],buffer[2],
			outbuffer[0],outbuffer[1],outbuffer[2],outbuffer[3]);
		jar += outbuffer[0];
		jar += outbuffer[1];
		jar += outbuffer[2];
		jar += outbuffer[3];
		fourcount +=1;
		if ( fourcount > linewidth ) {
			fourcount = 0;
			//jar += '\n';
		}
	}
	//jar += '\n';
}

/// @brief Given 3*N bytes of memory to fill, and a string containing 4*N characters, decode the
/// characters (interpreting 4 bytes of ASCII text as 3 binary bytes) and populate the block of memory.
/// @param[in] memory A pointer to the first byte of memory to fill.  It is assumed that we're filling a
/// contiguous block of memory.
/// @param[in] jar The string containing the characters that will be decoded and converted to bytes.
/// @param[in] maxbytes The maximum number of bytes to put into the memory pointed to by the "memory" pointer.
/// (i.e. The size of the array that we're filling).
/// @note Assumes memory already allocated!!!  There is no direct check for vector overflows, since this function
/// has no knowlege of what it is that it's filling (or how big the object is).  The function relies on maxbytes
/// to prevent overflows.
/// @returns The number of bytes filled.
/// @author Originally author unknown.  Revised in 2021 by Vikram K. Mulligan (vmulligan@flatironinstitute.org).
platform::Size
decode6bit(
	unsigned char * memory,
	std::string const & jar,
	platform::Size const maxbytes
) {
	//printf("-->%s\n",jar.c_str());
	const unsigned char *jarmemory = (unsigned char *)jar.c_str();
	platform::Size mempos = 0;
	//unsigned int memlength = (unsigned int)jar.length() * 3 / 4 + 1;
	//*memory = new unsigned char [memlength];
	for ( platform::Size i=0; i<jar.size(); ) {
		unsigned char inbuf[4] = {0,0,0,0};
		unsigned char ibuf=0;
		for ( ; ((i<jar.size())&&(ibuf<4)); i++ ) {
			inbuf[ibuf] = jarmemory[i];
			if ( inbuf[ibuf] < 32 ) continue;
			ibuf++;
		}
		unsigned short const goodbytes = std::min( static_cast<platform::Size>(3), maxbytes-mempos );

		decode_32_to_24(inbuf[0],inbuf[1],inbuf[2],inbuf[3], &((memory)[mempos]), goodbytes);

		mempos += goodbytes;
		//Stop filling bytes:
		if ( mempos > maxbytes ) {
			break;
		}
	}
	return mempos;
}
///
//////////////

}
