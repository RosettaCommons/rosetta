// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

namespace utility {

void encode6bit(const unsigned char* memory, unsigned int length, std::string &jar){
	jar = "";
	unsigned int i;
	unsigned int fourcount=0;
	unsigned int linewidth=15;
	for(i=0;i<length;){
		unsigned char buffer[3] = {0,0,0};
		unsigned char outbuffer[4] = {0,0,0,0};
		int ibuf=0;
		for(;((i<length)&&(ibuf<3));i++){
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
		if(fourcount > linewidth){
			fourcount = 0;
			//jar += '\n';
		}
	}
	//jar += '\n';
}

// assumes memory already allocated!!!!
int decode6bit(unsigned char* memory, const std::string &jar){
	//printf("-->%s\n",jar.c_str());
	const unsigned char *jarmemory = (unsigned char *)jar.c_str();
	unsigned int mempos = 0;
	//unsigned int memlength = (unsigned int)jar.length() * 3 / 4 + 1;
	//*memory = new unsigned char [memlength];
	for(size_t i=0;i<jar.size();){
		unsigned char inbuf[4] = {0,0,0,0};
		unsigned char ibuf=0;
		for(;((i<jar.size())&&(ibuf<4));i++){
			inbuf[ibuf] = jarmemory[i];
			if(inbuf[ibuf] < 32) continue;
			ibuf++;
		}

		decode_32_to_24(inbuf[0],inbuf[1],inbuf[2],inbuf[3],
		//                (*memory)[mempos+0],(*memory)[mempos+1],(*memory)[mempos+2]);
		                (memory)[mempos+0],(memory)[mempos+1],(memory)[mempos+2]);

		mempos += 3;
	}
	return mempos;
}
///
//////////////

}
