// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/VallChunkFilter.cc
/// @brief  defines two basic chunk filters based on pdb id
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

// package headers
#include <protocols/frag_picker/PdbIdChunkFilter.hh>

// mini
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>


// C++ stuff
#include <string>
#include <map>

#include <protocols/frag_picker/VallChunk.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {

static basic::Tracer trPdbFilter( "protocols.frag_picker.PdbIdChunkFilter" );

void PdbIdChunkFilter::load_pdb_id_from_file(std::string file_name) {

	utility::io::izstream data(file_name.c_str());
	trPdbFilter.Info << "read PDB ids from: " << file_name
		<< std::endl;
	if ( !data ) {
		utility_exit_with_message("[ERROR] Unable to open a file fith allowed PDB ids: "
			+ file_name);
	}
	std::string line, pdb_id;
	while ( getline(data, line) ) {
		if ( line.size() && line[0] == '#' ) continue;
		std::istringstream line_stream( line );
		line_stream >> pdb_id;
		add_pdb_id(pdb_id);

		if ( pdb_id.size() > 4 && pdb_id[4] == '_' ) {
			pdb_id[4] = 'A';
			add_pdb_id(pdb_id);
		}
	}
	data.close();
}

void PdbIdChunkFilter::show_pdb_ids(std::ostream& out) {

	core::Size cnt = 0;
	std::map<std::string, bool>::iterator iter;
	out << '\n';
	for ( iter = pdb_hash_.begin(); iter != pdb_hash_.end(); ++iter ) {
		out << iter->first << " ";
		cnt++;
		if ( cnt % 12 == 0 ) {
			out << '\n';
		}
	}
	out << std::endl;
}

bool DenyPdbIdFilter::test_chunk(VallChunkOP a_chunk) {

	std::string id = a_chunk->get_pdb_id() + a_chunk->get_chain_id();
	trPdbFilter.Debug << "Testing " << id << " ... ";
	if ( pdb_hash_.find(id) != pdb_hash_.end() ) {
		trPdbFilter.Debug << "FAILED" << std::endl;
		return false;
	}
	if ( pdb_hash_.find(a_chunk->get_pdb_id()) != pdb_hash_.end() ) {
		trPdbFilter.Debug << "FAILED" << std::endl;
		return false;
	}
	trPdbFilter.Debug << "OK" << std::endl;
	return true;
}

bool AllowPdbIdFilter::test_chunk(VallChunkOP a_chunk) {

	std::string id = a_chunk->get_pdb_id() + a_chunk->get_chain_id();
	trPdbFilter.Debug << "Testing " << id << " ... ";
	if ( pdb_hash_.find(id) != pdb_hash_.end() ) {
		trPdbFilter.Debug << "OK" << std::endl;
		return true;
	}
	if ( pdb_hash_.find(a_chunk->get_pdb_id()) != pdb_hash_.end() ) {
		trPdbFilter.Debug << "OK" << std::endl;
		return true;
	}
	trPdbFilter.Debug << "FAILED" << std::endl;
	return false;
}
} // frag_picker
} // protocols


