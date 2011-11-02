// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   protocols/frag_picker/FragmentCandidate.hh
/// @brief Something that might become a fragment if its scores will be good enough
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallProvider.hh>

#include <basic/Tracer.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>


// Utility
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <cctype>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {

static basic::Tracer trFragmentCandidate(
		"protocols.frag_picker.FragmentCandidate");

const std::string FragmentCandidate::unknown_pool_name_ = "UNKNOWN_POOL_NAME";


utility::vector1<FragmentCandidateOP> read_fragment_candidates(
		std::string file_name, VallProviderOP chunk_owner) {

	utility::vector1<FragmentCandidateOP> candidates;

	utility::io::izstream data(file_name.c_str());
	trFragmentCandidate.Info << "Reading fragments  from: " << file_name
			<< std::endl;

	if (!data)
		utility_exit_with_message(
				"[ERROR] Unable to open a file with fragments: " + file_name);

	std::string pdb_id = "";
	char chain_id;
	Size res_id = 0, qpos = 1, n_res = 0;
	std::string line;
	std::string tmp;
	Size n_frags = 0;
	while (data) {
		getline(data, line);
		Size found = line.find_first_not_of(" \t");
		if (found == std::string::npos) { // the line is empty!
			if (n_res > 0) {
				Size vpos = 0;
				VallChunkOP chunk = chunk_owner->find_chunk(pdb_id, chain_id,
						res_id);
				if (chunk != 0) {
					for (Size j = 1; j <= chunk->size(); ++j) {
						if (chunk->at(j)->resi() == res_id) {
							vpos = j;
							if (vpos > 99999)
								trFragmentCandidate.Warning
										<< "Supprisingly high residue id for a vall postion: "
										<< vpos << std::endl;
							break;
						}
					}
				} else {
					trFragmentCandidate.Warning
							<< "Can't find the following chunk in a vall: "
							<< pdb_id << " " << chain_id << " " << res_id
							<< std::endl;
					n_res = 0;
					pdb_id = "";
					continue;
				}
				if (vpos == 0) {
					trFragmentCandidate.Warning << "Can't find a residue: "
							<< res_id << " within a chunk" << std::endl;
					n_res = 0;
					pdb_id = "";
					continue;
				}
				FragmentCandidateOP c = new FragmentCandidate(qpos, vpos,
						chunk, n_res);
				candidates.push_back(c);
				++n_frags;
				n_res = 0;
				pdb_id = "";
			}
			continue;
		}
		if ((line.substr(0, 9) == "Position:") || (line.substr(0, 9)
				== "position:") || (line.substr(0, 9) == " Position")
				|| (line.substr(0, 9) == " position")) {
			std::istringstream line_stream(line);
			line_stream >> tmp >> qpos;
			if (n_frags > 0)
				trFragmentCandidate.Info << " ... " << n_frags << " found"
						<< std::endl;
			trFragmentCandidate.Info << "Reading fragments for a position: "
					<< qpos << " in a query";
			n_frags = 0;
			continue;
		}
		n_res++;
		if ((pdb_id.size() == 4) && (res_id > 0) && (res_id < 99999))
			continue;
		std::istringstream line_stream(line);
		line_stream >> pdb_id >> chain_id >> res_id;
		if (pdb_id.size() != 4) {
			trFragmentCandidate.Info << "Skiping strange PDB ID: " << pdb_id
					<< std::endl;
			pdb_id = "";
			continue;
		}
	}
	trFragmentCandidate.Info << " ... " << n_frags << " found"<< std::endl;
	trFragmentCandidate.Info << candidates.size()<<" fragments read from a file"<<std::endl;

	return candidates;
}

	/// @brief Prints fragment data, the output can be directly loaded to minirosetta
void FragmentCandidate::print_fragment(std::ostream& out) {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		bool if_ca_in_output = false;
		if (option[frags::write_ca_coordinates].user())
			if_ca_in_output = option[frags::write_ca_coordinates]();
		for (Size i = 1; i <= fragmentLength_; ++i) {
			VallResidueOP r = get_residue(i);
			char aa_upper( toupper(r->aa()) );
			out << " " << get_pdb_id() << " " << get_chain_id() << " " << I(5,r->resi())
					<<" " << aa_upper << " " << r->ss() << F(9, 3, r->phi())
					<< F(9, 3, r->psi()) << F(9, 3, r->omega());
			if(if_ca_in_output)
				out << F(9, 3,r->x()) << F(9, 3, r->y()) << F(9, 3, r->z());
			out << std::endl;
		}
	}

/* To be used in the nearest future...
ConstantLengthFragSet FragmentCandidate::create_frag_set(utility::vector1<std::pair<
                        FragmentCandidateOP, scores::FragmentScoreMapOP> >& selected_fragments) {

    std::stringstream s;
    for (Size fi = 1; fi <= out.size(); ++fi) {
        selected_fragments[i].first->print_fragment(output);
        s << std::endl;
    }
} */

bool FragmentCandidate::same_chain( FragmentCandidateCOP fr ) {
	if (get_pdb_id() != fr->get_pdb_id() || get_chain_id() != fr->get_chain_id()) return false;
	return true;
}

} // frag_picker
} // protocols

