// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/FragmentCandidate.hh
/// @brief Something that might become a fragment if its scores will be good enough
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallProvider.hh>

#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>

// core
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/fragment/util.hh>
#include <core/pose/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>

// Project headers
#include <core/scoring/rms_util.hh>
#include <numeric/model_quality/rms.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>


// Utility
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {

static THREAD_LOCAL basic::Tracer trFragmentCandidate(
	"protocols.frag_picker.FragmentCandidate");

const std::string FragmentCandidate::unknown_pool_name_ = "UNKNOWN_POOL_NAME";


utility::vector1<FragmentCandidateOP> read_fragment_candidates(
	std::string file_name, VallProviderOP chunk_owner, Size max_nfrags_per_pos) {

	utility::vector1<FragmentCandidateOP> candidates;

	utility::io::izstream data(file_name.c_str());
	trFragmentCandidate.Info << "Reading fragments  from: " << file_name
		<< std::endl;

	if ( !data ) {
		utility_exit_with_message(
			"[ERROR] Unable to open a file with fragments: " + file_name);
	}

	std::string pdb_id = "";
	char chain_id;
	Size res_id = 0, qpos = 1, n_res = 0;
	std::string line;
	std::string tmp;
	Size n_frags = 0;
	while ( data ) {
		getline(data, line);
		Size found = line.find_first_not_of(" \t");
		if ( found == std::string::npos ) { // the line is empty!
			if ( n_res > 0 ) {
				Size vpos = 0;
				VallChunkOP chunk = chunk_owner->find_chunk(pdb_id, chain_id,
					res_id);
				if ( chunk != nullptr ) {
					for ( Size j = 1; j <= chunk->size(); ++j ) {
						if ( chunk->at(j)->resi() == res_id ) {
							vpos = j;
							if ( vpos > 99999 ) {
								trFragmentCandidate.Warning
									<< "Supprisingly high residue id for a vall postion: "
									<< vpos << std::endl;
							}
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
				if ( vpos == 0 ) {
					trFragmentCandidate.Warning << "Can't find a residue: "
						<< res_id << " within a chunk" << std::endl;
					n_res = 0;
					pdb_id = "";
					continue;
				}
				FragmentCandidateOP c( new FragmentCandidate(qpos, vpos,
					chunk, n_res) );
				if ( n_frags < max_nfrags_per_pos ) candidates.push_back(c);
				++n_frags;
				n_res = 0;
				pdb_id = "";
			}
			continue;
		}
		if ( (line.substr(0, 9) == "Position:") || (line.substr(0, 9)
				== "position:") || (line.substr(0, 9) == " Position")
				|| (line.substr(0, 9) == " position") ) {
			std::istringstream line_stream(line);
			line_stream >> tmp >> qpos;
			if ( n_frags > 0 ) {
				trFragmentCandidate.Info << " ... " << n_frags << " found"
					<< std::endl;
			}
			trFragmentCandidate.Info << "Reading fragments for a position: "
				<< qpos << " in a query";
			n_frags = 0;
			continue;
		}
		n_res++;
		if ( (pdb_id.size() == 4) && (res_id > 0) && (res_id < 99999) ) {
			continue;
		}
		std::istringstream line_stream(line);
		line_stream >> pdb_id >> chain_id >> res_id;
		if ( pdb_id.size() != 4 ) {
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


void FragmentCandidate::print_fragment_index(std::ostream& out, bool vall_index_database_exists) {
	if ( vall_index_database_exists ) {
		VallResidueOP r = get_residue(1);
		out << r->key() << " " << fragmentLength_ << std::endl;
	} else {
		out << "0 " << fragmentLength_ << std::endl;
		for ( Size i = 1; i <= fragmentLength_; ++i ) {
			VallResidueOP r = get_residue(i);
			out << r->aa() << " " << r->ss() << " " << utility::trim(F(9, 3, r->phi())) << " " << F(9, 3, r->psi()) << " " << F(9, 3, r->omega()) << std::endl;
		}
	}
}

/// @brief Prints fragment data, the output can be directly loaded to minirosetta
void FragmentCandidate::print_fragment(std::ostream& out, scores::FragmentScoreMapOP sc, scores::FragmentScoreManagerOP ms) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	bool if_ca_in_output = false;
	//bool score_in_output = false;
	if ( option[frags::write_ca_coordinates].user() ) {
		//  std::cerr << "pf 2" << std::endl;
		if_ca_in_output = option[frags::write_ca_coordinates]();
	}
	// optionally add the total score as a comment line.
	// code for reading this score is in src/core/fragment/ConstantLengthFragSet.cc
	if ( option[frags::write_scores].user() && option[frags::write_scores]() && sc && ms ) {
		out << "# score " << F(9,3, ms->total_score(sc)) << std::endl;
	}
	for ( Size i = 1; i <= fragmentLength_; ++i ) {
		//  std::cerr << "pf 4 " << i << std::endl;
		VallResidueOP r = get_residue(i);
		if ( !r ) continue;
		//  std::cerr << "pf 5 " << i << std::endl;
		char aa_upper( toupper(r->aa()) );
		out << " " << get_pdb_id() << " " << get_chain_id() << " " << I(5,r->resi())
			<<" " << aa_upper << " " << r->ss() << F(9, 3, r->phi())
			<< F(9, 3, r->psi()) << F(9, 3, r->omega());
		if ( if_ca_in_output ) {
			out << F(9, 3,r->x()) << F(9, 3, r->y()) << F(9, 3, r->z());
		}
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

/// @brief Prints fragment data, the output can be directly loaded to minirosetta
void FragmentCandidate::print_fragment_seq(std::ostream& out) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	out << " " << get_pdb_id() << " " << get_chain_id();
	VallResidueOP r1 = get_residue(1);
	out << " " << I(5,r1->resi()) << " ";
	for ( Size i = 1; i <= fragmentLength_; ++i ) {
		VallResidueOP r = get_residue(i);
		char aa_upper( toupper(r->aa()) );
		out << aa_upper;
	}
}

/// @brief Prints fragment to silent struct
void FragmentCandidate::output_silent(core::io::silent::SilentFileData & sfd, std::string const & sequence, std::string const & silent_file_name, std::string const & tag, scores::FragmentScoreMapOP sc, scores::FragmentScoreManagerOP ms) {

	using namespace ObjexxFCL;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//if (option[frags::write_ca_coordinates].user())
	//    if_ca_in_output = option[frags::write_ca_coordinates]();

	pose::Pose pose;
	// make pose from frag
	utility::vector1<fragment::FragDataCOP> fragdata;
	fragdata.push_back(get_frag_data());
	fragment::make_pose_from_frags( pose, sequence, fragdata, false );

	if ( option[frags::score_output_silent]() ) {
		pose::Pose relax_pose = pose;

		core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();

		if ( !relax_pose.is_fullatom() ) {
			core::util::switch_to_residue_type_set( relax_pose, core::chemical::FA_STANDARD );
		}


		relax_pose.energies().clear();
		relax_pose.data().clear();

		//   Detect disulfides
		relax_pose.conformation().detect_disulfides();
		relax_pose.conformation().detect_bonds();

		utility::vector1< bool > needToRepack( relax_pose.total_residue(), true );
		core::pack::task::PackerTaskOP taskstd = core::pack::task::TaskFactory::create_packer_task( relax_pose );
		taskstd->restrict_to_repacking();
		taskstd->or_include_current(true);
		taskstd->restrict_to_residues( needToRepack );
		protocols::simple_moves::PackRotamersMover pack1( sfxn, taskstd );
		pack1.apply( relax_pose );

		// quick SC minimization
		core::optimization::AtomTreeMinimizer mzr;
		core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 1e-5, true, false );
		core::kinematics::MoveMapOP mm_min( new core::kinematics::MoveMap() );
		mm_min->set_bb( false );
		mm_min->set_chi( true );
		mzr.run( relax_pose, *mm_min, *sfxn, options );

		Real sc_min_score = (*sfxn)( relax_pose );
		core::pose::setPoseExtraScore( pose, "sc_min_score", sc_min_score );

		// RELAX fragment
		// setup relax protocol for pose
		protocols::relax::RelaxProtocolBaseOP sub_pose_relax_protocol = protocols::relax::generate_relax_from_cmd();
		kinematics::MoveMapOP mm = sub_pose_relax_protocol->get_movemap();
		sub_pose_relax_protocol->apply( relax_pose );

		// check relaxed pose
		Real relaxed_rmsd = scoring::CA_rmsd( relax_pose, pose );
		Real relaxed_score = (*sfxn)( relax_pose );

		core::pose::setPoseExtraScore( pose, "rlx_rms", relaxed_rmsd );
		core::pose::setPoseExtraScore( pose, "rlx_score", relaxed_score );
	}

	// add fragment info and scores
	VallResidueOP r = get_residue(1);
	core::pose::add_score_line_string( pose, "query_pos", string_of(get_first_index_in_query()) ); // query position
	core::pose::add_score_line_string( pose, "vall_pos", string_of(r->resi()) ); // vall position
	core::pose::add_score_line_string( pose, "pdbid", get_pdb_id() ); // pdbid
	core::pose::add_score_line_string( pose, "c", string_of(get_chain_id()) ); // chain id
	core::pose::add_score_line_string( pose, "ss", string_of(get_middle_ss()) ); // secondary structure

	bool if_quota = false;
	if ( sc->get_quota_score() < 999.98 ) if_quota = true;

	for ( Size i = 1; i <= sc->size(); i++ ) {
		scores::FragmentScoringMethodOP ms_i = ms->get_component(i);
		core::pose::setPoseExtraScore( pose, ms_i->get_score_name(), sc->at(i) );
	}
	if ( if_quota ) {
		core::pose::setPoseExtraScore( pose, "QUOTA_TOT", sc->get_quota_score());
		core::pose::setPoseExtraScore( pose, "TOTAL", ms->total_score(sc));
		core::pose::add_score_line_string( pose, "POOL_NAME", get_pool_name());
	} else {
		core::pose::setPoseExtraScore( pose, "TOTAL", ms->total_score(sc));
	}
	debug_assert ( key() > 0 );
	debug_assert ( key() < 40000000 );
	core::pose::add_score_line_string( pose, "FRAG_ID", string_of(key()));

	core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( pose );
	ss->fill_struct( pose, tag );
	sfd.write_silent_struct( *ss, silent_file_name );
}

bool FragmentCandidate::same_chain( FragmentCandidateCOP fr ) {
	if ( get_pdb_id() != fr->get_pdb_id() || get_chain_id() != fr->get_chain_id() ) return false;
	return true;
}

} // frag_picker
} // protocols

