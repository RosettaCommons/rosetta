// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @author Alex Ford (fordas@uw.edu)
//

#include <ndarray.h>
#include <utility/exit.hh>

#include <json.hpp>

#include <boost/format.hpp>
#include <boost/range.hpp>
#include <boost/range/combine.hpp>
#include <boost/range/adaptors.hpp>


#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/vector_tools.hh>
#include <protocols/indexed_structure_store/pose_utility.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/types.hh>

static basic::Tracer TR("core.indexed_structure_store.pose_utility");


namespace protocols { namespace indexed_structure_store {

ResidueEntry extract_residue_entry(core::conformation::Residue const & res) {
	using namespace core::id;
	ResidueEntry out;

	out.bb.phi = res.mainchain_torsion(phi_torsion);
	out.bb.psi = res.mainchain_torsion(psi_torsion);
	out.bb.omega = res.mainchain_torsion(omega_torsion);

	out.sc.aa = res.name1();

	int nchi = res.nchi();
	out.sc.chi1 = 1 <= nchi ? res.chi(1) : 0;
	out.sc.chi2 = 2 <= nchi ? res.chi(2) : 0;
	out.sc.chi3 = 3 <= nchi ? res.chi(3) : 0;
	out.sc.chi4 = 4 <= nchi ? res.chi(4) : 0;

	auto N_location = res.xyz("N");
	out.orient.N[0] = N_location.x();
	out.orient.N[1] = N_location.y();
	out.orient.N[2] = N_location.z();

	auto C_location = res.xyz("C");
	out.orient.C[0] = C_location.x();
	out.orient.C[1] = C_location.y();
	out.orient.C[2] = C_location.z();

	auto CA_location = res.xyz("CA");
	out.orient.CA[0] = CA_location.x();
	out.orient.CA[1] = CA_location.y();
	out.orient.CA[2] = CA_location.z();

	auto O_location = res.xyz("O");
	out.orient.O[0] = O_location.x();
	out.orient.O[1] = O_location.y();
	out.orient.O[2] = O_location.z();

	return out;
}

ndarray::Array<ResidueEntry, 1>
extract_residue_entries(
	core::pose::Pose const & pose,
	bool ignore_non_protein
) {
	std::vector<core::Size> rnum;
	for ( std::size_t r = 1; r <=pose.size(); ++r ) {
		if ( pose.residue(r).is_protein() ) {
			rnum.push_back(r);
		} else if ( !ignore_non_protein ) {
			utility_exit_with_message(boost::str(boost::format(
				"Non-protein residue in extract_residue_entries_from_pose residue number: %i residue: %s")
				% r % pose.residue(r)
				));
		}
	}

	std::vector<ResidueEntry> result = v_map<ResidueEntry>(
		rnum,
		[&](core::Size &n){ return extract_residue_entry(pose.residue(n)); }
	);

	for ( auto t: boost::combine(result, rnum) ) {
		boost::get<0>(t).residue_id = boost::get<1>(t);
		boost::get<0>(t).structure_id = 0;
	}

	for ( core::Size c = 0; c < pose.conformation().num_chains(); ++c ) {
		result[pose.conformation().chain_end(c+1) - 1].chain_ending = true;
	}

	ndarray::Array<ResidueEntry, 1> result_array(result.size());
	std::copy(result.begin(), result.end(), result_array.begin());

	return result_array;
}


void
apply_residue_entries_to_pose(
	ndarray::Array<ResidueEntry, 1> residue_entries,
	core::pose::Pose & pose, core::Size start_residue,
	bool apply_bb, bool apply_sidechain, bool apply_orient
) {
	// Unpack and validate entry length
	core::Size num_entries = residue_entries.getSize<0>();
	runtime_assert(num_entries + start_residue - 1 <= pose.total_residue());

	if ( apply_bb ) {
		for ( core::Size i = 0; i < num_entries; ++i ) {
			pose.set_phi(i + start_residue, residue_entries[i].bb.phi);
			pose.set_psi(i + start_residue, residue_entries[i].bb.psi);
			pose.set_omega(i + start_residue, residue_entries[i].bb.omega);

			if ( residue_entries[i].chain_ending ) {
				if ( std::find(
						pose.conformation().chain_endings().begin(),
						pose.conformation().chain_endings().end(),
						i+start_residue ) != pose.conformation().chain_endings().end() ) {}
				else if ( i+start_residue == pose.total_residue() ) {}
				else {
					utility_exit_with_message("Applied chain_ending residue entry to non-chain_ending pose residue.");
				}
			}
		}
	}

	if ( apply_orient ) {
		core::id::AtomID_Mask rebuild_mask(false);
		rebuild_mask.resize(pose.total_residue());
		for ( core::Size r = 1; r <= pose.total_residue(); ++r ) {
			rebuild_mask.resize(r, pose.residue(r).natoms());
		}

		for ( core::Size i = 0; i < num_entries; ++i ) {
			core::conformation::Residue new_residue(pose.residue(i + start_residue));

			// Orient current atom coordinates onto backbone stub, sets sidechain atomic coords
			new_residue.orient_onto_location(
				new_residue.atom_index("CA"),
				new_residue.atom_index("N"),
				new_residue.atom_index("C"),
				core::Vector(residue_entries[i].orient.CA.data()),
				core::Vector(residue_entries[i].orient.N.data()),
				core::Vector(residue_entries[i].orient.O.data())
			);

			// Directly set backbone coords
			new_residue.set_xyz(
				new_residue.atom_index("CA"),
				core::Vector(residue_entries[i].orient.CA.data())
			);

			new_residue.set_xyz(
				new_residue.atom_index("O"),
				core::Vector(residue_entries[i].orient.O.data())
			);

			new_residue.set_xyz(
				new_residue.atom_index("C"),
				core::Vector(residue_entries[i].orient.C.data())
			);

			new_residue.set_xyz(
				new_residue.atom_index("N"),
				core::Vector(residue_entries[i].orient.N.data())
			);

			// Check for not-a-number coordinates, which may occur during invalid orient operations
			runtime_assert(!std::isnan(new_residue.xyz(1).x()));

			pose.replace_residue(i + start_residue, new_residue, false);

			// Identify all non-oriented backbone atoms for rebuild
			for ( auto bb_atom : new_residue.all_bb_atoms() ) {
				rebuild_mask[core::id::AtomID(bb_atom, start_residue + i)] = true;
			}
			rebuild_mask[core::id::AtomID(new_residue.atom_index("N"), start_residue + i)] = false;
			rebuild_mask[core::id::AtomID(new_residue.atom_index("CA"), start_residue + i)] = false;
			rebuild_mask[core::id::AtomID(new_residue.atom_index("C"), start_residue + i)] = false;
			rebuild_mask[core::id::AtomID(new_residue.atom_index("O"), start_residue + i)] = false;
		}

		pose.conformation().fill_missing_atoms(rebuild_mask);
	}

	if ( apply_sidechain ) {
		for ( core::Size i = 0; i < num_entries; ++i ) {

			core::Size resi = i + start_residue;

			if ( pose.residue(resi).name1() != residue_entries[i].sc.aa ) {
				using namespace core::chemical;
				ResidueTypeCOP new_restype(core::pose::get_restype_for_pose(
					pose,
					name_from_aa(aa_from_oneletter_code(residue_entries[i].sc.aa)),
					pose.residue_type(resi).mode())
				);

				// Create the new residue and replace it
				core::conformation::ResidueOP new_res = core::conformation::ResidueFactory::create_residue(
					*new_restype, pose.residue( resi ),
					pose.conformation());
				// Make sure we retain as much info from the previous res as possible
				core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms(
					pose.residue( resi ),
					*new_res,
					pose.conformation(),
					true
				);
				pose.replace_residue( resi, *new_res, false );
			}

			runtime_assert(pose.residue(i + start_residue).name1() == residue_entries[i].sc.aa);
			core::Size nchi = pose.residue(i + start_residue).nchi();
			if ( 1 <= nchi ) {
				pose.set_chi(1, i + start_residue, residue_entries[i].sc.chi1);
			}
			if ( 2 <= nchi ) {
				pose.set_chi(2, i + start_residue, residue_entries[i].sc.chi2);
			}
			if ( 3 <= nchi ) {
				pose.set_chi(3, i + start_residue, residue_entries[i].sc.chi3);
			}
			if ( 4 <= nchi ) {
				pose.set_chi(4, i + start_residue, residue_entries[i].sc.chi4);
			}
		}
	}
}

core::pose::PoseOP
residue_entries_to_pose(
	ndarray::Array<protocols::indexed_structure_store::ResidueEntry, 1> residue_entries,
	std::string residue_type, bool auto_termini
) {

	core::pose::PoseOP work_pose = initial_pose_for_residues(
		residue_entries, residue_type, auto_termini);

	apply_residue_entries_to_pose(residue_entries, *work_pose, 1, true, true, true);

	return work_pose;
}

template <typename Range>
std::string to_str(Range range) {
	return std::string(boost::begin(range), boost::end(range));
}

template<typename Range, typename Op>
std::vector<size_t>
index_if(Range const & r, Op & pred) {
	std::vector<size_t> result;

	for ( auto & pair : r | boost::adaptors::indexed(0) ) {
		if ( pred(pair.value()) ) {
			result.push_back(pair.index());
		}
	}

	return result;
}

core::pose::PoseOP
initial_pose_for_residues(
	ndarray::Array<protocols::indexed_structure_store::ResidueEntry, 1> residue_entries,
	std::string residue_type, bool auto_termini
) {
	typedef std::pair<std::size_t, std::size_t> Span;
	std::vector<Span> chain_spans;
	size_t start = 0;
	for ( std::size_t r = 0; r < residue_entries.getSize<0>(); ++r ) {
		if ( residue_entries[r].chain_ending ) {
			chain_spans.push_back(Span(start, r + 1));
			start = r + 1;
		}
	}
	if ( start != residue_entries.getSize<0>() ) {
		chain_spans.push_back(Span(start, residue_entries.getSize<0>()));
	}

	core::pose::PoseOP work_pose(new core::pose::Pose());

	std::vector<char> chain_chars;
	for ( auto & entry: residue_entries[ndarray::view(chain_spans[0].first, chain_spans[0].second)] ) {
		chain_chars.push_back(entry.sc.aa);
	}

	std::string chain_seq(chain_chars.begin(), chain_chars.end());
	core::pose::make_pose_from_sequence(
		*work_pose, chain_seq, residue_type, auto_termini);

	for ( Span chain_span : boost::make_iterator_range(chain_spans.begin() + 1, chain_spans.end()) ) {
		std::vector<char> chain_chars2;
		for ( auto & entry: residue_entries[ndarray::view(chain_span.first, chain_span.second)] ) {
			chain_chars2.push_back(entry.sc.aa);
		}
		std::string chain_seq2(chain_chars2.begin(), chain_chars2.end());

		core::pose::Pose chain_pose;
		core::pose::make_pose_from_sequence(
			chain_pose, chain_seq2, residue_type, auto_termini);
		work_pose->append_pose_by_jump(chain_pose, 1);
	}

	work_pose->pdb_info(utility::pointer::make_shared< core::pose::PDBInfo >(*work_pose));

	for ( size_t i = 0; i < work_pose->total_residue(); ++i ) {
		work_pose->pdb_info()->number(i + 1, residue_entries[i].residue_id);
	}

	return work_pose;
}


} }
