// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constel/Primitives.cc
/// @brief Definition of primitive functions and classes used by the constel program.
/// @author jk
/// @author Andrea Bazzoli

#include <protocols/constel/Primitives.hh>
#include <protocols/constel/SingResCnlCrea.hh>
#include <protocols/constel/MasterFilter.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <utility/graph/Graph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/chemical/AA.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.hh>
#include <numeric/xyzVector.hh>
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>
#include <iomanip>
#include <sstream>
#include <string>

static basic::Tracer TR( "protocols.constel.Primitives" );

namespace protocols {
namespace constel {

using core::Real;
using core::Size;
using core::pose::Pose;

/// @brief Outputs all pair-constellations between a given pair of residues
///
/// @param[in] i pose index of the 1st residue.
/// @param[in] j pose index of the 2nd residue.
/// @param[in] pose_init the pose.
///
/// @remark Output files have the names and format specified in out_pair_constel()
///
void pair_constel_set_idx2(core::Size const i, core::Size const j, Pose const& pose_init) {

	// get info about either residue
	char aai = core::chemical::oneletter_code_from_aa(pose_init.aa(i));
	int i_pdb_number = pose_init.pdb_info()->number(i);
	std::string i_pdb_chain = pose_init.pdb_info()->chain(i);
	utility::vector1<char> allowable_primary_mutations =
		SingResCnlCrea::list_allowable_mutations(aai);

	char aaj = core::chemical::oneletter_code_from_aa(pose_init.aa(j));
	int j_pdb_number = pose_init.pdb_info()->number(j);
	std::string j_pdb_chain = pose_init.pdb_info()->chain(j);
	utility::vector1<char> allowable_secondary_mutations =
		SingResCnlCrea::list_allowable_mutations(aaj);

	utility::vector1<core::Size> cnl;
	cnl.push_back(i);
	cnl.push_back(j);

	// loop over all allowed combinations of mutations for the residue pair
	core::Size const N1MUT = allowable_primary_mutations.size();
	core::Size const N2MUT = allowable_secondary_mutations.size();
	for ( core::Size imut=1; imut <= N1MUT; ++imut ) {

		char aa_imut = allowable_primary_mutations.at(imut);
		Pose primary_mut_pose = pose_init;
		SingResCnlCrea::zero_occ_for_deleted_atoms( primary_mut_pose, i, aa_imut);

		for ( core::Size jmut=1; jmut <= N2MUT; ++jmut ) {

			char aa_jmut = allowable_secondary_mutations.at(jmut);
			Pose secondary_mut_pose = primary_mut_pose;
			SingResCnlCrea::zero_occ_for_deleted_atoms( secondary_mut_pose, j, aa_jmut);

			if ( MasterFilter::is_constel_valid( secondary_mut_pose, cnl ) ) {

				// print the atoms that would be removed by these mutations to a pdb file
				ResMut mut1(aai, aa_imut, i_pdb_chain, i_pdb_number, i);
				ResMut mut2(aaj, aa_jmut, j_pdb_chain, j_pdb_number, j);
				out_pair_constel(mut1, mut2, secondary_mut_pose);
			}
		}
	}
}


/// @brief concatenates into a string all the PDB ATOM records of a constellation
///
/// @param[in] ps pose containing the constellation
/// @param[in] residues_to_print pose indexes of the residues forming the constellation
///
std::string create_constel_records_from_sfr(
	Pose const & ps,
	utility::vector1< core::Size > const & residues_to_print ) {

	core::io::pose_to_sfr::PoseToStructFileRepConverter pose_to_sfr;
	pose_to_sfr.init_from_pose( ps, residues_to_print );

	utility::vector1< core::io::pdb::Record > records( core::io::pdb::create_records_from_sfr(
		*pose_to_sfr.sfr(), utility::pointer::make_shared< core::io::StructFileRepOptions >() ) );

	std::string pdb_contents;
	for ( core::Size i = 1, imax=records.size(); i <= imax; ++i ) {
		if ( records[ i ]["type"].value == "ATOM  " ) {
			pdb_contents += core::io::pdb::create_pdb_line_from_record( records[ i ] ) + "\n";
		}
	}

	return pdb_contents;
}



/// @brief Outputs to file a constellation obtained from mutating a pair of
///  residues.
///
/// @param[in] mut1 representation of the mutation of the first residue.
/// @param[in] mut2 representation of the mutation of the second residue.
/// @param[in] ps Rosetta pose that both residues belong to. In the pose,
///  the occupancy of atoms in either residue is that AFTER the mutation.
///
/// @remarks
///  1. The file name of the constellation has the following format,
///   SIIIIEC_siiiiec.pdb, where:
///  - S is the first residue's start amino acid type
///  - IIII is a four-digit field indicating the first residue's number in the
///   input PDB file
///  - E is the first residue's end amino acid type
///  - C is the first residue's chain in the input PDB file
///  - s is the second residue's start amino acid type
///  - iiii is a four-digit field indicating the second residue's number in the
///   input PDB file
///  - e is the second residue's end amino acid type
///  - c is the second residue's chain in the input PDB file
///
///  2. This function may set to zero the occupancy of additional atoms of the
///   residues forming the constellation..
///
void out_pair_constel(ResMut const& mut1, ResMut const& mut2, Pose& ps) {

	// prepare file name
	std::ostringstream outPDB_name;
	outPDB_name.fill('0');

	outPDB_name << "constel_" <<
		mut1.saa << std::setw(4) << mut1.pdbn << mut1.eaa << mut1.cid << "_" <<
		mut2.saa << std::setw(4) << mut2.pdbn << mut2.eaa << mut2.cid << ".pdb";

	// print header
	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(outPDB_name.str(), std::ios::out);
	outPDB_stream << "HEADER   " <<
		" 1st MUTATION: " << mut1.cid << ':' << mut1.saa << mut1.pdbn <<
		mut1.eaa <<
		" 2nd MUTATION: " << mut2.cid << ':' << mut2.saa << mut2.pdbn <<
		mut2.eaa << std::endl;

	// check whether atoms must be stripped off the constellation
	if ( SingResCnlCrea::stripped() ) {
		SingResCnlCrea::strip_atoms(ps, mut1.psn, mut1.eaa);
		SingResCnlCrea::strip_atoms(ps, mut2.psn, mut2.eaa);
	}

	// print constellation
	utility::vector1< core::Size > residues_to_print;
	residues_to_print.push_back(mut1.psn);
	residues_to_print.push_back(mut2.psn);
	std::string data = create_constel_records_from_sfr( ps, residues_to_print );
	outPDB_stream.write( data.c_str(), data.size() );
	outPDB_stream.close();
	outPDB_stream.clear();
}


/// @brief Outputs all triple-constellations among a given triple of residues
///
/// @param[in] i pose index of the 1st residue.
/// @param[in] j pose index of the 2nd residue.
/// @param[in] k pose index of the 3rd residue.
/// @param[in] pose_init the pose.
///
/// @remark Output files have the names and format specified in out_triple_constel()
///
void triple_constel_set_idx3(core::Size const i, core::Size const j, core::Size const k,
	Pose const& pose_init ) {

	using utility::vector1;

	// get info about each residue of the triple
	char aai = oneletter_code_from_aa(pose_init.aa(i));
	int i_pdb_number = pose_init.pdb_info()->number(i);
	std::string i_pdb_chain = pose_init.pdb_info()->chain(i);
	vector1<char> allowable_primary_mutations =
		SingResCnlCrea::list_allowable_mutations(aai);

	char aaj = oneletter_code_from_aa(pose_init.aa(j));
	int j_pdb_number = pose_init.pdb_info()->number(j);
	std::string j_pdb_chain = pose_init.pdb_info()->chain(j);
	vector1<char> allowable_secondary_mutations =
		SingResCnlCrea::list_allowable_mutations(aaj);

	char aak = oneletter_code_from_aa(pose_init.aa(k));
	int k_pdb_number = pose_init.pdb_info()->number(k);
	std::string k_pdb_chain = pose_init.pdb_info()->chain(k);
	vector1<char> allowable_tertiary_mutations =
		SingResCnlCrea::list_allowable_mutations(aak);

	vector1<core::Size> cnl;
	cnl.push_back(i);
	cnl.push_back(j);
	cnl.push_back(k);

	// loop over all allowed combinations of mutations for the residue triple
	core::Size const N1MUT = allowable_primary_mutations.size();
	core::Size const N2MUT = allowable_secondary_mutations.size();
	core::Size const N3MUT = allowable_tertiary_mutations.size();

	for ( core::Size imut=1; imut <= N1MUT; ++imut ) {

		char aa_imut = allowable_primary_mutations.at(imut);
		Pose primary_mut_pose = pose_init;
		SingResCnlCrea::zero_occ_for_deleted_atoms( primary_mut_pose, i, aa_imut);

		for ( core::Size jmut=1; jmut <= N2MUT; ++jmut ) {

			char aa_jmut = allowable_secondary_mutations.at(jmut);
			Pose secondary_mut_pose = primary_mut_pose;
			SingResCnlCrea::zero_occ_for_deleted_atoms( secondary_mut_pose, j, aa_jmut);

			for ( core::Size kmut=1; kmut <= N3MUT; ++kmut ) {

				char aa_kmut = allowable_tertiary_mutations.at(kmut);
				Pose tertiary_mut_pose = secondary_mut_pose;
				SingResCnlCrea::zero_occ_for_deleted_atoms( tertiary_mut_pose, k, aa_kmut);

				if ( MasterFilter::is_constel_valid( tertiary_mut_pose, cnl ) ) {

					// print the atoms that would be removed by these mutations to
					// a pdb file
					ResMut mut1(aai, aa_imut, i_pdb_chain, i_pdb_number, i);
					ResMut mut2(aaj, aa_jmut, j_pdb_chain, j_pdb_number, j);
					ResMut mut3(aak, aa_kmut, k_pdb_chain, k_pdb_number, k);
					out_triple_constel(mut1, mut2, mut3, tertiary_mut_pose);
				}
			}
		}
	}
}


/// @brief Outputs to file a constellation obtained from mutating a triple of
///  residues.
///
/// @param[in] mut1 representation of the mutation of the first residue.
/// @param[in] mut2 representation of the mutation of the second residue.
/// @param[in] mut3 representation of the mutation of the third residue.
/// @param[in] ps Rosetta pose that all three residues belong to. In the pose,
///  the occupancy of atoms in each residue is that AFTER the mutation.
///
/// @remarks
///  1. The file name of the constellation has the following format,
///   SIIIIEC_siiiiec_tjjjjfd.pdb, where:
///  - S is the first residue's start amino acid type
///  - IIII is a four-digit field indicating the first residue's number in the
///   input PDB file
///  - E is the first residue's end amino acid type
///  - C is the first residue's chain in the input PDB file
///  - s is the second residue's start amino acid type
///  - iiii is a four-digit field indicating the second residue's number in the
///   input PDB file
///  - e is the second residue's end amino acid type
///  - c is the second residue's chain in the input PDB file
///  - t is the third residue's start amino acid type
///  - jjjj is a four-digit field indicating the third residue's number in the
///   input PDB file
///  - f is the third residue's end amino acid type
///  - d is the third residue's chain in the input PDB file
///
void out_triple_constel(ResMut const& mut1, ResMut const& mut2,
	ResMut const& mut3, Pose& ps) {

	// prepare file name
	std::ostringstream outPDB_name;
	outPDB_name.fill('0');

	outPDB_name << "constel_" <<
		mut1.saa << std::setw(4) << mut1.pdbn << mut1.eaa << mut1.cid << "_" <<
		mut2.saa << std::setw(4) << mut2.pdbn << mut2.eaa << mut2.cid << "_" <<
		mut3.saa << std::setw(4) << mut3.pdbn << mut3.eaa << mut3.cid << ".pdb";

	// print header
	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(outPDB_name.str(), std::ios::out);
	outPDB_stream << "HEADER   " <<
		" 1st MUTATION: " << mut1.cid << ':' << mut1.saa << mut1.pdbn <<
		mut1.eaa <<
		" 2nd MUTATION: " << mut2.cid << ':' << mut2.saa << mut2.pdbn <<
		mut2.eaa <<
		" 3rd MUTATION: " << mut3.cid << ':' << mut3.saa << mut3.pdbn <<
		mut3.eaa << std::endl;

	// check whether atoms must be stripped off the constellation
	if ( SingResCnlCrea::stripped() ) {
		SingResCnlCrea::strip_atoms(ps, mut1.psn, mut1.eaa);
		SingResCnlCrea::strip_atoms(ps, mut2.psn, mut2.eaa);
		SingResCnlCrea::strip_atoms(ps, mut3.psn, mut3.eaa);
	}

	// print constellation
	utility::vector1< core::Size > residues_to_print;
	residues_to_print.push_back(mut1.psn);
	residues_to_print.push_back(mut2.psn);
	residues_to_print.push_back(mut3.psn);
	std::string data = create_constel_records_from_sfr( ps, residues_to_print );
	outPDB_stream.write( data.c_str(), data.size() );
	outPDB_stream.close();
	outPDB_stream.clear();
}


/// @brief common database for computation of hydrogen bonds
core::scoring::hbonds::HBondDatabaseCOP HBondCommon::hb_database;


/// @brief Initializes the data structures of this class.
///
void HBondCommon::init() {

	hb_database = core::scoring::hbonds::HBondDatabase::get_database
		(core::scoring::hbonds::HBondOptions().params_database_tag());
}


/// @brief Given a residue's moiety in a constellation, returns true if it forms
///  at least a hydrogen bond; returns false otherwise.
///
/// @param[in] ps pose to which the constellation belongs.
/// @param[in] cnl indexes in the pose of the residues that contribute atoms to
///  the constellation.
/// @param[in] im index in cnl of the residue providing the moiety in question.
/// @param[in] is_donor boolean flag indicating whether the moiety has to be
///  considered a donor (true) or an acceptor (false) of hydrogen bonds.
/// @param[in] hb_atoms names of the moiety's atoms whose hydrogen bonds are
///  relevant.
///
/// @remarks
///  Assumption: the atom names in 'hb_atoms' are consistent with the value of
///   'is_donor'.
///  Note: this function neglects hydrogen bonds formed by the moiety and any
///   atom in the side chain of a residue that contributes to the
///   constellation. Such a neglection is overrestrictive with respect to
///   those atoms in the side chain that are not part of the constellation.
///
bool HBondCommon::is_rmoi_hbonded(Pose const& ps,
	utility::vector1<core::Size> const& cnl, core::Size const im, bool const is_donor,
	utility::vector1<std::string> const& hb_atoms) {

	using core::scoring::hbonds::identify_hbonds_1way;
	using core::conformation::Residue;

	debug_assert( ps.energies().residue_neighbors_updated() );

	core::Size const pim = cnl[im];

	// get ensemble of potential hbond partners for moiety
	core::scoring::EnergyGraph const & energy_graph( ps.energies().energy_graph() );
	core::scoring::TenANeighborGraph const & tenA_neighbor_graph( ps.energies().tenA_neighbor_graph() );
	int const nnm = tenA_neighbor_graph.get_node( pim )->num_neighbors_counting_self_static();

	// loop over potential hbond partners and collect hbonds formed by moiety
	core::scoring::hbonds::HBondSet hb_set;
	hb_set.clear();

	Residue const& rm( ps.residue( pim ) );
	for ( utility::graph::Graph::EdgeListConstIter
			nit = energy_graph.get_node(pim)->const_edge_list_begin(),
			nite = energy_graph.get_node(pim)->const_edge_list_end();
			nit != nite; ++nit ) {

		core::Size const pin( (*nit)->get_other_ind(pim) );

		Residue const& rn( ps.residue( pin ) );
		int const nnn = tenA_neighbor_graph.get_node( pin )->num_neighbors_counting_self_static();

		core::Size const CSIZ = cnl.size();
		bool nei_in_cnl = false;
		for ( core::Size i=1; i<=CSIZ; ++i ) {
			if ( cnl[i] == pin ) {
				nei_in_cnl = true;
				break;
			}
		}

		if ( is_donor ) {

			if ( nei_in_cnl ) {
				identify_hbonds_1way( *HBondCommon::hb_database, rm, rn, nnm, nnn, false,
					true, false, true, true, hb_set);
			} else {
				identify_hbonds_1way( *HBondCommon::hb_database, rm, rn, nnm, nnn, false,
					true, false, true, false, hb_set);
			}
		} else {
			if ( nei_in_cnl ) {
				identify_hbonds_1way( *HBondCommon::hb_database, rn, rm, nnn, nnm, false,
					true, true, false, true, hb_set);
			} else {
				identify_hbonds_1way( *HBondCommon::hb_database, rn, rm, nnn, nnm, false,
					true, true, false, false, hb_set);
			}
		}
	}

	// check that at least one relevant hydrogen bond is present
	core::Size N_HB_ATOMS = hb_atoms.size();
	for ( core::Size i=1; i<=N_HB_ATOMS; ++i ) {
		if ( !rm.has( hb_atoms[i] ) ) {
			return false;
		}
	}

	for ( core::Size i=1; i<=N_HB_ATOMS; ++i ) {
		core::Size aidx = rm.atom_index( hb_atoms[i] );
		core::id::AtomID aid(aidx, pim);
		if ( hb_set.atom_hbonds(aid).size() ) {
			return true;
		}
	}

	return false;
}


/// @brief Tells whether atom 'low' is closer to atom 'tgt' than atom 'hi' is.
///
/// @param[in] low Cartesian coordinates of a first atom.
/// @param[in] hi Cartesian coordinates of a second atom.
/// @param[in] tgt Cartesian coordinates of a third, target atom.
///
/// @return true if the distance of 'low' to 'tgt' is less than the distance of
///  'hi' to 'tgt'; false otherwise.
///
bool OrientCommon::is_closer_to_tgt(numeric::xyzVector<Real> const& low,
	numeric::xyzVector<Real> const& hi, numeric::xyzVector<Real> const& tgt) {

	Real dl = tgt.distance_squared( low );
	Real dh = tgt.distance_squared( hi );

	return dl < dh;
}


/// @brief Records the presence of given amino acid types in a given
///  constellation.
///
/// @param[in] ps a pose containing a target constellation.
/// @param[in] cnl indexes in 'ps' of the residues contributing to the
///  constellation.
/// @param[in] aa_typs desired amino acid types for the residues contributing to
///   the constellation.
/// @param[out] aa_idxs vector that will map the desired amino acid types to the
///  residue indexes in 'cnl'.
///
/// @return false if at least one amino acid type in 'aa_typs' is absent from
///  the constellation. Returns true otherwise; in this case, the ith value of
///   'aa_idxs' is the least index in 'cnl' of a residue having the ith amino
///   acid type in 'aa_typs' (i=1,...,aa_typs.size()).
///
/// @remarks It is assumed that:
///  1. 'cnl', 'aa_typs', and 'aa_idxs' have the same size.
///   2.  different elements in 'aa_typs' have different values.
///
bool PresenceCommon::are_aa_pres(core::pose::Pose const& ps,
	utility::vector1<core::Size> const& cnl,
	utility::vector1<core::chemical::AA> const& aa_typs,
	utility::vector1<core::Size>& aa_idxs) {

	const core::Size N = aa_typs.size();
	for ( core::Size i=1; i<=N; ++i ) {

		core::chemical::AA tgt_aa = aa_typs[i];

		core::Size j;
		for ( j=1; j<=N; ++j ) {
			if ( ps.aa(cnl[j]) == tgt_aa ) {
				aa_idxs[i] = j;
				break;
			}
		}

		if ( j>N ) {
			return false;
		}
	}

	return true;
}


/// @brief Records the presence of given atoms in a given residue.
///
/// @param[in] res: a residue.
/// @param[in] anams vector of atom names to be found in the residue.
/// @param[out] aidxs vector that will map the desired atom names to their
///  indexes in the residue.
///
/// @return false if at least one atom with a desired name is absent from
///  the residue. Returns true otherwise; in this case, the ith element of
///  'aidxs' will hold the index in 'res' of the atom having the ith name
///  in 'anams' (i=1,...,anams.size()).
///
bool PresenceCommon::are_atoms_pres(core::conformation::Residue const& res,
	utility::vector1<std::string> const& anams,
	utility::vector1<core::Size>& aidxs) {

	const core::Size N = anams.size();
	for ( core::Size i=1; i<=N; ++i ) {
		if ( res.has(anams[i]) ) {
			aidxs[i] = res.atom_index(anams[i]);
		} else {
			return false;
		}
	}

	return true;
}

} // constel
} // protocols
