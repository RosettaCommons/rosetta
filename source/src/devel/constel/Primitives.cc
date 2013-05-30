// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief Definition of primitive functions and classes used by the constel program.
/// @author jk
/// @author Andrea Bazzoli

#include <devel/constel/Primitives.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/graph/Graph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/util.hh>
#include <core/chemical/AA.hh>
#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <numeric/xyzVector.hh>
#include <core/io/pdb/file_data.hh>
#include <basic/Tracer.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <iomanip>
#include <sstream>
#include <string>

using core::Real;


static basic::Tracer TR("devel.constel.Primitives");

namespace devel {
namespace constel {

///
/// @brief Returns the list of amino acid types that a given amino acid type can
/// 	be reduced to.
///
/// @param[in] starting_aa the given amino acid type.
///
utility::vector1<char> list_allowable_mutations( char const starting_aa ) {

	utility::vector1<char> allowed_list;

	if ( starting_aa == 'G' ) return allowed_list;

	// general cases
	// anything other than Gly can become Gly, anything other than Gly/Ala can become Ala
	allowed_list.push_back('G');
	if ( starting_aa != 'A' ) {
		allowed_list.push_back('A');
	}

	// special cases
	// Thr can become Ser
	if ( starting_aa == 'T' ) {
		allowed_list.push_back('S');
	}
	// Ile can become Val
	if ( starting_aa == 'I' ) {
		allowed_list.push_back('V');
	}
	// Tyr can become Phe and Leu
	if ( starting_aa == 'Y' ) {
		allowed_list.push_back('F');
		allowed_list.push_back('L');
	}
	// Phe can become Leu
	if ( starting_aa == 'F' ) {
		allowed_list.push_back('L');
	}
	// Trp can become Leu
	if ( starting_aa ==	'W' ) {
		allowed_list.push_back('L');
	}

	return allowed_list;
}


///
/// @brief Sets occupancy to zero for a residue's non-constellation atoms.
///
/// @details Given a residue to be reduced to a smaller amino acid type, sets
/// 	occupancy to zero for any atoms that don't need to be printed (those
/// 	forming the new residue too) and sets occupancy to 1 for any atoms
/// 	that *do* need to be printed (those forming the deleted constellation).
///
/// @param[out] pose pose to which the residue belongs.
/// @param[in] seqpos index of the residue in the pose.
/// @param[in] target_aa amino acid type that the residue has to be mutated
/// 	into.
///
void zero_occ_for_deleted_atoms(Pose & pose, core::Size seqpos,
	char const target_aa) {

	char const starting_aa = core::chemical::oneletter_code_from_aa( pose.aa(seqpos) );

	core::conformation::Residue const & rsd( pose.residue(seqpos) );

	// we never need to print backbone atoms or hydrogens
	for ( Size i=1; i<= rsd.natoms(); ++i ) {
		pose.pdb_info()->occupancy( seqpos, i, 1 ); // to overwrite (rare) zero occ.
		if ( rsd.atom_is_hydrogen(i) || rsd.atom_is_backbone(i) ) {
			pose.pdb_info()->occupancy( seqpos, i, 0 );
		}
	}

	// if it's a mutation to Gly, we're done
	if ( target_aa == 'G' ) return;

	// for anything other than a mutation to Gly, suppress the C-beta
	pose.pdb_info()->occupancy( seqpos, rsd.atom_index("CB"), 0. );

	// if it's a mutation to Ala, we're done
	if ( target_aa == 'A' ) return;

	// I-->V: suppress everything other than CD1
	if ( ( starting_aa == 'I' ) && ( target_aa == 'V' ) ) {
		Size atom_inx_to_keep = rsd.atom_index("CD1");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( i != atom_inx_to_keep )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
		return;
	}

	// T-->S: suppress everything other than CG2
	if ( ( starting_aa == 'T' ) && ( target_aa == 'S' ) ) {
		Size atom_inx_to_keep = rsd.atom_index("CG2");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( i != atom_inx_to_keep )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
		return;
	}

	// Y-->F: suppress everything other than OH
	if ( ( starting_aa == 'Y' ) && ( target_aa == 'F' ) ) {
    		Size atom_inx_to_keep = rsd.atom_index("OH");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( i != atom_inx_to_keep )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
    		return;
  	}

	// Y-->L: suppress everything other than OH, CE1, CE2, CZ
	if ( ( starting_aa == 'Y' ) && ( target_aa == 'L' ) ) {
		Size inx1 = rsd.atom_index("OH");
		Size inx2 = rsd.atom_index("CE1");
		Size inx3 = rsd.atom_index("CE2");
		Size inx4 = rsd.atom_index("CZ");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( (i != inx1) && (i != inx2) && (i != inx3) && (i !=inx4) )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
		return;
	}

	// F-->L: suppress everything other than CE1, CE2, CZ
	if ( ( starting_aa == 'F' ) && ( target_aa == 'L' ) ) {
		Size inx1 = rsd.atom_index("CE1");
		Size inx2 = rsd.atom_index("CE2");
		Size inx3 = rsd.atom_index("CZ");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( (i != inx1) && (i != inx2) && (i != inx3) )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
    		return;
  	}

	// W-->L: suppress everything other than NE1, CE2, CE3, CZ2, CZ3, CH2
	if ( ( starting_aa == 'W' ) && ( target_aa == 'L' ) ) {
		Size inx1 = rsd.atom_index("NE1");
		Size inx2 = rsd.atom_index("CE2");
		Size inx3 = rsd.atom_index("CE3");
		Size inx4 = rsd.atom_index("CZ2");
		Size inx5 = rsd.atom_index("CZ3");
		Size inx6 = rsd.atom_index("CH2");
		for ( Size i=1; i<= rsd.natoms(); ++i )
			if ( (i != inx1) && (i != inx2) && (i != inx3) && (i != inx4) && (i != inx5) && (i != inx6) )
				pose.pdb_info()->occupancy( seqpos, i, 0. );
		return;
	}

	TR << "DANGER DANGER - COULD NOT MAKE REQUESTED MUTATION!!" << std::endl;
	TR << "REQUESTED " << starting_aa << " TO " << target_aa << " BUT NO INFO ON WHAT ATOMS TO SUPPRESS....." << std::endl;
	exit(1);
}


///
/// @brief NeighTeller's constructor
///
/// @param[out] ref_pose pose containing the residues whose neighbors	are to be
/// 	identified.
///
/// @remarks For speed, a scorefxn containing only fa_atr is used.
///
NeighTeller::NeighTeller(Pose& ref_pose) :
  scorefxn(core::scoring::getScoreFunction()),
  fa_atr_weight(0.8), interaction_score_threshold(-0.3) {

	(*scorefxn)(ref_pose);
}


///
/// @brief Tells whether a probe residue is a neighbor of a target residue.
///
/// @param[in] tgt the target residue.
/// @param[in] prb the probe residue.
/// @ref_pose[in] pose both residues belong to.
///
/// @return true if prb is a neighbor of tgt; false otherwise.
///
bool NeighTeller::isneigh(core::conformation::Residue const & tgt,
                          core::conformation::Residue const & prb,
                          Pose const& ref_pose) {

	unweighted_emap.clear();

	core::Vector const tgt_centroid = core::scoring::compute_sc_centroid(tgt);
	core::Vector const prb_centroid = core::scoring::compute_sc_centroid(prb);
	Real const tgt_rad = core::scoring::compute_sc_radius( tgt, tgt_centroid );
	Real const prb_rad = core::scoring::compute_sc_radius( prb, prb_centroid );

	core::scoring::eval_scsc_sr2b_energies( tgt, prb, tgt_centroid, prb_centroid, tgt_rad, prb_rad, ref_pose, *scorefxn, unweighted_emap );

	Real const weighted_fa_atr = fa_atr_weight * unweighted_emap[ core::scoring::fa_atr ];

	return weighted_fa_atr < interaction_score_threshold;
}


///
/// @brief Returns the residue number of a residue in a pose.
///
/// @parm[in] pdbnum residue number of the residue in its PDB file.
/// @parm[in] pdbchn chain identifier of the residue in the PDB file.
/// @parm[in] ps pose that the residue has been loaded into.
///
core::Size get_pose_resnum(int const pdbnum, char const pdbchn, Pose& ps) {

	for ( core::Size j = 1; j <= ps.total_residue(); ++j )
		if ( ( ps.pdb_info()->chain(j) == pdbchn ) &&	(ps.pdb_info()->number(j) == pdbnum) )
			return j;

	// residue not found
	TR << "ERROR!! Could not find residue" << pdbnum << " and chain " << pdbchn << std::endl;
	exit(1);
}


///
/// @brief Creates the list of residues that are neighbors of a given target
/// 	residue.
///
/// @param[in] tgtnum residue number of the target residue in its pose.
/// @param[out] neighs boolean mask filled in by this function to represent the
/// 	list of neighbors. neighs[i] will be true if the ith residue in the pose
/// 	is a neighbor of the target; otherwise neighs[i] will be false.
/// @param[in] ps: pose that the target and its potential neighbors belong to.
///
void mk_neigh_list(core::Size const tgtnum, utility::vector1<bool>& neighs,
	Pose& ps) {

	NeighTeller nt(ps);

	core::conformation::Residue tgtres = ps.residue(tgtnum);

	for(Size i = 1; i <= ps.total_residue(); ++i)
		if(i != tgtnum)
			if(nt.isneigh(tgtres, ps.residue(i), ps))
				neighs.at(i)=true;
}


///
/// @brief Outputs to file a constellation obtained from mutating a pair of
/// 	residues.
///
/// @param[in] mut1 representation of the mutation of the first residue.
/// @param[in] mut2 representation of the mutation of the second residue.
/// @param[in] cslnum a number to identify the constellation.
/// @param[in] ps Rosetta pose that both residues belong to. In the pose,
/// 	the occupancy of atoms in either residue is that AFTER the mutation.
///
/// @remarks
/// 	1. The file names of constellations have the following format,
/// 		SIIIIEC_siiiiec.pdb, where:
/// 	- S is the first residue's start amino acid type
/// 	- IIII is a four-digit field indicating the first residue's number in the
/// 		input PDB file
/// 	- E is the first residue's end amino acid type
/// 	- C is the first residue's chain in the input PDB file
/// 	- s is the second residue's start amino acid type
/// 	- iiii is a four-digit field indicating the second residue's number in the
/// 		input PDB file
/// 	- e is the second residue's end amino acid type
/// 	- c is the second residue's chain in the input PDB file
///
void out_pair_constel(ResMut const& mut1, ResMut const& mut2, int const cslnum, Pose& ps) {

	std::ostringstream outPDB_name;
	outPDB_name.fill('0');

	outPDB_name	<<	"constel_"	<<	mut1.saa	<<	std::setw(4)	<< mut1.pdbn	<<	mut1.eaa	<<	mut1.cid	<<	"_"
  														<<	mut2.saa	<<	std::setw(4)	<< mut2.pdbn	<<	mut2.eaa	<<	mut2.cid	<<	".pdb";

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(outPDB_name.str(), std::ios::out);
	outPDB_stream << "HEADER   CONST NUM "	<<	cslnum	<<
                   " TARGET MUTATION: "	<<	mut1.cid	<<	':'	<<	mut1.saa	<<	mut1.pdbn	<<	mut1.eaa	<<
              "  SECONDARY_MUTATION: "	<<	mut2.cid	<<	':'	<<	mut2.saa	<<	mut2.pdbn <<	mut2.eaa	<<	std::endl;

	core::io::pdb::FileData fd;
	std::string data;
	utility::vector1< core::Size > residues_to_print;
	residues_to_print.push_back(mut1.psn);
	residues_to_print.push_back(mut2.psn);
	fd.init_from_pose( ps, residues_to_print );
	data = core::io::pdb::PDB_DReader::createPDBData(fd);
	outPDB_stream.write( data.c_str(), data.size() );
	outPDB_stream.close();
	outPDB_stream.clear();
}

/// @brief common database for computation of hydrogen bonds
core::scoring::hbonds::HBondDatabaseCOP HBondCommon::hb_database;


///
/// @brief Initializes the data structures of this class.
///
void HBondCommon::init() {

	hb_database = &(*core::scoring::hbonds::HBondDatabase::get_database
                 (core::scoring::hbonds::HBondOptions().params_database_tag()));
}


///
/// @brief Given a residue's moiety in a constellation, returns true if it forms
/// 	at least a hydrogen bond; returns false otherwise.
///
/// @param[in] ps pose to which the constellation belongs.
/// @param[in] cnl indexes in the pose of the residues that contribute atoms to
/// 	the constellation.
///	@param[in] im index in cnl of the residue providing the moiety in question.
/// @param[in] is_donor boolean flag indicating whether the moiety has to be
/// 	considered a donor (true) or an acceptor (false) of hydrogen bonds.
/// @param[in] hb_atoms names of the moiety's atoms whose hydrogen bonds are
/// 	relevant.
///
/// @remarks
///		Assumption: the atom names in 'hb_atoms' are consistent with the value of
/// 		'is_donor'.
/// 	Note: this function neglects hydrogen bonds formed by the moiety and any
/// 		atom in the side chain of a residue that contributes to the
/// 		constellation. Such a neglection is overrestrictive with respect to
/// 		those atoms in the side chain that are not part of the constellation.
///
bool HBondCommon::is_rmoi_hbonded(Pose const& ps,
	utility::vector1<Size> const& cnl, Size const im, bool const is_donor,
  utility::vector1<std::string> const& hb_atoms) {

	using core::scoring::hbonds::identify_hbonds_1way;
	using core::conformation::Residue;

	assert( ps.energies().residue_neighbors_updated() );

	Size const pim = cnl[im];

	// get ensemble of potential hbond partners for moiety
	core::scoring::EnergyGraph const & energy_graph( ps.energies().energy_graph() );
	core::scoring::TenANeighborGraph const & tenA_neighbor_graph( ps.energies().tenA_neighbor_graph() );
	int const nnm = tenA_neighbor_graph.get_node( pim )->num_neighbors_counting_self_static();

	// loop over potential hbond partners and collect hbonds formed by moiety
	core::scoring::hbonds::HBondSet hb_set;
	hb_set.clear();

	Residue const& rm( ps.residue( pim ) );
	for ( core::graph::Graph::EdgeListConstIter
	      nit = energy_graph.get_node(pim)->const_edge_list_begin(),
	      nite = energy_graph.get_node(pim)->const_edge_list_end();
	      nit != nite; ++nit ) {

		Size const pin( (*nit)->get_other_ind(pim) );

		Residue const& rn( ps.residue( pin ) );
		int const nnn = tenA_neighbor_graph.get_node( pin )->num_neighbors_counting_self_static();

		Size const CSIZ = cnl.size();
		bool nei_in_cnl = false;
		for(Size i=1; i<=CSIZ; ++i)
			if(cnl[i] == pin) {
				nei_in_cnl = true;
				break;
			}

		if(is_donor) {

			if(nei_in_cnl)
				identify_hbonds_1way( *HBondCommon::hb_database, rm, rn, nnm, nnn, false,
				                       true, false, true, true, hb_set);
			else
				identify_hbonds_1way( *HBondCommon::hb_database, rm, rn, nnm, nnn, false,
                               true, false, true, false, hb_set);
		}
		else {
			if(nei_in_cnl)
				identify_hbonds_1way( *HBondCommon::hb_database, rn, rm, nnn, nnm, false,
		                           true, true, false, true, hb_set);
			else
				identify_hbonds_1way( *HBondCommon::hb_database, rn, rm, nnn, nnm, false,
                               true, true, false, false, hb_set);
		}
	}

	// check that at least one relevant hydrogen bond is present
	Size N_HB_ATOMS = hb_atoms.size();
	for(Size i=1; i<=N_HB_ATOMS; ++i)
		if( !rm.has( hb_atoms[i] ) )
			return false;

	for(Size i=1; i<=N_HB_ATOMS; ++i) {
		Size aidx = rm.atom_index( hb_atoms[i] );
		core::id::AtomID aid(aidx, pim);
		if(hb_set.atom_hbonds(aid).size())
			return true;
	}

	return false;
}


///
/// @brief Tells whether atom 'low' is closer to atom 'tgt' than atom 'hi' is.
///
/// @param[in] low Cartesian coordinates of a first atom.
/// @param[in] hi Cartesian coordinates of a second atom.
/// @param[in] tgt Cartesian coordinates of a third, target atom.
///
/// @return true if the distance of 'low' to 'tgt' is less than the distance of
/// 	'hi' to 'tgt'; false otherwise.
///
bool OrientCommon::is_closer_to_tgt(numeric::xyzVector<Real> const& low,
	numeric::xyzVector<Real> const& hi, numeric::xyzVector<Real> const& tgt) {

	Real dl = tgt.distance_squared( low );
	Real dh = tgt.distance_squared( hi );

	return dl < dh;
}


///
/// @brief Records the presence of given amino acid types in a given
/// 	constellation.
///
/// @param[in] ps a pose containing a target constellation.
/// @param[in] cnl indexes in 'ps' of the residues contributing to the
/// 	constellation.
/// @param[in] aa_typs desired amino acid types for the residues contributing to
/// 	 the constellation.
/// @param[out] aa_idxs vector that will map the desired amino acid types to the
/// 	residue indexes in 'cnl'.
///
/// @return false if at least one amino acid type in 'aa_typs' is absent from
/// 	the	constellation. Returns true otherwise; in this case, the ith value of
/// 	 'aa_idxs' is the least index in 'cnl' of a residue having the ith amino
/// 	 acid type in 'aa_typs' (i=1,...,aa_typs.size()).
///
/// @remarks It is assumed that:
/// 	1. 'cnl', 'aa_typs', and 'aa_idxs' have the same size.
///   2.  different elements in 'aa_typs' have different values.
///
bool PresenceCommon::are_aa_pres(core::pose::Pose const& ps,
                                 utility::vector1<Size> const& cnl,
                                 utility::vector1<core::chemical::AA> const& aa_typs,
                                 utility::vector1<Size>& aa_idxs) {

	const Size N = aa_typs.size();
	for(Size i=1; i<=N; ++i) {

		core::chemical::AA tgt_aa = aa_typs[i];

		Size j;
		for(j=1; j<=N; ++j)
			if(ps.aa(cnl[j]) == tgt_aa) {
				aa_idxs[i] = j;
				break;
			}

		if(j>N)
			return false;
	}

	return true;
}


///
/// @brief Records the presence of given atoms in a given residue.
///
/// @param[in] res: a residue.
/// @param[in] anams vector of atom names to be found in the residue.
/// @param[out] aidxs vector that will map the desired atom names to their
/// 	indexes in the residue.
///
/// @return false if at least one atom with a desired name is absent from
/// 	the residue. Returns true otherwise; in this case, the ith element of
/// 	'aidxs' will hold the index in 'res' of the atom having the ith name
/// 	in 'anams' (i=1,...,anams.size()).
///
bool PresenceCommon::are_atoms_pres(core::conformation::Residue const& res,
                                    utility::vector1<std::string> const& anams,
                                    utility::vector1<Size>& aidxs) {

	const Size N = anams.size();
	for(Size i=1; i<=N; ++i)
		if(res.has(anams[i]))
			aidxs[i] = res.atom_index(anams[i]);
		else
			return false;

	return true;
}

} // constel
} // devel
