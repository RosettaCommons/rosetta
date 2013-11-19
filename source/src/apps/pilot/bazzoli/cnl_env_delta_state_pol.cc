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

///
/// @brief Prints the difference in hydrogen-bonding state of polar atoms
/// 	between	a wild-type structure containing a constellation and a mutant
/// 	structure without the constellation.
///
/// @details The difference is computed with respect to a fixed set of
/// 	environment residues, namely, residues which are neighbors of the
/// 	constellation.
///
/// @param[in] -s <PDBWT>, where <PDBWT> is the path to the PDB file containing
/// 	the wild-type structure.
///
/// @param[in] -mut <PDBMUT>, where <PDBMUT> is the path to the PDB file
/// 	containing the mutant structure.
///
/// @param[in] -cnl_resfile <CNL>, where <CNL> is the path to a file enumerating
/// 	the residues that form the constellation. The file has the following
/// 	format:
/// 	I1 C1\n
///   ...
///   IN CN\n ,
///   where Ii and Ci are the residue index and the chain ID, respectively, of
/// 	the ith residue forming the constellation (i=1,...,N).
///
/// @param[in] -extra_res_fa <PARAMS>. This parameter is required only when the
/// 	mutant structure contains a ligand for which a .params file is needed. In
/// 	such a case, <PARAMS> provides the path to the .params file.
///
/// @details The input constellation file requires that a chain identifier equal
/// 	to ' ' in the PDB file be represented by '_'.
///
/// @author Andrea Bazzoli (bazzoli@ku.edu)
///
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/graph/Graph.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <numeric/xyzVector.hh>
#include <core/chemical/AA.hh>
#include <utility/vector1.hh>
#include <protocols/neighbor/Neighborhood.hh>

#include <string>
#include <fstream>
#include <map>

using namespace basic::options::OptionKeys;
using namespace core::scoring::hbonds;
using core::pose::Pose;
using core::conformation::Residue;
using core::id::AtomID;
using core::Size;
using core::Real;
using utility::vector1;
using numeric::xyzVector;

static basic::Tracer TR( "apps.pilot.cnl_env_delta_state_pol" );


///
/// @brief: prints the identifiers of a set of residues
///
///	@param[in]: vec indexes of the residues in the pose
/// @param[in]: ps the pose
/// @param[in]: tr output tracer
///
/// @details: the ith output line contains the identifier of the residue
///   having the ith index in vec (i=1,...,vec.size()).
///
void print_res_ids(vector1<Size> vec, Pose const& ps, basic::Tracer& tr) {

	for(Size i=1; i<=vec.size(); ++i) {
		Size ri = vec[i];
		tr << ri << '(' << ps.pdb_info()->chain(ri) << ps.pdb_info()->number(ri) <<
			')' << std::endl;
	}
}


///
/// @brief Returns the residue number of a residue in a pose.
///
/// @parm[in] pdbnum residue number of the residue in its PDB file.
/// @parm[in] pdbchn chain identifier of the residue in the PDB file.
/// @parm[in] ps pose that the residue has been loaded into.
///
Size get_pose_resnum(int const pdbnum, char const pdbchn, Pose& ps) {

	for ( Size j = 1; j <= ps.total_residue(); ++j )
		if ( ( ps.pdb_info()->chain(j) == pdbchn ) && (ps.pdb_info()->number(j) == pdbnum) )
			return j;

	// residue not found
	TR << "ERROR!! Could not find residue" << pdbnum << " and chain " << pdbchn << std::endl;
	exit(1);
}


///
/// @brief returns the index in the residue of the atom designated to represent
/// 	a target atom for	hydrogen bonding.
///
/// @param[in] tgt index in the residue of the target atom
/// @param[in] res residue to which the target atom belongs
///
/// @details: the representative atom is used to make some atoms equivalent with
/// 	respect to hydrogen bonding. Currently:
/// 	- all hydrogen atoms bound to the same heavy atom are equivalent.
///   - Arg's 1HH1, 1HH2, 2HH1, and 2HH2 are all equivalent.
///   - Asp's OD1 and OD2 are equivalent.
///   - Glu's OE1 and OE2 are equivalent.
///
Size rep_hb_atom(Size const tgt, Residue const& res) {

	using core::chemical::aa_arg;
	using core::chemical::aa_asp;
	using core::chemical::aa_glu;

	if(res.atom_is_hydrogen(tgt)) {
		if(res.aa() == aa_arg)
			if((res.atom_name(tgt) == "1HH2") || (res.atom_name(tgt) == "2HH2"))
				return res.atom_index("NH1");

		return res.atom_base(tgt);
	}
	if(res.aa() == aa_asp) {
		if(res.atom_name(tgt) == " OD2")
			return res.atom_index("OD1");
	}
	else
		if(res.aa() == aa_glu)
			if(res.atom_name(tgt) == " OE2")
				return res.atom_index("OE1");

	return tgt;
}


///
/// @brief Tells whether two atoms form a hydrogen bond.
///
/// @param[in] don candidate donor atom.
/// @param[in] acc candidate acceptor atom.
/// @param[in] ps pose containing the atoms.
/// @param[in] database hbond database.
/// @param[in] hbond_set dummy hbond set (will not be filled in).
///
/// @return true if the two atoms form a hydrogen bond; false otherwise.
///
/// @details It is assumed that don is a polar hydrogen atom and acc may
/// 	accept hydrogen bonds.
///
bool is_hbond(AtomID const& don, AtomID const& acc, Pose const& ps,
	HBondDatabase const& database, HBondSet const& hbond_set) {

	Residue const& don_rsd = ps.residue(don.rsd());
	Residue const& acc_rsd = ps.residue(acc.rsd());

	Size const hatm = don.atomno();
	Size const aatm = acc.atomno();

	HBondDerivs derivs;

	xyzVector<Real> const& hatm_xyz(don_rsd.atom(hatm).xyz());

	Size const datm(don_rsd.atom_base(hatm));
	xyzVector<Real> const& datm_xyz(don_rsd.atom(datm).xyz());

	xyzVector<Real> const& aatm_xyz(acc_rsd.atom(aatm).xyz());

	if(hatm_xyz.distance_squared(aatm_xyz) <= MAX_R2) {

		Real unweighted_energy(0.0);
		HBEvalTuple hbe_type(datm, don_rsd, aatm, acc_rsd);

		int const base(acc_rsd.atom_base(aatm));
		int const base2(acc_rsd.abase2(aatm));
		assert(base2 > 0 && base != base2);

		core::scoring::hbonds::hb_energy_deriv(database,
			hbond_set.hbond_options(), hbe_type, datm_xyz, hatm_xyz, aatm_xyz,
			acc_rsd.atom(base).xyz(), acc_rsd.atom(base2).xyz(),
			unweighted_energy, false, derivs);

		if(unweighted_energy < MAX_HB_ENERGY)
			return true;
		else
			return false;
	}
	else return false;
}


///
/// @brief computes and stores the number of hydrogen bonds for each
/// 	representative atom of a residue.
///
/// @param[in] ridx index of the residue in the pose
/// @param[in] ps the pose
/// @param[in] database hbond database.
/// @param[in] hbond_set dummy hbond set (will not be filled in).
/// @param[out] dcounts storage for the donor atoms' hydrogen bond counts
/// @param[out] acounts storage for the acceptor atoms' hydrogen bond counts
///
/// @details For each representative atom r of donor atoms,
/// 	dcounts[AtomID(r, ridx)] will contain the total number of hydrogen bonds
/// 	made by the donor atoms it represents.
///
/// @details For each representative atom s of acceptor atoms,
/// 	acounts[AtomID(s, ridx)] will contain the total number of hydrogen bonds
/// 	made by the acceptor atoms it represents.
///
/// @details It is assumed that neither dcounts nor acounts contain, when
/// 	the function begins, any item for residue ridx.
///
void store_hb_counts(Size const ridx, Pose const& ps,
	HBondDatabase const& database, HBondSet const& hbond_set,
	std::map<AtomID, Size>& dcounts, std::map<AtomID, Size>& acounts) {

	Residue const& res = ps.residue(ridx);

	core::scoring::EnergyGraph const& energy_graph(ps.energies().energy_graph());
	core::graph::Graph::EdgeListConstIter NITB =
		energy_graph.get_node(ridx)->const_edge_list_begin();
	core::graph::Graph::EdgeListConstIter NITE =
		energy_graph.get_node(ridx)->const_edge_list_end();

	// collect hbonds for donors
	for(core::chemical::AtomIndices::const_iterator
		hnum = res.Hpos_polar().begin(), hnume = res.Hpos_polar().end();
		hnum != hnume; ++hnum) {

		Size hatm(*hnum);
		AtomID don(hatm, ridx);

		Size ratm = rep_hb_atom(hatm, res);
		AtomID rep(ratm, ridx);
		if(dcounts.find(rep) == dcounts.end())
			dcounts[rep] = 0;

		for(core::graph::Graph::EdgeListConstIter nit = NITB; nit != NITE; ++nit) {

			Size nidx((*nit)->get_other_ind(ridx));
			Residue const& ngb = ps.residue(nidx);

			for(core::chemical::AtomIndices::const_iterator
				anum = ngb.accpt_pos().begin(), anume = ngb.accpt_pos().end();
				anum != anume; ++anum) {

				Size aatm(*anum);
				AtomID acc(aatm, nidx);

				if(is_hbond(don, acc, ps, database, hbond_set))
					dcounts[rep]++;
			}
		}
	}

	// collect hbonds for acceptors
	for(core::chemical::AtomIndices::const_iterator
		anum = res.accpt_pos().begin(), anume = res.accpt_pos().end();
		anum != anume; ++anum) {

		Size aatm(*anum);
		AtomID acc(aatm, ridx);

		Size ratm = rep_hb_atom(aatm, res);
		AtomID rep(ratm, ridx);
		if(acounts.find(rep) == acounts.end())
			acounts[rep] =  0;

		for(core::graph::Graph::EdgeListConstIter nit = NITB; nit != NITE; ++nit) {

			Size nidx((*nit)->get_other_ind(ridx));
			Residue const& ngb = ps.residue(nidx);

			for(core::chemical::AtomIndices::const_iterator
				hnum = ngb.Hpos_polar().begin(), hnume = ngb.Hpos_polar().end();
				hnum != hnume; ++hnum) {

				Size hatm(*hnum);
				AtomID don(hatm, nidx);

				if(is_hbond(don, acc, ps, database, hbond_set))
					acounts[rep]++;
			}
		}
	}
}


///
/// @brief prints the contents of a storage of hydrogen bond counts
///
/// @param[in] counts the storage of hydrogen bond counts
/// @param[in] ps pose from which hydrogen bond counts were computed
///
/// @details The output can be seen as a sequnce of N blocks, where block i
/// 	reports the counts for the ith residue having at least one atom
/// 	contributing to counts (i=1,...,N, where N is the number of such
/// 	residues).
///   Within block i, line j reports the identifier and hydrogen bond count
///   for the residue's jth atom that contributes to counts (i=1,...,N;
/// 	j=1,...,M(i), where	M(i) is the number of such atoms).
///
void print_hb_counts(std::map<AtomID, Size>& counts, Pose const& ps) {

	for(Size i=1; i<=ps.total_residue(); ++i) {

		Residue const& r = ps.residue(i);

		for(Size j=1; j<=r.natoms(); ++j) {
			AtomID id(j, i);
			if(counts.find(id) != counts.end())
				TR << i << '(' << ps.pdb_info()->chain(i) << ps.pdb_info()->number(i) <<
					')' << r.atom_name(j) << ':' << counts[id] << std::endl;
		}
	}
}


///
/// @brief given two atoms, returns true if they have the same number of hydrogen
/// 	bonds; returns false otherwise.
///
/// @param[in] p1: AtomID and number of hydrogen bonds for the 1st atom
/// @param[in] p2: AtomID and number of hydrogen bonds for the 2nd atom
///
bool same_count(std::pair<AtomID, Size> const& p1,
	std::pair<AtomID, Size> const& p2) {

	return p1.second == p2.second;
}


/// definition of sameness between equivalent atoms
typedef bool (*SAME_ST_PTR)(std::pair<AtomID, Size> const& p1,
	std::pair<AtomID, Size> const& p2);

SAME_ST_PTR same_state = same_count;


///
/// @brief given the hydrogen-bond-count repositories of two equivalent sets of
/// 	atoms, prints the number of equivalent pairs of atoms that have the same
/// 	hydrogen-bond state.
///
/// @param[in] m1: hydrogen-bond-count repository for the 1st set of atoms.
/// @param[in] m2: hydrogen-bond-count repository for the 2nd set of atoms.
///
/// @details: it is assumed that m1 and m2 represent equivalent sets of atoms,
/// 	namely, for each AtomID x in m1 there exists an AtomID x in m2, and
/// 	vice versa.
///
/// @details: in particular, the function returns the number of equivalent
/// 	pairs of atoms (id1, id2) such that same_state((id1,c1), (id2,c2)) is
/// 	true,	where c1 and c2 are the counts associated with id1 and id2,
/// 	respectively.
///
Size num_same_state(std::map<AtomID, Size> m1, std::map<AtomID, Size> m2) {

	Size num_same = 0;
	for(std::map<AtomID, Size>::const_iterator it = m1.begin(), ite = m1.end();
		it != ite; ++it) {

		AtomID id = it->first;
		Size count = it->second;
		std::pair<AtomID, Size> p1(id, count);

		assert(m2.find(id) != m2.end());

		std::pair<AtomID, Size> p2(id, m2[id]);
		if(same_state(p1, p2))
			num_same++;
	}

	return num_same;
}


///
/// @brief copies a hydrogen-bond-count storage by keeping only contributions
/// 	from backbone	atoms.
///
/// @param[in] ps pose from which the storage was built
/// @param[in] src the storage (source)
/// @param[out] pruned space for the new, pruned storage (initially empty)
///
void prune_nobb(Pose const& ps,	std::map<AtomID, Size> const& src,
	std::map<AtomID, Size> &pruned) {

	for(std::map<AtomID, Size>::const_iterator it = src.begin(),
		ite = src.end(); it != ite; ++it) {

		AtomID id = it->first;
		Residue r = ps.residue(id.rsd());
		if(r.atom_is_backbone(id.atomno()))
			pruned[id] = it->second;
	}
}


OPT_KEY( String, cnl_resfile )
OPT_KEY( String, mut )
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                   MAIN                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char * argv [] )
{
	NEW_OPT( cnl_resfile, "set of residues forming the constellation", "cnl_resfile.txt" );
	NEW_OPT( mut, "PDB file containing the mutant structure", "" );

	devel::init(argc, argv);

	// create wild-type pose
	Pose ps;
	std::string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_pdb( ps, input_pdb_name );

	// create mutant pose
	Pose mut_ps;
	std::string const mut_pdb_name =
		basic::options::option[basic::options::OptionKeys::mut];
	core::import_pose::pose_from_pdb(mut_ps, mut_pdb_name);

	// score poses
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
	(*scorefxn)(ps);
	(*scorefxn)(mut_ps);

	// load constellation: cnl[i] contains the pose index of the ith residue in
	// the constellation input file (i=1,...,N, where N is the number of residues
	// in the file)
	std::string const cnlf = basic::options::option[cnl_resfile];
	std::ifstream cnlfs(cnlf.c_str());
	vector1<Size> cnl;
	int ri;
	char rc;
	while(cnlfs >> ri >> rc) {
		if(rc == '_')
			rc = ' ';
		cnl.push_back(get_pose_resnum(ri, rc, ps));
	}

	// find neighbors of the constellation in wild-type pose
	protocols::neighbor::Neighborhood n(cnl, ps, protocols::neighbor::in_ngbat_sphere);
	vector1<Size> ngb = n.get();
	TR << "neighbors of the constellation: " << std::endl;
	print_res_ids(ngb, ps, TR);
	TR << std::endl;

	// store hydrogen bond counts for wild-type neighbors
	HBondSet hb_set;
	hb_set.clear();
	HBondDatabase const& hb_db(
		*HBondDatabase::get_database(hb_set.hbond_options().params_database_tag()));

	std::map<AtomID, Size> dcounts;
	std::map<AtomID, Size> acounts;

	for(Size i=1; i<=ngb.size(); ++i)
		store_hb_counts(ngb[i], ps, hb_db, hb_set, dcounts, acounts);

	TR << "wild-type neighbor hbond counts: " << std::endl;
	print_hb_counts(dcounts, ps);
	TR << std::endl;
	print_hb_counts(acounts, ps);
	TR << std::endl;

	// store hydrogen bond counts for mutant neighbors
	std::map<AtomID, Size> mut_dcounts;
	std::map<AtomID, Size> mut_acounts;

	for(Size i=1; i<=ngb.size(); ++i)
		store_hb_counts(ngb[i], mut_ps, hb_db, hb_set, mut_dcounts, mut_acounts);

	TR << "mutant neighbor hbond counts: " << std::endl;
	print_hb_counts(mut_dcounts, mut_ps);
	TR << std::endl;
	print_hb_counts(mut_acounts, mut_ps);
	TR << std::endl;

	// store hydrogen bond counts for wild-type constellation backbone
	std::map<AtomID, Size> cnl_dcounts;
	std::map<AtomID, Size> cnl_acounts;
	std::map<AtomID, Size> bbcnl_dcounts;
	std::map<AtomID, Size> bbcnl_acounts;

	for(Size i=1; i<=cnl.size(); ++i)
		store_hb_counts(cnl[i], ps, hb_db, hb_set, cnl_dcounts, cnl_acounts);
	prune_nobb(ps, cnl_dcounts, bbcnl_dcounts);
	prune_nobb(ps, cnl_acounts, bbcnl_acounts);

	TR << "wild-type constellation backbone hbond counts: " << std::endl;
	print_hb_counts(bbcnl_dcounts, ps);
	TR << std::endl;
	print_hb_counts(bbcnl_acounts, ps);
	TR << std::endl;

	// store hydrogen bond counts for mutant constellation backbone
	std::map<AtomID, Size> cnl_mut_dcounts;
	std::map<AtomID, Size> cnl_mut_acounts;
	std::map<AtomID, Size> bbcnl_mut_dcounts;
	std::map<AtomID, Size> bbcnl_mut_acounts;

	for(Size i=1; i<=cnl.size(); ++i)
		store_hb_counts(cnl[i], mut_ps, hb_db, hb_set, cnl_mut_dcounts, cnl_mut_acounts);
	prune_nobb(mut_ps, cnl_mut_dcounts, bbcnl_mut_dcounts);
	prune_nobb(mut_ps, cnl_mut_acounts, bbcnl_mut_acounts);

	TR << "mutant constellation backbone hbond counts: " << std::endl;
	print_hb_counts(bbcnl_mut_dcounts, mut_ps);
	TR << std::endl;
	print_hb_counts(bbcnl_mut_acounts, mut_ps);
	TR << std::endl;

	Size nss =
		num_same_state(dcounts, mut_dcounts) +
		num_same_state(acounts, mut_acounts) +
		num_same_state(bbcnl_dcounts, bbcnl_mut_dcounts) +
		num_same_state(bbcnl_acounts, bbcnl_mut_acounts);

	Size tot =
		dcounts.size() +
		acounts.size() +
		bbcnl_dcounts.size() +
		bbcnl_acounts.size();

	TR << "fraction of same-state polar atoms: " << (nss/float(tot)) << std::endl;
}
