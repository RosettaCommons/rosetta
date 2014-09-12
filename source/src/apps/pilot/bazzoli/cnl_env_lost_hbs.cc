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
/// @brief Lists the atoms that form hydrogen bonds with a constellation of
/// 	residues--in a wild-type structure--but not with the ligand replacing the
///   constellation--in a mutant, holo structure.
///
/// @param[in] -s <PDBWT>, where <PDBWT> is the path to the PDB file containing
/// 	the wild-type structure.
///
/// @param[in] -mut <PDBMUT>, where <PDBMUT> is the path to the PDB file
/// 	containing the mutant structure.
///
/// @param[in] -extra_res_fa <PARAMS>, where <PARAMS> is the path to the params
/// 	file containing the Rosetta description of the ligand in the mutant
/// 	structure.
///
/// @param[in] -cnl_resfile <CNL>, where <CNL> is the path to a file enumerating
/// 	the residue-specific contributions to the constellation. The file has the
///   following	format:
/// 	I1 C1 NEW1\n
/// 	...
/// 	IN CN NEWN\n ,
/// 	where Ii, Ci, and NEWi are the residue index, chain ID, and new amino acid
/// 	type, respectively, of the ith residue forming the constellation
/// 	(i=1,...,N).
///
/// @details The input constellation file requires that a chain identifier equal
/// 	to ' ' in the PDB file be represented by '_'.
///
/// @details It is assumed that the ligand is the last residue in the mutant
/// 	pose.
///
/// @details Some atoms are considered to be equivalent regarding hydrogen bonding
/// 	and are therefore	represented	by a single atom name. See function
/// 	hb_rep_atom() for details.
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
using core::chemical::oneletter_code_from_aa;

static basic::Tracer TR( "apps.pilot.cnl_env_lost_hbs" );


// data structure holding the set of h-bond partners for each atom
typedef std::map<AtomID, vector1<AtomID> > HB_Partners;


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
/// 	a target atom for hydrogen bonding.
///
/// @param[in] tgt index in the residue of the target atom
/// @param[in] res residue to which the target atom belongs
///
/// @details: the representative atom is used to make some atoms equivalent with
/// 	respect to hydrogen bonding. Currently:
/// 	- all hydrogen atoms bound to the same heavy atom are equivalent.
///   	- Arg's 1HH1, 1HH2, 2HH1, and 2HH2 are all equivalent.
///   	- Asp's OD1 and OD2 are equivalent.
///   	- Glu's OE1 and OE2 are equivalent.
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
/// @brief prints an atom's PDB identifier
///
/// param[in] aid atom identifier in the pose
/// param[in] ps the pose
/// param[out] os output stream
///
void print_atom_pdbid(AtomID const& aid, Pose const& ps, std::ostream& os) {

	Size const ridx = aid.rsd();
	Size const aidx = aid.atomno();
	Residue const& res = ps.residue(ridx);
	os <<	oneletter_code_from_aa(res.aa()) <<
		ps.pdb_info()->number(ridx) <<
		ps.pdb_info()->chain(ridx) << '(' << ridx << ')' <<
		res.atom_name(aidx);
}


///
/// @brief prints the contents of an HB_Partners instance
///
/// @param[in] ptns the HB_Partners instance
/// @param[in] ps the relevant pose
/// @param[out] os output stream
///
/// @details The output can be seen as a sequence of N blocks, where N is the
/// 	number of keys (atoms whose partners are stored) in ptns. The ith block
/// 	lists the partners of the ith key in ptns (i=0,...,N-1). Within block i,
/// 	the jth	line contains the identifier of ptns[i][j] (i=0,...,N-1;
/// 	j=1,...,M(i), where M(i) is the number of partners of ptns[i]).
///
void print_hb_partners(HB_Partners const& ptns, Pose const& ps, std::ostream& os) {

	for(HB_Partners::const_iterator i=ptns.begin(),	END=ptns.end();
		i != END; ++i) {

		os << "H-bond partners of ";
		print_atom_pdbid(i->first, ps, os);
		os << ':' << std::endl;

		vector1<AtomID> const& pv = i->second;
		Size const NP = pv.size();
		for(Size j=1; j<=NP; ++j) {
			print_atom_pdbid(pv[j], ps, os);
			os << std::endl;
		}
	}
}


///
/// @brief records a donor atom as an hbond partner for each atom accepting
/// 	hbonds from it
///
/// @param[in] don the donor atom
/// @param[in] ps pose to which the atom and its partners belong
/// @param[in] database hbond database
/// @param[in] hbond_set dummy hbond set (will not be filled in)
/// @param[out] hbptns stores the donor atom as an hbond partner
///
/// @details For each atom x accepting an hbond from 'don', 'don' is added as a
/// 	a partner of R(x) in 'hbptns' (i.e., to vector hbptns[R(x)]), where R(x)
/// 	is the hbond representative atom of x.
///
void don_store_hbs(AtomID const& don, Pose const& ps,
	HBondDatabase const& database, HBondSet const& hbond_set, HB_Partners& hbptns) {

	Size const dri = don.rsd();
	core::scoring::EnergyGraph const& energy_graph(ps.energies().energy_graph());

	core::graph::Graph::EdgeListConstIter NITB =
		energy_graph.get_node(dri)->const_edge_list_begin();
	core::graph::Graph::EdgeListConstIter NITE =
		energy_graph.get_node(dri)->const_edge_list_end();

	for(core::graph::Graph::EdgeListConstIter nit = NITB; nit != NITE; ++nit) {

		Size nri((*nit)->get_other_ind(dri));
		Residue const& ngb = ps.residue(nri);

		for(core::chemical::AtomIndices::const_iterator
			anum = ngb.accpt_pos().begin(), anume = ngb.accpt_pos().end();
			anum != anume; ++anum) {

			Size aai(*anum);
			AtomID acc(aai, nri);

			if(is_hbond(don, acc, ps, database, hbond_set)) {
				AtomID repacc(rep_hb_atom(aai, ngb), nri);
				hbptns[repacc].push_back(don);
			}
		}
	}
}


///
/// @brief records an acceptor atom as an hbond partner for each atom donating
/// 	hbonds to it.
///
/// @param[in] acc the acceptor atom
/// @param[in] ps pose to which the atom and its partners belong
/// @param[in] database hbond database
/// @param[in] hbond_set dummy hbond set (will not be filled in)
/// @param[out] hbptns stores the acceptor atom as an hbond partner
///
/// @details For each atom x donating an hbond to 'acc', 'acc' is added as a
/// 	a partner of R(x) in 'hbptns' (i.e., to vector hbptns[R(x)]), where R(x)
/// 	is the hbond representative atom of x.
///
void acc_store_hbs(AtomID const& acc, Pose const& ps,
	HBondDatabase const& database, HBondSet const& hbond_set, HB_Partners& hbptns) {

	Size const ari = acc.rsd();
	core::scoring::EnergyGraph const& energy_graph(ps.energies().energy_graph());

	core::graph::Graph::EdgeListConstIter NITB =
		energy_graph.get_node(ari)->const_edge_list_begin();
	core::graph::Graph::EdgeListConstIter NITE =
		energy_graph.get_node(ari)->const_edge_list_end();

	for(core::graph::Graph::EdgeListConstIter nit = NITB; nit != NITE; ++nit) {

		Size nri((*nit)->get_other_ind(ari));
		Residue const& ngb = ps.residue(nri);

		for(core::chemical::AtomIndices::const_iterator
			dnum = ngb.Hpos_polar().begin(), dnume = ngb.Hpos_polar().end();
			dnum != dnume; ++dnum) {

			Size dai(*dnum);
			AtomID don(dai, nri);

			if(is_hbond(don, acc, ps, database, hbond_set)) {
				AtomID repdon(rep_hb_atom(dai, ngb), nri);
				hbptns[repdon].push_back(acc);
			}
		}
	}
}


///
/// @brief stores a residue's constellation atoms as hbond partners of the
/// 	atoms, if any, with which they form hydrogen bonds.
///
/// @param[in] ridx index of the residue in its pose
/// @param[in] newaa amino acid type to which the residue is mutated (to yield
/// 	the constellation)
/// @param[in] ps the pose
/// @param[in] database hbond database
/// @param[in] hbond_set dummy hbond set (will not be filled in)
/// @param[out] hbptns repository of hbond partners
///
/// @details it is assumed that, for all constellations but Thr->Ser, all polar
/// 	side-chain atoms are part of the constellation (this holds for the present
/// 	ensemble of donors and acceptors in Rosetta).
///
void cnlres_store_hbs(const Size ridx, const char newaa, Pose const& ps,
	HBondDatabase const& database, HBondSet const& hbond_set,
	HB_Partners& hbptns) {

	Residue const& res = ps.residue(ridx);

	if((ps.aa(ridx) == core::chemical::aa_thr) && (newaa == 'S'))
		return;

	// donors
	for(core::chemical::AtomIndices::const_iterator
		dnum = res.Hpos_polar().begin(), dnume = res.Hpos_polar().end();
		dnum != dnume; ++dnum) {

		Size const didx = *dnum;
		if(!(res.atom_is_backbone(didx))) {
			AtomID don(didx, ridx);
			don_store_hbs(don, ps, database, hbond_set, hbptns);
		}
	}

	// acceptors
	for(core::chemical::AtomIndices::const_iterator
		anum = res.accpt_pos().begin(), anume = res.accpt_pos().end();
		anum != anume; ++anum) {

		Size const aidx = *anum;
		if(!(res.atom_is_backbone(aidx))) {
			AtomID acc(aidx, ridx);
			acc_store_hbs(acc, ps, database, hbond_set, hbptns);
		}
	}
}

/// @brief data structure for a residue contributing to a constellation
struct CnlRes {

	Size ridx; // residue index in its pose
	char newaa; // amino acid type to which the residue is mutated

	CnlRes(Size i, char a) : ridx(i), newaa(a) {}
};


///
/// @brief stores a constellation's atoms as hbond partners of the
/// 	atoms, if any, with which they form hydrogen bonds.
///
/// @param[in] cnl the constellation
/// @param[in] ps pose to which the constellation belongs
/// @param[in] database hbond database
/// @param[in] hbond_set dummy hbond set (will not be filled in)
/// @param[out] hbptns repository of hbond partners
///
void cnl_store_hbs(vector1<CnlRes> const& cnl, Pose const& ps,
  HBondDatabase const& database, HBondSet const& hbond_set,
  HB_Partners& hbptns) {

	Size const N = cnl.size();
	for(Size i=1; i<=N; ++i)
		cnlres_store_hbs(cnl[i].ridx, cnl[i].newaa, ps, database, hbond_set,
			hbptns);
}


///
/// @brief prints a constellation's description
///
/// @param[in] cnl residue-specific contributions to the constellation
/// @param[in] ps pose to which the constellation belongs
/// @param[out] os output stream
///
/// @details the ith output line contains the description of the ith
/// 	residue contribution in cnl (i=1,...,N, where N is the size of cnl).
///
void cnl_print(vector1<CnlRes> const& cnl, Pose const& ps, std::ostream& os) {

	Size const N = cnl.size();
	for(Size i=1; i<=N; ++i) {
		Size const ridx = cnl[i].ridx;
		char const newaa = cnl[i].newaa;
		os << oneletter_code_from_aa(ps.aa(ridx)) << ps.pdb_info()->number(ridx) <<
			ps.pdb_info()->chain(ridx) <<	'(' << ridx << ')' << " --> " << newaa <<
			std::endl;
	}
}


///
/// @brief returns true if atom 'polat' belongs to constellation 'cnl'; returns
/// 	false otherwise.
///
/// @param[in] polat an hbond-representative polar atom
/// @param[in] cnl the constellation
/// @param[in] ps pose to which the atom and the constellation belong
///
/// @details it is assumed that, for all single-residue constellations but
/// 	Thr->Ser, all polar side-chain atoms are part of the constellation (this
/// 	holds for the present ensemble of donors and acceptors in Rosetta).
///
bool in_cnl(AtomID const& polat, vector1<CnlRes> const& cnl, Pose const& ps) {

	using core::chemical::aa_thr;

	Size const ridx = polat.rsd();
	Size const aidx = polat.atomno();
	Residue const& res = ps.residue(ridx);

	Size N = cnl.size();
	for(Size i=1; i<=N; ++i)
		if(ridx == cnl[i].ridx)
			if(!res.atom_is_backbone(aidx))
				if((res.aa() != aa_thr) || cnl[i].newaa != 'S')
					return true;

	return false;
}


//
// @brief deletes from a repository of hbond partners those atoms (keys) that
// 	belong to a constellation
//
// @param[out] hbptns the repository of hbond partners
// @param[in] cnl the constellation
// @param[in] ps pose from which the repository and the constellation were
// 	built.
//
void prune_away_cnl(HB_Partners& hbptns, vector1<CnlRes> const& cnl,
	Pose const& ps) {

	bool found_cnl;
	do {

		found_cnl = false;

		for(HB_Partners::iterator i=hbptns.begin(),	END=hbptns.end();
			i != END; ++i)
			if(in_cnl(i->first, cnl, ps)) {
				hbptns.erase(i);
				found_cnl = true;
				break;
			}
	}
	while(found_cnl == true);
}


///
/// @brief Detects the hbonds between an hbond-representative atom and a partner
/// 	residue
///
/// @param[in] rep the hbond-representative atom
/// @param[in] ptnidx index of the partner residue in its pose
/// @param[in] ps pose to which both atom 'rep' and residue 'ptnidx' belong
/// @param[in] database hbond database
/// @param[in] hbond_set dummy hbond set (will not be filled in)
/// @param[out] os output stream where hbonds are listed
///
/// @return number of hydrogen bonds between atom 'rep' (i.e., all atoms it
/// 	represents) and residue 'ptnidx'.
///
/// @details if 'rep' represents only acceptors or only donors,	the output can
/// 	be seen as a sequence of N blocks, where N is the number of atoms
/// 	represented by 'rep' that form at least one hydrogen bond with residue
/// 	'ptnidx'. The ith block lists the hydrogen bonds made by the ith such
/// 	represented atom, according to the order specified by Residue::accpt_pos()
/// 	(acceptors) or Residue::Hpos_polar() (donors) (i=0,...,N-1). Within block
/// 	i, the jth line denotes the jth polar atom in residue 'ptnidx' to form a
/// 	hydrogen bond with the represented atom	in question for block i
/// 	(i=0,...,N-1; j=0,...,M(i)-1, where M(i) is the number of such polar
/// 	atoms, and acceptors and donors in residue 'ptnidx' are also ordered
/// 	according to Residue::accpt_pos() and Residue::Hpos_polar(),
/// 	respectively).
///
/// @details if 'rep' represents both donors and acceptors, 2 instances of the
/// 	above output are produced, the first for the represented acceptors and the
/// 	second, immediately following, for the represented donors.
///
Size hbonds_to_ptnres(AtomID const& rep, Size ptnidx, Pose const& ps,
	HBondDatabase const& database, HBondSet const& hbond_set, std::ostream& os) {

	Size const rridx = rep.rsd();
	Size const raidx = rep.atomno();
	Residue const& rres = ps.residue(rridx);

	Residue const& pres = ps.residue(ptnidx);

	Size nbonds = 0;

	// atom 'rep' represents acceptors
	if(rres.heavyatom_is_an_acceptor(raidx)) {

		for(core::chemical::AtomIndices::const_iterator
			anum = rres.accpt_pos().begin(), anume = rres.accpt_pos().end();
			anum != anume; ++anum) {

			Size aai(*anum);
			if(rep_hb_atom(aai, rres) == raidx) {

				AtomID acc(aai, rridx);

				for(core::chemical::AtomIndices::const_iterator
					dnum = pres.Hpos_polar().begin(), dnume = pres.Hpos_polar().end();
					dnum != dnume; ++dnum) {

					AtomID don(*dnum, ptnidx);

					if(is_hbond(don, acc, ps, database, hbond_set)) {
						print_atom_pdbid(don, ps, os);
						os << std::endl;
						nbonds++;
					}
				}
			}
		}
	}

	// atom 'rep' represents donors
	if(rres.heavyatom_has_polar_hydrogens(raidx)) {

		for(core::chemical::AtomIndices::const_iterator
			dnum = rres.Hpos_polar().begin(), dnume = rres.Hpos_polar().end();
			dnum != dnume; ++dnum) {

			Size dai(*dnum);
			if(rep_hb_atom(dai, rres) == raidx) {

				AtomID don(dai, rridx);

				for(core::chemical::AtomIndices::const_iterator
					anum = pres.accpt_pos().begin(), anume = pres.accpt_pos().end();
					anum != anume; ++anum) {

					AtomID acc(*anum, ptnidx);

					if(is_hbond(don, acc, ps, database, hbond_set)) {
						print_atom_pdbid(acc, ps, os);
						os << std::endl;
						nbonds++;
					}
				}
			}
		}
	}

	return nbonds;
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

try {

	NEW_OPT( cnl_resfile, "set of residues forming the constellation", "" );
	NEW_OPT( mut, "PDB file containing the mutant structure", "");

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
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	(*scorefxn)(ps);
	(*scorefxn)(mut_ps);

	// load constellation: cnl[i] describes the contribution to the constellation
	// by the ith residue in the constellation file (i=1,...,N, where N is the
	// number of residues in the constellation file)
	std::string const cnlf =
		basic::options::option[basic::options::OptionKeys::cnl_resfile];
	std::ifstream cnlfs(cnlf.c_str());
	vector1<CnlRes> cnl;
	int ri;
	char rc;
	char newaa;
	while(cnlfs >> ri >> rc >> newaa) {
		if(rc == '_')
			rc = ' ';
		cnl.push_back(CnlRes(get_pose_resnum(ri, rc, ps), newaa));
	}

	std::cout << std::endl;
	std::cout << "### computing H-bond partners for the following constellation:  ";
	std::cout << std::endl;
	cnl_print(cnl, ps, std::cout);

	// for each hbond-representative atom involved in a hydrogen bond with the
	// constellation, build the set of its constellation partners. envptns[x][i]
	// contains the ith partner to be found for atom x, where (1) constellation
	// residues are evaluated in the order specified by vector cnl, (2) for each
	// constellation residue, donors are evaluated before acceptors, (3)
	// a constellation residue's donors are evaluated in the order specified by
	// Residue::Hpos_polar(), and (4) a constellation residue's acceptors are
	// evaluated in the order specified by Residue::accpt_pos().
	HB_Partners env_hb_ptns;
	HBondSet hb_set;
	hb_set.clear();
	HBondDatabase const& hb_db(
 		*HBondDatabase::get_database(hb_set.hbond_options().params_database_tag()));

	cnl_store_hbs(cnl, ps, hb_db, hb_set, env_hb_ptns);
	std::cout << std::endl;
	prune_away_cnl(env_hb_ptns, cnl, ps);

	std::cout << "### H-bonds to the constellation in the WT pose: "
		<< std::endl;
	print_hb_partners(env_hb_ptns, ps, std::cout);

	// find hbonds between constellation partners and ligand in the mutant pose.
	// The output is a list of N blocks, where N is the size of env_hb_ptns. The
	// ith block lists the hbond partners in the mutant pose of the constellation
	// partner atom at env_hb_ptnsp[i] (i=0,...,N-1). For each block, partners
	// are listed according to function hbonds_to_ptnres()
	std::cout << std::endl << "### conserved H-bonds to the ligand in the mutant pose: "
		<< std::endl;

	using core::chemical::aa_thr;
	using core::chemical::aa_ser;
	const Size LIGIDX = mut_ps.total_residue();
	Size lost_hbs = 0;
	for(HB_Partners::const_iterator i=env_hb_ptns.begin(), END=env_hb_ptns.end();
		i!=END; ++i ) {

		AtomID envat = i->first;

		Size wt_ridx = envat.rsd();
		char pdb_cid = ps.pdb_info()->chain(wt_ridx);
		int pdb_ridx = ps.pdb_info()->number(wt_ridx);

		Size mut_ridx = get_pose_resnum(pdb_ridx, pdb_cid, mut_ps);
		envat.rsd() = mut_ridx;

		// in case of a T->S mutation, change OG1 (if any) to OG
		if(ps.aa(wt_ridx) == aa_thr && mut_ps.aa(mut_ridx) == aa_ser)
			if(!ps.residue(wt_ridx).atom_is_backbone(envat.atomno()))
				envat.atomno() = mut_ps.residue(mut_ridx).atom_index("OG");

		std::cout << "H-bond partners of ";
		print_atom_pdbid(envat, mut_ps, std::cout);
		std::cout << ':' << std::endl;

		Size hbs = hbonds_to_ptnres(envat, LIGIDX, mut_ps, hb_db, hb_set, std::cout);
		lost_hbs += (i->second.size() - hbs);
	}

	std::cout << std::endl << "### number of H-bonds lost by the environment: " <<
		lost_hbs << std::endl;

} // try
catch ( utility::excn::EXCN_Base const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
}

}
