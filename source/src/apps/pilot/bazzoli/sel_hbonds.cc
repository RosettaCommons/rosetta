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

/// @brief Outputs selected hydrogen bonds from a pose.
///
/// @param[in] -s <PDB>, where <PDB> is the path to the PDB file containing the
/// 	pose.
///
/// @param[in] -bond_fil <BOFIL>, where <BOFIL> is the path to a file specifying
/// 	the hbonds to be searched for in the pose. The file has the following
/// 	format:
///
/// 	C1 I1 A1, c1 i1 a1\n
/// 	...
/// 	CM IM AM, cM iM aM\n
///
/// 	Here, M is the number of hydrogen bonds to be searched for. Ci, Ii, and Ai
/// 	indicate, respectively, the chain identifier, residue index, and name of
/// 	the DONOR atom in the ith  hbond to be searched for (i=1,...,M); ci,	ii,
/// 	and ai indicate, respectively, the chain identifier, residue index and
/// 	name of the ACCEPTOR atom in the ith hbond to be searched for (i=1,...,M).
///
/// @param[in] -atom_fil <ATFIL>, where <ATFIL> is the path to a file specifying
/// 	the atoms whose hbonds are to be searched for. The file has the following
/// 	format:
///
/// 	C1 I1 A1\n
/// 	...
/// 	CN IN AN\n
///
/// 	Here Ci, Ii, and Ai indicate, respectively, the chain identifier, residue
/// 	index, and name of the ith atom whose hbonds are to be searched for
/// 	(i=1,...,N, where N is the number of such atoms).
///
/// @details If both the -bond_fil and -atom_fil options are specified,	only
/// 	the former will be applied.
///
/// @details In the case of option -bond_fil, the ith output line indicates the
/// 	ith hbond from the input list (i=0,...,N-1, where N is the number of such
/// 	hbonds).
///
/// @details In the case of option -atom_fil, the output can be seen as a sequence
/// 	of S blocks, where the ith block lists the hbond identifiers of the ith
/// 	input atom that forms any hydrogen bonds (i=1,...,S).
///
/// @author Andrea Bazzoli (bazzoli@ku.edu)


#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/types.hh>
#include <devel/init.hh>
#include <numeric/xyzVector.hh>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>

static basic::Tracer TR( "apps.pilot.sel_hbonds.main" );

using core::Size;
using core::Real;
using core::id::AtomID;
using core::pose::Pose;
using std::string;
using core::conformation::Residue;
using numeric::xyzVector;
using std::setw;
using std::setfill;
using namespace core::scoring::hbonds;


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
/// @brief Tells whether a hydrogen atom may be an hbond donor.
///
/// @param[in] aid atom identifier in its pose.
/// @param[in] ps the pose.
///
/// @return true if the atom is a polar hydrogen bond; false otherwise.
///
bool is_donor_h(AtomID const& aid, Pose const& ps) {

	Residue const& rsd = ps.residue(aid.rsd());
	Size const HATM = aid.atomno();

	for(core::chemical::AtomIndices::const_iterator
		hnum = rsd.Hpos_polar().begin(), hnume = rsd.Hpos_polar().end();
		hnum != hnume; ++hnum) {

		Size const hatm(*hnum);
		if(hatm == HATM)
			return true;
	}

	return false;
}


///
/// @brief Tells whether an atom may be an hbond acceptor.
///
/// @param[in] aid atom identifier in its pose.
/// @param[in] ps the pose.
///
/// @return true if the atom may be an acceptor atom; false otherwise.
///
bool is_acceptor(AtomID const& aid, Pose const& ps) {

	Residue const& rsd = ps.residue(aid.rsd());
	Size const AATM = aid.atomno();

	for(core::chemical::AtomIndices::const_iterator
		anum = rsd.accpt_pos().begin(), anume = rsd.accpt_pos().end();
		anum != anume; ++anum) {

		Size const aatm(*anum);
		if(aatm == AATM)
			return true;
	}

	return false;
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
/// @brief Prints an hbond identifier
///
/// @param[in] ndon: name of the donor atom.
/// @param[in] idon: PDB index of the donor residue.
/// @param[in] cdon: PDB chain identifier of the donor residue.
/// @param[in] nacc: name of the acceptor atom.
/// @param[in] iacc: PDB index of the acceptor residue.
/// @param[in] cacc: PDB chain identifier of the acceptor residue.
///
void print_hbond(string const& ndon, int const idon, char const cdon,
	string const& nacc, int const iacc, char const cacc) {

	// strip blanks off atom names
	std::istringstream ndss(ndon);
	string ndon_;
	ndss >> ndon_;

	std::istringstream nass(nacc);
	string nacc_;
	nass >> nacc_;

	// print
	TR	<< setw(4) << setfill(' ') << ndon_ <<
		'(' << setw(4) << setfill('0') << idon << cdon << ')' <<
		" --- "
			<< setw(4) << setfill(' ') << nacc_ <<
		'(' << setw(4) << setfill('0') << iacc << cacc << ')' <<
		std::endl;
}


OPT_KEY( String, atom_fil )
OPT_KEY( String, bond_fil )
using basic::options::option;
using basic::options::OptionKeys::atom_fil;
using basic::options::OptionKeys::bond_fil;

int main( int argc, char * argv [] )
{
try {

	NEW_OPT( atom_fil, "selected hbonded atoms", "" );
	NEW_OPT( bond_fil, "selected hydrogen bonds", "" );

	devel::init(argc, argv);

	// create pose from PDB file
	Pose ps;
	string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_pdb( ps, input_pdb_name );

	// score pose
	core::scoring::ScoreFunctionOP scorefxn(
		core::scoring::get_score_function());
	(*scorefxn)(ps);

	// open input hbond file
	string const bf = option[bond_fil];
	string const af = option[atom_fil];
	string hbf;
	if(bf != "")
		hbf = bf;
	else if(af != "")
		hbf = af;
	else {
		TR << "Need option -atom_fil or option -bond_fil" << std::endl;
		return 0;
	}

	std::ifstream hbfs(hbf.c_str());
	if(!hbfs) {
		TR << "Can't find " << hbf << std::endl;
		return 0;
	}

	// initialize hbond auxiliary data structures
	HBondSet hb_set;
	hb_set.clear();
	HBondDatabase const & hb_db(
		*HBondDatabase::get_database(hb_set.hbond_options().params_database_tag()));
	core::scoring::EnergyGraph const& energy_graph(ps.energies().energy_graph());

	int const LINE_TERM_UB = 10;

	// option -bond_fil: read hbond identifiers and check hbond presence
	if(bf != "") {

		char cdon, cacc;
		int idon, iacc;
		string ndon, nacc;

		while(hbfs.get(cdon)) {

			// read hbond identifier
			hbfs >> idon >> ndon;
			ndon = ndon.substr(0, ndon.length() - 1);

			hbfs.ignore();
			hbfs.get(cacc);
			hbfs >> iacc >> nacc;

			hbfs.ignore(LINE_TERM_UB, '\n');

			// check if hbond is present
			Size p_idon = get_pose_resnum(idon, cdon, ps);
			Size p_iacc = get_pose_resnum(iacc, cacc, ps);

			for(core::graph::Graph::EdgeListConstIter
				nit = energy_graph.get_node(p_idon)->const_edge_list_begin(),
				nite = energy_graph.get_node(p_idon)->const_edge_list_end();
				nit != nite; ++nit ) {

				Size p_inb((*nit)->get_other_ind(p_idon));
				if(p_inb == p_iacc) {

					Size ai_don = ps.residue(p_idon).atom_index(ndon);
					Size ai_acc = ps.residue(p_iacc).atom_index(nacc);

					AtomID id_don(ai_don, p_idon);
					AtomID id_acc(ai_acc, p_iacc);

					if(is_donor_h(id_don, ps))
						if(is_acceptor(id_acc, ps))
							if(is_hbond(id_don, id_acc, ps, hb_db, hb_set))
								print_hbond(ndon, idon, cdon, nacc, iacc, cacc);
				}
			}
		}
	}

	// option -atom_fil: read atom identifiers and print their hbonds
	else {

		char c;
		int i;
		string n;

		while(hbfs.get(c)) {

			hbfs >> i >> n;
			hbfs.ignore(LINE_TERM_UB, '\n');

			Size pi = get_pose_resnum(i, c, ps);
			Size ai = ps.residue(pi).atom_index(n);
			AtomID aid(ai, pi);

			bool donor;
			if(!(donor = is_donor_h(aid, ps)))
				if(!is_acceptor(aid, ps))
					continue;

			// the ith output block lists the hbond identifiers of the ith neighbor
			// residue that forms hbonds with the atom (i=1,...,P, where P is the
			// number of such neighbor residues).
			for(core::graph::Graph::EdgeListConstIter
				nit = energy_graph.get_node(pi)->const_edge_list_begin(),
				nite = energy_graph.get_node(pi)->const_edge_list_end();
				nit != nite; ++nit ) {

				Size p_inb((*nit)->get_other_ind(pi));
				Residue const& rnb = ps.residue(p_inb);

				if(donor) {
					// the atom is a donor: the ith output line describes the hbond with
					// the neighbor residue's ith acceptor atom that does form an hbond
					// (i=1,...,Q, where Q is the number of such acceptor atoms in the
					// neighbor residue).
					for(core::chemical::AtomIndices::const_iterator
						anum = rnb.accpt_pos().begin(), anume = rnb.accpt_pos().end();
						anum != anume; ++anum) {

						Size aatm(*anum);
						AtomID nid(aatm, p_inb);

						if(is_hbond(aid, nid, ps, hb_db, hb_set))
							print_hbond(n, i, c, rnb.atom_name(aatm),
								ps.pdb_info()->number(p_inb), ps.pdb_info()->chain(p_inb));
					}
				}
				else
					// the atom is an acceptor: the ith output line describes the hbond
					// with the neighbor residue's ith donor atom that does form an hbond
					// (i=1,...,R, where R is the number of such donor atoms in the
					// neighbor residue).
					for(core::chemical::AtomIndices::const_iterator
						hnum = rnb.Hpos_polar().begin(), hnume = rnb.Hpos_polar().end();
						hnum != hnume; ++hnum) {

						Size hatm(*hnum);
						AtomID nid(hatm, p_inb);

						if(is_hbond(nid, aid, ps, hb_db, hb_set))
							print_hbond(rnb.atom_name(hatm), ps.pdb_info()->number(p_inb),
								ps.pdb_info()->chain(p_inb), n, i, c);
					}
			}
		}
	}

	return 0;

} // try
catch ( utility::excn::EXCN_Base const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}
} // main
