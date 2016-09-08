// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @brief Prints the fraction of polar atoms in the ligand that are satisfied
///  by a hydrogen bond.
///
/// @param[in] -s <PDBFIL>, where <PDBFIL> is the path to the PDB file
///  containing the protein+ligand complex.
///
/// @param[in] -extra_res_fa <PARAMS>, where <PARAMS> is the path to the .params
///  file describing the ligand.
///
/// @details It is assumed that the ligand is the last residue of the pose.
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
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <numeric/xyzVector.hh>
#include <core/chemical/AA.hh>
#include <utility/vector1.hh>

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

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.lig_polar_sat" );


/// @brief returns the index in the residue of the atom designated to represent
///  a target atom for hydrogen bonding.
///
/// @param[in] tgt index in the residue of the target atom
/// @param[in] res residue to which the target atom belongs
///
/// @details: the representative atom is used to make some atoms equivalent with
///  respect to hydrogen bonding. Currently:
///  - all hydrogen atoms bound to the same heavy atom are equivalent.
///
Size rep_hb_atom(Size const tgt, Residue const& res) {

	if ( res.atom_is_hydrogen(tgt) ) {
		return res.atom_base(tgt);
	} else {
		return tgt;
	}
}


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
///  accept hydrogen bonds.
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

	if ( hatm_xyz.distance_squared(aatm_xyz) <= MAX_R2 ) {

		Real unweighted_energy(0.0);
		HBEvalTuple hbe_type(datm, don_rsd, aatm, acc_rsd);

		int const base(acc_rsd.atom_base(aatm));
		int const base2(acc_rsd.abase2(aatm));
		assert(base2 > 0 && base != base2);

		core::scoring::hbonds::hb_energy_deriv(database,
			hbond_set.hbond_options(), hbe_type, datm_xyz, hatm_xyz, aatm_xyz,
			acc_rsd.atom(base).xyz(), acc_rsd.atom(base2).xyz(),
			unweighted_energy, false, derivs);

		if ( unweighted_energy < MAX_HB_ENERGY ) {
			return true;
		} else {
			return false;
		}
	} else return false;
}


/// @brief computes and stores the number of hydrogen bonds for each
///  representative atom of a residue.
///
/// @param[in] ridx index of the residue in the pose
/// @param[in] ps the pose
/// @param[in] database hbond database.
/// @param[in] hbond_set dummy hbond set (will not be filled in).
/// @param[out] dcounts storage for the donor atoms' hydrogen bond counts
/// @param[out] acounts storage for the acceptor atoms' hydrogen bond counts
///
/// @details For each representative atom r of donor atoms,
///  dcounts[AtomID(r, ridx)] will contain the total number of hydrogen bonds
///  made by the donor atoms it represents.
///
/// @details For each representative atom s of acceptor atoms,
///  acounts[AtomID(s, ridx)] will contain the total number of hydrogen bonds
///  made by the acceptor atoms it represents.
///
/// @details It is assumed that neither dcounts nor acounts contain, when
///  the function begins, any item for residue ridx.
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
	for ( core::chemical::AtomIndices::const_iterator
			hnum = res.Hpos_polar().begin(), hnume = res.Hpos_polar().end();
			hnum != hnume; ++hnum ) {

		Size hatm(*hnum);
		AtomID don(hatm, ridx);

		Size ratm = rep_hb_atom(hatm, res);
		AtomID rep(ratm, ridx);
		if ( dcounts.find(rep) == dcounts.end() ) {
			dcounts[rep] = 0;
		}

		for ( core::graph::Graph::EdgeListConstIter nit = NITB; nit != NITE; ++nit ) {

			Size nidx((*nit)->get_other_ind(ridx));
			Residue const& ngb = ps.residue(nidx);

			for ( core::chemical::AtomIndices::const_iterator
					anum = ngb.accpt_pos().begin(), anume = ngb.accpt_pos().end();
					anum != anume; ++anum ) {

				Size aatm(*anum);
				AtomID acc(aatm, nidx);

				if ( is_hbond(don, acc, ps, database, hbond_set) ) {
					dcounts[rep]++;
				}
			}
		}
	}

	// collect hbonds for acceptors
	for ( core::chemical::AtomIndices::const_iterator
			anum = res.accpt_pos().begin(), anume = res.accpt_pos().end();
			anum != anume; ++anum ) {

		Size aatm(*anum);
		AtomID acc(aatm, ridx);

		Size ratm = rep_hb_atom(aatm, res);
		AtomID rep(ratm, ridx);
		if ( acounts.find(rep) == acounts.end() ) {
			acounts[rep] =  0;
		}

		for ( core::graph::Graph::EdgeListConstIter nit = NITB; nit != NITE; ++nit ) {

			Size nidx((*nit)->get_other_ind(ridx));
			Residue const& ngb = ps.residue(nidx);

			for ( core::chemical::AtomIndices::const_iterator
					hnum = ngb.Hpos_polar().begin(), hnume = ngb.Hpos_polar().end();
					hnum != hnume; ++hnum ) {

				Size hatm(*hnum);
				AtomID don(hatm, nidx);

				if ( is_hbond(don, acc, ps, database, hbond_set) ) {
					acounts[rep]++;
				}
			}
		}
	}
}


/// @brief prints the contents of a storage of hydrogen bond counts
///
/// @param[in] counts the storage of hydrogen bond counts
/// @param[in] ps pose from which hydrogen bond counts were computed
///
/// @details The output can be seen as a sequnce of N blocks, where block i
///  reports the counts for the ith residue having at least one atom
///  contributing to counts (i=1,...,N, where N is the number of such
///  residues).
///   Within block i, line j reports the identifier and hydrogen bond count
///   for the residue's jth atom that contributes to counts (i=1,...,N;
///  j=1,...,M(i), where M(i) is the number of such atoms).
///
void print_hb_counts(std::map<AtomID, Size>& counts, Pose const& ps) {

	for ( Size i=1; i<=ps.size(); ++i ) {

		Residue const& r = ps.residue(i);

		for ( Size j=1; j<=r.natoms(); ++j ) {
			AtomID id(j, i);
			if ( counts.find(id) != counts.end() ) {
				TR << i << '(' << ps.pdb_info()->chain(i) << ps.pdb_info()->number(i) <<
					')' << r.atom_name(j) << ':' << counts[id] << std::endl;
			}
		}
	}
}


/// @brief given a hydrogen-bond-count storage, returns the number of atoms
///  which have at least a hydrogen bond.
///
/// @param[in] counts the storage
///
Size nsat(std::map<AtomID, Size> const& counts) {

	Size ns = 0;
	for ( std::map<AtomID, Size>::const_iterator it = counts.begin(),
			ite = counts.end(); it != ite; ++it ) {
		if ( it->second ) {
			ns++;
		}
	}

	return ns;
}


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                   MAIN                                     //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);

		// create pose
		Pose ps;
		std::string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_file( ps, input_pdb_name , core::import_pose::PDB_file);

		// score pose
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();

		(*scorefxn)(ps);

		// store hydrogen bond counts for ligand
		HBondSet hb_set;
		hb_set.clear();
		HBondDatabase const& hb_db(
			*HBondDatabase::get_database(hb_set.hbond_options().params_database_tag()));

		std::map<AtomID, Size> dcounts;
		std::map<AtomID, Size> acounts;

		store_hb_counts(ps.size(), ps, hb_db, hb_set, dcounts, acounts);

		TR << "hbond counts: " << std::endl;
		print_hb_counts(dcounts, ps);
		TR << std::endl;
		print_hb_counts(acounts, ps);
		TR << std::endl;

		Size nsd = nsat(dcounts);
		Size nsa = nsat(acounts);
		Size ns = nsd + nsa;

		Size tot = dcounts.size() + acounts.size();

		Size nu = tot - ns;
		TR << "number of unsatisfied polar atoms: " << nu << std::endl;

		float frac = tot ? ns / float(tot) : 1.0F;
		TR << "fraction of satisfied polar atoms: " << frac << std::endl;

	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
}
