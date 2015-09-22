// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps.pilot.kevin.khxtal_water_bunsat.cc
/// @brief
/// @details
/// @author Kevin Houlihan

//#include <iomanip>
#include <algorithm>
#include <cctype>
#include <limits>
#include <boost/unordered/unordered_map.hpp>

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>
//#include <numeric/xyzVector.hh>
#include <core/import_pose/import_pose.hh>

//#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondSet.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>

#include <core/id/AtomID_Map.hh>

//#include <core/graph/Graph.hh>
//#include <core/scoring/EnergyGraph.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

#include <core/scoring/hbonds/hbonds.hh>

//#include <core/scoring/sasa.hh>

//#include <devel/buns/BuriedUnsatHbondFilter2.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>

#include <core/types.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/options/keys/RealOptionKey.hh>
//#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>

#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <core/pose/util.hh>

//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static THREAD_LOCAL basic::Tracer TR( "buns_stats" );

basic::options::RealOptionKey const AHdist_opkey("AHdist");
basic::options::RealOptionKey const AHD_opkey("AHD");

using namespace core;
using namespace core::pose::metrics;
using namespace devel::vardist_solaccess;
using core::Distance;

void read_FileVector(const basic::options::FileVectorOptionKey& fv_key, std::vector<std::string>& str_vec) {
	if(!basic::options::option[fv_key].user()) return;
	std::string list_file = basic::options::option[fv_key]()[1];
	utility::io::izstream list;
	list.open(list_file,std::_S_in);
	while(!list.eof()) {
		std::string line;
		getline(list,line);
		if(line.size() > 0) str_vec.push_back(line);
	}
	list.close();
}


class hbond_strat {
	public:
	hbond_strat() {};
	virtual ~hbond_strat() {}

	virtual bool evaluate(const pose::Pose& pose, scoring::hbonds::HBondCOP h) const = 0;
};


class hbond_energy_strat : public hbond_strat {
	public:
	hbond_energy_strat(Real energy_cutoff) : energy_cutoff_(energy_cutoff) {}

	bool evaluate(const pose::Pose&, scoring::hbonds::HBondCOP h) const {
		return ( h->energy() * h->weight() <= energy_cutoff_ );
	}
	private:
	Real energy_cutoff_;
};


class hbond_geom_strat : public hbond_strat {
	public:
	hbond_geom_strat (
		Real AHdist_geom_eval_threshold,
		Real AHD_geom_eval_threshold
	) :
		AHdist_geom_eval_threshold_(AHdist_geom_eval_threshold),
		AHD_geom_eval_threshold_(AHD_geom_eval_threshold)
	{}

	bool evaluate (const pose::Pose& pose, scoring::hbonds::HBondCOP h) const {
		Real AHdist = h->get_HAdist(pose);
		Real AHD = h->get_AHDangle(pose);

		return (AHdist <= AHdist_geom_eval_threshold_ &&
		       AHD >= AHD_geom_eval_threshold_);
	}

	private:
	Real AHdist_geom_eval_threshold_;
	Real AHD_geom_eval_threshold_;
};


void
vsasa_bunsats(
	pose::Pose& pose,
	const Real AHdist_threshold,
	const hbond_strat& hb_eval,
	VarSolDistSasaCalculator& vsasa_calc,
	const Real burial_cutoff
) {
	pose.update_residue_neighbors();
	scoring::ScoreFunctionOP scorefxn = scoring::get_score_function();
	scorefxn->score(pose);
	pose.energies();

	const Real max_nbr_radius = pose::pose_max_nbr_radius(pose);
	//constexpr const Real max_covalent_hbond_length = 1.5;
	// max_nbr_radius considers heavyatoms,
	// distance from covalent bond to hydrogen + hbond interaction range is required
	const Real max_covalent_hbond_length = 1.5;
	const Real neighbor_cutoff = 2 * max_nbr_radius + max_covalent_hbond_length + AHdist_threshold;

	conformation::PointGraphOP pg = new conformation::PointGraph;
	conformation::residue_point_graph_from_conformation(pose.conformation(), *pg);
	conformation::find_neighbors<conformation::PointGraphVertexData,
		conformation::PointGraphEdgeData>(pg, neighbor_cutoff);

	id::AtomID_Map<Real> atom_sasa = vsasa_calc.calculate(pose);

	core::Size nres = pose.total_residue();

	// find hbonds and intersect residues
	scoring::hbonds::HBondSet hbond_set;

	const scoring::hbonds::HBondDatabaseCOP hb_database = scoring::hbonds::HBondDatabase::get_database();
	const scoring::TenANeighborGraph& tenA_neighbor_graph(pose.energies().tenA_neighbor_graph());

	for (core::Size lowerResNum = 1; lowerResNum <= nres; ++lowerResNum ) {
		const conformation::Residue& lowerRes(pose.residue(lowerResNum));
		const core::Size nbl = tenA_neighbor_graph.get_node(lowerResNum)->
			num_neighbors_counting_self_static();
		for (conformation::PointGraph::UpperEdgeListConstIter
				ue  = pg->get_vertex(lowerResNum).const_upper_edge_list_begin(),
				ue_end = pg->get_vertex(lowerResNum).const_upper_edge_list_end();
				ue != ue_end; ++ue ) {
			core::Size upperResNum = ue->upper_vertex();
			const conformation::Residue& upperRes(pose.residue(upperResNum));
			const core::Size nbu = tenA_neighbor_graph.get_node(upperResNum)->
				num_neighbors_counting_self_static();
			scoring::hbonds::identify_hbonds_1way_AHdist(*hb_database,
				lowerRes, upperRes, nbl, nbu, AHdist_threshold, hbond_set);
			scoring::hbonds::identify_hbonds_1way_AHdist(*hb_database,
				upperRes, lowerRes, nbu, nbl, AHdist_threshold, hbond_set);
		}
	}

	boost::unordered_map<std::string, Size> res_counts;
	res_counts["ALA"] = 0;
	res_counts["CYS"] = 0;
	res_counts["ASP"] = 0;
	res_counts["GLU"] = 0;
	res_counts["PHE"] = 0;
	res_counts["GLY"] = 0;
	res_counts["HIS"] = 0;
	res_counts["ILE"] = 0;
	res_counts["LYS"] = 0;
	res_counts["LEU"] = 0;
	res_counts["MET"] = 0;
	res_counts["PRO"] = 0;
	res_counts["GLN"] = 0;
	res_counts["ARG"] = 0;
	res_counts["SER"] = 0;
	res_counts["THR"] = 0;
	res_counts["VAL"] = 0;
	res_counts["TRP"] = 0;
	res_counts["TYR"] = 0;

	core::Size buns = 0;
	const pose::PDBInfo& pdb_info = *(pose.pdb_info());
	for (core::Size resNum = 1; resNum <= nres; ++resNum) {
		const conformation::Residue& res = pose.residue(resNum);
		const chemical::ResidueType& res_type = res.type();

		const std::string name3 = res_type.name3();
		if ( res_counts.find(name3) != res_counts.end() ) {
			res_counts[ name3 ]++;
		}

		for (core::Size atomNum = 1, natom = pose.residue(resNum).natoms(); atomNum <= natom; ++atomNum){
			if (!res.heavyatom_is_an_acceptor(atomNum) && !res.atom_is_polar_hydrogen(atomNum)) continue;
			if (pdb_info.temperature(resNum, atomNum) > 30) continue;
			if (pdb_info.occupancy(resNum, atomNum) < 1) continue;

			id::AtomID at(atomNum, resNum);
			const Real vsasa = atom_sasa[at];

			//const core::chemical::AtomType& atom_type = res_type.atom_type(atomNum);

			if (vsasa > burial_cutoff) continue;
			utility::vector1<scoring::hbonds::HBondCOP> hbonds = hbond_set.atom_hbonds(at, false /*include only allowed*/);
			bool hbonded = false;
			for (utility::vector1<scoring::hbonds::HBondCOP>::iterator h = hbonds.begin(), end = hbonds.end();
			     h != end; h++) {
				if (hb_eval.evaluate(pose, *h)){
					hbonded = true;
					break;
				}
			}
			if (!hbonded) {
				++buns;
			}
		}
	}
	for (boost::unordered_map<std::string, Size>::const_iterator rc = res_counts.begin(); rc != res_counts.end(); rc++ ) {
		std::cout << rc->first << ", " << rc->second << ", ";
	}
	std::cout << buns << std::endl;
}


void buns_stats(
	const std::vector<std::string>& pdb_list,
	const Real AHdist_geom_eval_threshold,
	const hbond_strat& hb_eval,
	VarSolDistSasaCalculator& vsasa_calc,
	const Real burial_cutoff
) {
	for (core::Size i = 0, size = pdb_list.size(); i != size; i++) {
		pose::Pose pose;
		import_pose::pose_from_pdb(pose, pdb_list[i]);
		vsasa_bunsats(pose,
		              AHdist_geom_eval_threshold,
		              hb_eval,
		              vsasa_calc,
		              burial_cutoff);
	}
}

int main(int argc, char* argv[])
{
	try {

	devel::init(argc, argv);

	basic::options::option.add(AHdist_opkey, "acceptor-hydrogen distance").def(2.8);
	basic::options::option.add(AHD_opkey, "acceptor-hydrogen-donor angle").def(90.0);

	const Real AHdist = basic::options::option[AHdist_opkey].value();
	const Real AHD = basic::options::option[AHD_opkey].value();

	const hbond_geom_strat hb_geom_strat1(AHdist, AHD);

	std::vector<std::string> pdb_list;
	read_FileVector(basic::options::OptionKeys::in::file::l, pdb_list);

	VarSolDistSasaCalculator vsasa_calc;
	const Real burial_cutoff = 0.0;

	buns_stats(pdb_list, AHdist, hb_geom_strat1, vsasa_calc, burial_cutoff);

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	return 0;
}
