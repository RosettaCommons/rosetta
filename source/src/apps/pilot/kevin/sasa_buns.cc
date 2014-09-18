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
//#include <unordered_map>

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

//#include <basic/MetricValue.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/id/AtomID_Map.hh>

//#include <core/graph/Graph.hh>
//#include <core/scoring/EnergyGraph.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

#include <core/scoring/hbonds/hbonds.hh>

#include <core/scoring/sasa.hh>

#include <devel/buns/BuriedUnsatHbondFilter2.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>

#include <core/types.hh> // Distance
//#include <boost/dynamic_bitset.hpp>
//#include <boost/zip_iterator.hpp>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>

#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <core/pose/util.hh>
//typedef std::vector<boost::dynamic_bitset<> > bitset

//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static thread_local basic::Tracer TR( "apps.pilot.kevin.buns" );

basic::options::FileVectorOptionKey const nat_list_opkey("nat_list");
basic::options::FileVectorOptionKey const rlx_list_opkey("rlx_list");

basic::options::RealOptionKey const AHdist_opkey("AHdist");
basic::options::RealOptionKey const AHD_opkey("AHD");
basic::options::RealOptionKey const probe_radius_opkey("probe_radius");

basic::options::BooleanOptionKey const use_varsoldist_sasa_calc( "vsasa" );
basic::options::BooleanOptionKey const energy_hbonds_opkey( "energy_hbonds" );

//enum BUNSCOMPARISON {NAT_LESS, NAT_SAME, NAT_MORE};

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
	//hbond_strat(const hbond_strat&) = default;
	//hbond_strat(hbond_strat&&) = default;
	//hbond_strat& operator=(const hbond_strat&) & = default;
	//hbond_strat& operator=(hbond_strat&&) & = default;
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
		Real AHdist_geom_eval_theeshold,
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


core::Size
sasa_bunsats(
	pose::Pose& pose,
	const Real AHdist_threshold,
	const hbond_strat& hb_eval,
	const Real probe_radius,
	const Real burial_cutoff
) {
	pose.update_residue_neighbors();
	scoring::ScoreFunctionOP scorefxn = scoring::get_score_function();
	scorefxn->score(pose);
	pose.energies();

	const Real max_nbr_radius = pose::pose_max_nbr_radius(pose);
	//constexpr const Real max_covalent_hbond_length = 1.5;
	const Real max_covalent_hbond_length = 1.5;
	const Real neighbor_cutoff = 2 * max_nbr_radius + max_covalent_hbond_length + AHdist_threshold;

	conformation::PointGraphOP pg = new conformation::PointGraph;
	conformation::residue_point_graph_from_conformation(pose.conformation(), *pg);
	conformation::find_neighbors<conformation::PointGraphVertexData,
		conformation::PointGraphEdgeData>(pg, neighbor_cutoff);

	id::AtomID_Map<Real> atom_sasa;

	if(basic::options::option[use_varsoldist_sasa_calc]) {
		VarSolDistSasaCalculator vsasa_calc;
		atom_sasa = vsasa_calc.calculate(pose);
	}
	else {
		utility::vector1< core::Real > residue_sasa;
		core::scoring::calc_per_atom_sasa(pose, atom_sasa, residue_sasa, probe_radius);
	}

	Size nres = pose.total_residue();

	// find hbonds and intersect residues
	scoring::hbonds::HBondSet hbond_set;

	const scoring::hbonds::HBondDatabaseCOP hb_database = scoring::hbonds::HBondDatabase::get_database();
	const scoring::TenANeighborGraph& tenA_neighbor_graph(pose.energies().tenA_neighbor_graph());

	for (Size lowerResNum = 1; lowerResNum <= nres; ++lowerResNum ) {
		const conformation::Residue& lowerRes(pose.residue(lowerResNum));
		const Size nbl = tenA_neighbor_graph.get_node(lowerResNum)->
			num_neighbors_counting_self_static();
		for (conformation::PointGraph::UpperEdgeListConstIter
				ue  = pg->get_vertex(lowerResNum).const_upper_edge_list_begin(),
				ue_end = pg->get_vertex(lowerResNum).const_upper_edge_list_end();
				ue != ue_end; ++ue ) {
			Size upperResNum = ue->upper_vertex();
			const conformation::Residue& upperRes(pose.residue(upperResNum));
			const Size nbu = tenA_neighbor_graph.get_node(upperResNum)->
				num_neighbors_counting_self_static();
			scoring::hbonds::identify_hbonds_1way_AHdist(*hb_database,
				lowerRes, upperRes, nbl, nbu, AHdist_threshold, hbond_set);
			scoring::hbonds::identify_hbonds_1way_AHdist(*hb_database,
				upperRes, lowerRes, nbu, nbl, AHdist_threshold, hbond_set);
			//rotamer_dots[lowerResNum]->intersect_residues(*rotamer_dots[upperResNum]);
		}
	}

	Size buns = 0;
	const pose::PDBInfo& pdb_info = *(pose.pdb_info());
	for (Size resNum = 1; resNum <= nres; ++resNum) {
		const conformation::Residue& res = pose.residue(resNum);
		//const chemical::ResidueType& res_type = res.type();
		for (Size atomNum = 1, natom = pose.residue(resNum).natoms(); atomNum <= natom; ++atomNum){
			if (!res.heavyatom_is_an_acceptor(atomNum) && !res.atom_is_polar_hydrogen(atomNum)) continue;
			if (pdb_info.temperature(resNum, atomNum) > 30) continue;
			if (pdb_info.occupancy(resNum, atomNum) < 1) continue;

			id::AtomID at(atomNum, resNum);
			const Real sasa = atom_sasa[at];

			//const core::chemical::AtomType& atom_type = res_type.atom_type(atomNum);

			if (sasa > burial_cutoff) continue;
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
	return buns;
}


int sasa_separation(const std::vector<std::string>& nat_list,
                      const std::vector<std::string>& rlx_list,
                      const Real AHdist_geom_eval_threshold,
                      const hbond_strat& hb_eval,
                      const Real probe,
                      const Real burial_cutoff) {
	int separation = 0;
	for (Size i = 1, size = nat_list.size(); i != size; i++) {
		pose::Pose pose;
		//poses.push_back(new pose::Pose);
		//pose::Pose& pose = *poses.back();
		//std::cout << i << std::endl;
		//std::cout << nat_list[i] << std::endl;
		import_pose::pose_from_pdb(pose, nat_list[i]);
		Size nat_buns = sasa_bunsats(pose,
		                             AHdist_geom_eval_threshold,
		                             hb_eval,
							         probe,
			                         burial_cutoff);
		import_pose::pose_from_pdb(pose, rlx_list[i]);
		Size rlx_buns = sasa_bunsats(pose,
		                             AHdist_geom_eval_threshold,
		                             hb_eval,
							         probe,
			                         burial_cutoff);
		//std::string nat_vs_rlx = nat_buns < rlx_buns ? "nat_less" : (nat_buns == rlx_buns ? "nat_same" : "nat_more");
		//std::cout << nat_buns << "," << rlx_buns << "," << nat_vs_rlx << "\n";
		if (nat_buns < rlx_buns) {
			separation++;
		}
		else if (nat_buns > rlx_buns) {
			separation--;
		}
	}
	return separation;
}


int main( int argc, char* argv[] ) {
	try {

	basic::options::option.add(nat_list_opkey, "list of crystal native pdbs");
	basic::options::option.add(rlx_list_opkey, "list of relaxed native pdbs");

	basic::options::option.add(AHdist_opkey, "acceptor-hydrogen distance").def(2.8);
	basic::options::option.add(AHD_opkey, "acceptor-hydrogen-donor angle").def(90.0);
	basic::options::option.add(probe_radius_opkey, "probe radius").def(1.4);

	basic::options::option.add( use_varsoldist_sasa_calc, "var sol d sasa calculator" ).def(false);

	basic::options::option.add( energy_hbonds_opkey, "hbonds by energy" ).def(false);

	devel::init(argc, argv);

	Real AHdist = basic::options::option[AHdist_opkey].value();
	Real AHD = basic::options::option[AHD_opkey].value();
	Real probe_radius = basic::options::option[probe_radius_opkey].value();
	const hbond_geom_strat hb_geom_strat1(AHdist, AHD);
	Real burial_cutoff = 0.0;

	hbond_strat* hb_strat;
	if(basic::options::option[energy_hbonds_opkey]) {
		hb_strat = new hbond_energy_strat(0.0);
	}
	else {
		hb_strat = new hbond_geom_strat(AHdist, AHD);
	}

	std::vector<std::string> nat_list;
	read_FileVector(nat_list_opkey, nat_list);
	std::vector<std::string> rlx_list;
	read_FileVector(rlx_list_opkey, rlx_list);

	int sep = sasa_separation(nat_list, rlx_list, AHdist, *hb_strat, probe_radius, burial_cutoff);
	delete hb_strat;

	TR << "separation: " << sep << std::endl;

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}
