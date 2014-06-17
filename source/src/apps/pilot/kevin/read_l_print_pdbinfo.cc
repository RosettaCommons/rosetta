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
#include <unordered_map>

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
//#include <core/pose/metrics/CalculatorFactory.hh>

#include <devel/buns/BuriedUnsatHbondFilter2.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>

#include <core/types.hh> // Distance
//#include <boost/dynamic_bitset.hpp>
#include <utility/options/keys/FileVectorOptionKey.hh>

#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <core/pose/util.hh>
//typedef std::vector<boost::dynamic_bitset<> > bitset

//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static basic::Tracer TR("apps.pilot.kevin.buns");

//basic::options::BooleanOptionKey const use_varsoldist_sasa_calc( "use_varsoldist_sasa_calc" );
//basic::options::BooleanOptionKey const water_dist_H( "water_dist_H_cutoff" );
//basic::options::BooleanOptionKey const water_dist_O( "water_dist_O_cutoff" );
basic::options::FileVectorOptionKey const nat_list_opkey("nat_list");
basic::options::FileVectorOptionKey const rlx_list_opkey("rlx_list");

using namespace core;
using namespace core::pose::metrics;
using namespace devel::vardist_solaccess;
using core::Distance;
using platform::Size;

static void read_FileVector(const basic::options::FileVectorOptionKey& fv_key, std::vector<std::string>& str_vec) {
	if(!basic::options::option[basic::options::OptionKeys::in::file::l].user()) return;
	std::string list_file = basic::options::option[basic::options::OptionKeys::in::file::l]()[1];
	utility::io::izstream list;
	list.open(list_file,std::_S_in);
	while(!list.eof()) {
		std::string line;
		getline(list,line);
		if(line.size() > 0) str_vec.push_back(line);
	}
	list.close();
}

///@brief
class xtal_water_bunsat : public protocols::moves::Mover {

private:
	static constexpr const Real AHdist_geom_eval_threshold = 2.8;
	static constexpr const Real AHD_geom_eval_threshold = 90;

	bool evaluate_hbond_geom(Real AHdist, Real AHD) {
		return(AHdist <= AHdist_geom_eval_threshold &&
		       AHD >= AHD_geom_eval_threshold);
	}

public:
	xtal_water_bunsat()
	{
	}
	~xtal_water_bunsat(){};
	
	std::string
	get_name() const { return "xtal_water_bunsat"; }

	void
	apply(core::pose::Pose & pose, ){

		pose.update_residue_neighbors();
		scoring::ScoreFunctionOP scorefxn = scoring::getScoreFunction();
		scorefxn->score(pose);
		pose.energies();

		const Real AHdist_threshold = 4.0;

		const Real max_nbr_radius = core::pose::pose_max_nbr_radius( pose );
		const Real max_covalent_hbond_length = 1.5;
		const Real neighbor_cutoff = 2 * max_nbr_radius + max_covalent_hbond_length + AHdist_threshold;

		conformation::PointGraphOP pg = new conformation::PointGraph;
		conformation::residue_point_graph_from_conformation( pose.conformation(), *pg );
		conformation::find_neighbors<core::conformation::PointGraphVertexData,
			core::conformation::PointGraphEdgeData>( pg, neighbor_cutoff );

		Size nres = pose.total_residue();

		utility::vector1< VarSolDRotamerDotsOP > rotamer_dots(nres);
		for (Size resNum = 1; resNum <= nres; ++resNum ) {
			rotamer_dots[resNum] = new VarSolDRotamerDots(
					new conformation::Residue( pose.residue(resNum) ), true);
			rotamer_dots[resNum]->increment_self_overlap();
		}

		// find hbonds and intersect residues
		scoring::hbonds::HBondSet hbond_set;
		const scoring::hbonds::HBondDatabaseCOP hb_database = scoring::hbonds::HBondDatabase::get_database();
		const scoring::TenANeighborGraph& tenA_neighbor_graph(pose.energies().tenA_neighbor_graph());

		for (Size lowerResNum = 1; lowerResNum < nres; ++lowerResNum ) {
		const conformation::Residue& lowerRes(pose.residue(lowerResNum));
		const Size nbl = tenA_neighbor_graph.get_node(lowerResNum)->
				num_neighbors_counting_self_static();
			for (conformation::PointGraph::UpperEdgeListConstIter
					ue  = pg->get_vertex(lowerResNum).const_upper_edge_list_begin(),
					ue_end = pg->get_vertex(lowerResNum).const_upper_edge_list_end();
					ue != ue_end; ++ue ) {
				Size upperResNum = ue->upper_vertex();
    	    	const conformation::Residue& upperRes(pose.residue(lowerResNum));
    	    	const Size nbu = tenA_neighbor_graph.get_node(upperResNum)->
					num_neighbors_counting_self_static();
    	        scoring::hbonds::identify_hbonds_1way_AHdist(*hb_database,
					lowerRes, upperRes, nbl, nbu, AHdist_threshold, hbond_set);
    	        scoring::hbonds::identify_hbonds_1way_AHdist(*hb_database,
					upperRes, lowerRes, nbu, nbl, AHdist_threshold, hbond_set);
				rotamer_dots[lowerResNum]->intersect_residues(*rotamer_dots[upperResNum]);
			}
		}

		id::AtomID_Map< Real > vsasa_map;
		core::pose::initialize_atomid_map(vsasa_map, pose, 0.0);
		for (Size resNum = 1; resNum <= nres; ++resNum) {
			for ( Size atomNum = 1, natom = pose.residue(resNum).natoms(); atomNum <= natom; ++atomNum ){
				id::AtomID at(atomNum, resNum);
				vsasa_map[at] = rotamer_dots[resNum]->msas_for_atom(atomNum);
			}
		}

		std::vector<id::AtomID> buns;

		const core::pose::PDBInfo& pdb_info = *(pose.pdb_info());
		for (Size resNum = 1; resNum < nres; ++resNum) {
			const core::conformation::Residue& res = pose.residue(resNum);
			const core::chemical::ResidueType& res_type = res.type();
			for (Size atomNum = 1, natom = pose.residue(resNum).natoms(); atomNum <= natom; ++atomNum){
				if (pdb_info.temperature(resNum, atomNum) > 30) continue;
				if (pdb_info.occupancy(resNum, atomNum) < 1) continue;

				id::AtomID at(atomNum, resNum);
				const Real vsasa = vsasa_map[at];
				utility::vector1<scoring::hbonds::HBondCOP> hbonds = hbond_set.atom_hbonds(at, false);

				for (auto h : hbonds) {
					Real AHdist = h->get_HAdist(pose);
					Real AHD = h->get_AHDangle(pose);
					if (evaluate_hbond_geom(AHdist, AHD)){
						h->show(pose, true, std::cout);
						break;
					}
				}

				//const core::chemical::AtomType& atom_type = res_type.atom_type(atomNum);
				//const char& raw_icode = pdb_info.icode(resNum);
				//utility::file::FileName filename(pdb_info.modeltag());
				//std::cout << filename.base() << ","
				//	<< pdb_info.chain(resNum)                             << ","
				//	<< pdb_info.number(resNum)                            << ","
				//	<< raw_icode                                          << ","
				//	<< res_type.name3()                                   << ","
				//	<< atom_type.name()                                   << ","
				//	<< atom_type.element()                                << ","
				//	<< vsasa                                              << "\n";
			}
		}

		return;
	}


private:

};

typedef utility::pointer::owning_ptr< xtal_water_bunsat > xtal_water_bunsatOP;

int main( int argc, char* argv[] )
{
	try {

	basic::options::option.add(nat_list_opkey, "nat_list");
	basic::options::option.add(rlx_list_opkey, "nat_list");

	devel::init(argc, argv);

	// Basic -l parsing snippet	

	//std::list<std::string> file_list;
	//if(basic::options::option[basic::options::OptionKeys::in::file::l].user()) {
	//	std::string list_file = basic::options::option[basic::options::OptionKeys::in::file::l]()[1];
	//	utility::io::izstream list;
	//	list.open(list_file,std::_S_in);
	//	while(!list.eof())
	//	{
	//		std::string line;
	//		getline(list,line);
	//		if(line.size() > 0) file_list.push_back(line);
	//	}
	//	list.close();
	//} else if(basic::options::option[basic::options::OptionKeys::in::file::s].user()) {
	//	std::string mol_file = basic::options::option[basic::options::OptionKeys::in::file::s]()[1];
	//	file_list.push_back(mol_file);
	//}



	std::vector<std::string> nat_list;
	read_FileVector(nat_list_opkey, nat_list);
	std::vector<std::string> rlx_list;
	read_FileVector(rlx_list_opkey, rlx_list);

	for (auto i : nat_list) {
		std::cout << i << std::endl;
	}
	for (auto i : rlx_list) {
		std::cout << i << std::endl;
	}

	return 0;

	//protocols::jd2::JobDistributor::get_instance()->go(new xtal_water_bunsat);
	pose::Pose pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	import_pose::pose_from_pdb(pose, option[in::file::s]()[1]);
	import_pose::pose_from_pdb(pose, option[in::file::s]()[1]);
	xtal_water_bunsat xwb;
	xwb.apply(pose);

	TR << "************************d**o**n**e**************************************" << std::endl;

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}
