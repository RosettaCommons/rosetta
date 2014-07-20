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

#ifdef USEMPI
#include <mpi.h>

//#include <iomanip>
#include <algorithm>
#include <cctype>
#include <limits>
//#include <unordered_map>

#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
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

#include <core/scoring/sasa.hh>

#include <devel/buns/BuriedUnsatHbondFilter2.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>

#include <core/types.hh> // Distance

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>

#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <core/pose/util.hh>

#include <nlopt.hpp>

//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static basic::Tracer TR("apps.pilot.kevin.buns");

basic::options::FileVectorOptionKey const nat_list_opkey("nat_list");
basic::options::FileVectorOptionKey const rlx_list_opkey("rlx_list");

basic::options::RealOptionKey const AHdist_opkey("AHdist");
basic::options::RealOptionKey const AHD_opkey("AHD");
basic::options::BooleanOptionKey const energy_hbonds_opkey("energy_hbonds");

basic::options::RealOptionKey const O_inner_opkey("O_inner");
basic::options::RealOptionKey const O_outer_opkey("O_outer");
basic::options::RealOptionKey const N_inner_opkey("N_inner");
basic::options::RealOptionKey const N_outer_opkey("N_outer");
basic::options::RealOptionKey const C_inner_opkey("C_inner");
basic::options::RealOptionKey const Hpol_inner_opkey("Hpol_inner");
basic::options::RealOptionKey const Hpol_outer_opkey("Hpol_outer");
basic::options::RealOptionKey const Hnonpol_inner_opkey("Hnonpol_inner");

using namespace core;
using namespace core::pose::metrics;
using namespace devel::vardist_solaccess;
using core::Distance;
using platform::Size;

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


//void
//fill_atom_sasa(id::AtomID_Map<Real>& atom_sasa) {
//    if (vsasa) {
//        VarSolDistSasaCalculator vsasa_calc;
//        atom_sasa = vsasa_calc.calculate(pose);
//    }
//    else {
//        utility::vector1< core::Real > residue_sasa;
//        core::scoring::calc_per_atom_sasa(pose, atom_sasa, residue_sasa, probe_radius);
//    }
//}


	core::Size
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
			//rotamer_dots[lowerResNum]->intersect_residues(*rotamer_dots[upperResNum]);
		}
	}

	core::Size buns = 0;
	const pose::PDB_Info& pdb_info = *(pose.pdb_info());
	for (core::Size resNum = 1; resNum <= nres; ++resNum) {
		const conformation::Residue& res = pose.residue(resNum);
		const chemical::ResidueType& res_type = res.type();
		for (core::Size atomNum = 1, natom = pose.residue(resNum).natoms(); atomNum <= natom; ++atomNum){
			if (!res.heavyatom_is_an_acceptor(atomNum) && !res.atom_is_polar_hydrogen(atomNum)) continue;
			if (pdb_info.temperature(resNum, atomNum) > 30) continue;
			if (pdb_info.occupancy(resNum, atomNum) < 1) continue;

			id::AtomID at(atomNum, resNum);
			const Real vsasa = atom_sasa[at];

			const core::chemical::AtomType& atom_type = res_type.atom_type(atomNum);
			//std::cout << atom_type.name() << "\n";
			//std::cout << "	vsasa: " << vsasa << "\n";

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
	return buns;
}


int vsasa_separation(
	const std::vector<std::string>& nat_list,
	const std::vector<std::string>& rlx_list,
	const Real AHdist_geom_eval_threshold,
	const hbond_strat& hb_eval,
	VarSolDistSasaCalculator& vsasa_calc,
	const Real burial_cutoff
) {
	int separation = 0;

	for (core::Size i = 0, size = nat_list.size(); i < size; i++) {
		pose::Pose pose;
		//poses.push_back(new pose::Pose);
		//pose::Pose& pose = *poses.back();
		//std::cout << i << std::endl;
		//std::cout << nat_list[i] << std::endl;
		import_pose::pose_from_pdb(pose, nat_list[i]);
		core::Size nat_buns = vsasa_bunsats(pose,
		                        AHdist_geom_eval_threshold,
		                        hb_eval,
					            vsasa_calc,
			                    burial_cutoff);
		import_pose::pose_from_pdb(pose, rlx_list[i]);
		core::Size rlx_buns = vsasa_bunsats(pose,
		                        AHdist_geom_eval_threshold,
		                        hb_eval,
					            vsasa_calc,
			                    burial_cutoff);
		//std::string nat_vs_rlx = nat_buns < rlx_buns ? "nat_less" : (nat_buns == rlx_buns ? "nat_same" : "nat_more");
		//std::cout << nat_buns << "," << rlx_buns << "," << nat_vs_rlx << "\n";
		if (nat_buns < rlx_buns) {
			separation--;
		}
		else if (nat_buns > rlx_buns) {
			separation++;
		}
	}

	//std::cout << "Slave batch separation: " << separation << std::endl;
	return separation;
}

double vsasa_objective_master(const std::vector<double> &x, std::vector<double> &grad, void* func_data) {
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD,&world_size);
	//int num_slaves = world_size - 1;

	std::vector<double> mpi_data = x;

	// broadcast vsasa parameters
	MPI_Bcast(&mpi_data[0], x.size(), MPI_DOUBLE, 0,
		MPI_COMM_WORLD);

	int numPdbs = *( reinterpret_cast<int*>(func_data) );

	int list_index = 0;
	int batch_size = 10;

	// send initial batch to every node
	for (int rank = 1; rank < world_size; rank++) {
		if (list_index >= numPdbs) {
			break;
		}
		MPI_Send(&list_index, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
		list_index += batch_size;
	}

	int separation = 0;
	int subseparation;
	MPI_Status status;

	// dispatch batches of structures
	while (list_index < numPdbs) {
		MPI_Recv(&subseparation, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		separation += subseparation;
		MPI_Send(&list_index, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		list_index += batch_size;
	}
	for (int rank = 1; rank < world_size; rank++) {
		MPI_Recv(&subseparation, 1, MPI_INT, rank, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		separation += subseparation;
		int obfunc_done = -1;
		MPI_Send(&obfunc_done, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
	}

	std::cout << "Master params evaluation:\n"
			  << "\t\tAHdist = " << x[0] << "\n"
			  << "\t\tAHD = " << x[1] << "\n"
			  << "\t\tO_inner = " << x[2] << "\n"
			  << "\t\tO_outer = " << x[3] << "\n"
			  << "\t\tN_inner = " << x[4] << "\n"
			  << "\t\tN_outer = " << x[5] << "\n"
			  << "\t\tC_inner = " << x[6] << "\n"
			  << "\t\tHpol_inner = " << x[7] << "\n"
			  << "\t\tHpol_outer = " << x[8] << "\n"
			  << "\t\tHnonpol_inner = " << x[9] << "\n"
			  << "\tseparation = " << separation << std::endl;

	return static_cast<double>(separation);
}

//int vsasa_objective_slave(std::vector<double> &x, std::vector<double> &grad, void *my_func_data) {
//
//	// receive vsasa parameters
//
//	// receive batches of structures, evaluate batch, send sub-separation, repeat
//
//}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	try {

	basic::options::option.add(nat_list_opkey, "nat_list");
	basic::options::option.add(rlx_list_opkey, "rlx_list");

	basic::options::option.add(AHdist_opkey, "acceptor-hydrogen distance").def(2.8);
	basic::options::option.add(AHD_opkey, "acceptor-hydrogen-donor angle").def(90.0);
	basic::options::option.add(energy_hbonds_opkey, "use hbonds determined by energy").def(false);

	basic::options::option.add(O_inner_opkey, "O element vsasa collision radius").def(2.6);
	basic::options::option.add(O_outer_opkey, "O element vsasa interaction radius").def(3.0);
	basic::options::option.add(N_inner_opkey, "N element vsasa collision radius").def(2.7);
	basic::options::option.add(N_outer_opkey, "N element vsasa interaction radius").def(3.1);
	basic::options::option.add(C_inner_opkey, "C element vsasa collision radius").def(3.1);
	basic::options::option.add(Hpol_inner_opkey, "Hpol vsasa collision radius").def(1.7);
	basic::options::option.add(Hpol_outer_opkey, "Hpol vsasa interaction radius").def(2.1);
	basic::options::option.add(Hnonpol_inner_opkey, "Hnonpol vsasa collision radius").def(2.4);

	devel::init(argc, argv);

	Real AHdist = basic::options::option[AHdist_opkey].value();
	Real AHD = basic::options::option[AHD_opkey].value();

	Real O_inner = basic::options::option[O_inner_opkey].value();
	Real O_outer = basic::options::option[O_outer_opkey].value();
	Real N_inner = basic::options::option[N_inner_opkey].value();
	Real N_outer = basic::options::option[N_outer_opkey].value();
	Real C_inner = basic::options::option[C_inner_opkey].value();
	Real Hpol_inner = basic::options::option[Hpol_inner_opkey].value();
	Real Hpol_outer = basic::options::option[Hpol_outer_opkey].value();
	Real Hnonpol_inner = basic::options::option[Hnonpol_inner_opkey].value();

	//const hbond_geom_strat hb_geom_strat1(AHdist, AHD);
	Real burial_cutoff = 0.0;

	//VarSolDistSasaCalculator vsasa_calc;

	//vsasa_calc.set_element_radii("O", O_inner, O_outer, 5);
	//vsasa_calc.set_element_radii("N", N_inner, N_outer, 5);
	//vsasa_calc.set_element_radii("C", O_inner, C_inner, 1);
	//vsasa_calc.set_element_radii("H", Hnonpol_inner, Hnonpol_inner, 1);
	//vsasa_calc.set_atom_type_radii("Hpol", Hpol_inner, Hpol_outer, 5);

	std::vector<std::string> nat_list;
	read_FileVector(nat_list_opkey, nat_list);
	std::vector<std::string> rlx_list;
	read_FileVector(rlx_list_opkey, rlx_list);

	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	if (world_rank == 0) {
		std::cout << "I'm the master" << std::endl;
		/* do master stuff */

		/* BOBYQA */
		nlopt::opt local_opt(nlopt::LN_BOBYQA, 10);
		nlopt::opt opt(nlopt::G_MLSL_LDS, 10);
		//nlopt::opt opt(nlopt::LN_BOBYQA, 10);

		opt.set_local_optimizer(local_opt);
		//opt.set_xtol_abs1(0.01);

		std::cout << "Master setting reasonable bounds on objective function" << std::endl;

		std::cout << "Master setting lowerBounds" << std::endl;
	  std::vector<double> lb(10, 0.1) ;
		lb[0] = 2.5;
		lb[1] = 80;
	  opt.set_lower_bounds(lb);
		lb[4] = 1.0;
		lb[5] = 1.0;
		lb[6] = 1.0;
		// C_inner
		lb[7] = 1.0;
		lb[8] = 0.1;
		lb[9] = 1.0;
		lb[10] = 0.4;
	  opt.set_lower_bounds(lb);

		std::cout << "Master setting upperBounds" << std::endl;
	  std::vector<double> ub(10, 5.0);
		ub[1] = 180;
		opt.set_upper_bounds(ub);

		std::cout << "Master about to set objective function" << std::endl;
		size_t numPdbs = nat_list.size();
		opt.set_min_objective(&vsasa_objective_master, &numPdbs );

		std::cout << "Master has set objective function" << std::endl;
		//std::vector<double> x = {AHdist, O_inner, O_outer, N_inner, N_outer, C_inner, Hpol_inner, Hpol_outer, Hnonpol_inner};
		double tmp[] = {AHdist, AHD, O_inner, O_outer, N_inner, N_outer, C_inner, Hpol_inner, Hpol_outer, Hnonpol_inner};
		std::vector<double> x (tmp, tmp + sizeof(tmp) / sizeof(double));
		double minf;
		std::cout << "Master calling optimize()" << std::endl;
		nlopt::result result = opt.optimize(x, minf);

	    if (result < 0) {
	        printf("nlopt failed!\n");
	    }
	    else {
			std::cout << "Found minimum.\n"
			          << "AHdist = " << x[0] << "\n"
			          << "AHD = " << x[1] << "\n"
			          << "O_inner = " << x[2] << "\n"
			          << "O_outer = " << x[3] << "\n"
			          << "N_inner = " << x[4] << "\n"
			          << "N_outer = " << x[5] << "\n"
			          << "C_inner = " << x[6] << "\n"
			          << "Hpol_inner = " << x[7] << "\n"
			          << "Hpol_outer = " << x[8] << "\n"
			          << "Hnonpol_inner = " << x[9] << "\n"
			          << "separation = " << minf << std::endl;
	    }

		// send terminate signal
		x = std::vector<double>(10,0);
		MPI_Bcast(&x[0], x.size(), MPI_DOUBLE, 0,
			MPI_COMM_WORLD);
	}
	else {
		/* do slave stuff */
		std::cout << "Slaving away" << std::endl;
		while (true) {
			double x[10];
			// receive vsasa params or terminate signal
			MPI_Bcast(&x[0], 10, MPI_DOUBLE, 0,
				MPI_COMM_WORLD);
			// check for terminate signal
			bool allZero = true;
			for (size_t i = 0; i<10; i++) {
				if (x[i] != 0) {
					allZero = false;
					break;
				}
			}
			if (allZero) {
				// shut down
				break;
			}

			// set vsasa parameters from x
			AHdist = x[0];
			AHD = x[1];
			O_inner = x[2];
			O_outer = x[3];
			N_inner = x[4];
			N_outer = x[5];
			C_inner = x[6];
			Hpol_inner = x[7];
			Hpol_outer = x[8];
			Hnonpol_inner = x[9];

			const hbond_geom_strat hb_geom_strat1(AHdist, AHD);

			VarSolDistSasaCalculator vsasa_calc;
			vsasa_calc.set_element_radii("O", O_inner, O_outer, 5);
			vsasa_calc.set_element_radii("N", N_inner, N_outer, 5);
			vsasa_calc.set_element_radii("C", O_inner, C_inner, 1);
			vsasa_calc.set_element_radii("H", Hnonpol_inner, Hnonpol_inner, 1);
			vsasa_calc.set_atom_type_radii("Hpol", Hpol_inner, Hpol_outer, 5);

			// do work
			// receive either pdbs to work on or signal that the objective function is done
			int list_index = 0;
			size_t batch_size = 10;
			size_t numPdbs = nat_list.size();
			MPI_Status status;
			while (true) {
				MPI_Recv(&list_index, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				if (list_index == -1) {
					break;
				}
				std::vector<std::string>::iterator last = nat_list.begin() + list_index + batch_size;
				if (last > nat_list.end()) {
					last = nat_list.end();
				}
				std::vector<std::string> nat_sublist(nat_list.begin() + list_index, last);
				std::vector<std::string> rlx_sublist(rlx_list.begin() + list_index, rlx_list.begin() + list_index + batch_size);
				// calculate subseparation of list_index to list_index + batch_size
				int sep = vsasa_separation(nat_sublist, rlx_sublist, AHdist, hb_geom_strat1, vsasa_calc, burial_cutoff);
				MPI_Send(&sep, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			}
		}
	}

	MPI_Finalize();

	//int sep = vsasa_separation(nat_list, rlx_list, AHdist, hb_geom_strat1, vsasa_calc, burial_cutoff);

	//std::cout << "separation: " << sep << std::endl;

	} catch (utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}

#else

// no non-mpi version
int main(int argc, char* argv[]) {
	return 0;
}

#endif
