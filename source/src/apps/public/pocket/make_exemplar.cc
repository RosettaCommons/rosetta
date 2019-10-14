// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>
#include <algorithm>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/NoOutputJobOutputter.hh>

#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>



// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>
#include <utility/io/mpistream.hh>
#include <core/kinematics/MoveMap.hh>
#include <numeric/random/random.fwd.hh>

//Numeric Headers
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>

//Protocol Headers
#include <protocols/pockets/PocketGrid.hh>
//#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;
using namespace conformation;
using namespace core::pose::datacache;
using namespace core::id;
using namespace protocols::rigid;
using namespace protocols::simple_moves;
using namespace protocols;
using namespace protocols::moves;


OPT_KEY( String, central_relax_pdb_num )
//OPT_KEY( Integer, pocket_num_angles )
OPT_KEY( Integer, max_atoms )
OPT_KEY( Integer, min_atoms)

static basic::Tracer TR( "apps.public.make_exemplar.main", basic::t_debug );

//This Mover creates an exemplar
class ExemplarMover : public moves::Mover {

public:
	ExemplarMover();

	~ExemplarMover() override;

	MoverOP clone() const override;
	MoverOP fresh_instance() const override;

	void apply( Pose & pose ) override;
	std::string get_name() const override;
	void test_move( Pose & pose ) override
	{
		apply(pose);
	}
};

ExemplarMover::ExemplarMover() :
	Mover( "benchmark" )
{

}

ExemplarMover::~ExemplarMover() = default;

MoverOP ExemplarMover::clone() const {
	return utility::pointer::make_shared< ExemplarMover >( *this );
}

MoverOP ExemplarMover::fresh_instance() const {
	return utility::pointer::make_shared< ExemplarMover >();
}

void ExemplarMover::apply( Pose & pose ) {
	using namespace pose;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	core::Size const max_size = option[ max_atoms ];
	core::Size const min_size = option[ min_atoms ];
	core::Size const angles = option[ basic::options::OptionKeys::pocket_grid::pocket_num_angles ];
	numeric::xyzVector<core::Size> best_exemplar(2,1,0);
	std::string const resid_c = option[central_relax_pdb_num];
	std::string resid_tag = resid_c;
	//replace comma with dash in file output name
	std::replace( resid_tag.begin(), resid_tag.end(), ',', '-');
	std::string best_exemplar_name = get_current_tag().substr(0,get_current_tag().find_last_of("_")) +".pdb."+ resid_tag + ".exemplar.pdb";
	std::string top_exemplar_name = "";
	std::string best_grid_name = get_current_tag().substr(0,get_current_tag().find_last_of("_")) +".pdb."+ resid_tag + ".grid.pdb";

	std::vector< conformation::ResidueCOP > residues = protocols::pockets::PocketGrid::getRelaxResidues(pose, resid_c);

	for ( Size i=1; i<=angles; i++ ) {
		protocols::pockets::PocketGrid comparison_pg( residues );
		if ( resid_c.length() ) {
			resid_tag = resid_c;
			while ( true ) {
				std::size_t fpos( resid_tag.find(','));
				if ( fpos == std::string::npos ) break;
				resid_tag[fpos]='-';
			}

			if ( i>1 ) {
				comparison_pg.randomAngle();
			} else {
				comparison_pg.zeroAngle();
			}
			std::string i_str = std::to_string(i);
			std::string out_exfname = get_current_tag().substr(0,get_current_tag().find_last_of("_")) +"_pdb_" + i_str +"." + resid_tag + ".exemplar.pdb";
			comparison_pg.autoexpanding_pocket_eval( residues, pose );
			core::Size clusts = comparison_pg.GetExemplarNumClusters();
			std::string clust_string = std::to_string(clusts);
			std::cout << "File: " << out_exfname << std::endl;
			std::cout <<"Number of Clusters: "<< clust_string << std::endl;
			core::Size num_atoms = comparison_pg.GetExemplarNumAtoms();
			std::string num_atoms_string = std::to_string(num_atoms);
			std::cout <<"Number of atoms: " << num_atoms_string << std::endl;
			core::Size num_polar_atoms = comparison_pg.GetExemplarNumPolarAtoms();
			std::string num_polar_atoms_string = std::to_string(num_polar_atoms);
			std::cout <<"Number of polar atoms: " << num_polar_atoms_string <<std::endl;

			numeric::xyzVector<core::Size> current_exemplar(clusts, num_atoms, num_polar_atoms);
			std::string current_exemplar_name = out_exfname;
			if ( min_size <= num_atoms && num_atoms <= max_size && current_exemplar[0] == 1 ) {
				core::Real current_exemplar_polar = (core::Real)current_exemplar[2];
				core::Real best_exemplar_polar = (core::Real)best_exemplar[2];
				core::Real current_exemplar_atoms = (core::Real)current_exemplar[1];
				core::Real best_exemplar_atoms = (core::Real)best_exemplar[1];

				core::Real best_exemplar_ratio = best_exemplar_polar/best_exemplar_atoms;
				core::Real current_exemplar_ratio = current_exemplar_polar/current_exemplar_atoms;

				std::cout<<"Best ratio prior to current exemplar: " << std::to_string(best_exemplar_ratio) << std::endl;
				std::cout<<"Current exemplar ratio: " <<std::to_string(current_exemplar_ratio) << std::endl;

				if ( current_exemplar_ratio > best_exemplar_ratio ) {
					best_exemplar=current_exemplar;
					top_exemplar_name = current_exemplar_name;
					comparison_pg.dumpExemplarToFile( best_exemplar_name.c_str());
					if ( option[ OptionKeys::pocket_grid::pocket_dump_pdbs ]() ) {
						comparison_pg.dumpGridToFile(best_grid_name.c_str());
					}
				} else {
					continue;
				}
			} else {
				std::cout<<"Current exemplar not in size range or has 0 clusters" << std::endl;
			}


		}
	}
	std::cout<<"Name of best exemplar up to this point: " << top_exemplar_name << std::endl;
	std::cout<<"Number of clusters in best exemplar: " << std::to_string(best_exemplar[0]) << std::endl;
	std::cout<<"Number of atoms in best exemplar: " << std::to_string(best_exemplar[1]) << std::endl;
	std::cout<<"Number of polar atoms in best exemplar: "<<std::to_string(best_exemplar[2]) << std::endl;
}

std::string ExemplarMover::get_name() const {
	return "ExemplarMover";
}


/// General testing code
int main( int argc, char * argv [] ) {

	try{
		using namespace core;
		using namespace protocols::moves;
		using namespace scoring;
		using namespace basic::options;
		using namespace protocols::jobdist;
		using namespace basic::options::OptionKeys;
		using protocols::moves::MoverOP;

		NEW_OPT( central_relax_pdb_num, "target residue(s)", "-1");
		//NEW_OPT( pocket_num_angles, "number of times to rotate the grid", 1);
		NEW_OPT( max_atoms, "maximum number of atoms that the exemplar can have", 100);
		NEW_OPT( min_atoms, "minimum number of atoms that the exemplar can have", 1);
		// APL NOTE: Tracers cannot be written to before devel::init gets called. TR << "Calling init" << std::endl;
		//initializes Rosetta functions
		jd2::register_options();
		option.add_relevant( OptionKeys::in::file::fullatom );
		option.add_relevant( OptionKeys::in::file::movemap );
		devel::init(argc, argv);
		MoverOP protocol( new ExemplarMover() );
		protocols::jd2::JobDistributor::get_instance()->go( protocol,  utility::pointer::make_shared< jd2::NoOutputJobOutputter >() );
	}
catch (utility::excn::Exception const & e ) {
	//std::cerr << "caught exception " << e.msg() << std::endl;
	e.display();
	return -1;
}
	return 0;
}

