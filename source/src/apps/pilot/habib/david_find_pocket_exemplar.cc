// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author David Johnson

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

// Protocol Headers
#include <devel/init.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <protocols/pockets/PocketExemplarMultifunc.hh>
#include <core/optimization/ParticleSwarmMinimizer.hh>
#include <basic/options/option_macros.hh>

// Utility Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

#include <basic/options/option_macros.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>


#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>


using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

OPT_KEY( Integer, num_runs )
OPT_KEY( Integer, num_particles )
OPT_KEY( Integer, num_carbons )
OPT_KEY( Real, c_radius )
OPT_KEY( Real, rep_weight )


int main( int argc, char * argv [] ) {

	try {

		NEW_OPT( num_runs, "no. of runs for PSO", 100 );
		NEW_OPT( num_particles, "no. of particles for PSO", 100 );
		NEW_OPT( num_carbons, "Number of carbon atoms in exemplar", 0 );
		NEW_OPT( c_radius, "Radius of carbon atoms", 1.4 );
		NEW_OPT( rep_weight, "Penalty for clash with protein", 10.0 );
		devel::init(argc, argv);

		int particle_size = option[ num_particles ];
		int run_size = option[ num_runs ];
		core::Size nCarbons = option [ num_carbons ];
		core::Real c_rad = option [ c_radius ];
		core::Real repW = option [ rep_weight ];
		std::string const input_pdb_name ( basic::options::start_file() );
		std::string const resid (option[ OptionKeys::pocket_grid::central_relax_pdb_num ]  );

		if ( nCarbons == 0 ) {
			std::cerr<<"Must specify the number of carbon atoms for the exemplar\n";
			exit(10);
		}

		//these store the max and min values from for the vector of variables to optimize
		utility::vector1<core::Real> p_min(3*nCarbons);
		utility::vector1<core::Real> p_max(3*nCarbons);

		//declare an array of aprticles
		ParticleOPs particles;

		//A pocket multifunction for exemplar optimization
		protocols::pockets::PocketExemplarMultifunc pm(input_pdb_name, resid, c_rad, repW, p_min,p_max);

		std::cout<<std::endl<<p_min[1]<<" "<<p_max[1]<<" "<<p_min[2]<<" "<<p_max[2]<<" "<<p_min[3]<<" "<<p_max[3]<<std::endl;
		std::cout<<std::endl<<p_min[4]<<" "<<p_max[4]<<" "<<p_min[5]<<" "<<p_max[5]<<" "<<p_min[6]<<" "<<p_max[3]<<std::endl;
		//Make a minimizer with the constraints of the variables to optimize
		core::optimization::ParticleSwarmMinimizer pso(p_min, p_max);

		//run the PSO, get the particles
		particles = pso.run(run_size, pm, particle_size);

		//Particles are sorted based on fitness, so first is global optima
		ParticleOP p = particles[1];
		core::optimization::Particle parti(*p);
		core::Real fit_best = -(parti.fitness_pbest());
		utility::vector1<core::Real> best_vars(3*nCarbons);
		best_vars = parti.pbest();
		std::cout<<"Optimal Z: "<<fit_best<<std::endl;

		utility::io::ozstream outPDB_stream;
		outPDB_stream.open("test.pdb", std::ios::out);
		int counter=1;
		for ( int i = 1; i<(int)best_vars.size(); i+=3 ) {
			std::string concatenated_pdb_info;
			concatenated_pdb_info += "HETATM";
			std::stringstream  tmp;
			tmp<<counter;
			if ( counter<10 ) concatenated_pdb_info += "    ";
			else if ( counter<100 ) concatenated_pdb_info += "   ";
			else if ( counter<1000 ) concatenated_pdb_info += "  ";
			else if ( counter<10000 ) concatenated_pdb_info += " ";
			else concatenated_pdb_info += "";
			concatenated_pdb_info += tmp.str()+"  ";
			concatenated_pdb_info += " C  TMP A";
			tmp.str(std::string());
			tmp<<"1";
			concatenated_pdb_info += "   ";
			concatenated_pdb_info += tmp.str()+"  ";
			tmp.str(std::string());
			tmp<<"  "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<best_vars[i]<<std::setw(8)<<best_vars[i+1]<<std::setw(8)<<best_vars[i+2];
			tmp << "  1.00  2.03           C";
			tmp << std::endl;
			concatenated_pdb_info += tmp.str();
			counter++;
			outPDB_stream<<concatenated_pdb_info;
		}
		outPDB_stream.close();
		outPDB_stream.clear();
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
