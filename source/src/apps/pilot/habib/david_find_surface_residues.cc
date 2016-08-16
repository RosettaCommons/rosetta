// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk + dj

#include <iostream>
#include <iomanip>
#include <algorithm>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <basic/MetricValue.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/rms_util.hh>

#include <core/pose/metrics/CalculatorFactory.hh>
#include <numeric/random/random.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/scoring/dssp/Dssp.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <string>
#include <ObjexxFCL/string.functions.hh>

#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );


OPT_KEY( Integer, max_residues )
OPT_KEY( Real, min_sasa )
OPT_KEY( String, contact_list )

//set to store pdb info keys
std::vector<std::string> surface;
std::set <std::string> interface;

//stores resid of the ligand residue
void register_metrics() {

	core::pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );

}

void
setup_secstruct_dssp( pose::Pose & pose )
{
	core::scoring::dssp::Dssp dssp( pose );
	ObjexxFCL::FArray1D_char dssp_secstruct( pose.total_residue() );
	dssp.dssp_reduced( dssp_secstruct );
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		pose.set_secstruct(i,  dssp_secstruct(i) );
	}

}


/// General testing code
int
main( int argc, char * argv [] )
{
	try {

		std::vector<std::string> surface;
		std::set <std::string> interface;
		NEW_OPT( max_residues, "Maximum number of residues to return", 1000);
		NEW_OPT( min_sasa, "Minimum SASA to classify a residue as surface", 10);
		NEW_OPT ( contact_list, "File name for optional list of contact residues to check","");

		using namespace core;
		using namespace core::scoring;

		devel::init(argc, argv);
		register_metrics();
		pose::Pose pose;
		std::string const input_pdb_name ( basic::options::start_file() );
		core::import_pose::pose_from_file( pose, input_pdb_name , core::import_pose::PDB_file);

		std::string const cfilename = option[ contact_list ];
		if ( cfilename != "" ) {
			std::ifstream ifs(cfilename.c_str(), std::ifstream::in);
			if ( !ifs.is_open() ) {
				std::cout<< "Error opening contact list file "<<cfilename<<std::endl;
				return -100;
			}
			//ifb.open (cfilename,std::ios::in);
			//std::ostream ios(&ifb);
			std::string intres;
			while ( ifs.good() ) {
				ifs >> intres;
				interface.insert(intres);
			}
		}

		core::Size sasa_threshold = option [ min_sasa ];
		core::Size max_resi = option [ max_residues ];
		basic::MetricValue< utility::vector1< Real > > resisasa;
		pose.metric( "sasa", "residue_sasa", resisasa );
		setup_secstruct_dssp(pose);

		for ( Size i = 1; i <= pose.total_residue(); ++i ) {
			int close = 0;
			for ( Size j = 1; j <= pose.total_residue(); ++j ) {
				if ( i == j ) continue;
				//make sure it's not within 12A of the interface
				std::ostringstream residuestream;
				residuestream << pose.pdb_info()->chain(j) << pose.pdb_info()->number(j);
				std::string res_id = residuestream.str();
				if ( interface.find(res_id) != interface.end() ) {
					if ( pose.residue(i).xyz( pose.residue(i).nbr_atom() ).distance( pose.residue(j).xyz( pose.residue(j).nbr_atom() ) ) <= 12 ) {
						close = 1;
						break;
					}
				}
				//make sure it's not within 12A of a *meric interface
				if ( pose.pdb_info()->chain(i) == pose.pdb_info()->chain(j) ) continue;
				if ( pose.residue(i).xyz( pose.residue(i).nbr_atom() ).distance( pose.residue(j).xyz( pose.residue(j).nbr_atom() ) ) <= 12 ) {
					close = 1;
				}

			}
			if ( !close ) {
				if ( resisasa.value()[i]>= sasa_threshold && pose.secstruct(i) != 'L' ) {
					std::ostringstream residuestream;
					residuestream << pose.pdb_info()->chain(i) << pose.pdb_info()->number(i);
					std::string res_id = residuestream.str();
					if ( interface.find(res_id) == interface.end() ) {
						//std::cout<<res_id<<std::endl;
						surface.push_back(res_id);
					}
				}
			}
		}

		std::filebuf fb;
		std::stringstream filename;
		filename<<option[ OptionKeys::out::output_tag ]()<<".randomsurface";
		fb.open (filename.str().c_str(),std::ios::out);
		std::ostream os(&fb);


		if ( max_resi < surface.size() ) {
			std::set<std::string> indeces;
			for ( core::Size i = 0; i < max_resi; i++ ) {
				int r=(int) (numeric::random::uniform() * surface.size());
				std::string rname(surface[r]);
				if ( indeces.find(surface[r]) == indeces.end() ) {
					indeces.insert(surface[r]);

					//remove nearest neighbors
					core::Size pos = 0;
					for ( Size j = 1; j <= pose.total_residue(); ++j ) {
						std::ostringstream residuestream;
						residuestream << pose.pdb_info()->chain(j) << pose.pdb_info()->number(j);
						std::string res_id = residuestream.str();
						if ( !res_id.compare(surface[r]) ) {
							pos = j;
							break;
						}
					}
					for ( std::vector<std::string>::iterator it2=surface.begin(); it2 != surface.end(); it2++ ) {
						core::Size(pos2) = 0;
						for ( Size j = 1; j <= pose.total_residue(); ++j ) {
							std::ostringstream residuestream;
							residuestream << pose.pdb_info()->chain(j) << pose.pdb_info()->number(j);
							std::string res_id = residuestream.str();
							if ( !res_id.compare(*it2) ) {
								pos2 = j;
								break;
							}
						}
						if ( pos == pos2 ) continue;
						if ( pose.residue(pos).xyz( pose.residue(pos).nbr_atom() ).distance( pose.residue(pos2).xyz( pose.residue(pos2).nbr_atom() ) ) <= 12 ) {
							if ( it2 != surface.begin() ) --it2;
							surface.erase(it2);
							if ( max_resi == surface.size() ) max_resi--;
						}

					}
				} else {
					i--;
				}
			}
			//sort (indeces.begin(), indeces.end());
			for ( std::set<std::string>::iterator it=indeces.begin(); it != indeces.end(); it++ ) {
				os<<*it<<std::endl;
			}
		} else {
			for ( std::vector<std::string>::iterator it=surface.begin(); it != surface.end(); it++ ) {
				os<<*it<<std::endl;
			}
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}
