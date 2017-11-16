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
#include <sstream>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <basic/MetricValue.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

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

static basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );


OPT_KEY( String, input_TMalignment )

//set to store pdb info keys
std::vector<std::string> surface;
std::set <std::string> interface;


/// General testing code
int
main( int argc, char * argv [] )
{
	try {

		std::vector<std::string> surface;
		std::set <std::string> interface;
		NEW_OPT ( input_TMalignment, "File name specifying the alignment from the last 3 lines of output of TMalign","");

		using namespace core;
		using namespace core::scoring;

		devel::init(argc, argv);
		pose::Pose pose;
		std::string const input_pdb_name ( basic::options::start_file() );
		core::import_pose::pose_from_file( pose, input_pdb_name );

		std::string const aln_file = option[ input_TMalignment ];
		if ( aln_file.length() == 0 ) {
			std::cerr<<"Error: No TMalignment file specified"<<std::endl;
			return -102;
		}
		std::ifstream ifs(aln_file.c_str(), std::ifstream::in);
		if ( !ifs.is_open() ) {
			std::cerr<< "Error opening TMalignment file "<<std::endl;
			return -100;
		}

		std::string target_aln;
		if ( ifs.good() ) {
			ifs >> target_aln;
			std::cout<<target_aln<<std::endl;
		} else {
			std::cerr<< "Error reading TMalignment file"<<std::endl;
			return -101;
		}
		std::string musashi_aln;
		if ( ifs.good() ) {
			ifs.ignore(1000, '\n');
			ifs.ignore(1000, '\n');
		} else {
			std::cerr<< "Error reading TMalignment file"<<std::endl;
			return -101;
		}
		if ( ifs.good() ) {
			ifs >> musashi_aln;
			std::cout<<musashi_aln<<std::endl;
		} else {
			std::cerr<< "Error reading TMalignment file"<<std::endl;
			return -101;
		}

		char *target_caln = new char[target_aln.length()+1];
		char *musashi_caln = new char[musashi_aln.length()+1];
		std::strcpy(target_caln, target_aln.c_str());
		std::strcpy(musashi_caln, musashi_aln.c_str());

		int pos1=-1,pos2=-1,pos3=-1;
		int pos=1;
		int tpos=1;
		for ( int i=0; i< (int)musashi_aln.length(); ++i ) {
			if ( pos == 4 && pos1== -1 ) {
				if ( musashi_caln[i]!='F' ) {
					std::cerr<<"Error: alignment should be a F, not an "<<musashi_caln[i]<<std::endl;
					return -105;
				}
				pos1=tpos;
				if ( target_caln[i]=='-' ) {
					std::cerr<<"Error: alignment at B1-2 is a gap"<<std::endl;
					return -108;
				}

			}

			if ( pos == 44 && pos2== -1 ) {
				if ( musashi_caln[i]!='F' ) {
					std::cerr<<"Error: alignment should be a F, not an "<<musashi_caln[i]<<std::endl;
					return -106;
				}
				pos2=tpos;
				if ( target_caln[i]=='-' ) {
					std::cerr<<"Error: alignment at B3-3 is a gap"<<std::endl;
					return -109;
				}
			}

			if ( pos == 46 && pos3== -1 ) {
				if ( musashi_caln[i]!='F' ) {
					std::cerr<<"Error: alignment should be a F, not an "<<musashi_caln[i]<<std::endl;
					return -107;
				}
				pos3=tpos;
				if ( target_caln[i]=='-' ) {
					std::cerr<<"Error: alignment at B3-5 is a gap"<<std::endl;
					return -110;
				}
			}

			if ( musashi_caln[i]!='-' ) {
				pos++;
			}
			if ( target_caln[i]!='-' ) {
				tpos++;
			}
		}

		int res1=-1,res2=-1,res3=-1;
		pos=1;
		for ( int i=1; i <= (int)pose.size(); i++ ) {
			if ( pose.residue(i).is_protein() ) {
				if ( pos==pos1 && res1==-1 ) {
					res1=i;
				}
				if ( pos==pos2 && res2==-1 ) {
					res2=i;
				}
				if ( pos==pos3 && res3==-1 ) {
					res3=i;
				}
				pos++;
			}
		}

		std::cout<<"Residue 1 "<<pose.pdb_info()->chain(res1) << pose.pdb_info()->number(res1)<<std::endl;
		std::cout<<"Residue 2 "<<pose.pdb_info()->chain(res2) << pose.pdb_info()->number(res2)<<std::endl;
		std::cout<<"Residue 3 "<<pose.pdb_info()->chain(res3) << pose.pdb_info()->number(res3)<<std::endl;

		scoring::ScoreFunctionOP scorefxn( get_score_function() );
		(*scorefxn)(pose);

		Size rpos1=0, rpos2=0;
		Real score1=1, score2=1;

		EnergyGraph & energy_graph(pose.energies().energy_graph());
		for ( int i=1; i <= (int)pose.size(); i++ ) {

			if ( pose.residue(i).is_RNA()||pose.residue(i).is_DNA() ) {

				Size cpos=0;
				Real cscore=1;
				for ( utility::graph::Graph::EdgeListIter
						iru  = energy_graph.get_node( i )->lower_edge_list_begin(),
						irue = energy_graph.get_node( i )->lower_edge_list_end();
						iru != irue; ++iru ) {
					EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
					Size const j( edge->get_first_node_ind() );
					if ( !((int)j==res1 || (int)j==res2 || (int)j==res3) ) {
						continue;
					}

					// the pair energies cached in the link
					EnergyMap const & emap( edge->fill_energy_map());
					Real const attr( emap[ fa_atr ] );
					//TR<<"\n"<<j<<": "<<attr<<"\n";
					if ( attr < -1. ) {
						if ( pose.residue(j).is_protein() ) {
							if ( attr < cscore ) {
								cpos=i;
								cscore=attr;
							}
							std::cout<<"Score "<<i<<" "<<pose.pdb_info()->chain(i) << pose.pdb_info()->number(i)<<" "<<j<<" "<<pose.pdb_info()->chain(j) << pose.pdb_info()->number(j)<<" "<<attr<<std::endl;
						}
					}
				}
				if ( cscore < score1 ) {
					rpos2=rpos1;
					score2=score1;
					rpos1=cpos;
					score1=cscore;
				} else if ( cscore<score2 ) {
					rpos2=cpos;
					score2=cscore;
				}
			}
		}

		for ( int i=1; i <= (int)pose.size(); i++ ) {
			Size cpos=0;
			Real cscore=1;
			//std::cout<<i<<" "<<res1<<" "<<res2<<" "<<res3<<std::endl;
			if ( pose.residue(i).is_protein() ) {
				if ( !(i == res1 || i == res2 ||i == res3) ) {
					continue;
				}

				for ( utility::graph::Graph::EdgeListIter
						iru  = energy_graph.get_node( i )->lower_edge_list_begin(),
						irue = energy_graph.get_node( i )->lower_edge_list_end();
						iru != irue; ++iru ) {
					EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
					Size const j( edge->get_first_node_ind() );
					if ( !(pose.residue(j).is_RNA() || pose.residue(j).is_DNA()) ) {
						continue;
					}

					// the pair energies cached in the link
					EnergyMap const & emap( edge->fill_energy_map());
					Real const attr( emap[ fa_atr ] );
					//TR<<"\n"<<j<<": "<<attr<<"\n";
					if ( attr < -1. ) {
						if ( attr < cscore ) {
							cpos=j;
							cscore=attr;
						}
						std::cout<<"Score "<<i<<" "<<pose.pdb_info()->chain(i) << pose.pdb_info()->number(i)<<" "<<j<<" "<<pose.pdb_info()->chain(j) << pose.pdb_info()->number(j)<<" "<<attr<<std::endl;
					}
				}
			}
			if ( cscore < score1 ) {
				rpos2=rpos1;
				score2=score1;
				rpos1=cpos;
				score1=cscore;
			} else if ( cscore<score2 ) {
				rpos2=cpos;
				score2=cscore;
			}
		}

		if ( rpos1>0 || rpos2>0 ) {
			if ( rpos1<rpos2 ) {
				if ( rpos1>1 ) {
					pose.conformation().delete_residue_range_slow(1,rpos1-1);
				}
				if ( rpos2>rpos1+1 ) {
					pose.conformation().delete_residue_range_slow(2, rpos2-rpos1);
				}
				if ( pose.size()>2 ) {
					pose.conformation().delete_residue_range_slow(3, pose.size());
				}
			} else {
				if ( rpos2>1 ) {
					pose.conformation().delete_residue_range_slow(1,rpos2-1);
				}
				if ( rpos1>rpos2+1 ) {
					pose.conformation().delete_residue_range_slow(2, rpos1-rpos2);
				}
				if ( pose.size()>2 ) {
					pose.conformation().delete_residue_range_slow(3, pose.size());
				}

			}
			if ( pose.size()>0 ) {
				std::stringstream dnpfilename;
				if ( pose.size()>1 ) {
					dnpfilename<<"dinucleotide_pair.pdb";
				} else {
					dnpfilename<<"single_nucleotide.pdb";
				}
				std::filebuf fb;
				fb.open(dnpfilename.str().c_str(), std::ios::out);
				std::ostream os(&fb);
				//pose.dump_pdb("test.pdb");

				Size acount=1;
				for ( Size i=1; i<= pose.size(); ++i ) {
					core::conformation::Residue const & rsd( pose.conformation().residue(i) );
					for ( Size j=1; j<= rsd.natoms(); ++j ) {
						if ( rsd.atom_name(j).compare(" P  ")&&rsd.atom_name(j).compare(" OP2")&&rsd.atom_name(j).compare(" OP1")&&rsd.atom_name(j).compare(" O5'")&&rsd.atom_name(j).compare(" C5'")&&rsd.atom_name(j).compare(" C4'")&&rsd.atom_name(j).compare(" O4'")&&rsd.atom_name(j).compare(" C3'")&&rsd.atom_name(j).compare(" O3'")&&rsd.atom_name(j).compare(" C2'")&&rsd.atom_name(j).compare(" H5'")&&rsd.atom_name(j).compare("H5''")&&rsd.atom_name(j).compare(" H4'")&&rsd.atom_name(j).compare(" H3'")&&rsd.atom_name(j).compare("H2''")&&rsd.atom_name(j).compare(" H2'")&&rsd.atom_name(j).compare("HO2'")&&rsd.atom_name(j).compare(" H1'")&&rsd.atom_name(j).compare(" O2'") ) {
							os<<"ATOM  "<<std::setw(5)<<acount;
							os<<" "<<rsd.atom_name(j)<<rsd.name3()<<"  A"<<std::setw(4)<<i;
							os<<"    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<rsd.atom(j).xyz()(1)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<rsd.atom(j).xyz()(2)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<rsd.atom(j).xyz()(3);
							os<<"  1.00  0.00           "<<rsd.atom_type(j).element()<<std::endl;
							acount++;
						}
					}
				}
				fb.close();
			}

		}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}



