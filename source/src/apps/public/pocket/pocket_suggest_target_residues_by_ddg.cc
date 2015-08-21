// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk + dj

#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/scoring/rms_util.hh>

// Numeric Headers
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

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

OPT_KEY( String, ddg_list )
OPT_KEY( String, target_chain_list )
//OPT_KEY( Integer, num_angles )

static basic::Tracer TR( "apps.pilot.david_choose_target_residues.main" );
core::Real interface_residue_ddg (core::pose::Pose const & pose, core::Size resno);
bool is_interface_residue (const char chain, const int resno);


//set to store pdb info keys
std::set <std::string> interface;
std::set <std::string> interface_resi;
std::set <char> chains;

//stores resid of the ligand residue

/// General testing code
int
main( int argc, char * argv [] )
{
	try {



		utility::vector1< std::string > targets;
		utility::vector1< core::Real > scores;
		int target_count=1;

		NEW_OPT ( ddg_list, "File name from Robetta DDG calculations","");
		NEW_OPT ( target_chain_list, "Comma separated list of chains to use for target residue selection","");
		//NEW_OPT ( num_angles, "Number of different pose angles to measure score at", 1);
		devel::init(argc, argv);
		pose::Pose input_pose;
		//int angles = option[ num_angles ];
		//if (angles <1){
		//  fprintf (stderr, "Error: invalid number of angles.  Must be greather than 0\n");
		//  return -1;
		//}
		//read in pdb file from command line
		std::string const input_pdb_name ( basic::options::start_file() );
		core::import_pose::pose_from_pdb( input_pose, input_pdb_name );
		core::Real ddg=0;
		core::Real dim = 12.;

		std::string const dfilename = option[ ddg_list ];
		std::string const chain_list = option[ target_chain_list ];
		if ( dfilename != "" ) {
			std::ifstream ifs(dfilename.c_str(), std::ifstream::in);
			if ( !ifs.is_open() ) {
				std::cout<< "Error opening ddg list file "<<dfilename<<std::endl;
				return -100;
			}
			//ifb.open (dfilename,std::ios::in);
			//std::ostream ios(&ifb);
			std::string intres;
			while ( ifs.good() ) {
				std::getline(ifs, intres);
				interface.insert(intres);
			}
		}
		char * token;
		char * str = new char[chain_list.length()+1];
		strcpy (str, chain_list.c_str());
		token = strtok (str, ",");
		while ( token != NULL ) {
			if ( strlen(token)>1 ) {
				std::cout << "Chains shouls only be one character" << std::endl;
				return -1;
			}
			chains.insert(token[0]);
			token = strtok (NULL, ",");
		}
		if ( chains.empty()/*size() == 0*/ ) {
			std::cout << "No target chains specified. use -target_chain_list A[,B,C...]"<<std::endl;
			return -1;
		}



		//std::filebuf fb,fb2;
		//std::stringstream filename, filename2;
		//if (!option[ OptionKeys::out::output_tag ]().empty()){
		//filename<<option[ OptionKeys::out::output_tag ]()<<".pscore";
		//filename2<<option[ OptionKeys::out::output_tag ]()<<".lpscore";
		//}else{
		//filename<<basic::options::start_file()<<".pscore";
		//filename2<<basic::options::start_file()<<".lpscore";
		//}
		//fb.open (filename.str().c_str(),std::ios::out);
		//fb2.open (filename2.str().c_str(),std::ios::out);
		//std::ostream os(&fb);
		//std::ostream os2(&fb2);


		scoring::ScoreFunctionOP scorefxn( get_score_function() );
		(*scorefxn)(input_pose);
		EnergyGraph & energy_graph(input_pose.energies().energy_graph());

		for ( int j = 1, resnum = input_pose.total_residue(); j <= resnum; ++j ) {
			if ( chains.count(input_pose.pdb_info()->chain(j)) ==0 ) {
				for ( graph::Graph::EdgeListIter
						iru  = energy_graph.get_node( j )->lower_edge_list_begin(),
						irue = energy_graph.get_node( j )->lower_edge_list_end();
						iru != irue; ++iru ) {
					EnergyEdge * edge( static_cast< EnergyEdge *> (*iru) );
					Size const tpos( edge->get_first_node_ind() );
					if ( chains.count(input_pose.pdb_info()->chain(tpos)) ==0 ) continue;
					// the pair energies cached in the link
					EnergyMap const & emap( edge->fill_energy_map());
					Real const attr( emap[ fa_atr ] );
					//TR<<"\n"<<j<<": "<<attr<<"\n";
					if ( attr < 0 ) {
						// create string id to store in set
						std::ostringstream residuestream;
						residuestream << input_pose.pdb_info()->chain(tpos) << input_pose.pdb_info()->number(tpos);
						std::string res_id = residuestream.str();
						if ( interface_resi.count(res_id) ==0 ) {
							interface_resi.insert(res_id);
						}
					}
				}
			}
		}

		for ( int j = 1, resnum = input_pose.total_residue(); j <= resnum; ++j ) {
			bool iface_ala = false;
			if ( input_pose.residue(j).name1() == 'A' ) {
				std::ostringstream residuestream;

				residuestream << input_pose.pdb_info()->chain(j) << input_pose.pdb_info()->number(j);
				std::string res_id = residuestream.str();
				if ( interface_resi.count(res_id) > 0 ) {
					iface_ala=true;
				}

			}


			if ( iface_ala || is_interface_residue(input_pose.pdb_info()->chain(j), input_pose.pdb_info()->number(j)) ) {
				core::Real tx1=0,ty1=0,tz1=0;
				int count1=0;


				core::conformation::Residue const & rsd1 = input_pose.residue(j);
				for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {
					if ( !rsd1.atom_is_backbone(i) && !rsd1.atom_is_hydrogen(i) ) {
						tx1+=rsd1.atom(i).xyz()(1);
						ty1+=rsd1.atom(i).xyz()(2);
						tz1+=rsd1.atom(i).xyz()(3);
						count1++;
					}
				}
				tx1+=rsd1.xyz("CA")(1);
				ty1+=rsd1.xyz("CA")(2);
				tz1+=rsd1.xyz("CA")(3);
				count1++;
				tx1/=count1;
				ty1/=count1;
				tz1/=count1;

				for ( int k = j+1, resnum = input_pose.total_residue(); k <= resnum; ++k ) {
					iface_ala = false;
					if ( input_pose.residue(k).name1() == 'A' ) {
						std::ostringstream residuestream;

						residuestream << input_pose.pdb_info()->chain(k) << input_pose.pdb_info()->number(k);
						std::string res_id = residuestream.str();
						if ( interface_resi.count(res_id) > 0 ) {
							iface_ala=true;
						}

					}

					if ( iface_ala || is_interface_residue(input_pose.pdb_info()->chain(k), input_pose.pdb_info()->number(k)) ) {
						core::Real tx2=0,ty2=0,tz2=0;
						int count2=0;
						core::conformation::Residue const & rsd2 = input_pose.residue(k);
						for ( Size i = 1, i_end = rsd2.nheavyatoms(); i <= i_end; ++i ) {
							if ( !rsd2.atom_is_backbone(i) && !rsd2.atom_is_hydrogen(i) ) {
								tx2+=rsd2.atom(i).xyz()(1);
								ty2+=rsd2.atom(i).xyz()(2);
								tz2+=rsd2.atom(i).xyz()(3);
								count2++;
							}
						}
						tx2+=rsd2.xyz("CA")(1);
						ty2+=rsd2.xyz("CA")(2);
						tz2+=rsd2.xyz("CA")(3);
						count2++;
						tx2/=count2;
						ty2/=count2;
						tz2/=count2;

						double ires_dist = sqrt (pow (tx1 - tx2, 2) + pow (ty1 - ty2, 2) + pow (tz1-tz2,2));
						//if ( (input_pose.pdb_info()->number(j) == 100 || input_pose.pdb_info()->number(k) == 100) && (input_pose.pdb_info()->number(j) == 146 || input_pose.pdb_info()->number(k) == 146)){
						//std::cout<<ires_dist<<"\t"<<tx1<<" "<<ty1<<" "<<tz1<<" "<<tx2<<" "<<ty2<<" "<<tz2<<std::endl;
						//}
						//double ires_dist = sqrt( pow(input_pose.residue(j).xyz("CA")(1)-input_pose.residue(k).xyz("CA")(1),2) + pow(input_pose.residue(j).xyz("CA")(2)-input_pose.residue(k).xyz("CA")(2),2) + pow(input_pose.residue(j).xyz("CA")(3)-input_pose.residue(k).xyz("CA")(3),2) );
						if ( ires_dist > 15. ) continue;
						if ( input_pose.pdb_info()->number(j) - input_pose.pdb_info()->number(k)<6 && input_pose.pdb_info()->number(k) - input_pose.pdb_info()->number(j)<6 ) continue;


						core::Real tx=0,ty=0,tz =0;
						core::Real cminX=tx1,cminY=ty1,cminZ=tz1,cmaxX=tx1,cmaxY=ty1,cmaxZ=tz1;
						if ( cminX > tx2 ) cminX = tx2;
						if ( cminY > ty2 ) cminY = ty2;
						if ( cminZ > tz2 ) cminZ = tz2;
						if ( cmaxX < tx2 ) cmaxX = tx2;
						if ( cmaxY < ty2 ) cmaxY = ty2;
						if ( cmaxZ < tz2 ) cmaxZ = tz2;
						// AMW: cppcheck notes the immediate reassignment...
						// Commenting out the meaningless statements and moving declaration down to make the relationship more clear.
						//core::Real dimX = dim - (cmaxX-cminX);
						//core::Real dimY = dim - (cmaxY-cminY);
						//core::Real dimZ = dim - (cmaxZ-cminZ);
						core::Real dimX = dim;
						core::Real dimY = dim;
						core::Real dimZ = dim;
						//if (dimX<4.) dimX=4.;
						//if (dimY<4.) dimY=4.;
						//if (dimZ<4.) dimZ=4.;



						int count=0;
						//core::conformation::Residue const & rsd1 = input_pose.residue(j);
						//core::conformation::Residue const & rsd2 = input_pose.residue(k);
						for ( Size i = 1, i_end = rsd1.nheavyatoms(); i <= i_end; ++i ) {
							if ( !rsd1.atom_is_backbone(i) && !rsd1.atom_is_hydrogen(i) ) {
								tx+=rsd1.atom(i).xyz()(1);
								ty+=rsd1.atom(i).xyz()(2);
								tz+=rsd1.atom(i).xyz()(3);
								count++;
							}
						}
						tx+=rsd1.xyz("CA")(1);
						ty+=rsd1.xyz("CA")(2);
						tz+=rsd1.xyz("CA")(3);
						count++;
						for ( Size i = 1, i_end = rsd2.nheavyatoms(); i <= i_end; ++i ) {
							if ( !rsd2.atom_is_backbone(i) && !rsd2.atom_is_hydrogen(i) ) {
								tx+=rsd2.atom(i).xyz()(1);
								ty+=rsd2.atom(i).xyz()(2);
								tz+=rsd2.atom(i).xyz()(3);
								count++;
							}
						}
						tx+=rsd2.xyz("CA")(1);
						ty+=rsd2.xyz("CA")(2);
						tz+=rsd2.xyz("CA")(3);
						count++;
						tx/=count;
						ty/=count;
						tz/=count;

						int gridResCount=0;
						core::Real gridDDG = 0.;
						for ( int i = 1, resnum = input_pose.total_residue(); i <= resnum; ++i ) {
							ddg=interface_residue_ddg(input_pose, i);
							core::Real minX=10000.;
							core::Real maxX=-10000.;
							core::Real minY=10000.;
							core::Real maxY=-10000.;
							core::Real minZ=10000.;
							core::Real maxZ=-10000.;
							if ( ddg>0 ) {
								core::conformation::Residue const & rsd3 = input_pose.residue(i);
								for ( Size a = 1, a_end = rsd3.nheavyatoms(); a <= a_end; ++a ) {
									if ( !rsd3.atom_is_backbone(a) && !rsd3.atom_is_hydrogen(a) ) {
										if ( minX > rsd3.atom(a).xyz()(1) ) minX = rsd3.atom(a).xyz()(1);
										if ( minY > rsd3.atom(a).xyz()(2) ) minY = rsd3.atom(a).xyz()(2);
										if ( minZ > rsd3.atom(a).xyz()(3) ) minZ = rsd3.atom(a).xyz()(3);
										if ( maxX < rsd3.atom(a).xyz()(1) ) maxX = rsd3.atom(a).xyz()(1);
										if ( maxY < rsd3.atom(a).xyz()(2) ) maxY = rsd3.atom(a).xyz()(2);
										if ( maxZ < rsd3.atom(a).xyz()(3) ) maxZ = rsd3.atom(a).xyz()(3);
									}
								}
								//if ( (input_pose.pdb_info()->number(j) == 97 || input_pose.pdb_info()->number(k) == 97) && (input_pose.pdb_info()->number(j) == 143 || input_pose.pdb_info()->number(k) == 143)){
								//std::cout<<input_pose.pdb_info()->number(i)<<" "<<minX<<" "<<minY<<" "<<minZ<<" - "<<maxX<<" "<<maxY<<" "<<maxZ<<" - "<<cminX-dimX<<" "<<cminY-dimY<<" "<<cminZ-dimZ<<" - "<<cmaxX+dimZ<<" "<<cmaxY+dimY<<" "<<cmaxZ+dimZ<<std::endl;
								//}

								if ( (minX > (cminX - dimX)) && (minY > (cminY - dimY)) && (minZ > (cminZ - dimZ)) && (maxX < (cmaxX + dimX)) && (maxY < (cmaxY + dimY)) && (maxZ < (cmaxZ + dimZ)) ) {
									//if ( (input_pose.pdb_info()->number(j) == 97 || input_pose.pdb_info()->number(k) == 97) && (input_pose.pdb_info()->number(j) == 143 || input_pose.pdb_info()->number(k) == 143)){
									//std::cout<<input_pose.pdb_info()->number(i)<<" "<<ddg<<std::endl;
									//}
									gridResCount++;
									gridDDG += ddg;
								}

							}
						}
						scores.push_back(gridDDG);
						std::stringstream tmpstream;
						tmpstream<<input_pose.pdb_info()->number(j)<<":"<<input_pose.pdb_info()->chain(j)<<","<<input_pose.pdb_info()->number(k)<<":"<<input_pose.pdb_info()->chain(k);
						targets.push_back(tmpstream.str());
						//scores[target_count]=gridDDG;
						//targets[target_count]<<input_pose.pdb_info()->number(j)<<":"<<input_pose.pdb_info()->chain(j)<<","<<input_pose.pdb_info()->number(k)<<":"<<input_pose.pdb_info()->chain(k);
						target_count++;
					}
				}
			}

		}

		for ( int j = target_count-1; j>1; --j ) {
			int minpos=1;
			core::Real minscore=scores[1];
			for ( int i=1; i<=j; ++i ) {
				if ( minscore>scores[i] ) {
					minscore=scores[i];
					minpos=i;
				}
			}
			if ( minpos ==j ) continue;
			core::Real tmpscore = scores[j];
			std::string tmptarget;
			tmptarget = targets[j];
			scores[j] = scores[minpos];
			targets[j] = targets[minpos];
			scores[minpos]=tmpscore;
			targets[minpos] = tmptarget;
		}
		std::ofstream outfile;
		outfile.open( "targets.txt");

		for ( int i =1; i < target_count; ++i ) {
			outfile<<scores[i]<<"\t"<<targets[i]<<std::endl;
		}
		outfile.close();


		//fb.close();
		//fb2.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}

core::Real interface_residue_ddg (core::pose::Pose const & pose, core::Size resno){
	core::Real out=0;
	for ( std::set<std::string>::iterator it=interface.begin(); it!=interface.end();
			++it ) {
		char * token;
		char * str = new char[it->length()+1];
		strcpy (str, it->c_str());
		token = strtok (str, " ");
		int token_count=0;
		while ( token != NULL ) {
			if ( token_count==0 ) {
				if ( pose.pdb_info()->number(resno) == atoi(token) ) {
					token_count++;
					token = strtok (NULL, " ");
				} else { break; }
			} else if ( token_count == 1 ) {
				if ( pose.pdb_info()->chain(resno) == token[0] ) {
					token_count++;
					token = strtok (NULL, " ");
				} else { break; }
			} else if ( token_count == 5 ) {
				out=(core::Real)atof(token);
				if ( out < 0 ) out=0;
				break;
			} else {
				token_count++;
				token = strtok (NULL, " ");
			}
		}
		delete[] str;
	}
	return out;
}

bool is_interface_residue (const char chain, const int resno){
	bool out=false;
	if ( chains.count(chain)==0 ) return out;
	for ( std::set<std::string>::iterator it=interface.begin(); it!=interface.end();
			++it ) {
		char * token;
		char * str = new char[it->length()+1];
		strcpy (str, it->c_str());
		token = strtok (str, " ");
		int token_count=0;
		while ( token != NULL ) {
			if ( token_count==0 ) {
				if ( resno == atoi(token) ) {
					token_count++;
					token = strtok (NULL, " ");
				} else { break; }
			} else if ( token_count == 1 ) {
				if ( chain == token[0] ) {
					delete[] str;
					return true;
				} else { break; }
			}
		}
		delete[] str;
	}
	return out;
}

