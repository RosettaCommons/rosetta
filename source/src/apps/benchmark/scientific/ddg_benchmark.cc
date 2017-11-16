// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Liz Kellogg ekellogg@u.washington.edu

#include <core/types.hh>

#include <core/chemical/AA.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>

#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <core/pack/task/ResfileReader.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>

// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector0.hh>


using basic::Error;
using basic::Warning;


using namespace core;
using namespace basic;
using namespace scoring;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

typedef std::vector<Real> ddGs;
static basic::Tracer TR( "apps.pilot.yiliu.ddg" );

///////////////////////////////////////////////////////////////////////////////
// YAML helper function
std::ostream & writeYamlValue(std::ostream & S, std::string name, core::Real value)
{
	S << "'" << name << "' : " << value << ", ";
	return S;
}

///////////////////////////////////////////////////////////////////////////////
// YAML helper function
std::ostream & writeYamlValue(std::ostream & S, std::string name, bool value)
{
	std::string sv;
	if ( value ) sv = "True";
	else sv = "False";

	S << "'" << name << "' : " << sv << ", ";
	return S;
}


////////////////////////// ddG Functions ////////////////////////////////////////

Real
sum(ddGs &scores_to_sum)
{
	Real sum=0;
	for ( Size i =0; i<scores_to_sum.size(); i++ ) {
		sum+=scores_to_sum[i];
	}
	return sum;
}

Real
average( utility::vector1<Real> &scores_to_average)
{
	Real sum = 0;
	for ( Size i =1; i<=scores_to_average.size(); i++ ) {
		sum+=scores_to_average[i];
	}
	return (sum/scores_to_average.size());
}

Size
store_energies( ObjexxFCL::FArray2D< Real > &two_d_e_arrays,
	scoring::ScoreFunction &s,
	pose::Pose &p, Size next_index ,  Size size_to_expect)
{
	s(p); //score the pose
	//all this to determine how many non-zero weights there are
	Size num_score_components = 0;
	EnergyMap::const_iterator it = s.weights().begin();
	while ( it != (s.weights()).end() ) {
		++it;
		if ( *it != 0 ) {
			num_score_components++;
		}
	}


	two_d_e_arrays.dimension(num_score_components,size_to_expect);

	Size current_score_component=0;
	Size j =1;
	for ( EnergyMap::const_iterator i = (s.weights()).begin(); i != s.weights().end(); ++i ) {
		//get score component of pose, then store in next slot of two_d_e_arrays

		current_score_component++;
		if ( *i != 0 ) {
			two_d_e_arrays(j++,next_index)=(*i) *
				((p.energies()).total_energies())[ScoreType(current_score_component)];
		}
	}
	return 0;
}

Size
average_score_components( ObjexxFCL::FArray2D< Real > &scores_to_average,
	utility::vector1<Real> &averaged_scores )
{
	averaged_scores = utility::vector1<Real>(scores_to_average.u1()-scores_to_average.l1()+1);
	for ( Size i = Size(scores_to_average.l1()); i <= Size(scores_to_average.u1()); i++ ) {
		Real sum_score_component = 0;
		for ( Size j = Size(scores_to_average.l2()); j <= Size(scores_to_average.u2()); j++ ) {
			sum_score_component += scores_to_average(i,j);
		}
		averaged_scores[i]=
			sum_score_component/(scores_to_average.u2()-scores_to_average.l2()+1);
	}
	return 0;
}

////////////////// Correlation Functions /////////////////////////////////
//obtain the correlation coefficient for two array.
Real correlation_coefficient(utility::vector1<Real> sim, utility::vector1<Real> exp){
	if ( sim.size()!=exp.size() ) {
		TR << "size does not match" << std::endl;
		//error, quit
		exit(1);
	}
	Real xx=0, yy=0, xy=0, x=0,y=0;
	//assuming the array start from 1 as in default of C/C++
	for ( Size i=1; i<=sim.size(); i++ ) {
		TR << "sim[i]: " << sim[i] << std::endl;
		TR << "exp[i]: " << exp[i] << std::endl;
		xx += sim[i]*sim[i];
		yy += exp[i]*exp[i];
		xy += sim[i]*exp[i];
		x+=sim[i];
		y+=sim[i];
	}
	xx/=sim.size();
	yy/=sim.size();
	xy/=sim.size();
	x/=sim.size();
	y/=sim.size();
	xy = xy-x*y;
	xx = xx-x*x;
	yy = yy-y*y;
	return xy/sqrt(xx*yy);
}

void correlation(std::string ddg_out){
	using namespace std;

	std::string const sim_filename( ddg_out.c_str() );
	std::string  const exp_filename( "./2lzm_ddg_exp.txt" );
	std::ifstream sim_data( sim_filename.c_str() );
	std::ifstream exp_data( exp_filename.c_str() );

	std::string sim_data_line,exp_data_line;
	char tmp[100],wt[100], res[100], mut[100], ddg[100];
	utility::vector1<Real> ddg_sim, ddg_exp;
	utility::vector1<std::string> mut_sim, mut_exp;
	std::string pos_info;

	/* get the ddg_exp vector */
	if ( !exp_data.good() ) {
		std::cout << "Unable to open exp data!" <<std::endl;
	} else {
		while ( getline(exp_data, exp_data_line) ) {
			sscanf(exp_data_line.c_str(),"%s %s %s %s %s", tmp, wt, res, mut, ddg);
			pos_info = string(wt) + string(" ") + string(res) + string(" ") + string(mut);
			mut_exp.push_back(pos_info);
			ddg_exp.push_back(atof(ddg));
		}
	}

	if ( !sim_data.good() ) {
		std::cout << "Unable to open exp data!" <<std::endl;
	} else {
		/* get the ddg_sim vector */
		while ( getline( sim_data,sim_data_line) ) {
			sscanf(sim_data_line.c_str(),"%s %s %s %s", wt, res, mut, ddg);
			pos_info = string(wt) + string(" ") + string(res) + string(" ") + string(mut);
			mut_sim.push_back(pos_info);
			ddg_sim.push_back(atof(ddg));
		}
	}

	Real correlation_result;
	correlation_result = correlation_coefficient(ddg_sim, ddg_exp);

	std::string results_fname( ".results.log" );
	std::ofstream staResult( results_fname.c_str() );
	if ( !staResult ) {
		TR.Error << "Can not open file " << results_fname;
	} else {
		char correlation_out[100];
		sprintf(correlation_out,"%0.2f",correlation_result);
		staResult << "The correlation coefficient for monomer ddG protocol on the t4-lysozyme test set: " << correlation_out << std::endl;
	}

	std::string yaml_fname( ".results.yaml" );
	std::ofstream yaml( yaml_fname.c_str() );
	if ( !yaml ) {
		TR.Error << "Can not open file " << yaml_fname;
	} else {
		yaml << "{ ";
		writeYamlValue(yaml, "CorrelationCoefficient", correlation_result);
		writeYamlValue(yaml, "_isTestPassed", correlation_result > 0.81 );
		yaml << "}\n";
	}
}


///////////////// Main ////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{

		using namespace pose;
		using namespace scoring;
		using namespace conformation;

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::pack::task;


		// setup random numbers and options
		devel::init(argc, argv);

		// read the pose
		pose::Pose pose;
		import_pose::pose_from_file( pose, basic::options::start_file() , core::import_pose::PDB_file); // gets filename from -s option

		// this is the numbering system relevant for the resfiles (currently... ie not pdb resnums)

		//store all energies of 20 repacked poses
		utility::vector1<Real> wt_averaged_score_components; // to be filled in and averaged later


		std::string weight_file = option[ OptionKeys::ddg::weight_file ]();
		ScoreFunctionOP score_structure_scorefxn(ScoreFunctionFactory::create_score_function(weight_file));

		Size num_iterations = option[ OptionKeys::ddg::iterations ]();
		//initialize output options.
		//debug output?
		bool debug_output = option[ OptionKeys::ddg::debug_output ]();
		if ( debug_output ) {
			std::cout << "weights being used: " <<
				score_structure_scorefxn->weights() << "\n";
		}

		//dump repacked pdbs?
		bool dump_pdbs = option[ OptionKeys::ddg::dump_pdbs ]();

		//output ddgs into what file?
		std::string ddg_out = option[ OptionKeys::ddg::out ]();
		std::ofstream ddg_output(ddg_out.c_str(), std::ios_base::app);
		if ( !ddg_output ) {
			std::cout << "having trouble opening output file for dumping predicted ddgs"
				<< ddg_out << std::endl;
			utility::exit(EXIT_FAILURE, __FILE__, __LINE__);
		}

		// initialize the scoring function stuff
		(*score_structure_scorefxn)(pose);
		/// Now handled automatically.  score_structure_scorefxn->accumulate_residue_total_energies(pose);

		//initialize vector for storing ddg energy components
		utility::vector1< ddGs> delta_delta_energy_components;
		utility::vector1< std::string > delta_delta_G_prefix;

		pack::task::PackerTaskOP storage_task(pack::task::TaskFactory::create_packer_task(pose));

		storage_task->initialize_from_command_line();
		pack::task::parse_resfile(pose, *storage_task);
		storage_task->or_include_current(true);

		//write out information to repack_native logfile
		//this eliminates unnecessary computations
		std::ostringstream native_logfile;
		native_logfile << "native.logfile";

		//initialize dG_wildtype which is final score
		Real dG_wildtype;

		std::ifstream fh(native_logfile.str().c_str(),std::ios::in);
		//check existence of file and for specific key-word
		bool dG_native_calculated=false;
		bool averaged_wt_scores_calculated=false;
		if ( fh.is_open() ) {
			if ( debug_output ) { std::cout <<
				"native logfile exists! checking for dG of native"
				<< std::endl;
			}
			std::string line;
			std::string result="";
			std::string averaged_score_components="";

			while ( fh >> line ) {
				//check for key-word dG_native:

				if ( line.compare("dG_wildtype:") == 0 ) {
					fh >> result;
				}
				if ( !result.empty() ) {
					if ( debug_output ) {
						std::cout << "result string is not empty.\n";
					}
					dG_native_calculated=true;
				}

				if ( line.compare("averaged_score_components:") == 0 ) {
					while ( fh >> averaged_score_components ) {
						//convert to double and store in wt_averaged_scores array
						std::istringstream convert_to_double(averaged_score_components);
						Real e_component;
						convert_to_double >> e_component;
						wt_averaged_score_components.push_back(e_component);
					}
					averaged_wt_scores_calculated=true;
				}

				if ( dG_native_calculated ) {
					if ( debug_output ) {
						std::cout << "converting to double\n";
					}
					std::istringstream convert_to_double(result);
					convert_to_double >> dG_wildtype;
					if ( debug_output ) {
						std::cout <<
							"final value saved to variable dG_wildtype is: "
							<< dG_wildtype << std::endl;
					}
					if ( dG_native_calculated && averaged_wt_scores_calculated ) {
						break;
					}
				}
			}
		}
		//store all energies of 20 repacked poses
		ObjexxFCL::FArray2D<Real> wt_score_components;

		if ( !dG_native_calculated ) {

			if ( debug_output ) {
				std::cout << "dG_wildtype hasn't been calculated yet! " <<
					"starting to calculate dG for native structure\n";
			}

			//start filehandler for recording the native structure logfile
			std::ofstream record_trajectories;
			record_trajectories.open(native_logfile.str().c_str());

			//start usual calculation of DG for native structure
			utility::vector1<Real> free_energy_wildtype(num_iterations,-999.999);

			//measure dG of input structure by repacking 20 times and taking the average score

			pack::task::PackerTaskOP repack_native(pack::task::TaskFactory::create_packer_task(pose));
			repack_native->restrict_to_repacking();
			for ( Size j =1; j<=pose.size(); j++ ) {
				//by default use ex1 and ex2
				repack_native->nonconst_residue_task(j).or_ex1(true);
				repack_native->nonconst_residue_task(j).or_ex2(true);
			}
			for ( Size i=1; i<=num_iterations; i++ ) { //repack for 20 cycles

				pose::Pose temporary_pose = pose;

				std::ostringstream q;
				q << i;
				//initialize packertask and prevent design at all positions but allow repacking.
				Real start_score_dG_wildtype( (*score_structure_scorefxn)( temporary_pose ) );
				if ( debug_output ) {
					record_trajectories << "round: "
						<< q.str() << " score before repacking wildtype "
						<< start_score_dG_wildtype << " \n"
						<< "start packing pose round: " << q.str() << std::endl;
				}//output
				pack::pack_rotamers(temporary_pose,(*score_structure_scorefxn),repack_native);
				if ( debug_output ) {
					record_trajectories << "end packing pose round: " << q.str() << std::endl;
				}//output

				//store final score
				free_energy_wildtype[i]=(*score_structure_scorefxn)(temporary_pose);
				Real final_score_dG_wildtype((*score_structure_scorefxn)(temporary_pose));

				//store components of energy breakdown in 2D array for easy averaging later on
				store_energies(wt_score_components, (*score_structure_scorefxn), temporary_pose,i,num_iterations);
				//end store components


				if ( debug_output ) {
					record_trajectories << "round: " << q.str() <<
						" score after repacking wildtype " <<
						final_score_dG_wildtype << ' ' <<
						temporary_pose.energies().total_energies().weighted_string_of( score_structure_scorefxn->weights() )
						<< std::endl;
				} //debug
				//output the repacked native for this round
				if ( dump_pdbs ) {
					std::string out_path;
					if ( option[out::path::pdb].active() ) {
						out_path = option[out::path::pdb].value();
					} else {
						out_path = option[out::path::pdb].default_value();
					}
					std::string dump_repacked_wt = out_path + "repacked_wt_round_" + q.str() + ".pdb";
					temporary_pose.dump_pdb(dump_repacked_wt);
				}
			}
			//average scores
			dG_wildtype = average(free_energy_wildtype);
			average_score_components(wt_score_components,wt_averaged_score_components);

			//output averages to file
			//===ddg_output << "WILDTYPE AVERAGED ENERGY COMPONENTS\n";
			//output energy component headers
			Size score_component =0;
			std::string header="";
			for ( EnergyMap::const_iterator i = (score_structure_scorefxn->weights()).begin();
					i != score_structure_scorefxn->weights().end(); ++i ) {
				score_component++;
				if ( *i != 0 ) {
					header = header + name_from_score_type(ScoreType(score_component)) +  " ";
				}
			}

			//record the dG of the repacked wild-type structure to save computer time later
			record_trajectories << "dG_wildtype: " << dG_wildtype << std::endl;
			record_trajectories << "averaged_score_components: ";
			//end output energy component headers
			record_trajectories << std::endl; //stored in case we start running another parallel job
			//end output averages to file

			record_trajectories.close(); //close the file after you're done.
		}
		//somehow create individual packertasks for each mutation.
		//residues to store information on which residues to mutate, and which amino acid to mutate to
		utility::vector1<bool> residues_to_mutate(pose.size(),false);
		utility::vector1<bool> residues_to_repack(pose.size(),true);
		utility::vector1<bool> aminoacids_to_mutate(20,false);

		//for each new point_mutant_task, create a new packertask object, copy the relevent information over
		//then create the mutant, score, and spit out.
		for ( Size i =1; i<=pose.size(); i++ ) {

			if ( storage_task->design_residue(i) ) {
				for ( ResidueLevelTask::ResidueTypeCOPListConstIter aa_iter(storage_task->residue_task(i).allowed_residue_types_begin()),
						aa_end(storage_task->residue_task(i).allowed_residue_types_end());
						aa_iter != aa_end; ++aa_iter ) {

					if ( debug_output ) {
						std::cout << " mutating residue " << i << " from " <<
							core::chemical::name_from_aa(pose.aa(i)) << " to amino acid: " <<
							core::chemical::name_from_aa((*aa_iter)->aa()) << std::endl;
					}
					char mutant_aa = core::chemical::oneletter_code_from_aa((*aa_iter)->aa());
					char wildtype_aa = core::chemical::oneletter_code_from_aa(pose.aa(i));
					std::string str_wildtype_aa;
					str_wildtype_aa.assign(&wildtype_aa,1);
					std::string str_mutant_aa;
					str_mutant_aa.assign(&mutant_aa,1);

					// create a log file for each structure for each mutation
					// only create the file if not already existing, otherwise die
					std::ofstream record_mutant_trajectories;
					std::ostringstream filename;
					filename << wildtype_aa << "_" << i << "_" << mutant_aa << ".logfile";

					std::ifstream fh( filename.str().c_str(), std::ios::in );
					if ( fh.is_open() ) {
						std::cout << "[DEBUG] file exists. "<< filename.str() << " skipping this mutation" << std::endl;
						continue;
					} else {
						std::cout << "[DEBUG] file does not exist. "<< filename.str() << ". Creating logfile now.\n";
						fh.close();
						record_mutant_trajectories.open( filename.str().c_str()  );
						if ( debug_output ) {
							record_mutant_trajectories << "beginning logfile" <<std::endl;
						}

					}
					//  end

					utility::vector1<Real> free_energy_mutants(num_iterations,-999.999);
					//storage for energy components of each of 20 repacks
					FArray2D<Real> mutant_energy_components = FArray2D<Real>();
					for ( Size k = 1; k <= num_iterations; k++ ) { //do this for 20 cycles
						//restrict the amino acids as specified by PIKAA
						//initialize new objects for the point-mutant
						utility::vector1<bool> restrict_to_aa = aminoacids_to_mutate;
						restrict_to_aa[(*aa_iter)->aa()]=true;

						//initialize packer and copy pose
						pose::Pose temporary_pose = pose;
						pack::task::PackerTaskOP point_mutant_packer_task(pack::task::TaskFactory::create_packer_task(temporary_pose));
						utility::vector1<bool> point_mutant = residues_to_mutate;
						point_mutant[i]=true;

						//debug output
						Real start_score_dG_mutant( (*score_structure_scorefxn)( temporary_pose ) ); //debug ek

						if ( debug_output ) {
							record_mutant_trajectories << "round: " <<
								k << " score before mutation " <<
								start_score_dG_mutant << ' ' << std::endl;
						}//debug

						//restrict to repacking for everything but the mutant
						for ( Size j=1; j <=temporary_pose.size(); j++ ) {
							if ( j!=i ) {
								point_mutant_packer_task->nonconst_residue_task(j).restrict_to_repacking();
								//always increase levels of rotamer sampling to ex1 ex2
								point_mutant_packer_task->nonconst_residue_task(j).or_ex1(true);
								point_mutant_packer_task->nonconst_residue_task(j).or_ex2(true);
							}
						}

						//restrict to the amino acid specified in the resfile for the residue-to-mutate
						point_mutant_packer_task->nonconst_residue_task(i).restrict_absent_canonical_aas(restrict_to_aa);
						point_mutant_packer_task->nonconst_residue_task(i).or_ex1(true);
						point_mutant_packer_task->nonconst_residue_task(i).or_ex2(true);
						//then do the mutation
						if ( debug_output ) {
							record_mutant_trajectories <<
								"start packing pose round: " << k <<
								std::endl;
						}

						pack::pack_rotamers(temporary_pose,(*score_structure_scorefxn),point_mutant_packer_task);

						if ( debug_output ) {
							record_mutant_trajectories <<
								"end packing pose round: " << k <<
								std::endl;
						}

						Real const final_score( (*score_structure_scorefxn)( temporary_pose ) );
						//store scores
						store_energies(mutant_energy_components, (*score_structure_scorefxn),temporary_pose, k,num_iterations);

						record_mutant_trajectories <<
							"score after mutation: residue " << wildtype_aa << " "
							<< i  << " " << mutant_aa << " "
							<< final_score << ' ' <<
							temporary_pose.energies().total_energies().weighted_string_of( score_structure_scorefxn->weights() ) << std::endl;

						//output the mutated pdb
						std::ostringstream q;
						q << k;
						std::ostringstream o;
						o << i;
						record_mutant_trajectories << "round " << q.str()
							<< " mutate " << wildtype_aa
							<< " " << o.str() << " " <<
							mutant_aa << " "  << (final_score) << std::endl;

						free_energy_mutants[k]=final_score;

						if ( dump_pdbs ) {
							std::string out_path;
							if ( option[out::path::pdb].active() ) {
								out_path = option[out::path::pdb].value();
							} else {
								out_path = option[out::path::pdb].default_value();
							}
							std::string output_pdb = out_path + "mut_" +
								str_wildtype_aa + "_" + o.str() + "_" +
								str_mutant_aa + "round_" + q.str() + ".pdb";
							temporary_pose.dump_pdb(output_pdb);
						}
					}

					std::ostringstream resnum;
					resnum << i;
					//average scores
					Real end_dG = average(free_energy_mutants);
					utility::vector1<Real> averaged_score_components;
					average_score_components(mutant_energy_components,averaged_score_components);

					//output averages to file
					//end output averages to file
					//store differences in energy
					ddGs delta_delta_G;
					for ( Size k =1; k<= averaged_score_components.size(); k++ ) {
						if ( debug_output ) {
							std::cout << "mutant score component at index " << k
								<< " is: " << averaged_score_components[k]
								<< " " << "wt score component at index " << k
								<< " is: " << wt_averaged_score_components[k] << std::endl;
						}
						Real delta_e = averaged_score_components[k]-wt_averaged_score_components[k];

						delta_delta_G.push_back(delta_e);
						//would have been easier if i just did it with FArrays
					}

					delta_delta_energy_components.push_back( delta_delta_G );
					delta_delta_G_prefix.push_back(str_wildtype_aa + " " + resnum.str() + " " + str_mutant_aa);
					std::ostringstream r;
					r << i;
					record_mutant_trajectories << "mutate " << wildtype_aa <<
						" " << r.str() << " " << mutant_aa  << " wildtype_dG is: "
						<< dG_wildtype << " and mutant_dG is: "
						<< end_dG << " ddG is: " << (end_dG-dG_wildtype)
						<< std::endl;
					record_mutant_trajectories.close();
				}
			}
		}

		//now output delta energy components for all mutations
		for ( Size c=1; c<=delta_delta_G_prefix.size(); c++ ) {
			ddg_output << delta_delta_G_prefix[c] << " ";
			ddGs ddg_score_components = delta_delta_energy_components[c];
			ddg_output << F(9,3,sum(ddg_score_components)) << " ";
			ddg_output << std::endl;
		}
		//end
		std::string const sim_filename( ddg_out.c_str() );
		std::string  const exp_filename( "2lzm_ddg_exp.txt" );
		std::ifstream sim_data( sim_filename.c_str() );
		std::ifstream exp_data( exp_filename.c_str() );
		correlation(ddg_out);
		exit(0);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

