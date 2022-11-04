// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
// /// @file apps/public/cl_complex_rescore
// /// @brief Application to rescore poses from RosettaDock for use with differential CL. In a manuscript to be published in 2022
// /// @author Zachary Drake (drake.463@osu.edu)
//
// // In Covalent Labeling (CL) Mass Spectrometry, proteins are exposed to labeling agents such as
// // hydroxyl radicals or DEPC which will label particular residues depending on solvent
// // accessbility and reactivity. The degree of modifiation of particular residues can be
// // assessed using Mass Spectrometry. When comparing the modification of labeled residues
// // between the unbound and bound forms of complex, the proximity of labeled residues to the
// // binding interface of the complex can be determined. This application uses user input
// // differential CL data to assign score to a given model.
//
// // To use this application, input structures from RosettaDock and CL Modification.

#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/option.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <utility/io/ozstream.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <iostream>
#include <core/pose/chains_util.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>



static basic::Tracer TR( "apps.pilot.zdrake.cl_complex_rescore" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

//local options
basic::options::StringOptionKey const interface( "interface" );
basic::options::StringOptionKey const input_file( "cl_file" );


void read_in_cl_data(
	utility::vector1< std::pair< core::Size, std::pair< core::Real, utility::vector1<char> > > > & input_cl_data,
	std::string const & input_fil
){
	utility::io::izstream input(input_fil);
	if ( !input.good() ) {
		std::string const msg( "Error opening file: " + input_fil);
		utility_exit_with_message( msg );
	}
	std::string line;
	//Import residue information from cl data file
	while ( getline(input,line) ) {
		if ( line.substr(0,1) == "#" ) continue;
		std::istringstream ss(line);
		core::Size resi;
		core::Real mon_rate;
		core::Real com_rate;
		std::string chains;

		ss >> resi >> mon_rate >> com_rate >> chains;

		//Calculate modification percentage from monomer and complex rates

		core::Real mod_per = 100.0 * ((mon_rate - com_rate)/(mon_rate));

		runtime_assert_string_msg( !(ss.fail() || ss.bad()), "Error in cl_complex_rescore::cl_file Could not parse line\"" + line + "\".");

		utility::vector1< char > temp_v(chains.begin(), chains.end());

		input_cl_data.push_back(make_pair(resi, make_pair( mod_per, temp_v) ) );
	}
}
void read_in_pdbs( utility::vector1<core::pose::PoseOP> & poses) {
	core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
	//Import each pose
	while ( input.has_another_pose() ) {
		core::pose::PoseOP mypose( new core::pose::Pose);
		input.fill_pose( *mypose );
		poses.push_back(mypose);
	}
	TR << "Number of poses entered: " << poses.size() << std::endl;
}

void cl_score(
	utility::vector1<core::pose::PoseOP> const & poses,
	utility::vector1< std::pair < core::Size, std::pair < core::Real, utility::vector1<char> > > > const & input_cl_data,
	utility::vector1< std::pair <std::string, core::Real >> & model_penalties,
	utility::vector1< std::pair < std::string, core::Real >> & cl_scores,
	std::string const & partner1,
	std::string const & partner2
){
	utility::vector1<char> partner1_chains(partner1.begin(), partner1.end());
	utility::vector1<char> partner2_chains(partner2.begin(), partner2.end());
	for ( auto const & pose : poses ) {
		//Iterate through each pose
		core::pose::Pose mypose = *pose;
		core::pose::PDBInfo const & pdb_info = *(mypose.pdb_info());
		utility::vector1<std::pair < core::Size, core::Real > > interface_distances;
		for ( core::Size lres = 1; lres <= input_cl_data.size(); ++lres ) {
			//Iterate through each labeled residue to find minimum interface distance
			core::Real min_dist = 5000.0;
			for ( core::Size lchain =1; lchain <= input_cl_data[lres].second.second.size(); ++lchain ) {
				if ( std::find(partner1_chains.begin(), partner1_chains.end(), input_cl_data[lres].second.second[lchain]) != partner1_chains.end() ) {
					utility::vector1<core::Size> lchain_res = core::pose::get_resnums_for_chain(mypose, input_cl_data[lres].second.second[lchain]);
					for ( core::Size res1 = 1; res1 <= lchain_res.size(); ++res1 ) {
						core::Size first_res_id = lchain_res[1]-1;
						if ( input_cl_data[lres].first == mypose.residue(res1).seqpos()-first_res_id && pdb_info.chain(res1) == input_cl_data[lres].second.second[lchain] ) {
							for ( core::Size res2 = 1; res2 <= mypose.size(); ++res2 ) {
								if ( pdb_info.chain(res2) != input_cl_data[lres].second.second[lchain] && std::find(partner2_chains.begin(), partner2_chains.end(), pdb_info.chain(res2)) != partner2_chains.end() ) {
									for ( core::Size latom = 1; latom <= mypose.residue(res1).nheavyatoms(); ++latom ) {
										if ( !mypose.residue(res1).atom_type(latom).is_virtual() ) {
											for ( core::Size atom = 1; atom <= mypose.residue(res2).nheavyatoms(); ++atom ) {
												if ( !mypose.residue(res2).atom_type(atom).is_virtual() ) {
													core::Real temp_dist = mypose.residue(res1).xyz(latom).distance(mypose.residue(res2).xyz(atom));
													if ( temp_dist < min_dist ) {
														min_dist = temp_dist;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
				if ( std::find(partner2_chains.begin(), partner2_chains.end(), input_cl_data[lres].second.second[lchain]) != partner2_chains.end() ) {
					utility::vector1<core::Size> lchain_res = core::pose::get_resnums_for_chain(mypose, input_cl_data[lres].second.second[lchain]);
					for ( core::Size res1 = 1; res1 <= lchain_res.size(); ++res1 ) {
						core::Size first_res_id = lchain_res[1]-1;
						if ( input_cl_data[lres].first == mypose.residue(res1).seqpos()-first_res_id && pdb_info.chain(res1) == input_cl_data[lres].second.second[lchain] ) {
							for ( core::Size res2 = 1; res2 <= mypose.size(); ++res2 ) {
								if ( pdb_info.chain(res2) != input_cl_data[lres].second.second[lchain] && std::find(partner1_chains.begin(), partner1_chains.end(), pdb_info.chain(res2)) != partner1_chains.end() ) {
									for ( core::Size latom = 1; latom <= mypose.residue(res1).nheavyatoms(); ++latom ) {
										if ( !mypose.residue(res1).atom_type(latom).is_virtual() ) {
											for ( core::Size atom = 1; atom <= mypose.residue(res2).nheavyatoms(); ++atom ) {
												if ( !mypose.residue(res2).atom_type(atom).is_virtual() ) {
													core::Real temp_dist = mypose.residue(res1).xyz(latom).distance(mypose.residue(res2).xyz(atom));
													if ( temp_dist < min_dist ) {
														min_dist = temp_dist;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			interface_distances.push_back(std::pair < core::Size, core::Real>(input_cl_data[lres].first, min_dist));
		}

		//The slope and intercept values were obtained from the interface distance and modification percentage relationship described in Drake 2022 manuscript. The linear parameters are used to predict modification percentage of model residues given their interface distances.
		core::Real slope = -2.06661005172;
		core::Real intercept = 46.2673485645;
		//A and B are functional parameters which have been otpimized for the sigmoidal penalty function described in Drake 2022 manuscript.
		core::Real A = 1.88;
		core::Real B = 38.0;
		core::Real model_penalty = 0.0;
		//calculate predicted modification change and model penalty
		for ( core::Size j = 1; j <= interface_distances.size(); ++j ) {
			if ( interface_distances[j].first == input_cl_data[j].first ) {
				core::Real temp_pred_mod = slope*interface_distances[j].second + intercept;
				core::Real temp_difference = fabs(input_cl_data[j].second.first - temp_pred_mod);
				model_penalty += 1.0 - ((1.0)/(1.0 + std::exp(A*(temp_difference-B))));

			}
		}
		model_penalties.push_back(std::pair < std::string, core::Real>(pdb_info.name(),model_penalty));
	}
	core::Real max_penalty = 0.0;
	//determine max penalty across all models
	for ( core::Size j = 1; j <= model_penalties.size(); ++j ) {
		if ( model_penalties[j].second > max_penalty ) {
			max_penalty = model_penalties[j].second;
		}
	}
	//normalize each model penalty with max penalty
	for ( core::Size j = 1; j <= model_penalties.size(); ++j ) {
		core::Real cl_score = 65.0*(model_penalties[j].second/max_penalty);
		cl_scores.push_back( std::pair < std::string, core::Real>(model_penalties[j].first,cl_score));
	}
}

void output_results (
	utility::vector1< std::pair < std::string, core::Real >> const & model_penalties,
	utility::vector1< std::pair < std::string, core::Real >> const & cl_scores
) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string outfile = "cl_score.out";
	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o]();
	utility::io::ozstream outz( outfile.c_str() );
	outz << "Model\t\t\tRaw_Penalty\t\tWeighted_CL_ScoreTerm" << std::endl;
	for ( core::Size i=1; i <= cl_scores.size(); ++i ) {
		outz << cl_scores[i].first << "\t\t" << model_penalties[i].second << "\t\t" << cl_scores[i].second << std::endl;
	}
	outz.close();
	outz.clear();
}

int
main(int argc, char * argv [])
{

	try{

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring;
		using utility::file::FileName;
		option.add( interface, "interface definition" ).def("A_B");
		option.add( input_file, "cl_file" ).def("");
		devel::init( argc, argv );

		//Read in variables from command line
		std::string intf( option[ interface ].value() );
		std::string input_fil( option[ input_file ].value() );

		if ( ! intf.find('_') ) {
			utility_exit_with_message("Unrecognized interface: " + intf + " must have side1 and side2, ex: A_BC or A_B");
		}
		utility::vector1<std::string> partner_chains = utility::string_split(intf, '_');
		std::string dock_partner1 = partner_chains[1];
		std::string dock_partner2 = partner_chains[2];

		TR << "Interface: " << intf << std::endl;
		TR << "Partner 1: " << dock_partner1 << std::endl;
		TR << "Partner 2: " << dock_partner2 << std::endl;
		//Import CL data
		utility::vector1< std::pair< core::Size, std::pair < core::Real, utility::vector1<char> > > > input_cl_data;
		read_in_cl_data(input_cl_data,input_fil);

		//Import Poses
		utility::vector1<core::pose::PoseOP> poses;
		read_in_pdbs( poses );

		//Claculate scores for each pose
		utility::vector1< std::pair< std::string, core::Real >> cl_scores;
		utility::vector1< std::pair< std::string, core::Real >> model_penalties;
		cl_score( poses, input_cl_data, model_penalties, cl_scores, dock_partner1, dock_partner2);

		//Output results
		output_results(model_penalties, cl_scores);

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}

