// -*- mode:c++;tab-width:2;incdent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/analysis/SID_ERMS_prediction
/// @brief Application to predict SID ERMS data
/// @author Justin Seffernick (seffernick.9@osu.edu)

#include <basic/options/option.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/chains_util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <utility/io/ozstream.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>
#include <numeric/random/uniform.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <cmath>

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <protocols/pose_metric_calculators/SaltBridgeCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricContainer.hh>

#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>


static basic::Tracer TR( "apps.public.analysis.SID_ERMS_prediction" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

//local options
basic::options::FileOptionKey const complex_type( "complex_type" );
basic::options::FileOptionKey const ERMS( "ERMS" );
basic::options::BooleanOptionKey const RMSE( "RMSE");
basic::options::FileOptionKey const B_vals( "B_vals" );
basic::options::RealOptionKey const steepness( "steepness" );
basic::options::RealOptionKey const breakage_cutoff( "breakage_cutoff" );

//calculate probability
core::Real calc_prob(core::Real x,core::Real a, core::Real b) {
	return -1.0 / (1.0 + exp(a*(x-b)) ) + 1;
}

//read in the complex type: subunits and connectivities (nodes and edges)
void read_complex_type(core::Size &n_chains, utility::vector1<char> &nodes, utility::vector1<utility::vector1<std::string>> &edges) {
	std::string complex_type_filename;
	complex_type_filename = option[ complex_type ]();
	utility::io::izstream input(complex_type_filename);

	//error if invalid file
	if ( !input.good() ) {
		std::string const msg( "Error opening file: " + complex_type_filename );
		utility_exit_with_message( msg );
	}

	//read file
	std::string line;
	core::Size count = 0;
	while ( getline(input,line) ) {
		std::istringstream ss(line);
		utility::vector1< std::string > inputs = utility::split_whitespace(ss.str());
		//read first line, which contains chain IDs
		if ( count==0 ) {
			n_chains = inputs.size();
			TR << "Complex type file has " << n_chains << " chains." << std::endl;
			for ( core::Size i=1; i<=n_chains; i++ ) {
				if ( inputs[i].size() != 1 ) {
					std::string const msg( "Chain input incorrectly: " + inputs[i]);
					utility_exit_with_message( msg );
				}
				char node_curr;
				node_curr = inputs[i][0];
				nodes.push_back(node_curr);
				TR << "Chain " << i << ": " << node_curr << std::endl;
			}
		} else { //read in subsequent lines, which contain the interfaces
			core::Size n_edges = inputs.size(); //(ss.str().length()+1)/4;
			if ( n_edges == 0 ) {
				std::string const msg( "Interface input incorrectly");
				utility_exit_with_message( msg );
			}
			utility::vector1<std::string> edges_curr;
			for ( core::Size i=1; i<=n_edges; i++ ) {
				std::string edge = inputs[i];
				//check edge format
				if ( edge.size()!=3 || edge[1]!='_' ) {
					std::string const msg( "Interface input incorrectly: " + edge);
					utility_exit_with_message( msg );
				}

				edges_curr.push_back(edge);
			}
			TR << "Interface type " << count << " has " << n_edges << " symmetric interfaces. Interface " << edges_curr[1] << " is used." << std::endl;
			edges.push_back(edges_curr);
		}
		count++;
	}
	//check to make sure all chains are involved in at least one interface
	for ( core::Size i=1; i<=nodes.size(); i++ ) {
		char node_curr = nodes[i];
		bool node_found(false);
		for ( core::Size j=1; j<=edges.size(); j++ ) {
			for ( core::Size k=1; k<=edges[j].size(); k++ ) {
				if ( node_curr==edges[j][k][0] || node_curr==edges[j][k][2] ) {
					node_found = true;
				}
			}
		}
		if ( !node_found ) {
			std::string msg( "Complex type file incorrect. Chain ");
			msg.push_back(node_curr);
			msg+=" not in an interface. All chains must participate in at least one interface";
			utility_exit_with_message( msg );
		}
	}
}


//make sure intensities sum to 1
void check_intensities( const utility::vector1<utility::vector1<core::Real>> &ERMS_data ) {
	for ( core::Size i=1; i<=ERMS_data.size(); i++ ) {
		core::Real sum = 0;
		for ( core::Size j=1; j<=ERMS_data[i].size(); j++ ) {
			sum += ERMS_data[i][j];
		}
		if ( sum < 0.95 || sum > 1.05 ) {
			std::string const msg( "Error in ERMS input. Intensities at each acceleration energy must add to 1.");
			utility_exit_with_message( msg );
		}
	}
}

//read in acceleration energies or full ERMS
void read_ERMS(utility::vector1<core::Real> &ACE, utility::vector1<utility::vector1<core::Real>> &ERMS_read, const core::Size n_chains) {

	bool input_ERMS = false;
	if ( option[ RMSE ].user() ) {
		input_ERMS = true;
	} else {
		TR.Warning << "Acceleration energies read in, but not ERMS values." << std::endl;
	}
	std::string ERMS_filename;
	ERMS_filename = option[ ERMS ]();
	utility::io::izstream input(ERMS_filename);

	//error if invalid file
	if ( !input.good() ) {
		std::string const msg( "Error opening file: " + ERMS_filename );
		utility_exit_with_message( msg );
	}

	//read file
	std::string line;
	while ( getline(input,line) ) {
		std::istringstream ss(line);
		core::Real ACE_current;

		ss >> ACE_current;

		if ( ACE_current < 0 ) {
			std::string const msg( "Error in acceleration energy input. All acceleration energies must be positive." );
			utility_exit_with_message( msg );
		}

		utility::vector1<core::Real> ERMS_current;
		//read in ERMS data
		if ( input_ERMS ) {
			for ( core::Size i=1; i<=n_chains; i++ ) {
				core::Real ERMS_val_current;
				ss >> ERMS_val_current;

				if ( ERMS_val_current < 0 or ERMS_val_current > 1 ) {
					std::string const msg( "Error in experimental intensity input. All intensities should be between 0 and 1." );
					utility_exit_with_message( msg );
				}

				ERMS_current.push_back(ERMS_val_current);
			}
			ERMS_read.push_back(ERMS_current);
		}

		runtime_assert_string_msg( !(ss.fail() || ss.bad()), "Error in SID_ERMS_prediction:  Could not parse line \"" + line + "\"." );
		ACE.push_back(ACE_current);
	}
	check_intensities(ERMS_read);
}

//read in B values if given via file input
void read_B_vals(utility::vector1<core::Real> &B) {
	std::string B_filename;
	B_filename = option[ B_vals ]();
	utility::io::izstream input(B_filename);

	//error if invalid file
	if ( !input.good() ) {
		std::string const msg( "Error opening file: " + B_filename );
		utility_exit_with_message( msg );
	}

	std::stringstream err_msg;
	err_msg << "Incorrect number of B values given, must match the number of interface types." << std::endl;

	//read file
	std::string line;
	core::Size count = 1;
	while ( getline(input,line) ) {
		std::istringstream ss(line);

		if ( count <= B.size() ) {
			ss >> B[count];
			count ++;
		} else {
			utility_exit_with_message(err_msg.str()); //too many B values input
		}
		runtime_assert_string_msg( !(ss.fail() || ss.bad()), "Error in SID_ERMS_prediction:  Could not parse line \"" + line + "\"." );
	}
	if ( count-1 != B.size() ) {
		utility_exit_with_message(err_msg.str()); //too few B values input
	}
}

//check to make sure interfaces are symmetric
void check_interface_symmetry(const core::pose::PoseOP &pose_check, const utility::vector1<std::string> &edges_check, core::Real &dSASA, core::Real &PRE) {
	core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();

	//calculate values for first interface
	std::string interface_1(edges_check[1]);
	protocols::analysis::InterfaceAnalyzerMoverOP IAM = utility::pointer::make_shared< protocols::analysis::InterfaceAnalyzerMover >(interface_1, true, core::scoring::ScoreFunctionFactory::create_score_function("ref2015")/*scorefxn_*/, false/*compute_packstat_*/, false/*pack_together_*/, false/*pack_separated_*/);
	IAM->apply(*pose_check);
	protocols::analysis::InterfaceData int_data = IAM->get_all_data();
	utility::vector1< core::Real > dSASA_vec = int_data.dSASA;
	core::Real dSASA_1 = dSASA_vec[1];//get_all_per_residue_data
	dSASA = dSASA_1; //set actual value

	protocols::analysis::PerResidueInterfaceData int_data_per = IAM->get_all_per_residue_data();
	utility::vector1< core::Real > PRE_vec = int_data_per.regional_avg_per_residue_energy_int;
	PRE = PRE_vec[1]; //set actual value

	//calculate size for remaining interfaces and compare to first. Should be within 10%.
	for ( core::Size j=2; j<=edges_check.size(); j++ ) {
		std::string interface_curr(edges_check[j]);
		protocols::analysis::InterfaceAnalyzerMoverOP IAM = utility::pointer::make_shared< protocols::analysis::InterfaceAnalyzerMover >(interface_curr, true, core::scoring::ScoreFunctionFactory::create_score_function("ref2015")/*scorefxn_*/, false/*compute_packstat_*/, false/*pack_together_*/, false/*pack_separated_*/);
		IAM->apply(*pose_check);
		protocols::analysis::InterfaceData int_data = IAM->get_all_data();
		utility::vector1< core::Real > dSASA_vec = int_data.dSASA;
		core::Real dSASA_curr = dSASA_vec[1];

		//warn if size difference is greater than 10% compared to first interface.
		core::Real percent_diff = std::abs(dSASA_1-dSASA_curr)/(dSASA_1)*100.0;
		if ( percent_diff > 10 ) {
			TR.Warning << "Interfaces " << interface_1 << " and " << interface_curr << " input as symmetric, but size varies by " << percent_diff << " %. Based on the PDB file, there is a high likelihood that the complex type file is incorrect. It is recommended that you adjust the complex type file to remove the symmetry."  << std::endl;
		}
	}
}

//calculate B values from PDB (and check for disulfide bond if dimer to shift steepness)
void calc_B_values(utility::vector1<core::Real>  &B, const utility::vector1<utility::vector1<std::string>> &edges, core::Real &A, const core::Size n_chains ) {
	core::Real w_SA(0.488748);
	core::Real w_PRE(-457.753);
	core::Real w_int(-1488.57);
	if ( option[ breakage_cutoff ].user() ) {
		w_SA = 0.464145;
		w_PRE = -589.803;
		w_int = -1834.75;
	}
	//read in PDB
	utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ]();
	if ( filenames.size() == 1 ) {
		std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
	} else {
		std::stringstream err_msg;
		err_msg << "You didn’t provide a single PDB file with the -in::file::s option" << std::endl;
		utility_exit_with_message(err_msg.str());
	}
	core::pose::PoseOP pose = core::import_pose::pose_from_file( filenames[1] );

	//read in a vector of chains present in the input PDB. To be used to check against input interfaces
	utility::vector1< char > chains;
	for ( core::Size i=1; i<=pose->num_chains(); i++ ) {
		chains.push_back(core::pose::get_chain_from_chain_id(i, *pose));
	}

	//loop through the interfaces
	for ( core::Size i=1; i<=edges.size(); i++ ) {

		//check to make sure both chains in the interface are in the input PDB
		std::stringstream err_msg;
		if ( std::find(chains.begin(), chains.end(), edges[i][1][0]) == chains.end() ) { //check first chain
			err_msg << "Chain " << edges[i][1][0] << " not in PDB." << std::endl;
			utility_exit_with_message(err_msg.str());
		}
		if ( std::find(chains.begin(), chains.end(), edges[i][1][2]) == chains.end() ) { //check second chain
			err_msg << "Chain " << edges[i][1][2] << " not in PDB." << std::endl;
			utility_exit_with_message(err_msg.str());
		}

		std::string interface(edges[i][1]);

		//Initialize IA mover
		protocols::analysis::InterfaceAnalyzerMoverOP IAM = utility::pointer::make_shared< protocols::analysis::InterfaceAnalyzerMover >(interface, true, core::scoring::ScoreFunctionFactory::create_score_function("ref2015")/*scorefxn_*/, false/*compute_packstat_*/, false/*pack_together_*/, false/*pack_separated_*/);

		core::Real dSASA;
		core::Real PRE;
		//check interface symmetry and (if valid), calculate dSASA and PRE for first interface
		check_interface_symmetry(pose, edges[i], dSASA, PRE);

		B[i] = w_SA*dSASA + w_PRE*PRE + w_int;
		TR << "Interface " << edges[i][1] << ": B = " << B[i] << std::endl;
	}

	//if dimer, check for disulfide bonds to use different A
	if ( n_chains==2 ) {
		bool is_DS = false;
		//loop through the two chains
		for ( core::Size i=1; i<=2; i++ ) {
			core::pose::PoseOP current_chain(pose->split_by_chain(i));
			core::Real num_res = current_chain->size();
			core::Size DS = 0; //intrasubunit disulfide bond calculation
			for ( core::Size j=1; j<=num_res; j++ ) {
				for ( core::Size k=1; k<=num_res; k++ ) {
					if ( current_chain->residue(j).is_bonded(k) && ! current_chain->residue(j).is_polymer_bonded(k) ) {
						DS++;
					}
				}
			}
			if ( DS>0 ) {
				is_DS = true;
			}
		}
		//if there is a disulfide bond, use different steepness
		if ( is_DS && !option[ steepness ].user() ) {
			A = 0.01;
			TR << "Disulfide bond detected for dimer. A = " << A << std::endl;
		}
	}
}

//function to simulate ERMS
void simulate_ERMS( utility::vector1<utility::vector1<core::Real>> &ERMS_prediction, const utility::vector1<core::Real> &ACE, const utility::vector1<core::Real> &B, const core::Size n_chains, const utility::vector1<utility::vector1<std::string>> &edges, const core::Real A, const core::Real breakage_cut, const utility::vector1<char> &nodes) {
	typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS > Graph;
	core::Size n_sims = 1000;

	TR << "Steepness = " << A << std::endl;

	//loop through the acceleration energies
	for ( core::Size i=1; i<=ACE.size(); i++ ) {
		core::Real X = ACE[i];
		utility::vector1<core::Real> PB; //probability of interface breakage
		//loop through each B value to calculate the breakage probabilities
		for ( core::Size j=1; j<=B.size(); j++ ) {
			PB.push_back(calc_prob(X, A, B[j]));
		}

		utility::vector1<utility::vector1<core::Real>> n_complexes_all;
		//loop through each simulation run
		for ( core::Size j=1; j<=n_sims; j++ ) {
			Graph G(n_chains);
			//loop through each interface type
			for ( core::Size k=1; k<=B.size(); k++ ) {
				//loop through each interface of that type
				for ( core::Size l=1; l<=edges[k].size(); l++ ) {
					core::Real random_num = numeric::random::uniform();
					if ( random_num > PB[k] ) {
						int node_0(nodes.index(edges[k][l][0])-1);
						if ( node_0==-1 ) {
							std::stringstream err_msg;
							err_msg << "Chain " << edges[k][l][0] << " does not match input chains." << std::endl;
							utility_exit_with_message(err_msg.str());
						}
						int node_1(nodes.index(edges[k][l][2])-1);
						if ( node_1==-1 ) {
							std::stringstream err_msg;
							err_msg << "Chain " << edges[k][l][2] << " does not match input chains." << std::endl;
							utility_exit_with_message(err_msg.str());
						}
						add_edge(node_0, node_1, G);
					}
				}
			}
			utility::vector1<core::Size> c(num_vertices(G));
			core::Size num = connected_components(G, make_iterator_property_map(c.begin(), get(boost::vertex_index, G), c[1]));
			utility::vector1<core::Size>::iterator k;
			utility::vector1<core::Size> n_components(num, 0);
			for ( k = c.begin(); k != c.end(); ++k ) {
				n_components[*k+1]++;
			}
			utility::vector1<core::Real> n_complexes;
			for ( core::Size k=1; k<=n_chains; k++ ) {
				n_complexes.push_back(count(n_components.begin(), n_components.end(), k));
			}
			n_complexes_all.push_back(n_complexes);
		}
		//calculate the average complexes
		utility::vector1<core::Real> avg_complexes(n_chains, 0.0);
		for ( core::Size j=1; j<=n_sims; j++ ) {
			for ( core::Size k=1; k<=n_chains; k++ ) {
				avg_complexes[k] += n_complexes_all[j][k];
			}
		}
		for ( core::Size j=1; j<=n_chains; j++ ) {
			avg_complexes[j] = avg_complexes[j]/n_sims;
		}
		//normalize average complexes
		if ( X==0 or X <= breakage_cut ) { //breakage not yet allowed, all precursor
			for ( core::Size j=1; j<=n_chains; j++ ) {
				if ( j==n_chains ) {
					ERMS_prediction[j].push_back(1.0);
				} else {
					ERMS_prediction[j].push_back(0.0);
				}
			}
		} else {
			for ( core::Size j=1; j<=n_chains; j++ ) {
				ERMS_prediction[j].push_back(avg_complexes[j]*(j)/n_chains);
			}
		}
	}
}

void calc_RMSE(const utility::vector1<utility::vector1<core::Real>> & ERMS_prediction, const utility::vector1<utility::vector1<core::Real>> &ERMS, const core::Size n_chains, const core::Size n_ACE) {
	core::Real RMSE_val(0.0);
	core::Size count(0);
	for ( core::Size i=1; i<=n_ACE; i++ ) {
		for ( core::Size j=1; j<=n_chains; j++ ) {
			RMSE_val += pow(ERMS_prediction[n_chains-(j-1)][i]-ERMS[i][j],2);
			count +=1;
		}
	}

	RMSE_val = pow(RMSE_val/count, 0.5);
	TR << "RMSE = " << RMSE_val << std::endl;
}

void output_results (const utility::vector1<core::Real> &ACE, const utility::vector1<utility::vector1<core::Real>> &ERMS_prediction) {
	std::ostringstream out;
	for ( core::Size i=1; i<=ACE.size(); i++ ) {
		out << ACE[i] << "\t";
		for ( core::Size j=1; j<=ERMS_prediction.size(); j++ ) {
			out << ERMS_prediction[ERMS_prediction.size()-(j-1)][i] << "\t";
		}
		out << std::endl;
	}
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	std::string outfile("ERMS_prediction.tsv");
	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();
	utility::io::ozstream outz( outfile.c_str() );
	outz << out.str();
	outz.close();
	outz.clear();
}

int
main( int argc, char * argv [] )
{

	try{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring;

		//local options
		option.add( complex_type, "TSV file containing the information on complex type (nodes and edges)" ).def("comp.tsv");
		option.add( ERMS, "TSV file containing the ERMS data or acceleration energies" ).def("ERMS.tsv");
		option.add( RMSE, "Calculate RMSE?" ).def(false);
		option.add( B_vals, "TSV file containing the B values (eV) for the interfaces" ).def("B.tsv");
		option.add( steepness, "Steepness value (A) as input by the user" );
		option.add( breakage_cutoff, "Breakage cutoff for simulation (eV)" ).def(0.0);

		devel::init( argc, argv );

		//default steepness (A) parameter
		core::Real A_(0.0025);

		//breakage cutoff
		core::Real breakage_cutoff_(0.0);
		if ( option[ breakage_cutoff ].user() ) {
			breakage_cutoff_ = option[ breakage_cutoff ]();
			if ( breakage_cutoff_ < 0 ) {
				std::stringstream err_msg;
				err_msg << "Breakage cutoff must be greater than zero. Value of " << breakage_cutoff_ << " eV given." << std::endl;
				utility_exit_with_message(err_msg.str());
			}
			A_ = 0.002; //default steepness when using breakage cutoff
		}

		//get steepness parameter from command line if given
		if ( option[ steepness ].user() ) {
			A_ = option[ steepness ]();
			TR.Warning << "Custom steepness of " << A_ << " input. Not recommended." << std::endl;
		}

		//read in complex type (nodes and edges)
		core::Size n_chains_;
		utility::vector1<char> nodes_;
		utility::vector1<utility::vector1<std::string>> edges_;
		if ( option[ complex_type ].user() ) {
			read_complex_type(n_chains_, nodes_, edges_);
		} else {
			std::stringstream err_msg;
			err_msg << "Must input the complex type using the complex_type option." << std::endl;
			utility_exit_with_message(err_msg.str());
		}
		TR << "Complex type read in" << std::endl;

		//set up vector of B values. Length is number of different interface types
		utility::vector1<core::Real> B_(edges_.size(), 0.0);
		if ( option[ B_vals ].user() ) {
			//read in B values if given by file input
			read_B_vals(B_);
			TR << "B values read in from file" << std::endl;
			if ( option[ basic::options::OptionKeys::in::file::s ].user() ) {
				TR.Warning << "PDB file input not used. B values read from file instead. To calculate B from structure, remove the flag (B_vals)." << std::endl;
			}
		} else {
			if ( !option[ basic::options::OptionKeys::in::file::s ].user() ) {
				std::stringstream err_msg;
				err_msg << "Must input either file containing B values (using flag B_vals) or PDB file to calculate B values for simulation." << std::endl;
				utility_exit_with_message(err_msg.str());
			}

			//B value calculation from structure
			calc_B_values(B_, edges_, A_, n_chains_); //read in PDB and calculate B values for each the first in each interface type
			TR << "B values read in from structure" << std::endl;
		}

		//read in ERMS data: acceleration energies and ERMS data (if applicable)
		utility::vector1<core::Real> ACE_; //Acceleration energies, read from file or input automatically
		utility::vector1<utility::vector1<core::Real>> ERMS_;
		if ( option[ ERMS ].user() ) {
			read_ERMS(ACE_, ERMS_, n_chains_);
			TR << "Acceleration energies and ERMS read in from file" << std::endl;
		} else {
			//if no acceleration energies are input, the range is set from 0 to 2x max B, with 10 steps in between
			core::Real step_size = (*max_element(B_.begin(), B_.end()))*2/10;
			ACE_.push_back(0.0);
			for ( core::Size i=1; i<=10; i++ ) {
				ACE_.push_back(i*step_size);
			}
			TR << "Acceleration energies set up based on max B" << std::endl;
		}


		//simulate ERMS
		utility::vector1<utility::vector1<core::Real>> ERMS_prediction_(n_chains_);
		simulate_ERMS(ERMS_prediction_, ACE_, B_, n_chains_, edges_, A_, breakage_cutoff_, nodes_);
		TR << "ERMS simulation complete" << std::endl;

		//calculte RMSE (if applicable)
		if ( option[ RMSE ].user() && option[ ERMS ].user() ) {
			calc_RMSE(ERMS_prediction_, ERMS_, n_chains_, ERMS_.size());
		}
		if ( option[ RMSE ].user() && !option[ ERMS ].user() ) {
			TR.Warning << "RMSE cannot be calculated without ERMS data input." << std::endl;
		}

		//output results to file
		output_results(ACE_, ERMS_prediction_);
		TR << "Output results complete." << std::endl;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
