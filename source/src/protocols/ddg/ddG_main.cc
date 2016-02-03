// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Liz Kellogg ekellogg@u.washington.edu

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <protocols/scoring/Interface.hh>
#include <core/conformation/Residue.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>


#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/ddg.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/database/open.hh>


#include <numeric/xyzVector.hh>

#include <fstream>
#include <iostream>
#include <sstream>

// C++ headers
#include <cstdlib>
#include <string>
#include <protocols/ddg/ddGMover.hh>
#include <protocols/ddg/ddG_main.hh>

#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>

using basic::T;
using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "pilot_apps.fix_bb_monomer_ddg" );


//numeric::random::RandomGenerator RG(54324); // <- Magic number, do not change it!!!

using namespace core;
using namespace scoring;

typedef utility::vector1<core::chemical::AA> mutations;
typedef utility::vector1<double> ddgs;

namespace protocols {
void
print_ddgs(std::string ddg_out,
	std::string label,
	ddgs delta_e_components,
	ddgs /*mut_avg_components*/,
	double total_ddgs
){
	std::ofstream ddg_output(ddg_out.c_str(), std::ios_base::app);
	if ( !ddg_output ) {
		TR << "having trouble opening output file for dumping predicted ddgs"
			<< ddg_out << std::endl;
		utility::exit(EXIT_FAILURE, __FILE__, __LINE__);
	}

	//  utility::vector1<std::string> scorefxn_header = get_wt_score.get_scorefunction_header();

	//  ddg_output << "ddG: description total ";
	//  for(Size i =1; i <=scorefxn_header.size();i++){
	//    ddg_output << scorefxn_header[i] << " ";
	//  }
	if ( label.compare("") != 0 ) {
		using namespace ObjexxFCL::format;
		ddg_output << "ddG: " << label << " " << F(9,3,total_ddgs) << " ";
		for ( Size m=1; m<=delta_e_components.size(); m++ ) {
			ddg_output << F(9,3,delta_e_components[m]) << " ";
		}
		ddg_output << "\n";
	}

	ddg_output << std::endl;
}

void
read_in_mutations(utility::vector1<mutations> & res_to_mut, std::string filename, pose::Pose & pose){
	std::ifstream inputstream;
	inputstream.open(filename.c_str());
	if ( inputstream.is_open() ) {
		int total; std::string total_keyword;
		inputstream >> total_keyword;
		assert(total_keyword.compare("total") == 0);

		inputstream >> total; //keep for cross-checking
		while ( !inputstream.eof() ) {
			mutations current_mutation(pose.total_residue(),core::chemical::aa_unk);
			int num_mutations;
			inputstream >> num_mutations;
			while ( num_mutations > 0 ) {
				char wt; int resnum; char mut;
				inputstream >> wt >> resnum >> mut;
				TR << "wt is " << wt << " resnum is " << resnum << " and mut is " << mut << std::endl;
				assert(pose.residue(resnum).name1() == wt);
				core::chemical::AA mutation= core::chemical::aa_from_oneletter_code(mut);
				current_mutation[resnum]=mutation;
				num_mutations--; total--;
			}
			TR << "end reading mutations for this" << std::endl;
			if ( num_mutations < 0 ) {
				TR.Error << "number of mutations mismatch! num_mutations < 0" << std::endl;
				return;
			} else {
				res_to_mut.push_back(current_mutation);
			}
		}
		if ( total < 0 ) {
			TR.Error << "total number of mutations mismatch! total < 0" << std::endl;
			return;
		}
	}
}

int
ddG_main()
{
	using namespace pose;
	using namespace scoring;
	using namespace conformation;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pack::task;
	using namespace protocols::moves;


	OPT(ddg::weight_file);
	OPT(ddg::iterations);
	OPT(ddg::debug_output);
	OPT(ddg::dump_pdbs);
	OPT(ddg::out);
	OPT(ddg::interface_ddg);
	OPT(ddg::opt_radius);
	// OPT(score::weights);
	// OPT(score::patch);
	OPT(in::file::s);

	// setup random numbers and options
	//  devel::init(argc, argv);

	// read the pose
	pose::Pose pose;
	core::import_pose::pose_from_file( pose, basic::options::start_file() , core::import_pose::PDB_file); // gets filename from -s option

	std::string weight_file = option[ OptionKeys::ddg::weight_file ]();
	basic::options::option[ score::fa_max_dis ](9.0); //set fa_max_dis before scorefunction is created!
	ScoreFunctionOP score_structure_scorefxn(ScoreFunctionFactory::create_score_function(weight_file));

	ScoreFunctionOP minimize_sfxn;
	if ( basic::options::option[OptionKeys::ddg::minimization_scorefunction].user() && basic::options::option[OptionKeys::ddg::minimization_patch].user() ) {
		minimize_sfxn=ScoreFunctionFactory::create_score_function(basic::options::option[OptionKeys::ddg::minimization_scorefunction](),
			basic::options::option[OptionKeys::ddg::minimization_patch]());
	} else if ( basic::options::option[OptionKeys::ddg::minimization_scorefunction].user() ) {
		minimize_sfxn=ScoreFunctionFactory::create_score_function(basic::options::option[OptionKeys::ddg::minimization_scorefunction]());
	} else {
		minimize_sfxn=get_score_function();
	}

	int num_iterations = option[ OptionKeys::ddg::iterations ]();
	bool opt_nbrs = false;
	double cutoff = -1;
	if ( basic::options::option[ OptionKeys::ddg::opt_radius].user() ) {
		opt_nbrs = true;
		cutoff = basic::options::option[ OptionKeys::ddg::opt_radius ]();
	} else if ( basic::options::option[OptionKeys::ddg::local_opt_only]() ) {
		opt_nbrs = true;
		cutoff = 8.0; //default cutoff
	}

	//initialize output options.
	//debug output?
	bool debug_output = option[ OptionKeys::ddg::debug_output ]();
	if ( debug_output ) {
		TR << "weights being used: " <<
			score_structure_scorefxn->weights() << "\n";
	}

	//dump repacked pdbs?
	bool dump_pdbs = option[ OptionKeys::ddg::dump_pdbs ]();

	//output ddgs into what file?
	std::string ddg_out = option[ OptionKeys::ddg::out ]();


	//interface mode?
	bool interface_ddg = option[ OptionKeys::ddg::interface_ddg ]();

	//minimize after repacking?
	bool min_cst = option[OptionKeys::ddg::min_cst]();

	//take mean or min energy as predicted ddg?
	bool mean = option[OptionKeys::ddg::mean]();
	bool min = option[OptionKeys::ddg::min]();

	ObjexxFCL::FArray2D<double> wt_scores(20,num_iterations);

	utility::vector1<core::chemical::AA> all_unk(pose.total_residue(),core::chemical::aa_unk);

	utility::vector1<double> wt_averaged_score_components;
	utility::vector1<ddgs> delta_energy_components;
	utility::vector1<double> total_ddgs;
	utility::vector1<ddgs> mutant_averaged_score_components;
	utility::vector1<std::string> delta_delta_G_label;

	ddg::ddGMover get_wt_score(score_structure_scorefxn,minimize_sfxn,all_unk);

	if ( basic::options::option[OptionKeys::ddg::mut_only].user() && basic::options::option[OptionKeys::ddg::mut_only]() ) {}
	else {


		get_wt_score.set_min_cst(min_cst);
		get_wt_score.set_min(min);
		get_wt_score.set_mean(mean);
		if ( !opt_nbrs ) {
			get_wt_score.restrict_to_nbrs(opt_nbrs);
			get_wt_score.neighbor_cutoff(cutoff);
			get_wt_score.num_iterations(num_iterations);
			get_wt_score.dump_pdbs(dump_pdbs);
			get_wt_score.is_interface_ddg(interface_ddg);
			get_wt_score.debug_output(debug_output);
			get_wt_score.num_iterations(num_iterations);
			get_wt_score.residues_to_mutate(all_unk);
			get_wt_score.apply(pose);
			wt_averaged_score_components=get_wt_score.get_wt_averaged_score_components();
		}

		if ( basic::options::option[OptionKeys::ddg::wt_only].user() && basic::options::option[OptionKeys::ddg::wt_only]() ) {
			std::cout <<  "optimizing only wt structure. exiting" << std::endl;
			return 1; //early exit
		}
	}
	if ( option[ OptionKeys::ddg::mut_file ].user() ) { //check if mutfile is specified
		TR << "reading in mutfile" << std::endl;
		std::string filename = option[OptionKeys::ddg::mut_file]();
		utility::vector1<mutations> res_to_mut;
		read_in_mutations(res_to_mut, filename,pose);
		TR << "size of res_to_mut is: " << res_to_mut.size() << std::endl;
		//initialize wildtype scores
		for ( Size i=1;  i <= res_to_mut.size(); i++ ) {
			utility::vector1<core::chemical::AA> residues_to_mutate = res_to_mut[i];
			bool mutation_defined = false; //to check if any mutation was specified
			for ( Size m =1; m<= residues_to_mutate.size(); m++ ) {
				if ( residues_to_mutate[m] != core::chemical::aa_unk ) {
					mutation_defined=true;
				}
			}
			if ( mutation_defined ) {
				ddg::ddGMover point_mutation(score_structure_scorefxn,minimize_sfxn,residues_to_mutate);
				point_mutation.set_min_cst(min_cst);
				point_mutation.set_min(min);
				point_mutation.set_mean(mean);
				if ( !opt_nbrs && get_wt_score.is_wt_calc_complete() &&
						!(basic::options::option[OptionKeys::ddg::mut_only].user() && basic::options::option[OptionKeys::ddg::mut_only]()) ) {
					TR << "testing if wt calc is complete. should be complete!" << std::endl;
					point_mutation.wt_score_components(get_wt_score.wt_score_components());
				}
				point_mutation.restrict_to_nbrs(opt_nbrs);
				point_mutation.neighbor_cutoff(cutoff);
				point_mutation.dump_pdbs(dump_pdbs);
				point_mutation.debug_output(debug_output);
				point_mutation.num_iterations(num_iterations);
				point_mutation.apply(pose);
				delta_delta_G_label.push_back(point_mutation.mutation_label(pose));
				TR << "mutation label for this round is " << point_mutation.mutation_label(pose) << std::endl;
				if ( point_mutation.is_wt_calc_complete() &&
						point_mutation.is_mutant_calc_complete() ) {
					//TR << " both calculations are complete so start storing info!" << std::endl;
					//output everything or store everything for output later
					delta_energy_components.push_back(point_mutation.get_delta_energy_components());
					mutant_averaged_score_components.push_back(point_mutation.get_mutant_averaged_score_components());
					total_ddgs.push_back(point_mutation.ddG());
					print_ddgs(ddg_out,
						point_mutation.mutation_label(pose),
						point_mutation.get_delta_energy_components(),
						point_mutation.get_mutant_averaged_score_components(),
						point_mutation.ddG());
				}
			}
		}
	}

	if ( option[packing::resfile].user() ) { //check is resfile is specified
		pack::task::PackerTaskOP storage_task(pack::task::TaskFactory::create_packer_task(pose));

		storage_task->initialize_from_command_line();
		pack::task::parse_resfile(pose, *storage_task);
		storage_task->or_include_current(true);

		for ( Size i =1; i<=pose.total_residue(); i++ ) {
			if ( storage_task->design_residue(i) ) {
				for ( ResidueLevelTask::ResidueTypeCOPListConstIter aa_iter(storage_task->residue_task(i).allowed_residue_types_begin()),
						aa_end(storage_task->residue_task(i).allowed_residue_types_end());
						aa_iter != aa_end; ++aa_iter ) {
					utility::vector1<core::chemical::AA> residues_to_mutate(pose.total_residue(),core::chemical::aa_unk);
					residues_to_mutate[i]=((*aa_iter)->aa());
					if ( residues_to_mutate[i] != core::chemical::aa_unk ) {
						ddg::ddGMover point_mutation(score_structure_scorefxn,minimize_sfxn,residues_to_mutate);
						point_mutation.set_min_cst(min_cst);
						point_mutation.set_min(min);
						point_mutation.set_mean(mean);
						//initialize wildtype scores
						if ( !opt_nbrs && get_wt_score.is_wt_calc_complete() ) {
							TR << "testing if wt calc is complete. should be complete!" << std::endl;
							point_mutation.wt_score_components(get_wt_score.wt_score_components());
						}
						point_mutation.restrict_to_nbrs(opt_nbrs);
						point_mutation.neighbor_cutoff(cutoff);
						point_mutation.dump_pdbs(dump_pdbs);
						point_mutation.debug_output(debug_output);
						point_mutation.num_iterations(num_iterations);
						point_mutation.apply(pose);
						delta_delta_G_label.push_back(point_mutation.mutation_label(pose));
						if ( point_mutation.is_wt_calc_complete() &&
								point_mutation.is_mutant_calc_complete() ) {
							//TR << " both calculations are complete so start storing info!" << std::endl;
							//output everything
							delta_energy_components.push_back(point_mutation.get_delta_energy_components());
							mutant_averaged_score_components.push_back(point_mutation.get_mutant_averaged_score_components());
							total_ddgs.push_back(point_mutation.ddG());
							//output information to file
							print_ddgs(ddg_out,
								point_mutation.mutation_label(pose),
								point_mutation.get_delta_energy_components(),
								point_mutation.get_mutant_averaged_score_components(),
								point_mutation.ddG());
						} else {
							std::cout << "cannot output ddgs. either mut or wt arrays are incomplete\ni hope you don't mind!\n" << std::endl;
						}
					}
				}
			}
		}
	}

	//INTERFACE MODE
	if ( interface_ddg ) {

		//TR << "[DEBUG]: now  in interface mode"<< std::endl;
		//detect interface residues
		using namespace core;
		using namespace core::conformation;
		using namespace core::chemical;

		utility::vector1<core::chemical::AA> residues_to_mutate;
		//set up interface object
		protocols::scoring::Interface protein_interface(1);
		protein_interface.distance(10.0);
		protein_interface.calculate(pose);
		protein_interface.print(pose);

		//debug statement
		for ( Size i =1; i<=pose.total_residue(); i++ ) {
			if ( protein_interface.is_interface(i) ) {
				TR << "[DEBUG]:" << i << " is in the interface " << std::endl;
			}
		}
		//debug statement end

		for ( Size i =1; i<=pose.total_residue(); i++ ) {
			if ( protein_interface.is_interface(i) ) { //is interface residue
				for ( Size j =1; j <= 20 ; j++ ) { //iterate through all amino acids
					residues_to_mutate = all_unk; //mutate each interface residue one at a time
					core::chemical::AA curr_aa = (core::chemical::AA)j;
					if ( curr_aa != pose.aa(i) && (pose.aa(i) != aa_unk)/*this hopefully will never happen?*/ ) {
						residues_to_mutate[i]=curr_aa;
						ddg::ddGMover interface_mutation(score_structure_scorefxn,minimize_sfxn,residues_to_mutate);
						interface_mutation.set_min_cst(min_cst);
						interface_mutation.is_interface_ddg(interface_ddg);
						interface_mutation.set_min(min);
						interface_mutation.set_mean(mean);
						if ( get_wt_score.is_wt_calc_complete() ) {
							TR << "testing if wt calc is complete. should be complete!" << std::endl;
							interface_mutation.wt_score_components(get_wt_score.wt_score_components());
							interface_mutation.wt_unbound_score_components(get_wt_score.wt_unbound_score_components());
						}

						if ( dump_pdbs ) interface_mutation.dump_pdbs(dump_pdbs);
						if ( debug_output ) interface_mutation.debug_output(debug_output);

						interface_mutation.num_iterations(num_iterations);
						interface_mutation.apply(pose);
						delta_delta_G_label.push_back(interface_mutation.mutation_label(pose));
						TR << "mutation label for this round is " << interface_mutation.mutation_label(pose) << std::endl;
						if ( interface_mutation.is_wt_calc_complete() &&
								interface_mutation.is_mutant_calc_complete() ) {

							delta_energy_components.push_back(interface_mutation.get_delta_energy_components());
							mutant_averaged_score_components.push_back(interface_mutation.get_mutant_averaged_score_components());
							total_ddgs.push_back(interface_mutation.ddG());
							print_ddgs(ddg_out,
								interface_mutation.mutation_label(pose),
								interface_mutation.get_delta_energy_components(),
								interface_mutation.get_mutant_averaged_score_components(),
								interface_mutation.ddG());
							TR << "interface mutation is complete and ddg is: " << interface_mutation.ddG() << std::endl;
						}
					}
				}//iterate through all amino acids
			}
		}
	}

	//format and output all the stored information
	//  utility::vector1<std::string> scorefxn_header = get_wt_score.get_scorefunction_header();

	return 1;
}
} //namespace protocols

