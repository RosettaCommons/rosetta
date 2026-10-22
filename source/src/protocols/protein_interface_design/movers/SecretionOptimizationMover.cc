// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_interface_design/movers/SecretionOptimizationMover.cc
/// @author John Wang (jyjwang@uw.edu)

// Unit Headers
#include <protocols/protein_interface_design/movers/SecretionOptimizationMover.hh>
#include <protocols/protein_interface_design/movers/SecretionOptimizationMoverCreator.hh>

#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <numeric/xyz.functions.hh>
#include <core/pose/extra_pose_info_util.hh>


// Packer Headers
#include <core/conformation/Residue.hh>
#include <core/pack/task/ResidueLevelTask.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/chains_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

// TaskOp Headers
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>

#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/score_filters/ScoreTypeFilter.hh>
#include <protocols/simple_moves/DumpPdb.hh>
#include <protocols/jd2/util.hh>

#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <map>
#include <algorithm>
#include <unordered_map>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <boost/functional/hash.hpp>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {
static basic::Tracer TR( "protocols.protein_interface_design.movers.SecretionOptimizationMover" );
//@details when the Degreaser is initialized, the scorefunction is cloned, not used directly.
SecretionOptimizationMover::SecretionOptimizationMover():
	scorefxn_( /* NULL */ )
{
	SecretionOptimizationMover::make_table();
}

SecretionOptimizationMover::SecretionOptimizationMover(core::scoring::ScoreFunctionCOP scorefxn):
	scorefxn_(scorefxn->clone())
{
	SecretionOptimizationMover::make_table();
}

// Tag Parsing
void SecretionOptimizationMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap  & data
)
{
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	set_max_mutations( tag->getOption<core::Size>("max_mutations",3) );
	set_dG_ins_threshold( tag->getOption<core::Real>("dG_ins_threshold",2.7) );
	set_score_tolerance( tag->getOption<core::Real>("score_tolerance",0.0) );
	set_aas_allowed( tag->getOption<std::string>("aas_allowed","DEKRQN") );
	set_chain_to_degrease( tag->getOption<std::string>("chain_to_degrease","A") );
	repack_shell( tag->getOption<core::Real>("repack_shell",8.0) );
	set_dG_tolerance( tag->getOption<core::Real>("dG_tolerance",0.27) );
	set_dump_multis( tag->getOption<bool>("dump_multis",false) );
	set_dump_singles( tag->getOption<bool>("dump_singles",false) );
	set_dump_dir( tag->getOption<std::string>("dump_dir",".") );
	use_largest_ddg( tag->getOption<bool>("largest_ddg",true) );
	use_smallest_dscore( tag->getOption<bool>("smallest_dscore",false) );
	// Construct score function.
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();

}

//@brief find designable residues, then allow packing in a certain size shell around each designable residue. based on DesignAroundApply TaskOp by Neil King
void
SecretionOptimizationMover::DesignAroundApply( core::pose::Pose const & pose, utility::vector1<mutt> & targets) const
{
	using namespace core::pack::task::operation;
	core::pack::task::TaskFactoryOP task_factory_(new core::pack::task::TaskFactory);

	std::set<core::Size> focus_residues;
	RestrictAbsentCanonicalAAS set_focus;
	RestrictResidueToRepacking repack;
	PreventRepacking no_repack;
	focus_residues.clear();

	//go through focus residue as designated by targets. each mutt object has a resnum and set of aa's to design to.
	for ( core::Size i = 1; i <= targets.size(); ++i ) {
		focus_residues.insert(targets[i].resnum);
		set_focus.include_residue(targets[i].resnum);
		set_focus.keep_aas(targets[i].aa);
		set_focus.apply(pose, *task_);
	}
	//scan along a (in the ASU [if symmetric] or chain [if asymmetric]) to be packed. if it's NOT a focus residue, then see if it's near a reflection [if symmetric] of a focus residue. if so, pack it.
	core::Size all_residues = pose.total_residue();
	core::Size iter_a_residues = all_residues;
	bool symmetric = core::pose::symmetry::is_symmetric(pose);
	if ( symmetric ) {
		core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);

		all_residues = sym_info->num_total_residues_without_pseudo();
		iter_a_residues = sym_info->num_independent_residues();
	}
	core::Real const repack_shell_squared(repack_shell() * repack_shell());
	for ( core::Size a = 1; a <= iter_a_residues; ++a ) {
		//if it's a focus residue, don't touch it
		if ( !pose.residue(a).is_protein() ) {
			continue;
		}
		if ( focus_residues.count(a) ) {
			continue;
		}
		bool allow_packing(false);
		for ( core::Size b = 1; b <= all_residues; ++b ) {
			core::Size b_check = b;
			if ( symmetric ) {
				core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
				b_check = sym_info->chi_follows(b);
			}
			if ( focus_residues.find(b) != focus_residues.end() || focus_residues.find(b_check) != focus_residues.end() ) {
				core::Real const distance_sq( pose.residue(a).xyz(pose.residue(a).nbr_atom()).distance_squared(pose.residue(b).xyz(pose.residue(b).nbr_atom())) );
				if ( distance_sq <= repack_shell_squared ) {
					allow_packing=true;
					break;
				}
			}
		}//for b
		if ( allow_packing ) {
			task_->nonconst_residue_task(a).restrict_to_repacking();
			repack.include_residue(a);
		} else {
			task_->nonconst_residue_task(a).prevent_repacking();
			no_repack.include_residue(a);
		}
	}// for a
	repack.apply(pose, *task_);
	no_repack.apply(pose, *task_);
	task_factory_->modify_task(pose, task_);
}

//@brief call and execute the packer after task has been modified
void SecretionOptimizationMover::packeruser( core::pose::Pose & pose) const
{
	//task_ is being updated by DesignAroundapply. we'll make a new SymPackRotamersMover (because it handles symmetric and asymmetric poses) and pack given the task_
	protocols::minimization_packing::PackRotamersMoverOP pack_it_( utility::pointer::make_shared< protocols::minimization_packing::PackRotamersMover >(scorefxn_, task_, 1));
	//pack_it_->task_factory(task_factory());
	pack_it_->apply(pose);
}

//@brief Build a table, the hard way
void
SecretionOptimizationMover::make_table()
{
	coeffs["A"]["a0"] = 1.27e-01;
	coeffs["A"]["a1"] = 2.15e-02;
	coeffs["A"]["a2"] = 0;
	coeffs["A"]["a3"] = 0;
	coeffs["A"]["a4"] =0;
	coeffs["C"]["a0"] = -7.65e-02;
	coeffs["C"]["a1"] = 9.94e-02;
	coeffs["C"]["a2"] = 0;
	coeffs["C"]["a3"] = 0;
	coeffs["C"]["a4"] =0;
	coeffs["D"]["a0"] = 1.79e00;
	coeffs["D"]["a1"] = 1.73e-02;
	coeffs["D"]["a2"] = 0;
	coeffs["D"]["a3"] = 0;
	coeffs["D"]["a4"] =0;
	coeffs["E"]["a0"] = 1.42e00;
	coeffs["E"]["a1"] = 8.94e-03;
	coeffs["E"]["a2"] = 0;
	coeffs["E"]["a3"] = 0;
	coeffs["E"]["a4"] =0;
	coeffs["F"]["a0"] = -2.77e-01;
	coeffs["F"]["a1"] = 1.03e-03;
	coeffs["F"]["a2"] = 0;
	coeffs["F"]["a3"] = 0;
	coeffs["F"]["a4"] =0;
	coeffs["G"]["a0"] = 4.81e-01;
	coeffs["G"]["a1"] = 4.72e-03;
	coeffs["G"]["a2"] = 0;
	coeffs["G"]["a3"] = 0;
	coeffs["G"]["a4"] =0;
	coeffs["H"]["a0"] = 1.20e00;
	coeffs["H"]["a1"] = 8.01e-03;
	coeffs["H"]["a2"] = 0;
	coeffs["H"]["a3"] = 0;
	coeffs["H"]["a4"] =0;
	coeffs["I"]["a0"] = -4.60e-01;
	coeffs["I"]["a1"] = 1.81e-02;
	coeffs["I"]["a2"] = 0;
	coeffs["I"]["a3"] = 0;
	coeffs["I"]["a4"] =0;
	coeffs["K"]["a0"] = 1.85e00;
	coeffs["K"]["a1"] = 2.18e-02;
	coeffs["K"]["a2"] = 0;
	coeffs["K"]["a3"] = 0;
	coeffs["K"]["a4"] =0;
	coeffs["L"]["a0"] = -4.28e-01;
	coeffs["L"]["a1"] = 2.38e-03;
	coeffs["L"]["a2"] = 0;
	coeffs["L"]["a3"] = 0;
	coeffs["L"]["a4"] =0;
	coeffs["M"]["a0"] = -7.75e-02;
	coeffs["M"]["a1"] = 9.84e-02;
	coeffs["M"]["a2"] = 0;
	coeffs["M"]["a3"] = 0;
	coeffs["M"]["a4"] =0;
	coeffs["N"]["a0"] = 1.33e00;
	coeffs["N"]["a1"] = 9.24e-03;
	coeffs["N"]["a2"] = 0;
	coeffs["N"]["a3"] = 0;
	coeffs["N"]["a4"] =0;
	coeffs["P"]["a0"] = 1.09e00;
	coeffs["P"]["a1"] = 1.01e-02;
	coeffs["P"]["a2"] = 0;
	coeffs["P"]["a3"] = 0;
	coeffs["P"]["a4"] =0;
	coeffs["Q"]["a0"] = 1.33e00;
	coeffs["Q"]["a1"] = 1.12e-02;
	coeffs["Q"]["a2"] = 0;
	coeffs["Q"]["a3"] = 0;
	coeffs["Q"]["a4"] =0;
	coeffs["R"]["a0"] = 1.65e00;
	coeffs["R"]["a1"] = 5.12e-02;
	coeffs["R"]["a2"] = 0;
	coeffs["R"]["a3"] = 0;
	coeffs["R"]["a4"] =0;
	coeffs["S"]["a0"] = 7.02e-01;
	coeffs["S"]["a1"] = 7.77e-03;
	coeffs["S"]["a2"] = 0;
	coeffs["S"]["a3"] = 0;
	coeffs["S"]["a4"] =0;
	coeffs["T"]["a0"] = 5.27e-01;
	coeffs["T"]["a1"] = 3.12e-02;
	coeffs["T"]["a2"] = 0;
	coeffs["T"]["a3"] = 0;
	coeffs["T"]["a4"] =0;
	coeffs["V"]["a0"] = -2.45e-01;
	coeffs["V"]["a1"] = 9.79e-02;
	coeffs["V"]["a2"] = 0;
	coeffs["V"]["a3"] = 0;
	coeffs["V"]["a4"] =0;
	coeffs["W"]["a0"] = 2.91e-01;
	coeffs["W"]["a1"] = 1.89e-02;
	coeffs["W"]["a2"] = -5.48e-01;
	coeffs["W"]["a3"] = 9.30e-02;
	coeffs["W"]["a4"] =6.47e00;
	coeffs["Y"]["a0"] = 6.28e-01;
	coeffs["Y"]["a1"] = 1.04e-02;
	coeffs["Y"]["a2"] = -5.74e-01;
	coeffs["Y"]["a3"] = 9.48e-02;
	coeffs["Y"]["a4"] =6.92e00;
}

//@brief Given a window size and a start position, return the window at that position
std::string
SecretionOptimizationMover::get_window_at_index( core::pose::Pose & pose, core::Size & wlen, core::Size & start )
{
	std::string seq = pose.sequence();
	return seq.substr(start-1,wlen);//thanks @mlnance for saving a line here
}

// @brief slightly cleaner way to calculate position-specific dG_ins contribution
core::Real SecretionOptimizationMover::dG_term_for_residue_at_position(core::Real & position, std::string & residue)
{
	core::Real a0 = coeffs[residue]["a0"]; //call coefficients from table
	core::Real a1 = coeffs[residue]["a1"]; //" "
	core::Real a2 = coeffs[residue]["a2"]; //" "
	core::Real a3 = coeffs[residue]["a3"]; //" "
	core::Real a4 = coeffs[residue]["a4"]; //" "

	return a0*exp(-a1*position*position) + a2*exp(-(a3*(position-a4)*(position-a4)))+a2*exp(-a3*(position+a4)*(position+a4)); //equation 4 in Hessa et. al 2007
}
core::Real SecretionOptimizationMover::normalize_position_in_helix_given_length(core::Size & position, core::Real & helix_length) //I needed to use core::Real for helix length because i needed to multiply them elsewhere
{
	return 9.0*((position/(helix_length-1.0))*2.0 - 1.0); //for a 19-residue window, you should just get -9 to +9 for position
}
//@brief Calculate the dG_insertion for a given window, equations directly from Hessa et. al 2007
core::Real
SecretionOptimizationMover::dG_ins_for_window( std::string & window )
{
	core::Real dG_tot = 0; //Real for the rounded dG_ins
	core::Real dG_aa = 0; //core::Reals for more precise numbers
	core::Real mu_tot = 0;
	core::Real length_term = 0;
	core::Real mu0 = 0;
	core::Real mu1 = 0;
	core::Real k = window.length();

	for ( core::Size j = 0; j < window.length(); j++ ) {
		core::Real i = normalize_position_in_helix_given_length(j, k);  // formula for calculating (symmetric) position within k-residue segment
		// below equation (6) in Hessa et. al 2007
		std::string AA = window.substr(j,1); //AA identity at that position
		core::Real dG_aa_i = dG_term_for_residue_at_position(i, AA);  // calculate the individual residue contribution using another method (for clarity)
		mu0 = mu0 + dG_aa_i*sin(100.0*(j+1)*numeric::constants::d::pi/180); //hydrophobic moment component calculation, equation 2
		mu1 = mu1 + dG_aa_i*cos(100.0*(j+1)*numeric::constants::d::pi/180); //hydrophobic moment component calculation, equation 2

		dG_aa = dG_aa + dG_aa_i; //direct summation of energetic terms, equation 3
	}

	mu_tot = 0.27*sqrt(mu0*mu0 + mu1*mu1); //hydrophobic moment term, equation 2
	length_term = 9.29 - 0.645*k + 0.00822*k*k; //length corrected equation 1, parameter values in supplementary table 2
	dG_tot = dG_aa + mu_tot + length_term; //total predicted dG_ insertion (equation 1)

	return dG_tot;
}


// @brief default destructor
SecretionOptimizationMover::~SecretionOptimizationMover() = default;

// @brief apply function. kinda sloppy right now but i'll make more classes soon(tm)
void
SecretionOptimizationMover::apply( core::pose::Pose & pose)
{
	runtime_assert_string_msg(task_factory_ !=nullptr, "A task needs to be non-null. Should only be an error if you've specifically deleted all Tasks.");
	//check decision making set by user
	//defaults to making the largest dG_ins change at the lowest dG_ins region in the protein
	if ( largest_ddg_ && smallest_dscore_ ) {
		TR.Warning << "You've enabled highest ddG and lowest score. Can't guarantee that those are the same pose, defaulting to highest ddG" << std::endl;
		largest_ddg_ = true;
		smallest_dscore_ = false;
	}
	if ( !largest_ddg_ && !smallest_dscore_ ) {
		TR.Warning << "You've disabled highest ddG and lowest score. You have to pick one, defaulting to highest ddG" << std::endl;
		largest_ddg_ = true;
	}
	//setup for score calculations
	core::scoring::ScoreFunctionOP sct_sfxn = scorefxn_->clone();
	protocols::score_filters::ScoreTypeFilter const energy_filter(sct_sfxn, core::scoring::total_score, 0);
	core::Real old_score = energy_filter.compute(pose);
	//report_test regurgitates the user-designated parameters to the tracer
	//and fills in found_regions_ so that we can operate on those regions
	report_test(TR, pose);
	if ( !found_regions_.size() ) {
		TR << "FOUND NO TRANSMEMBRANE REGIONS WITH GIVEN PARAMETERS" << std::endl;
		return;
	}
	//could be a threshold thing or could just be a protein that doesn't need degreasing, exits if there's nothing to do
	utility::vector1<tm_region> outputs;
	core::pose::Pose pose_clean = pose;
	core::pose::Pose pose_smallest_dscore_ = pose;
	core::pose::Pose pose_largest_ddg_for_lowest_dg_ = pose;
	core::Real global_lowest_dscore = score_tolerance_;
	core::Real global_lowest_dg = found_regions_[1].dG_ins;
	core::Real global_highest_ddg = 0;
	//initialize these private variables (a little clunky with poses)
	//we need to keep track of the stats on the original pose (pose_clean)
	//we also need some global parameters to decision-make at the end
	//
	//this loop iterates over each detected (putative) transmembrane domain in the protein
	for ( core::Size REGIONS = 1; REGIONS <= found_regions_.size(); ++REGIONS ) {
		core::Size i = 1;
		mutids_.clear();
		found_regions_[REGIONS].score = old_score; //stats for the original pose

		outputs.push_back(found_regions_[REGIONS]);  //outputs is a vector of mutts (the same type
		//of struct as found_regions_) that we use to
		//keep track of the possible outputs
		while ( i <= window_size_ ) //in each found region, there's a window of (default 19) residues to attempt to perturb
				{
			resids_.clear();
			resids_.push_back(found_regions_[REGIONS].possible_mutants[i]);  //possible mutants can shrink, so we clear and restore them
			//each time we try to mutate a certain position
			pose = pose_clean;
			char orig_res = pose_clean.residue(resids_[1].resnum).name1();
			task_ = task_factory_->create_task_and_apply_taskoperations(pose);  //inherit the task factory (in case the user has supplied one)
			DesignAroundApply(pose, resids_);     //designate designable / repackable residues with DesignArroundApply
			//task_factory_->modify_task(pose, task_);    //legacy - adding this here was necessary to actually apply a taskoperation
			char oldres = pose.residue(resids_[1].resnum).name1();  //separate tracker for what the residue could be designed to over loops
			//orig_res is from pose_clean and won't change
			//oldres could change as we try mutations at a particular position
			TR << "working on " << resids_[1].resnum << oldres << std::endl;
			TR << "allowed residues: " << resids_[1].aa << std::endl;
			packeruser(pose); // actual designing / repacking
			char newres = pose.residue(resids_[1].resnum).name1(); //evaluation of post-mutation + repacking metrics
			std::string wind = get_window_at_index(pose, window_size_, found_regions_[REGIONS].index); //calling the same functions as before but for the new sequence
			core::Real newdG = dG_ins_for_window(wind);

			core::Real new_score = energy_filter.compute(pose);

			if ( new_score - old_score > score_tolerance_ ) {
				TR << "new score with mutation " << oldres << resids_[1].resnum << newres << " is too high. throwing out this pose." << std::endl;  //we outright reject any position that cannot accommodate
				//a packer-decided mutation
				++i;
			} else if ( newdG - found_regions_[REGIONS].dG_ins < dG_tolerance_ ) {  //having a perturbation that fits score-wise but not dG_ins
				//wise is the tricky case, so we further iterate
				//to make sure that we don't miss out on
				//possible mutations that could be better than
				//for example, an A -> S mutation
				if ( newres == oldres ) {
					if ( newres == orig_res ) {
						TR << "new dG not high enough, throwing out original residue " << newres << std::endl;
						found_regions_[REGIONS].possible_mutants[i].aa.erase(std::remove(found_regions_[REGIONS].possible_mutants[i].aa.begin(),found_regions_[REGIONS].possible_mutants[i].aa.end(),newres), found_regions_[REGIONS].possible_mutants[i].aa.end()); //from the string of available mutants "possible_mutants"
						//we use the found_regions iterator in order to delete just one character
						//such that the publically-available "possible_mutants" does not include the residue that failed
						//this is specifically for the case that the mover tries
						//to put the residue that the pose started with
						if ( !resids_[1].aa.size() ) {
							TR << "we've run out of possible mutants. moving on from this position" << std::endl;
							++i;
							continue; //if we've run out of possible mutants, then we move on from this position
						}
						if ( resids_[1].aa.find(newres) == resids_[1].aa.npos ) {
							++i;
							continue;
							TR <<"SOMEHOW the packer put "<<newres<<" in there but it WASN'T EVEN ALLOWED!!!! skipping this position..." << std::endl;
							//this was happening quite a bit early on but should be resolved
						}
						continue;
					}
				}
				TR << "new dG not high enough. throwing out this pose and trying without "<<newres<<std::endl;
				found_regions_[REGIONS].possible_mutants[i].aa.erase(std::remove(found_regions_[REGIONS].possible_mutants[i].aa.begin(),found_regions_[REGIONS].possible_mutants[i].aa.end(),newres), found_regions_[REGIONS].possible_mutants[i].aa.end());
				//deleting the tried residue from the pool before retrying (as above)
			} else {
				TR << "got a hit! " << std::endl;
				TR << "old sequence was: " << found_regions_[REGIONS].sequence << '\t';
				TR << "new sequence is: " << wind << std::endl;
				TR << "old dG_ins was: " << found_regions_[REGIONS].dG_ins << '\t';
				TR << "new dG_ins is: " << newdG << std::endl;
				TR << "old score was: " << old_score << '\t';
				TR << "new score is: " << new_score << std::endl;
				TR << "this pose with mutation "<< oldres <<resids_[1].resnum<< newres<<" meets requirements. storing for later." << std::endl;
				//save metrics of the hit
				resids_[1].aa = newres;
				resids_[1].ddG_ins = newdG - found_regions_[REGIONS].dG_ins;
				resids_[1].dscore = new_score - old_score;
				if ( resids_[1].dscore < global_lowest_dscore ) {
					TR << "this pose is now the lowest dscore pose" << std::endl;
					pose_smallest_dscore_ = pose;
					global_lowest_dscore = resids_[1].dscore;
				}
				if ( found_regions_[REGIONS].dG_ins <= global_lowest_dg && resids_[1].ddG_ins > global_highest_ddg ) {
					TR << "this pose is now the greatest dG_ins change at the region of lowest dG_ins" << std::endl;
					global_lowest_dg = found_regions_[REGIONS].dG_ins;
					pose_largest_ddg_for_lowest_dg_ = pose;
					global_highest_ddg = resids_[1].ddG_ins;
				}
				if ( resids_[1].ddG_ins > 0 ) {
					mutids_.push_back(resids_[1]);
					outputs.push_back({found_regions_[REGIONS].index, wind, newdG, new_score, resids_});
				}
				//save the outputs for later
				if ( dump_singles_ ) {
					std::string name = protocols::jd2::current_output_name();
					name = name + "_" + oldres + std::to_string(resids_[1].resnum) + newres + ".pdb";
					pose.dump_pdb(dump_dir_+"/"+name);
				}
				//optional pdb dumping, appends the mutation to the job name
				++i;
			}
		}
		if ( !mutids_.size() ) {
			TR << "FOUND NO CANDIDATE MUTANTS (NOT EVEN ONE!) FOR THIS TM REGION GIVEN YOUR PARAMETERS" << std::endl;
			continue;
		}
		//mutant combiner - this thing is 'dumb' and just tries the top three candidates
		//that were previously found. good for cases where score tolerance is high
		//or if making a bunch of mutations doesn't significantly perturb the structure (e.g. some natural proteins)
		//this is more or less a copy of the above code, except multiple mutations are allowed at once
		//there may be a more clever way of doing this
		pose = pose_clean;
		task_ = task_factory_->create_task_and_apply_taskoperations(pose);
		std::sort(mutids_.rbegin(),mutids_.rend()); // RBEGIN AND REND

		utility::vector1<mutt>::const_iterator first = mutids_.begin();
		utility::vector1<mutt>::const_iterator last;

		if ( mutids_.size() >= max_mutations_ ) {
			last = mutids_.begin() + max_mutations_;
		} else {
			last = mutids_.end();
		}

		utility::vector1<mutt> to_mut(first, last);

		mutids_ = to_mut;

		DesignAroundApply(pose, mutids_);
		packeruser(pose);

		std::string wind = get_window_at_index(pose, window_size_, found_regions_[REGIONS].index);
		core::Real newdG = dG_ins_for_window(wind);
		core::Real new_score = energy_filter.compute(pose);

		if ( new_score - old_score > score_tolerance_ ) {
			TR << "new score for combination mutant is too high. throwing out this pose." << std::endl;
		} else if ( newdG - found_regions_[REGIONS].dG_ins < dG_tolerance_ ) {
			TR << "new dG for combination mutant not high enough. throwing out this pose." << std::endl;
		} else {
			TR << "combination mutant found !" << std::endl;
			TR << "old sequence was: " << found_regions_[REGIONS].sequence << '\t';
			TR << "new sequence is: " << wind << std::endl;
			TR << "old dG_ins was: " << found_regions_[REGIONS].dG_ins << '\t';
			TR << "new dG_ins is: " << newdG << std::endl;
			TR << "old score was: " << old_score << '\t';
			TR << "new score is: " << new_score << std::endl;
			if ( new_score - old_score < global_lowest_dscore ) {
				pose_smallest_dscore_ = pose;
				global_lowest_dscore = new_score - old_score;
			}
			if ( found_regions_[REGIONS].dG_ins <= global_lowest_dg && newdG - found_regions_[REGIONS].dG_ins > global_highest_ddg ) {
				pose_largest_ddg_for_lowest_dg_ = pose;
				global_highest_ddg = newdG - found_regions_[REGIONS].dG_ins;
			}
			outputs.push_back({found_regions_[REGIONS].index, wind, newdG, new_score, mutids_});
			if ( dump_multis_ ) {
				std::string name = protocols::jd2::current_output_name();
				name = name + "_mut" + std::to_string(found_regions_[REGIONS].index) + ".pdb";
				pose.dump_pdb(dump_dir_+"/"+name);
			}
		}// else
	} //for REGIONS

	if ( largest_ddg_ ) {
		TR << "using the pose with the largest delta(dG_ins) at the lowest identified dG_ins" << std::endl;
		pose = pose_largest_ddg_for_lowest_dg_;
	}
	if ( smallest_dscore_ ) {
		TR << "using hte pose with the smallest delta(score)" << std::endl;
		pose = pose_smallest_dscore_;
	}
	//once we've examined each found tm region, we can write the output metrics to the PDB
	std::string header = "index\t sequence\t   dG_ins\t    score\t  ddG_ins\t   dscore\t  mutants\n";
	core::pose::setPoseExtraScore(pose, "DEGREASER000", header);
	for ( core::Size o = 1; o <= outputs.size(); ++o ) {
		std::string mutants = " \t";
		if ( !outputs[o].possible_mutants.size() ) {
			TR << "no possible mutants for this output line!" << std::endl;
			continue;
		}
		for ( core::Size a = 1; a <= outputs[o].possible_mutants.size(); ++a ) {
			char orig_res = pose_clean.residue(outputs[o].possible_mutants[a].resnum).name1();
			char mut_res = outputs[o].possible_mutants[a].aa[0];
			if ( mut_res != orig_res && (outputs[o].possible_mutants.size() <= max_mutations_) ) {
				if ( a > 1 ) {
					mutants += ", ";
				}
				mutants += orig_res;
				mutants += std::to_string(outputs[o].possible_mutants[a].resnum);
				mutants += mut_res;
			}
		}
		core::Size ind = outputs[o].index;
		std::string wind = get_window_at_index(pose_clean, window_size_, ind);
		core::Real total_ddG_ins = outputs[o].dG_ins - dG_ins_for_window(wind);
		core::Real total_dscore = outputs[o].score - old_score;
		std::string output_string = std::to_string(outputs[o].index) + '\t' + outputs[o].sequence + '\t' + ObjexxFCL::format::F (9,3,outputs[o].dG_ins) + '\t' + ObjexxFCL::format::F(9,2,outputs[o].score) + '\t' + ObjexxFCL::format::F(9,3,total_ddG_ins) + '\t' + ObjexxFCL::format::F(9,3,total_dscore) + mutants + '\n';
		std::stringstream o_str;
		o_str << std::setfill('0') << std::setw(3) << o;
		core::pose::setPoseExtraScore( pose, "DEGREASER"+o_str.str(), output_string);
	}
	TR.flush();
	//move on with the pose that was decided with the parameters given
}

//@brief Report on the input parameters and the found hydrophobic regions
void SecretionOptimizationMover::report_test( std::ostream & out, core::pose::Pose & pose )
{
	using namespace core::scoring::methods;
	using ObjexxFCL::format::F;
	using ObjexxFCL::format::LJ;
	out << "\n--------------------------------------\n" ;
	out << " You set:                             \n" ;
	out << " Max mutations " << max_mutations_ << '\n' ;
	out << " dG threshold " << dG_ins_threshold_<< '\n';
	out << " designing to " << keep_ << '\n';
	out << " Score tolerance " << score_tolerance_ << '\n' ;
	out << "--------------------------------------" <<std::endl ;
	found_regions_ = SecretionOptimizationMover::find_tm_regions(TR, pose);
	out << " hooray! " << std::endl;
}

//@brief This identifies TM regions in the pose
utility::vector1<SecretionOptimizationMover::tm_region>
SecretionOptimizationMover::find_tm_regions( std::ostream & out, core::pose::Pose & pose)
{
	//it'll return output, which is a subset of all of the tm energies in the pose
	utility::vector1<tm_region> output;
	//energies will have ALL of the tm regions in the sequence
	utility::vector1<tm_region> energies;
	//there was a SymmetryInfoCOP here but it is no longer necessary
	utility::vector1<core::Size> chain_resnums = core::pose::get_resnums_for_chain(pose,chain_to_degrease_);
	core::Size scan_start = chain_resnums.front();
	core::Size scan_end = chain_resnums.back() - window_size_ + 1;
	out << "Degreasing chain " << chain_to_degrease_ << " from indices " << std::to_string(scan_start) << " to " << std::to_string(scan_end) << std::endl;

	//i'm not good at c++ so i instantiate a temporary tm_region struct just to make life a little easier
	tm_region temp;
	core::SSize temp_counter = 1;
	for ( core::Size q = scan_start; q <= scan_end; q++ ) {
		temp.index = q;
		temp.sequence = get_window_at_index(pose,window_size_,q);
		temp.dG_ins = dG_ins_for_window(temp.sequence);
		energies.push_back({temp.index, temp.sequence, temp.dG_ins, 0, temp.possible_mutants});
		out << " Sequence " << energies[temp_counter].sequence << " at index " << energies[temp_counter].index <<  " has dG_ins_pred of: " << energies[temp_counter].dG_ins << std::endl;
		temp_counter++;
		//TR << " the coefficient a0 for Ala is " << coeffs["A"]["a0"];
	}

	utility::vector1<tm_region> energies_sorted = energies;
	std::sort(energies_sorted.begin(), energies_sorted.end());

	temp_counter = 1;
	for ( core::Size q = scan_start; q <= scan_end; q++ ) {//i just scan the top half of the scoring regions for local minima. there might be a better way to do this.
		core::SSize i = temp_counter;
		core::Size iplus = i + 1;
		if ( iplus > energies_sorted.size() ) {
			iplus = energies_sorted.size();
		}
		core::SSize iminus = i - 1;
		if ( iminus < 1 ) {
			iminus = 1;
		}
		core::Size iplus2 = i + 2;
		if ( iplus2 > energies_sorted.size() ) {
			iplus2 = energies_sorted.size();
		}
		core::SSize iminus2 = i - 2;
		if ( iminus2 < 1 ) {
			iminus2 = 1;
		}
		out << "i: " << i << ",iplus: " << iplus << ",iminus: " << iminus << ",iplus2: " << iplus2 << ",iminus2: " << iminus2;
		if ( energies_sorted[iplus].index > scan_end ) {
			iplus = i;
		}
		if ( energies_sorted[iplus2].index > scan_end ) {
			iplus2 = i;
		}
		if ( energies_sorted[iminus].index < scan_start ) {
			iminus = i;
		}
		if ( energies_sorted[iminus2].index < scan_start ) {
			iminus2 = i;
		}
		out << "i_indices is sorted out at " << i;
		if ( energies_sorted[i].dG_ins <= energies_sorted[iplus].dG_ins && energies_sorted[i].dG_ins <= energies_sorted[iminus].dG_ins && energies_sorted[i].dG_ins <= dG_ins_threshold_ && energies_sorted[i].dG_ins <= energies_sorted[iplus2].dG_ins && energies_sorted[i].dG_ins <= energies_sorted[iminus2].dG_ins && energies_sorted[i].dG_ins <= dG_ins_threshold_ ) {
			out << "Found a local minimum at index " << energies_sorted[i].index << " with sequence "<< energies_sorted[i].sequence <<  " with dG_ins_pred of: " << energies_sorted[i].dG_ins << std::endl;
			for ( core::Size j = energies_sorted[i].index; j < energies_sorted[i].index + window_size_ + 1 ; ++j ) {
				energies_sorted[i].possible_mutants.push_back({j, keep_, 0, 0});
			}
			output.push_back({energies_sorted[i].index,energies_sorted[i].sequence,energies_sorted[i].dG_ins, 0, energies_sorted[i].possible_mutants});
		}  //this is a little complicated. energies_sorted and output are vectors of tm_region.
		//tm_region also contains vectors of mutts as their possible_mutants.
		//for now, just fill in keep_ (designable residues) for the possible mutants of each output tm_region.
		temp_counter++;
	}
	return output;
}

protocols::moves::MoverOP
SecretionOptimizationMover::clone() const
{
	return (protocols::moves::MoverOP) protocols::moves::MoverOP( new SecretionOptimizationMover( *this ) );
}

protocols::moves::MoverOP
SecretionOptimizationMover::fresh_instance() const
{
	return (protocols::moves::MoverOP) protocols::moves::MoverOP( new SecretionOptimizationMover );
}

//utility::tag::XMLSchemaComplexTypeGeneratorOP
//SecretionOptimizationMover::define_schema()
void SecretionOptimizationMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::rosetta_scripts::attributes_for_parse_task_operations(attlist);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"max_mutations", xsct_non_negative_integer,
		"maximum number of point mutants to jam into each particular tm region. it should take the X highest mutants in terms of ddG_ins effect.",
		"3");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"aas_allowed", xs_string,
		"one-letter amino acid specification (non-delimited) for allowable residues",
		"DEKRQN");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"chain_to_degrease", xs_string,
		"single-character chain to degrease, defaulted to chain A",
		"A");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"dG_ins_threshold", xs_decimal,
		"What counts as a tm region? +2.7 seems high, but some of these problem regions seem to have dG_ins greater than 2.0. 2.7 was determined kind of empirically.",
		"2.7");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"dG_tolerance", xs_decimal,
		"how much of an effect does a particular mutant need to have? 0.27 was determined as a sort of happy compromise.",
		"0.27");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"score_tolerance", xs_decimal,
		"how much can the total score of the pose change due to a mutation? I actually found that many IMPROVE the score of poses, so it's set to 0 now. +1.0 only gave me 20% more hits, +3.0 only 40%. This can also depend on your application. Set it to something really high to get as many outputs as possible. +45 REU seems to be tolerable for binder design purposes.",
		"0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"dump_multis", xsct_rosetta_bool,
		"if you want the PDBs for each of the combo mutants, set this to true. writes directly to disk, and so avoid if running large parallel jobs.",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"dump_singles", xsct_rosetta_bool,
		"if you want every single point mutant along the way, set this to true. writes directly to disk, and so avoid if running large parallel jobs.",
		"false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"dump_dir", xs_string,
		"directory into which your dumped pdbs are dumped. by default it's where you run your script (maybe not a good idea ?)",
		"."
	);

	attlist + XMLSchemaAttribute::attribute_w_default(
		"repack_shell", xs_decimal,
		"distance around designable residues to repack. 8.0 seems like a standard distance.",
		"8.0");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"largest_ddg", xsct_rosetta_bool,
		"output the pose with the highest ddG_ins - enabled by default",
		"true");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"smallest_dscore", xsct_rosetta_bool,
		"output the pose with the lowest score - dsiabled by default",
		"false");
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This is the SecretionOptimizationMover, aka the Degreaser. It uses a model from Hessa et al. 2007 to identify regions of low transmembrane insertion potential (dG_ins,pred, or just dG in this Mover) and design them away while maintaining structural stability near the mutation site. This Mover is semi-exhaustive and allows the user to set parameters such as how much dG_ins,pred change is necessary or how much score difference is tolerable. With all default options, a maximum of three mutations is made to a pose in the region of lowest initial dG_ins,pred such that the change in dG_ins,pred is largest. Other options enable more complex usage of this Mover. See preprint with experimental results available at https://www.biorxiv.org/content/10.1101/2022.08.04.502842v1", attlist);
}

void SecretionOptimizationMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SecretionOptimizationMover::provide_xml_schema(xsd);
}

//what's in a name?
std::string
SecretionOptimizationMover::get_name() const {
	return "SecretionOptimizationMover";
}

std::string SecretionOptimizationMover::mover_name() {
	return "SecretionOptimizationMover";
}

std::string SecretionOptimizationMoverCreator::keyname() const {
	return SecretionOptimizationMover::mover_name();
}

protocols::moves::MoverOP
SecretionOptimizationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SecretionOptimizationMover);
}
//for unit testing
void SecretionOptimizationMover::test_move( core::pose::Pose & pose) {
	SecretionOptimizationMover::apply(pose);
}
} // movers
} // protein_interface_design
} // protocols

