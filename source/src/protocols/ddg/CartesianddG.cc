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
/// @author Brandon Frenz
/// @author Frank DiMaio
/// @author Hahnbeom Park

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>

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
#include <numeric/xyz.functions.hh>

#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <protocols/ddg/ddGMover.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/constraint_movers/AddConstraintsToCurrentConformationMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/hybridization/util.hh>
#include <protocols/hybridization/CartesianSampler.hh>
#include <protocols/ddg/CartesianddG.hh>

//Auto Headers
#include <utility/json_utilities.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <basic/Tracer.hh>

//Auto Headers
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ddg {
namespace CartesianddG {

using namespace core;

static basic::Tracer TR("protocols.CartesianddG");

#ifdef _NLOHMANN_JSON_ENABLED_

utility::vector1<core::Size>
MutationSet::get_prolines(){
	utility::vector1<core::Size> prolines;
	for ( core::Size i=1; i<=mutations_.size(); i++ ) {
		if ( mutations_[i] == core::chemical::aa_pro ) {
			prolines.push_back(resnums_[i]);
		}
	}
	return prolines;
}

utility::vector1<std::pair<core::Size,std::string>>
MutationSet::get_mutation_pairs(){
	runtime_assert(mutations_.size() == resnums_.size());
	utility::vector1<std::pair<core::Size,std::string>> mutpairs;
	for ( core::Size i=1; i<=resnums_.size(); i++ ) {
		std::pair<core::Size,std::string> resmut = std::make_pair(resnums_[i],core::chemical::name_from_aa(mutations_[i]));
		mutpairs.push_back(resmut);
	}
	return mutpairs;
}

nlohmann::json
MutationSet::to_json(core::pose::Pose & pose){
    add_wildtypes(pose);
    nlohmann::json mutationset;
    utility::vector1<nlohmann::json> mutations;
    for( core::Size i=1; i<=resnums_.size(); i++){
        nlohmann::json mutdata;
        char mut = core::chemical::oneletter_code_from_aa(mutations_[i]);
        char wt = wild_types_[i];
        core::Size pos = resnums_[i];
        mutdata["mut"] = utility::to_string(mut);
        mutdata["wt"] = utility::to_string(wt);
        mutdata["pos"] = pos;
        mutations.push_back(mutdata);
    }
    return mutations;
}

//bool
//MutationSet::swaps_charge(utility::vector1<core::chemical::AA> native_aas){
//  //load residue classes
//  utility::vector1<core::chemial::AA> positive;
//  utility::vector1<core::chemial::AA> negative;
//  positive.push_back(core::chemical::lys);
//  positive.push_back(core::chemical::arg);
//  negative.push_back(core::chemical::asp);
//  negative.push_back(core::chemical::glu);
//  runtime_assert(native_aas.size() == mutations_.size());
//  for( core::Size i=1; i<=native_aas.size(); i++){
//
//  }
//
//
//}

bool
MutationSet::is_converged(
	const core::Size n_to_check,
	const core::Real cutoff)
{
	std::sort(scores_.begin(), scores_.end());
	core::Size n_converged = 1;
	if ( scores_.size() < n_to_check ) {
		return false;
	}
	core::Real best_score = scores_[1];
	for ( core::Size i=2; i<=scores_.size(); i++ ) {
		if ( std::abs(best_score-scores_[i]) < cutoff ) n_converged++;
		if ( n_to_check == n_converged ) {
			return true;
		}
	}
	return false;
}

std::string
MutationSet::generate_tag(){
	//generate tag
	std::ostringstream tag;
	tag << "MUT";
	for ( core::Size i=1; i <= resnums_.size(); i++ ) {
		tag << "_" << resnums_[i] << mutations_[i];
	}
	return tag.str();
}

utility::vector1<core::Size>
find_neighbors(
	MutationSet mutations,
	core::pose::Pose const & pose,
	Real const heavyatom_distance_threshold )
{
	utility::vector1<core::Size> neighbors;

	for ( core::Size i : mutations.get_resnums() ) {
		core::conformation::Residue const & rsd1( pose.residue(i) );
		for ( core::Size j=1; j<=pose.size(); j++ ) {
			if ( j == i or mutations.get_resnums().has_value(j) ) continue;
			core::conformation::Residue const & rsd2( pose.residue(j) );
			if ( rsd1.nbr_atom_xyz().distance_squared( rsd2.nbr_atom_xyz() ) <=
					numeric::square( rsd1.nbr_radius() + rsd2.nbr_radius() + heavyatom_distance_threshold ) ) {
				if ( !neighbors.has_value(j) ) {
					neighbors.push_back(j);
				}
			}
		}
	}
	return neighbors;
}

utility::vector1<core::Size>
find_neighbors_directional(
	MutationSet mutations,
	core::pose::Pose const & pose,
	core::Real const K = 8.0
) {

	// interface parameters
	//   if (angle_degrees > K*exp(b*distance_angstrom) we have an interface interaction
	//   reduce K to increase # detected interface residues
	//      max dist = (1/b)*ln(180/k)
	//   at K= 10, maxdist = 10.32
	//         12             9.67
	//         14             9.12
	//         16             8.64
	core::Real b=0.28;

	// make pose polyA
	core::pose::Pose pose_working = pose;
	utility::vector1< Size > protein_residues;
	for ( core::Size i=1; i<=pose.size(); i++ ) {
		if ( pose.residue( i ).is_protein() ) {
			protein_residues.push_back( i );
		}
	}
	protocols::toolbox::pose_manipulation::construct_poly_XXX_pose( "ALA", pose_working, protein_residues, false, false, false );

	utility::vector1<core::Size> neighbors;

	for ( core::Size resnum : mutations.get_resnums() ) {
		neighbors.push_back(resnum);
		conformation::Residue const & rsd1( pose_working.residue( resnum ) );
		for ( core::Size i=1; i<=pose.size(); i++ ) {
			conformation::Residue const & rsd2( pose_working.residue( i ) );
			core::Real dist = (rsd1.atom("CB").xyz() - rsd2.atom("CB").xyz()).length();
			core::Real angle1 = numeric::angle_degrees(rsd1.atom("CA").xyz(), rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz() ) ;
			core::Real angle2 = numeric::angle_degrees(rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz(), rsd2.atom("CA").xyz() ) ;

			core::Real angle_tgt = K*exp(b*dist);

			if ( angle_tgt < 180 && angle1 > angle_tgt && angle2 > angle_tgt ) {
				neighbors.push_back(i);
			}
		}
	}
	return neighbors;

}


/// @brief The input file is a list of mutation blocks.  Usually, this will be a set of point mutations.
/// where each "block" lists a single mutation.  However, it is possible to specify multiple mutations
/// together in a single block.
///
/// The file format is:
/// "total N"
/// followed by N blocks, where each block is
/// "M"
/// specifying followed by M lines of wt/resid/mutaa triples
/// "wtaa resid mutaa"
/// N, M and resid are all supposed to be integers.
/// wtaa, and mutaa are supposed to be 1-letter amino acid codes.
utility::vector1<MutationSet>
read_in_mutations(
	const std::string filename,
	core::pose::Pose & pose)
{
	core::Size const n_iters = basic::options::option[ basic::options::OptionKeys::ddg::iterations ].value();
	std::ifstream inputstream;
	inputstream.open(filename.c_str());
	utility::vector1<MutationSet> mutation_sets;
	if ( inputstream.is_open() ) {
		int total;
		std::string total_keyword;
		inputstream >> total_keyword;
		debug_assert(total_keyword.compare("total") == 0);

		inputstream >> total; //keep for cross-checking
		TR << " total read as " << total << std::endl;
		while ( !inputstream.eof() && total>0 ) {
			utility::vector1<core::Size> resnums;
			utility::vector1<core::chemical::AA> mutations;
			int num_mutations;
			inputstream >> num_mutations;
			runtime_assert(num_mutations>0);
			while ( num_mutations > 0 ) {
				char wt; int resnum; char mut;
				inputstream >> wt >> resnum >> mut;
				TR.Debug << "wt is " << wt << " resnum is " << resnum << " and mut is " << mut << std::endl;
				runtime_assert( pose.residue(resnum).name1() == wt ); /// APL -- never use regular asserts when it comes to user input
				runtime_assert( core::chemical::oneletter_code_specifies_aa( mut ) ); /// APL -- input should specify the 1-letter code for an amino acid.
				core::chemical::AA mutation = core::chemical::aa_from_oneletter_code( mut );
				resnums.push_back(resnum);
				mutations.push_back(mutation);
				num_mutations--;
			}
			total--;
			MutationSet current_set(resnums,mutations,n_iters);
			if ( num_mutations < 0 ) {
				TR.Error << "number of mutations mismatch! num_mutations < 0" << std::endl;
				runtime_assert(num_mutations==0);
			} else {
				mutation_sets.push_back(current_set);
			}
		}
		if ( total < 0 ) {
			TR.Error << "total number of mutations mismatch! total < 0" << std::endl;
			runtime_assert(total>0);
		}
	}
	return mutation_sets;
}

void
read_existing_json(utility::vector1<MutationSet> & existing_mutsets, const std::string filename, nlohmann::json & results_json, const core::Size iters){
    if( !utility::file::file_exists(filename) ) return;
    std::ifstream i(filename);
    i >> results_json;
    for( nlohmann::json mutationset : results_json ){
       utility::vector1<core::Size> resnums;
       utility::vector1<core::chemical::AA> mutations;
       for( nlohmann::json mutation : mutationset["mutations"] ){
            if( mutation["mut"] == "WT" ){
                continue;
            }
            core::Size pos = mutation["pos"];
            resnums.push_back(pos);
            std::string mutc = mutation["mut"].get<std::string>();
            //mutation["mut"].get<std::string>().c_str()(1,mutc);
		    core::chemical::AA mut = core::chemical::aa_from_oneletter_code( mutc[0] );
            mutations.push_back(mut);
       }
       MutationSet mutset(resnums,mutations,iters);
       for( MutationSet & existing_mutset : existing_mutsets ){
           if( existing_mutset.is_same(mutset) ){
               existing_mutset.subtract_iterations(1);
               core::Real score = mutationset["scores"]["total"];
               existing_mutset.add_score(score);
           }
       }
    }
}

void
sample_fragments(
	core::pose::Pose & pose,
	MutationSet & mutations,
	core::scoring::ScoreFunctionOP sf,
	const core::Size bbnbrs,
	const core::Size ncycles)
{

	core::pose::Pose best_pose = pose;
	core::pose::Pose work_pose = pose;
	core::Real best_score = 1e9;
	core::scoring::ScoreFunctionOP cen_sf = core::scoring::get_score_function(false);

	for ( std::pair<core::Size, core::fragment::FragSetOP> fragment : mutations.get_fragments() ) {
		work_pose = best_pose;
		utility::vector1<core::fragment::FragSetOP> cart_frags;
		cart_frags.push_back(fragment.second);

		//Setup cart samplers
		protocols::hybridization::CartesianSampler cart_sampler( cart_frags );
		cart_sampler.set_scorefunction(cen_sf);
		cart_sampler.set_fa_scorefunction(sf);
		cart_sampler.set_strategy("user");
		core::select::residue_selector::ResidueSelectorOP selector( new core::select::residue_selector::ResidueIndexSelector( fragment.first ) );
		cart_sampler.set_userpos(selector);

		//core::Size tag = 1;
		for ( core::fragment::ConstFrameIterator iter = fragment.second->begin(); iter != fragment.second->end(); ++iter ) {
			work_pose = best_pose;
			core::fragment::Frame frame = **iter;
			frame.shift_to( fragment.first );
			utility::vector1<core::Size> frag_resis;
			core::Size lower = std::max(core::Size(1),fragment.first-bbnbrs);
			core::Size upper = std::min(fragment.first+frame.length()+bbnbrs,pose.size());
			for ( core::Size i=lower; i<=upper; i++ ) {
				frag_resis.push_back(i);
			}
			for ( core::Size i=1; i<=ncycles; i++ ) {
				cart_sampler.apply_frame(work_pose,frame);
				min_pack_min_element(work_pose,frag_resis,sf);
				core::Real score = (*sf)(work_pose);
				if ( score < best_score ) {
					best_score = score;
					best_pose = work_pose;
				}
			}
		}
	}
	const core::Real original_score = (*sf)(pose);
	if ( original_score > best_score ) {
		TR << " accepting sampler best score is " << best_score << std::endl;
		pose = best_pose;
	}
}

void
min_pack_min_element(
	core::pose::Pose & pose,
	utility::vector1<core::Size> min_resis,
	core::scoring::ScoreFunctionOP sfxn)
{

	core::kinematics::MoveMap mm;
	for ( core::Size i : min_resis ) {
		mm.set_chi( i, true );
		mm.set_bb( i, true );
	}

	// setup minimizer
	core::optimization::CartesianMinimizerOP minimizer( new core::optimization::CartesianMinimizer );
	core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 1e-4, true, false, false);
	options.max_iter(200);
	//run
	minimizer->run(pose,mm,*sfxn,options);

	pack::task::PackerTaskOP packer_task(pack::task::TaskFactory::create_packer_task(pose));
	packer_task->restrict_to_repacking();

	//Block repacking on unchanged residues
	for ( core::Size i=1; i<=pose.size(); i++ ) {
		if ( !min_resis.has_value(i) ) {
			packer_task->nonconst_residue_task(i).prevent_repacking();
		}
	}
	//Repack
	protocols::minimization_packing::PackRotamersMoverOP packer(new protocols::minimization_packing::PackRotamersMover( sfxn, packer_task ));
	packer->apply(pose);
	//Minimize
	minimizer->run(pose,mm,*sfxn,options);

}

void
optimize_structure(
	MutationSet mutations,
	core::scoring::ScoreFunctionOP fa_scorefxn,
	core::pose::Pose & pose,
	utility::vector1<core::Size> neighbors,
	const bool flex_bb,
	const bool cartesian,
	const core::Size bbnbrs)
{

	protocols::relax::FastRelax fastrelax( fa_scorefxn, 0 );
	fastrelax.cartesian( cartesian );

	core::kinematics::MoveMapOP movemap(new core::kinematics::MoveMap);
	movemap->set_bb( false );
	movemap->set_chi( false );
	movemap->set_jump( false );

	for ( core::Size i : neighbors ) {
		movemap->set_chi( i, true );
		if ( flex_bb ) {
			movemap->set_bb( i, true );
		}
	}
	for ( core::Size i : mutations.get_resnums() ) {
		for ( core::Size j=std::max(i-bbnbrs,core::Size(1)); j<=std::min(i+bbnbrs,pose.size()); j++ ) {
			movemap->set_bb( i, true );
		}
	}

	fastrelax.set_movemap( movemap );
	fastrelax.apply(pose);
	pose.remove_constraints();

}

void
optimize_native(
	utility::vector1<MutationSet> mutationsets,
	core::pose::Pose & pose,
	core::scoring::ScoreFunctionOP fa_scorefxn,
	const core::Size bbnbrs,
	const bool cartesian)
{
	utility::vector1<core::Size> all_muts;
	for ( MutationSet mutset : mutationsets ) {
		for ( core::Size i : mutset.get_resnums() ) {
			all_muts.push_back(i);
		}
	}
	utility::vector1<core::chemical::AA> alanines;
	core::chemical::AA ala = core::chemical::aa_from_oneletter_code( 'A' );
	for ( core::Size i=1; i<=all_muts.size(); i++ ) {
		alanines.push_back(ala);
	}

	core::Real cutoff = 6.0;
	core::Size iterations = 1; //Doesn't get used
	MutationSet newmuts(all_muts,alanines,iterations);
	utility::vector1<core::Size> neighbors = find_neighbors(newmuts,pose,cutoff);
	optimize_structure(newmuts,fa_scorefxn,pose,neighbors,bbnbrs,cartesian);
}

void
extract_element(
	core::pose::Pose& context_pose,
	core::pose::Pose &peptide_pose,
	MutationSet mutations,
	const core::Size neighbors_to_extract)
{

	utility::vector1<core::Size> peptide_resis;

	for ( core::Size i : mutations.get_resnums() ) {
		core::Size lower = std::max(i-neighbors_to_extract,core::Size(1));
		core::Size upper = std::min(i+neighbors_to_extract,context_pose.size());
		for ( core::Size j=lower; j<=upper; j++ ) {
			peptide_resis.push_back(j);
		}
	}

	if ( peptide_resis.size() == 0 ) return;
	for ( core::Size i=context_pose.size(); i>=1; i-- ) {
		if ( peptide_resis.has_value(i) ) {
			context_pose.delete_residue_slow(i);
		}
	}
	for ( core::Size i=peptide_pose.size(); i>=1; i-- ) {
		if ( !peptide_resis.has_value(i) ) {
			peptide_pose.delete_residue_slow(i);
		}
	}
}

core::Size
find_interface_jump(core::pose::Pose & pose, core::Size interface_ddg){

	//fd try to be smart .. look for interchain jump
	if ( interface_ddg > pose.num_jump() ) {
		for ( int i=1; i<=(int)pose.fold_tree().num_jump(); ++i ) {
			kinematics::Edge jump_i = pose.fold_tree().jump_edge( i ) ;
			if ( pose.pdb_info()->chain(jump_i.start()) != pose.pdb_info()->chain(jump_i.stop()) ) {
				runtime_assert( interface_ddg > pose.num_jump() );
				interface_ddg = i;
			}
		}
		TR << "Autosetting interface jump to " << interface_ddg << std::endl;
	}
	runtime_assert ( interface_ddg <= pose.num_jump() );

	return interface_ddg;
}

void
mutate_sequence( std::string & sequence, MutationSet mutations){

	utility::vector1<core::Size> resnums = mutations.get_resnums();
	utility::vector1<core::chemical::AA> mutated_aas = mutations.get_mutations();
	for ( core::Size i=1; i<=resnums.size(); i++ ) {
		sequence[resnums[i]-1] = core::chemical::oneletter_code_from_aa(mutated_aas[i]);
	}
}

utility::vector1<core::Size>
involves_prolines(
	core::pose::Pose & pose,
	MutationSet mutations)
{
	utility::vector1<core::Size> prolines = mutations.get_prolines();
	for ( core::Size resnum : mutations.get_resnums() ) {
		if ( pose.residue(resnum).name3() == "PRO" ) {
			prolines.push_back(resnum);
		}
	}
	return prolines;
}

void
mutate_pose(
	core::pose::Pose & pose,
	MutationSet mutations,
	core::scoring::ScoreFunctionOP score_fxn)
{

	pack::task::PackerTaskOP packer_task(pack::task::TaskFactory::create_packer_task(pose));
	for ( core::Size i=1; i<=mutations.get_resnums().size(); i++ ) {
		utility::vector1< bool > allowable( 20, false );
		core::Size resnum = mutations.get_resnum(i);
		core::chemical::AA aa = mutations.get_aa(i);
		allowable[int(aa)] = true;
		packer_task->nonconst_residue_task(resnum).restrict_absent_canonical_aas(allowable);
	}
	//Block repacking on unchanged residues
	for ( core::Size i=1; i<=pose.size(); i++ ) {
		if ( !mutations.get_resnums().has_value(i) ) {
			packer_task->nonconst_residue_task(i).prevent_repacking();
		}
	}
	//Repack
	protocols::minimization_packing::PackRotamersMoverOP packer(new protocols::minimization_packing::PackRotamersMover( score_fxn, packer_task ));
	packer->apply(pose);
}

void
pick_fragments(
	core::pose::Pose & pose,
	utility::vector1<MutationSet> & mutationsets,
	const core::Size frag_nbrs)
{
	if ( pose.size() < 4 ) return;//We cannot do fragment insertion for poses smaller than 4 residues.

	for ( MutationSet& mutations : mutationsets ) {
		utility::vector1<core::Size> prolines = involves_prolines(pose,mutations);
		if ( prolines.size() > 0 ) {
			std::string sequence = pose.sequence();
			mutate_sequence(sequence,mutations);
			for ( int proline : prolines ) {
				core::Size lower = core::Size(std::max(int(proline-frag_nbrs),int(1)));
				core::Size upper = std::min(proline+frag_nbrs,pose.size());
				core::Size fragsize = std::max(upper-lower+1,core::Size(4)); //Fragments of at least length 4 are required
				if ( lower+fragsize > pose.size() ) lower = pose.size()-fragsize+1;
				std::string subseq = sequence.substr(lower-1,fragsize);
				core::fragment::FragSetOP frags = protocols::hybridization::create_fragment_set_no_ssbias(subseq,
					fragsize, basic::options::option[basic::options::OptionKeys::ddg::nfrags].value(), 'D');
				mutations.add_fragments(lower,frags);
			}
		}
	}
}

void subtract_iterations(
	utility::vector1<MutationSet> & mutationsets,
	utility::vector1<std::pair<core::Size,std::string>> finished_residues,
	const core::Real score)
{
	for ( MutationSet & mutations: mutationsets ) {
		utility::vector1<std::pair<core::Size,std::string>> mutpairs = mutations.get_mutation_pairs();
		std::sort(finished_residues.begin(), finished_residues.end());
		std::sort(mutpairs.begin(), mutpairs.end());
		if ( finished_residues == mutpairs ) {
			mutations.subtract_iterations(1);
			mutations.add_score(score);
			TR << " adding score " << score << std::endl;
		}
	}
}

void
read_existing(
 const std::string filename,
 utility::vector1<MutationSet> & mutationsets)
{

	//output file, mark if rerun
	std::map< std::string, Size > before_jump_done_list;
	std::map< std::string, Size > after_jump_done_list;
	utility::io::ozstream ofp;
	utility::file::FileName parsefn(filename);
	std::string ofn = parsefn.base()+".ddg";
	if ( utility::file::file_exists(ofn) ) {
		utility::io::izstream ifp;
		ifp.open(ofn);
		if ( !ifp.good() ) utility_exit_with_message( "can't open cluster file " + ofn );
		while ( !ifp.eof() ) {
			std::string line;
			utility::io::getline( ifp, line );
			//parse the line
			std::istringstream istr(line);
			std::string cat, rd, mut;
			core::Real score;
			istr >> cat >> rd >> mut >> score;
			mut = mut.substr(0, mut.size()-1); //Remove ending colon.

			//skip if no keyword
			const size_t pos1 = rd.find("Round");
			const size_t pos2 = rd.find(':');
			if ( std::string::npos == pos1 || std::string::npos == pos2 ) continue;
			//std::string nrdstr = rd.substr(pos1+5, pos2-5);
			//Size nrd = std::atoi(nrdstr.c_str());

			char delimiter = '_';
			size_t pos = 0;
			std::string aa_resnum;
			utility::vector1<std::pair<core::Size,std::string>> mutated_resis;
			bool done = false;
			if ( mut.find(delimiter) == std::string::npos ) done = true;
			while ( !done ) {
				if ( (pos = mut.find(delimiter)) == std::string::npos ) done = true;
				aa_resnum = mut.substr(0,pos);
				mut.erase(0, pos+1);
				if ( aa_resnum == "MUT" ) continue;
				std::string name3 = aa_resnum.substr(aa_resnum.length()-3,aa_resnum.length());
				std::string str_num = aa_resnum.substr(0,aa_resnum.length()-3);
				core::Size resnum;
				resnum = std::stoi(str_num);
				std::pair<core::Size,std::string> pair = std::make_pair(resnum,name3);
				mutated_resis.push_back(pair);
			}
			subtract_iterations(mutationsets,mutated_resis,score);
		}
		ifp.close();
	}
}

core::Real
extracted_score(
	core::pose::Pose & pose,
	MutationSet mutations,
	core::scoring::ScoreFunctionOP score_fxn,
	const core::Size n_nbrs_to_extract)
{

	//Declare score types to be used in binding calculations
	utility::vector1<core::scoring::ScoreType> binding_scoretypes;
	binding_scoretypes.push_back(core::scoring::fa_elec);
	binding_scoretypes.push_back(core::scoring::fa_sol);
	binding_scoretypes.push_back(core::scoring::hbond_sr_bb);
	binding_scoretypes.push_back(core::scoring::hbond_lr_bb);
	binding_scoretypes.push_back(core::scoring::hbond_sc);

	//Setup score functions
	core::scoring::ScoreFunctionOP local_sfxn = score_fxn->clone();
	for ( core::scoring::ScoreType st : binding_scoretypes ) {
		local_sfxn->set_weight(st,0);
	}

	core::scoring::ScoreFunctionOP binding_sfxn = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
	for ( core::scoring::ScoreType st : binding_scoretypes ) {
		binding_sfxn->set_weight(st,score_fxn->get_weight(st));
	}

	//Extract elements
	core::pose::Pose context = pose;
	core::pose::Pose peptide = pose;
	extract_element(context,peptide,mutations,n_nbrs_to_extract);

	//calculate score
	core::Real base_score = (*local_sfxn)(pose);
	core::Real total_binding = (*binding_sfxn)(pose);
	core::Real context_binding = (*binding_sfxn)(context);
	core::Real peptide_binding = (*binding_sfxn)(peptide);
	core::Real binding_energy = total_binding-context_binding-peptide_binding;
	core::Real final_score = binding_energy + base_score;
	return final_score;

}

void
write_json(const std::string filename, nlohmann::json results_json){
	utility::io::ozstream ofp;
    ofp.open(filename);
    ofp << std::setw(4) << results_json << std::endl;
    ofp.close();

}

nlohmann::json
get_scores_as_json(
    core::pose::Pose & pose,
    core::scoring::ScoreFunctionOP score_fxn,
    core::Real total_score)
{

    nlohmann::json scores;
    for( core::scoring::ScoreType st : score_fxn->get_nonzero_weighted_scoretypes() ){
        std::string st_name = core::scoring::name_from_score_type(st);
        core::Real score = pose.energies().total_energies()[st];
        scores[st_name] = score;
    }
    scores["total"] = total_score;
    return scores;
}

void
run(core::pose::Pose & pose){

	//Json is required to use this protocol.
    nlohmann::json results_json;
    utility::vector1<nlohmann::json> mutset;
    results_json = mutset;

	bool fd_mode = basic::options::option[ basic::options::OptionKeys::ddg::fd_mode ].value();
	const core::Real restricted_neighbors_k = 6.0;
	const core::Real permissive_neighbors_k = 8.0;
	core::Real cutoff = fd_mode? permissive_neighbors_k : restricted_neighbors_k;
	if ( basic::options::option[ basic::options::OptionKeys::ddg::opt_radius].user() ) {
		cutoff = basic::options::option[ basic::options::OptionKeys::ddg::opt_radius ].value();
	}
	//Parse command line option
	Size bbnbrs = basic::options::option[ basic::options::OptionKeys::ddg::bbnbrs ].value();
	core::Size frag_nbrs = basic::options::option[ basic::options::OptionKeys::ddg::frag_nbrs ].value();
	core::Size ncycles = basic::options::option[ basic::options::OptionKeys::ddg::ntrials ].value();
	bool force_iterations = basic::options::option[ basic::options::OptionKeys::ddg::force_iterations ].value();
	bool flex_bb = basic::options::option[ basic::options::OptionKeys::ddg::flex_bb ].value();
	bool cartesian = basic::options::option[ basic::options::OptionKeys::ddg::cartesian ].value();

	core::Size extract_nbrs = basic::options::option[ basic::options::OptionKeys::ddg::extract_element_nbrs].value();

	if( basic::options::option[basic::options::OptionKeys::ddg::optimize_proline].value() == true ){
	TR << " -frag_nbrs must be assigned in order to optimize proline using fragments. 4 is a good value." << std::endl;
	runtime_assert( basic::options::option[basic::options::OptionKeys::ddg::frag_nbrs].value() != 0 );
	}

	core::scoring::ScoreFunctionOP score_fxn;
	score_fxn = core::scoring::get_score_function();

	core::Size interface_ddg = find_interface_jump(pose, basic::options::option[basic::options::OptionKeys::ddg::interface_ddg].value());
    TR << " interface ddg is " << interface_ddg << std::endl;

	if ( !basic::options::option[ basic::options::OptionKeys::ddg::mut_file ].user() ) {
		utility_exit_with_message("Option -ddg::mut_file must be set.");
	}

	//Read the in the target mutations and the completed jobs.
	TR.Debug << "reading in mutfile" << std::endl;
	std::string filename = basic::options::option[basic::options::OptionKeys::ddg::mut_file].value();
	utility::vector1<MutationSet> mutationsets = read_in_mutations( filename, pose);

	utility::io::ozstream ofp;
	utility::file::FileName parsefn(filename);
	std::string ofn = parsefn.base()+".ddg";
    const std::string jsonout = parsefn.base()+".json";
	if ( utility::file::file_exists(ofn) ) {
		//append in previous file
		ofp.open_append(ofn);
	} else {
		//create new file
		ofp.open(ofn);
	}
    if( basic::options::option[ basic::options::OptionKeys::ddg::json ].value() ){
	    core::Size const n_iters = basic::options::option[ basic::options::OptionKeys::ddg::iterations ].value();
	    read_existing_json( mutationsets, jsonout, results_json, n_iters);
    }else{
        read_existing( filename, mutationsets );
    }


	//Pick Fragments if they are being used to to optimize proline.
	if ( basic::options::option[basic::options::OptionKeys::ddg::optimize_proline].value() ) {
		pick_fragments(pose, mutationsets, frag_nbrs);
	}

	//optimize_native(mutationsets,pose,score_fxn,bbnbrs,cartesian);
	core::Real native_score = (*score_fxn)(pose);
	std::string tag  = "WT";
    nlohmann::json wt_results;
    utility::vector1<nlohmann::json> mutations_json;
    nlohmann::json mutation;
    mutation["mut"] = "WT";
    mutations_json.push_back(mutation);
    wt_results["mutations"] = mutations_json;
    wt_results["scores"] = get_scores_as_json(pose,score_fxn,native_score);
	ofp << "COMPLEX:   Round1: " << tag << ": " << ObjexxFCL::format::F(9,3,native_score) << " "
		<< pose.energies().total_energies().weighted_string_of( score_fxn->weights() ) << std::endl;
    results_json.push_back(wt_results);
    if( basic::options::option[ basic::options::OptionKeys::ddg::json ].value() ){
        write_json(jsonout,results_json);
    }


	for ( MutationSet mutations : mutationsets ) {
		if ( mutations.iterations() == 0 ) continue;
        nlohmann::json mutationset_results;
        nlohmann::json mutations_json = mutations.to_json(pose);
        mutationset_results["mutations"] = mutations_json;
		core::pose::Pose work_pose(pose);
		mutate_pose(work_pose,mutations,score_fxn);
		utility::vector1<core::Size> neighbors;
		if ( fd_mode ) {
			neighbors = find_neighbors_directional(mutations,work_pose,cutoff);
		} else {
			neighbors = find_neighbors(mutations,work_pose,cutoff);
		}

		//This is where the real work gets done
		for ( core::Size i=1; i<= mutations.iterations(); i++ ) {
			if ( !force_iterations ) {
			    core::Size n_converged = basic::options::option[ basic::options::OptionKeys::ddg::n_converged ].value();
			    core::Real score_cutoff = basic::options::option[ basic::options::OptionKeys::ddg::score_cutoff ].value();
				if ( mutations.is_converged(n_converged,score_cutoff) ) {
					TR << " is converged!" << std::endl;
					break;
				}
			}
			core::pose::Pose local_pose = work_pose;
			sample_fragments(local_pose,mutations,score_fxn,bbnbrs,ncycles);
			optimize_structure(mutations,score_fxn,local_pose,neighbors,flex_bb,bbnbrs,cartesian);
			core::Real score = 0;
			core::Size round = basic::options::option[ basic::options::OptionKeys::ddg::iterations ].value()-mutations.iterations()+i;
			if ( extract_nbrs != 0 ) {
				score = extracted_score(local_pose,mutations,score_fxn,extract_nbrs);
				core::Real native_score = extracted_score(pose,mutations,score_fxn,extract_nbrs);
				core::Real ddg = score-native_score;
				ofp << "COMPLEX_DDG:   Round" << utility::to_string(round) << ": " << mutations.generate_tag() << ": " << ObjexxFCL::format::F(9,3,ddg) << std::endl;
			} else {
				score = (*score_fxn)(local_pose);
			}
			mutations.add_score(score);
			//Write Results
			ofp << "COMPLEX:   Round" << utility::to_string(round) << ": " << mutations.generate_tag() << ": " << ObjexxFCL::format::F(9,3,score) << " "
				<< local_pose.energies().total_energies().weighted_string_of( score_fxn->weights() ) << std::endl;
            mutationset_results["scores"] = get_scores_as_json(local_pose,score_fxn,score);
            results_json.push_back(mutationset_results);
            if( basic::options::option[ basic::options::OptionKeys::ddg::json ].value() ){
                write_json(jsonout,results_json);
            }

			//dump pdbs
			if ( basic::options::option[basic::options::OptionKeys::ddg::dump_pdbs].value() ) {
				std::ostringstream dump_fn;
				dump_fn << mutations.generate_tag() << "_bj" << round << ".pdb";
				local_pose.dump_pdb(dump_fn.str());
			}
			if ( interface_ddg > 0 ) {
				Size rb_jump(interface_ddg);
				protocols::rigid::RigidBodyTransMoverOP separate_partners( new protocols::rigid::RigidBodyTransMover( local_pose, rb_jump ) );
				separate_partners->step_size(1000.0);
				separate_partners->apply(local_pose);
                core::Real final_score = (*score_fxn)( local_pose );
				ofp << "APART:     Round" << utility::to_string(round) << ": " << mutations.generate_tag() << ": " << ObjexxFCL::format::F(9,3,final_score) << " "
					<< local_pose.energies().total_energies().weighted_string_of( score_fxn->weights() ) << std::endl;
			}
		}
		//interface mode, seperate and score
		if ( interface_ddg > 0 ) {
			Size rb_jump(interface_ddg);
			protocols::rigid::RigidBodyTransMoverOP separate_partners( new protocols::rigid::RigidBodyTransMover( work_pose, rb_jump ) );
			separate_partners->step_size(1000.0);
			separate_partners->apply(work_pose);

			//repack or not?
			//seperate partners energy
		    for ( core::Size i=1; i<= mutations.iterations(); i++ ) {
                core::pose::Pose local_pose(work_pose);
			    optimize_structure(mutations,score_fxn,local_pose,neighbors,flex_bb,bbnbrs,cartesian);
			    core::Size round = basic::options::option[ basic::options::OptionKeys::ddg::iterations ].value()-mutations.iterations()+i;

				//output
				Real const final_score( (*score_fxn)( local_pose ) );
				if ( basic::options::option[basic::options::OptionKeys::ddg::dump_pdbs].value() ) {
					std::ostringstream dump_fn;
					dump_fn << mutations.generate_tag() << "_aj" << utility::to_string(round) << ".pdb";
					local_pose.dump_pdb(dump_fn.str());
				}
				ofp << "OPT_APART: Round" << utility::to_string(round) << ": " << mutations.generate_tag() << ": " << ObjexxFCL::format::F(9,3,final_score) << " "
					<< local_pose.energies().total_energies().weighted_string_of( score_fxn->weights() ) << std::endl;
			}
		}
	}
	ofp.close();
}

#endif //_NLOHMANN_JSON_ENABLED_
}//CartesianddG
}//ddg
}//protocols
