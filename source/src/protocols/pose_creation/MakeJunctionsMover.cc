// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/pose_creation/MakeJunctionsMover.cc
/// @brief Reads in a file that contains a descriptions of the junctions. Then goes through that file making 1 junction at a time.
/// @detailed
/// @author TJ Brunette (tjbrunette@gmail.com)


#include <protocols/pose_creation/MakeJunctionsMoverCreator.hh>
#include <protocols/pose_creation/MakeJunctionsMover.hh>

#include <protocols/pose_creation/MergePDBatOverlapMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/pose_creation/RepeatPropagationMover.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>


#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/util.hh>
#include <core/select/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/ConstraintSet.hh>


#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

#include <iostream>
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

//#include <basic/datacache/BasicDataCache.hh>
//#include <basic/datacache/CacheableString.hh>
//#include <basic/datacache/DataCache.hh>
//#include <basic/datacache/DataMap.hh>

#include <basic/Tracer.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <map>
#include <queue>

using namespace core;
using namespace protocols::simple_moves;

namespace protocols {
namespace pose_creation {

using utility::vector1;
static basic::Tracer TR( "protocols.pose_creation.MakeJunctionsMover" );


std::string MakeJunctionsMoverCreator::keyname() const
{
	return MakeJunctionsMover::mover_name();
}

protocols::moves::MoverOP
MakeJunctionsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MakeJunctionsMover );
}

void MakeJunctionsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MakeJunctionsMover::provide_xml_schema( xsd );
}


std::string
MakeJunctionsMover::mover_name()
{
	return "MakeJunctionsMover";
}

MakeJunctionsMover::MakeJunctionsMover()
: moves::Mover("MakeJunctionsMover")
{
}

std::string
MakeJunctionsMover::get_name() const {
	return MakeJunctionsMover::mover_name();
}

moves::MoverOP
MakeJunctionsMover::clone() const
{
	return moves::MoverOP( new MakeJunctionsMover( *this ) );
}

moves::MoverOP
MakeJunctionsMover::fresh_instance() const
{
	return moves::MoverOP( new MakeJunctionsMover );
}

MakeJunctionsMover::Design MakeJunctionsMover::line_to_design(std::string line){
	utility::vector1<std::string> jobs = utility::string_split(line,'|',std::string());
	std::stack<std::string> n_term_stack;
	std::stack<std::string> c_term_stack;
	std::stack<std::string> c_term_backward;
	std::string start_struct;
	bool in_c_term=false;
	utility::vector1<std::string> name_plus_junk = utility::string_split(jobs[1],' ',std::string());
	std::string name =  name_plus_junk[1];
	for ( Size ii=2; ii<=jobs.size(); ++ii ) {
		if ( jobs[ii].substr(0,4)=="True" ) {
			start_struct = jobs[ii].substr(0,jobs[ii].size());
			in_c_term = true;
		}
		if ( jobs[ii].substr(0,5)=="False" ) {
			if ( !in_c_term ) {
				n_term_stack.push(jobs[ii]);
			}
			if ( in_c_term ) {
				c_term_backward.push(jobs[ii]);
			}
		}
	}
	while ( !c_term_backward.empty() ) {
		c_term_stack.push(c_term_backward.top());
		c_term_backward.pop();
	}
	MakeJunctionsMover::Design design_tmp(name,n_term_stack,c_term_stack,start_struct);
	return(design_tmp);
}

std::queue<MakeJunctionsMover::Design> MakeJunctionsMover::read_in_designs(){
	std::queue<MakeJunctionsMover::Design> designs;
	utility::io::izstream stream(designs_fn_);
	std::string line;
	while ( getline(stream,line) ) {
		MakeJunctionsMover::Design design_tmp = line_to_design(line);
		designs.push(design_tmp);
	}
	return(designs);
}

void MakeJunctionsMover::parse_attach_description(std::string attach_description,std::string & pdb_location,char & chain, Size & n_term_trim, Size & c_term_trim, Size & n_repeats,Size &n_term_attach_length, Size &c_term_attach_length,std::string & n_term_seq, std::string & c_term_seq){
	//TR << "attach_description" << attach_description << std::endl;
	utility::vector1<std::string> command = utility::string_split(attach_description,',');
	pdb_location = command[2];
	if ( command[3].size() == 0 ) {
		std::stringstream err;
		err << "\nUser must input chain name" << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,err.str());
	} else {
		chain = command[3].at(0);
	}
	n_term_trim = std::stoi(command[4]);
	c_term_trim = std::stoi(command[5]);
	n_repeats = std::stoi(command[6]);
	n_term_attach_length = std::stoi(command[7]);
	c_term_attach_length = std::stoi(command[8]);
	if ( command[9].size()==0 ) {
		n_term_seq = "";
	} else {
		n_term_seq = command[9];
	}
	if ( command[10].size()==0 ) {
		c_term_seq = "";
	} else {
		c_term_seq = command[10];
	}
}



void MakeJunctionsMover::generate_start_pose(core::pose::Pose & pose, core::pose::Pose & background_pose, std::string attach_description){
	Size n_term_trim,c_term_trim,n_term_attach_length,c_term_attach_length,n_repeats;
	std::string pdb_location,n_term_seq,c_term_seq;
	parse_attach_description(attach_description,pdb_location,chain_,n_term_trim,c_term_trim,n_repeats,n_term_attach_length,c_term_attach_length,n_term_seq,c_term_seq);
	core::pose::Pose start_pose = get_and_cache_pdb(pdb_location);
	sfxn_->score(start_pose);
	if ( !has_chain(chain_,start_pose) ) {
		std::stringstream err;
		err << "\nCan't find chain name" << chain_ << "in pdb"  << attach_description << std::endl;
		throw CREATE_EXCEPTION(utility::excn::BadInput,err.str());
	}
	Size chain_id =  get_chain_id_from_chain(chain_,start_pose);
	core::pose::PoseOP desired_chain = start_pose.split_by_chain(chain_id);
	pose = *desired_chain;
	if ( pose.size() == start_pose.size() ) {
		background_pose=*(core::pose::PoseOP( new Pose()));
	} else {  //delete chain from full_pose. Last step will be pasting the original chain back ina
		start_pose.conformation().delete_residue_range_slow( start_pose.conformation().chain_begin( chain_id), start_pose.conformation().chain_end( chain_id ) );
		background_pose = start_pose;
	}
	trim_pose(pose,n_term_trim,c_term_trim);
	if ( n_repeats>0 ) {
		RepeatPropagationMoverOP propagateOP(new RepeatPropagationMover(n_repeats));
		propagateOP->apply(pose);
	}
	sfxn_->score(pose); //necessary to run before assigning sequence in case re-assigning sequence
	assign_seq(pose,n_term_seq,c_term_seq);

}


bool MakeJunctionsMover::attach_next_part(core::pose::Pose & pose, std::string attach_termini, std::string attach_description){
	Size n_term_trim,c_term_trim,n_term_attach_length,c_term_attach_length,n_repeats;
	std::string pdb_location,n_term_seq,c_term_seq;
	parse_attach_description(attach_description,pdb_location,chain_,n_term_trim,c_term_trim,n_repeats,n_term_attach_length,c_term_attach_length,n_term_seq,c_term_seq);
	core::pose::Pose attach_pose = get_and_cache_pdb(pdb_location);
	core::Real pose_score = sfxn_->score(pose);
	core::Real attach_score = sfxn_->score(attach_pose);
	std::cout << "score before attach " << pose_score <<"," << attach_score << std::endl;
	//assumes only chain A in attach pose.
	trim_pose(attach_pose,n_term_trim,c_term_trim);
	if ( n_repeats>0 ) {
		RepeatPropagationMoverOP propagateOP(new RepeatPropagationMover(n_repeats));
		propagateOP->apply(attach_pose);
		//if(relax_during_build_)
	}
	attach_score = sfxn_->score(attach_pose);
	std::cout << "score before assigning seq " << attach_score << std::endl;
	std::cout << "fullatom:" << attach_pose.is_fullatom() <<"," << pose.is_fullatom() << std::endl;
	assign_seq(attach_pose,n_term_seq,c_term_seq);
	attach_score = sfxn_->score(attach_pose);
	std::cout << "score after assigning seq " << attach_score << std::endl;
	core::scoring::ScoreFunctionOP tmp_sfxn = sfxn_->clone();
	MergePDBatOverlapMoverOP mergePDBOP(new MergePDBatOverlapMover(tmp_sfxn));
	Size attachment_length=0;
	if ( attach_termini=="n_term" ) {
		attachment_length=c_term_attach_length;
	} else {
		attachment_length=n_term_attach_length;
	}
	bool success = mergePDBOP->makeJunctions_apply(pose,attach_pose,attachment_length,junction_rmsd_thresh_,attach_termini);
	pose_score = sfxn_->score(pose);
	std::cout << "score after attach" << pose_score << std::endl;
	if ( !success ) {
		junction_failure_set_.insert(attach_description);
	}
	return(success);
}

void MakeJunctionsMover::trim_pose(core::pose::Pose & pose, Size n_term_trim,Size c_term_trim){

	if ( n_term_trim>0 ) {
		//trim residues from n_term
		Size start_location = 1;
		Size end_location = n_term_trim;
		pose.conformation().delete_residue_range_slow(start_location,end_location);
		renumber_pdbinfo_based_on_conf_chains(pose,true,false,false,false);
	}
	if ( c_term_trim>0 ) {
		//delete c_term
		Size start_location = pose.size()-c_term_trim+1;
		Size end_location = pose.size();
		pose.conformation().delete_residue_range_slow(start_location,end_location);
		renumber_pdbinfo_based_on_conf_chains(pose,true,false,false,false);
	}
}

void MakeJunctionsMover::assign_seq(core::pose::Pose & pose, std::string n_term_seq,std::string c_term_seq){
	using namespace core::chemical;
	using namespace core::pack::task;
	std::string pose_seq = pose.sequence();
	utility::vector1< bool > overlap_and_neighbors( pose.size(), false);
	bool residue_mutated = false;
	if ( n_term_seq.size()>0 ) {
		std::string current_n_term_seq = pose_seq.substr(0,n_term_seq.size());
		if ( n_term_seq != current_n_term_seq ) {
			simple_moves::MutateResidueOP mutation_mover;
			for ( Size ii=1; ii<=n_term_seq.size(); ii++ ) {
				AA my_aa =aa_from_oneletter_code(n_term_seq.at(ii-1));
				mutation_mover = simple_moves::MutateResidueOP ( new simple_moves::MutateResidue (
					ii, //position
					my_aa//residue
					) );
				overlap_and_neighbors[ii]=true;
				residue_mutated = true;
				mutation_mover->apply( pose );
			}
		}
	}
	if ( c_term_seq.size()>0 ) {
		std::string current_c_term_seq = pose_seq.substr(pose.size()-c_term_seq.size(),c_term_seq.size());
		if ( c_term_seq != current_c_term_seq ) {
			simple_moves::MutateResidueOP mutation_mover;
			for ( Size ii=1; ii<=c_term_seq.size(); ii++ ) {
				AA my_aa =aa_from_oneletter_code(c_term_seq.at(ii-1));
				mutation_mover = simple_moves::MutateResidueOP ( new simple_moves::MutateResidue (
					pose.size()-c_term_seq.size()+ii, //position
					my_aa//residue
					) );
				overlap_and_neighbors[pose.size()-c_term_seq.size()+ii]=true;
				residue_mutated = true;

				mutation_mover->apply( pose );
			}
		}
	}
	//relax the residues into good positions
	using namespace core::scoring;
	using namespace core::optimization;
	if ( residue_mutated ) {
		Size packing_range = 5;
		core::select::fill_neighbor_residues(pose, overlap_and_neighbors, packing_range);
		optimization::MinimizerOptions minopt( "lbfgs_armijo_nonmonotone", 0.02, true, false, false );
		kinematics::MoveMap mm_loc;
		mm_loc.set_jump( false ); mm_loc.set_bb( false ); mm_loc.set_chi( true );
		TR << "mutations allow residue to relax" << std::endl;
		for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
			if ( overlap_and_neighbors[ resnum ] ) { //should only minimize residues that have been allowed to be minimized.
				TR << resnum << ",";
				mm_loc.set_chi( resnum, true);
			}
		}
		TR << std::endl;
		minopt.max_iter( 100 );
		AtomTreeMinimizer minimizer;
		minimizer.run( pose, mm_loc, *sfxn_, minopt );
	}
}


core::pose::Pose MakeJunctionsMover::get_and_cache_pdb(std::string pdb_location){
	if ( pdb_cache_bool_ ) {
		if ( pdb_cache_.find(pdb_location) != pdb_cache_.end() ) {
			return(pdb_cache_[pdb_location]);
		} else {
			core::pose::Pose pose;
			core::import_pose::pose_from_file(pose,*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ), pdb_location , core::import_pose::PDB_file);
			if ( pdb_cache_max_size_<pdb_cache_.size() ) {
				pdb_cache_.insert(std::pair<std::string,core::pose::Pose>(pdb_location,pose));
			}
			return(pose);
		}
	}
	core::pose::Pose pose;
	core::import_pose::pose_from_file(pose, pdb_location , core::import_pose::PDB_file);
	return(pose);
}
bool MakeJunctionsMover::check_all_junctions_good(MakeJunctionsMover::Design design){
	std::queue<std::string> tmp_queue;
	while ( !design.n_term_stack.empty() ) {
		std::string item = design.n_term_stack.top();
		if ( junction_failure_set_.find(item) != junction_failure_set_.end() ) {
			return(false);
		}
		tmp_queue.push(design.n_term_stack.top());
		design.n_term_stack.pop();
	}
	while ( !tmp_queue.empty() ) {
		design.n_term_stack.push(tmp_queue.front());
		tmp_queue.pop();
	}
	while ( !design.c_term_stack.empty() ) {
		std::string item = design.c_term_stack.top();
		if ( junction_failure_set_.find(item) != junction_failure_set_.end() ) {
			return(false);
		}
		tmp_queue.push(design.c_term_stack.top());
		design.c_term_stack.pop();
	}
	while ( !tmp_queue.empty() ) {
		design.c_term_stack.push(tmp_queue.front());
		tmp_queue.pop();
	}
	return(true);
}


bool MakeJunctionsMover::make_pose_from_design(MakeJunctionsMover::Design design,core::pose::PoseOP & return_pose){
	using namespace protocols::symmetry;
	using namespace basic::options;
	if ( !check_all_junctions_good(design) ) {
		return(false);
	}
	core::pose::Pose background_pose;
	core::pose::Pose pose;
	generate_start_pose(pose,background_pose,design.start_struct);
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	job->inner_job_nonconst()->optional_output_name(design.name);  //such a hack. Will only work with jd2.
	bool success = true;
	while ( (!design.n_term_stack.empty()) && success ) {
		std::string command = design.n_term_stack.top();
		success = attach_next_part(pose,"n_term",command);
		design.n_term_stack.pop();
	}
	while ( (!design.c_term_stack.empty())&& success ) {
		std::string command = design.c_term_stack.top();
		success = attach_next_part(pose,"c_term",command);
		design.c_term_stack.pop();
	}
	if ( background_pose.size()>0 ) {
		append_pose_to_pose(pose,background_pose);
	}
	renumber_pdbinfo_based_on_conf_chains(pose);
	return_pose = core::pose::PoseOP( new Pose( pose ) );
	return(success);
}


void
MakeJunctionsMover::apply( core::pose::Pose & pose )
{
	design_q_ = read_in_designs();
	core::pose::PoseOP tmpPoseOP=get_additional_output();
	if ( tmpPoseOP!=NULL ) {
		pose=*tmpPoseOP;
	}
}

core::pose::PoseOP MakeJunctionsMover::get_additional_output(){
	bool new_pose_generated = false;
	core::pose::PoseOP return_pose;
	while ( !new_pose_generated ) {
		if ( design_q_.empty() ) {
			set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
			return nullptr;
		}
		new_pose_generated = make_pose_from_design(design_q_.front(),return_pose);
		design_q_.pop();
	}
	return_pose->remove_constraints();
	sfxn_->score(*return_pose);
	return(return_pose);
}

void MakeJunctionsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ){
	designs_fn_ = tag->getOption< std::string >( "designs" ,"designs.txt" );
	pdb_cache_bool_ = tag->getOption< bool >( "pdb_cache" ,true);
	pdb_cache_max_size_ = tag->getOption< Size> ("pdb_size_size",2000);
	junction_rmsd_thresh_ = tag->getOption<Real>("junction_rmsd_thresh", 1.5);
	//relax_during_build_ = tag->getOption<bool>("relax_during_build",true);
	if ( tag->hasOption("scorefxn") ) {
		std::string const scorefxn_key( tag->getOption<std::string>("scorefxn") );
		if ( datamap.has( "scorefxns", scorefxn_key ) ) {
			sfxn_ = datamap.get_ptr< core::scoring::ScoreFunction >( "scorefxns", scorefxn_key );
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,"ScoreFunction " + scorefxn_key + " not found in basic::datacache::DataMap.");
		}
	}
}
void MakeJunctionsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "pdb_cache", xsct_rosetta_bool, "caches the junctions/DHR as required during runtime", "true" )
		+ XMLSchemaAttribute::attribute_w_default("pdb_cache_max_size",xsct_non_negative_integer,"max number of pdbs to save in cache.","2000")
		+ XMLSchemaAttribute::required_attribute( "designs", xs_string, "File with one design per line" )
		+ XMLSchemaAttribute::required_attribute( "scorefxn", xs_string, "Score function used for packing and design." )
		//+ XMLSchemaAttribute::attribute_w_default( "relax_during_build", xsct_rosetta_bool, "relaxes structure during build", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "junction_rmsd_thresh",xsct_real,"Not all junctions are perfectly identical structurally. in the first junction. This sets a threshold for allowed difference.","1.5");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Makes junctions from file", attlist );
}


} //simple_moves
} //protocols
