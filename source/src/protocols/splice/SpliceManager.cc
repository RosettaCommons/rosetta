// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/splice/SpliceManager.cc
/// @brief
/// @author Gideon Lapidoth (glapidoth@gmail.com)

// Unit headers
#include <protocols/splice/SpliceManager.hh>
#include <protocols/splice/Splice.hh>
#include <core/sequence/SequenceProfile.hh>
#include <boost/foreach.hpp>
#include <protocols/splice/util.hh>
#include <protocols/splice/SpliceSegment.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
#include <core/pose/util.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <sstream>
#include <utility/string_util.hh>
#include <dirent.h>
#include <protocols/toolbox/superimpose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <map>
#include <protocols/moves/Mover.hh>
#include <core/sequence/SequenceProfile.hh>
#include <utility/exit.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <protocols/splice/RBInMover.hh>
#include <protocols/splice/RBOutMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <boost/algorithm/string.hpp>
#include <core/kinematics/MoveMap.hh>



namespace protocols {
namespace splice {

static  basic::Tracer TR( "protocols.splice.SpliceManager" );
static  basic::Tracer TR_constraints( "protocols.splice.Splice_constraints" );

using namespace core::sequence;
using namespace std;

SpliceManager::SpliceManager():template_from_res_(0),template_to_res_(0),pose_from_res_(0),pose_to_res_(0),tail_seg_(""),residue_diff_(0),chain_num_(1), use_sequence_profiles_(false),segment_type_(""),add_sequence_constraints_only_ (false), profile_weight_away_from_interface_(1.0),
	rb_sensitive_(0),cut_site_(0),mm_(nullptr)
{
	dofs_.clear();
	splice_segments_.clear();
	pdb_segments_.clear();
}


SpliceManager::~SpliceManager() {}
///@brief set the fold tree around start/stop/cut sites.

void SpliceManager::fold_tree(core::pose::Pose & pose,utility::vector1<std::array<int, 3>> positions) const {
	core::kinematics::FoldTree ft;
	ft.clear();
	for ( auto i: positions ) {
		TR<<i[0]<<","<<i[1]<<","<<i[2]<<std::endl;
		ft.add_edge(i[0],i[1],i[2]);
	}
	TR<<ft<<std::endl;
	ft.reorder(anchor_res_);
	//ft.delete_extra_vertices();
	ft.delete_self_edges();
	ft.check_fold_tree();
	TR << "Previous ft: " << pose.fold_tree() << std::endl;
	pose.fold_tree(ft);
	TR << "Current ft: " << pose.fold_tree() << std::endl;
}


//set the phi/psi/omega angles of the segment in the pose. The angles are taken from either the torsion DB (in splice in) of the source PDB (in Splice out)
void
SpliceManager::set_BB_dofs(core::pose::Pose & pose){
	for ( core::Size i = 0; i < dofs().size(); ++i ) {
		core::Size const pose_resi(pose_from_res() + i);
		TR<<"Previous phi/psi/omega at resi: "<<pose_resi<<" "<<pose.phi( pose_resi )<<'/'<<pose.psi( pose_resi )<<'/'<<pose.omega( pose_resi )<<std::endl;
		pose.set_phi(pose_resi, dofs()[i + 1].phi());
		pose.set_psi(pose_resi, dofs()[i + 1].psi());
		pose.set_omega(pose_resi, dofs()[i + 1].omega());
		//pose.dump_pdb( "dump"+ utility::to_string( pose_resi ) + ".pdb" );
		TR<<"requested phi/psi/omega: "<<dofs()[ i + 1 ].phi()<<'/'<<dofs()[i+1].psi()<<'/'<<dofs()[i+1].omega()<<std::endl;
		TR<<"After change, phi/psi/omega: "<< pose_resi<<' '<<pose.phi( pose_resi )<<'/'<<pose.psi( pose_resi )<<'/'<<pose.omega( pose_resi )<<std::endl;
	}
}//set_BB_dofs

//for segments with no PSSM (e.g. H3 in antibodies) creates a pssm from the sequence using the BLOSUM62 profile and adds the profile to the Splice mover segment/pssm map
void SpliceManager::generate_profile_from_seq(std::map< std::string/*1AHW*/, std::string/*L1.1*/ >  segment_seq_map, std::string segmentName, std::string sourcename ,SpliceSegmentOP SpliceSegment){
	std::string segment_seq="";
	if ( segment_seq_map.size() ) {
		std::map<std::string,std::string>::iterator it;
		it=segment_seq_map.find(sourcename);
		if ( it !=segment_seq_map.end() ) {
			segment_seq=it->second;
		} else {
			utility_exit_with_message("Could not find pdb name: "+sourcename+" in "+segmentName+" sequence file\n");
		}
	}
	SequenceProfileOP segment_profile( new SequenceProfile );
	//TR<<"sequence of H3:"<<H3_seq<<std::endl;
	if ( segment_seq!="" ) {
		core::sequence::Sequence seq;
		seq.sequence(segment_seq);
		segment_profile->generate_from_sequence(seq);
		//TR<<"Size of H3 profile is: "<<H3_profile->size()<<std::endl;
	}
	SpliceSegment->add_profile(sourcename,segment_profile);
}

core::sequence::SequenceProfileOP
SpliceManager::generate_sequence_profile(core::pose::Pose & pose) {
	if ( !use_sequence_profiles_ ) {
		return NULL;
	}
	using namespace core::sequence;
	using namespace std;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	map<string, string> const comments = core::pose::get_all_comments(pose);
	if ( comments.size() < 1 ) { /// SJF changed from <3 22Jul14
		utility_exit_with_message(
			"Please check comments field in the pdb file (header= ##Begin comments##), could not find any comments");
	}

	if ( !add_sequence_constraints_only_ ) { ///If only doing sequence constraints then don't add to pose comments source name
		TR << "The current segment is: " << segment_type() << " and the source pdb is " << source_pdb_name() << std::endl;
		core::pose::add_comment(pose, "segment_" + segment_type(), source_pdb_name());//change correct association between current loop and pdb file
	}
	protocols::splice::load_pdb_segments_from_pose_comments(pose,pdb_segments_); // get segment name and pdb association from comments in pdb file
	TR << "There are " << pdb_segments().size() << " PSSM segments" << std::endl;

	runtime_assert(pdb_segments().size()); //This assert is in place to make sure that the pdb file has the correct comments, otherwise this function will fail

	utility::vector1<SequenceProfileOP> profile_vector;

	profile_vector.clear(); //this vector holds all the pdb segment profiless

	for ( std::string const segment_type: segment_names_ordered_ ) { //<- Start of PDB segment iterator
		TR<<"segment_type: "<<segment_type<<std::endl;
		TR<<"Map size: "<<splice_segments_[ segment_type ]->pdb_to_profile_map().size()<<std::endl;
		if ( splice_segments_[ segment_type ]->pdb_profile(pdb_segments_[segment_type])==0 ) {
			utility_exit_with_message(" could not find the source pdb name: "+ pdb_segments_[segment_type]+ ", in pdb_profile_match file."+segment_type+" or PSSM file is missing\n");
		}
		TR<<"reading profile:"<< pdb_segments_[segment_type]<<std::endl;
		profile_vector.push_back( splice_segments_[ segment_type ]->pdb_profile( pdb_segments_[segment_type] ));
	} // <- End of PDB segment iterator
	TR << "The size of the profile vector is: " << profile_vector.size() << std::endl;

	///Before upwweighting constraint we check that the PSSMs are properly aligned by making sure
	///that the PSSM score of the falnking segments of the current designed segment agree with the identity of the aa (i.e if
	/// we are designing L1 then we would expect that segments Frm.light1 and Frm.light2 have concensus aa identities)


	std::map<std::string,std::string> pose_comments=get_all_comments(pose);
	std::string H3_pdb_source=pose_comments["segment_H3"];
	//TR<<"H3_pdb_source: "<<H3_pdb_source<<std::endl;
	std::string H3_seq="";
	//TR<<"Size of pdb_to_H3_seq_map_: "<<pdb_to_H3_seq_map_.size()<<std::endl;
	//TR<<"H3 sequence: "<<H3_seq<<std::endl;
	return concatenate_profiles(profile_vector, segment_names_ordered_,H3_seq);
}//generate sequence profiles

void
SpliceManager::add_sequence_constraints(core::pose::Pose & pose) {
	if ( !add_sequence_constraints_only_ ) { ///If only doing sequence constraints then don't add to pose comments source name
		TR << "The current segment is: " << segment_type() << " and the source pdb is " << source_pdb_name() << std::endl;
		core::pose::add_comment(pose, "segment_" + segment_type(), source_pdb_name()); //change correct association between current loop and pdb file
	}

	if ( !use_sequence_profiles() ) {
		TR << TR.Red << "WARNING! Not using sequence constraints"""<<TR.Reset << std::endl;
		return;
	}
	using namespace core::scoring::constraints;

	/// first remove existing sequence constraints
	TR << "Removing existing sequence profile constraints from pose" << std::endl;
	ConstraintCOPs constraints(pose.constraint_set()->get_all_constraints());
	TR << "Total number of constraints at start: " << constraints.size() << std::endl;
	core::Size cst_num(0);
	for ( ConstraintCOP const c: constraints ) {
		if ( ( c->type() == "SequenceProfile" ) and (c->residues()[1]>=pose.conformation().chain_begin(chain_num())) and (c->residues()[1]<=pose.conformation().chain_end(chain_num())) ) { //only remove profile sequence constraints of pose that is being splice out
			pose.remove_constraint( c );
			cst_num++;
		}
	}
	TR << "Removed a total of " << cst_num << " sequence constraints." << std::endl;
	TR << "After removal the total number of constraints is: " << pose.constraint_set()->get_all_constraints().size()<< std::endl;
	/// then impose new sequence constraints
	core::sequence::SequenceProfileOP seqprof(generate_sequence_profile(pose));
	TR << "chain_num: " << chain_num() << std::endl;
	TR << "Chain length/seqprof size: " << pose.conformation().chain_end(chain_num()) - pose.conformation().chain_begin(chain_num()) + 1 << ", "<< seqprof->size() - 1 << std::endl;
	/*   std::string pdb_dump_fname_("after_splice.pdb");
	std::ofstream out( pdb_dump_fname_.c_str() );
	pose.dump_pdb(out); //Testi*/

	runtime_assert( seqprof->size() - 1 == pose.conformation().chain_end(chain_num()) - pose.conformation().chain_begin(chain_num()) + 1); //Note that the minus 1 after seqprof size is because seqprof size is always +1 to the actual size. Do not chnage this!!
	cst_num = 0;
	TR << "Up-weighting sequence constraints " << std::endl;

	//If pose has more than one chain the sequence profile mapping needs to be modified accordingly
	core::id::SequenceMappingOP smap( new core::id::SequenceMapping() );
	for ( core::Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ) {
		if ( (seqpos >= pose.conformation().chain_begin(chain_num()))
				and (seqpos <= pose.conformation().chain_end(chain_num())) ) {
			smap->push_back(seqpos - pose.conformation().chain_begin(chain_num()) + 1);
			//if position is in chain >1 then it should have a sequence profile, the mapping should look like this: 0|0|0|1|2|3|
			// So the first positions ('0') does not have a seqeuence constraint
		} else {
			smap->push_back(0);
		}
	}
	check_sequence_profile(pose,smap,seqprof );
	//smap->show();//uncomment to see mapping in tracer

	utility::vector1<core::Size> const non_upweighted_residues(protocols::splice::find_residues_on_chain1_inside_interface(pose, chain_num()));
	for ( core::Size seqpos = pose.conformation().chain_begin(chain_num()); seqpos <= pose.conformation().chain_end(chain_num()); ++seqpos ) {
		using namespace core::scoring::constraints;

		SequenceProfileConstraintOP spc(new SequenceProfileConstraint(pose, seqpos, seqprof, smap));
		//TR << "sepos=" << spc->seqpos() << std::endl;
		if ( std::find(non_upweighted_residues.begin(), non_upweighted_residues.end(), seqpos)
				== non_upweighted_residues.end() ) { //seqpos not in interface so upweight
			spc->weight(profile_weight_away_from_interface());
		}
		//TR << "Now adding constraint to aa: " << seqpos << pose.aa(seqpos)<< std::endl;
		//TR<<"The sequence profile row for this residue is: "<<seqprof->prof_row(seqpos-)<<std::endl;
		pose.add_constraint(spc);
		cst_num++;
		//TR << "Current constraints size: "<< pose.constraint_set()->get_all_constraints().size() << std::endl;
	}
	TR << "Added a total of " << cst_num << " sequence constraints." << std::endl;
	TR << "Now the pose has a total of " << pose.constraint_set()->get_all_constraints().size() << " constraints"
		<< std::endl;

	/// just checking that the scorefxn has upweighted res_type_constraint
	core::Real const score_weight(scorefxn()->get_weight(core::scoring::res_type_constraint));
	TR << "res_type_constraint weight is set to " << score_weight << std::endl;
	if ( score_weight <= 0.001 ) {
		TR
			<< "Warning! res_type_constraint weight is low, even though I've just added sequence constraints to the pose! These sequence constraints will have no effect. This could be an ERROR"
			<< std::endl;
	}
}

void
SpliceManager::rb_adjust_template(core::pose::Pose const & pose) const {
	if ( !rb_sensitive() ) {
		return;
	}
	TR << "CNDEBUG Now applying the Rigid body orientation of pose to template pose " << std::endl;

	//template_pose_->dump_pdb(mover_name_ + "_template_should_be_pose_aligned_BEFORE_ALIGN.pdb" );
	//pose.dump_pdb( mover_name_ + "_pose_should_be_template_aligned_BEFORE_ALIGN.pdb" );

	RBOutMover rbo;
	//template_pose_->dump_pdb( "pre_jump_template.pdb" );
	core::pose::Pose copy_pose(pose); /// I don't want the fold tree to change on pose...
	core::pose::Pose chainA = *pose.split_by_chain(1);
	core::kinematics::Jump const pose_jump = rbo.get_disulf_jump(copy_pose, chainA);

	RBInMover rbi;
	rbi.set_fold_tree(*template_pose_);
	TR << "The RBO will be applied on template with foldtree " << (*template_pose_).fold_tree() << std::endl;
	template_pose_->set_jump(1, pose_jump);
	/*if (debug_){
	template_pose_->dump_pdb(mover_name_ + "_template_should_be_pose_aligned.pdb");
	pose.dump_pdb( mover_name_ + "_pose_should_be_template_aligned.pdb");
	// runtime_assert( 0 );
	}*/
}

//The from_res and to_res are given by the template pose. The input pose could be different than the template, so we adjust the from_res and to_res to the numbering
//according to the input pose
void
SpliceManager::adjust_stem_positions_by_template(core::pose::Pose const & pose) {
	if ( !(adjust_stem_positions_by_template_) ) {
		return;
	}
	if ( template_pose() == nullptr ) {
		return;
	}
	rb_adjust_template(pose);
	if ( template_from_res() && template_to_res() && tail_seg()=="" ) {
		pose_from_res(protocols::rosetta_scripts::find_nearest_res(pose, *template_pose_, template_from_res(), 0/*chain*/));
		pose_to_res (protocols::rosetta_scripts::find_nearest_res(pose, *template_pose_, template_to_res(), 0/*chain*/));
		TR<<"pose_from_res():"<<pose_from_res()<<std::endl;
		TR<<"pose_to_res():"<<pose_to_res()<<std::endl;
		runtime_assert(template_from_res()>0);
		runtime_assert(template_to_res()>0);
		return;
	} else if ( template_from_res() && tail_seg()!="" ) {
		core::conformation::Conformation const & conf(pose.conformation());
		if ( boost::iequals(tail_seg(), "c") ) {
			pose_to_res(conf.chain_end(1));
			pose_from_res(protocols::rosetta_scripts::find_nearest_res(pose, *template_pose_, template_from_res(), 0/*chain*/));

		} else if ( boost::iequals(tail_seg(), "n") ) {
			pose_from_res(1);
			pose_to_res(protocols::rosetta_scripts::find_nearest_res(pose, *template_pose_, template_to_res(), 0/*chain*/));
		}
		TR<<"pose_from_res():"<<pose_from_res()<<std::endl;
		TR<<"pose_to_res():"<<pose_to_res()<<std::endl;
		runtime_assert(template_from_res()>0);
		runtime_assert(template_to_res()>0);
		return;
	}
}

void
SpliceManager::update_pose_stem_positions() {
	pose_to_res(pose_to_res() +residue_diff());
}

void SpliceManager::check_sequence_profile(core::pose::Pose & pose, core::id::SequenceMappingOP smap, core::sequence::SequenceProfileOP seqprof){
	for ( core::Size row = 1; row < seqprof->size(); row++ ) { //go over all the PSSM sements provided by the user
		utility::vector1< core::Size > cur_prof_row = seqprof->prof_row(row);


		auto it = std::find_if(std::begin(cur_prof_row), std::end(cur_prof_row), [](int pos){return pos >= 8;});
		core::Size prof_AA = std::distance(std::begin(cur_prof_row), it)+1;

		if ( prof_AA<=20 ) {
			core::Size pose_pos = smap->get_corresponding_residue_in_current(row);
			if ( string(1,pose.residue(pose_pos).name1())!=seqprof->alphabet()[prof_AA] ) {
				std::string seqpos, pssm_AA,pose_AA;
				std::ostringstream convert;
				convert << pose_pos; // insert the textual representation of 'Number' in the characters in the stream
				seqpos = convert.str();
				convert.str("");
				convert.clear();
				convert<<seqprof->alphabet()[prof_AA];
				pose_AA = convert.str();
				TR<<TR.Red<<" PSSM and pose might be misaligned, position " +  seqpos + " should be a "+pose_AA+ "\n"<<TR.Reset<<std::endl;
			}

		}
	} //end pssm segment for
}

void SpliceManager::parse_tags(TagCOP const tag,basic::datacache::DataMap &data){
	design_shell(tag->getOption<core::Real>("design_shell", 6.0)); //Added by gideonla May13,
	repack_shell(tag->getOption<core::Real>("repack_shell", 8.0)); //Added by gideonla May13,
	mover_name_=tag->getOption<std::string>("name");//for debugging purposes
	debug_=tag->getOption<bool>("debug", false);
	rb_sensitive(tag->getOption<bool>("rb_sensitive", false));
	chain_num(tag->getOption<core::Size>("chain_num", 1));


	template_file(tag->getOption<std::string>("template_file", ""));
	if ( template_file_ != "" ) { /// using a template file to determine from_res() to_res()
		if ( data.has("poses", template_file_) ) {
			template_pose_ = data.get_ptr<core::pose::Pose>("poses", template_file_);
			TR << "using template pdb from datamap" << std::endl;
		} else  {
			template_pose_ = core::pose::PoseOP( new core::pose::Pose );
			template_pose_ = core::import_pose::pose_from_file(template_file_);
			data.add("poses", template_file_, template_pose_);
			TR << "loading template_pose from " << template_file_ << std::endl;
		}
	}
	if ( tag->hasOption("task_operations") ) {
		task_factory(protocols::rosetta_scripts::parse_task_operations(tag->getOption<std::string>("task_operations"),data));
	}
	dbase_file_name(tag->getOption<std::string>("torsion_database", ""));
	thread_original_sequence(tag->getOption<bool>("thread_original_sequence", false)); // Will not do sequence design on the inserted segments and save the sequence from the source PDB/Torsion db

}

void SpliceManager::modify_pdb_segments_with_current_segment(std::string const & pdb_name) {
	pdb_segments()[segment_type()] = pdb_name;
}

void
SpliceManager::chainbreak_check( core::pose::Pose const & pose , core::Real const tolerance , bool fail_retry_if_found , bool crash_if_found){

	core::Size cuts = 0;
	core::Size chain_id = chain_num();
	core::Real max_bl = 1.33+tolerance; // The tolerance should be no less than 0.13
	core::Real min_bl = 1.33-tolerance;
	if ( chain_id == 0 ) {
		chain_id = 1;
	}
	TR << "Will check peptide bond lengths between " << pose.conformation().chain_begin( chain_id ) << " to " << pose.conformation().chain_end( chain_id )-1 << std::endl;
	TR<< "bond length tolerance value is:"<<tolerance<<std::endl;
	for ( core::Size resj = pose.conformation().chain_begin( chain_id ); resj <= pose.conformation().chain_end( chain_id )-1; ++resj ) {
		core::Real const distance = pose.residue( resj + 1 ).xyz( "N" ).distance(pose.residue( resj ).xyz( "C" ));
		if ( distance > max_bl || distance < min_bl ) {
			TR<<"distance is: "<<distance<<" max_bl"<<max_bl<<std::endl;
			TR << "The distance from " << resj << " to " << resj+1 << " is " << distance << std::endl;
			++cuts;
		}
		if ( cuts > allowed_cuts_ ) {
			TR << TR.Red << "WARNING! There are more than "<<allowed_cuts_<<" cutpoints"<<TR.Reset << std::endl;
			if ( crash_if_found ) {
				TR <<TR.Red<< "Dumping pdb for inspection: " << TR.Underline<< mover_name()+ "_more_than_one_cutpoint.pdb" << std::endl;
				pose.dump_pdb(mover_name()+ "_more_than_one_cutpoint.pdb" );
				utility_exit_with_message("ERROR: There are more than one cutpoint \n");
			} else if ( fail_retry_if_found ) {
				//protocols::moves::Mover::set_last_move_status(protocols::moves::FAIL_RETRY);
				return;
			}
		}
	}
}

void SpliceManager::add_coordinate_constraints(core::pose::Pose & pose, core::pose::Pose const & source_pose, core::Size from, core::Size to, core::Size anchor, std::string atom_type, core::pack::task::PackerTaskOP task) {
	// if( !superimposed() ){
	//  TR<<"The source pose is not superimposed on the template, so we're not using coordinate restraints, only dihedrals"<<std::endl;
	//  return;
	// }
	core::scoring::constraints::ConstraintOPs cst;
	TR_constraints << "pose anchor residue : " <<anchor << std::endl;
	core::Size anchor_source = protocols::rosetta_scripts::find_nearest_res(source_pose, pose, anchor, 0);
	runtime_assert( anchor_source );
	TR_constraints << "closest residue to anchor residue on source is : " <<anchor_source << std::endl;
	int res_diff = anchor_source - anchor;
	TR_constraints << "res diff in coordinate constarint is: " << res_diff << std::endl;
	core::Size const fixed_res(anchor);
	TR_constraints << "Anchor residue for the coordinate constraint is " << fixed_res << std::endl;
	TR_constraints << "Current pose "<<atom_type<< " xyz coordinate/source pdb "<<atom_type<< " xyz coordinate:" << std::endl;
	core::id::AtomID const anchor_atom(core::id::AtomID(pose.residue(fixed_res).atom_index("CA"), fixed_res));
	for ( core::Size i = from; i <= to; ++i ) {
		if ( atom_type=="CB" && ((pose.residue(i).name3()=="GLY")||(source_pose.residue(i + res_diff).name3()=="GLY")) ) { //Gly doesn't have CB so we should just skip
			TR_constraints<<"Found GLY! Not applying CB constraints here"<<std::endl;
			continue;
		}

		if ( atom_type=="CG" ) {
			if ( pose.fold_tree().is_cutpoint(i) ) {
				TR<<i<<"is cutpoint skipping"<<std::endl;
				continue;//don't want to add CG constraints to cut point residues. This causes problems to CCD.
			}
			TR_constraints<<"Allowed aa's for residue "<<i<<" are: ";
			std::list< core::chemical::ResidueTypeCOP > allowed_aas =task->residue_task(i).allowed_residue_types();
			std::vector<char> allowed_aas_names;
			allowed_aas_names.clear();
			for ( std::list< core::chemical::ResidueTypeCOP >::const_iterator restype = allowed_aas.begin(); restype != allowed_aas.end(); ++restype ) {
				TR_constraints<<(*restype )->name1()<<",";
				allowed_aas_names.push_back((*restype )->name1());
			}
			TR_constraints<<std::endl;
			std::vector<char> aromatic_and_his;
			aromatic_and_his.clear();
			aromatic_and_his.push_back('Y');aromatic_and_his.push_back('F');aromatic_and_his.push_back('W');aromatic_and_his.push_back('H');
			std::vector<char> intersect;
			std::sort(allowed_aas_names.begin(), allowed_aas_names.end());
			std::vector<char>::iterator it;
			it = std::unique(allowed_aas_names.begin(), allowed_aas_names.end());
			allowed_aas_names.resize( std::distance(allowed_aas_names.begin(),it) );
			std::sort(aromatic_and_his.begin(), aromatic_and_his.end());
			std::set_intersection(allowed_aas_names.begin(),allowed_aas_names.end(),aromatic_and_his.begin(),aromatic_and_his.end(),back_inserter(intersect));
			TR<<"size of intersect is"<<intersect.size()<<std::endl;
			if ( (allowed_aas_names.size()==1) && (pose.residue(i).name3()!="GLY")&&(pose.residue(i).name3()==source_pose.residue(i+res_diff).name3()) ) { //Apply when we have strict conservation
				//change chi to match the source pdb
				for ( core::Size chi=1; chi<source_pose.residue(i + res_diff).nchi(); chi++ ) {
					core::Real const chi_curr=source_pose.chi(chi,i + res_diff);
					pose.set_chi(chi,i,chi_curr);
				}
				for ( core::Size atmnum=source_pose.residue(i + res_diff).first_sidechain_atom(); atmnum<=source_pose.residue(i + res_diff).natoms(); atmnum++ ) {
					if ( source_pose.residue(i + res_diff).atom_is_hydrogen(atmnum) ) { //Don't want to add constraints to H atoms
						continue;
					}
					if ( pose.residue(i).atom_name(atmnum).compare("CB") == 0 ) { //Already added CB constratins
						TR<<"FOUND CB ATOM IN :"<<pose.residue(i).atom_name(atmnum)<<std::endl;
						continue;
					}
					core::scoring::func::FuncOP coor_cont_fun( new core::scoring::func::HarmonicFunc(0.0, 1) );
					TR_constraints<<"Applying constraints to atom:"<<pose.residue(i).atom_name(atmnum)<<",of residue "<<pose.residue(i).name3()<<i;
					TR_constraints<<"Taking xyz coordinates from atom:"<<source_pose.residue(i + res_diff).xyz(atmnum)[0]<<"(x),"<<source_pose.residue(i + res_diff).xyz(atmnum)[1]<<"(y),"<<source_pose.residue(i + res_diff).xyz(atmnum)[2]<<"(z), pose atom xyz: "<<pose.residue(i).xyz(atmnum)[0]<<std::endl;
					cst.push_back(core::scoring::constraints::ConstraintOP( new core::scoring::constraints::CoordinateConstraint(core::id::AtomID(atmnum, i),anchor_atom, source_pose.residue(i + res_diff).xyz(atmnum), coor_cont_fun) ));
					//pose.add_constraints(cst);
				}//for atom_it
				//using namespace std;
				//string String = static_cast<ostringstream*>( &(ostringstream() << i) )->str();

			} else if ( (allowed_aas_names.size()>0)&&(intersect.size()==allowed_aas_names.size())&&(std::find(aromatic_and_his.begin(), aromatic_and_his.end(),source_pose.residue(i + res_diff).name1())!=aromatic_and_his.end()) ) { //fi allowed_aas_name//if same size then all allowed identities are in the aromatic vector
				for ( core::Size chi=1; chi<source_pose.residue(i + res_diff).nchi(); chi++ ) {
					core::Real const chi_curr=source_pose.chi(chi,i + res_diff);
					pose.set_chi(chi,i,chi_curr);
					utility::vector1< core::Size > const chiAtoms=pose.residue(i).chi_atoms(chi);//get chi atoms from chi angle so we can apply coordiante constraints
					for ( core::Size chiAtom=1; chiAtom<=4; chiAtom++ ) {
						core::scoring::func::FuncOP coor_cont_fun( new core::scoring::func::HarmonicFunc(0.0, 1) );
						cst.push_back(core::scoring::constraints::ConstraintOP( new core::scoring::constraints::CoordinateConstraint(core::id::AtomID(chiAtoms[chiAtom], i),anchor_atom, source_pose.residue(i + res_diff).xyz(chiAtoms[chiAtom]), coor_cont_fun) ));
						TR<<"Applying constraints to chi_atom:"<<chiAtoms[chiAtom]<<pose.residue(i).atom_name(chiAtoms[chiAtom])<<",of residue "<<pose.residue(i).name3()<<i<<std::endl;
						//pose.add_constraints(cst);
					}//for chiAtom
				}//for chi
				// using namespace std;
				// string String = static_cast<ostringstream*>( &(ostringstream() << i) )->str();
			} //else fi
			continue;
		}//fi CG_atom
		core::scoring::func::FuncOP coor_cont_fun( new core::scoring::func::HarmonicFunc(0.0, 1) );
		cst.push_back(core::scoring::constraints::ConstraintOP( new core::scoring::constraints::CoordinateConstraint(core::id::AtomID(pose.residue(i).atom_index(atom_type), i),
			anchor_atom, source_pose.residue(i + res_diff).atom(atom_type).xyz(), coor_cont_fun) ));
		//Print xyz coor of current pose CA atoms vs. source pose
		TR_constraints << i << pose.aa(i) << " " << pose.residue(i).atom(atom_type).xyz()[0] << "," << pose.residue(i).atom(atom_type).xyz()[1]
			<< "," << pose.residue(i).atom(atom_type).xyz()[2] << " / " << i + res_diff << source_pose.aa(i + res_diff) << " "
			<< source_pose.residue(i + res_diff).atom(atom_type).xyz()[0] << "," << source_pose.residue(i + res_diff).atom(atom_type).xyz()[1] << "," << source_pose.residue(i + res_diff).atom(atom_type).xyz()[2]
			<< std::endl;
		//pose.add_constraints(cst);
	}//for
	pose.add_constraints(cst);
	protocols::splice::report_coordinate_constraints(pose);
}//func

void SpliceManager::add_dihedral_constraints(core::pose::Pose & pose, core::pose::Pose const & source_pose, core::Size start_res_on_pose, core::Size end_res_on_pose) {
	int residue_diff;
	core::Size nearest_res_on_source_pdb(0);
	TR<<"anchor_res:"<<anchor_res_<<std::endl;
	nearest_res_on_source_pdb=(protocols::rosetta_scripts::find_nearest_res(source_pose,pose,anchor_res_,0));

	TR<<"nearest_res_on_source_pdb:"<<nearest_res_on_source_pdb<<std::endl;
	/* if (boost::iequals(tail_seg(), "n"))
	nearest_res_on_source_pdb=(protocols::rosetta_scripts::find_nearest_res(source_pose,pose,end_res_on_pose,0));//in tail segment the first res in the n-ter tail will probably won't be aligned
	else
	nearest_res_on_source_pdb=(protocols::rosetta_scripts::find_nearest_res(source_pose,pose,start_res_on_pose,0));
	runtime_assert(nearest_res_on_source_pdb); //if find_nearest_res fails it return 0
	//To compare the angles we keep track of the corresponding residues in the template pose
	runtime_assert(nearest_res_on_source_pdb); //if find_nearest_res fails it return 0*/
	residue_diff =nearest_res_on_source_pdb- anchor_res_ ;
	TR_constraints<<"Residue diff is:"<<residue_diff<<std::endl;
	TR_constraints << "Applying dihedral constraints to pose:" << std::endl;
	core::scoring::constraints::ConstraintOPs csts; //will hold dihedral constraints
	//Set up constraints for the phi angle
	for ( core::Size i = start_res_on_pose + residue_diff; i < end_res_on_pose + residue_diff; ++i ) {

		core::id::AtomID phi_resi_n(source_pose.residue_type(i).atom_index("N"), i);
		core::id::AtomID pose_phi_resi_n(pose.residue_type(i-residue_diff).atom_index("N"), i - residue_diff);
		numeric::xyzVector<core::Real> xyz_Ni = source_pose.residue(i).atom("N").xyz();

		core::id::AtomID phi_resj_c(source_pose.residue_type(i - 1).atom_index("C"), i - 1);
		core::id::AtomID pose_phi_resj_c(pose.residue_type(i - 1-residue_diff).atom_index("C"), i - 1 -residue_diff);
		numeric::xyzVector<core::Real> xyz_Cj = source_pose.residue(i - 1).atom("C").xyz();

		core::id::AtomID phi_resi_co(source_pose.residue_type(i).atom_index("C"), i);
		core::id::AtomID pose_phi_resi_co(pose.residue_type(i-residue_diff).atom_index("C"), i-residue_diff);
		numeric::xyzVector<core::Real> xyz_Ci = source_pose.residue(i).atom("C").xyz();

		core::id::AtomID phi_resi_ca(source_pose.residue_type(i).atom_index("CA"), i);
		core::id::AtomID pose_phi_resi_ca(pose.residue_type(i-residue_diff).atom_index("CA"), i-residue_diff);
		numeric::xyzVector<core::Real> xyz_Cai = source_pose.residue(i).atom("CA").xyz();

		TR_constraints << "Phi: " << i - residue_diff << pose.aa(i - residue_diff) << ":" << pose.phi(i - residue_diff)
			<< " / " << i << source_pose.aa(i) << ":" << numeric::dihedral_degrees(xyz_Cj, xyz_Ni, xyz_Cai, xyz_Ci)
			<< std::endl;
		core::scoring::func::FuncOP di_const_func_phi( new core::scoring::func::CircularHarmonicFunc(
			(source_pose.phi(i) * numeric::constants::d::pi_2) / 360, 1) );
		csts.push_back(
			core::scoring::constraints::ConstraintOP( new core::scoring::constraints::DihedralConstraint(pose_phi_resj_c, pose_phi_resi_n, pose_phi_resi_ca, pose_phi_resi_co,
			di_const_func_phi) ));

		//Set up constraints for the psi angle
		core::id::AtomID psi_resi_n(source_pose.residue_type(i).atom_index("N"), i);
		core::id::AtomID pose_psi_resi_n(pose.residue_type(i-residue_diff).atom_index("N"), i-residue_diff);
		xyz_Ni = source_pose.residue(i).atom("N").xyz();

		core::id::AtomID psi_resj_n(source_pose.residue_type(i + 1).atom_index("N"), i + 1);
		core::id::AtomID pose_psi_resj_n(pose.residue_type(i + 1-residue_diff).atom_index("N"), i + 1-residue_diff);
		numeric::xyzVector<core::Real> xyz_Nj = source_pose.residue(i + 1).atom("N").xyz();

		core::id::AtomID psi_resi_co(source_pose.residue_type(i).atom_index("C"), i);
		core::id::AtomID pose_psi_resi_co(pose.residue_type(i-residue_diff).atom_index("C"), i-residue_diff);
		xyz_Ci = source_pose.residue(i).atom("C").xyz();

		core::id::AtomID psi_resi_ca(source_pose.residue_type(i).atom_index("CA"), i);
		core::id::AtomID pose_psi_resi_ca(pose.residue_type(i-residue_diff).atom_index("CA"), i-residue_diff);
		xyz_Cai = source_pose.residue(i).atom("CA").xyz();

		//for each residue the ideal angle is taken from the source pdb
		core::scoring::func::FuncOP di_const_func_psi( new core::scoring::func::CircularHarmonicFunc(
			(source_pose.psi(i) * numeric::constants::d::pi_2) / 360, 1) );
		csts.push_back(
			core::scoring::constraints::ConstraintOP( new core::scoring::constraints::DihedralConstraint(pose_psi_resi_n, pose_psi_resi_ca, pose_psi_resi_co, pose_psi_resj_n,
			di_const_func_psi) ));
		TR_constraints << "Psi: " << i - residue_diff << pose.aa(i - residue_diff) << ":" << pose.psi(i - residue_diff)
			<< " / " << i << source_pose.aa(i) << ":" << numeric::dihedral_degrees(xyz_Ni, xyz_Cai, xyz_Ci, xyz_Nj)
			<< std::endl;
		//Set up constraints for the omega angle
		core::id::AtomID omega_resj_n(source_pose.residue_type(i + 1).atom_index("N"), i + 1);
		core::id::AtomID pose_omega_resj_n(pose.residue_type(i + 1-residue_diff).atom_index("N"), i + 1-residue_diff);
		xyz_Ni = source_pose.residue(i + 1).atom("N").xyz();

		core::id::AtomID omega_resi_ca(source_pose.residue_type(i).atom_index("CA"), i);
		core::id::AtomID pose_omega_resi_ca(pose.residue_type(i-residue_diff).atom_index("CA"), i-residue_diff);
		xyz_Cai = source_pose.residue(i).atom("CA").xyz();

		core::id::AtomID omega_resi_co(source_pose.residue_type(i).atom_index("C"), i);
		core::id::AtomID pose_omega_resi_co(pose.residue_type(i-residue_diff).atom_index("C"), i-residue_diff);
		xyz_Ci = source_pose.residue(i).atom("C").xyz();

		core::id::AtomID omega_resj_ca(source_pose.residue_type(i + 1).atom_index("CA"), i + 1);
		core::id::AtomID pose_omega_resj_ca(pose.residue_type(i + 1-residue_diff).atom_index("CA"), i + 1-residue_diff);
		numeric::xyzVector<core::Real> xyz_Caj = source_pose.residue(i + 1).atom("CA").xyz();
		TR_constraints << "omega: " << i - residue_diff << pose.aa(i - residue_diff) << ":" << pose.omega(i - residue_diff)
			<< " / " << i << source_pose.aa(i) << ":" << numeric::dihedral_degrees(xyz_Cai, xyz_Ci, xyz_Nj, xyz_Caj)
			<< std::endl;
		//for each residue the ideal angle is taken from the "donor" pdb
		core::scoring::func::FuncOP di_const_func_omega( new core::scoring::func::CircularHarmonicFunc(
			(source_pose.omega(i) * numeric::constants::d::pi_2) / 360, 1) );
		csts.push_back(
			core::scoring::constraints::ConstraintOP( new core::scoring::constraints::DihedralConstraint(pose_omega_resi_ca, pose_omega_resi_co, pose_omega_resj_n, pose_omega_resj_ca,
			di_const_func_omega) ));
	}
	core::Real const score_weight(scorefxn()->get_weight(core::scoring::dihedral_constraint));
	TR_constraints << "dihedral_constraint weight is set to " << score_weight << std::endl;
	//scorefxn()->show(pose);
	pose.add_constraints(csts);
	//pose.dump_pdb("at_end_of_dihedral_const.pdb");
}

void SpliceManager::parse_segments(utility::vector1<TagCOP> const & sub_tags,TagCOP const tag, basic::datacache::DataMap &data){
	typedef utility::vector1<std::string> StringVec;
	bool check_segment = false;
	if ( sub_tags.empty() ) {
		return;
	}
	for ( TagCOP const sub_tag: sub_tags ) {
		TR<<sub_tag->getName()<<std::endl;
		if ( sub_tag->getName() == "Segments" ) {

			segment_names_ordered().clear(); //This string vector holds all the segment names inserted by the user, the pssm will be constructed in this order so user must make sure the order is correct
			check_segment = false;//This is set to false unless current segment appears in the segment list
			use_sequence_profiles(true);
			profile_weight_away_from_interface( tag->getOption< core::Real >( "profile_weight_away_from_interface", 1.0 ) );
			segment_type(sub_tag->getOption< std::string >( "current_segment" ));

			TR<<"reading segments in "<<tag->getName()<<std::endl;

			utility::vector1< TagCOP > const segment_tags( sub_tag->getTags() );
			for ( TagCOP const segment_tag: segment_tags ) {
				SpliceSegmentOP splice_segment( new SpliceSegment );
				std::string const segment_name( segment_tag->getOption< std::string >( "name" )); //get name of segment from xml
				std::string const pdb_profile_match( segment_tag->getOption< std::string >( "pdb_profile_match" ) );// get name of pdb profile match, this file contains all the matching between pdb name and sub segment name, i.e L1.1,L1.2 etc
				std::string const profiles_str( segment_tag->getOption< std::string >( "profiles" ) );
				StringVec const profile_name_pairs( utility::string_split( profiles_str, ',' ) );
				if ( data.has("segments", segment_name) ) {
					splice_segments().insert( std::pair< std::string, SpliceSegmentOP >( segment_name, data.get_ptr<SpliceSegment>("segments", segment_name)  ));
					TR << "using segment "<<segment_name<<" from datamap" << std::endl;
					segment_names_ordered().push_back(segment_name);
					if ( segment_name.compare(segment_type()) == 0 ) {
						check_segment=true;
					}
					continue;
				}
				// TR<<"Now working on segment:"<<segment_name<<std::endl;
				for ( std::string const s: profile_name_pairs ) {
					std::size_t found=s.find(':');
					if ( found==std::string::npos ) {
						continue; //if no : seperator found in string this is an error and we skip this entry
					}
					StringVec const profile_name_file_name( utility::string_split( s, ':' ) );
					//TR<<"pssm file:"<<profile_name_file_name[ 2 ]<<",segment name:"<<profile_name_file_name[ 1 ]<<std::endl;

					splice_segment->read_profile( profile_name_file_name[ 2 ], profile_name_file_name[ 1 ] );
				}

				splice_segment->read_pdb_profile( pdb_profile_match );
				//TR<<"the segment name is: "<<segment_name<<std::endl;
				if ( segment_name.compare(segment_type()) == 0 ) { //make sure that the current segment being spliced exists in the pssm segments
					check_segment=true;
				}
				splice_segments().insert( std::pair< std::string, SpliceSegmentOP >( segment_name, splice_segment ) );
				segment_names_ordered().push_back(segment_name);
				data.add("segments", segment_name, splice_segment);
			} //foreach segment_tag
		} // fi Segments
		if ( tag->hasOption("use_sequence_profile") ) {
			add_sequence_constraints_only(tag->getOption<bool>("use_sequence_profile"));
		}
		if ( !check_segment && !segment_names_ordered().empty() ) { //sanity check to make sure the current segment (what is being spliced) is also in the list of segments
			utility_exit_with_message("Segment " + segment_type() + " was not found in the list of segemnts. Check XML file\n");
		}
		if ( segment_names_ordered().empty() ) { //If splicing segment but not using sequence profile then turn off "use_sequence_profiles"
			use_sequence_profiles(false);
		}
	}
} //fi (sub_tags!=NULL)


} //splice
} //protocols
