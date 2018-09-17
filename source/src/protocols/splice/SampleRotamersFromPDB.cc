/*
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/splice/SampleRotamersFromPDB.cc
/// @brief  Given source PDBs limits the rotamer sampling to those found in equivalent positions in the PDB
/// @detail At initialization the taskoperation reads a database with allowed rotamers for each position
/// and changes the Rotmer set accordingly
/// @author Gideon Lapidoth ( glapidoth@gmail.com )
*/

/// @brief samples rotamers from a pre-computed db
/// @author Gideon Lapidoth (glapidoth@gmail.com)

// Unit Headers
#include <core/pose/extra_pose_info_util.hh>
#include <protocols/splice/SampleRotamersFromPDB.hh>
#include <protocols/splice/SampleRotamersFromPDBCreator.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>
#include <sstream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <map>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/splice/util.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/util.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include "SampleRotamersFromPDB.hh"
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/PoseResidueTypeSet.hh>
#include <boost/algorithm/string.hpp> //using this to parse up the rotamer db



#include <core/pack/task/operation/task_op_schemas.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


#ifdef WIN32
#include <core/graph/Graph.hh>
#endif

namespace protocols {
namespace splice {
typedef utility::vector1<utility::vector1<core::conformation::ResidueOP> > res_matrix;
static  basic::Tracer TR("core.pack.rotamer_set.SampleRotamersFromPDB_RotamerSetOperation");

SampleRotamersFromPDB_RotamerSetOperation::SampleRotamersFromPDB_RotamerSetOperation() :
	RotamerSetOperation(),
	add_rotamer_(false),
	SampleAtAlignedpositions_(),
	debug_(false),
	ccd_(1),
	db_file_(""),
	resi_vec_()
{
}

SampleRotamersFromPDB_RotamerSetOperation::SampleRotamersFromPDB_RotamerSetOperation(bool add_rot, utility::vector1<
	core::Size> SampleAtAlignedpositions,bool d,bool ccd, std::string file_name) :
	RotamerSetOperation(),
	add_rotamer_(add_rot),
	SampleAtAlignedpositions_(SampleAtAlignedpositions),
	debug_(d),
	ccd_(ccd),
	db_file_(file_name),
	resi_vec_() {
}

SampleRotamersFromPDB_RotamerSetOperation::~SampleRotamersFromPDB_RotamerSetOperation() {
}

///////////////////////////// HELPER FUNCTIONs ///////////////////////////
//@ breif This function prints the chi values of a given residue
void printChi(core::conformation::Residue const res) {
	utility::vector1<core::Real> chi = res.chi();
	TR << res.name3() << ":" << chi << std::endl;
}

//To avoid redundancy in the rotamer database I will not store Rotamers that are identical
bool
is_identical_rotamer( ROT const existing_res,ROT const new_res )
{
	//TR<<"Searching for indentical residue"<<std::endl;
	bool match = true;
	if ( existing_res.chi_vec.size() != new_res.chi_vec.size() || existing_res.AA != new_res.AA  ) {
		return false;
	} else {
		for ( core::Size i = 1; i<= existing_res.chi_vec.size(); ++i ) {
			if ( std::abs( existing_res.chi_vec[i] - new_res.chi_vec[i]) >= 5 ) {
				match = false;
			}
		}
	}
	// if (match) TR<<ex/isting_res.chi_vec<<":"<<new_res.chi_vec<<std::endl;
	return match;
}
core::pack::rotamer_set::RotamerSetOperationOP SampleRotamersFromPDB_RotamerSetOperation::clone() const {
	return core::pack::rotamer_set::RotamerSetOperationOP(new SampleRotamersFromPDB_RotamerSetOperation(*this));
}


void SampleRotamersFromPDB_RotamerSetOperation::add_pose(core::pose::PoseCOP pose) {
	poses_.push_back(pose);
}

/*void  SampleRotamersFromPDB_RotamerSetOperation::add_cutpoint_variant(core::conformation::ResidueOP db_rot, core::conformation::ResidueOP pose_rot){
if (pose_rot->has_variant_type(core::chemical::CUTPOINT_LOWER)){
add_variant_type_to_pose_residue(pose, CUTPOINT_LOWER, cutpoint );
}
}*/
bool SampleRotamersFromPDB_RotamerSetOperation::is_residue_allowed(core::chemical::ResidueType const & restype, core::pack::task::ResidueLevelTask const & rtask) {
	bool allowed = false;
	for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter j = rtask.allowed_residue_types_begin(),
			j_end = rtask.allowed_residue_types_end(); j != j_end; ++j ) {
		//TR<<restype.name3()<<"=="<<(**j).name3()<<std::endl;
		if ( restype.name3() == (**j).name3() ) {
			allowed = true;
			//TR << "Residue type from source is: " << restype << std::endl;
			break;
		}
	}
	return allowed;
}

void SampleRotamersFromPDB_RotamerSetOperation::add_rotamer_to_rotamer_set(core::pose::Pose const & pose, core::pack::rotamer_set::RotamerSet & rotamer_set, core::conformation::ResidueOP cur_rot) {
	core::Size seqnum = (core::Size) rotamer_set.resid();
	TR<<"Adding Rotamer!"<<std::endl;
	TR<<"rotamer_set size before addition:"<<rotamer_set.num_rotamers()<<std::endl;
	core::conformation::ResidueOP cur_res = cur_rot->clone();
	core::conformation::Residue const & existing_residue = pose.residue(seqnum);
	cur_res->place(existing_residue, pose.conformation());
	cur_res->seqpos(seqnum);
	cur_res->chain(existing_residue.chain());
	cur_res->copy_residue_connections_from(existing_residue);
	rotamer_set.add_rotamer_into_existing_group(*cur_res);
	TR<<"rotamer_set size after addition:"<<rotamer_set.num_rotamers()<<std::endl;
}
//@brief rather than filtering rotamers that match those in the db file I can rescore them to be favoured by the packer
void SampleRotamersFromPDB_RotamerSetOperation::add_rotamer_constraints(Pose & pose, core::Size seqnum, core::conformation::ResidueOP cur_rot) {
	core::pack::dunbrack::RotamerConstraintOP constraint(new core::pack::dunbrack::RotamerConstraint(pose, seqnum));
	constraint->add_residue(*cur_rot);
	pose.add_constraint(constraint);
}

void SampleRotamersFromPDB_RotamerSetOperation::initialize_from_command_line() {
	using namespace basic::options;
	if ( !option[OptionKeys::packing::unboundrot].active() ) {
		return;
	}
	for ( Size i = 1; i <= option[OptionKeys::packing::unboundrot]().size(); ++i ) {
		std::string filename = option[OptionKeys::packing::unboundrot]()[i].name();
		TR << "Adding rotamers from " << filename << std::endl;
		core::pose::PoseOP pose(new core::pose::Pose());
		//core::import_pose::pose_from_pdb( *pose, filename );
		pose=core::import_pose::pose_from_file(filename);
		this->add_pose(pose);
	}
}
void SampleRotamersFromPDB_RotamerSetOperation::copy_rotamer_matrix(rot_matrix const & rm) {
	//resi_vec_.resize(rm.size());
	for ( core::Size resi_num = 1; resi_num <= rm.size(); ++resi_num ) {
		resi_vec_.push_back(rm[resi_num]);
	}
}
// given a rotamer db file create a matrix
void SampleRotamersFromPDB_RotamerSetOperation::fill_rotamer_matrix_from_db_file() {
	TR << "Filling rotamers from DB!" << std::endl;
	utility::vector1<utility::vector1<std::string> > resi_vec;
	//resi_vec.resize()
	std::string line;
	std::ifstream db_in;
	db_in.open(db_file_.c_str(), std::ios::in);

	if ( db_in.is_open() ) {
		core::Size line_num = 1;

		while ( (getline(db_in, line)) && (!line.empty()) ) {
			utility::vector1<std::string> Rots;
			boost::algorithm::split(Rots, line, boost::is_any_of(","), boost::token_compress_on);
			//TR<<line<<std::endl;
			Rots.erase(Rots.begin());   //remove residue number
			Rots.pop_back();   //remove trailing newline
			//TR<<Rots.back()<<std::endl;
			resi_vec.push_back(Rots);
			//TR<<resi_vec[line_num][2]<<std::endl;
			line_num++;
		}
		db_in.close();
	} else {
		utility_exit_with_message("Cannot open rotamer file: " + db_file_ + "\n");
	}

	//////////Create_Rotamers_from_text///////////
	using namespace core::chemical;
	using namespace core::conformation;
	resi_vec_.resize(resi_vec.size());
	for ( core::Size resi_num = 1; resi_num <= resi_vec_.size(); ++resi_num ) {
		for ( core::Size rot = 1; rot <= resi_vec[resi_num].size(); rot++ ) {
			utility::vector1<std::string> Rots;
			boost::split(Rots, resi_vec[resi_num][rot], boost::is_any_of(" "), boost::token_compress_on);
			//TR<<Rots[1]<<std::endl;
			// ResidueTypeSet const & residue_set(pose.residue(1).residue_type_set());
			ROT new_res;
			//ResidueOP new_res = ResidueFactory::create_residue(residue_set.name_map(name_from_aa(aa_from_oneletter_code(*Rots[1].c_str()))));
			new_res.AA = *Rots[1].c_str();
			Rots.erase(Rots.begin());
			utility::vector1<core::Real> Rots_real;
			for ( core::Size st = 1; st <= Rots.size(); st++ ) {
				new_res.chi_vec.push_back(std::atof(Rots[st].c_str()));
			}
			//new_res->set_all_chi(Rots_real);
			bool res_exsits_in_db = false;
			for ( core::Size i=1; i<resi_vec_[resi_num].size(); i++ ) {
				res_exsits_in_db = is_identical_rotamer(resi_vec_[resi_num][i],new_res);
				if ( res_exsits_in_db ) break;
				//TR<<"Found identical residue"<<std::endl;
			}//for

			if ( !res_exsits_in_db ) resi_vec_[resi_num].push_back(new_res);

		}//for rot
	}//for resi_num
}
void SampleRotamersFromPDB_RotamerSetOperation::alter_rotamer_set_from_db(core::pose::Pose const & pose, core::pack::task::PackerTask const & ptask, core::pack::rotamer_set::RotamerSet & rotamer_set) {
	using namespace core;

	//TR<<"Size of SampleAtAlignedpositions_ is: "<<SampleAtAlignedpositions_.size()<<std::endl;
	Size const seqnum = (Size) rotamer_set.resid();  //seqnum holds the the position of the current rotamer vector
	if ( pose.residue(seqnum).chain() > 1 ) {
		return; //I can only apply this to chain 1
	}
	//if residue is a cut point don't change rotamer vector, it causes problem with CCD
	//TR<<"CCD val is:"<<ccd_<<std::endl;
	if ( (ccd_)&&(pose.residue(seqnum).has_variant_type( core::chemical::CUTPOINT_UPPER)||pose.residue(seqnum).has_variant_type( core::chemical::CUTPOINT_LOWER)) ) {
		TR<<"Rseidue is cut"<<std::endl;
		return;
	}
	TR << "Now looking at seqnum: " << seqnum << std::endl;
	core::pack::task::ResidueLevelTask const & rtask = ptask.residue_task(seqnum);
	// core::chemical::ResidueType const & restype = pose.residue_type(seqnum);
	utility::vector1<bool> rotamers_to_delete( rotamer_set.num_rotamers(), true ); //by default I am deleting all the rotamers
	TR<<"Allowed residues: ";
	for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter j = rtask.allowed_residue_types_begin(),
			j_end = rtask.allowed_residue_types_end(); j != j_end; ++j ) {
		TR << (**j).name1()<<",";
	}
	TR<<std::endl;
	TR << "Original size of Rosetta calc. rotamer vector " << rotamer_set.num_rotamers() << std::endl;
	//of the rotmer vector of the specific position
	Size irot(0), num_to_delete(0);

	TR << "Size of rotamer vector from db file: " << resi_vec_[seqnum].size() << std::endl;

	std::list< core::conformation::ResidueOP > new_rotamers;

	runtime_assert(pose.split_by_chain(1)->total_residue() == resi_vec_.size());//assume we only applying to chain 1
	utility::vector1< core::conformation::ResidueOP > temp_res_vec;
	for ( Size i = 1; i <= resi_vec_[seqnum].size(); ++i ) {    //go over all rotamers in DB rotamer vector
		core::conformation::ResidueOP cur_rot = ROT2res(resi_vec_[seqnum][i],pose);    //rotamer from db
		//if residue is not allowed skip iteration
		if ( !(is_residue_allowed(cur_rot->type(),rtask)) ) {
			continue;
		}
		//TR<<"Residue is allowed"<<std::endl;
		core::Size count = 0;
		irot = 0;
		for ( Rotamers::const_iterator it = rotamer_set.begin(), ite = rotamer_set.end(); it != ite; ++it ) {
			count++;
			irot++;    //moving irot to the first position of the rotamers_to_delete vecotr1
			core::conformation::ResidueOP rop(*it); //rotamer from rotamer set of pose
			if ( rop->is_similar_rotamer(*cur_rot) ) {
				if ( debug_ ) {
					TR << "The DB rotamer for residue" << seqnum << " is: ";
				}
				if ( debug_ ) {
					printChi(*cur_rot);
				}
				if ( debug_ ) {
					TR << "The PDB rotamer is: ";
				}
				if ( debug_ ) {
					printChi(*rop);
				}
				if ( debug_ ) {
					TR << " Found similar Rot! at rotamer vector postion " << irot << std::endl;
				}
				rotamers_to_delete[irot]=false;
			}
		}
		temp_res_vec.push_back(cur_rot);
		if ( add_rotamer_ ) {
			//Add pose rotamer to rotamer vector
			new_rotamers.push_back( cur_rot );
		}
	} //for resi_vec_
	num_to_delete = std::count(rotamers_to_delete.begin(), rotamers_to_delete.end(), true);
	TR << "number of rotamers to delete is :" << num_to_delete << std::endl;
	//After going over all rotamers now append PDB rotmaers to rotamer vector
	//flo nov 2010: if all the rotamers in the rotamer set are to be deleted,
	//and there is no input rotamer at that position, this will cause a program exit
	//if this is the case (rare), it's probably best to not delete any rotamers
	TR << rotamers_to_delete << std::endl;
	if ( num_to_delete == rotamer_set.num_rotamers() )  {

		TR.Error << "Did not find similar rotamers in DB rotamer vector " << rotamer_set.resid()  << std::endl;
		if ( !add_rotamer_ ) {
			return;//if adding rotamer I will never be in a probelm where I'm delteing all the rotamers in the vector //gideonla,may17
		}
		/*
		if (temp_res_vec.size()==0){
		std::cerr <<"Not deleteing any rotamers" <<std::endl;
		return;
		}
		TR << "SIZE of temp_res_vec: "<< temp_res_vec.size()<<std::endl;
		for ( core::Size i =1; i<=temp_res_vec.size(); ++i ) {
		this->add_rotamer_to_rotamer_set(pose, rotamer_set, temp_res_vec[i]);
		rotamers_to_delete.push_back(false);
		}
		std::cerr << ", Using only DB vector rotamers!" << std::endl;
		TR << rotamers_to_delete << std::endl;
		*/
	}
	for ( auto const & rot : new_rotamers ) {
		// Note: while the rotamer_set is being added to, you do not want to ask for (even to print out)
		// any of its data, as that will trigger an O(N) integration of the new rotamers into the set;
		// If you were to do this in every iteration of this loop, you would pay O(N^2) for something
		// you do not want or need.
		add_rotamer_to_rotamer_set( pose, rotamer_set, rot );
	}
	utility::vector1<bool> rotamers_to_delete2( rotamer_set.num_rotamers(), true ); //by default I am deleting all the rotamer
	for ( Size i = 1; i <= resi_vec_[seqnum].size(); ++i ) {    //go over all rotamers in DB rotamer vector
		core::conformation::ResidueOP cur_rot = ROT2res(resi_vec_[seqnum][i],pose);    //rotamer from db
		//if residue is not allowed skip iteration
		if ( !(is_residue_allowed(cur_rot->type(),rtask)) ) {
			continue;
		}
		//TR<<"Residue is allowed"<<std::endl;
		core::Size count = 0;
		irot = 0;
		for ( Rotamers::const_iterator it = rotamer_set.begin(), ite = rotamer_set.end(); it != ite; ++it ) {
			count++;
			irot++;    //moving irot to the first position of the rotamers_to_delete vecotr1
			core::conformation::ResidueOP rop(*it); //rotamer from rotamer set of pose
			if ( rop->is_similar_rotamer(*cur_rot) ) {
				rotamers_to_delete2[irot]=false;
			}
		}
	} //for resi_vec_
	Size num_to_delete2( std::count(rotamers_to_delete2.begin(), rotamers_to_delete2.end(), true));
	TR<<"rotamers_to_delete size:"<<num_to_delete2<<std::endl;
	TR<<"rotamer_set size:"<<rotamer_set.num_rotamers()<<std::endl;
	if ( (rotamer_set.num_rotamers() ==  num_to_delete2) ) {
		return;
	}
	rotamer_set.drop_rotamers(rotamers_to_delete2);
	return;

}
void SampleRotamersFromPDB_RotamerSetOperation::alter_rotamer_set_from_pdb(core::pose::Pose const & pose, core::pack::task::PackerTask const & ptask, core::pack::rotamer_set::RotamerSet & rotamer_set) {
	using namespace core;

	//TR<<"Size of SampleAtAlignedpositions_ is: "<<SampleAtAlignedpositions_.size()<<std::endl;
	Size const seqnum = (Size) rotamer_set.resid();  //seqnum holds the the position of the current rotamer vector
	if ( pose.residue(seqnum).chain() > 1 ) {
		return; //I can only apply this to chain 1
	}
	if ( (ccd_)&&(pose.residue(seqnum).has_variant_type( core::chemical::CUTPOINT_UPPER)||pose.residue(seqnum).has_variant_type( core::chemical::CUTPOINT_LOWER)) ) {
		TR<<"Rseidue is cut"<<std::endl;
		return;
	}
	if ( debug_ ) {
		TR << "Now looking at seqnum: " << seqnum << std::endl;
	}
	utility::vector1<bool> rotamers_to_delete;
	rotamers_to_delete.resize(rotamer_set.num_rotamers()); //resizing the size of the bool vector to be the size of the rotamer_set
	std::fill(rotamers_to_delete.begin(), rotamers_to_delete.end(), true); //by default I am deleting all the rotamers
	if ( debug_ ) {
		TR << "Original size of rotamer vector " << rotamer_set.num_rotamers() << std::endl;
	}
	core::pack::task::ResidueLevelTask const & rtask = ptask.residue_task(seqnum);
	//of the rotmer vector of the specific position
	Size irot(0), num_to_delete(0);

	for ( Size i = 1; i <= poses_.size(); ++i ) {
		core::pose::Pose const & ubr_pose = *(poses_[i]);
		//check the aligned positions on the pose compared to the refrence
		core::Size nearest_on_pose=0;
		core::Size nearest_on_source=0;

		bool found_aligned_postion_match = false;
		for ( core::Size const ref_pose_pos: SampleAtAlignedpositions_ ) {
			TR <<"Checking ref pos: "<<ref_pose_pos<<std::endl;
			if ( ref_pose_pos>ubr_pose.total_residue() ) continue; //if the refernece position if larger than the pose then just skip this posistion,Gideon11mar15
			nearest_on_pose = protocols::rosetta_scripts::find_nearest_res(pose,ubr_pose, ref_pose_pos, 0);
			TR<<"nearest_on_pose:"<<nearest_on_pose<<std::endl;
			if ( nearest_on_pose==seqnum ) {
				found_aligned_postion_match=true;
				nearest_on_source = ref_pose_pos;
			}
		}
		//TR<<"The value of found_aligned_postion_match ="<<found_aligned_postion_match<<std::endl;
		if ( !found_aligned_postion_match and SampleAtAlignedpositions_.size() ) {
			return;
		}
		core::Size const nearest_on_ubr_pose(protocols::rosetta_scripts::find_nearest_res(ubr_pose, pose, seqnum, 0)); //find the closestset residue between the pose and ubr_pose
		if ( nearest_on_ubr_pose == 0 ) {
			continue; //if did not find on nearest pose from pose list then there is no point
		}
		if ( debug_ ) {
			TR << "Nearest residue on " << ubr_pose.pdb_info()->name() << " is " << nearest_on_ubr_pose << std::endl;
		}
		core::chemical::ResidueType const & restype = ubr_pose.residue_type(nearest_on_source);
		TR<<"Allowed residues: ";
		for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter j = rtask.allowed_residue_types_begin(),
				j_end = rtask.allowed_residue_types_end(); j != j_end; ++j ) {
			TR << (**j).name1()<<",";
		}
		TR<<std::endl;
		TR<<restype<<std::endl;
		TR<<is_residue_allowed(restype, rtask)<<std::endl;
		if ( is_residue_allowed(restype, rtask) ) {
			core::conformation::ResidueOP cur_rot;
			core::Size count = 0;
			irot = 0;
			for ( Rotamers::const_iterator it = rotamer_set.begin(), ite = rotamer_set.end(); it != ite; ++it ) {
				count++;
				irot++; //moving irot to the first position of the rotamers_to_delete vecotr1
				cur_rot = ubr_pose.residue(nearest_on_ubr_pose).clone();
				core::conformation::ResidueOP rop(*it);
				if ( rop->is_similar_rotamer(*cur_rot) ) {
					TR << "The PDB rotamer is: ";
					printChi(*cur_rot);
					TR << "The DB rotamer is: ";
					printChi(*rop);
					TR << " Found similar Rot! at rotamer vector postion " << irot << std::endl;
					rotamers_to_delete[irot] = false;
				}
			}
			if ( add_rotamer_ ) {
				//Add pose rotamer to rotamer vector
				this->add_rotamer_to_rotamer_set(pose, rotamer_set, cur_rot);
				rotamers_to_delete.push_back(false);
			}
		} else { //fi type_is_allowed
			TR << "Residue is not allowed" << std::endl;
		}
	} //for(Size i = 1; i <= poses_.size(); ++i)
	num_to_delete = std::count(rotamers_to_delete.begin(), rotamers_to_delete.end(), true);
	TR << "number of rotamers for current res before deletion:" << rotamer_set.num_rotamers() << std::endl;
	TR << "number of rotamers to delete is :" << num_to_delete << std::endl;
	//After going over all rotamers now append PDB rotmaers to rotamer vector
	//flo nov 2010: if all the rotamers in the rotamer set are to be deleted,
	//and there is no input rotamer at that position, this will cause a program exit
	//if this is the case (rare), it's probably best to not delete any rotamers
	TR << rotamers_to_delete << std::endl;
	if ( (num_to_delete == rotamer_set.num_rotamers()) && (rotamer_set.id_for_current_rotamer() == 0) ) {
		//std::cerr << "shit condition at position " << rotamer_set.resid() << ", not deleting any of the " << rotamer_set.num_rotamers() << " rotamers." << std::endl;
		return;
	}
	rotamer_set.drop_rotamers(rotamers_to_delete);
	return;

}
void SampleRotamersFromPDB_RotamerSetOperation::alter_rotamer_set(core::pose::Pose const & pose, core::scoring::ScoreFunction const & /*sfxn*/, core::pack::task::PackerTask const & ptask, utility::graph::GraphCOP /*packer_neighbor_graph*/, core::pack::rotamer_set::RotamerSet & rotamer_set) {
	//TR<<"Size of resi_vec_: "<<resi_vec_.size()<<std::endl;
	//TR << "The address of resi_vec_: " << &resi_vec_ << std::endl;
	if ( poses_.size() > 0 ) {
		alter_rotamer_set_from_pdb(pose, ptask, rotamer_set);
	} else {
		alter_rotamer_set_from_db(pose, ptask, rotamer_set);
	}
} // alter_rotamer_set

core::conformation::ResidueOP SampleRotamersFromPDB_RotamerSetOperation::ROT2res(ROT rot, core::pose::Pose const & pose){
	using namespace core::chemical;
	using namespace core::conformation;


	PoseResidueTypeSetOP residue_set= pose.conformation().modifiable_residue_type_set_for_conf( core::chemical::FULL_ATOM_t );
	ResidueOP new_res = ResidueFactory::create_residue(residue_set->name_map(name_from_aa(aa_from_oneletter_code(rot.AA[0]))));
	new_res->set_all_chi(rot.chi_vec);
	return new_res;
}
///////////////////////////////////////////////////////////////////////////////////////////
core::pack::task::operation::TaskOperationOP SampleRotamersFromPDBCreator::create_task_operation() const {
	return core::pack::task::operation::TaskOperationOP(new SampleRotamersFromPDB);
}

/// @brief default constructor
SampleRotamersFromPDB::SampleRotamersFromPDB() :
	TaskOperation() {

}

/// @brief destructor
SampleRotamersFromPDB::~SampleRotamersFromPDB() {
}

/// @brief clone
core::pack::task::operation::TaskOperationOP SampleRotamersFromPDB::clone() const {
	return core::pack::task::operation::TaskOperationOP(new SampleRotamersFromPDB(*this));
}

void SampleRotamersFromPDBCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SampleRotamersFromPDB::provide_xml_schema( xsd );
}

std::string SampleRotamersFromPDBCreator::keyname() const
{
	return SampleRotamersFromPDB::keyname();
}

/// @brief
void SampleRotamersFromPDB::apply(Pose const & pose, PackerTask & task) const {
	//TR<<ROTdb_segments_["frm1"]->pdb_profile("2VC5")[1][1]->name()<<std::endl;
	//TR<<ROTdb_segments_.find("frm1")->second->pdb_profile("2VC5")[1][1]->name()<<std::endl;
	SampleRotamersFromPDB_RotamerSetOperationOP rso(new SampleRotamersFromPDB_RotamerSetOperation(add_rotamer_, SampleAtAlignedpositions(), debug_, ccd_,db_fname_));


	//So we have 3 different methods to fill the rotamer matrix :1) from a single db file (in case we are applying this to a single strcture),
	//2) In the case of splice we can combine different rotamer db files
	//3) copy rotmaers from input pdbs
	if ( db_fname_ != "" ) {
		rso->fill_rotamer_matrix_from_db_file();
		task.append_rotamerset_operation(rso);
	} else if ( segment_names_ordered_.size() > 0 ) {
		//if segment subtags are specfied in the xml
		rso->copy_rotamer_matrix(this->combine_rot_dbs(pose));
		task.append_rotamerset_operation(rso);
	} else {
		rso->initialize_from_command_line();
	}
	task.append_rotamerset_operation(rso);
	//if(debug_) pose.dump_pdb("sample_rotamer.pdb");
}

void SampleRotamersFromPDB::parse_tag(TagCOP tag, DataMap &) {
	db_fname_ = tag->getOption < std::string > ("db_file_name", "");
	add_rotamer_ = tag->getOption<bool>("add_rotamer", false);//if set to true rotamer from db file or pdb are added to pose residue vector
	std::string const SampleAtAlignedpositions(tag->getOption < std::string > ("aligned_positions", ""));
	//parse SampleAtAlignedpositions
	utility::vector1<std::string> const split_reslist(utility::string_split(SampleAtAlignedpositions, ','));
	for ( std::string const &res_str: split_reslist ) {
		if ( res_str=="" ) break; //with this no residue numbers given an empty string causes an error here. Gideon,11Mar15
		using namespace std;
		TR<<"res_str:"<<res_str<<std::endl;
		int num;
		istringstream ( res_str ) >> num;
		TR<<"num is:"<<num<<std::endl;
		SampleAtAlignedpositions_.push_back( num );
	}
	debug_ = tag->getOption<bool>("debug", false); //if set to true throws more output
	ccd_ = tag->getOption<bool>("ccd", false); //if set to true than cutpoint residues are not modified

	utility::vector1<TagCOP> const sub_tags(tag->getTags());
	core::pose::PoseCOP pose(protocols::jd2::get_current_jobs_starting_pose());
	for ( TagCOP sub_tag: sub_tags ) {
		if ( sub_tag->getName() == "Segments" ) {
			if ( db_fname_!="" ) {
				utility_exit_with_message( "Cannot use \"db_file_name\" with sub tags \"Segments\"");
			}
			utility::vector1< TagCOP > const segment_tags( sub_tag->getTags() );
			for ( TagCOP segment_tag: segment_tags ) {
				RotLibdbOP RotLib_segment( new RotLibdb );
				std::string const segment_name( segment_tag->getOption< std::string >( "name" )  ); //get name of segment from xml
				std::string const pdb_profile_match( segment_tag->getOption< std::string >( "pdb_profile_match" ) );// get name of pdb profile match, this file contains all the matching between pdb name and sub segment name, i.e L1.1,L1.2 etc
				std::string const profiles_str( segment_tag->getOption< std::string >( "rot_lib" ) );
				typedef utility::vector1<std::string> StringVec;
				StringVec const profile_name_pairs( utility::string_split( profiles_str, ',' ) );

				// TR<<"Now working on segment:"<<segment_name<<std::endl;
				for ( std::string const & s: profile_name_pairs ) {
					StringVec const profile_name_file_name( utility::string_split( s, ':' ) );
					if ( debug_ ) TR<<"Rotamer DB file: "<<profile_name_file_name[ 2 ]<<",segment name: "<<profile_name_file_name[ 1 ]<<std::endl;

					RotLib_segment->add_segment_fname(profile_name_file_name[ 2 ], profile_name_file_name[ 1 ] );
				}
				RotLib_segment->read_pdb_profile( pdb_profile_match );
				/*TR<<"the segment name is: "<<segment_name<<std::endl;
				if (segment_name.compare(segment_type_) == 0) {
				check_segment=true;
				}*/
				ROTdb_segments_.insert( std::pair< std::string, RotLibdbOP >( segment_name, RotLib_segment ) );
				segment_names_ordered_.push_back(segment_name);
			}
			//foreach segment_tag
		} else { //fi
			utility_exit_with_message( "SampleRotamersFromPDB subtag not recognized: " + sub_tag->getName() );
		}
	} //foreach

	//TR<<ROTdb_segments_["frm1"]->pdb_profile("2VC5")[1][1]->name()<<std::endl;
}


//@brief using the pose comments the create a single rotamer matrix
rot_matrix SampleRotamersFromPDB::combine_rot_dbs(core::pose::Pose const & pose) const {
	using namespace std;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::jd2::JobDistributor;

	//TR<<ROTdb_segments_.find("frm1")->second->pdb_profile("2VC5")[1][1]->name()<<std::endl;


	std::map< std::string/*which segment (L1,L2...)*/, std::string/*pdb name*/ > pdb_segments_;
	rot_matrix profile_vector_;
	protocols::splice::load_pdb_segments_from_pose_comments(pose, pdb_segments_); // get segment name and pdb accosiation from comments in pdb file
	if ( pdb_segments_.size() != segment_names_ordered_.size() ) {
		utility_exit_with_message(
			"The number of segments in the input pose comments does not match the number of segments in the XML\n");
	} else {
		TR << "There are " << pdb_segments_.size() << " ROTlib segments" << std::endl;
	}

	runtime_assert(pdb_segments_.size()); //This assert is in place to make sure that the pdb file has the correct comments, otherwise this function will fail
	for ( std::string const & segment_type: segment_names_ordered_ ) { //<- Start of PDB segment iterator
		TR<<"segment_type: "<<segment_type<<std::endl;
		if ( ROTdb_segments_.find(segment_type )-> second->get_segment_fname(pdb_segments_[segment_type]).empty() ) {
			utility_exit_with_message(" could not find the source pdb name: "+ pdb_segments_[segment_type]+ ", in the pdb_profile_match \n");
		}
		TR<<"reading profile:"<< pdb_segments_[segment_type]<<std::endl;
		//profile_vector.push_back(ROTdb_segments_[ segment_type ]->pdb_profile( pdb_segments_[segment_type] )[1]);
		RotLibdb RotLib;
		rot_matrix tmp_rot_mat=RotLib.fill_rotamer_matrix_from_db_file(ROTdb_segments_.find(segment_type )-> second->get_segment_fname( pdb_segments_[segment_type] ));
		concatenate_rot_matrix(profile_vector_,tmp_rot_mat );
		TR << "The size of the profile vector is: " << profile_vector_.size() << std::endl;
		//TR<<ROTdb_segments_.find(segment_type )-> second->pdb_profile( pdb_segments_[segment_type] )[1][15]->name()<<std::endl;
	} // <- End of PDB segment iterator

	if ( profile_vector_.size() != pose.conformation().chain_end(1)  ) {
		utility_exit_with_message(
			"The number of residues in the first chain of the pose does not match the number of rotamer vectors in the rotamer db \n");
	}
	return profile_vector_;
}
void SampleRotamersFromPDB::concatenate_rot_matrix(rot_matrix & a, rot_matrix const & b) const{
	//for(std::vector<utility::vector1< core::conformation::ResidueOP > >::iterator it = b.begin(); it != b.end(); ++it) {
	for ( utility::vector1< ROT > const & i: b ) {
		TR<<i[1].AA<<i[1].chi_vec<<std::endl;;

		a.push_back(i);
	}
}

std::string Srfp_complex_type_name_for_subsubtag( std::string const & foo ) {
	return "subsubtag_Srfp_" + foo + "_type";
}

std::string Srfp_complex_type_name_for_subtag( std::string const & foo ) {
	return "subtag_Srfp_" + foo + "_type";
}
void SampleRotamersFromPDB::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using namespace core::pack::task::operation;
	XMLSchemaSimpleSubelementList subelements;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default(  "add_rotamer", xsct_rosetta_bool, "Which chain identifier" ,"1")
		+ XMLSchemaAttribute::attribute_w_default(  "debug", xsct_rosetta_bool, "make output more verbose" ,"0" )
		+ XMLSchemaAttribute::attribute_w_default(  "ccd", xsct_rosetta_bool, "change behavior if running this with CCD, default it true" ,"1" )
		+ XMLSchemaAttribute(  "aligned_positions", xsct_int_cslist, "which positions are aligned" );


	// The "Segments" subtag
	AttributeList segments_subtag_attlist;


	// The "segment" sub-subtag"
	AttributeList subtag_segments_subtag_attlist;
	subtag_segments_subtag_attlist + XMLSchemaAttribute( "pdb_profile_match", xs_string, "text file matching pdb name to profile" )
		+ XMLSchemaAttribute( "rot_lib", xs_string, "XRW TO DO" );

	XMLSchemaComplexTypeGenerator segment_subsubtag_gen;
	segment_subsubtag_gen.complex_type_naming_func( & Srfp_complex_type_name_for_subsubtag )
		.element_name( "Segment" )
		.description( "individual segment tag" )
		.add_attributes( subtag_segments_subtag_attlist )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList subsubelements;
	subsubelements.add_already_defined_subelement( "Segment", Srfp_complex_type_name_for_subsubtag/*, 0*/ );

	XMLSchemaComplexTypeGenerator segments_subtag_gen;
	segments_subtag_gen.complex_type_naming_func( & Srfp_complex_type_name_for_subtag )
		.element_name( "Segments" )
		.description( "Wrapper for multiple segments tags" )
		.add_attributes( segments_subtag_attlist )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subsubelements )
		.write_complex_type_to_schema( xsd );

	subelements.add_already_defined_subelement( "Segments", Srfp_complex_type_name_for_subtag/*, 0*/ );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.element_name( keyname() )
		.description( "Restrict rotamer selctions to ones derived from natural PDBs" )
		.complex_type_naming_func( & complex_type_name_for_task_op )
		.add_attributes( attlist )
		.set_subelements_repeatable( subelements )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );


}
///////////////////////////////////////////////////////////////////ROTLIBDB/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//@brief used for checking backbone segment compatibilty
RotLibdb::RotLibdb() :
	pdb_names_("") {
}
void RotLibdb::read_profile(std::string const & file_name, std::string const & segment_name) {
	rot_matrix res_mat;
	utility::file::FileName const fname(file_name);
	//TR<<"Segment name: "<<segment_name<<",reading sequence profile from "<<file_name<<std::endl;
	res_mat = this->fill_rotamer_matrix_from_db_file(file_name);
	//TR<< "SpliceSegmentsSeqProf:"<< new_seqprof->prof_row(1)<<std::endl;
	//TR<< "SpliceSegmentspprobabilty:"<< new_seqprof->probability_row(1)<<std::endl;
	rot_matrix_profile_.insert(std::pair<std::string, rot_matrix>(segment_name, res_mat));
}
void RotLibdb::add_segment_fname(std::string const & file_name, std::string const & segment_name) {
	ROTdb_segments_fileName_.insert(std::pair<std::string, std::string>(segment_name,file_name));
}
rot_matrix RotLibdb::fill_rotamer_matrix_from_db_file(std::string fname) {
	utility::vector1<utility::vector1<std::string> > resi_vec;
	rot_matrix resi_mat;
	//resi_vec.resize()
	std::string line;
	std::ifstream db_in;
	db_in.open(fname.c_str(), std::ios::in);

	if ( db_in.is_open() ) {
		core::Size line_num = 1;

		while ( (getline(db_in, line)) && (!line.empty()) ) {
			utility::vector1<std::string> Rots;
			boost::split(Rots, line, boost::is_any_of(","), boost::token_compress_on);
			//TR<<line<<std::endl;
			Rots.erase(Rots.begin());      //remove residue number
			Rots.pop_back();      //remove trailing newline
			//TR<<Rots.back()<<std::endl;
			resi_vec.push_back(Rots);
			//TR<<resi_vec[line_num][2]<<std::endl;
			line_num++;
		}
		db_in.close();
	} else {
		utility_exit_with_message("Cannot open rotamer db file: " + fname + "\n");
	}
	//////////Create_Rotamers_from_text///////////
	using namespace core::chemical;
	using namespace core::conformation;
	resi_mat.resize(resi_vec.size());
	for ( core::Size resi_num = 1; resi_num <= resi_mat.size(); ++resi_num ) {
		for ( core::Size rot = 1; rot <= resi_vec[resi_num].size(); rot++ ) {
			utility::vector1<std::string> Rots;
			boost::split(Rots, resi_vec[resi_num][rot], boost::is_any_of(" "), boost::token_compress_on);
			//TR<<Rots[1]<<std::endl;
			ROT new_res;;
			//ResidueOP new_res = ResidueFactory::create_residue(residue_set.name_map(name_from_aa(aa_from_oneletter_code(*Rots[1].c_str()))));
			new_res.AA = *Rots[1].c_str();
			Rots.erase(Rots.begin());
			utility::vector1<core::Real> Rots_real;
			for ( core::Size st = 1; st <= Rots.size(); st++ ) {
				new_res.chi_vec.push_back(std::atof(Rots[st].c_str()));
			}

			bool res_exsits_in_db = false;
			for ( core::Size i=1; i<resi_mat[resi_num].size(); i++ ) {
				res_exsits_in_db = is_identical_rotamer(resi_mat[resi_num][i],new_res);
				if ( res_exsits_in_db ) break;
				//TR<<"Found identical residue"<<std::endl;
			}//for

			if ( !res_exsits_in_db ) resi_mat[resi_num].push_back(new_res);
		}
	}      //for res_num
	return resi_mat;
}
void RotLibdb::read_pdb_profile(std::string const & file_name) {
	utility::io::izstream data(file_name);
	if ( !data ) {
		utility_exit_with_message("File not found " + file_name);
	}
	//TR<<"Loading pdb profile pairs from file "<<file_name<<std::endl;
	std::string line;
	while ( getline(data, line) ) {
		std::istringstream line_stream(line);
		while ( !line_stream.eof() ) {
			std::string pdb, profile;
			line_stream >> pdb >> profile;
			pdb_to_profile_map_.insert(std::pair<std::string, std::string>(pdb, profile));
			//TR<<"Loading pdb-profile pair: "<<pdb<<" "<<profile<<std::endl;
		}
	}

}
rot_matrix RotLibdb::pdb_profile(std::string const & pdb_name) {
	//TR<<"size of sequence_Profile_ is:"<<sequence_profile_.size()<<std::endl;
	//TR<<"pdb to profile map is:"<<pdb_name<<":"<<pdb_to_profile_map_[pdb_name]<<std::endl;
	return rot_matrix_profile_[pdb_to_profile_map_[pdb_name]];
}

std::string RotLibdb::get_segment_fname(std::string const & segment_name) {

	return ROTdb_segments_fileName_[segment_name];
}


} // splice
} // protocols

