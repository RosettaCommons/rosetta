// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RotLibOut.cc
/// @brief generate a database of rotamers from a list of pdbs
/// @author Gideon Lapidoth (glapidoth@gmail.com)

// Unit headers
#include <protocols/splice/RotLibOut.hh>
#include <protocols/splice/RotLibOutCreator.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
//using basic::T;
using basic::Error;
using basic::Warning;
static  basic::Tracer TR( "protocols.simple_moves.RotLibOut" );
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <fstream>
#include <vector>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <core/io/pdb/build_pose_as_is.hh>

#include <core/pose/PDBInfo.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <numeric/xyzVector.hh>
#include <algorithm>
#include <stdio.h>
#include <ctype.h>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace splice {


using namespace::protocols;

std::string
RotLibOutCreator::keyname() const
{
	return RotLibOutCreator::mover_name();
}

protocols::moves::MoverOP
RotLibOutCreator::create_mover() const {
	return protocols::moves::MoverOP( new RotLibOut );
}

std::string
RotLibOutCreator::mover_name()
{
	return "RotLibOut";
}

RotLibOut::RotLibOut(): moves::Mover("RotLibOut"), min_dist_(4.0)
{}

RotLibOut::~RotLibOut(){}

void RotLibOut::set_min_dist( core::Real min_dist){
	min_dist_ = min_dist;
}

core::Real RotLibOut::get_min_dist(){
	return min_dist_;
}

bool atom_dist_with_zero( core::pose::Pose & pose, core::pose::Pose & hit_pose, core::Size cur_target, core::Size cur_nearest ,core::Size min_dist) {
	if ( cur_nearest == 0 ) return false;
	//TR<<"cur_target="<<cur_target<<std::endl;
	//TR<<"cur_nearest="<<cur_nearest<<std::endl;
	//TR<<hit_pose.residue(1)<<std::endl;
	TR<<" template residue "<<cur_target<<" phi psi are: "<<pose.phi(cur_target)<<" "<<pose.psi(cur_target)<<std::endl;
	TR<<" hit residue "<<cur_nearest<<" phi psi are: "<<hit_pose.phi(cur_nearest)<<" "<<hit_pose.psi(cur_nearest)<<std::endl;
	TR<<"Distance between CA atoms of template residue "<<cur_target<<" and target residue "<<cur_nearest<<" is: "<<pose.residue( cur_target ).xyz( "CA" ).distance( hit_pose.residue( cur_nearest ).xyz( "CA" )) <<std::endl;
	TR<<"Distance between N atoms of template residue "<<cur_target<<" and target residue "<<cur_nearest<<" is: "<<pose.residue( cur_target ).xyz( "N" ).distance( hit_pose.residue( cur_nearest ).xyz( "N" )) <<std::endl;
	TR<<"Distance between O atoms of template residue "<<cur_target<<" and target residue "<<cur_nearest<<" is: "<<pose.residue( cur_target ).xyz( "O" ).distance( hit_pose.residue( cur_nearest ).xyz( "O" )) <<std::endl;
	TR<<"Distance between C atoms of template residue "<<cur_target<<" and target residue "<<cur_nearest<<" is: "<<pose.residue( cur_target ).xyz( "C" ).distance( hit_pose.residue( cur_nearest ).xyz( "C" )) <<std::endl;

	return ((pose.residue( cur_target ).xyz( "CA" ).distance( hit_pose.residue( cur_nearest ).xyz( "CA" )) < min_dist )  and
		(pose.residue( cur_target ).xyz( "N" ).distance( hit_pose.residue( cur_nearest ).xyz( "N" )) < min_dist ) and
		(pose.residue( cur_target ).xyz( "O" ).distance( hit_pose.residue( cur_nearest ).xyz( "O" )) < min_dist ) and
		(pose.residue( cur_target ).xyz( "C" ).distance( hit_pose.residue( cur_nearest ).xyz( "C" )) < min_dist ));
}

void
RotLibOut::find_matching_res( core::pose::Pose & pose, core::pose::Pose & hit_pose){
	// core::Size cur_target = 1; // the current position on the target examined
	// core::Size cur_nearest = rosetta_scripts::nearest_residue(pose with cur_target, hit_pose);
	// core::Size prev_nearest = cur_nearest;
	// resi_vec[ cur_target ].push_back( hit_pose[ cur_nearest ] );
	// ++cur_target;

	core::Size prev_hit_pose=0;
	core::Size cur_target=0, cur_nearest=0;//uninitialized error, Gideon Jun16
	std::string base_filename = hit_pose.pdb_info()->name().substr(hit_pose.pdb_info()->name().find_last_of("/\\") + 1);
	for ( core::Size j = 1; j <= pose.total_residue(); ++j ) {

		cur_nearest = rosetta_scripts::find_nearest_res(hit_pose,pose, j,0 /*search all chains*/);
		if ( cur_nearest != 0 ) {
			//cur_target = j + 1;
			//cur_nearest = cur_nearest + 1;
			cur_target = j;
			TR<<"cur_target="<<cur_target<<std::endl;
			TR<<"cur_nearest="<<cur_nearest<<std::endl;
			//prev_target_pose = cur_target;
			prev_hit_pose = cur_nearest;
			break;
		}
	}

	TR << "finished for, found cur_target " << cur_target << " cur_nearest " << cur_nearest << std::endl;

	while ( cur_target <= pose.total_residue() ) {
		core::Size current_frag_length=0;
		utility::vector1 <utility::vector1< core::conformation::ResidueOP > > resi_vec_tmp;
		resi_vec_tmp.resize(pose.total_residue());
		std::map< std::string, utility::vector1< std::string > > AA_vec_tmp;
		AA_vec_tmp[base_filename].resize(pose.total_residue());

		while ( (cur_nearest<=hit_pose.total_residue()) && (atom_dist_with_zero(pose, hit_pose, cur_target, cur_nearest,min_dist_))  ) {
			//TR<<"cur_nearest "<<cur_nearest<<std::endl;
			core::conformation::ResidueOP cur_nearest_res = hit_pose.residue(cur_nearest).clone();
			//resi_vec[cur_target].push_back( cur_nearest_res );
			resi_vec_tmp[cur_target].push_back( cur_nearest_res );
			//AA_vec[base_filename][cur_target] = hit_pose.residue(cur_nearest).name1();
			AA_vec_tmp[base_filename][cur_target] = hit_pose.residue(cur_nearest).name1();
			++cur_target;
			++cur_nearest;
			++current_frag_length;
			if ( cur_nearest>hit_pose.total_residue()||cur_target>pose.total_residue() ) break;
		}
		if ( current_frag_length>=min_frag_length() ) {
			//TR<<"current_frag_length="<<current_frag_length<<std::endl;
			for ( core::Size i=cur_target-current_frag_length; i<cur_target; i++ ) {
				//TR<<i<<std::endl;
				//TR<<AA_vec_tmp[base_filename][i]<<std::endl;
				resi_vec[i].push_back( resi_vec_tmp[i][1] );
				AA_vec[base_filename][i] = AA_vec_tmp[base_filename][i];
			}
		}//fi (current_frag_length
		//For debugging//
		TR<<"cur_nearest:" <<cur_nearest<<std::endl;
		TR<<"prev_hit_pose:" <<prev_hit_pose<<std::endl;
		if ( ((int)cur_nearest-(int)prev_hit_pose>1)&&(prev_hit_pose>0)&&(cur_nearest>0) &&dump_pdb() ) {
			TR<<"dumping pdb"<<std::endl;
			core::pose::Pose aligned_frag( hit_pose, prev_hit_pose , cur_nearest);
			//aligned_frag->delete_residue_range_slow(cur_nearest+1,aligned_frag->total_residue());
			std::string seqpos1,seqpos2;
			std::ostringstream convert;
			convert << cur_nearest; // insert the textual representation of 'Number' in the characters in the stream
			seqpos2 = convert.str();
			convert.str("");
			convert.clear();
			convert << prev_hit_pose; // insert the textual representation of 'Number' in the characters in the stream
			seqpos1 = convert.str();
			aligned_frag.dump_pdb(seqpos1+"_"+seqpos2+"_test.pdb");
			//prev_target_pose = cur_target;
		}
		prev_hit_pose  = cur_nearest;

		if ( cur_nearest>hit_pose.total_residue()||cur_target>pose.total_residue() ) break;
		//  erase_residue(n-ter, cur_nearest,hit_pose);
		//  int delta=cur_nearest - prev_hit_pose-1;
		//  TR<<"prev_hit_pose:"<<prev_hit_pose<<std::endl;
		//  prev_hit_pose = cur_nearest-delta;
		//  TR<<"delta="<<delta<<std::endl;
		//if (cur_nearest>0) hit_pose.conformation().delete_residue_range_slow( 1,  std::max((int) cur_nearest-1,1)); //remove residues we have gone over

		cur_nearest = rosetta_scripts::find_nearest_res(hit_pose, pose,cur_target,0);
		TR<<"cur_target:"<<cur_target<<std::endl;
		TR<<"cur_nearest:"<<cur_nearest<<std::endl;
		if ( ( cur_nearest == 0 ) || !(atom_dist_with_zero(pose, hit_pose, cur_target, cur_nearest,min_dist_)) ) ++ cur_target; // if a match on target is not found increment the current res on the target
	}//while
	//if the alignmnet is empty (i.e. not letters  just "-") then don't keep alignmnet
	bool good_aln=false;
	for ( core::Size res = 1; res <= AA_vec[base_filename].size() ; ++res ) {
		if ( isalpha(AA_vec[base_filename][res][0]) ) {
			good_aln = true;
		}
	}
	if ( !good_aln ) AA_vec.erase (base_filename);  //if seqeunce vectir is empty remove entry from map
}

//This function as not thread safe as only one thread can write to a file
void RotLibOut::print_lib_to_file() {
	std::ofstream db_file;
	std::stringstream seq_stringstream;
	db_file.open(Rotamer_db_file().c_str(), std::ios::app);
	for ( core::Size res_i = 1; res_i <= resi_vec.size(); ++res_i ) {//go over all the residues in the input pose
		seq_stringstream<<res_i<<" ,";
		for ( core::Size rot_i = 1; rot_i <= resi_vec[ res_i ].size(); ++rot_i ) {//
			seq_stringstream << resi_vec[ res_i ][ rot_i ]->name1()<<" ";
			for ( utility::vector1< core::Real >::iterator it = resi_vec[ res_i ][ rot_i ]->chi().begin(); it != resi_vec[ res_i ][ rot_i ]->chi().end(); ++it ) {
				seq_stringstream<<*it<<" ";
			}
			seq_stringstream<<",";
		}
		seq_stringstream<<std::endl;;
	}
	db_file << seq_stringstream.str()<<std::endl;
	db_file.close();

	std::ofstream seq_aln_file;
	seq_stringstream.str(std::string());
	seq_aln_file.open(seq_aln_file_.c_str(), std::ios::app);
	for ( std::map< std::string, utility::vector1< std::string > >::iterator it = AA_vec.begin(); it != AA_vec.end(); ++it ) {
		seq_stringstream<<">" << it->first << std::endl;
		for ( core::Size res = 1; res <= it->second.size() ; ++res ) {
			seq_stringstream<<it->second[res];
		}
		seq_stringstream<<std::endl;
		//TR<<seq_stringstream.str()<<std::endl;
	}
	seq_aln_file << seq_stringstream.str()<<std::endl;
	seq_aln_file.close();
}


void //print to screen rotamer vector
RotLibOut::print(){

}
void
RotLibOut::apply( Pose & pose )
{
	using namespace basic::options;
	TR<<"min_dist is set to:" <<min_dist_<<std::endl;
	resi_vec.resize(pose.total_residue());

	pdb_name_list = option[ OptionKeys::packing::unboundrot ]() ;
	for ( core::Size i = 1; i <= pdb_name_list.size(); ++i ) {
		TR<< "working on PDB " << pdb_name_list[ i ] << std::endl;
		// std::string base_filename = pdb_name_list[ i ].substr(pdb_name_list[ i ].find_last_of("/\\") + 1);
		core::pose::PoseOP hit_pose( new core::pose::Pose() );
		core::io::pdb::build_pose_from_pdb_as_is( *hit_pose, pdb_name_list[ i ] );
		std::string base_filename = hit_pose->pdb_info()->name().substr(hit_pose->pdb_info()->name().find_last_of("/\\") + 1);
		AA_vec[base_filename].resize(pose.total_residue());
		//TR<<base_filename<<std::endl;
		std::fill(AA_vec[base_filename].begin(), AA_vec[base_filename].end(), "-");//initialize arrays to "-" (gaps)
		this->find_matching_res( pose, *hit_pose);
	}
	this->print_lib_to_file();
}

std::string
RotLibOut::get_name() const {
	return RotLibOutCreator::mover_name();
}

std::string
RotLibOut::mover_name()  {
	return "RotLibOut";
}

moves::MoverOP
RotLibOut::clone() const
{
	return moves::MoverOP( new RotLibOut( *this ) );
}

moves::MoverOP
RotLibOut::fresh_instance() const
{
	return moves::MoverOP( new RotLibOut );
}

void
RotLibOut::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	set_min_dist(tag->getOption< core::Real >( "min_dist", 2.0 ));
	Rotamer_db_file(tag->getOption< std::string >( "rotamer_db_filename", "Rotamer_db.dat" ));
	seq_aln_file(tag->getOption< std::string >( "sequence_alignment_filename", "seqeunce_aln.fa" ));
	dump_pdb(tag->getOption< bool >( "dump_pdb", "0" ));// if set to true dump aligned fragments between template and target
	min_frag_length(tag->getOption< core::Size >( "min_frag_len", 3 ) );// if set to true dump aligned fragments between template and target
}
void RotLibOut::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "min_frag_len", xs_string, "Set minimum fragment length to use in alignmnet" )
		+ XMLSchemaAttribute( "jump_dbase_fname", xs_string, "filename of jump db file" )
		+ XMLSchemaAttribute::attribute_w_default( "jump_from_foldtree", xsct_rosetta_bool, "boolean- whether or not use jump from foldtree", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "min_dist", xsct_non_negative_integer, "Minimum distance between aligned residues", "3" )
		+ XMLSchemaAttribute::attribute_w_default( "dump_pdb", xsct_rosetta_bool, "for debuging, whether to dump pdbs during run", "false" )
		+ XMLSchemaAttribute( "rotamer_db_filename", xs_string, "Path to save rotamer db file")
		+ XMLSchemaAttribute( "sequence_alignment_filename", xs_string, "Path to save rotamer db file");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}


void RotLibOutCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RotLibOut::provide_xml_schema( xsd );
}
} // splice
} // protocols

