// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file 
/// @author Eva-Maria Strauch (evas01@u.washington.edu)
/// @brief based on Bruno's CA cst generator
// Unit headers
#include <protocols/seeded_abinitio/CAcstGenerator.hh>
#include <protocols/seeded_abinitio/CAcstGeneratorCreator.hh>

//#include <protocols/protein_interface_design/util.hh>
#include <utility/string_util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <boost/foreach.hpp>

#include <core/types.hh>

#define foreach BOOST_FOREACH

#include <basic/Tracer.hh>
#include <protocols/moves/DataMap.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

//Auto Headers
#include <core/pose/util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/import_pose.hh>


//protocols
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/loops_main.hh>

//constraints
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/Func.hh>

//options
#include <basic/options/option.hh>
#include <basic/options/keys/fold_from_loops.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


namespace protocols {
namespace seeded_abinitio {

using namespace core;
using namespace scoring::constraints;	
using namespace protocols::moves;
	
static basic::Tracer TR( "protocols.seeded_abinitio.CAcstGenerator" );
//static basic::Tracer TR_debug( "DEBUG.CAcstGenerator" );
	
	std::string
	CAcstGeneratorCreator::keyname() const
	{
		return CAcstGeneratorCreator::mover_name();
	}
	
	protocols::moves::MoverOP
	CAcstGeneratorCreator::create_mover() const {
		return new CAcstGenerator;
	}
	
	std::string
	CAcstGeneratorCreator::mover_name()
	{
		return "CAcstGenerator";
	}
	
	CAcstGenerator::~CAcstGenerator() {}

	CAcstGenerator::CAcstGenerator() :
	protocols::moves::Mover( CAcstGeneratorCreator::mover_name() )
	{
		stddev_ = 3.0;
	}
	
	
bool 
is_cut(	utility::vector1<core::Size> & cut_points,
				core::Size & residue){
	
	bool res_cut = false;
	
	for ( core::Size it = cut_points[1]; it <= cut_points[ cut_points.size() ]; it++ ){
		if (it == residue)
			res_cut = true; 
		}
	return res_cut;
}


///this is still not bulletproof the numbering of the actual seeds in case it is either the template used or the input	
void add_dist_constraints( 					
						 				pose::Pose & pose,
									 	pose::PoseOP & pose_of_int,
										core::Size start_relevant_chain,
										core::scoring::constraints::ConstraintSetOP & cst,
										protocols::loops::Loops & seeds,// referring to template if one is given! or just individual chain numbering
										protocols::loops::Loops & clear_area,
									 	utility::vector1< core::Size > cut_points,
										bool add_cst_seed,
					 					core::Real stddev,
										core::Size seq_separation
										){
	
	using namespace scoring::constraints;
	using namespace id;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	TR<<"stddev for harmonic constraints: " << stddev <<std::endl;
	
	if( cut_points.empty() )
		TR<<"there are no cut point registered" <<std::endl;		
	
	TR.Debug << "start_relevant_chain: " << start_relevant_chain << std::endl;
	TR.Debug << "seeds: " << seeds << std::endl;
	TR.Debug << "clear_area " << clear_area << std::endl;
	
	for( Size i=1;i<=cut_points.size(); i++) {
		TR <<"cutpoints: " << cut_points[i] << std::endl;
	}
	
	//adjust cutpoints first to relevant numbering
	for( Size i = 1 ; i <= cut_points.size(); ++i ){
		TR.Debug <<"old cutpoint: "<< cut_points[i] <<std::endl;
		cut_points[i] = cut_points[i] - (start_relevant_chain -1);
		TR.Debug<<"adjusted cutpoint: "<< cut_points[i] <<std::endl;
	}
	
	TR.Debug  << "adjustied cutpoints:" << std::endl;

  for( Size i=1;i<=cut_points.size(); i++) {
    TR <<"cutpoints: " << cut_points[i] << std::endl;
  }


	for ( Size pos = 1; pos <= pose_of_int->total_residue(); pos++ ) {
			for (Size pos_2 = 1; pos_2 <=pose_of_int->total_residue(); pos_2++   ){
				
				//if residues are part of the loop, set to true					
				bool res_is_loop = seeds.is_loop_residue( pos );
				bool res2_is_loop = seeds.is_loop_residue( pos_2 );
				
				//std::cout<<"pos "<< pos <<" "<< res_is_loop << ", pos2: "<<" "<<pos_2 <<" " << res2_is_loop <<std::endl;
				
				bool cut_point_pos = false;
				bool cut_point_pos_2 = false;
				bool res_cst_free = false;
		
				//add the cst of the segment that will be replaced by the seed(s)
				if ( add_cst_seed ){
					res_is_loop = false;
					res2_is_loop = false;
				}
				
				//user specified area that is allowed to float
				if ( clear_area.size() > 0 ){ 
					if ( clear_area.is_loop_residue(pos) || clear_area.is_loop_residue(pos_2) ){
						res_cst_free = true;
					}
				}
				
				//mark cut points
				if ( cut_points.size() != 0 ) {
					cut_point_pos = is_cut(cut_points, pos );
					cut_point_pos_2 = is_cut(cut_points, pos_2 );
				}
				
				//avoiding doubling of constraints
				Size seq_sep = 0;
				
				if (pos > pos_2)
					seq_sep = pos - pos_2;
				else 
					seq_sep = 0; //to avoid repetition
					//seq_sep = pos_2 - pos;
				
				if ( seq_sep >= seq_separation  ){
					if ( !res_is_loop && !res2_is_loop ) {
						if( !cut_point_pos && !cut_point_pos_2 ){
							//if ( !res_is_loop_neighbor && !res2_is_loop_neighbor) {
							
							if (!res_cst_free){
		 						core::conformation::Residue res_pos = pose_of_int->residue(pos);
								core::conformation::Residue res_pos_2 = pose_of_int->residue(pos_2);
								
								Real const distance_ca( res_pos.xyz( res_pos.atom_index("CA") ).distance( res_pos_2.xyz( res_pos_2.atom_index("CA") )));
								TR.Debug <<"template: " << pos_2 <<" "<< pos << " "<<distance_ca <<std::endl;
								TR.Debug <<"updated: "<<pos_2 + start_relevant_chain -1 <<" " << pos + start_relevant_chain - 1<< " " <<distance_ca << std::endl;
								
								//need to adjust numbering to the current input pose!
								core::conformation::Residue res_in_pose = pose.residue( pos + start_relevant_chain - 1 );
								core::conformation::Residue res_in_pose2= pose.residue( pos_2 + start_relevant_chain - 1);
								
								cst->add_constraint( new AtomPairConstraint ( AtomID(res_in_pose.atom_index("CA"), pos + start_relevant_chain - 1), AtomID(res_in_pose2.atom_index("CA"),pos_2 + start_relevant_chain - 1), new HarmonicFunc( distance_ca, stddev ) ) );
								}
							//}
						}
					}
				}//end seq_sep
			}
		}
}

	
void
CAcstGenerator::apply( pose::Pose & pose ){

	using namespace core;

	utility::vector1< Size > cutpoints = pose.fold_tree().cutpoints();
  
	for( Size i=1;i<=cutpoints.size(); i++) {
    TR.Debug<<"cutpoints: " << cutpoints[i] << std::endl;
  }

	TR.Debug << "foldtree: " << pose.fold_tree();

	pose::PoseOP donor_poseOP;
	Size start_recipient_chain;
	ca_cst_ =  new core::scoring::constraints::ConstraintSet();
	
	//this all has to be checked during run time since parse time input will be inaccurate if the pose has been grown before
	if( from_chain_ == 0 ){
		TR<<"user did not specify for which chain constraints should be derrived, defaulting to all chains"<<std::endl;
		if( template_presence_ ){
			donor_poseOP = *template_pdb_ ;
			TR<<"derriving CA distance constraints from the template pdb, since a template was given"<<std::endl;
		}
		else 
			donor_poseOP = pose;
	}
	
	//statements to prevent bogus from happening
	if( to_chain_ == 0 && template_presence_ ){
		if (pose.total_residue() != template_pdb_->total_residue() )
			utility_exit_with_message("chain(s) to derrive constraints form and chain(s) to apply it to, do NOT have the same residue number. NOT supported." );
	}
	
	if( from_chain_ != 0 ){
		if( template_presence_ ){
			donor_poseOP = new pose::Pose( template_pdb_->split_by_chain( from_chain_ ) );
			TR<<"derriving CA distance constraints from the template pdb"<<std::endl;
		}
		else{
			donor_poseOP = new pose::Pose( pose.split_by_chain( from_chain_ ) );
		}
	}
	
	if( to_chain_ != 0 ){
		start_recipient_chain = pose.conformation().chain_begin( to_chain_ );
		TR<<"adding constraints starting from index residue: " <<start_recipient_chain <<std::endl;	
	}
											
	if( from_chain_ != 0 && to_chain_ != 0 ){
		TR<<"donor chain length: " << donor_poseOP->total_residue() <<", recipient chain length: " << pose.conformation().chain_end( to_chain_ ) - pose.conformation().chain_begin( to_chain_) +1 <<std::endl;
		if( donor_poseOP->total_residue() != (pose.conformation().chain_end( to_chain_ ) - pose.conformation().chain_begin( to_chain_) + 1) )
			utility_exit_with_message("donor pose and recipient chain do not have the same residue numbers, check your template or input pdb");
	}

	
	add_dist_constraints( pose, donor_poseOP , start_recipient_chain , ca_cst_, all_seeds_ , clear_seeds_ , cutpoints, add_cst_seed_ , stddev_, seq_separation_ );
	
	if( replace_ ){
		TR<<"replacing all constraints with newly generated constraint set" <<std::endl;
		pose.constraint_set( ca_cst_ );
	}
	else {
		utility_exit_with_message("ADDing new constraints to pose, is currently not supported, try just replacing"); //<<std::endl;
	}
	TR.flush();
	//pose.dump_pdb("constraints.pdb");

}
	
std::string
CAcstGenerator::get_name() const {
		return CAcstGeneratorCreator::mover_name();
		}
	
	
void
CAcstGenerator::parse_my_tag( TagPtr const tag,
							  DataMap & /*data*/,
							  protocols::filters::Filters_map const & /*filters*/,
							  Movers_map const &,
							  Pose const & pose){
	
	TR<<"CAcstGenerator has been invoked"<<std::endl;
	
	/// the constraints can be either derrived from a chain of the input pose or from a template pose
	/// If a template pose is provided, the constraints will be derrived from that pdb, defaulting to chain 1 of it
	/// otherwise, the user can specify for which chain of the input pose constraints will be derrived
	/// alternatively, the user can specify for which chain of the input pdb constraints shoudl be derrived
	/// currently it only allows one chain at a time or all chains....
	
	stddev_ = tag->getOption<core::Real>( "stddev", 3.0);/// this needs to be in an constructor too
	
	TR<<"setting constraint standard deviation to "<< stddev_<< std::endl;
	
	if( tag->hasOption( "template_pdb" ) ){
		std::string const template_pdb_fname( tag->getOption< std::string >( "template_pdb" ));
		template_pdb_ =  new core::pose::Pose ;
		core::import_pose::pose_from_pdb( *template_pdb_, template_pdb_fname );
		TR<<"read in a template pdb with " <<template_pdb_->total_residue() <<"residues"<<std::endl;
		template_presence_ = true;
	}

	add_cst_seed_ = tag->getOption< bool >("add_cst_seed", 0 ); ///header
	
	replace_ = tag->getOption< bool >("replace", 1 );
	
	seq_separation_ = tag->getOption< core::Size >( "seq_separation", 6 );

	//parsing branch tags
	utility::vector0< TagPtr > const branch_tags( tag->getTags() );
	
	foreach( TagPtr const btag, branch_tags ){
		//parse the pdb of interest, which is either the template or the input pdb depending on the users specificiation
		if( template_presence_ )
			curr_pose_ = template_pdb_;
			
		if( btag->getName() == "Seeds" ) { //need an assertion for the presence of these or at least for the option file
			//needs some assertions to avoid bogus input
			std::string const beginS( btag->getOption<std::string>( "begin" ) );
			std::string const endS( btag->getOption<std::string>( "end" ) );
			core::Size const begin( protocols::rosetta_scripts::parse_resnum( beginS, *curr_pose_ ) );
			core::Size const end( protocols::rosetta_scripts::parse_resnum( endS, *curr_pose_ ) );
			
			TR.Debug <<"parsing seeds: \n"<< begin <<" and " << end <<std::endl; 
			TR.Debug <<"seeds: "<< all_seeds_ <<std::endl;
			
			all_seeds_.add_loop( begin , end , 0, 0, false );
			
		}//end seed tags
		
		if( btag->getName() == "Clear_cst_segment" ) { //need an assertion for the presence of these or at least for the option file
			
			std::string const begin_str( btag->getOption<std::string>( "begin" ) );
			std::string const end_str( btag->getOption<std::string>( "end" ) );
			core::Size const begin( protocols::rosetta_scripts::parse_resnum( begin_str, *curr_pose_ ) );
			core::Size const end( protocols::rosetta_scripts::parse_resnum( end_str, *curr_pose_ ) );
			clear_seeds_.add_loop( begin , end , 0, 0, false );
			
		}//end seed tags
	}//end branch tags
	
	///could be eventually a vector of chains if desired
	///but currently just the simple version...
	
	if( tag->hasOption( "from_chain" ) ){
		from_chain_ = tag->getOption< core::Size >( "from_chain", 1 );
		TR<<"chain to derrive constraints from: "<< from_chain_ <<std::endl; 
	}
	
	if( tag->hasOption( "to_chain" ) ){
		to_chain_ = tag->getOption< core::Size >( "to_chain", 1 );
		TR<<"chain to apply constraints to: "<< to_chain_ <<std::endl; 
	}
								 
	}//end parse my tag							 
							 
	}//CAcstGenerator
}//protocol
