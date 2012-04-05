// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody2/Ab_GraftOneCDR_Mover.cc
/// @brief grafts a cdr onto the template of an antibody framework
/// @detailed
/// @author Jianqing Xu (xubest@gmail.com)

#include <protocols/antibody2/Ab_GraftOneCDR_Mover.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/antibody2/Ab_TemplateInfo.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <protocols/idealize/IdealizeMover.hh>


#include <core/id/AtomID_Map.hh>
#include <core/pose/util.tmpl.hh>

#include <basic/Tracer.hh>


static basic::Tracer TRG("protocols.antibody2.Ab_GraftOneCDR_Mover");

namespace protocols {
namespace antibody2 {
using namespace core;



Ab_GraftOneCDR_Mover::Ab_GraftOneCDR_Mover(){}

Ab_GraftOneCDR_Mover::Ab_GraftOneCDR_Mover(std::string cdr_name, 
                                           Size query_start, 
                                           Size query_end, 
                                           scoring::ScoreFunctionOP scorefxn ) : Mover( "Ab_GraftOneCDR_Mover" )
{
    scorefxn_ = scorefxn;
	query_start_ = query_start;
	query_end_   = query_end;
	set_default( cdr_name );
    std::string const path = basic::options::option[ basic::options::OptionKeys::in::path::path ]()[1];
    TRG << "Reading in template: " << path << cdr_name << ".pdb " << std::endl;
    import_pose::pose_from_pdb( template_pose_, path + cdr_name + ".pdb" );
    
} //default constructor

    
    

Ab_GraftOneCDR_Mover::Ab_GraftOneCDR_Mover( std::string cdr_name, 
                                           AntibodyInfoOP ab_info, 
                                           Ab_TemplateInfoOP ab_t_info, 
                                           scoring::ScoreFunctionOP scorefxn ) : Mover( "Ab_GraftOneCDR_Mover" )
{
    scorefxn_ = scorefxn;
    query_start_ = ab_info->get_CDR_loop(cdr_name)->start();
	query_end_   = ab_info->get_CDR_loop(cdr_name)->stop();
	set_default(cdr_name);
    template_pose_ = ab_t_info->get_one_template_pose(cdr_name) ;
    TRG<< "template_pose_ "<<std::endl;
    TRG<< template_pose_<<std::endl;
}
    
    
    
    

// Ab_GraftOneCDR_Mover default destructor
Ab_GraftOneCDR_Mover::~Ab_GraftOneCDR_Mover() {}


    
    

void Ab_GraftOneCDR_Mover::set_default( std::string template_name )
{

/*        
    //idealize the loop
    idealize::IdealizeMover idealizer;
    idealizer.fast( false );
    idealizer.apply( template_pose_ );
    template_pose_.dump_pdb("idealized.pdb");
  */  
    

	template_name_  = template_name;
    flank_size_ = 4 ;//^^^^^^^^^^  
    TRG<<"flank_size: "<<flank_size_<<std::endl;
    // JQX:  the default value of flank_size_ is equle to 4, meaning there are 4 residues
    //       on the C-ter and N-ter of the actual loop
    //       However, based on the old R2 antibody code, only 3 residues on each stem 
} 

    
    
    
    
    




void Ab_GraftOneCDR_Mover::apply( pose::Pose & pose_in )
{

    TRG<<"flank_size: "<<flank_size_<<std::endl;

    TRG<<"Start to Graft CDRs   "<< template_name_  <<" ............"<<std::endl;
    Size const nres( pose_in.total_residue() ); // Total residues
    Size query_size = ( query_end_ - query_start_ )+1;




	// create a sub pose with  4 flanking residues on either side of CDR loop
//        TRG<<"query_start="<<query_start_<<std::endl;
//        TRG<<"query_end="<<query_end_<<std::endl;
//        TRG<<"flank_size="<<flank_size_<<std::endl;
//        TRG<<"query_size="<<query_size<<std::endl;
//        TRG<<"truncated_pose will be from "<<query_start_-flank_size_<<" to "<<query_end_+flank_size_<<std::endl;
    pose::Pose truncated_pose( pose_in, query_start_-flank_size_, query_end_+flank_size_ );

    // Just want to make life a litter easier
    //TRG<<"**************&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*************"<<std::endl;
    //TRG<<pose_in.sequence()<<std::endl;
    //TRG<<"  the template sequence:  "<<template_pose_.sequence()<<std::endl;
    //TRG<<"trucated_query sequence:  "<<truncated_pose.sequence()<<std::endl;
    //TRG<<"**************&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*************"<<std::endl;
    
    
	// create atom map for superimposing 2 flanking resiudes
    id::AtomID_Map< id::AtomID > atom_map;
    pose::initialize_atomid_map( atom_map, template_pose_, id::BOGUS_ATOM_ID );


    
    //   ****AAAAAAAAAAAAAAAAAAA****  the template pose should have 4 residues each side
    //    @@@AAAAAAAAAAAAAAAAAAA@@@   the real alignment only based on 3 residues each side
    
    for( Size start_stem = 2; start_stem <= flank_size_; ++start_stem ) { 
    /// starting from the 2nd residue in the l1-3.pdb, H1-3.pdb
        Size const ref_stem ( start_stem  );
        for( Size j=1; j <= 4; j++ ) {    /// four backbone heavy atoms
//          TRG<<"j="<<j<<"  start_stem_in_template_pose_="<<start_stem<<"  ref_stem_in_truncated_pose_="<<ref_stem<<std::endl;
            id::AtomID const id1( j, start_stem );
            id::AtomID const id2( j, ref_stem );
            atom_map[ id1 ] = id2;
        }
    }

	// start at the end of the actual loop
    for( Size end_stem = query_size+flank_size_+1; end_stem <= query_size+flank_size_+flank_size_-1; ++end_stem ) { 
        Size const ref_stem ( end_stem);  
//        if(template_name_ == "h3") Size const ref_stem(end_stem+1);
        for( Size j=1; j <= 4; j++ ) {    /// four backbone heavy atoms
//            TRG<<"j="<<j<<"  end_stem_in_template_pose_="<<end_stem<<"  ref_stem_in_truncated_pose_="<<ref_stem<<std::endl;
            id::AtomID const id1( j, end_stem );
            id::AtomID const id2( j, ref_stem );
            atom_map[ id1 ] = id2;
        }
    }


    scoring::superimpose_pose( template_pose_, truncated_pose, atom_map );
    template_pose_.dump_pdb(template_name_);

    // TODO:
    // JQX:
    // when copying the backbone heavy atom coordinates, you copy two more residues on each side
    // pay attention to this, it is "2", not "3" here. "3" was used to superimpose like the code did above. 
    // No reason for that, just to match R2_antibody. 
    // This may be the reason why you need hfr to reasign the phi, psi, or omega angle of h3_start-1
    // when doing H3 loop modeling. But this is still weird, because you need to take care
    // of 2 residues on each side then, not just one more on h3_start-1 or +1 position using hfr.pdb. 
    // Maybe I should remove this number 2, just copy the loop itself !!!!!
    
    for ( Size i=query_start_-2; i <= query_end_+2; ++i ) {
        Size template_num ( i - (query_start_-5) );  // this "5" is the default in the l1.pdb, you have 4 residues before the 1st one
//        TRG<<" i="<<i<<"    template_num="<<template_num<<std::endl;
        conformation::Residue const & source_rsd( pose_in.residue( i ) );
        conformation::Residue const & target_rsd( template_pose_.residue( template_num) );
//        TRG<<source_rsd<<std::endl;
//        TRG<<target_rsd<<std::endl;
        
        Size const natoms( source_rsd.natoms() );

        bool any_missing( false );
        id::AtomID_Mask missing( false );
        // dimension the missing-atom mask
        pose::initialize_atomid_map( missing, pose_in );

        if( source_rsd.name() != target_rsd.name() ) pose_in.set_secstruct( i, 'X' );
        for ( Size j=1; j<= natoms; ++j ) {
            std::string const & atom_name( source_rsd.atom_name(j) );
            if ( target_rsd.has( atom_name ) ) {
                pose_in.set_xyz( id::AtomID( source_rsd.atom_index( atom_name),i ), target_rsd.xyz( atom_name ) );
                //make sure the source_rsd in pose_in only has 4 backbone atoms
                //my script should take care of it
            } 
            else {
                TRG<<"watch out !!! missing "<<atom_name << "   !!!!!" <<std::endl;
                any_missing = true;
                missing[ id::AtomID( pose_in.residue_type(i).atom_index(source_rsd.atom_name(j)), i ) ] = true;
            }
        }

        if ( any_missing ) {
            pose_in.conformation().fill_missing_atoms( missing );
        }
    }

    pose_in.dump_pdb(template_name_+"_graft");

    
} // Ab_GraftOneCDR_Mover::apply


    
    
    
    
    
    
    
    
    
    
    
    
    
std::string Ab_GraftOneCDR_Mover::get_name() const { return "Ab_GraftOneCDR_Mover"; }
    
// copy ctor
Ab_GraftOneCDR_Mover::Ab_GraftOneCDR_Mover( Ab_GraftOneCDR_Mover const & rhs ) {
    initForEqualOperatorAndCopyConstructor(*this, rhs);
}
    
///@brief assignment operator
Ab_GraftOneCDR_Mover & Ab_GraftOneCDR_Mover::operator=( Ab_GraftOneCDR_Mover const & rhs ){
    //abort self-assignment
    if (this == &rhs) return *this;
    Mover::operator=(rhs);
    initForEqualOperatorAndCopyConstructor(*this, rhs);
    return *this;
}
    
void Ab_GraftOneCDR_Mover::initForEqualOperatorAndCopyConstructor(Ab_GraftOneCDR_Mover & lhs, Ab_GraftOneCDR_Mover const & rhs) {
        
}
    
    
    
    

}  // namespace antibody2
}  // namespace protocols
