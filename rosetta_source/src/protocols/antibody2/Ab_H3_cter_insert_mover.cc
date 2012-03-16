// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/Ab_H3_cter_insert_mover.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#include <protocols/antibody2/Ab_H3_cter_insert_mover.hh>

#include <protocols/loops/loop_closure/ccd/CcdLoopClosureMover.hh>
#include <protocols/loops/loops_main.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/antibody2/Ab_Info.hh>


//#include <protocols/antibody2/CDRH3Modeler2.hh>
//TODO: JQX: the reason I used this CDRH3Modeler2 head file 
// is because of the "simple_one_loop_fold_tree" function in CDRH3Modeler2 file.
// Planning to make a independent file to save this function, since it is widely used
#include <protocols/antibody2/Ab_util.hh>



static numeric::random::RandomGenerator RG(19810606);

static basic::Tracer TR("protocols.antibody2.Ab_H3_cter_insert_mover");
using namespace core;



namespace protocols {
namespace antibody2 {
    
    
    
    
        
Ab_H3_cter_insert_mover::Ab_H3_cter_insert_mover(antibody2::Ab_Info & ab_info) : Mover()   
{
    init(ab_info, false, false);     
}
        
        
        
Ab_H3_cter_insert_mover::Ab_H3_cter_insert_mover(antibody2::Ab_Info & ab_info, bool camelid ) : Mover()
{
    user_defined_ = true;
            
    init(ab_info, camelid, false);
}
        
        
        
        
// Ab_H3_cter_insert_mover default destructor
Ab_H3_cter_insert_mover::~Ab_H3_cter_insert_mover() {}
        
        
        
        
        
        
        
void Ab_H3_cter_insert_mover::init(antibody2::Ab_Info & ab_info, bool camelid, bool benchmark)
{
    Mover::type( "Ab_H3_cter_insert_mover" );
            
    set_default();
            
    if ( user_defined_ ) {
        is_camelid_ = camelid;
        benchmark_  = benchmark;
    }
            
//    setup_objects();
    
    
    read_H3_cter_fragment(ab_info, camelid ) ;
    
    
    
}
        
    
    
    
        
void Ab_H3_cter_insert_mover::set_default(){
    is_camelid_ = false;
    benchmark_ = false;
}
        
        
//void Ab_H3_cter_insert_mover::setup_objects(){
//}
        
        
        
//void Ab_H3_cter_insert_mover::finalize_setup(core::pose::Pose & pose ){
        
//}
        
        
        
void Ab_H3_cter_insert_mover::apply(pose::Pose & pose)    
{        
//    finalize_setup(pose);
            
            
        
    antibody_modeling_insert_ter(pose) ;
            
}
        
        
        
        
std::string Ab_H3_cter_insert_mover::get_name() const {
    return "Ab_H3_cter_insert_mover";
}
        
        
        
        
        
        
        
        
void Ab_H3_cter_insert_mover::read_H3_cter_fragment(
                                antibody2::Ab_Info & ab_info,
                                bool is_camelid ) 
{
    using namespace fragment;
            
    std::string const path = basic::options::option[ basic::options::OptionKeys::in::path::path ]()[1];
            
    TR <<  "H3M Reading CDR H3 C-ter Fragments" << std::endl;
            
    bool is_kinked( ab_info.is_kinked() );
    bool is_extended( ab_info.is_extended() );
            
    // extract single letter aa codes for the chopped loop residues
    Size cdr_h3_size = ( ab_info.get_CDR_loop("h3")->stop()-1 -
                         ab_info.get_CDR_loop("h3")->start() ) + 1;
    utility::vector1< char > aa_1name;
    for( Size ii = ab_info.get_CDR_loop("h3")->start() - 2;
              ii <= ( ab_info.get_CDR_loop("h3")->start() - 2 ) + cdr_h3_size + 3; ++ii )
        aa_1name.push_back( ab_info.get_Fv_sequence()[ii] );
            
    // used only when no length & kink match are found
    utility::vector1< FragData > H3_base_library_seq_kink;
            
    // used only when no (length & kink) or (length & seq) are found
    utility::vector1< FragData > H3_base_library_kink;
            

    // file is read in from where other contraints are supposed to exist
    if( is_camelid )
        H3_ter_library_filename_ = path+"camelid_H3_CTERM";
    else
        H3_ter_library_filename_ = path+"H3_CTERM";
            
    // Read the file defined by command line option
    utility::io::izstream H3_ter_library_stream( H3_ter_library_filename_ );
            
    // Check to see if file exists
    if ( !H3_ter_library_stream ) {
        TR << "[Error]: Could not open H3 base library file: "
           << H3_ter_library_filename_ << std::endl
           << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
        std::exit( EXIT_FAILURE );
    }
            
    std::string pdb_name;
    std::string res_no;
    char res_name;
    Real phi(0.0);
    Real psi(0.0);
    Real omega(0.0);
    Size H3_length(0);
    Real resolution(0.0);
    std::string base_type;
            
    Size pdb_H3_length = cdr_h3_size;
    Size h3_base_frag_size( is_camelid ? 6 : 4 ); // changed from 4:6
    bool end_not_reached(true);
    while(end_not_reached){
        bool seq_match( true );
        bool kink_match( false );
                
        FragData f;
        f.set_valid( true );
                
        for ( Size i = 1; i <= h3_base_frag_size; ++i ) {
            H3_ter_library_stream >> pdb_name
                                  >> res_no
                                  >> res_name
                                  >> omega
                                  >> phi
                                  >> psi
                                  >> H3_length
                                  >> resolution
                                  >> base_type
                                  >> std::skipws;
            
            if(H3_ter_library_stream.eof()) {
                end_not_reached = false;
                break;
            }
            if( res_name != aa_1name[aa_1name.size() - 5 + i] ){
                seq_match = false;
            }
                    
            utility::pointer::owning_ptr< BBTorsionSRFD > res_torsions(
                        new BBTorsionSRFD( 3, 'L', res_name ) ); // 3 protein torsions
            
            // ugly numbers 1-3, but pose.set_phi also uses explicit numbers
            res_torsions->set_torsion   ( 1, phi   );
            res_torsions->set_torsion   ( 2, psi   );
            res_torsions->set_torsion   ( 3, omega );
            res_torsions->set_secstruct ( 'L' );
            f.add_residue( res_torsions );
        }
        if( !is_camelid ) {
            if( is_kinked && base_type == "KINK" )
                kink_match = true;
            else if( is_extended && base_type == "EXTENDED" )
                kink_match = true;
            else if( !is_kinked && !is_extended && base_type == "NEUTRAL" )
                kink_match = true;
        }
        else {
            if( is_extended && base_type == "EXTENDED" )
                kink_match = true;
            else if( ( is_kinked && base_type == "KINK" ) ||
                    ( is_kinked && base_type == "EXTENDED" ) )
                kink_match = true;
            else if( !is_kinked && !is_extended )
                kink_match = true;
        }
        if( is_camelid && end_not_reached && kink_match ){
            H3_base_library_.push_back( f );
        }
        else if( end_not_reached && (H3_length==pdb_H3_length) && kink_match ){
            H3_base_library_.push_back( f );
        }
        if( end_not_reached && seq_match && kink_match ){
            H3_base_library_seq_kink.push_back( f );
        }
        if( end_not_reached && kink_match  ){
            H3_base_library_kink.push_back( f );
        }
    }
            
    H3_ter_library_stream.close();
    H3_ter_library_stream.clear();
            
    // if no match found based on sequence and kink match criterion
    // then choose based on size and kink match criterion
    // if still no match, then choose based only on kink
    if( H3_base_library_.size() == 0 ) {
        H3_base_library_ = H3_base_library_seq_kink;
    }
    if( H3_base_library_.size() == 0 ) {
        H3_base_library_ = H3_base_library_kink;
    }
            
    TR <<  "H3M Finished reading CDR H3 C-ter Fragments" << std::endl;
            
    return;
}
        
    
    
    
    
    
    
    
    
    
    
void Ab_H3_cter_insert_mover::antibody_modeling_insert_ter( core::pose::Pose & pose) 
{
            
    TR <<  "H3M Inserting CDR H3 C-ter Fragments" << std::endl;
            
    // Storing initial fold tree
    kinematics::FoldTree const input_tree( pose.fold_tree() );
            
    Size loop_begin(0), loop_end(0), cutpoint(0), random_H3_ter(0);
    //utility::vector1< fragment::FragData >::const_iterator H3_ter;
    fragment::FragData f;
            
    loop_begin = ab_info_->get_CDR_loop("h3")->start();
    cutpoint   = ab_info_->get_CDR_loop("h3")->start() + 1;
    random_H3_ter = RG.random_range( 1, H3_base_library_.size() );
    //H3_ter = H3_base_library_.begin();
            
    loop_end = ab_info_->get_CDR_loop("h3")->stop();
            
    loops::Loop cdr_h3( loop_begin, loop_end, cutpoint,	0, true );
    simple_one_loop_fold_tree( pose, cdr_h3 );
            
    // choosing a base randomly
    //H3_ter = H3_ter + random_H3_ter;
    f = H3_base_library_[ random_H3_ter ];
            
    //			pose::Pose start_pose = ab_info_.Fv;
    //inserting base dihedrals
    Size cter_insertion_pos( is_camelid_ ? 4 : 2 );
    if( (ab_info_->get_CDR_loop("h3")->stop()-1 - cter_insertion_pos) <=
        ab_info_->get_CDR_loop("h3")->start() )
        TR << "H3 LOOP IS TOO SHORT: CAN NOT USE N-TERM INFORMATION"
            << std::endl;
    else {
        // H3_ter->apply(...);
        f.apply( pose, ab_info_->get_CDR_loop("h3")->stop()-1 -
                cter_insertion_pos, ab_info_->get_CDR_loop("h3")->stop() );
    }
            
    // Restoring pose fold tree
    pose.fold_tree( input_tree );
            
    TR <<  "H3M Finished Inserting CDR H3 C-ter Fragments" << std::endl;
            
    return;
} 
        
        
        
        
        
}//antibody2
}//protocols





