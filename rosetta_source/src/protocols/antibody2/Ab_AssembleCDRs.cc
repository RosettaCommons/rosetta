// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @author Jianqing Xu ( xubest@gmail.com )

#include <protocols/jobdist/JobDistributors.hh> // SJF Keep first for mpi

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueSelector.hh>

#include <core/chemical/VariantType.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
//#include <basic/options/keys/antibody2.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/prof.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include <basic/Tracer.hh>


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
using namespace ObjexxFCL::fmt;

#include <protocols/jd2/ScoreMap.hh>

#include <protocols/antibody2/Ab_Info.hh>
#include <protocols/antibody2/Ab_TemplateInfo.hh>

#include <protocols/antibody2/Ab_AssembleCDRs.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>


#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/moves/TrialMover.hh>



#include <protocols/antibody2/Ab_GraftOneCDR_Mover.hh>
#include <protocols/antibody2/Ab_CloseOneCDR_Mover.hh>
#include <protocols/antibody2/Ab_RelaxCDRs_Mover.hh>



//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody2.Ab_AssembleCDRs");
using namespace core;

namespace protocols {
namespace antibody2 {

// default constructor
Ab_AssembleCDRs::Ab_AssembleCDRs() : Mover() {
	user_defined_ = false;
	init();
}

// default destructor
Ab_AssembleCDRs::~Ab_AssembleCDRs() {}

//clone
protocols::moves::MoverOP
Ab_AssembleCDRs::clone() const {
	return( new Ab_AssembleCDRs() );
}

    
    
    
    
    
    
    
void Ab_AssembleCDRs::init() {
	Mover::type( "Ab_AssembleCDRs" );

	// setup all the booleans with default values
	// they will get overwritten by the options and/or passed values
    
	set_default();
	register_options();
	init_from_options();

    
//	if ( ab_model_score() == NULL ) { //<- use this if we want to pass in score functions
		// score functions
		scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "standard", "score12" );
		scorefxn_->set_weight( core::scoring::chainbreak, 1.0 );
		scorefxn_->set_weight( core::scoring::overlap_chainbreak, 10./3. );
		scorefxn_->set_weight( core::scoring::atom_pair_constraint, 1.00 );
//	}

	setup_objects();

}

    
    
    
    
    
    
void Ab_AssembleCDRs::set_default()
{
	TR <<  "Setting up default settings to all FALSE" << std::endl;
	graft_l1_  = false;
	graft_l2_  = false;
	graft_l3_  = false;
	graft_h1_  = false;
	graft_h2_  = false;
	graft_h3_  = false;
	benchmark_ = false;
	camelid_   = false;
	camelid_constraints_ = false;

}

    
    
    
    
    
    
void Ab_AssembleCDRs::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::antibody::camelid );
	option.add_relevant( OptionKeys::antibody::camelid_constraints );
	option.add_relevant( OptionKeys::antibody::graft_l1 );
	option.add_relevant( OptionKeys::antibody::graft_l2 );
	option.add_relevant( OptionKeys::antibody::graft_l3 );
	option.add_relevant( OptionKeys::antibody::graft_h1 );
	option.add_relevant( OptionKeys::antibody::graft_h2 );
	option.add_relevant( OptionKeys::antibody::graft_h3 );
	option.add_relevant( OptionKeys::constraints::cst_weight );
	option.add_relevant( OptionKeys::run::benchmark );
	option.add_relevant( OptionKeys::in::file::native );
}

    
    
    
    
void Ab_AssembleCDRs::init_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	TR <<  "Reading Options" << std::endl;

	if ( option[ OptionKeys::antibody::graft_l1 ].user() )
                set_graft_l1( option[ OptionKeys::antibody::graft_l1 ]() );
	if ( option[ OptionKeys::antibody::graft_l2 ].user() )
                set_graft_l2( option[ OptionKeys::antibody::graft_l2 ]() );
	if ( option[ OptionKeys::antibody::graft_l3 ].user() )
                set_graft_l3( option[ OptionKeys::antibody::graft_l3 ]() );
	if ( option[ OptionKeys::antibody::graft_h1 ].user() )
                set_graft_h1( option[ OptionKeys::antibody::graft_h1 ]() );
	if ( option[ OptionKeys::antibody::graft_h2 ].user() )
                set_graft_h2( option[ OptionKeys::antibody::graft_h2 ]() );
	if ( option[ OptionKeys::antibody::graft_h3 ].user() )
                set_graft_h3( option[ OptionKeys::antibody::graft_h3 ]() );
	if ( option[ OptionKeys::antibody::camelid ].user() )
                set_camelid( option[ OptionKeys::antibody::camelid ]() );
	if ( option[ OptionKeys::antibody::camelid_constraints ].user() )
                set_camelid_constraints( option[ OptionKeys::antibody::camelid_constraints ]() );
	if ( option[ OptionKeys::run::benchmark ].user() )
                set_benchmark( option[ OptionKeys::run::benchmark ]() );
    
    if ( option[ OptionKeys::constraints::cst_weight ].user() )
        set_cst_weight( option[ OptionKeys::constraints::cst_weight ]() );


	//set native pose if asked for
	if ( option[ OptionKeys::in::file::native ].user() ) {
		core::pose::PoseOP native_pose = new core::pose::Pose();
		core::import_pose::pose_from_pdb( *native_pose, option[ OptionKeys::in::file::native ]() );
		set_native_pose( native_pose );
	}
	else{
		set_native_pose(NULL);
	}
    
    
	cst_weight_ = option[ OptionKeys::constraints::cst_weight ]();
    
	if( camelid_ ) {
		graft_l1_ = false;
		graft_l2_ = false;
		graft_l3_ = false;
	}
    
    grafts_.insert( std::pair< std::string, bool >("l1", graft_l1_) );
    grafts_.insert( std::pair< std::string, bool >("l2", graft_l2_) );
    grafts_.insert( std::pair< std::string, bool >("l3", graft_l3_) );
    grafts_.insert( std::pair< std::string, bool >("h1", graft_h1_) );
    grafts_.insert( std::pair< std::string, bool >("h2", graft_h2_) );
    grafts_.insert( std::pair< std::string, bool >("h3", graft_h3_) );
    


}

    
    
    
    
    
    
    
    
void Ab_AssembleCDRs::setup_objects() {
    ab_info_ = NULL;
    ab_t_info_ = NULL;
    
    graft_sequence_ = NULL;
    relax_sequence_= NULL;
    packer_ = NULL;
    pymol_=NULL;
    
    scorefxn_ = NULL;

	sync_objects_with_flags();
}
    
    
    
void Ab_AssembleCDRs::sync_objects_with_flags() {

	using namespace protocols::moves;


	flags_and_objects_are_in_sync_ = true;
	first_apply_with_current_setup_ = true;
}


    
    
    
    
    

    

void Ab_AssembleCDRs::finalize_setup( pose::Pose & frame_pose ) {
    TR<<" finalize_setup ............."<<std::endl;
	TR<<"AAAAAAAA     cst_weight: "<<cst_weight_<<std::endl;

	// check for native and input pose
	if ( !get_input_pose() ) {
		pose::PoseOP input_pose = new pose::Pose(frame_pose);  //JQX: QUESTION: why owning pointer here
		set_input_pose( input_pose );   // JQX: pass the input_pose to the mover.input_pose_
	}


	pose::PoseOP native_pose;
	if ( !get_native_pose() ) {
		TR << "Danger Will Robinson! Native is an impostor!" << std::endl;
        TR << "   'native_pose' is just a copy of the 'input_pose'    " << std::endl;
        TR << "    since you didn't sepcifiy the native pdb name"<<std::endl;
		native_pose = new pose::Pose(frame_pose);
	} else {
		native_pose = new pose::Pose( *get_native_pose() );
	}

    // JQX: this is the secondary structure from the native pose
	pose::set_ss_from_phipsi( *native_pose ); 

	set_native_pose( native_pose ); // pass the native pose to the mover.native_pose_

    
    ab_info_   =  new Ab_Info(frame_pose, camelid_);
    ab_t_info_ =  new Ab_TemplateInfo(graft_l1_, graft_l2_, graft_l3_,
                                      graft_h1_, graft_h2_, graft_h3_, camelid_);
    
    graft_sequence_ = new moves::SequenceMover();
    relax_sequence_ = new moves::SequenceMover();
    pymol_ = new moves::PyMolMover();
    
    
    TR<<" Checking Ab_Info object: "<<std::endl<<*ab_info_<<std::endl<<std::endl;
    TR<<" Checking Ab_TemplateInfo object: "<<std::endl<<*ab_t_info_<<std::endl<<std::endl;

    for ( GraftMap::const_iterator it = grafts_.begin(); it != grafts_.end(); ++it ) {
        if ( it->second ) {
            TR << "Creating movers for " << it->first << std::endl;
            TR << "                  start (chothia): "<<ab_info_->get_CDR_loop(it->first)->start()<<std::endl;
            TR << "                   stop (chothia): "<<ab_info_->get_CDR_loop(it->first)->stop()<<std::endl;
            
            Ab_GraftOneCDR_MoverOP graft_one_cdr = new Ab_GraftOneCDR_Mover( it->first, ab_info_, ab_t_info_, scorefxn_) ;
            graft_one_cdr->enable_benchmark_mode( benchmark_ );
            graft_sequence_->add_mover( graft_one_cdr);
            //              graft_sequence_->add_mover( pymol_ );
            
            
            /*
             Ab_CloseOneCDR_MoverOP closeone( new Ab_CloseOneCDR_Mover( ab_info.get_loop(it->first)->start(),
             ab_info.get_loop(it->first)->stop()   )     );
             closeone->enable_benchmark_mode( benchmark_ );
             closeone->set_pymol( pymol_ );
             graft_sequence_->add_mover( closeone );
             graft_sequence_->add_mover( pymol_ );
             
            
            
            Ab_RelaxCDRs_MoverOP rlx_one_cdr = new Ab_RelaxCDRs_Mover(  ab_info_->get_CDR_loop(it->first)->start(),
                                                                        ab_info_->get_CDR_loop(it->first)->stop()    );
            
            rlx_one_cdr->enable_benchmark_mode( benchmark_ );
            relax_sequence_->add_mover( rlx_one_cdr );
            relax_sequence_->add_mover( pymol_ );
             */
        }
    }
    
    // Exact match Aroop's old code in Rosetta 2:
    // graft all CDRs by superimpose stems, then pack the whole new pose
    
    
    set_packer_default(frame_pose, true /* include_current */)  ;
    graft_sequence_->add_mover(packer_);


}



    
    
//APPLY
void Ab_AssembleCDRs::apply( pose::Pose & frame_pose ) {

    using namespace chemical;
    using namespace id;
    using namespace scoring;
    using namespace core::scoring::constraints;
    using namespace protocols::moves;

//  I assume the pose is from the job distributor, which can take the -s flag to get the pose
    // the below test proves that the inital secstruct is all "L"!
/*    TR<<"JQX:    this is the 1st time that the 'pose' is used in the code: "<<std::endl;
    TR<<pose<<std::endl;
    for ( Size i = 1; i <= pose.total_residue(); ++i ) {
            TR<<"JQX:   residue: "<<i<<"       secstruct: "<<pose.secstruct(i)<<std::endl;
    }
    exit(-1);   */
    
    TR<<" in the apply function "<<std::endl;
    
    
    protocols::moves::PyMolMover pymol;
    if ( !flags_and_objects_are_in_sync_ ){ 
       sync_objects_with_flags(); 
    }
    
    if ( first_apply_with_current_setup_ ){ 
        finalize_setup(frame_pose);  
        first_apply_with_current_setup_=false; 
    }



	basic::prof_reset();
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	// utility::exit( EXIT_FAILURE, __FILE__, __LINE__);

	pose::set_ss_from_phipsi( frame_pose );
    
	pymol.apply( frame_pose );
	pymol.send_energy( frame_pose );

	// display constraints and return
	if( camelid_constraints_ ) {
		display_constraint_residues( frame_pose );
		return;
	}

    
    
    
    Size nres = frame_pose.total_residue();
    
    // Storing secondary structure
    utility::vector1<char> secondary_struct_storage;
    for( Size i = 1; i <= nres; i++ ) {
        secondary_struct_storage.push_back( frame_pose.secstruct( i ) );
        //        TR<<"JQX:   residue: "<<i<<"       secstruct: "<<frame_pose.secstruct(i)<<std::endl;
    }
    
    graft_sequence_->apply( frame_pose );
    
    frame_pose.dump_pdb("finish_grafting_and_packing.pdb");

    // Recover secondary structures
    for( Size i = 1; i <= nres; i++ ) frame_pose.set_secstruct( i, secondary_struct_storage[ i ] );
    
    // relax optimized CDR grafted regions
    relax_sequence_->apply( frame_pose );
    
    
    // Recover secondary structures
    for( Size i = 1; i <= nres; i++ ) frame_pose.set_secstruct( i, secondary_struct_storage[ i ] );
    
    // align pose to native pose
    pose::Pose native_pose;
    if( get_native_pose() ) native_pose = *get_native_pose();
    else                    native_pose = frame_pose;
    
    Ab_Info native_ab( native_pose, camelid_ );
    
    
    ab_info_->align_to_native( frame_pose, native_ab, native_pose );
    
    
    basic::prof_show();


}// end apply








std::string Ab_AssembleCDRs::get_name() const {
	return "Ab_AssembleCDRs";
}

void Ab_AssembleCDRs::set_packer_default(pose::Pose & pose, bool include_current) {
    //set up packer
    pack::task::PackerTaskOP task;
    task = pack::task::TaskFactory::create_packer_task( pose );
    task->restrict_to_repacking();
    task->or_include_current( include_current );
    packer_ = new simple_moves::PackRotamersMover( scorefxn_, task );
} // Ab_GraftCDRs_Mover set_packer_default





void Ab_AssembleCDRs::display_constraint_residues( core::pose::Pose & pose ) {		
    // Detecting di-sulfide bond

    Size H1_Cys(0), H3_Cys(0);

    if(      pose.residue( pose.pdb_info()->pdb2pose('H',32 ) ).name3() == "CYS" )
        H1_Cys = pose.pdb_info()->pdb2pose( 'H', 32 );
    else if( pose.residue( pose.pdb_info()->pdb2pose('H',33 ) ).name3() == "CYS" )
        H1_Cys = pose.pdb_info()->pdb2pose( 'H', 33 );

    for( Size ii = ab_info_->get_CDR_loop("h3")->start(); ii <= ab_info_->get_CDR_loop("h3")->stop(); ii++ )
        
        if( pose.residue(ii).name3() == "CYS" ) H3_Cys = ii;

    if( ( H1_Cys != 0 ) && ( H3_Cys != 0 ) )
        TR << "CONSTRAINTS: "<< "AtomPair CA " << H1_Cys << " CA " << H3_Cys
                << " BOUNDED 4.0 6.1 0.6 BOND; mean 5.6 sd 0.6" << std::endl;

    // Specifying extended kink

    Size hfr_46(0), h3_closest(0);
    hfr_46 = pose.pdb_info()->pdb2pose( 'H', 46 );
    if( ab_info_->is_extended() ) h3_closest = ab_info_->get_CDR_loop("h3")->stop() - 5;
    if( h3_closest != 0 )
        TR << "CONSTRAINTS: " << "AtomPair CA " << hfr_46 << " CA " << h3_closest
            << " BOUNDED 6.5 9.1 0.7 DISTANCE; mean 8.0 sd 0.7" << std::endl;

    return;
} // display_constraint_residues

    
    
    
    
    
    
    
    
    
    
/// @details  Show the complete setup of the docking protocol
void Ab_AssembleCDRs::show( std::ostream & out ) {
    if ( !flags_and_objects_are_in_sync_ ){
        sync_objects_with_flags();
    }
    out << *this;
}
    
std::ostream & operator<<(std::ostream& out, const Ab_AssembleCDRs & ab_m_2 )
{
    using namespace ObjexxFCL::fmt;
        
    // All output will be 80 characters - 80 is a nice number, don't you think?
    std::string line_marker = "///";
    out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
    out << line_marker << A( 47, "Rosetta 3 Antibody Modeler" ) << space( 27 ) << line_marker << std::endl;
    out << line_marker << A( 47, "     CDR Assembling       " ) << space( 27 ) << line_marker << std::endl;
    out << line_marker << space( 74 ) << line_marker << std::endl;
        
    // Display the state of the low resolution docking protocol that will be used
    out << line_marker << " Graft_l1:  " << ab_m_2.graft_l1_<<std::endl;
    out << line_marker << " Graft_l2:  " << ab_m_2.graft_l2_<<std::endl;
    out << line_marker << " Graft_l3:  " << ab_m_2.graft_l3_<<std::endl;
    out << line_marker << " Graft_h1:  " << ab_m_2.graft_h1_<<std::endl;
    out << line_marker << " Graft_h2:  " << ab_m_2.graft_h2_<<std::endl;
    out << line_marker << " Graft_h3:  " << ab_m_2.graft_h3_<<std::endl;
    out << line_marker << "  camelid:  " << ab_m_2.camelid_ <<std::endl;
        
        
    // Display the state of the low resolution docking protocol that will be used
        
        
    // Close the box I have drawn
    out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
    return out;
}
    

    



    
    
    
    

} // end antibody2
} // end protocols






