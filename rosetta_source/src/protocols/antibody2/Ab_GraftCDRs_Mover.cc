// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file antibody2/Ab_GraftCDRs_Mover.cc
/// @brief grafts a cdr onto the template of an antibody framework
/// @detailed
/// @author Jianqing Xu (xubest@gmail.com)


// Rosetta Headers
#include <protocols/antibody2/Ab_GraftCDRs_Mover.hh>


#include <core/conformation/Conformation.hh>
#include <core/id/types.hh>
#include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>

#include <basic/Tracer.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/antibody2/Ab_Info.hh>
//#include <protocols/loops/LoopMover.fwd.hh>
//#include <protocols/loops/LoopMover.hh>
#include <protocols/antibody2/CDRH3Modeler2.hh>
#include <protocols/antibody2/Ab_GraftOneCDR_Mover.fwd.hh>
#include <protocols/antibody2/Ab_GraftOneCDR_Mover.hh>
#include <protocols/antibody2/Ab_CloseOneCDR_Mover.fwd.hh>
#include <protocols/antibody2/Ab_CloseOneCDR_Mover.hh>
#include <protocols/antibody2/Ab_RelaxCDRs_Mover.fwd.hh>
#include <protocols/antibody2/Ab_RelaxCDRs_Mover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>



#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <utility/vector1.hh>

#ifdef WIN32
#include <protocols/antibody2/Ab_TemplateInfo.hh>
#endif


static basic::Tracer TR("protocols.antibody2.Ab_GraftCDRs_Mover");



namespace protocols {
namespace antibody2 {
using namespace core;

Ab_GraftCDRs_Mover::Ab_GraftCDRs_Mover() : moves::Mover()
{
	user_defined_ = false;
	init( false, false, false, false, false, false, false, false );
} // Ab_GraftCDRs_Mover default constructor

    
Ab_GraftCDRs_Mover::Ab_GraftCDRs_Mover(bool l1,bool l2,bool l3,bool h1,bool h2,bool h3,bool camelid,bool benchmark ) : Mover() {
	user_defined_ = true;
	init(l1,l2,l3,h1,h2,h3,camelid,benchmark );
}


Ab_GraftCDRs_Mover::Ab_GraftCDRs_Mover(Ab_InfoCOP ab_info, Ab_TemplateInfoCOP ab_template): Mover(){
    user_defined_ = true;
//    init(l1,l2,l3,h1,h2,h3,camelid,benchmark );
}
    
    
    
    
// Ab_GraftCDRs_Mover default destructor
Ab_GraftCDRs_Mover::~Ab_GraftCDRs_Mover() {}

    
    
    
    
    
    
void Ab_GraftCDRs_Mover::init(bool l1,bool l2,bool l3,bool h1,bool h2,bool h3,bool camelid,bool benchmark){
	Mover::type("Ab_GraftCDRs_Mover");

	// setup all the booleans with default values
	// they will get overwritten by the options and/or passed values
    
	set_default();
    
//	register_options();
//	init_from_options();
	if ( user_defined_ ) {
		graft_l1_ = l1;
		graft_l2_ = l2;
		graft_l3_ = l3;
		graft_h1_ = h1;
		graft_h2_ = h2;
		graft_h3_ = h3;
		camelid_ = camelid;
		benchmark_ = benchmark;
        TR<<"User defined values: "<<std::endl;
        TR<<"graft_l1_="<<graft_l1_<<std::endl;
        TR<<"graft_l2_="<<graft_l2_<<std::endl;
        TR<<"graft_l3_="<<graft_l3_<<std::endl;
        TR<<"graft_h1_="<<graft_h1_<<std::endl;
        TR<<"graft_h2_="<<graft_h2_<<std::endl;
        TR<<"graft_h3_="<<graft_h3_<<std::endl;
        TR<<"camelid_=" <<camelid_ <<std::endl;
        TR<<"benchmark_="<<benchmark_<<std::endl;
	}

	// ensure that if camelid is set, no l-loops are set to true
	if ( camelid_ ) graft_l1_=graft_l2_=graft_l3_=false;

	grafts_.insert( std::pair< std::string, bool >("l1", graft_l1_) );
	grafts_.insert( std::pair< std::string, bool >("l2", graft_l2_) );
	grafts_.insert( std::pair< std::string, bool >("l3", graft_l3_) );
	grafts_.insert( std::pair< std::string, bool >("h1", graft_h1_) );
	grafts_.insert( std::pair< std::string, bool >("h2", graft_h2_) );
	grafts_.insert( std::pair< std::string, bool >("h3", graft_h3_) );

	// set up objects based on the boolean values defined above
//	setup_objects();
}

    
    
    
    
    
    
void
Ab_GraftCDRs_Mover::set_default() {
	TR <<  "Setting up default settings, setting everything to false" << std::endl;
	graft_l1_=graft_l2_=graft_l3_=graft_h1_=graft_h2_=graft_h3_=false;
    camelid_  = false;

	// High resolution scores
	scorefxn_ = scoring::ScoreFunctionFactory::create_score_function( scoring::STANDARD_WTS );

	first_apply_with_current_setup_ = true;
} // Ab_GraftCDRs_Mover set_default

    
    
    
    
    
    
void Ab_GraftCDRs_Mover::finalize_setup(core::pose::Pose & frame_pose, Ab_Info & ab_info )
{
    graft_sequence_ = new moves::SequenceMover();
    relax_sequence_ = new moves::SequenceMover();
    pymol_ = new moves::PyMolMover();

    
    Size nres = frame_pose.total_residue();
    
    for ( GraftMap::const_iterator it = grafts_.begin(); it != grafts_.end(); ++it ) {
        if ( it->second ) {
                    TR << "Creating movers for " << it->first << std::endl;
                    TR << "                  start (chothia): "<<ab_info.get_loop(it->first)->start()<<std::endl;
                    TR << "                   stop (chothia): "<<ab_info.get_loop(it->first)->stop()<<std::endl;
            
                Ab_GraftOneCDR_MoverOP graftone ( new Ab_GraftOneCDR_Mover( ab_info.get_loop(it->first)->start(), 
                                                              ab_info.get_loop(it->first)->stop(), 
                                                              it->first, scorefxn_ )   );

            
                graftone->enable_benchmark_mode( benchmark_ );
                graft_sequence_->add_mover( graftone );
//              graft_sequence_->add_mover( pymol_ );
                
            
                /*
                Ab_CloseOneCDR_MoverOP closeone( new Ab_CloseOneCDR_Mover( ab_info.get_loop(it->first)->start(),
                                                             ab_info.get_loop(it->first)->stop()   )     );
                closeone->enable_benchmark_mode( benchmark_ );
                closeone->set_pymol( pymol_ );
                graft_sequence_->add_mover( closeone );
                graft_sequence_->add_mover( pymol_ );
                 */
            

                Ab_RelaxCDRs_MoverOP rlx_one_loop(new Ab_RelaxCDRs_Mover( ab_info.get_loop(it->first)->start(),
                                                              ab_info.get_loop(it->first)->stop()   )    );
                rlx_one_loop->enable_benchmark_mode( benchmark_ );
                relax_sequence_->add_mover( rlx_one_loop );
                relax_sequence_->add_mover( pymol_ );
        }
    }
    
    // Exact match Aroop's old code in Rosetta 2:
    // graft all CDRs by superimpose stems, then pack the whole new pose
    
    
    set_packer_default(frame_pose, true /* include_current */)  ;
    graft_sequence_->add_mover(packer_);  
    
    
}

    
    
    
    
    
    
    
void Ab_GraftCDRs_Mover::apply( pose::Pose & frame_pose )
{
	TR <<  "Grafting designated CDRs" << std::endl;


	Ab_Info ab_info( frame_pose, camelid_ );

	if ( first_apply_with_current_setup_ ){ 
	    finalize_setup(frame_pose, ab_info ); 
	    first_apply_with_current_setup_ = false; 
	}

	Size nres = frame_pose.total_residue();

	// Storing secondary structure
	utility::vector1<char> secondary_struct_storage;
	for( Size i = 1; i <= nres; i++ ) {
        secondary_struct_storage.push_back( frame_pose.secstruct( i ) );
        TR<<"JQX:   residue: "<<i<<"       secstruct: "<<frame_pose.secstruct(i)<<std::endl;
    }

	graft_sequence_->apply( frame_pose );
    
        frame_pose.dump_pdb("finish_grafting_and_packing.pdb");
        exit(-1); //JQX

	if ( !graft_h3_ ) {
		TR << "Extending CDR H3" << std::endl;
		loops::Loop cdr_h3( ab_info.get_loop("h3")->start()-1, 
                            ab_info.get_loop("h3")->stop(), 
                            ab_info.get_loop("h3")->start(), 0, false );
		simple_one_loop_fold_tree( frame_pose, cdr_h3);

		// silly hack to make extended loops work
		loops::LoopsOP loop_list = new loops::Loops();
		loop_list->add_loop( cdr_h3 );
/* Commented out by BDW with JX's consent
		loops::LoopMoverOP my_loop_move = new loops::LoopMover( loop_list );
		my_loop_move->set_extended_torsions( frame_pose, cdr_h3 );
*/
	}


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


	ab_info.align_to_native( frame_pose, native_ab, native_pose );
} // Ab_GraftCDRs_Mover::apply()

    
    
    
    
    
    
    
    
    
    
    
    
std::string
Ab_GraftCDRs_Mover::get_name() const { return "Ab_GraftCDRs_Mover"; }

// copy ctor
Ab_GraftCDRs_Mover::Ab_GraftCDRs_Mover( Ab_GraftCDRs_Mover const & rhs ) {
    initForEqualOperatorAndCopyConstructor(*this, rhs);
}

///@brief assignment operator
Ab_GraftCDRs_Mover & Ab_GraftCDRs_Mover::operator=( Ab_GraftCDRs_Mover const & rhs ){
    //abort self-assignment
    if (this == &rhs) return *this;
    Mover::operator=(rhs);
    initForEqualOperatorAndCopyConstructor(*this, rhs);
    return *this;
}

void Ab_GraftCDRs_Mover::initForEqualOperatorAndCopyConstructor(Ab_GraftCDRs_Mover & lhs, Ab_GraftCDRs_Mover const & rhs)
{
    lhs.graft_l1_ = rhs.graft_l1_;
    lhs.graft_l2_ = rhs.graft_l2_;
    lhs.graft_l3_ = rhs.graft_l3_;
    lhs.graft_h1_ = rhs.graft_h1_;
    lhs.graft_h2_ = rhs.graft_h2_;
    lhs.graft_h3_ = rhs.graft_h3_;
    
    lhs.grafts_ = rhs.grafts_;
    
    lhs.benchmark_ = rhs.benchmark_;
    lhs.camelid_ = rhs.camelid_;
    lhs.user_defined_ = rhs.user_defined_;
    lhs.first_apply_with_current_setup_ = rhs.first_apply_with_current_setup_;
    lhs.graft_sequence_ = new moves::SequenceMover( *( rhs.graft_sequence_ ) );
    lhs.relax_sequence_ = new moves::SequenceMover( *( rhs.relax_sequence_ ) );
    lhs.packer_ = new simple_moves::PackRotamersMover( *( rhs.packer_ ) );
    lhs.pymol_ = new moves::PyMolMover( *( rhs.pymol_ ) );
    lhs.scorefxn_ = new core::scoring::ScoreFunction( *( rhs.scorefxn_ ) );
}

    
    
    
    
    
    
    
    
    
void Ab_GraftCDRs_Mover::set_packer_default(pose::Pose & pose, bool include_current) {

	//set up packer
	pack::task::PackerTaskOP task;
	task = pack::task::TaskFactory::create_packer_task( pose );
	task->restrict_to_repacking();
	task->or_include_current( include_current );
	packer_ = new simple_moves::PackRotamersMover( scorefxn_, task );

} // Ab_GraftCDRs_Mover set_packer_default

    
    
    
    
    
    
    
    
    
    
    
    
    
    
void Ab_GraftCDRs_Mover::relax_optimized_CDR_grafts( pose::Pose & pose ) {
	Size loop_begin(0), loop_end(0);
	bool detect_flag( false );
	for( Size ii = 1; ii <= pose.total_residue(); ii++ ) {
		if( (pose.secstruct(ii) == 'Y') && !detect_flag ) {
			loop_begin = ii;
			detect_flag = true;
		}
		if( (pose.secstruct(ii) != 'Y') && detect_flag ) {
			loop_end = ii - 1;
			detect_flag = false;
		}
		if((detect_flag == false) && (loop_begin != 0) && (loop_end != 0 )) {
			Ab_RelaxCDRs_MoverOP rlx_one_loop(new Ab_RelaxCDRs_Mover( loop_begin, loop_end));
			rlx_one_loop->enable_benchmark_mode( benchmark_ );
			rlx_one_loop->apply( pose );
			loop_begin = 0;
			loop_end = 0;
		}
	} // for ii <= nres
} // Ab_GraftCDRs_Mover::relax_optimized_CDR_grafts





}  // namespace antibody2
}  // namespace protocols
