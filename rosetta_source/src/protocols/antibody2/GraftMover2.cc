// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file antibody2/GraftMover2.cc
/// @brief grafts a cdr onto the template of an antibody framework
/// @detailed
/// @author Jianqing Xu (xubest@gmail.com)


// Rosetta Headers
#include <protocols/antibody2/GraftMover2.hh>


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

#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/loops/LoopMover.fwd.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/antibody2/CDRH3Modeler2.hh>
#include <protocols/antibody2/GraftOneMover.fwd.hh>
#include <protocols/antibody2/GraftOneMover.hh>
#include <protocols/antibody2/CloseOneMover.fwd.hh>
#include <protocols/antibody2/CloseOneMover.hh>
#include <protocols/antibody2/LoopRlxMover.fwd.hh>
#include <protocols/antibody2/LoopRlxMover.hh>
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


static basic::Tracer TR("protocols.antibody2.GraftMover2");



namespace protocols {
namespace antibody2 {
using namespace core;

GraftMover2::GraftMover2() : protocols::moves::Mover()
{
	user_defined_ = false;
	init( false, false, false, false, false, false, false, false );
} // GraftMover2 default constructor

GraftMover2::GraftMover2( bool l1, bool l2, bool l3, bool h1, bool h2, bool h3, bool camelid, bool benchmark ) : Mover() {
	std::cout<<"I am here 1.1"<<std::endl;
	user_defined_ = true;
	init( l1, l2, l3, h1, h2, h3, camelid, benchmark );
}

// GraftMover2 default destructor
GraftMover2::~GraftMover2() {}

void GraftMover2::init( bool l1, bool l2, bool l3, bool h1, bool h2, bool h3, bool camelid, bool benchmark ) {
	Mover::type( "GraftMover2" );

	std::cout<<"I am here 1.1.1"<<std::endl;
	std::cout<<"l1="<<l1<<std::endl;
	std::cout<<"l2="<<l2<<std::endl;
	std::cout<<"l3="<<l3<<std::endl;
	std::cout<<"h1="<<h1<<std::endl;
	std::cout<<"h2="<<h2<<std::endl;
	std::cout<<"h3="<<h3<<std::endl;
	std::cout<<"camelid="<<camelid<<std::endl;
	std::cout<<"benchmark="<<benchmark<<std::endl;
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
GraftMover2::set_default() {
	std::cout<<"I am here 1.1.1.1"<<std::endl;
	TR <<  "Setting up default settings" << std::endl;
	graft_l1_ = false;
	graft_l2_ = false;
	graft_l3_ = false;
	graft_h1_ = false;
	graft_h2_ = false;
	graft_h3_ = false;
	camelid_ = false;

	// High resolution scores
	scorefxn_ = scoring::ScoreFunctionFactory::create_score_function( scoring::STANDARD_WTS );

	first_apply_with_current_setup_ = true;
} // GraftMover2 set_default

void GraftMover2::finalize_setup( AntibodyInfo & ab_info )
{
	graft_sequence_ = new protocols::moves::SequenceMover();
	relax_sequence_ = new protocols::moves::SequenceMover();
	pymol_ = new protocols::moves::PyMolMover();

	for ( GraftMap::const_iterator it = grafts_.begin(); it != grafts_.end(); ++it ) {
		if ( it->second ) {
			TR << "Creating movers for " << it->first << std::endl;
			protocols::antibody2::GraftOneMoverOP graftone ( new protocols::antibody2::GraftOneMover( ab_info.get_loop(it->first)->start(), ab_info.get_loop(it->first)->stop(), it->first, scorefxn_ ));
			graftone->enable_benchmark_mode( benchmark_ );
			graft_sequence_->add_mover( graftone );
			graft_sequence_->add_mover( pymol_ );

			CloseOneMoverOP closeone( new CloseOneMover( ab_info.get_loop(it->first)->start(),ab_info.get_loop(it->first)->stop() ) );
			closeone->enable_benchmark_mode( benchmark_ );
			closeone->set_pymol( pymol_ );
			graft_sequence_->add_mover( closeone );
			graft_sequence_->add_mover( pymol_ );

			LoopRlxMoverOP rlx_one_loop(new LoopRlxMover( ab_info.get_loop(it->first)->start(),ab_info.get_loop(it->first)->stop() ) );
			rlx_one_loop->enable_benchmark_mode( benchmark_ );
			relax_sequence_->add_mover( rlx_one_loop );
			relax_sequence_->add_mover( pymol_ );
		}
	}
}

void GraftMover2::apply( pose::Pose & pose )
{
	TR <<  "Grafting designated CDRs" << std::endl;

	std::cout<<"I am here 7.1"<<std::endl;
	AntibodyInfo ab_info( pose, camelid_ );
	std::cout<<"I am here 7.2"<<std::endl;
	if ( first_apply_with_current_setup_ ){ finalize_setup( ab_info ); first_apply_with_current_setup_ = false; }
	std::cout<<"I am here 7.3"<<std::endl;
	Size nres = pose.total_residue();

	// Storing secondary structure
	utility::vector1<char> secondary_struct_storage;
	for( Size i = 1; i <= nres; i++ ) secondary_struct_storage.push_back( pose.secstruct( i ) );

	std::cout<<"I am here 7.4"<<std::endl;
	graft_sequence_->apply( pose );

	std::cout<<"I am here 7.5"<<std::endl;
	if ( !graft_h3_ ) {
		TR << "Extending CDR H3" << std::endl;
		loops::Loop cdr_h3( ab_info.get_loop("h3")->start()-1, ab_info.get_loop("h3")->stop(), ab_info.get_loop("h3")->start(), 0, false );
		simple_one_loop_fold_tree( pose, cdr_h3);

		// silly hack to make extended loops work
		loops::Loops loop_list;
		loop_list.add_loop( cdr_h3 );

		loops::LoopMoverOP my_loop_move( new loops::LoopMover( loop_list ) );
		my_loop_move->set_extended_torsions( pose, cdr_h3 );
	}

	std::cout<<"I am here 7.6"<<std::endl;
	// Recover secondary structures
	for( Size i = 1; i <= nres; i++ ) pose.set_secstruct( i, secondary_struct_storage[ i ] );

	// relax optimized CDR grafted regions
	relax_sequence_->apply( pose );

	std::cout<<"I am here 7.7"<<std::endl;
	// Recover secondary structures
	for( Size i = 1; i <= nres; i++ ) pose.set_secstruct( i, secondary_struct_storage[ i ] );

	// align pose to native pose
	pose::Pose native_pose;
	if( get_native_pose() )
		native_pose = *get_native_pose();
	else
		native_pose = pose;
	antibody2::AntibodyInfo native_ab( native_pose, camelid_ );
	std::cout<<"I am here 7.8"<<std::endl;

	ab_info.align_to_native( pose, native_ab, native_pose );
} // GraftMover2::apply()

std::string
GraftMover2::get_name() const { return "GraftMover2"; }

// copy ctor
GraftMover2::GraftMover2( GraftMover2 const & rhs ) {
    initForEqualOperatorAndCopyConstructor(*this, rhs);
}

///@brief assignment operator
GraftMover2 & GraftMover2::operator=( GraftMover2 const & rhs ){
    //abort self-assignment
    if (this == &rhs) return *this;
    Mover::operator=(rhs);
    initForEqualOperatorAndCopyConstructor(*this, rhs);
    return *this;
}

void GraftMover2::initForEqualOperatorAndCopyConstructor(GraftMover2 & lhs, GraftMover2 const & rhs)
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
    lhs.graft_sequence_ = new protocols::moves::SequenceMover( *( rhs.graft_sequence_ ) );
    lhs.relax_sequence_ = new protocols::moves::SequenceMover( *( rhs.relax_sequence_ ) );
    lhs.packer_ = new protocols::simple_moves::PackRotamersMover( *( rhs.packer_ ) );
    lhs.pymol_ = new protocols::moves::PyMolMover( *( rhs.pymol_ ) );
    lhs.scorefxn_ = new core::scoring::ScoreFunction( *( rhs.scorefxn_ ) );
}

void GraftMover2::set_packer_default(
	pose::Pose & pose,
	bool include_current) {

	//set up packer
	pack::task::PackerTaskOP task;
	task = pack::task::TaskFactory::create_packer_task( pose );
	task->restrict_to_repacking();
	task->or_include_current( include_current );
	packer_ = new protocols::simple_moves::PackRotamersMover( scorefxn_, task );

} // GraftMover2 set_packer_default

void GraftMover2::relax_optimized_CDR_grafts( pose::Pose & pose ) {
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
			LoopRlxMoverOP rlx_one_loop(new LoopRlxMover( loop_begin, loop_end));
			rlx_one_loop->enable_benchmark_mode( benchmark_ );
			rlx_one_loop->apply( pose );
			loop_begin = 0;
			loop_end = 0;
		}
	} // for ii <= nres
} // GraftMover2::relax_optimized_CDR_grafts





}  // namespace antibody2
}  // namespace protocols
