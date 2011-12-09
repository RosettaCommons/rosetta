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

#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/types.hh>
#include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/FoldTree.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>


#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/CcdLoopClosureMover.hh>
#include <protocols/loops/LoopMover.fwd.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/antibody2/CDRH3Modeler2.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <core/pose/util.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


static basic::Tracer TR("protocols.antibody2.GraftMover2");
static basic::Tracer TRG("protocols.antibody2.GraftOneMover");
static basic::Tracer TRC("protocols.antibody2.CloseOneMover");
static basic::Tracer TRL("protocols.antibody2.LoopRlxMover");

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
			GraftOneMoverOP graftone ( new GraftOneMover( ab_info.get_loop(it->first)->start(), ab_info.get_loop(it->first)->stop(), it->first, scorefxn_ ));
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
GraftMover2::get_name() const {
	return "GraftMover2";
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
















//**************************************************************************************
//**************************************************************************************
//**************************************************************************************
//**************************************************************************************
//**************************************************************************************








GraftOneMover::GraftOneMover(
	Size query_start,
	Size query_end,
	std::string template_name,
	scoring::ScoreFunctionOP scorefxn ) : Mover( "GraftOneMover" ), scorefxn_( scorefxn )
{
	query_start_ = query_start;
	query_end_ = query_end;
	set_default( template_name );
} // GraftOneMover default constructor

void GraftOneMover::set_default( std::string template_name )
{
	std::string const path = basic::options::option[ basic::options::OptionKeys::in::path::path ]()[1];
	TRG << "Reading in template: " << path << template_name << ".pdb " << std::endl;
	core::import_pose::pose_from_pdb( template_pose_, path + template_name + ".pdb" );
	antibody2::AntibodyInfo ab_info( template_pose_, template_name );
	template_start_ = ab_info.current_start;
	template_end_ = ab_info.current_end;
	template_name_ = template_name;
} // GraftOneMover::set_default

std::string
GraftOneMover::get_name() const {
	return "GraftOneMover";
}





void GraftOneMover::apply( pose::Pose & pose_in )
{
	TRG<<"step1        I am here 7.4.1"<<std::endl;
	Size const nres( pose_in.total_residue() ); // Total residues
	Size query_size = ( query_end_ - query_start_ )+1;
	Size flank_size ( 5 );

	//set up packer
	pack::task::PackerTaskOP task;
	task = pack::task::TaskFactory::create_packer_task( pose_in );
	task->restrict_to_repacking();
	task->or_include_current( false );
	utility::vector1< bool > allow_repack( nres, false );
	for ( Size i = query_start_; i <= query_end_; ++i ) allow_repack[ i ] = true;
	task->restrict_to_residues( allow_repack );
	protocols::simple_moves::PackRotamersMoverOP packer = new protocols::simple_moves::PackRotamersMover( scorefxn_, task );

	// create a sub pose with  5 flanking residues on either side of CDR loop
	pose::Pose truncated_pose( pose_in, query_start_-flank_size, query_end_+flank_size );
	truncated_pose.conformation().delete_residue_range_slow( flank_size+1, ( query_size ) + flank_size );

	// create atom map for superimposing 2 flanking resiudes
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, template_pose_, id::BOGUS_ATOM_ID );

	for( Size start_stem = 1; start_stem < flank_size; ++start_stem ) {
		Size const ref_stem ( start_stem+1 ); // remember there are 5 flanking residues
		for( Size j=1; j <= 4; j++ ) {
			id::AtomID const id1( j, start_stem );
			id::AtomID const id2( j, ref_stem );
			atom_map[ id1 ] = id2;
		}
	}

	// start at the end of the actual loop
	for( Size end_stem = query_size+flank_size; end_stem < query_size+flank_size+4; ++end_stem ) {
		Size const ref_stem ( end_stem-query_size+1 );
		for( Size j=1; j <= 4; j++ ) {
			id::AtomID const id1( j, end_stem );
			id::AtomID const id2( j, ref_stem );
			atom_map[ id1 ] = id2;
		}
	}
	scoring::superimpose_pose( template_pose_, truncated_pose, atom_map );
	template_pose_.dump_pdb(template_name_);

	for ( Size i=query_start_; i <= query_end_; ++i ) {
		Size template_num ( i - (query_start_-5) );
		core::conformation::Residue const & source_rsd( pose_in.residue( i ) );
		core::conformation::Residue const & target_rsd( template_pose_.residue( template_num) );

		Size const natoms( source_rsd.natoms() );

		bool any_missing( false );
		id::AtomID_Mask missing( false );
		// dimension the missing-atom mask
		core::pose::initialize_atomid_map( missing, pose_in );

		if( source_rsd.name() != target_rsd.name() )
			pose_in.set_secstruct( i, 'X' );
		for ( Size j=1; j<= natoms; ++j ) {
			std::string const & atom_name( source_rsd.atom_name(j) );
			if ( target_rsd.has( atom_name ) ) {
				pose_in.set_xyz( id::AtomID( source_rsd.atom_index( atom_name),i ), target_rsd.xyz( atom_name ) );
			} else {
				any_missing = true;
				missing[ id::AtomID( pose_in.residue_type(i).atom_index(source_rsd.atom_name(j)), i ) ] = true;
			}
		}

		if ( any_missing ) {
			pose_in.conformation().fill_missing_atoms( missing );
		}
	}
	packer->apply( pose_in );
	pose_in.dump_pdb(template_name_+"_graft");
} // GraftOneMover::apply













//**************************************************************************************
//**************************************************************************************
//**************************************************************************************
//**************************************************************************************
//**************************************************************************************


CloseOneMover::CloseOneMover(
	Size query_start,
	Size query_end ) : Mover( "CloseOneMover" )
{
	set_default();
	cdr_loop_start_ = query_start;
	cdr_loop_end_ = query_end;
	loop_start_ = query_start-flanking_residues_;
	loop_end_ = query_end+flanking_residues_;
} // CloseOneMover default constructor

void CloseOneMover::set_default()
{
	allowed_separation_ = 1.9;
	flanking_residues_ = 5; // default 5;
	movemap_ = new kinematics::MoveMap();
	movemap_->set_chi( false );
	movemap_->set_bb( false );
	pymol_ = new protocols::moves::PyMolMover();
} // CloseOneMover::set_default

void CloseOneMover::apply( pose::Pose & pose_in )
{
	TRC<<"step 2         I am here 7.4.2"<<std::endl;
	Size const N ( 1 ); // N atom
	Size const C ( 3 ); // C atom

	// Coordinates of the C and N atoms at stem
	numeric::xyzVector_float peptide_C, peptide_N;
	// N-terminal
	peptide_C = pose_in.residue( cdr_loop_start_ - 1 ).xyz( C );
	peptide_N = pose_in.residue( cdr_loop_start_ ).xyz( N );

	// C-terminal
	peptide_C = pose_in.residue( cdr_loop_end_ ).xyz( C );
	peptide_N = pose_in.residue( cdr_loop_end_ + 1 ).xyz( N );

	// calculate separation at ends to see if it needs to be closed
//	Real nter_separation=distance(peptide_C, peptide_N);
//	Real cter_separation=distance(peptide_C, peptide_N);
	Real nter_separation=peptide_C.distance(peptide_N);
	Real cter_separation=peptide_C.distance(peptide_N);

	// save the starting foldtree
	core::kinematics::FoldTree f( pose_in.fold_tree() );

	// setup movemap to only loop residues
	utility::vector1< bool> allow_bb_move( pose_in.total_residue(), false );
	for ( Size i=loop_start_; i<= loop_end_; ++i )
		allow_bb_move[ i ] = true;
	movemap_->set_bb( allow_bb_move );
	movemap_->set_jump( 1, false );

	pymol_->apply( pose_in );

	if( nter_separation > allowed_separation_ ) {
		loops::Loop one_loop( loop_start_, cdr_loop_start_, cdr_loop_start_-1, 0, false );
		simple_one_loop_fold_tree( pose_in, one_loop );
		loops::CcdMoverOP ccd_moves = new loops::CcdMover( one_loop, movemap_ );
		ccd_moves->apply( pose_in );
		pymol_->apply( pose_in );
	}

	if( cter_separation > allowed_separation_ ) {
		loops::Loop one_loop( cdr_loop_end_, loop_end_, cdr_loop_end_+1, 0, false );
		simple_one_loop_fold_tree( pose_in, one_loop );
		loops::CcdMoverOP ccd_moves = new loops::CcdMover( one_loop, movemap_ );
		ccd_moves->apply( pose_in );
		pymol_->apply( pose_in );
	}

	Real separation = 0.00;
	for( Size ii = loop_start_; ii <= loop_end_; ii++ ) {
		peptide_C = pose_in.residue( ii ).xyz( C );
		peptide_N = pose_in.residue( ii + 1 ).xyz( N );
		separation=peptide_C.distance(peptide_N);
//		separation=distance(peptide_C, peptide_N);
		if( separation > allowed_separation_ ) {
			Size cutpoint = ii;
			loops::Loop one_loop( loop_start_, loop_end_, cutpoint, 0, false );
			loops::CcdMoverOP ccd_moves = new loops::CcdMover( one_loop, movemap_ );
			ccd_moves->apply( pose_in );
			pymol_->apply( pose_in );
		}
	} // for

	// reset to original foldtree
	pose_in.fold_tree( f );
} // CloseOneMover::apply

std::string
CloseOneMover::get_name() const
{
	return "CloseOneMover";
}













//**************************************************************************************
//**************************************************************************************
//**************************************************************************************
//**************************************************************************************
//**************************************************************************************

LoopRlxMover::LoopRlxMover(
	Size query_start,
	Size query_end ) : Mover( "LoopRlxMover" )
{
	loop_start_ = query_start;
	loop_end_ = query_end;
	set_default();
} // LoopRlxMover default constructor

void LoopRlxMover::set_default()
{
	highres_scorefxn_ = scoring::ScoreFunctionFactory::create_score_function( "standard", "score12" );
	highres_scorefxn_->set_weight( scoring::chainbreak, 1.0 );
	highres_scorefxn_->set_weight( scoring::overlap_chainbreak, 10./3. );

	movemap_= new kinematics::MoveMap();
	movemap_->set_chi( false );
	movemap_->set_bb( false );

	mc_ = new moves::MonteCarlo( *highres_scorefxn_, 0.8 );
	tf_ = new core::pack::task::TaskFactory;

	Real const init_temp( 2.0 );
	Real const last_temp( 0.5 );
	inner_cycles_ = (loop_end_ - loop_start_) + 1;
	outer_cycles_ = 2;
	gamma_ = std::pow( (last_temp/init_temp), (1.0/inner_cycles_));
	temperature_ = init_temp;
} // LoopRlxMover::set_default

void LoopRlxMover::setup_objects( pose::Pose & pose )
{
	using namespace protocols::moves;
	using namespace protocols::loops;

	//setting MoveMap
	utility::vector1< bool> allow_bb_move( pose.total_residue(), false );
	for( Size ii = loop_start_; ii <= loop_end_; ii++ )
		allow_bb_move[ ii ] = true;
	movemap_->set_bb( allow_bb_move );
	movemap_->set_jump( 1, false );

	Size loop_size = ( loop_end_ - loop_start_ ) + 1;
	Size cutpoint = loop_start_ + int(loop_size/2);

	one_loop_ = new Loop( loop_start_, loop_end_, cutpoint, 0, false );
	simple_one_loop_fold_tree( pose, *one_loop_ );

	// set cutpoint variants for correct chainbreak scoring
	if( !pose.residue( cutpoint ).is_upper_terminus() ) {
		if( !pose.residue( cutpoint ).has_variant_type( chemical::CUTPOINT_LOWER ) )
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, cutpoint );
		if( !pose.residue( cutpoint + 1 ).has_variant_type( chemical::CUTPOINT_UPPER ) )
			core::pose::add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, cutpoint + 1 );
	}

	Real min_tolerance = 0.001;
	if( benchmark_ ) min_tolerance = 1.0;
	std::string min_type = "dfpmin_armijo_nonmonotone";
	bool nb_list = true;
	loop_min_mover_ = new simple_moves::MinMover( movemap_, highres_scorefxn_, min_type, min_tolerance, nb_list );

	// more params
	Size n_small_moves ( numeric::max(Size(5), Size(loop_size/2)) );
	if( benchmark_ ) {
		n_small_moves = 1;
		inner_cycles_ = 1;
		outer_cycles_ = 1;
	}

	Real high_move_temp = 2.00;
	// minimize amplitude of moves if correct parameter is set
    simple_moves::BackboneMoverOP small_mover = new simple_moves::SmallMover( movemap_, high_move_temp, n_small_moves );
    simple_moves::BackboneMoverOP shear_mover = new simple_moves::ShearMover( movemap_, high_move_temp, n_small_moves );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 5.0 );
	small_mover->angle_max( 'L', 6.0 );
	shear_mover->angle_max( 'H', 2.0 );
	shear_mover->angle_max( 'E', 5.0 );
	shear_mover->angle_max( 'L', 6.0 );

	CcdMoverOP ccd_moves = new CcdMover( *one_loop_, movemap_ );
	RepeatMoverOP ccd_cycle = new RepeatMover(ccd_moves, n_small_moves);

	wiggle_loop_ = new protocols::moves::SequenceMover();
	wiggle_loop_->add_mover( small_mover );
	wiggle_loop_->add_mover( shear_mover );
	wiggle_loop_->add_mover( ccd_cycle );
	wiggle_loop_->add_mover( loop_min_mover_ );
}


///////////////////////////////////////////////////////////////////////////
/// @begin LoopRlxMover::apply
///
/// @brief relaxes the region specified
///
/// @detailed This is all done in high resolution.Hence there are no rigid
///           body moves relative to the docking partners. Only small moves
///           are carried out here to see if there are better fits.
///           Repacking is carried out extensively after each move.
///
/// @param[in] pose, loop begin position, loop end position
///
/// @global_read none
///
/// @global_write none
///
/// @remarks
///
/// @references
///
/// @authors Aroop 05/12/2010
///
/// @last_modified 05/12/2010
///////////////////////////////////////////////////////////////////////////
void LoopRlxMover::apply( pose::Pose & pose_in )
{
	using namespace protocols;
	using namespace protocols::loops;
	using namespace protocols::moves;
	using namespace pack;
	using namespace pack::task;
	using namespace pack::task::operation;

	TRL<<"I am here 7.4.3"<<std::endl;
	TR << "LoopRlxMover: Apply" << std::endl;

	if ( !pose_in.is_fullatom() )
		utility_exit_with_message("Fullatom poses only");

	setup_objects( pose_in );

	// storing starting fold tree
	kinematics::FoldTree tree_in( pose_in.fold_tree() );

	utility::vector1< bool> allow_repack( pose_in.total_residue(), false );
	select_loop_residues( pose_in, *one_loop_, true /*include_neighbors*/, allow_repack);
	movemap_->set_chi( allow_repack );

    simple_moves::PackRotamersMoverOP loop_repack = new simple_moves::PackRotamersMover(highres_scorefxn_);
	setup_packer_task( pose_in );
	( *highres_scorefxn_ )( pose_in );
	tf_->push_back( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) );
	loop_repack->task_factory(tf_);
	loop_repack->apply( pose_in );

	// rotamer trials - select loop with new neighbors after packing occurs
	select_loop_residues( pose_in, *one_loop_, true /*include_neighbors*/, allow_repack);
	movemap_->set_chi( allow_repack );
	setup_packer_task( pose_in );
	( *highres_scorefxn_ )( pose_in );
	tf_->push_back( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) );
    simple_moves::RotamerTrialsMoverOP pack_rottrial = new simple_moves::RotamerTrialsMover( highres_scorefxn_, tf_ );

	pack_rottrial->apply( pose_in );

	mc_->reset( pose_in ); // monte carlo reset

	// outer cycle
	for(Size i = 1; i <= outer_cycles_; i++) {
		mc_->recover_low( pose_in );

		// inner cycle
		for ( Size j = 1; j <= inner_cycles_; j++ ) {
			temperature_ *= gamma_;
			mc_->set_temperature( temperature_ );
			wiggle_loop_->apply( pose_in );

			// rotamer trials- select loop with new neighbors
			select_loop_residues( pose_in, *one_loop_, true /*include_neighbors*/, allow_repack);
			movemap_->set_chi( allow_repack );
			setup_packer_task( pose_in );
			( *highres_scorefxn_ )( pose_in );
			tf_->push_back( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) );
            simple_moves::RotamerTrialsMoverOP pack_rottrial = new simple_moves::RotamerTrialsMover( highres_scorefxn_, tf_ );
			pack_rottrial->apply( pose_in );
			mc_->boltzmann( pose_in );

			if ( numeric::mod(j,Size(20))==0 || j==inner_cycles_ ) {
				// repack trial
				loop_repack = new simple_moves::PackRotamersMover( highres_scorefxn_ );
				setup_packer_task( pose_in );
				( *highres_scorefxn_ )( pose_in );
				tf_->push_back( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) );
				loop_repack->task_factory( tf_ );
				loop_repack->apply( pose_in );
				mc_->boltzmann( pose_in );
			}
		} // inner cycles
	} // outer cycles
	mc_->recover_low( pose_in );

	// minimize
	if( !benchmark_ )
		loop_min_mover_->apply( pose_in );

	// Restoring pose stuff
	pose_in.fold_tree( tree_in ); // Tree

	TR << "LoopRlxMover: Finished Apply" << std::endl;
} // LoopRlxMover::apply

std::string
LoopRlxMover::get_name() const
{
	return "LoopRlxMover";
}
void
LoopRlxMover::setup_packer_task(
	pose::Pose & pose_in )
{
	using namespace pack::task;
	using namespace pack::task::operation;
	tf_ = new TaskFactory;

	TR << "LoopRlxMover Setting Up Packer Task" << std::endl;

	tf_->push_back( new OperateOnCertainResidues( new PreventRepackingRLT, new ResidueLacksProperty("PROTEIN") ) );
	tf_->push_back( new InitializeFromCommandline );
	tf_->push_back( new IncludeCurrent );
	tf_->push_back( new RestrictToRepacking );
	tf_->push_back( new NoRepackDisulfides );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	pack::rotamer_set::UnboundRotamersOperationOP unboundrot =
		new pack::rotamer_set::UnboundRotamersOperation();
	unboundrot->initialize_from_command_line();
	operation::AppendRotamerSetOP unboundrot_operation =
		new operation::AppendRotamerSet( unboundrot );
	tf_->push_back( unboundrot_operation );
	// adds scoring bonuses for the "unbound" rotamers, if any
	core::pack::dunbrack::load_unboundrot( pose_in );

	TR << "LoopRlxMover Done: Setting Up Packer Task" << std::endl;

} // LoopRlxMover::setup_packer_task



}  // namespace antibody2
}  // namespace protocols
