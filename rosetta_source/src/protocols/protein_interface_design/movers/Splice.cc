// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/Splice.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/Splice.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/protein_interface_design/movers/SpliceCreator.hh>
#include <utility/string_util.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/AA.hh>
#include <protocols/protein_interface_design/filters/TorsionFilter.hh>
#include <boost/foreach.hpp>
#include <protocols/rigid/RB_geometry.hh>
#define foreach BOOST_FOREACH
// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/protein_interface_design/movers/AddChainBreak.hh>
#include <utility/io/izstream.hh>
#include <iostream>
#include <sstream>
#include <algorithm>
//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <protocols/toolbox/task_operations/ThreadSequenceOperation.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <numeric/xyzVector.hh>
#include <protocols/loops/FoldTreeFromLoopsWrapper.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <numeric/random/random.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

namespace protocols {
namespace protein_interface_design {
namespace movers {

static basic::Tracer TR( "protocols.protein_interface_design.movers.Splice" );
static basic::Tracer TR_ccd( "protocols.protein_interface_design.movers.Splice_ccd" );
static numeric::random::RandomGenerator RG( 78289 );
std::string
SpliceCreator::keyname() const
{
	return SpliceCreator::mover_name();
}

protocols::moves::MoverOP
SpliceCreator::create_mover() const {
	return new Splice;
}

std::string
SpliceCreator::mover_name()
{
	return "Splice";
}

Splice::Splice() :
	Mover( SpliceCreator::mover_name() ),
	from_res_( 0 ),
	to_res_( 0 ),
	saved_from_res_( 0 ),
	saved_to_res_( 0 ),
	source_pdb_( "" ),
	ccd_( true ),
	rms_cutoff_( 999999 ),
	res_move_( 4 ),
	randomize_cut_( false ),
	task_factory_( NULL ),
	torsion_database_fname_( "" ),
	database_entry_( 0 ),
	template_file_( "" ),
	poly_ala_( true ),
	equal_length_( false ),
	template_pose_( NULL ),
	start_pose_( NULL ),
	saved_fold_tree_( NULL ),
	design_( false )
{
	torsion_database_.clear();
}


Splice::~Splice() {}
/* this doesnt quite work, and is superseded by copy_stretch
void
remove_cutpoint_variants_in_interval( core::pose::Pose & pose, core::Size const from_res, core::Size const to_res ){
	using namespace core::chemical;
	using namespace core::pose;

	const core::kinematics::FoldTree& tree(pose.fold_tree());
	for( core::Size resi = from_res + 1; resi<to_res; ++resi ) {
//		if (!tree.is_cutpoint(resi) || resi >= (pose.total_residue() - 1))
//			continue;
		if( pose.residue( resi ).is_lower_terminus() ){
			core::conformation::idealize_position( resi, pose.conformation() );
			remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, resi );
			remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, resi + 1 );
			TR<<"removing cut at "<<resi<<std::endl;
		}
	}
}
*/

/// @brief copy a stretch of aligned phi-psi dofs from source to target. No repacking no nothing.
void
copy_stretch( core::pose::Pose & target, core::pose::Pose const & source, core::Size const from_res, core::Size const to_res ){
	using namespace core::pose;
	using namespace protocols::rosetta_scripts;
	using namespace core::chemical;
	core::Size const from_nearest_on_source( find_nearest_res( source, target, from_res ) );
	core::Size const to_nearest_on_source( find_nearest_res( source, target, to_res ) );
TR<<"target: "<<from_res<<" "<<to_res<<" source: "<<from_nearest_on_source<<" "<<to_nearest_on_source<<std::endl;
	runtime_assert( from_nearest_on_source && to_nearest_on_source );
	core::Size const residue_diff( to_nearest_on_source - from_nearest_on_source - (to_res - from_res ));
	protocols::protein_interface_design::movers::LoopLengthChange llc;
	llc.loop_start( from_res );
	llc.loop_end( to_res );
	llc.delta( residue_diff );
	llc.apply( target );

	target.copy_segment( to_nearest_on_source - from_nearest_on_source + 1, source, from_res, from_nearest_on_source );
	target.update_residue_neighbors();
}

void
Splice::apply( core::pose::Pose & pose )
{
	using namespace protocols::rosetta_scripts;
	using core::chemical::DISULFIDE;

	TR<<"Starting splice apply"<<std::endl;
	save_values();

/// from_res() and to_res() can be determined directly on the tag, through a taskfactory, or through a template file. If through a template file,
/// we start by translating from_res/to_res from the template file to the in coming pose as in the following paragraph
	if( template_file_ != "" ){ /// using a template file to determine from_res() to_res()
		core::Size template_from_res( 0 ), template_to_res( 0 );

		if( from_res() && to_res() ){
			template_from_res = find_nearest_res( pose, *template_pose_, from_res());
			template_to_res   = find_nearest_res( pose, *template_pose_, to_res()  );
			runtime_assert( template_from_res );
			runtime_assert( template_to_res );
		}

		from_res( template_from_res );
		to_res( template_to_res );
	}// fi template_file != ""

	core::pose::Pose const in_pose_copy( pose );
	pose.conformation().detect_disulfides(); // just in case; but I think it's unnecessary


/// from_res/to_res can also be determined through task factory, by identifying the first and last residues that are allowed to design in this tf
	if( torsion_database_fname_ == "" && from_res() == 0 && to_res() == 0 ){/// set the splice site dynamically according to the task factory
		utility::vector1< core::Size > designable( protocols::rosetta_scripts::residue_packer_states( pose, task_factory(), true/*designable*/, false/*packable*/ ) );
		std::sort( designable.begin(), designable.end() );
		from_res( designable[ 1 ] );
		to_res( designable[ designable.size() ] );
	}
	core::pose::Pose source_pose;
	core::Size nearest_to_from( 0 ), nearest_to_to( 0 ), residue_diff( 0 ); // residues on source_pose that are nearest to from_res and to_res; what is the difference in residue numbers between incoming pose and source pose
	ResidueBBDofs dofs; /// used to store the torsion/resid dofs from any of the input files
	dofs.clear();
	core::Size cut_site( 0 );
	if( torsion_database_fname_ == "" ){ // read dofs from source pose rather than database
		core::import_pose::pose_from_pdb( source_pose, source_pdb_ );
		nearest_to_from = find_nearest_res( source_pose, pose, from_res() );
		nearest_to_to = find_nearest_res( source_pose, pose, to_res() );
		residue_diff = nearest_to_to - nearest_to_from - ( to_res() - from_res() );
		if( nearest_to_from == 0 || nearest_to_to == 0 ){
			TR<<"nearest_to_from: "<<nearest_to_from<<" nearest_to_to: "<<nearest_to_to<<". Failing"<<std::endl;
			set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
			retrieve_values();
			return;
		}
		for( core::Size i = nearest_to_from; i <= nearest_to_to; ++i ){
		  if( source_pose.residue( i ).has_variant_type( DISULFIDE ) ){/// in future, using disulfides would be a great boon as it rigidifies loops.
				TR<<"Residue "<<i<<" is a disulfide. Failing"<<std::endl;
				set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
				retrieve_values();
				return;
			}
/// Feed the source_pose dofs into the BBDofs array
			BBDofs residue_dofs;
			residue_dofs.resid( i ); /// resid is probably never used
			residue_dofs.phi( source_pose.phi( i ) );
			residue_dofs.psi( source_pose.psi( i ) );
			residue_dofs.omega( source_pose.omega( i ) );

			core::Size const nearest_on_target( find_nearest_res( pose, source_pose, i ) );

/// convert 3let residue code to 1let code
			std::stringstream ss; std::string s;
			ss << source_pose.residue( i ).name1();
			ss >> s;
			residue_dofs.resn( s );

			dofs.push_back( residue_dofs );
		}// for i nearest_to_from..nearest_to_to
		cut_site = dofs.cut_site() ? dofs.cut_site() + from_res() - 1: to_res(); // isn't this always going to be to_res()? I think so...
	}// fi torsion_database_fname==NULL
	else{/// read from dbase
		core::Size dbase_entry( database_entry() );
		if( dbase_entry == 0 ){/// randomize entry
			core::Size rand_trials( 0 );
			core::Size pose_residues( 0 );/// number of residues in loop on incoming pose
			do{/// choose random dbase entry. If equal_length is on then make sure the random entry has the same length
				dbase_entry = (core::Size) ( RG.uniform() * torsion_database_.size() + 1);
				dofs = torsion_database_[ dbase_entry ];
				rand_trials++;
				core::Size const nearest_to_entry_start_on_pose( find_nearest_res( pose, *template_pose_, dofs.start_loop() ) );
				core::Size const nearest_to_entry_stop_on_pose( find_nearest_res( pose, *template_pose_, dofs.stop_loop() ) );
			  pose_residues = nearest_to_entry_stop_on_pose - nearest_to_entry_start_on_pose + 1;
			} while( equal_length() && dofs.size() != pose_residues && rand_trials < torsion_database_.size() * 10/*prevent infinite loops*/ );
			if( rand_trials >=  torsion_database_.size() * 10 ){
				TR<<"Loop of appropriate length not found in database. Returning"<<std::endl;
				retrieve_values();
				return;
			}
			else
				TR<<"Randomly chose dbase entry "<<dbase_entry<<" with loop start: "<<dofs.start_loop()<<" loop_stop: "<<dofs.stop_loop()<<" cut: "<<dofs.cut_site()<<" total_residues: "<<dofs.size()<<std::endl;
		}
		else /// don't randomize entry, pick the one stated on the tag.
			dofs = torsion_database_[ dbase_entry ];
		foreach( BBDofs & resdofs, dofs ){/// transform 3-letter code to 1-letter code
			using namespace core::chemical;
			if( resdofs.resn() == "CYD" ){// at one point it would be a good idea to use disfulfides rather than bail out on them...; I think disulfided cysteins wouldn't be written as CYD. This requires something more clever...
				TR<<"Residue "<<resdofs.resid()<<" is a disulfide. Failing"<<std::endl;
				set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
				retrieve_values();
				return;
			}
			std::stringstream ss; std::string s;
			ss << oneletter_code_from_aa( aa_from_name( resdofs.resn() ) );
			ss >> s;
			resdofs.resn( s );
		}/// foreach resdof
		nearest_to_to = dofs.size(); /// nearest_to_to and nearest_to_from are used below to compute the difference in residue numbers...
		nearest_to_from = 1;
 /// set from_res/to_res/cut_site on the incoming pose
		if( template_file_ != "" ){/// according to the template pose
			from_res( find_nearest_res( pose, *template_pose_, dofs.start_loop() ) );
			to_res( find_nearest_res( pose, *template_pose_, dofs.stop_loop() ) );
//			to_res( from_res() + dofs.size() -1);
			runtime_assert( from_res() );
			runtime_assert( to_res() );
			cut_site = dofs.cut_site() - dofs.start_loop() + from_res();
		} // fi template_file != ""
		else{/// according to the dofs array (taken from the dbase)
			from_res( dofs.start_loop() );
			to_res( dofs.stop_loop() );
			cut_site = dofs.cut_site();
			runtime_assert( from_res() && to_res() && cut_site );
		}
		residue_diff = dofs.size() - ( dofs.stop_loop() - dofs.start_loop()  + 1 );
  }// read from dbase
	TR<<"From res: "<<from_res()<<" to_res: "<<to_res()<<std::endl;
	runtime_assert( to_res() > from_res() );
	if( saved_fold_tree_ )/// is saved_fold_tree_ being used?
		pose.fold_tree( *saved_fold_tree_ );

/// The database is computed with respect to the template pose, so before applying dofs from the dbase it's important to make that stretch identical to
/// the template. from_res() and to_res() were previously computed to be with respect to the incoming pose, so within this subroutine the refer to pose rather
/// than template_pose (this is a bit confusing, but it works!)
	copy_stretch( pose, *template_pose_, from_res(), to_res() );

	using namespace utility;
/// randomize_cut() should not be invoked with a database entry, b/c the dbase already specified the cut sites.
/// this is important b/c nearest_to_from/nearest_to_to are degenerate if the dbase is used.
	if( randomize_cut() ){
/// choose cutsite randomly within loop residues on the loop (no 2ary structure)
		core::scoring::dssp::Dssp dssp( source_pose );
		dssp.dssp_reduced(); // switch to simplified H E L notation
		std::vector< core::Size > loop_positions_in_source;
		loop_positions_in_source.clear();
		TR<<"DSSP of source segment: ";
		for( core::Size i = nearest_to_from; i <= std::min( nearest_to_to, to_res() - from_res() + nearest_to_from ); ++i ){
			if( dssp.get_dssp_secstruct( i ) == 'L' )
				loop_positions_in_source.push_back( i );
			TR<<dssp.get_dssp_secstruct( i );
		}
		TR<<std::endl;
		cut_site = loop_positions_in_source[ (core::Size) ( RG.uniform() * loop_positions_in_source.size()) ] - nearest_to_from + from_res();
		TR<<"Cut placed at: "<<cut_site<<std::endl;
  }// fi randomize_cut

	fold_tree( pose, from_res(), to_res(), cut_site );/// the fold_tree routine will actually set the fold tree to surround the loop
/// change the loop length
	protocols::protein_interface_design::movers::LoopLengthChange llc;
	llc.loop_start( from_res() );
	llc.loop_end( cut_site + residue_diff < from_res() ? to_res() : cut_site );
	llc.delta( residue_diff );
	llc.apply( pose );

/// set torsions
	core::Size const total_residue_new( dofs.size() );
	TR<<"Changing dofs for resi: ";
	for( core::Size i = 0; i < total_residue_new; ++i ){
		core::Size const pose_resi( from_res() + i );
		pose.set_phi( pose_resi, dofs[ i + 1 ].phi() );
		pose.set_psi( pose_resi, dofs[ i + 1 ].psi() );
		pose.set_omega( pose_resi, dofs[ i + 1 ].omega() );
		TR<<pose_resi<<",";
//		char f[ 10];
//		sprintf( f, "b%d.pdb", i );
	}

	pose.conformation().detect_disulfides();/// probably unnecessary
	TR<<std::endl;
	std::string threaded_seq( "" );/// will be all ALA except for Pro/Gly on source pose and matching identities on source pose
/// Now decide on residue identities: Alanine throughout except when the template pose has Gly, Pro or a residue that is the same as that in the original pose
	for( core::Size i = 0; i < total_residue_new; ++i ){
		core::Size const pose_resi( from_res() + i );
		std::string const dofs_resn( dofs[ i + 1 ].resn() );
		runtime_assert( dofs_resn.length() == 1 );
		if( design() ){ // all non pro/gly residues in template are allowed to design
			if( dofs_resn == "G" || dofs_resn == "P" )
				threaded_seq += dofs_resn;
			else
				threaded_seq += 'x';
			continue;
		}
		core::Size const nearest_in_copy( find_nearest_res( in_pose_copy, pose, pose_resi ) );
		if( ( nearest_in_copy > 0 && dofs_resn[ 0 ] == in_pose_copy.residue( nearest_in_copy ).name1() )  || dofs_resn == "G" || dofs_resn == "P" )
			threaded_seq += dofs_resn;
		else{
			if( poly_ala() )
				threaded_seq += "A";
			else{
				char orig_residue( 0 );
				if( nearest_in_copy )
					orig_residue = in_pose_copy.residue( nearest_in_copy ).name1();
				if( orig_residue == 0 || orig_residue == 'G' || orig_residue == 'P' )
					threaded_seq += 'x'; // residues that were originally Gly/Pro can be designed now
				else
					threaded_seq += ' '; // only repack
			}
		}
	}

	using namespace protocols::toolbox::task_operations;
	using namespace core::pack::task;
	ThreadSequenceOperationOP tso = new ThreadSequenceOperation;
	tso->target_sequence( threaded_seq );
	tso->start_res( from_res() );
	tso->allow_design_around( false );
	TR<<"Threading sequence: "<<threaded_seq<<" starting from "<<from_res()<<std::endl;
	TaskFactoryOP tf = new TaskFactory;
	tf->push_back( new operation::InitializeFromCommandline );
	tf->push_back( new operation::NoRepackDisulfides );
	tf->push_back( tso );
	DesignAroundOperationOP dao = new DesignAroundOperation;
	dao->design_shell( 1.0 ); // threaded sequence operation needs to design, and will restrict design to the loop only
	dao->repack_shell( 8.0 );
	for( core::Size i = from_res() - 1; i <= from_res() + total_residue_new + 1; ++i ){
		if( !pose.residue( i ).has_variant_type( DISULFIDE ) )
			dao->include_residue( i );
	}
	tf->push_back( dao );
	operation::RestrictAbsentCanonicalAASOP racaas = new operation::RestrictAbsentCanonicalAAS;
	racaas->keep_aas( "ADEFIKLMNQRSTVWY" ); /// disallow pro/gly/cys/his
	racaas->include_residue( 0 ); /// restrict all residues
	tf->push_back( racaas);

	protocols::protein_interface_design::movers::AddChainBreak acb;
	acb.resnum( utility::to_string( cut_site + residue_diff ) );
	acb.find_automatically( false );
	acb.change_foldtree( false );
	acb.apply( pose );
	TR<<"Adding chainbreak at: "<<cut_site<<std::endl;

	if( ccd() ){
		using namespace protocols::loops;
		Loop loop( std::max( (core::Size) 2, from_res() - 6 )/*start*/, std::min( pose.total_residue()-1, to_res() + 6 )/*stop*/, cut_site/*cut*/ );
		LoopsOP loops = new Loops();
		loops->push_back( loop );

/// Set ccd to minimize 4 residues at each loop terminus including the first residue of the loop. This way,
/// the torsion in the loop are maintained. Allow repacking around the loop.
/// If disulfide occurs in the range that is allowed to minimize, adjust that region to not include disulf
		core::scoring::ScoreFunctionOP scorefxn_local( scorefxn()->clone() );/// in case you want to modify the scorefxn. Currently not used
		protocols::loops::loop_mover::refine::LoopMover_Refine_CCD ccd_mover( loops, scorefxn_local );
		ccd_mover.temp_initial( 1.5 );
		ccd_mover.temp_final( 0.5 );
		core::kinematics::MoveMapOP mm;
		mm = new core::kinematics::MoveMap;
		mm->set_chi( false ); mm->set_bb( false ); mm->set_jump( false );
	/// First look for disulfides. Those should never be moved.
		core::Size disulfn( 0 ), disulfc( 0 );
		for( core::Size i = from_res() - 3; i <= from_res(); ++i ){
			if( pose.residue( i ).has_variant_type( DISULFIDE ) ){
				disulfn = i;
			}
		}
		for( core::Size i = from_res() + total_residue_new - 1; i <= from_res() + total_residue_new + 2; ++i ){
			if( pose.residue( i ).has_variant_type( DISULFIDE ) ){
				disulfc = i;
				break;
			}
		}
		core::Size const startn( disulfn > 0 ? disulfn + 1 : from_res() - 3 );
		core::Size const startc( disulfc > 0 ? disulfc - 6 : from_res() + total_residue_new - ( res_move() - 3 ) );
		for( core::Size i = startn; i <= startn + res_move() - 1; ++i ){
			mm->set_chi( i, true );
			mm->set_bb( i, true );
		}
		for( core::Size i = startc; i <= startc + res_move() - 1; ++i ){
			mm->set_chi( i, true );
			mm->set_bb( i, true );
		}
		ccd_mover.set_task_factory( tf );
		ccd_mover.move_map( mm );
		ccd_mover.apply( pose );

/// following ccd, compute rmsd to source loop to ensure that you haven't moved too much. This is pretty decent filter
		if( torsion_database_fname_ == "" ){ // no use computing rms if coming from a database (no coordinates)
			core::Real rms( 0 );
			for( core::Size i = 0; i <= total_residue_new - 1; ++i ){
				core::Real const dist( pose.residue( from_res() + i ).xyz( "CA" ).distance( source_pose.residue( nearest_to_from+ i ).xyz("CA" ) ) );
				rms += dist;
			}
			core::Real const average_rms( rms / total_residue_new );
			TR<<"Average distance of spliced segment to original: "<< average_rms<<std::endl;
			if( average_rms >= rms_cutoff() ){
				TR<<"Failing."<<std::endl;
				set_last_move_status( protocols::moves::FAIL_RETRY );
				retrieve_values();
				return;
			}
		}
/// tell us what the torsions of the new (closed) loop are. This is used for dbase construction. At one point, might be a good idea to make the mover
/// output the dofs directly to a dbase file rather than to a log file.
		TaskFactoryOP tf_dofs = new TaskFactory;
		DesignAroundOperationOP dao_dofs = new DesignAroundOperation;
		for( core::Size i = startn; i <= startc + res_move() - 1; ++i )
			dao_dofs->include_residue( i );
		dao_dofs->design_shell( 0 );/// only include the loop residues
		tf_dofs->push_back( dao_dofs );
		protocols::protein_interface_design::filters::Torsion torsion;
		torsion.task_factory( tf_dofs );
		torsion.task_factory_set( true );
		torsion.apply( pose );
		core::Size const stop_on_template( startc + res_move() - 1 - residue_diff );
		TR_ccd << "start, stop, cut: "<<startn<<" "<<stop_on_template<<" "<<cut_site<<std::endl; /// used for the dbase
	}// fi ccd
	else{ // if no ccd, still need to thread sequence
		PackerTaskOP ptask = tf()->create_task_and_apply_taskoperations( pose );
		protocols::simple_moves::PackRotamersMover prm( scorefxn(), ptask );

		prm.apply( pose );
	}
	saved_fold_tree_ = new core::kinematics::FoldTree( pose.fold_tree() );
	retrieve_values();
}

/// splice apply might change the from_res/to_res internals since they sometimes refer to the template file. If that happens, we want the values to
/// revert to their original values before the end of the apply function (so retrieve_values) below must be called before return.
void
Splice::save_values(){
	saved_from_res_ = from_res();
	saved_to_res_ = to_res();
}

void
Splice::retrieve_values(){
	from_res( saved_from_res_ );
	to_res( saved_to_res_ );
}

std::string
Splice::get_name() const {
	return SpliceCreator::mover_name();
}

void
Splice::parse_my_tag( TagPtr const tag, protocols::moves::DataMap &data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose )
{
	start_pose_ = new core::pose::Pose( pose );
	runtime_assert( tag->hasOption( "task_operations" ) != (tag->hasOption( "from_res" ) || tag->hasOption( "to_res" ) ) || tag->hasOption( "torsion_database" ) ); // it makes no sense to activate both taskoperations and from_res/to_res.
	runtime_assert( tag->hasOption( "torsion_database" ) != tag->hasOption( "source_pdb" ) );
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	if( !tag->hasOption( "task_operations" ) ){
		from_res( protocols::rosetta_scripts::parse_resnum( tag->getOption< std::string >( "from_res", "0" ), pose ) );
		to_res( protocols::rosetta_scripts::parse_resnum( tag->getOption< std::string >( "to_res", "0" ), pose ) );
	}
	if( tag->hasOption( "torsion_database" ) ){
		torsion_database_fname( tag->getOption< std::string >( "torsion_database" ) );
		database_entry( tag->getOption< core::Size >( "database_entry", 0 ) );
		read_torsion_database();
		TR<<"torsion_database: "<<torsion_database_fname()<<" ";
		if( database_entry() == 0 )
			TR<<" database entry will be randomly picked at run time. ";
		else{
			TR<<" database_entry: "<<database_entry()<<" ";
			runtime_assert( database_entry() <= torsion_database_.size() );
		}
	}
	else
		source_pdb( tag->getOption< std::string >( "source_pdb" ) );
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	ccd( tag->getOption< bool >( "ccd", 1 ) );
	rms_cutoff( tag->getOption< core::Real >( "rms_cutoff", 999999 ) );
	runtime_assert( !(tag->hasOption( "torsion_database" ) && tag->hasOption( "rms_cutoff" )) ); // torsion database doesn't specify coordinates so no point in computing rms
	res_move( tag->getOption< core::Size >( "res_move", 4 ) );
	randomize_cut( tag->getOption< bool >( "randomize_cut", false ) );
	runtime_assert( tag->hasOption( "randomize_cut" ) && tag->hasOption( "source_pose" ) || !tag->hasOption( "source_pose" ) );
	template_file( tag->getOption< std::string >( "template_file", "" ) );
	equal_length( tag->getOption< bool >( "equal_length", false ) );
	poly_ala( tag->getOption< bool >( "thread_ala", true ) );

	if( template_file_ != "" ){ /// using a template file to determine from_res() to_res()
		if( data.has( "poses", template_file_ ) ){
			template_pose_ = data.get< core::pose::Pose * >( "poses", template_file_ );
			TR<<"using template pdb from datamap"<<std::endl;
		}
		else if( tag->hasOption( "template_file" ) ){
			template_pose_ = new core::pose::Pose;
			core::import_pose::pose_from_pdb( *template_pose_, template_file_ );
			data.add( "poses", template_file_, template_pose_ );
			TR<<"loading template_pose from "<<template_file_<<std::endl;
		}
	}
	else
		template_pose_ = new core::pose::Pose( pose );

	design( tag->getOption< bool >( "design", false ) );
	TR<<"from_res: "<<from_res()<<" to_res: "<<to_res()<<" randomize_cut: "<<randomize_cut()<<" source_pdb: "<<source_pdb()<<" ccd: "<<ccd()<<" rms_cutoff: "<<rms_cutoff()<<" res_move: "<<res_move()<<" template_file: "<<template_file()<<std::endl;
}

protocols::moves::MoverOP
Splice::clone() const {
    return( protocols::moves::MoverOP( new Splice( *this ) ));
}

void
Splice::scorefxn( core::scoring::ScoreFunctionOP sf ){
	scorefxn_ = sf;
}

core::scoring::ScoreFunctionOP
Splice::scorefxn() const{
	return scorefxn_;
}

core::pack::task::TaskFactoryOP
Splice::task_factory() const{ return task_factory_; }

void
Splice::task_factory( core::pack::task::TaskFactoryOP tf ){ task_factory_ = tf; }


/// the torsion dbase should have the following structure:
/// each line represents a single loop. Each four values represent <phi> <psi> <omega> <3-let resid>; the last entry in a line represents <loop start> <loop stop> <cut site> cut; where cut signifies that this is the loop designator
void
Splice::read_torsion_database(){
	using namespace std;

  utility::io::izstream data( torsion_database_fname_ );
  if ( !data ) {
    TR << "cannot open torsion database " << torsion_database_fname_ << std::endl;
		set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
		return;
  }
  std::string line;
  while( getline( data, line ) ) {
    std::istringstream line_stream( line );
		ResidueBBDofs bbdof_entry;
		bbdof_entry.clear();
  	while( !line_stream.eof() ){
			core::Real phi, psi, omega;
			std::string resn;
			line_stream >> phi >> psi >> omega >> resn;
			if( resn == "cut" || resn == "CUT" ){ /// ugly!!!
				bbdof_entry.start_loop( (core::Size ) phi );
				bbdof_entry.stop_loop( (core::Size ) psi );
				bbdof_entry.cut_site( (core::Size ) omega );
			}
			else
				bbdof_entry.push_back( BBDofs( 0/*resid*/, phi, psi, omega, resn ) ); /// resid may one day be used. Currently it isn't
		}
		torsion_database_.push_back( bbdof_entry );
	}
	TR<<"Finished reading torsion database with "<<torsion_database_.size()<<" entries"<<std::endl;
}

///@brief set the fold tree around start/stop/cut sites.
/// presently makes a simple fold tree, but at one point may be a more complicated function to include two poses
void
Splice::fold_tree( core::pose::Pose & pose, core::Size const start, core::Size const stop, core::Size const cut ) const{
	using namespace protocols::loops;
	core::conformation::Conformation const & conf( pose.conformation() );
	core::Size const s1 = std::max( (core::Size) 2, start - 6 );
	core::Size const s2 = std::min( conf.chain_end( 1 ) - 1, stop + 6 );
	if( conf.num_chains() == 1 ){
		FoldTreeFromLoops ffl;
		Loop loop( s1, s2, cut/*cut*/ );
		LoopsOP loops = new Loops();
		loops->push_back( loop );
		ffl.loops( loops );
		ffl.apply( pose );
		return;
	}
	core::Size from_res( 0 ), to_res( 0 );
	for( core::Size resi = conf.chain_begin( 1 ); resi <= conf.chain_end( 1 ); ++resi ){
		 if( pose.residue( resi ).has_variant_type( core::chemical::DISULFIDE ) ){
			from_res = resi;
			break;
		}
	}
	to_res = protocols::geometry::residue_center_of_mass( pose, conf.chain_begin( 2 ), conf.chain_end( 2 ) );
	core::kinematics::FoldTree ft;
	ft.clear();
	core::Size const first = std::min( s1, from_res );
	core::Size const second = std::max( s1, from_res );
//	ft.add_edge( 1, from_res, -1 );
//	ft.add_edge( from_res, s1, -1 );
	ft.add_edge( 1, s1, -1 );
	ft.add_edge( s1, cut, -1 );
	ft.add_edge( s2, cut + 1, -1 );
	ft.add_edge( s1, s2, 1 );
	ft.add_edge( s2, conf.chain_end( 1 ), -1 );
	ft.add_edge( 1, conf.chain_begin( 2 ), 2 );
	ft.add_edge( conf.chain_begin( 2 ), conf.chain_end( 2 ), -1 );
//	ft.add_edge( from_res, to_res, 2 );
//	ft.add_edge( to_res, conf.chain_begin( 2 ), -1 );
//	ft.add_edge( to_res, conf.chain_end( 2 ), -1 );
	ft.reorder(1);
	TR<<"Previous ft: "<<pose.fold_tree()<<std::endl;
	pose.fold_tree( ft );
	TR<<"Current ft: "<<pose.fold_tree()<<std::endl;
}

BBDofs::~BBDofs() {}

ResidueBBDofs::~ResidueBBDofs() {}

} //movers
} //protein_interface_design
} //protocols
