// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/fldsgn/BluePrintBDR.cc
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu ),  Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/fldsgn/BluePrintBDR.hh>
#include <protocols/fldsgn/BluePrintBDRCreator.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentInsert.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/methods/pose_mod.hh>

// project headers
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <protocols/denovo_design/constraints/FileConstraintGenerator.hh>
#include <protocols/forge/constraints/NtoCConstraintGenerator.hh>
#include <protocols/forge/constraints/InvrotTreeRCG.hh>
#include <protocols/fldsgn/SheetConstraintGenerator.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/toolbox/match_enzdes_util/AlignPoseToInvrotTreeMover.hh>
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTree.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/tag/Tag.hh>

// C++ headers
#include <utility>

// boost

#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace fldsgn {


static THREAD_LOCAL basic::Tracer TR( "protocols.fldsgn.BluePrintBDR" );

std::string
BluePrintBDRCreator::keyname() const
{
	return BluePrintBDRCreator::mover_name();
}

protocols::moves::MoverOP
BluePrintBDRCreator::create_mover() const {
	return protocols::moves::MoverOP( new BluePrintBDR );
}

std::string
BluePrintBDRCreator::mover_name()
{
	return "BluePrintBDR";
}


/// @brief default constructor
BluePrintBDR::BluePrintBDR() :
	Super( "BluePrintBDR" ),
	blueprint_( /* NULL */ ),
	sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "fldsgn_cen" ) ),
	loop_mover_str_( "RemodelLoopMover" ),
	use_fullmer_( false ),
	num_fragpick_( 200 ),
	use_sequence_bias_( false ),
	use_abego_bias_( false ),
	max_linear_chainbreak_( 0.07 ),
	initialized_( false ),
	ss_from_blueprint_( true ),
	constraints_NtoC_( -1.0 ),
	constraints_sheet_( -1.0 ),
	constraint_file_( "" ),
	dump_pdb_when_fail_( "" ),
	rmdl_attempts_( 1 ),
	use_poly_val_( true ),
	tell_vlb_to_not_touch_fold_tree_(false),
	invrot_tree_(/* NULL */),
	enzcst_io_(/* NULL */)
{
	rcgs_.clear();
}

/// @brief value constructor
BluePrintBDR::BluePrintBDR( String const & filename, bool const ss_from_blueprint ) :
	Super( "BluePrintBDR" ),
	sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "fldsgn_cen" ) ),
	loop_mover_str_( "RemodelLoopMover" ),
	use_fullmer_( false ),
	num_fragpick_( 200 ),
	use_sequence_bias_( false ),
	use_abego_bias_( false ),
	max_linear_chainbreak_( 0.07 ),
	initialized_( false ),
	ss_from_blueprint_( ss_from_blueprint ),
	constraints_NtoC_( -1.0 ),
	constraints_sheet_( -1.0 ),
	constraint_file_( "" ),
	dump_pdb_when_fail_( "" ),
	rmdl_attempts_( 1 ),
	use_poly_val_( true ),
	tell_vlb_to_not_touch_fold_tree_(false),
	invrot_tree_(/* NULL */),
	enzcst_io_(/* NULL */)
{
	set_blueprint( filename );
	rcgs_.clear();
}

/// @brief value constructor
BluePrintBDR::BluePrintBDR( BluePrintOP const & blueprintOP, bool const ss_from_blueprint ) :
	Super( "BluePrintBDR" ),
	blueprint_( blueprintOP ),
	sfx_( core::scoring::ScoreFunctionFactory::create_score_function( "fldsgn_cen" ) ),
	loop_mover_str_( "RemodelLoopMover" ),
	use_fullmer_( false ),
	num_fragpick_( 200 ),
	use_sequence_bias_( false ),
	use_abego_bias_( false ),
	max_linear_chainbreak_( 0.07 ),
	initialized_( false ),
	ss_from_blueprint_( ss_from_blueprint ),
	constraints_NtoC_( -1.0 ),
	constraints_sheet_( -1.0 ),
	constraint_file_( "" ),
	dump_pdb_when_fail_( "" ),
	rmdl_attempts_( 1 ),
	use_poly_val_( true ),
	tell_vlb_to_not_touch_fold_tree_(false),
	invrot_tree_(/* NULL */),
	enzcst_io_(/* NULL */)
{
	rcgs_.clear();
}

/// @Brief copy constructor
BluePrintBDR::BluePrintBDR( BluePrintBDR const & rval ) :
	//utility::pointer::ReferenceCount(),
	Super( rval ),
	blueprint_( rval.blueprint_ ),
	manager_( rval.manager_ ),
	sfx_( rval.sfx_ ),
	loop_mover_str_( rval.loop_mover_str_ ),
	use_fullmer_( rval.use_fullmer_ ),
	num_fragpick_( rval.num_fragpick_ ),
	use_sequence_bias_( rval.use_sequence_bias_ ),
	use_abego_bias_( rval.use_abego_bias_ ),
	max_linear_chainbreak_( rval.max_linear_chainbreak_ ),
	initialized_( rval.initialized_ ),
	ss_from_blueprint_( rval.ss_from_blueprint_ ),
	constraints_NtoC_( rval.constraints_NtoC_ ),
	constraints_sheet_( rval.constraints_sheet_ ),
	constraint_file_( rval.constraint_file_ ),
	dump_pdb_when_fail_( rval.dump_pdb_when_fail_ ),
	rmdl_attempts_( rval.rmdl_attempts_ ),
	use_poly_val_( rval.use_poly_val_ ),
	tell_vlb_to_not_touch_fold_tree_( rval.tell_vlb_to_not_touch_fold_tree_),
	invrot_tree_(rval.invrot_tree_),
	enzcst_io_(rval.enzcst_io_),
	rcgs_( rval.rcgs_ )
{
	if ( rval.vlb_.get() ) {
		vlb_ = VarLengthBuildOP( new VarLengthBuild( *rval.vlb_ ) );
	}
}

/// @brief default destructor
BluePrintBDR::~BluePrintBDR() = default;

/// @brief clone this object
BluePrintBDR::MoverOP
BluePrintBDR::clone() const
{
	return BluePrintBDR::MoverOP( new BluePrintBDR( *this ) );
}

/// @brief create this type of object
BluePrintBDR::MoverOP
BluePrintBDR::fresh_instance() const
{
	return BluePrintBDR::MoverOP( new BluePrintBDR() );
}

/// @brief the centroid level score function, default "remodel_cen"
BluePrintBDR::ScoreFunction const &
BluePrintBDR::scorefunction() const
{
	return *sfx_;
}

/// @brief add instruction to the manager of this BluePrintBDR (no copy)
/// @param[in] bi BuildInstruction
void
BluePrintBDR::add_instruction( BuildInstructionOP bi )
{
	manager_.add( bi );

	// additional instruction means we'll need a new re-init the VLB, so
	// go ahead and drop the existing one
	vlb_.reset();
}


/// @brief create directed dependency between two instructions
void
BluePrintBDR::create_directed_dependency(
	BuildInstructionOP u,
	BuildInstructionOP v
)
{
	manager_.create_directed_dependency( u, v );
}


/// @brief set the centroid level score function
void
BluePrintBDR::scorefunction( ScoreFunction const & sfx )
{
	sfx_ = sfx.clone();
}

/// @brief set the centroid level score function
void
BluePrintBDR::scorefunction( ScoreFunctionOP sfx )
{
	sfx_ = sfx->clone();
}

/// @brief use blueprint
void
BluePrintBDR::set_blueprint( String const & filename )
{
	blueprint_ = BluePrintOP( new BluePrint( filename ) );
}

/// @brief use blueprint
void
BluePrintBDR::set_blueprint( BluePrintOP const & blp )
{
	blueprint_ = blp;
}

/// @brief set constraint between N- and C- terminal residues, need to set the weight
void
BluePrintBDR::set_constraints_NtoC( Real const & weight )
{
	constraints_NtoC_ = weight;
}

/// @brief set constraint file
void
BluePrintBDR::set_constraint_file( String const & constraint_file )
{
	constraint_file_ = constraint_file;
}

/// @brief set list of remodel constraint generators
void
BluePrintBDR::set_rcgs( utility::vector1< protocols::forge::remodel::RemodelConstraintGeneratorOP > const & rcgs )
{
	rcgs_ = rcgs;
}

/// @brief dump pdb when this protocol failed
void
BluePrintBDR::dump_pdb_when_fail( String const & dump_pdb_when_fail )
{
	dump_pdb_when_fail_ = dump_pdb_when_fail;
}

/// @brief set instruction by blueprint
bool
BluePrintBDR::set_instruction_blueprint( Pose const & pose )
{
	using protocols::forge::build::Interval;
	using protocols::forge::build::SegmentRebuild;
	using protocols::forge::build::SegmentInsert;

	bool flag( false ), insert( false );
	String aa, ss, insert_name, ins_sec;
	Size left( 0 ), right( 0 ), count( 0 ), insnum( 0 );
	Pose insert_pose;
	for ( Size i=1; i<=blueprint_->total_residue(); i++ ) {

		if ( blueprint_->resnum( i ) != 0 ) {
			count++;
			if ( count > pose.size() ) {
				TR.Error << "Residue number in blueprint file is more than that of pose!,  pose/blueprint= "
					<< pose.size() << "/" << count << std::endl;
				return false;
			}
		}

		//if( blueprint_->buildtype( i ) != '.' && !flag ) {  // found R or X at first
		if ( blueprint_->buildtype( i ) == 'R' && blueprint_->buildtype( i ) != 'I' && !flag ) {  // found R at first

			if ( count == 0 || i==1 ) { // N-terminal extensione
				left = 1;
			} else {
				left = blueprint_->resnum( i-1 )+1; // insert in the middle of sequence
				if ( left > pose.size() ) {  // C-terminal extension
					left = pose.size();
				}
			}
			flag = true;

			// } else if( blueprint_->buildtype( i ) == '.' && flag ){ // add instruction
		} else if ( blueprint_->buildtype( i ) != 'R' && blueprint_->buildtype( i ) != 'I' && flag ) { // add instruction

			right = blueprint_->resnum( i )-1;
			if ( right == 0 ) { // at the N-termianl
				right = 1;
			} else if ( right < left ) { // insert residues into the successive residue numbers
				right = left;
			}
			runtime_assert( right <= pose.size() );

			if ( insert ) {
				if ( ins_sec.length() == 1 ) {
					TR << "Secondary structure of insert pose will be given by Dssp" << std::endl;
					Dssp dssp( pose );
					dssp.insert_ss_into_pose( insert_pose );
				} else {
					runtime_assert( insert_pose.size() == ins_sec.length() );
					for ( Size j=1; j<=insert_pose.size(); j++ ) {
						insert_pose.set_secstruct( j, ins_sec[ j-1 ] );
					}
				}
				TR << "SegmentInsert left " << left << ", right: " << right
					<< ", ss: " << ss << ", aa:" << aa << ", pdb:" << insert_name << std::endl;
				add_instruction( BuildInstructionOP( new SegmentInsert( Interval( left, right ), ss, aa, insert_pose ) ) );
			} else {
				TR << "SegmentRebuild left: " << left << ", right: " << right << ", ss: " << ss << ", aa:" << aa << std::endl;
				add_instruction( BuildInstructionOP( new SegmentRebuild( Interval( left, right ), ss, aa ) ) );
			}

			flag = false;
			insert = false;
			aa = "";
			ss = "";

		} //blurprint_->buildtype()

		if ( flag ) {
			if ( blueprint_->buildtype( i ) == 'I' && !insert ) {
				insert = true;
				insnum ++;
				ins_sec = blueprint_->secstruct( i );
				insert_name = blueprint_->insertion( insnum );
				core::import_pose::pose_from_file( insert_pose, insert_name , core::import_pose::PDB_file);
				aa += '^';
				ss += '^';
			} else if ( blueprint_->buildtype( i ) == 'I' ) {
				ins_sec += blueprint_->secstruct( i );
			} else {
				aa += blueprint_->sequence( i );
				ss += blueprint_->secstruct( i );
			}
		} // flag

	} // blueprint_->size()

	if ( flag ) {
		if ( blueprint_->resnum( blueprint_->total_residue() ) == 0 ) {
			right = pose.size();
		} else {
			right = blueprint_->resnum( blueprint_->total_residue() );
			runtime_assert( right <= pose.size() );
		}
		add_instruction( BuildInstructionOP( new SegmentRebuild( Interval( left, right ), ss, aa ) ) );
		TR << "SegmentRebuild left: " << left << ", right: " << right << ", ss: " << ss << ", aa:" << aa << std::endl;
	}

	if ( instruction_size() < 1 ) {
		TR << "There is no instruction in blueprint. " << std::endl;
		return false;
	}

	return true;
} // set_build_instruction

void
BluePrintBDR::setup_invrot_tree_in_vlb( VarLengthBuild & vlb, Pose & pose  ) const
{

	toolbox::match_enzdes_util::AllowedSeqposForGeomCstOP allowed_seqpos( new protocols::toolbox::match_enzdes_util::AllowedSeqposForGeomCst() );
	//stupid: apparently we have to make a copy of the pose on the heap for
	//the initialization of allowed_seqpos to work
	core::pose::PoseOP posecopy( new core::pose::Pose( pose ) );
	allowed_seqpos->initialize_from_command_line( posecopy ); //this could be moved somewhere else to only initialize once, but probably not that important
	toolbox::match_enzdes_util::AlignPoseToInvrotTreeMoverOP setup_align_pose( new toolbox::match_enzdes_util::AlignPoseToInvrotTreeMover( invrot_tree_, allowed_seqpos ) );
	setup_align_pose->set_add_target_to_pose( true );
	setup_align_pose->set_geomcst_for_superposition_from_enz_io( enzcst_io_);

	toolbox::match_enzdes_util::AlignPoseToInvrotTreeMoverOP run_align_pose( new toolbox::match_enzdes_util::AlignPoseToInvrotTreeMover( invrot_tree_, allowed_seqpos ) );
	run_align_pose->set_geomcst_for_superposition_from_enz_io( enzcst_io_);

	forge::constraints::InvrotTreeRCGOP invrot_rcg( new forge::constraints::InvrotTreeRCG( invrot_tree_, allowed_seqpos ) );
	vlb.add_rcg( invrot_rcg );
	vlb.loop_mover_fold_tree_constant( true ); //we're taking care of the fold tree through the above align movers
	vlb.clear_setup_movers(); //safety
	vlb.add_setup_mover( setup_align_pose );

	vlb.clear_user_provided_movers();
	vlb.add_user_provided_mover( run_align_pose );

	if ( use_abego_bias_ ) {
		utility::vector1< std::string > abego_to_use( blueprint_->abego() );
		abego_to_use.push_back("X"); //kinda hacky, assuming we're only adding one ligand to pose. Ideally we'd query InvrotTree for the number of targets
		vlb_->set_abego( abego_to_use );
	}
}


/// @brief apply defined moves to given Pose
void
BluePrintBDR::apply( Pose & pose )
{
	using protocols::moves::MS_SUCCESS;
	using protocols::moves::FAIL_DO_NOT_RETRY;
	using protocols::moves::FAIL_BAD_INPUT;
	using protocols::moves::FAIL_RETRY;

	//set instruction by blueprint file and store secondary structure information
	if ( !blueprint_ && !initialized_ ) {
		TR << "You need to set a blueprint file" << std::endl;
		set_last_move_status( FAIL_BAD_INPUT );
		return;
	} else if ( !initialized_ ) {
		if ( ! set_instruction_blueprint( pose ) ) {
			set_last_move_status( FAIL_BAD_INPUT );
			return;
		}
		// assign secondary structure   ... probably, this is bug, we have to set SS info everytime
		// if( !ss_from_blueprint_ ){
		// Dssp dssp( pose );
		// dssp.insert_ss_into_pose( pose );
		//  }else{
		// insert_ss_into_pose( pose );
		// }
		//initial_pose_ = new Pose( pose );
		initialized_ = true;
	}

	// assign secondary structure
	if ( !ss_from_blueprint_ ) {
		Dssp dssp( pose );
		dssp.insert_ss_into_pose( pose );
	} else {
		blueprint_->insert_ss_into_pose( pose );
	}

	// do centroid build
	if ( !centroid_build( pose ) ) { // build failed
		set_last_move_status( FAIL_RETRY );

		// dump pdb when failed
		if ( dump_pdb_when_fail_ != "" ) {
			pose.dump_pdb( dump_pdb_when_fail_ );
		}
		return;
	}

	// if we've gotten to this point, then the structure has been built properly
	set_last_move_status( MS_SUCCESS );
	// output score
	// sfx_->show( TR, pose );

	//fpd reinitialize PDBinfo
	if ( !pose.pdb_info() || pose.pdb_info()->obsolete() ) {
		pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo(pose, true) ) );
	}
}


std::string
BluePrintBDR::get_name() const {
	return BluePrintBDRCreator::mover_name();
}

/// @brief run the centroid level build stage
/// @return true if loop closed, false otherwise
bool BluePrintBDR::centroid_build(
	Pose & pose
)
{
	using core::scoring::ScoreFunctionOP;
	using core::scoring::ScoreFunctionFactory;
	using protocols::moves::MS_SUCCESS;

	using core::util::switch_to_residue_type_set;
	using protocols::forge::methods::restore_residues;
	using protocols::forge::constraints::NtoCConstraintGenerator;
	using protocols::forge::constraints::NtoCConstraintGeneratorOP;
	using protocols::denovo_design::constraints::FileConstraintGenerator;
	using protocols::denovo_design::constraints::FileConstraintGeneratorOP;
	//using protocols::forge::constraints::SheetConstraintsRCG;
	//using protocols::forge::constraints::SheetConstraintsRCGOP;

	using protocols::toolbox::pose_manipulation::construct_poly_XXX_pose;


	// safety, clear the energies object
	pose.energies().clear();

	// make backup Pose for transferring sidechains
	Pose archive_pose = pose;
	Pose modified_archive_pose = archive_pose;
	manager_.modify( modified_archive_pose );

	// ensure modified_archive_pose is completely full-atom, otherwise mismatch
	// will occur when restoring sidechains at the end of the procedure
	bool mod_ap_is_full_atom = true;
	for ( Size i = 1, ie = modified_archive_pose.size(); mod_ap_is_full_atom && i != ie; ++i ) {
		mod_ap_is_full_atom &= ( modified_archive_pose.residue( i ).residue_type_set()->category() == core::chemical::FULL_ATOM_t );
	}

	if ( !mod_ap_is_full_atom ) {
		core::util::switch_to_residue_type_set( modified_archive_pose, core::chemical::FULL_ATOM_t );
	}

	if ( use_poly_val_ ) {
		// flip to poly-ala-gly-pro-disulf pose
		utility::vector1< Size > protein_residues;
		for ( Size i = 1, ie = pose.size(); i <= ie; ++i ) {
			if ( pose.residue( i ).is_protein() ) {
				protein_residues.push_back( i );
			}
		}
		construct_poly_XXX_pose( "VAL", pose, protein_residues, false, true, false );
	}
	/////////////////////////////////////////////////////////////////////////////////////

	// Run VLB to build the new section, if no segments have been added/deleted
	// we use the same VLB so that fragment caching works properly
	if ( !vlb_.get() ) {
		vlb_ = VarLengthBuildOP( new VarLengthBuild( manager_ ) );
	}

	// set weight of constraints
	if ( constraints_NtoC_ > 0.0 || constraints_sheet_ > 0.0 ) {
		Real cst_weight( sfx_->get_weight( core::scoring::atom_pair_constraint ) );
		runtime_assert( cst_weight > 0.0 );
	}

	if ( constraints_NtoC_ > 0.0 ) {
		NtoCConstraintGeneratorOP rcg( new NtoCConstraintGenerator );
		vlb_->add_rcg( rcg );
	}

	if ( constraints_sheet_ > 0.0 ) {
		SheetConstraintGeneratorOP cg( new SheetConstraintGenerator );
		cg->initialize_from_blueprint( blueprint_ );
		cg->set_weight( constraints_sheet_ );
		protocols::forge::remodel::GenericRemodelConstraintGeneratorOP rcg(
			new protocols::forge::remodel::GenericRemodelConstraintGenerator( "bdr", cg ) );
		vlb_->add_rcg( rcg );
	}

	if ( constraint_file_ != "" ) {
		FileConstraintGeneratorOP cst( new FileConstraintGenerator( constraint_file_ ) );
		protocols::forge::remodel::GenericRemodelConstraintGeneratorOP rcg(
			new protocols::forge::remodel::GenericRemodelConstraintGenerator( "bdr_cstfile", cst ) );
		vlb_->add_rcg( rcg );
	}

	// TL: add user-specified RCGS
	for ( core::Size i=1; i<=rcgs_.size(); ++i ) {
		vlb_->add_rcg( rcgs_[i] );
		rcgs_[i]->init( pose );
		TR << "Initialized an RCG called " << rcgs_[i]->get_name() << std::endl;
	}

	vlb_->scorefunction( sfx_ );
	vlb_->vall_memory_usage( protocols::forge::components::VLB_VallMemoryUsage::CLEAR_IF_CACHING_FRAGMENTS );
	vlb_->use_fullmer( use_fullmer_ );
	vlb_->max_linear_chainbreak( max_linear_chainbreak_ );
	vlb_->loop_mover_str( loop_mover_str_ );
	vlb_->num_fragpick( num_fragpick_ );

	if ( use_abego_bias_ ) {
		vlb_->set_abego( blueprint_->abego() );
	}

	if ( tell_vlb_to_not_touch_fold_tree_ ) {
		vlb_->loop_mover_fold_tree_constant( true );
	}

	if ( invrot_tree_ ) {
		this->setup_invrot_tree_in_vlb( *vlb_, pose );
	}


	if ( use_sequence_bias_ ) {
		vlb_->original_sequence( archive_pose.sequence() );
	}

	core::scoring::constraints::ConstraintSet cst ( *pose.constraint_set() );
	if ( TR.visible() ) cst.show(TR);

	vlb_->apply( pose );

	vlb_->clear_rcgs();

	if ( vlb_->get_last_move_status() == MS_SUCCESS ) {

		// record the used manager w/ all mapping info
		manager_ = vlb_->manager();

		// safety, clear all the energies before restoring full-atom residues and
		// scoring
		pose.energies().clear();

		if ( !core::pose::symmetry::is_symmetric( pose ) ) {
			// Swap back original sidechains.  At the moment this is a two step process
			// in case any sidechains from SegmentInsert and the like that aren't in the
			// original archive pose need to be transferred.
			restore_residues( modified_archive_pose, pose );
			restore_residues( manager_.original2modified(), archive_pose, pose );
		} else {
			TR << "Side-chains won't be swapped back to original. " << std::endl;
		}

		return true; // loop closed

	} else {

		pose = archive_pose;

	}

	return false; // false if loop not closed
}

/// @brief parse xml
void
BluePrintBDR::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	Movers_map const & movers,
	Pose const & )
{
	String const blueprint( tag->getOption<std::string>( "blueprint", "" ) );
	if ( blueprint == "" ) {
		TR << "No input of blueprint file ! " << std::endl;
		runtime_assert( false );
	}
	set_blueprint( blueprint );

	// set secondary structure using blueprint
	ss_from_blueprint_ = tag->getOption<bool>( "ss_from_blueprint", 1 );

	// set scorefxn
	String const sfxn ( tag->getOption<String>( "scorefxn", "" ) );
	if ( sfxn != "" ) {
		sfx_ = data.get_ptr<ScoreFunction>( "scorefxns", sfxn );
		TR << "score function, " << sfxn << ", is used. " << std::endl;
	}

	// pick fragment using sequence information
	use_sequence_bias_ = tag->getOption<bool>( "use_sequence_bias", 0 );

	// pick fragment using abego torsion information
	use_abego_bias_ = tag->getOption<bool>( "use_abego_bias", 0 );

	// constraint N- and C- terminal
	constraints_NtoC_ = tag->getOption<Real>( "constraints_NtoC", -1.0 );

	// constraints between Ca atoms in sheet
	constraints_sheet_ = tag->getOption<Real>( "constraints_sheet", -1.0 );

	// set constraint file
	constraint_file_ = tag->getOption<String>( "constraint_file", "" );

	// dump pdb when the protocol fail
	dump_pdb_when_fail_ = tag->getOption<String>( "dump_pdb_when_fail", "" );

	// loop mover for rebuilding loops
	loop_mover_str_ = tag->getOption<String>( "loop_mover", "RemodelLoopMover" );

	// number of allowed_closure_attempts_ of RemodelLoopMover
	rmdl_attempts_ = tag->getOption<Size>( "rmdl_attempts", 1 );

	// entire sequence except for rebuilding parts become poly-Val ( default true )
	use_poly_val_ = tag->getOption<bool>( "use_poly_val", 1 );

	// Use specified constraint generator movers
	// these are called from VLB after the residues are added, but before the actual fragment insertions take place
	utility::vector1< std::string > const mover_names( utility::string_split( tag->getOption< std::string >( "constraint_generators", "" ), ',' ) );
	for ( core::Size i=1; i<=mover_names.size(); ++i ) {
		if ( mover_names[i] == "" ) continue;
		protocols::moves::MoverOP mover = protocols::rosetta_scripts::parse_mover( mover_names[i], movers );
		// check to make sure the mover provided is a RemodelConstraintGenerator and if so, add it to the list
		assert( utility::pointer::dynamic_pointer_cast< protocols::forge::remodel::RemodelConstraintGenerator >( mover ) );
		protocols::forge::remodel::RemodelConstraintGeneratorOP rcg = utility::pointer::static_pointer_cast< protocols::forge::remodel::RemodelConstraintGenerator >( mover );
		rcgs_.push_back( rcg );
		TR << "Added RCG " << mover_names[i] << std::endl;
	}

	if ( tag->hasOption("keep_fold_tree") ) {
		tell_vlb_to_not_touch_fold_tree_ = tag->getOption<bool>( "keep_fold_tree", false );
	}

	//in case we'ref folding up around a ligand
	if ( tag->hasOption("invrot_tree") ) {
		String cstfilename = tag->getOption<String>( "invrot_tree", "");
		toolbox::match_enzdes_util::EnzConstraintIOOP enzcst_io( new toolbox::match_enzdes_util::EnzConstraintIO( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ) );
		enzcst_io->read_enzyme_cstfile( cstfilename );
		invrot_tree_ = protocols::toolbox::match_enzdes_util::InvrotTreeOP( new protocols::toolbox::match_enzdes_util::TheozymeInvrotTree( enzcst_io ) );
		invrot_tree_->generate_targets_and_inverse_rotamers();
		enzcst_io_ = enzcst_io;

		//this also means that we'd like the constraint score terms turned on
		if ( sfx_->has_zero_weight( core::scoring::coordinate_constraint ) ) sfx_->set_weight( core::scoring::coordinate_constraint, 1.0 );
		if ( sfx_->has_zero_weight( core::scoring::atom_pair_constraint ) ) sfx_->set_weight( core::scoring::atom_pair_constraint, 1.0 );
		if ( sfx_->has_zero_weight( core::scoring::angle_constraint ) ) sfx_->set_weight( core::scoring::angle_constraint, 1.0 );
		if ( sfx_->has_zero_weight( core::scoring::dihedral_constraint ) ) sfx_->set_weight( core::scoring::dihedral_constraint, 1.0 );
		if ( sfx_->has_zero_weight( core::scoring::backbone_stub_constraint ) ) sfx_->set_weight( core::scoring::backbone_stub_constraint, 1.0 );

	}
}


} // namespace fldsgn
} // namespace protocols
