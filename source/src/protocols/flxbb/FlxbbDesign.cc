// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/flxbb/FlxbbDesign.cc
/// @brief perform cycles of design and relax with filter
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit headers
#include <protocols/flxbb/FlxbbDesign.hh>
#include <protocols/flxbb/FlxbbDesignCreator.hh>

// Package Headers
#include <protocols/parser/BluePrint.hh>
#include <protocols/flxbb/DesignTask.hh>
#include <protocols/flxbb/FilterStructs.hh>
#include <protocols/flxbb/utility.hh>
#include <protocols/rosetta_scripts/util.hh>

// Project headers

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <protocols/moves/Mover.hh>

// option key includes
#include <basic/options/keys/flxbb.OptionKeys.gen.hh>

// Utility headers
#include <utility>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// Parser headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/pose_creation/MakePolyXMover.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/task_operations/RestrictToMoveMapChiOperation.hh>
#include <utility/vector0.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.flxbb.FlxbbDesign" );

using namespace core;
using namespace core::pack::task::operation;
using namespace protocols::flxbb;
using namespace basic::options;
using namespace basic::options::OptionKeys;


namespace protocols {
namespace flxbb {

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// XRW TEMP std::string
// XRW TEMP FlxbbDesignCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return FlxbbDesign::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FlxbbDesignCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new FlxbbDesign );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FlxbbDesign::mover_name()
// XRW TEMP {
// XRW TEMP  return "FlxbbDesign";
// XRW TEMP }

/// @brief default constructor
FlxbbDesignPack::FlxbbDesignPack() :
	protocols::minimization_packing::PackRotamersMover( std::string("FlxbbDesignPack") )
{}

FlxbbDesignPack::FlxbbDesignPack(
	ScoreFunctionCOP scorefxn,
	PackerTaskCOP task,
	FilterStructsOP filter ):
	protocols::minimization_packing::PackRotamersMover( scorefxn, task ),
	filter_(std::move( filter ))
{}

/// @brief destructor
FlxbbDesignPack::~FlxbbDesignPack()= default;

/// @brief main operation
void
FlxbbDesignPack::apply( pose::Pose & pose )
{
	pose.update_residue_neighbors();
	this->setup( pose );

	if ( ! filter_ ) {
		this->run( pose );
	} else {
		filter_->reset( pose );
		while ( filter_->filter_on() ) {
			this->run( pose );
			filter_->apply( pose );
		}
		pose = *filter_->get_bestpose();
	}
}

std::string
FlxbbDesignPack::get_name() const {
	return "FlxbbDesignPack";
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
FlxbbDesign::FlxbbDesign() :
	Mover( "FlxbbDesign" ),
	scorefxn_design_( scoring::get_score_function() ),
	scorefxn_relax_ ( scoring::get_score_function() ),
	nflxbb_cycles_( 3 ),
	layer_mode_( "" ),
	use_origseq_for_not_dsgned_layer_( true ),
	no_relax_( false ),
	no_design_( false ),
	filter_during_design_( /* NULL */ ),
	blueprint_( /* NULL */ ),
	resfile_( "" ),
	constraints_sheet_( -1.0 ),
	constraints_NtoC_( -1.0 ),
	movemap_from_blueprint_( false ),
	movemap_( /* NULL */ ),
	use_fast_relax_( true ),
	clear_all_residues_( false ),
	relax_constraint_to_design_( false ),
	limit_aroma_chi2_( true )
{
	design_taskset_.clear();
	task_operations_.clear();
	read_options();
}

/// @brief value constructor
FlxbbDesign::FlxbbDesign(
	ScoreFunctionOP const sfxnd,
	ScoreFunctionOP const sfxnr,
	Size const ncycle,
	String const & layer_mode,
	bool const use_origseq_for_not_dsgned_layer,
	bool const no_relax
) :
	Mover( "FlxbbDesign" ),
	scorefxn_design_ ( sfxnd ),
	scorefxn_relax_ ( sfxnr ),
	nflxbb_cycles_( ncycle ),
	layer_mode_( layer_mode ),
	use_origseq_for_not_dsgned_layer_( use_origseq_for_not_dsgned_layer ),
	no_relax_( no_relax ),
	no_design_( false ),
	filter_during_design_( /* NULL */ ),
	blueprint_( /* NULL */ ),
	resfile_( "" ),
	constraints_sheet_( -1.0 ),
	constraints_NtoC_( -1.0 ),
	movemap_from_blueprint_( false ),
	movemap_( /* NULL */ ),
	use_fast_relax_( true ),
	clear_all_residues_( false ),
	relax_constraint_to_design_( false ),
	limit_aroma_chi2_( true )
{
	design_taskset_.clear();
	task_operations_.clear();
	read_options();
}

/// @brief copy constructor
FlxbbDesign::FlxbbDesign( FlxbbDesign const & ) = default;

/// @brief destructor
FlxbbDesign::~FlxbbDesign()= default;

/// @brief clone this object
FlxbbDesign::MoverOP FlxbbDesign::clone() const
{
	return FlxbbDesign::MoverOP( new FlxbbDesign( *this ) );
}

/// @brief create this type of object
FlxbbDesign::MoverOP FlxbbDesign::fresh_instance() const
{
	return FlxbbDesign::MoverOP( new FlxbbDesign() );
}

/// @brief initialize setups
void FlxbbDesign::read_options()
{
	if ( option[ OptionKeys::flxbb::ncycle ].user() ) nflxbb_cycles_ = option[ OptionKeys::flxbb::ncycle ]();
	if ( option[ OptionKeys::flxbb::blueprint ].user() ) blueprint_ = BluePrintOP( new BluePrint( option[ OptionKeys::flxbb::blueprint ]().name() ) );
	if ( option[ OptionKeys::flxbb::layer::layer ].user() ) layer_mode_ = option[ OptionKeys::flxbb::layer::layer ]();
	if ( option[ OptionKeys::flxbb::constraints_sheet ].user() ) constraints_sheet_ =  option[ OptionKeys::flxbb::constraints_sheet ]();
	if ( option[ OptionKeys::flxbb::constraints_NtoC ].user() ) constraints_NtoC_ = option[ OptionKeys::flxbb::constraints_NtoC ]();
	if ( option[ OptionKeys::flxbb::movemap_from_blueprint ] ) movemap_from_blueprint_ = true;

	Size filter_trial( 0 );
	String filter_type( "" );
	if ( option[ OptionKeys::flxbb::filter_trial ].user() ) filter_trial = option[ OptionKeys::flxbb::filter_trial ]();
	if ( option[ OptionKeys::flxbb::filter_type ].user() ) filter_type = option[ OptionKeys::flxbb::filter_type ]();
	if ( filter_trial > 0 && filter_type != "" ) {
		initialize_filter( filter_trial, filter_type );
	}
}

/// @brief
void FlxbbDesign::initialize_filter( Size const filter_trial, String const & filter_type ){
	if ( filter_type != "" ) {
		if ( filter_type == "packstat" ) {
			filter_during_design_ = FilterStructsOP( new FilterStructs_Packstat( filter_trial ) );
		} else {
			TR.Error << filter_type << " does not exists as filter name " << std::endl;
			utility_exit_with_message("Filter type can't be handled by FlxbbDesign.");
		}
	}
}

/// @details registering of options that are relevant for FlxbbDesign
void FlxbbDesign::register_options()
{
	option.add_relevant( OptionKeys::flxbb::view );
	option.add_relevant( OptionKeys::flxbb::ncycle );
	option.add_relevant( OptionKeys::flxbb::constraints_sheet );
	option.add_relevant( OptionKeys::flxbb::constraints_NtoC );
	option.add_relevant( OptionKeys::flxbb::blueprint );
	option.add_relevant( OptionKeys::flxbb::filter_trial );
	option.add_relevant( OptionKeys::flxbb::filter_type );
	option.add_relevant( OptionKeys::flxbb::layer::layer );
}

/// @brief set scorefxn for fixbb design
void FlxbbDesign::set_scorefxn_design( ScoreFunctionOP const scorefxn )
{
	scorefxn_design_ = scorefxn;
}

/// @brief set scorefxn for relax
void FlxbbDesign::set_scorefxn_relax( ScoreFunctionOP const scorefxn )
{
	scorefxn_relax_ = scorefxn;
}

/// @brief the number of cycles of fixbb and relax
void FlxbbDesign::set_ncycles( Size const ncycles )
{
	nflxbb_cycles_ = ncycles;
}

/// @brief set layer mode
void FlxbbDesign::set_layer_mode( String const & layer_mode )
{
	layer_mode_ = layer_mode;
}

/// @brief use original sequence for not designed region in layer_mode,
/// otherwise the residues of the region turned into Ala ( default true )
void FlxbbDesign::use_origseq_for_not_dsgned_layer( bool const use )
{
	use_origseq_for_not_dsgned_layer_ = use;
}

/// @brief relax is not performed after design (default false)
void FlxbbDesign::no_relax( bool const no_relax )
{
	no_relax_ = no_relax;
}

/// @brief design protocol will not be performed (default false)
void FlxbbDesign::no_design( bool const no_design )
{
	no_design_ = no_design;
}

/// @brief set BlueprintOP
void FlxbbDesign::set_blueprint( BluePrintOP const blueprint )
{
	blueprint_ = blueprint;
}

/// @brief set weight of constraints_sheet which constrains between Ca atoms in beta-sheets
/// if this weight > 0.0, the constraints is applied ( default = -1.0 )
void FlxbbDesign::set_weight_constraints_sheet( Real const value )
{
	constraints_sheet_ = value;
}

/// @brief set weight of constraints_NtoC which constrain between Ca atoms of C- and N-terminal
/// if this weight > 0.0, the constraints is applied ( default = -1.0 )
void FlxbbDesign::set_weight_constraints_NtoC( Real const value )
{
	constraints_NtoC_ = value;
}

/// @brief
void FlxbbDesign::movemap_from_blueprint( bool const value )
{
	movemap_from_blueprint_ = value;
}

/// @brief set movemap factory for relax
void FlxbbDesign::set_movemap_factory( core::select::movemap::MoveMapFactoryCOP movemap_factory )
{
	movemap_factory_ = movemap_factory;
}
/// @brief set movemap for relax
void FlxbbDesign::set_movemap( MoveMapOP const movemap )
{
	movemap_ = movemap;
}

/// @brief clear DesignTaskSet
void FlxbbDesign::set_design_taskset( DesignTaskSet const & design_taskset )
{
	design_taskset_ = design_taskset;
}

/// @brief set DesignTask
void FlxbbDesign::add_design_task( DesignTaskOP const design_task )
{
	design_taskset_.push_back( design_task );
}

/// @brief filtering during design
void FlxbbDesign::set_filter_during_design( FilterStructsOP const filter_during_design )
{
	filter_during_design_ = filter_during_design;
}

/// @brief clear DesignTaskSet
void FlxbbDesign::clear_design_taskset()
{
	design_taskset_.clear();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief build design_task_set
DesignTaskSet
FlxbbDesign::build_design_taskset( Pose const & pose )
{
	using core::pose::PoseOP;
	using core::pack::task::PackerTaskOP;
	using core::pack::task::TaskFactory;
	using core::pack::task::TaskFactoryOP;
	using protocols::moves::MoverOP;
	using protocols::relax::ClassicRelax;
	using protocols::relax::FastRelax;
	using protocols::relax::RelaxProtocolBaseOP;
	using protocols::flxbb::FilterStructs_Packstat;
	using protocols::flxbb::FilterStructs_PackstatOP;
	using protocols::flxbb::FilterStructs_TotalCharge;
	using protocols::flxbb::FilterStructs_TotalChargeOP;
	using protocols::task_operations::LimitAromaChi2Operation;
	using protocols::task_operations::RestrictToMoveMapChiOperation;
	DesignTaskSet dts;

	// create relax mover
	RelaxProtocolBaseOP rlx_mover;
	if ( use_fast_relax_ ) {
		rlx_mover = RelaxProtocolBaseOP( new FastRelax( scorefxn_relax_ ) );
	} else {
		rlx_mover = RelaxProtocolBaseOP( new ClassicRelax( scorefxn_relax_ ) );
	}
	if ( movemap_ ) {
		rlx_mover->set_movemap( movemap_ ); // movemap can be controled by blueprint
	} else if ( movemap_factory_ ) {
		rlx_mover->set_movemap( movemap_factory_->create_movemap_from_pose( pose ) );
	}
	if ( relax_constraint_to_design_ ) rlx_mover->constrain_relax_to_start_coords( true );
	if ( limit_aroma_chi2_ ) { // default true
		TaskFactoryOP tf( new TaskFactory );
		//jadolfbr - FastRelax and ClassicRelax now completely respects the tf you give it.
		//Nobu - Note that ClassicRelax ignored the set task factory before.
		tf->push_back(TaskOperationCOP( new InitializeFromCommandline ));
		tf->push_back(TaskOperationCOP( new RestrictToRepacking ) );
		tf->push_back( TaskOperationCOP( new LimitAromaChi2Operation ) );
		if ( movemap_ ) {
			tf->push_back( TaskOperationCOP( new RestrictToMoveMapChiOperation(movemap_) ) );
		} else if ( movemap_factory_ ) {
			tf->push_back( TaskOperationCOP( new RestrictToMoveMapChiOperation( movemap_factory_->create_movemap_from_pose( pose ) ) ) );
		}
		rlx_mover->set_task_factory( tf );
	}

	if ( layer_mode_ == "" ) {
		dts.push_back( DesignTaskOP( new DesignTask_Normal( nflxbb_cycles_, scorefxn_design_, rlx_mover, filter_during_design_ ) ) );
	} else {

		utility::vector1< String > layers( utility::string_split( layer_mode_, '_' ) );
		if ( layers.size() == 1 && layer_mode_ == "normal" ) {

			rlx_mover = RelaxProtocolBaseOP( new FastRelax( scorefxn_relax_ ) );
			FilterStructs_PackstatOP filter1( new FilterStructs_Packstat( pose, 10 ) );
			FilterStructs_TotalChargeOP filter2( new FilterStructs_TotalCharge( pose ) );

			// 1st stage: 2 cycles of the fixbb designinig only core and boundary and relax
			dts.push_back( DesignTaskOP( new DesignTask_Layer( true, true, false, false,
				2, scorefxn_design_, rlx_mover, filter1 ) ));
			// 2nd stage: 2 cycles of design all and relax
			dts.push_back( DesignTaskOP( new DesignTask_Layer( true, true, true, false,
				2, scorefxn_design_, rlx_mover, filter1 ) ));
			// 3rd stage: design only surface without relax and filter sequence of which totalcharge is non-zero
			dts.push_back( DesignTaskOP( new DesignTask_Layer( false, false, true, true,
				1, scorefxn_design_, nullptr, filter2 ) ));

		} else if ( layers.size() == 1 && layer_mode_ == "all" ) {
			dts.push_back( DesignTaskOP( new DesignTask_Layer( true, true, true, false,
				nflxbb_cycles_, scorefxn_design_, rlx_mover, filter_during_design_ ) ));
		} else if ( layers.size() >= 1 ) {
			bool dcore( false );
			bool boundary( false );
			bool surface( false );
			for ( utility::vector1< String >::const_iterator iter = layers.begin(); iter != layers.end() ; ++iter ) {
				String layer(*iter);
				if ( layer == "core" ) {
					dcore = true;
				} else if ( layer == "surface" ) {
					surface = true;
				} else if ( layer == "boundary" ) {
					boundary = true;
				} else {
					TR << "Error!, wrong specification of layer_mode " << layer << std::endl;
					TR << "Every layers are designed. " << std::endl;
					dcore = true;
					surface = true;
					boundary = true;
				}
			} // utility::vector1
			dts.push_back( DesignTaskOP( new DesignTask_Layer( dcore, boundary, surface,
				use_origseq_for_not_dsgned_layer_,
				nflxbb_cycles_, scorefxn_design_, rlx_mover, filter_during_design_ ) ));
		}
	} // if ( layer_mode_ == "" )

	// set the filename of resfile, which is given in parser, to all the DesignTask
	if ( resfile_ != "" ) {
		for ( DesignTaskSet::const_iterator it=dts.begin(); it!=dts.end(); ++it ) {
			(*it)->set_resfile( resfile_ );
		}
	}

	// set additional task operations
	if ( ! task_operations_.empty() ) {
		for ( DesignTaskSet::const_iterator it=dts.begin(); it!=dts.end(); ++it ) {
			(*it)->add_task_operations( task_operations_ );
		}
	}

	// Exclude aromatic chi2 rotamers, of which angles are around 0
	if ( limit_aroma_chi2_ ) { // default true
		for ( DesignTaskSet::const_iterator it=dts.begin(); it!=dts.end(); ++it ) {
			(*it)->add_task_operation( TaskOperationOP( new LimitAromaChi2Operation ) );
		}
	}

	return dts;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief mover apply
void FlxbbDesign::apply( pose::Pose & pose )
{
	using core::util::switch_to_residue_type_set;
	using core::pack::task::PackerTaskOP;
	using core::pack::task::TaskFactory;
	using core::scoring::constraints::ConstraintSet;
	using core::scoring::constraints::ConstraintSetOP;
	using protocols::flxbb::DesignTaskSet;
	using protocols::flxbb::DesignTask_Normal;
	using protocols::moves::MoverOP;
	using protocols::relax::FastRelax;
	using protocols::relax::ClassicRelax;
	using protocols::pose_creation::MakePolyXMover;

	// set pose to fullatom
	if ( ! pose.is_fullatom() ) {
		core::util::switch_to_residue_type_set( pose, core::chemical::FULL_ATOM_t );
	}

	// set movemap from blueprint for relax
	if ( movemap_from_blueprint_ ) {
		if ( movemap_ || movemap_factory_ ) {
			TR << "Movemap will be overrided by the definition of movemap in blueprint " << std::endl;
		}
		runtime_assert( blueprint_ != nullptr );
		runtime_assert( pose.size() == blueprint_->total_residue() );
		movemap_ = MoveMapOP( new core::kinematics::MoveMap );
		blueprint_->set_movemap( movemap_ );
	}

	// set constraints
	ConstraintSetOP cstset( new ConstraintSet );

	// set weight of constraints
	if ( constraints_NtoC_ > 0.0 || constraints_sheet_ > 0.0 ) {
		Real cst_weight( scorefxn_relax_->get_weight( core::scoring::atom_pair_constraint ) );
		runtime_assert( cst_weight > 0.0 );
	}

	// constraints in beta-sheet
	if ( constraints_sheet_ > 0.0 ) {
		if ( blueprint_ ) {
			cstset->add_constraints( constraints_sheet( pose, blueprint_, constraints_sheet_ ) );
		} else {
			cstset->add_constraints( constraints_sheet( pose, constraints_sheet_ ) );
		}
	}
	// constraints between N and C
	if ( constraints_NtoC_ > 0.0 ) {
		cstset->add_constraints( constraints_NtoC( pose, constraints_NtoC_ ) );
	}
	// attach constraints to pose
	pose.add_constraints( cstset->get_all_constraints() );

	// setup design_taskset
	if ( design_taskset_.empty() ) {
		design_taskset_ = build_design_taskset( pose );
	}

	// make pose to all ala
	// alert !! we might want to have the functionality to keep_disulfide or not
	if ( clear_all_residues_ ) {
		MakePolyXMover bap( "ALA", false/*kee_pro*/, true /*keep_gly*/, false /*keep_disulfide_cys*/ );
		bap.apply( pose );
	}

	// run
	Size num_task( 0 );
	for ( auto design_task : design_taskset_ ) {

		num_task ++;
		for ( Size i=1 ; i<=design_task->ncycle() ; i++ ) {

			TR << "current_cycle/total_cycle: " << i << "/" << design_task->ncycle() << " in DesignTask: " << num_task << std::endl;

			if ( ! no_design_ ) {
				// set packertask
				PackerTaskOP task( TaskFactory::create_packer_task( pose ));
				design_task->setup( pose, task );
				design_task->dump_packertask( TR );

				// run design
				runtime_assert( design_task->scorefxn() != nullptr );
				FlxbbDesignPack pack( design_task->scorefxn(), design_task->packertask(), design_task->filter_structs() );
				pack.apply( pose );

				TR << "Designed sequence: " << pose.sequence() << std::endl;
				TR << "Score after design by fixbb " << std::endl;
				scorefxn_design_->show( TR, pose );
				TR.flush();
			} else {
				TR << "Design step was skipped. ";
			}

			// run mover
			if ( design_task->mover() && !no_relax_ ) {
				design_task->mover()->apply( pose );
				TR << "Score after mover, " << design_task->mover()->type() << std::endl;
				scorefxn_relax_->show( TR, pose );
				TR.flush();
			}

		} // iteration of i: design_task->ncycle()
	} // iteration of it: DesignTaskSet


	// detach constraints from pose
	pose.remove_constraints( cstset->get_all_constraints() );

	// calc packstat
	core::Real packscore;
	packscore = core::scoring::packstat::compute_packing_score( pose );
	setPoseExtraScore( pose, "pack_stat", packscore );

} // FlxbbDesign::apply

// XRW TEMP std::string
// XRW TEMP FlxbbDesign::get_name() const {
// XRW TEMP  return FlxbbDesign::mover_name();
// XRW TEMP }

/// @brief parse xml
void
FlxbbDesign::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	scorefxn_design_ = protocols::rosetta_scripts::parse_score_function( tag, "sfxn_design", data );
	scorefxn_relax_ = protocols::rosetta_scripts::parse_score_function( tag, "sfxn_relax", data );

	std::string const blueprint( tag->getOption<std::string>( "blueprint", "" ) );
	if ( blueprint != "" ) {
		blueprint_ = BluePrintOP( new BluePrint( blueprint ) );
	}

	// the number of cycles of fixbb and relax
	nflxbb_cycles_ = tag->getOption<Size>( "ncycles", 3 );

	// perform fixbb in layer mode: core, boundary, and surfacee
	layer_mode_ = tag->getOption<String>( "layer_mode", "" );

	// use original sequence for not designed region in layer_mode
	use_origseq_for_not_dsgned_layer_ = tag->getOption<bool>( "use_original_seq", true );

	// set filter
	auto filter_trial = tag->getOption<Size>( "filter_trial", 10 );
	String filter_type = tag->getOption<String>( "filter_type", "packstat" );
	if ( filter_trial > 0 && filter_type != "" ) {
		initialize_filter( filter_trial, filter_type );
	}

	// do not perform relax after design
	no_relax_ = tag->getOption<bool>( "no_relax", false );
	movemap_factory_ = protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data );

	// do not perform design
	no_design_ = tag->getOption<bool>( "no_design", false );

	/// clear amino acid sequence at the very beginning of FlxbbDesign
	clear_all_residues_ = tag->getOption<bool>( "clear_all_residues", false );

	// resfile for fixbb
	resfile_ = tag->getOption<String>( "resfile", "" );

	// constraint N- and C- terminal
	constraints_NtoC_ = tag->getOption<Real>( "constraints_NtoC", -1.0 );

	// constraint N- and C- terminal
	constraints_sheet_ = tag->getOption<Real>( "constraints_sheet", -1.0 );

	// constraint of backbone to fixbb-designed structure during relax
	relax_constraint_to_design_ = tag->getOption<bool>( "constraints_to_backbone", false );

	// movemap for relax
	movemap_from_blueprint_ = tag->getOption<bool>( "movemap_from_blueprint", false );
	if ( movemap_from_blueprint_ ) {
		TR << "Movemap is defined based on blueprint " << std::endl;
		runtime_assert( blueprint != "" );
	}

	// read task operations
	task_operations_ = protocols::rosetta_scripts::get_task_operations( tag, data );

	/// use fast relax
	use_fast_relax_ = tag->getOption<bool>( "fast_relax", true );

	// Exclude aromatic chi2 rotamers, of which angles are around 0
	limit_aroma_chi2_ = tag->getOption<bool>( "limit_aroma_chi2", true );

	/// read options  Alert!! if there are options, settings of parser will be overrided.
	read_options();

}

std::string FlxbbDesign::get_name() const {
	return mover_name();
}

std::string FlxbbDesign::mover_name() {
	return "FlxbbDesign";
}

void FlxbbDesign::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// AMW TODO
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function( attlist, "sfxn_design" );
	rosetta_scripts::attributes_for_parse_score_function( attlist, "sfxn_relax" );
	attlist + XMLSchemaAttribute( "blueprint", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "ncycles", xsct_non_negative_integer, "XRW TO DO", "3" )
		+ XMLSchemaAttribute( "layer_mode", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "use_original_seq", xsct_rosetta_bool, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_trial", xsct_non_negative_integer, "XRW TO DO", "10" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_type", xs_string, "XRW TO DO", "packstat" )
		+ XMLSchemaAttribute::attribute_w_default( "no_relax", xsct_rosetta_bool, "XRW TO DO", "0" );

	XMLSchemaSimpleSubelementList ssl;
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy( xsd, ssl );
	attlist + XMLSchemaAttribute::attribute_w_default( "no_design", xsct_rosetta_bool, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "clear_all_residues", xsct_rosetta_bool, "XRW TO DO", "0" )
		+ XMLSchemaAttribute( "resfile", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "constraints_NtoC", xsct_real, "XRW TO DO", "-1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "constraints_sheet", xsct_real, "XRW TO DO", "-1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "constraints_to_backbone", xsct_rosetta_bool, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "movemap_from_blueprint", xsct_rosetta_bool, "XRW TO DO", "0" );

	rosetta_scripts::attributes_for_parse_task_operations( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default( "fast_relax", xsct_rosetta_bool, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "limit_aroma_chi2", xsct_rosetta_bool, "XRW TO DO", "1" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, ssl );
}

std::string FlxbbDesignCreator::keyname() const {
	return FlxbbDesign::mover_name();
}

protocols::moves::MoverOP
FlxbbDesignCreator::create_mover() const {
	return protocols::moves::MoverOP( new FlxbbDesign );
}

void FlxbbDesignCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FlxbbDesign::provide_xml_schema( xsd );
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


} // namespace of flxbb
} // namespace of protocols
