// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov

// Unit headers
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/PackRotamersMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSetsFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>



namespace protocols {
namespace minimization_packing {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;

using basic::Warning;
using basic::t_warning;
static basic::Tracer TR( "protocols.minimization_packing.PackRotamersMover" );

// PackRotamersMover

std::string
PackRotamersMoverCreator::keyname() const
{
	return PackRotamersMover::mover_name();
}

protocols::moves::MoverOP
PackRotamersMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PackRotamersMover );
}

void PackRotamersMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PackRotamersMover::provide_xml_schema( xsd );
}

//std::string
//PackRotamersMoverCreator::mover_name()
//{
// return "PackRotamersMover";
//}

PackRotamersMover::PackRotamersMover() :
	protocols::moves::Mover("PackRotamersMover"),
	scorefxn_(/* 0 */),
	task_(/* 0 */),
	nloop_( 1 ), // temporary -- overwritten by the value on the command line
	task_factory_(/* 0 */),
	rotamer_sets_( nullptr ),
	ig_(/* 0 */)
{
	initialize_from_options( basic::options::option );
}

PackRotamersMover::PackRotamersMover(
	utility::options::OptionCollection const & options
) :
	protocols::moves::Mover("PackRotamersMover"),
	scorefxn_(/* 0 */),
	task_(/* 0 */),
	nloop_( 1 ), // temporary -- overwritten by the value on the command line
	task_factory_(/* 0 */),
	rotamer_sets_( nullptr ),
	ig_(/* 0 */)
{
	initialize_from_options( options );
}

PackRotamersMover::PackRotamersMover( std::string const & type_name ) :
	protocols::moves::Mover( type_name ),
	scorefxn_(/* 0 */),
	task_(/* 0 */),
	nloop_( 1 ), // temporary -- overwritten by the value on the command line
	task_factory_(/* 0 */),
	rotamer_sets_( nullptr ),
	ig_(/* 0 */)
{
	initialize_from_options( basic::options::option );
}

// constructors with arguments
PackRotamersMover::PackRotamersMover(
	ScoreFunctionCOP scorefxn,
	PackerTaskCOP task,
	Size nloop
) :
	protocols::moves::Mover("PackRotamersMover"),
	scorefxn_(std::move( scorefxn )),
	task_(std::move( task )),
	nloop_( nloop ),
	task_factory_(/* 0 */),
	rotamer_sets_( nullptr ),
	ig_(/* 0 */)
{}

PackRotamersMover::~PackRotamersMover()= default;

PackRotamersMover::PackRotamersMover( PackRotamersMover const & other ) :
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( other )
{
	scorefxn_ = other.score_function();
	task_ = other.task();
	nloop_ = other.nloop();
	task_factory_ = other.task_factory();
	rotamer_sets_ = RotamerSetsOP( nullptr );
	ig_.reset();
}

void
PackRotamersMover::apply( Pose & pose )
{
	// get rotamers, energies
	// also performs lazy initialization of ScoreFunction, PackerTask
	this->setup( pose );

	core::PackerEnergy best_energy(0.);
	Pose best_pose;
	best_pose = pose;
	for ( Size run(1); run <= nloop_; ++run ) {
		// run SimAnnealer
		core::PackerEnergy packer_energy( this->run( pose ) );
		// Real const score( scorefxn_( pose ) ); another option for deciding which is the 'best' result
		if ( run == 1 || packer_energy < best_energy ) {
			best_pose = pose;
			best_energy = packer_energy;
		}
	}
	if ( nloop_ > 1 ) pose = best_pose;

	//guaruntees proper scoring if this mover is used as a protocol (as in fixbb)
	(*scorefxn_)(pose);
	/// Now handled automatically.  scorefxn_->accumulate_residue_total_energies(pose);
}

std::string
PackRotamersMover::get_name() const { return mover_name(); }

std::string
PackRotamersMoverCreator::mover_name() {
	return "PackRotamersMover";
}

std::string
PackRotamersMover::mover_name() {
	return PackRotamersMoverCreator::mover_name();
}

void
PackRotamersMover::show(std::ostream & output) const
{
	Mover::show(output);
	if ( score_function() != nullptr ) {
		output << "Score function: " << score_function()->get_name() << std::endl;
	} else { output << "Score function: none" << std::endl; }
}

/// @brief when the PackerTask was not generated locally, verify compatibility with pose
/// @details the pose residue types must be equivalent to the ones used to generate the ResidueLevelTasks, because of the way that prevent_repacking and its associated flags work
bool
PackRotamersMover::task_is_valid( Pose const & pose ) const
{
	if ( task_->total_residue() != pose.size() ) return false;
	for ( Size i(1); i <= pose.size(); ++i ) {
		chemical::ResidueTypeCOP r = pose.residue_type_ptr(i);
		if ( ! task_->residue_task(i).is_original_type( r ) ) return false;
	}
	return true;
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
PackRotamersMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose
)
{
	if ( tag->hasOption("nloop") ) {
		nloop_ = tag->getOption<Size>("nloop",1);
		runtime_assert( nloop_ > 0 );
	}
	parse_score_function( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );
}

/// @brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
void
PackRotamersMover::parse_score_function(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	ScoreFunctionOP new_score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	if ( new_score_function == nullptr ) return;
	score_function( new_score_function );
}

/// @brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
void
PackRotamersMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == nullptr ) return;
	task_factory( new_task_factory );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PackRotamersMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new PackRotamersMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PackRotamersMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::minimization_packing::PackRotamersMover( *this ) );
}

/// @brief get rotamers, energies. Also performs lazy initialization of ScoreFunction, PackerTask.
void PackRotamersMover::setup( Pose & pose )
{
	// jec update_residue_neighbors() required to update EnergyGraph (ensures graph_state == GOOD) when calling Interface.cc
	pose.update_residue_neighbors();
	// guarantee of valid ScoreFunction and PackerTask postponed until now
	if ( scorefxn_ == nullptr ) {
		TR.Warning << "undefined ScoreFunction -- creating a default one" << std::endl;
		scorefxn_ = get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	}

	// if present, task_factory_ always overrides/regenerates task_
	if ( task_factory_ != nullptr ) {
		task_ = task_factory_->create_task_and_apply_taskoperations( pose );
	} else if ( task_ == nullptr ) {
		TR.Warning << "undefined PackerTask -- creating a default one" << std::endl;
		task_ = TaskFactory::create_packer_task( pose );
	} else {
		runtime_assert( task_is_valid( pose ) );
	}

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		task_ = core::pack::make_symmetric_task( pose, task_ );
	}
	// in case PackerTask was not generated locally, verify compatibility with pose

	note_packertask_settings( pose );

	rotamer_sets_ = core::pack::rotamer_set::RotamerSetsFactory::create_rotamer_sets( pose );

	pack_rotamers_setup( pose, *scorefxn_, task_, rotamer_sets_, ig_ );
}

core::PackerEnergy PackRotamersMover::run( Pose & pose, utility::vector0< int > rot_to_pack ) const
{
	return pack_rotamers_run( pose, task_, rotamer_sets_, ig_, rot_to_pack );
}

/// @brief note PackerTask's packable and designable residues as string info
void PackRotamersMover::note_packertask_settings( Pose const & pose )
{
	std::ostringstream packable, designable;
	packable << "REMARK PackingRes";
	designable << "REMARK DesignRes";

	for ( Size i(1), end( task_->total_residue() ); i <= end; ++i ) {
		ResidueLevelTask const & rtask( task_->residue_task(i) );
		if ( rtask.being_designed() ) {
			designable << ", ";
			if ( pose.pdb_info() ) {
				designable << pose.pdb_info()->number(i) << " " << pose.pdb_info()->chain(i);
			} else {
				designable << i;
			}
		} else if ( rtask.being_packed() ) {
			packable << ", ";
			if ( pose.pdb_info() ) {
				packable << pose.pdb_info()->number(i) << " " << pose.pdb_info()->chain(i);
			} else {
				packable << i;
			}
		}
	}
	info().clear();
	info().push_back( packable.str() );
	info().push_back( designable.str() );
}

// setters
void PackRotamersMover::score_function( ScoreFunctionCOP sf )
{
	runtime_assert( sf != nullptr );
	scorefxn_ = sf;
}

void PackRotamersMover::task( task::PackerTaskCOP t ) { task_ = t; }

void PackRotamersMover::task_factory( TaskFactoryCOP tf )
{
	runtime_assert( tf != nullptr );
	task_factory_ = tf;
}

void PackRotamersMover::nloop( Size nloop_in ) { nloop_ = nloop_in; }

// accessors
ScoreFunctionCOP PackRotamersMover::score_function() const { return scorefxn_; }
PackerTaskCOP PackRotamersMover::task() const { return task_; }
TaskFactoryCOP PackRotamersMover::task_factory() const { return task_factory_; }
rotamer_set::RotamerSetsCOP PackRotamersMover::rotamer_sets() const { return rotamer_sets_; }


utility::tag::XMLSchemaComplexTypeGeneratorOP
PackRotamersMover::complex_type_generator_for_pack_rotamers_mover( utility::tag::XMLSchemaDefinition & )
{

	using namespace utility::tag;
	AttributeList attributes;

	attributes + XMLSchemaAttribute::attribute_w_default(  "nloop", xsct_non_negative_integer, "Equivalent to \"-ndruns\"."
		"Number of complete packing runs before an output (best score) is produced.",  "1"  );
	rosetta_scripts::attributes_for_parse_score_function( attributes );
	rosetta_scripts::attributes_for_parse_task_operations( attributes );

	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );
	ct_gen->complex_type_naming_func( & moves::complex_type_name_for_mover )
		.add_attributes( attributes )
		.add_optional_name_attribute();
	return ct_gen;
}



void PackRotamersMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = complex_type_generator_for_pack_rotamers_mover( xsd );
	ct_gen->element_name( mover_name() )
		.description( "Repacks sidechains with user-supplied options, including TaskOperations." )
		.write_complex_type_to_schema( xsd );
}

void PackRotamersMover::list_options_read( utility::options::OptionKeyList & opts )
{
	using namespace basic::options::OptionKeys;
	opts + packing::ndruns;
}

void
PackRotamersMover::initialize_from_options( utility::options::OptionCollection const & options )
{
	nloop( options[ basic::options::OptionKeys::packing::ndruns ] );
}

std::ostream &operator<< (std::ostream &os, PackRotamersMover const &mover)
{
	mover.show(os);
	return os;
}

} // moves
} // protocols

