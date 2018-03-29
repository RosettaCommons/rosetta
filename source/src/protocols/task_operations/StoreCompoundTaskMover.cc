// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/StoreCompoundTaskMover.cc
/// @brief  Combine tasks using boolean logic for residues that are packable or designable,
/// assign new packing behavior to residues that match or do not match the specified criteria,
/// and store the resulting task in the pose's cacheable data.
/// TODO (Jacob): Add option to combine allowed amino acid sets.
/// @author Jacob Bale (balej@uw.edu) (Much of this code was adapted from the CompoundStatementFilter, StoreTaskMover, and rosetta_scripts/util.cc)

// Unit Headers
#include <protocols/task_operations/StoreCompoundTaskMover.hh>
#include <protocols/task_operations/StoreCompoundTaskMoverCreator.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pack/make_symmetric_task.hh>
#include <protocols/task_operations/STMStoredTask.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <ObjexxFCL/format.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.task_operations.StoreCompoundTaskMover" );

namespace protocols {
namespace task_operations {

// @brief default constructor
StoreCompoundTaskMover::StoreCompoundTaskMover() = default;

// @brief destructor
StoreCompoundTaskMover::~StoreCompoundTaskMover() = default;

void StoreCompoundTaskMover::task_clear() { compound_task_.clear(); }
StoreCompoundTaskMover::task_iterator StoreCompoundTaskMover::task_begin() { return( compound_task_.begin() ); }
StoreCompoundTaskMover::const_task_iterator StoreCompoundTaskMover::task_begin() const { return( compound_task_.begin() ); }
StoreCompoundTaskMover::task_iterator StoreCompoundTaskMover::task_end() { return( compound_task_.end() ); }
StoreCompoundTaskMover::const_task_iterator StoreCompoundTaskMover::task_end() const { return( compound_task_.end() ); }

void StoreCompoundTaskMover::factory_clear() { compound_factory_.clear(); }
StoreCompoundTaskMover::factory_iterator StoreCompoundTaskMover::factory_begin() { return( compound_factory_.begin() ); }
StoreCompoundTaskMover::const_factory_iterator StoreCompoundTaskMover::factory_begin() const { return( compound_factory_.begin() ); }
StoreCompoundTaskMover::factory_iterator StoreCompoundTaskMover::factory_end() { return( compound_factory_.end() ); }
StoreCompoundTaskMover::const_factory_iterator StoreCompoundTaskMover::factory_end() const { return( compound_factory_.end() ); }

void StoreCompoundTaskMover::invert( bool const inv ) { invert_ = inv; }
void StoreCompoundTaskMover::verbose( bool const verb ) { verbose_ = verb; }
void StoreCompoundTaskMover::overwrite( bool const ow ) { overwrite_ = ow; }
void StoreCompoundTaskMover::task_name( std::string const & tn ) { task_name_ = tn; }
void StoreCompoundTaskMover::mode( std::string const & md ) { mode_ = md; }
void StoreCompoundTaskMover::true_behavior( std::string const & tb ) { true_behavior_ = tb; }
void StoreCompoundTaskMover::false_behavior( std::string const & fb ) { false_behavior_ = fb; }

void
StoreCompoundTaskMover::CompoundPackableTask( core::pose::Pose const & pose, core::Size & total_residue, core::pack::task::PackerTaskOP & task)
{
	std::string select_true_pos(task_name_+": select true_positions, resi ");
	for ( core::Size resi=1; resi<=total_residue; ++resi ) {

		bool value( true );

		for ( StoreCompoundTaskMover::const_task_iterator it=compound_task_.begin(); it!=compound_task_.end(); ++it ) {
			if ( it - compound_task_.begin() == 0 ) {
				// first logical op may only be NOT
				// ANDNOT and ORNOT are also treated as NOT (with a warning)
				value = it->first->being_packed( resi );
				if ( it->second == NOT ) value = !value;
				if ( it->second == ORNOT ) {
					TR.Warning << "StoreCompoundTaskMover treating operator ORNOT as NOT" << std::endl;
					value = !value;
				}
				if ( it->second == ANDNOT ) {
					TR.Warning << "StoreCompoundTaskMover treating operator ANDNOT as NOT" << std::endl;
					value = !value;
				}
			} else {
				switch( it->second  ) {
				case ( AND ) : value = value && it->first->being_packed( resi ); break;
				case ( OR  ) : value = value || it->first->being_packed( resi ); break;
				case ( XOR ) : value = value ^ it->first->being_packed( resi ); break;
				case ( ORNOT ) : value = value || !it->first->being_packed( resi ); break;
				case ( ANDNOT ) : value = value && !it->first->being_packed( resi ); break;
				case ( NOR ) : value = !( value || it->first->being_packed( resi ) ); break;
				case (NAND ) : value = !( value && it->first->being_packed( resi ) ); break;
				case (NOT ) :
					TR.Warning << "StoreCompoundTaskMover treating operator NOT as ANDNOT" << std::endl;
					value = value && !it->first->being_packed( resi );
					break;
				}
			}
		}
		if ( invert_ ) value = !value;
		if ( !value ) {
			if ( false_behavior_ == "prevent_repacking" ) {
				task->nonconst_residue_task(resi).prevent_repacking();
			} else if ( false_behavior_ == "restrict_to_repacking" ) {
				task->nonconst_residue_task(resi).restrict_to_repacking();
			}
		} else {
			if ( true_behavior_ == "prevent_repacking" ) {
				task->nonconst_residue_task(resi).prevent_repacking();
			} else if ( true_behavior_ == "restrict_to_repacking" ) {
				task->nonconst_residue_task(resi).restrict_to_repacking();
			}
			core::Size output_resi = resi;
			if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
				output_resi = pose.pdb_info()->number( resi );
			}
			select_true_pos.append(ObjexxFCL::string_of(output_resi) + "+");
		}
	}
	if ( verbose_ ) {
		TR << select_true_pos << std::endl;
	}
}

void
StoreCompoundTaskMover::CompoundDesignableTask( core::pose::Pose const & pose, core::Size & total_residue, core::pack::task::PackerTaskOP & task)
{
	std::string select_true_pos(task_name_+": select true_positions, resi ");
	for ( core::Size resi=1; resi<=total_residue; ++resi ) {

		bool value( true );

		for ( StoreCompoundTaskMover::const_task_iterator it=compound_task_.begin(); it!=compound_task_.end(); ++it ) {
			if ( it - compound_task_.begin() == 0 ) {
				// first logical op may only be NOT
				// ANDNOT and ORNOT are also treated as NOT (with a warning)
				value = it->first->being_designed( resi );
				if ( it->second == NOT ) value = !value;
				if ( it->second == ORNOT ) {
					TR.Warning << "StoreCompoundTaskMover treating operator ORNOT as NOT" << std::endl;
					value = !value;
				}
				if ( it->second == ANDNOT ) {
					TR.Warning << "StoreCompoundTaskMover treating operator ANDNOT as NOT" << std::endl;
					value = !value;
				}
			} else {
				switch( it->second  ) {
				case ( AND ) : value = value && it->first->being_designed( resi ); break;
				case ( OR  ) : value = value || it->first->being_designed( resi ); break;
				case ( XOR ) : value = value ^ it->first->being_designed( resi ); break;
				case ( ORNOT ) : value = value || !it->first->being_designed( resi ); break;
				case ( ANDNOT ) : value = value && !it->first->being_designed( resi ); break;
				case ( NOR ) : value = !( value || it->first->being_designed( resi ) ); break;
				case (NAND ) : value = !( value && it->first->being_designed( resi ) ); break;
				case (NOT ) :
					TR.Warning << "StoreCompoundTaskMover treating operator NOT as ANDNOT" << std::endl;
					value = value && !it->first->being_designed( resi );
					break;
				}
			}
		}
		if ( invert_ ) value = !value;
		if ( !value ) {
			if ( false_behavior_ == "prevent_repacking" ) {
				task->nonconst_residue_task(resi).prevent_repacking();
			} else if ( false_behavior_ == "restrict_to_repacking" ) {
				task->nonconst_residue_task(resi).restrict_to_repacking();
			}
		} else {
			if ( true_behavior_ == "prevent_repacking" ) {
				task->nonconst_residue_task(resi).prevent_repacking();
			} else if ( true_behavior_ == "restrict_to_repacking" ) {
				task->nonconst_residue_task(resi).restrict_to_repacking();
			}
			core::Size output_resi = resi;
			if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
				output_resi = pose.pdb_info()->number( resi );
			}
			select_true_pos.append(ObjexxFCL::string_of(output_resi) + "+");
		}
	}
	if ( verbose_ ) {
		TR << select_true_pos << std::endl;
	}
}

void
StoreCompoundTaskMover::apply( core::pose::Pose & pose )
{

	// Only consider the residues in the asymmetric unit if the pose is symmetric
	core::Size total_residue;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
		total_residue = symm_info->num_independent_residues();
	} else {
		total_residue = pose.size();
	}

	// Loop over the task factory, boolean operation pairs, apply the task operations to the pose, and store these new task, operator pairs
	// Note: This is performed here rather than in the parse_my_tag function so that the PackerTask is created using the pose present at
	// runtime rather than during parsing.
	for ( StoreCompoundTaskMover::const_factory_iterator it=compound_factory_.begin(); it!=compound_factory_.end(); ++it ) {
		std::pair< core::pack::task::PackerTaskOP, boolean_operations > task_pair;
		task_pair.second = it->second;
		core::pack::task::PackerTaskOP new_packer_task = it->first->create_task_and_apply_taskoperations( pose );
		task_pair.first = new_packer_task->clone(); //clone?
		runtime_assert( new_packer_task != nullptr );
		compound_task_.push_back( task_pair );
	}

	// Create blank task to be modified to yield the new compound task.
	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );

	////////////////////////////////////////////////////////////////////////////////////////////////
	/// For each of the possible modes:
	/// 1) loop over all of the independent residues
	/// 2) for each residue
	///  1) loop through the task/boolean_operation pairs
	///  2) determine the final state of each residue (ie, packable, designable, designable with amino acids TWER... etc)
	///  3) apply this state to the new compound task
	////////////////////////////////////////////////////////////////////////////////////////////////

	if ( mode_ == "packable" ) {
		CompoundPackableTask( pose, total_residue, task );
	} else if ( mode_ == "designable" ) {
		CompoundDesignableTask( pose, total_residue, task );
	}

	/* TODO (Jacob): Add option to combine aa_sets.
	Does a getter function exist to get a length-20 vector of bools for allowed_aas for a given residue?
	Does a function exist to check if an amino acid (given by either one or three letter aa code) is allowed at a given position?
	Add addition loop over allowed_aa vector for each residue, storing value in a vector of bools, which is then used to set
	the new set of allowed amino acids via restrict_absent_canonical_aas( final vector of bools )
	} else if( mode_ == "allowed_aas" ) {
	// utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, false )
	// for(ResidueLevelTask::ResidueTypeCOPListConstIter aa_iter(storage_task->residue_task(i).allowed_residue_types_begin()),
	// aa_end(storage_task->residue_task(i).allowed_residue_types_end());
	// aa_iter != aa_end; ++aa_iter){
	// residues_to_mutate[i]=((*aa_iter)->aa());

	}
	*/
	//////////////////////////////////////////////////////////////////////////////////////////////////
	/// Store the new compound task in the pose's cacheable data.
	/////////////////////////////////////////////////////////////////////////////////////////////////
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		core::pack::make_symmetric_PackerTask_by_truncation(pose, task); // Does this need to be fixed or omitted?
	}
	// If the pose doesn't have STM_STORED_TASK data, put a blank STMStoredTask in there.
	if ( !pose.data().has( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) {
		protocols::task_operations::STMStoredTaskOP blank_tasks( new protocols::task_operations::STMStoredTask() );
		pose.data().set( core::pose::datacache::CacheableDataType::STM_STORED_TASKS, blank_tasks );
	}
	// Grab a reference to the data
	protocols::task_operations::STMStoredTask & stored_tasks = *( utility::pointer::static_pointer_cast< protocols::task_operations::STMStoredTask > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::STM_STORED_TASKS ) ) );
	// If you haven't set overwrite to true and your task name already exists, fail. Otherwise, put the task you've made into the data cache.
	if ( overwrite_ || !stored_tasks.has_task(task_name_) ) {
		stored_tasks.set_task( task, task_name_ );
	} else {
		utility_exit_with_message("A stored task with the name " + task_name_ + " already exists; you must set overwrite flag to true to overwrite." );
	}
	compound_task_.clear();
}

void
StoreCompoundTaskMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/ )
{

	using StringVec = utility::vector1<std::string>;

	TR<<"StoreCompoundTask: "<<tag->getName()<<std::endl;
	task_name_ = tag->getOption< std::string >( "task_name", "" );
	true_behavior_ = tag->getOption<std::string>( "true_behavior", "" );
	if ( !( (true_behavior_ == "prevent_repacking") || (true_behavior_ == "restrict_to_repacking") || (true_behavior_ == "") ) ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Error: true_behavior in tag is undefined." );
	}
	false_behavior_ = tag->getOption<std::string>( "false_behavior", "prevent_repacking" );
	if ( !( (false_behavior_ == "prevent_repacking") || (false_behavior_ == "restrict_to_repacking") || (false_behavior_ == "") ) ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Error: false_behavior in tag is undefined." );
	}
	mode_ = tag->getOption<std::string>( "mode", "packable" );
	if ( !( (mode_ == "packable") || (mode_ == "designable") ) ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Error: mode in tag is undefined." );
	}
	invert_ = tag->getOption<bool>( "invert", false );
	verbose_ = tag->getOption<bool>( "verbose", false );
	overwrite_ = tag->getOption< bool >( "overwrite", false );

	/////////////////////////////////////////////////////////////////////////////////////////////////
	/// Loop through all user-provided subtags (ex. < AND task_name="bbi" />) and put these into a
	/// vector of (TaskFactoryOP, boolean_operation) pairs.
	/// Note: Do not apply tasks to pose until runtime
	/////////////////////////////////////////////////////////////////////////////////////////////////
	for ( TagCOP cmp_tag_ptr : tag->getTags() ) {
		std::string const operation( cmp_tag_ptr->getName() );
		std::pair< core::pack::task::TaskFactoryOP, boolean_operations > factory_pair;
		if ( operation == "AND" ) factory_pair.second = AND;
		else if ( operation == "OR" ) factory_pair.second = OR;
		else if ( operation == "XOR" ) factory_pair.second = XOR;
		else if ( operation == "NOR" ) factory_pair.second = NOR;
		else if ( operation == "NAND" ) factory_pair.second = NAND;
		else if ( operation == "ORNOT" ) factory_pair.second = ORNOT;
		else if ( operation == "ANDNOT" ) factory_pair.second = ANDNOT;
		else if ( operation == "NOT" ) factory_pair.second = NOT;
		else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Error: Boolean operation in tag is undefined." );
		}

		core::pack::task::TaskFactoryOP new_task_factory( new core::pack::task::TaskFactory );
		std::string const t_o_val( cmp_tag_ptr->getOption<std::string>("task_operations") );

		TR<<"Defined with operator: "<<operation<<" and tasks: "<<t_o_val<<std::endl;

		StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
		for ( auto const & t_o_key : t_o_keys ) {
			if ( data_map.has( "task_operations", t_o_key ) ) {
				new_task_factory->push_back( data_map.get_ptr< core::pack::task::operation::TaskOperation >( "task_operations", t_o_key ) );
			} else {
				utility_exit_with_message("TaskOperation " + t_o_key + " not found in basic::datacache::DataMap.");
			}
		}
		factory_pair.first = new_task_factory->clone(); //clone?
		runtime_assert( new_task_factory != nullptr );
		compound_factory_.push_back( factory_pair );
	}
}

// @brief Identification
// XRW TEMP std::string StoreCompoundTaskMoverCreator::keyname() const { return StoreCompoundTaskMover::mover_name(); }
// XRW TEMP std::string StoreCompoundTaskMover::mover_name() { return "StoreCompoundTaskMover"; }
// XRW TEMP std::string StoreCompoundTaskMover::get_name() const { return "StoreCompoundTaskMover"; }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP StoreCompoundTaskMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new StoreCompoundTaskMover );
// XRW TEMP }

protocols::moves::MoverOP
StoreCompoundTaskMover::clone() const {
	return protocols::moves::MoverOP( new StoreCompoundTaskMover( *this ) );
}

protocols::moves::MoverOP
StoreCompoundTaskMover::fresh_instance() const {
	return protocols::moves::MoverOP( new StoreCompoundTaskMover );
}

std::string StoreCompoundTaskMover::get_name() const {
	return mover_name();
}

std::string StoreCompoundTaskMover::mover_name() {
	return "StoreCompoundTaskMover";
}

void StoreCompoundTaskMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaRestriction compound_task_behavior;
	compound_task_behavior.name( "compound_task_behavior" );
	compound_task_behavior.base_type( xs_string );
	compound_task_behavior.add_restriction( xsr_enumeration, "" );
	compound_task_behavior.add_restriction( xsr_enumeration, "prevent_repacking" );
	compound_task_behavior.add_restriction( xsr_enumeration, "restrict_to_repacking" );
	xsd.add_top_level_element( compound_task_behavior );


	XMLSchemaRestriction compound_task_mode;
	compound_task_mode.name( "compound_task_mode" );
	compound_task_mode.base_type( xs_string );
	compound_task_mode.add_restriction( xsr_enumeration, "packable" );
	compound_task_mode.add_restriction( xsr_enumeration, "designable" );
	xsd.add_top_level_element( compound_task_mode );


	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "task_name", xsct_pose_cached_task_operation, "The name by which the task operation to be cached will be identified after it is cached in the datacache of a Pose object." )
		//I'm calling these defaults since the empty string is a valid setting
		+ XMLSchemaAttribute::attribute_w_default( "true_behavior", "compound_task_behavior", "XRW TO DO", "")
		+ XMLSchemaAttribute::attribute_w_default( "false_behavior", "compound_task_behavior", "XRW TO DO",  "")
		+ XMLSchemaAttribute::attribute_w_default( "mode", "compound_task_mode", "XRW TO DO", "packable" )
		+ XMLSchemaAttribute::attribute_w_default( "invert", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "verbose", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "overwrite", xsct_rosetta_bool, "XRW TO DO", "false");

	AttributeList subtag_attributes;
	//All subtags have the same possible attributes
	protocols::rosetta_scripts::attributes_for_parse_task_operations( subtag_attributes );

	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_simple_subelement( "AND", subtag_attributes, "XRW_TODO" )
		.add_simple_subelement( "OR", subtag_attributes, "XRW_TODO" )
		.add_simple_subelement( "XOR", subtag_attributes, "XRW_TODO" )
		.add_simple_subelement( "NOR", subtag_attributes, "XRW_TODO" )
		.add_simple_subelement( "NAND", subtag_attributes, "XRW_TODO" )
		.add_simple_subelement( "ORNOT", subtag_attributes, "XRW_TODO" )
		.add_simple_subelement( "ANDNOT", subtag_attributes, "XRW_TODO" )
		.add_simple_subelement( "NOT", subtag_attributes, "XRW_TODO" );

	// TO DO: perhaps this is not the right function to call? -- also, delete this comment
	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW_TODO", attlist, subelements );
}

std::string StoreCompoundTaskMoverCreator::keyname() const {
	return StoreCompoundTaskMover::mover_name();
}

protocols::moves::MoverOP
StoreCompoundTaskMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new StoreCompoundTaskMover );
}

void StoreCompoundTaskMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	StoreCompoundTaskMover::provide_xml_schema( xsd );
}


} // task_operations
} // protocols

