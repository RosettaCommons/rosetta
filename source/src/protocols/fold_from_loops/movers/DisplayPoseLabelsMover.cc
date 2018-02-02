// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/DisplayPoseLabelsMover.cc
/// @brief Prints all the Labels of the Pose
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/movers/DisplayPoseLabelsMover.hh>
#include <protocols/fold_from_loops/movers/DisplayPoseLabelsMoverCreator.hh>
#include <protocols/fold_from_loops/utils/utils.hh>
#include <protocols/fold_from_loops/selectors/ConstraintResidueSelector.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>
#include <core/kinematics/MoveMap.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/movemap/util.hh>
#include <core/pose/selection.hh>
#include <core/kinematics/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/import_pose/import_pose.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.fold_from_loops.DisplayPoseLabelsMover" );

namespace protocols {
namespace fold_from_loops {
namespace movers {

DisplayPoseLabelsMover::DisplayPoseLabelsMover():
	protocols::moves::Mover( mover_name() ),
	title_width_( default_title_width() ),
	use_dssp_( default_use_dssp() ),
	movemap_factory_( /* NULL */ ),
	tasks_( /* NULL */ ),
	write_( default_write() )
{}

DisplayPoseLabelsMover::~DisplayPoseLabelsMover()= default;

void
DisplayPoseLabelsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	title_width( tag->getOption< core::Size >( "title_width", default_title_width() ) );
	use_dssp( tag->getOption< bool >( "use_dssp", default_use_dssp() ) );
	write( tag->getOption< bool >( "write", default_write() ) );
	movemap_factory( core::select::movemap::parse_movemap_factory( tag, data, "movemap_factory" ) );
	if ( tag->hasOption("task_operations") ) {
		tasks( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	}
}

void
DisplayPoseLabelsMover::apply( core::pose::Pose & pose )
{
	using namespace core::select::residue_selector;


	// SEQUENCE AND STRUCTURE
	std::string sequence  = pose.sequence();
	std::string structure = pose.secstruct();
	if ( use_dssp_ and structure.find_first_not_of("L")==std::string::npos ) {
		protocols::moves::DsspMover dssp;
		core::pose::Pose copy_pose = pose;
		dssp.apply( copy_pose );
		structure = copy_pose.secstruct();
	}
	print_data("SEQUENCE",  sequence,  true);
	print_data("STRUCTURE", structure, true);

	// LABELS
	std::map< std::string, std::string > labels = find_labels( pose );
	for ( auto & label : labels ) {
		print_data(label.first,  label.second,  false);
	}

	// CONSTRAINTS
	ResidueSelectorOP constraints( new selectors::ConstraintResidueSelector );
	// std::string constraint_repr = represent_residue_selector( constraints->apply( pose ) );
	print_data("CONSTRAINTS", represent_residue_selector( constraints->apply( pose ) ),  false);
	// bool constraint_show = constraint_repr.find("*")!=std::string::npos;

	// MoveMap
	if ( movemap_factory_ ) {
		print_movemap( pose );
	}

	// TaskOperators
	if ( tasks_ ) {
		print_task_operators( pose );
	}

	core::Size margin = 15;

	TR << std::left << std::setw( margin ) << "FOLDTREE";
	simple_visualize_fold_tree( pose.fold_tree(), TR );
	TR << std::endl;

	TR << pose.fold_tree() << std::endl;
	TR << pose.annotated_sequence() << std::endl;
	pose.pdb_info()->show(TR);

	if ( write_ ) {
		add_labels_as_remark( pose );
	}
}

bool
DisplayPoseLabelsMover::is_nubinitio_tree( core::kinematics::FoldTree const & fold_tree ) {
	core::Size root=0;
	for ( auto jump: fold_tree.get_jump_edges() ) {
		if ( root == 0 ) {
			root = jump.start();
		} else {
			if ( (core::Size) jump.start() != root ) {
				return false;
			}
		}
	}
	return true;
}

void
DisplayPoseLabelsMover::simple_visualize_fold_tree( core::kinematics::FoldTree const & fold_tree, std::ostream& out ) {
	if ( is_nubinitio_tree( fold_tree ) ) {
		for ( core::Size pos = 1; pos <= fold_tree.nres(); pos++ ) {
			bool special( false );
			if ( fold_tree.is_jump_point( pos ) ) {
				if ( fold_tree.is_root( pos ) ) {
					out << "R";
					special = true;
				} else {
					for ( core::Size jnr=1; jnr<=fold_tree.num_jump(); jnr++ ) {
						if ( (core::Size) fold_tree.jump_edge( jnr ).start() == pos || (core::Size)fold_tree.jump_edge( jnr ).stop() == pos ) {
							if ( special ) out << "/";
							out << jnr;
							special = true;
						}
					}
				}
			}
			if ( fold_tree.is_cutpoint( pos ) ) {
				if ( special ) out << "/";
				out << "C";
				special = true;
			}
			if ( !special ) out << "*";
		}
		out << std::endl;
	} else {
		core::kinematics::simple_visualize_fold_tree( fold_tree, out );
	}
}

void
DisplayPoseLabelsMover::print_movemap( core::pose::Pose const & pose )
{
	core::kinematics::MoveMapOP movemap = movemap_from_pose( pose );
	std::string bb  = "";
	std::string chi = "";
	for ( core::Size i = 1; i <= pose.size(); i++ ) {
		if ( movemap->get_bb ( i ) )  bb.append( "*" );
		else                bb.append( "." );
		if ( movemap->get_chi ( i ) ) chi.append( "*" );
		else              chi.append( "." );
	}
	print_data( "MVMP_BB",  bb, true );
	print_data( "MVMP_CHI", chi, true );
}

void
DisplayPoseLabelsMover::print_task_operators( core::pose::Pose const & pose )
{
	using namespace core::select::residue_selector;
	core::pack::task::PackerTaskOP task = tasks_->create_task_and_apply_taskoperations( pose );
	print_data("DESIGN", represent_residue_selector( task->designing_residues() ), false );
	print_data("REPACK", represent_residue_selector( task->repacking_residues() ), false );
}

std::map< std::string, std::string >
DisplayPoseLabelsMover::find_labels( core::pose::Pose const & pose ) const
{
	using namespace core::select::residue_selector;

	std::set< std::string > labels;
	for ( core::Size i = 1; i <= pose.size(); i++ ) {
		utility::vector1< std::string > rlabs = pose.pdb_info()->get_reslabels( i );
		for ( auto rlab : rlabs ) {
			labels.emplace( rlab );
		}
	}
	std::map< std::string, std::string > strlabels;
	ResiduePDBInfoHasLabelSelector residue_label;
	for ( auto label : labels ) {
		residue_label.set_label( label );
		strlabels[label] = represent_residue_selector( residue_label.apply( pose ) );
	}
	return strlabels;
}

std::map< std::string, utility::vector1< core::Size > >
DisplayPoseLabelsMover::find_labels2( core::pose::Pose const & pose ) const
{
	using namespace core::select::residue_selector;

	std::set< std::string > labels;
	for ( core::Size i = 1; i <= pose.size(); i++ ) {
		utility::vector1< std::string > rlabs = pose.pdb_info()->get_reslabels( i );
		for ( auto rlab : rlabs ) {
			labels.emplace( rlab );
		}
	}
	std::map< std::string, utility::vector1< core::Size > > strlabels;
	ResiduePDBInfoHasLabelSelector residue_label;
	for ( auto label : labels ) {
		residue_label.set_label( label );
		strlabels[label] = selection_positions( residue_label.apply( pose ) );
	}
	return strlabels;
}

void
DisplayPoseLabelsMover::add_labels_as_remark( core::pose::Pose & pose ) const
{
	using namespace core::select::residue_selector;
	std::map< std::string, utility::vector1< core::Size > > labels_raw = find_labels2( pose );
	utility::vector1< std::string > labels;
	for ( auto lab : labels_raw ) {
		auto first = std::begin(lab.second), last = std::end(lab.second);
		std::stringbuf buffer;
		std::ostream os (&buffer);
		os << lab.first << ":";
		while ( true )
				{
			auto mid = std::adjacent_find(first, last,
				[](int x, int y){ return x + 1 != y; });
			if ( mid == last ) {
				break;
			}
			if ( first == mid ) {
				os << *first << ",";
			} else {
				os << *first << "-" << *mid << ",";
			}
			first = ++mid;
		}
		if ( first == --last ) {
			os << *first;
		} else {
			os << *first << "-" << *last;
		}
		labels.push_back( buffer.str() );
	}
	std::ostringstream imploded;
	const char* const delim = ";";
	std::copy( labels.begin(), labels.end(), std::ostream_iterator<std::string>(imploded, delim) );
	std::string result_labs = imploded.str();
	result_labs.pop_back();
	core::pose::add_comment( pose, "LABELS", result_labs );
}

void
DisplayPoseLabelsMover::print_data( std::string id, std::string content, bool force )
{
	if ( force or content.find("*")!=std::string::npos ) {
		TR << std::left << std::setw( title_width_ ) << id;
		TR << content << std::endl;
	}
}

void
DisplayPoseLabelsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "title_width", xs_string,
		"Width assigned to print the label names to garantee alignment of the representation.",
		std::to_string( default_title_width() ) )
		+ XMLSchemaAttribute::attribute_w_default( "use_dssp", xsct_rosetta_bool,
		"Use DSSP if the Pose has no Secondary Structure Assigned (is all loop).", std::to_string( default_use_dssp() ) )
		+ XMLSchemaAttribute::attribute_w_default( "write", xsct_rosetta_bool,
		"Write the LABELS as a REMARK in the output structure or silent file "
		"(overwrites previous REMARKS with the title LABELS).", std::to_string( default_write() ) );

	protocols::rosetta_scripts::attributes_for_parse_task_operations_w_factory(attlist);

	core::select::movemap::attributes_for_parse_movemap_factory(attlist, "movemap_factory",
		"MoveMapFactory to print aligned to the residues");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Print the Pose's labels", attlist );
}

protocols::moves::MoverOP
DisplayPoseLabelsMover::clone() const
{
	return protocols::moves::MoverOP( new DisplayPoseLabelsMover( *this ) );
}

protocols::moves::MoverOP
DisplayPoseLabelsMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new DisplayPoseLabelsMover );
}

std::string DisplayPoseLabelsMover::get_name() const {
	return mover_name();
}

std::string DisplayPoseLabelsMover::mover_name() {
	return "DisplayPoseLabelsMover";
}

std::string DisplayPoseLabelsMoverCreator::keyname() const {
	return DisplayPoseLabelsMover::mover_name();
}

protocols::moves::MoverOP
DisplayPoseLabelsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DisplayPoseLabelsMover );
}

void DisplayPoseLabelsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DisplayPoseLabelsMover::provide_xml_schema( xsd );
}

}
} //protocols
} //fold_from_loops
