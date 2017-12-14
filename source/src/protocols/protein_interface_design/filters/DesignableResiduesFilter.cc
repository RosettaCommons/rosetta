// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/DesignableResiduesFilter.hh>
#include <protocols/protein_interface_design/filters/DesignableResiduesFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

// JBB 120425
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/pose/symmetry/util.hh>
#include <ObjexxFCL/format.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace filters {

static basic::Tracer TR( "protocols.protein_interface_design.filters.DesignableResiduesFilter" );
static basic::Tracer TR_pymol( "protocols.protein_interface_design.filters.DesignableResiduesFilter_pymol" );

/// @brief default ctor
DesignableResiduesFilter::DesignableResiduesFilter() :
	parent( "DesignableResidues" ),
	task_factory_( /* NULL */ ),
	lower_threshold_( 0 ),
	upper_threshold_( 1000 ),
	packable_( false ),
	designable_( true )
{}

core::Size
DesignableResiduesFilter::lower_threshold() const{
	return lower_threshold_;
}

core::Size
DesignableResiduesFilter::upper_threshold() const{
	return upper_threshold_;
}

bool
DesignableResiduesFilter::packable() const{
	return packable_;
}

bool
DesignableResiduesFilter::designable() const{
	return designable_;
}

void
DesignableResiduesFilter::lower_threshold( core::Size const l ){
	lower_threshold_ = l;
}

void
DesignableResiduesFilter::upper_threshold( core::Size const u ){
	upper_threshold_ = u;
}

void
DesignableResiduesFilter::designable( bool const d ){
	designable_ = d;
}

void
DesignableResiduesFilter::packable( bool const p ){
	packable_ = p;
}

core::pack::task::TaskFactoryOP
DesignableResiduesFilter::task_factory() const
{
	return task_factory_;
}

void
DesignableResiduesFilter::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

bool
DesignableResiduesFilter::apply(core::pose::Pose const & ) const
{
	return true;
}

core::Size
DesignableResiduesFilter::compute( core::pose::Pose const & pose ) const{
	runtime_assert( task_factory() != nullptr );
	runtime_assert( packable() || designable() );
	core::pack::task::PackerTaskCOP packer_task( task_factory()->create_task_and_apply_taskoperations( pose ) );
	core::Size total_residue;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
		total_residue = symm_info->num_independent_residues();
	} else {
		total_residue = pose.size();
	}
	core::Size design_pos = 0;
	if ( designable() ) {
		std::string select_design_pos("select design_positions, resi ");
		TR<<"Designable residues:"<<std::endl;
		bool first_pass( true );
		for ( core::Size resi=1; resi<=total_residue; ++resi ) {
			if ( packer_task->being_designed( resi ) ) {
				if ( !first_pass ) {
					select_design_pos.append( "+" );
				}
				TR<<pose.residue( resi ).name3()<<" "<< pose.pdb_info()->number( resi )<<pose.pdb_info()->chain( resi )<<std::endl;
				design_pos++;
				select_design_pos.append(ObjexxFCL::string_of(resi));
				first_pass = false;
			}
		}
		TR<<"Number of design positions: "<<design_pos<<std::endl;
		TR_pymol<<select_design_pos<<std::endl;
	}
	core::Size packable_pos = 0;
	if ( packable() ) {
		std::string select_packable_pos("select repackable_positions, resi ");
		TR<<"Repackable residues:"<<std::endl;
		for ( core::Size resi=1; resi<=total_residue; ++resi ) {
			if ( packer_task->being_packed( resi ) ) {
				TR<<pose.residue( resi ).name3()<<" "<<pose.pdb_info()->number( resi )<<pose.pdb_info()->chain( resi )<<std::endl;
				packable_pos++;
				select_packable_pos.append(ObjexxFCL::string_of(resi) + "+");
			}
		}
		TR<<"Number of repackable positions: "<<packable_pos<<std::endl;
		TR<<select_packable_pos<<std::endl;
	}
	if ( designable() && !packable() ) {
		return( design_pos );
	} else {
		return( packable_pos);
	}
}

core::Real
DesignableResiduesFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Size design_pos(compute( pose ));
	return( design_pos );
}

void
DesignableResiduesFilter::report( std::ostream & , core::pose::Pose const & ) const
{
}

void
DesignableResiduesFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	TR << "DesignableResiduesFilter"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	lower_threshold( tag->getOption< core::Size >( "lower_cutoff", 0 ) );
	upper_threshold( tag->getOption< core::Size >( "upper_cutoff", 1000 ) );
	packable( tag->getOption< bool >( "packable", false ) );
	designable( tag->getOption< bool >( "designable", true ) );
	runtime_assert( designable() || packable() );
	TR<<"with options designable: "<<designable()<<", repackable: "<<packable()<<", lower_cutoff: "<<lower_threshold()<<", and upper_cutoff: "<<upper_threshold()<<std::endl;
}

protocols::filters::FilterOP
DesignableResiduesFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new DesignableResiduesFilter() );
}

DesignableResiduesFilter::~DesignableResiduesFilter()= default;

protocols::filters::FilterOP
DesignableResiduesFilter::clone() const{
	return protocols::filters::FilterOP( new DesignableResiduesFilter( *this ) );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP DesignableResiduesFilterCreator::create_filter() const { return protocols::filters::FilterOP( new DesignableResiduesFilter ); }

// XRW TEMP std::string
// XRW TEMP DesignableResiduesFilterCreator::keyname() const { return "DesignableResidues"; }

std::string DesignableResiduesFilter::name() const {
	return class_name();
}

std::string DesignableResiduesFilter::class_name() {
	return "DesignableResidues";
}

void DesignableResiduesFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "lower_cutoff", xsct_non_negative_integer, "This option presently has no effect.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "upper_cutoff", xsct_non_negative_integer, "This option presently has no effect.", "1000" )
		+ XMLSchemaAttribute::attribute_w_default( "packable", xsct_rosetta_bool, "Determine all packable positions", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "designable", xsct_rosetta_bool, "Determine all designable positions", "0" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Filters based on minimum and maximum number of designable residues allowed; useful for automatic interface detection", attlist );
}

std::string DesignableResiduesFilterCreator::keyname() const {
	return DesignableResiduesFilter::class_name();
}

protocols::filters::FilterOP
DesignableResiduesFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new DesignableResiduesFilter );
}

void DesignableResiduesFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DesignableResiduesFilter::provide_xml_schema( xsd );
}


} // filters
} // protein_interface_design
} // protocols
