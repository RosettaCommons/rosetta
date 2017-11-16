// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Jacob Bale (balej@uw.edu)
#include <protocols/simple_filters/MutationsFilter.hh>
#include <protocols/simple_filters/MutationsFilterCreator.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/symmetry/util.hh>
#include <protocols/jd2/util.hh>
#include <ObjexxFCL/format.hh>
#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.MutationsFilter" );

/// @brief default ctor
MutationsFilter::MutationsFilter() :
	parent( "Mutations" ),
	task_factory_( /* NULL */ ),
	reference_pose_( /* NULL */ ),
	rate_threshold_( 0.0 ),
	mutation_threshold_( 100 ),
	mutations_( 0 ),
	verbose_( 0 ),
	packable_( 0 )
{}

core::pack::task::TaskFactoryOP
MutationsFilter::task_factory() const
{
	return task_factory_;
}

void
MutationsFilter::task_factory( core::pack::task::TaskFactoryOP task_factory )
{
	task_factory_ = task_factory;
}

core::Real
MutationsFilter::rate_threshold() const
{
	return( rate_threshold_ );
}

void
MutationsFilter::rate_threshold( core::Real const rate )
{
	runtime_assert( rate_threshold_ >= 0 && rate_threshold_ <= 1.0 );
	rate_threshold_ = rate;
}

core::pose::PoseCOP
MutationsFilter::reference_pose() const
{
	return reference_pose_;
}

void
MutationsFilter::reference_pose( core::pose::PoseCOP pose )
{
	reference_pose_ = pose;
}

void
MutationsFilter::reference_pose( core::pose::Pose const & pose )
{
	reference_pose_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( pose ) ) );
}

core::Size
MutationsFilter::mutation_threshold() const
{
	return( mutation_threshold_ );
}

void
MutationsFilter::mutation_threshold( core::Size const mut )
{
	mutation_threshold_ = mut;
}

bool
MutationsFilter::mutations() const
{
	return( mutations_ );
}

void
MutationsFilter::mutations( bool const muts )
{
	mutations_ = muts;
}

bool
MutationsFilter::verbose() const
{
	return( verbose_ );
}

void
MutationsFilter::verbose( bool const verb )
{
	verbose_ = verb;
}

bool
MutationsFilter::packable() const
{
	return( packable_ );
}

void
MutationsFilter::packable( bool const pack )
{
	packable_ = pack;
}

bool
MutationsFilter::write2pdb() const
{
	return( write2pdb_ );
}

void
MutationsFilter::write2pdb( bool const write )
{
	write2pdb_ = write;
}

bool
MutationsFilter::apply(core::pose::Pose const & pose ) const
{
	if ( mutations_ ) {
		core::Size const num_mutations( (core::Size) compute( pose, false ) );
		TR<<"The designed pose possesses "<<num_mutations<<" compared to the reference pose. ";
		if ( num_mutations <= mutation_threshold_ ) {
			TR<<"Success."<<std::endl;
			return true;
		} else {
			TR<<"Failing."<<std::endl;
			return false;
		}
	} else {
		core::Real const recovery_rate( compute( pose, false ) );
		TR<<"Sequence recovery rate evaluates to "<<recovery_rate<<". ";
		if ( recovery_rate <= rate_threshold_ ) {
			TR<<"Failing."<<std::endl;
			return false;
		} else {
			TR<<"Success."<<std::endl;
			return true;
		}
	}
}

core::Real
MutationsFilter::compute( core::pose::Pose const & pose, bool const & write ) const{
	runtime_assert( task_factory() != nullptr );
	runtime_assert( reference_pose() != nullptr );
	core::Size total_residue_ref;
	core::pose::Pose asym_ref_pose;
	if ( core::pose::symmetry::is_symmetric( *reference_pose() ) ) {
		core::pose::symmetry::extract_asymmetric_unit( *reference_pose(), asym_ref_pose);
		for ( core::Size i = 1; i <= asym_ref_pose.size(); ++i ) {
			if ( asym_ref_pose.residue_type(i).name() == "VRT" ) {
				asym_ref_pose.conformation().delete_residue_slow(asym_ref_pose.size());
			}
		}
		total_residue_ref = asym_ref_pose.size();
	} else {
		total_residue_ref = reference_pose()->size();
		asym_ref_pose = *reference_pose();
	}
	core::Size total_residue;
	core::pose::Pose asym_pose;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::pose::symmetry::extract_asymmetric_unit(pose, asym_pose);
		for ( core::Size i = 1; i <= asym_pose.size(); ++i ) {
			if ( asym_pose.residue_type(i).name() == "VRT" ) {
				asym_pose.conformation().delete_residue_slow(asym_pose.size());
			}
		}
		total_residue = asym_pose.size();
	} else {
		total_residue = pose.size();
		asym_pose = pose;
	}
	if ( total_residue_ref != total_residue ) {
		utility_exit_with_message( "Reference pose and current pose have a different number of residues" );
	}
	core::pack::task::PackerTaskOP packer_task( task_factory_->create_task_and_apply_taskoperations( pose ) );
	core::Size resi_count( 0 );
	core::Size mutation_count( 0 );
	//core::Size output_resi;
	std::map< core::Size, std::string > res_names1;
	std::map< core::Size, std::string > res_names2;
	for ( core::Size resi=1; resi<=total_residue; ++resi ) {
		if ( packer_task->being_designed( resi ) || (packer_task->being_packed( resi) && packable_) ) {
			resi_count++;
			res_names1.insert( std::make_pair( resi, asym_ref_pose.residue(resi).name3() ));
			res_names2.insert( std::make_pair( resi, asym_pose.residue(resi).name3() ));
			if ( asym_ref_pose.residue(resi).name3() != asym_pose.residue(resi).name3() ) {
				mutation_count++;
				if ( verbose_ ) {
					if ( basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
						TR << asym_ref_pose.residue(resi).name3() << resi << asym_pose.residue(resi).name3() << std::endl;
					} else {
						TR << asym_ref_pose.residue(resi).name3() << pose.pdb_info()->number( resi ) << asym_pose.residue(resi).name3() << std::endl;
					}
				}
			}
		}
	}
	// AMW: cppcheck notes that if this can be zero, then we divide by zero by it in the next line
	if ( !resi_count ) {
		TR.Warning << "No designable residues identified in pose. Are you sure you have set the correct task operations?"<<std::endl;
	}
	core::Real const rate( 1.0 - (core::Real) mutation_count / resi_count );
	TR<<"Your design mover mutated "<<mutation_count<<" positions out of "<<resi_count<<" designable positions. Sequence recovery is: "<<rate<<std::endl;
	if ( write ) {
		write_to_pdb( pose, res_names1, res_names2 );
	}
	if ( mutations_ ) {
		return( (core::Real) mutation_count );
	} else {
		return( rate );
	}
}

/// @brief Add each mutation to the output pdb if desired
void
MutationsFilter::write_to_pdb( core::pose::Pose const & pose, std::map< core::Size, std::string > const & res_names1, std::map< core::Size, std::string > const & res_names2 ) const {

	std::string user_name = this->get_user_defined_name();
	auto it_name1 = res_names1.begin();
	auto it_name2 = res_names2.begin();
	core::Size output_resi = it_name1->first;
	while ( it_name1 != res_names1.end() ) {
		if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
			output_resi = pose.pdb_info()->number( it_name1->first );
		}
		std::string output_string = "MutationsFilter " + user_name + ": " + it_name2->second + ObjexxFCL::string_of( output_resi ) + it_name1->second;
		protocols::jd2::add_string_to_current_job( output_string );
		++it_name1; ++it_name2;
	}

}

core::Real
MutationsFilter::report_sm( core::pose::Pose const & pose ) const
{
	return( compute( pose, false ) );
}

void
MutationsFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out<<"MutationsFilter returns "<<compute( pose, write2pdb() )<<std::endl;
}

void
MutationsFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{
	TR << "MutationsFilter"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	rate_threshold( tag->getOption< core::Real >( "rate_threshold", 0.0 ) );
	mutation_threshold( tag->getOption< core::Size >( "mutation_threshold", 100 ) );
	mutations( tag->getOption< bool >( "report_mutations", 0 ) );
	verbose( tag->getOption< bool >( "verbose", 0 ) );
	packable( tag->getOption< bool >( "packable", 0 ) );
	write2pdb( tag->getOption< bool >( "write2pdb", 0 ) );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ in::file::native ].user() ) {
		std::string const reference_pdb = option[ in::file::native ]();
		core::pose::PoseOP temp_pose( new core::pose::Pose );
		core::import_pose::pose_from_file( *temp_pose, reference_pdb , core::import_pose::PDB_file);
		reference_pose( temp_pose );
		TR<<"Using native pdb "<<reference_pdb<<" as reference.";
	} else {
		TR<<"Using starting pdb as reference. You could use -in::file::native to specify a different pdb for reference";
		reference_pose( pose );
	}
	TR<<std::endl;
}

protocols::filters::FilterOP
MutationsFilter::fresh_instance() const{
	return protocols::filters::FilterOP( new MutationsFilter() );
}

MutationsFilter::~MutationsFilter()= default;


protocols::filters::FilterOP
MutationsFilter::clone() const{
	return protocols::filters::FilterOP( new MutationsFilter( *this ) );
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP MutationsFilterCreator::create_filter() const { return protocols::filters::FilterOP( new MutationsFilter ); }

// XRW TEMP std::string
// XRW TEMP MutationsFilterCreator::keyname() const { return "Mutations"; }

std::string MutationsFilter::name() const {
	return class_name();
}

std::string MutationsFilter::class_name() {
	return "Mutations";
}

void MutationsFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default("rate_threshold", xsct_real, "Lower cutoff for the acceptable recovery rate for a passing design. Will fail if actual rate is below this threshold.", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("mutation_threshold", xsct_non_negative_integer, "Upper cutoff for the number of mutations for an acceptable design. Only matters if report_mutations is set to true.", "100")
		+ XMLSchemaAttribute::attribute_w_default("report_mutations", xsct_rosetta_bool, "Defaults to false. If set to true, then will act as a filter for the number of mutations rather than the rate.", "0")
		+ XMLSchemaAttribute::attribute_w_default("verbose", xsct_rosetta_bool, "Defaults to false. If set to true, then will output the mutated positions and identities to the tracer.", "0")
		+ XMLSchemaAttribute::attribute_w_default("packable", xsct_rosetta_bool, "Defaults to false. If set to true, then will also consider mutations at packable positions in addition to designable positions.", "0")
		+ XMLSchemaAttribute::attribute_w_default("write2pdb", xsct_rosetta_bool, "Defaults to false. If set to true, then will output the mutated positions and identities to the output pdb.", "0");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Determines mutated residues in current pose as compared to a reference pose", attlist );
}

std::string MutationsFilterCreator::keyname() const {
	return MutationsFilter::class_name();
}

protocols::filters::FilterOP
MutationsFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new MutationsFilter );
}

void MutationsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MutationsFilter::provide_xml_schema( xsd );
}



} // simple_filters
} // protocols
