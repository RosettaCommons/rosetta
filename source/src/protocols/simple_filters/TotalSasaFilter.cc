// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/TotalSasaFilter.cc
/// @brief
/// @details (Based on InterfaceSasaFilter)
/// @author Rocco Moretti (rmoretti@u.washington.edu)

#include <protocols/simple_filters/TotalSasaFilter.hh>
#include <protocols/simple_filters/TotalSasaFilterCreator.hh>

#include <core/pose/Pose.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <basic/MetricValue.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <core/types.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.TotalSasaFilter" );




TotalSasaFilter::TotalSasaFilter() :
	Filter( "TotalSasa" ),
	lower_threshold_( 0.0 ),
	upper_threshold_(100000000000.0),
	hydrophobic_( false ),
	polar_( false ),
	report_per_residue_sasa_(false)
{}

TotalSasaFilter::TotalSasaFilter( core::Real const lower_threshold, bool const hydrophobic/*=false*/, bool const polar/*=false*/, core::Real upper_threshold, bool per_residue_sasa/*=false*/ ) :
	Filter( "TotalSasa" ),
	lower_threshold_( lower_threshold ),
	upper_threshold_(upper_threshold),
	hydrophobic_( hydrophobic ),
	polar_( polar ),
	report_per_residue_sasa_(per_residue_sasa)
{}

TotalSasaFilter::~TotalSasaFilter()= default;

filters::FilterOP
TotalSasaFilter::clone() const{
	return utility::pointer::make_shared< TotalSasaFilter >( *this );
}

filters::FilterOP
TotalSasaFilter::fresh_instance() const{
	return utility::pointer::make_shared< TotalSasaFilter >();
}

core::pack::task::TaskFactoryOP
TotalSasaFilter::task_factory() {
	return taskfactory_;
}
void
TotalSasaFilter::task_factory(core::pack::task::TaskFactoryOP task_factory) {
	taskfactory_ = task_factory;
}

void
TotalSasaFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	lower_threshold_ = tag->getOption<core::Real>( "threshold", 800.0 );
	upper_threshold_ = tag->getOption<core::Real>( "upper_threshold", 1000000.0 );

	hydrophobic_ = tag->getOption<bool>( "hydrophobic", false );
	polar_ = tag->getOption<bool>( "polar", false );
	report_per_residue_sasa_ = tag->getOption<bool>( "report_per_residue_sasa", false );

	if ( tag->hasOption("task_operations") ) {
		task_factory( protocols::rosetta_scripts::parse_task_operations(tag, data) );
	} else {
		task_factory( nullptr ); // No task factory - use all residues.
	}

	if ( polar_ && hydrophobic_ ) {
		TR.Error << "Polar and hydrophobic both flags specified in TotalSasa filter: " << tag << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Polar and hydrophobic flags specified in TotalSasa filter." );
	}

	TR.Debug << "Parsed TotalSasa Filter: <TotalSasa" <<
		" threshold=" << lower_threshold_ <<
		" upper_threshold=" << upper_threshold_ <<
		" hydrophobic=" << hydrophobic_ <<
		" polar=" << polar_ <<
		"report_per_residue_sasa_=" << report_per_residue_sasa_ <<
		" />" << std::endl;
}

bool
TotalSasaFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ) );

	TR<<"sasa is "<<sasa<<". ";
	if ( sasa >= lower_threshold_ && sasa <= upper_threshold_ ) {
		TR<<"passing." <<std::endl;
		return true;
	} else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
TotalSasaFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ));
	out<<"Sasa= "<< sasa<<'\n';

	if ( report_per_residue_sasa_ ) {

		std::string threshold_string = "The following residues have sasa within the desired threshold range: \n";

		basic::MetricValue< utility::vector1< core::Real > > residue_sasa;
		pose.metric( "sasa", "residue_sasa", residue_sasa ); //this does not represent a new calculation since pose metric calculators are smart enough to figure it out

		runtime_assert( pose.size() == (residue_sasa.value()).size() );
		for ( core::Size i = 1; i<=pose.size(); ++i ) {
			std::string res_chain = pose.pdb_info()->chain(i);
			int res_pdbnum = pose.pdb_info()->number(i);
			core::Real this_sasa = residue_sasa.value()[i];
			out << pose.residue( i ).name3() << res_pdbnum << " " << res_chain << " : " << this_sasa << '\n';
			if ( (this_sasa > lower_threshold_) && (this_sasa < upper_threshold_) ) {
				threshold_string = threshold_string + utility::to_string(res_pdbnum) + "+";
			}
		}
		threshold_string = threshold_string + "\n";
		out << threshold_string;
	}
}

core::Real
TotalSasaFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const sasa( compute( pose ));
	return( sasa );
}

core::Real
TotalSasaFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::pose::metrics;
	using basic::MetricValue;
	using namespace core;
	using namespace protocols::moves;

	if ( !core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "sasa" ) ) {
		utility_exit_with_message("Must define core::pose::metrics::simple_calculators::SasaCalculatorLegacy with name 'sasa' in order to use TotalSasaFilter.");
	}

	runtime_assert( ! (hydrophobic_ && polar_) );
	if ( !hydrophobic_ && !polar_ && !taskfactory_ ) {
		// Complete Sasa
		MetricValue< core::Real > mv_sasa;

		pose.metric( "sasa", "total_sasa", mv_sasa);
		core::Real const total_sasa( mv_sasa.value() );
		TR.Debug << "Sasa: " << total_sasa << std::endl;

		return( total_sasa );
	} else if ( !hydrophobic_ && !polar_ ) {
		// Per residue Sasa
		debug_assert( taskfactory_ );
		core::Real sasa( 0.0 );
		core::pack::task::PackerTaskOP task( taskfactory_->create_task_and_apply_taskoperations(pose) );
		MetricValue< utility::vector1< core::Real > > residue_sasa;
		pose.metric( "sasa", "residue_sasa", residue_sasa );
		for ( core::Size ii(1); ii <= pose.size(); ++ii ) {
			if ( task->being_packed(ii) ) {
				sasa += residue_sasa.value()[ii];
			}
		}
		return sasa;
	} else {
		//Split sasa
		core::pack::task::PackerTaskOP task;
		if ( taskfactory_ ) {
			task = taskfactory_->create_task_and_apply_taskoperations(pose);
		}
		MetricValue< id::AtomID_Map< core::Real > > atom_sasa;
		pose.metric( "sasa", "atom_sasa", atom_sasa );
		core::Real polar_sasa( 0.0 ), hydrophobic_sasa( 0.0 );
		for ( core::Size pos(1); pos<=pose.size(); ++pos ) {
			core::Real pos_polar_sasa( 0.0 ), pos_hydrophobic_sasa( 0.0 );

			if ( task && ! task->being_packed(pos) ) { continue; }
			for ( core::Size atomi( 1 ); atomi <= atom_sasa.value().n_atom( pos ); ++atomi ) {
				core::Real const atomi_sasa( atom_sasa.value()( pos, atomi ) );
				core::conformation::Residue const pos_rsd( pose.residue( pos ) );
				core::chemical::AtomType const atom_type( pos_rsd.atom_type( atomi ) );
				bool const is_polar( atom_type.is_donor() || atom_type.is_acceptor() || atom_type.is_polar_hydrogen() );
				if ( is_polar ) {
					polar_sasa += atomi_sasa;
					pos_polar_sasa += atomi_sasa;
				} else {
					hydrophobic_sasa += atomi_sasa;
					pos_hydrophobic_sasa += atomi_sasa;
				}
			}

			if ( report_per_residue_sasa_ && (polar_ || hydrophobic_) ) {
				std::string res_chain = pose.pdb_info()->chain(pos);
				int res_pdbnum = pose.pdb_info()->number(pos);

				TR << pose.residue( pos ).name3() << res_pdbnum << " " << res_chain << " ";
				if ( polar_ ) { TR << "POLAR SASA : " << pos_polar_sasa << std::endl; }
				if ( hydrophobic_ ) { TR << "HYDROPHOBIC SASA : " << pos_hydrophobic_sasa << std::endl; }
			}


		}


		if ( hydrophobic_ ) return hydrophobic_sasa;
		if ( polar_ ) return polar_sasa;
	}
	utility_exit_with_message( "Execution should never have reached this point." );
	return( 0 );
}

std::string TotalSasaFilter::name() const {
	return class_name();
}

std::string TotalSasaFilter::class_name() {
	return "TotalSasa";
}

void TotalSasaFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "If it is **higher** than threshold, it passes. " , "800" )
		+ XMLSchemaAttribute::attribute_w_default( "upper_threshold" , xsct_real , "Fails if it is above the upper_threshold." , "1000000" )
		+ XMLSchemaAttribute::attribute_w_default( "hydrophobic" , xsct_rosetta_bool , "Compute hydrophobic-only SASA." , "false" )
		// Hydrophobic/polar are computed by discriminating each atom into polar (acceptor/donor or polar hydrogen)
		// or hydrophobic (all else) and summing the SASA over each category.
		+ XMLSchemaAttribute::attribute_w_default( "polar" , xsct_rosetta_bool , "Compute polar_only SASA." , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "report_per_residue_sasa" , xsct_rosetta_bool , "Add the per-residue SASA to the tracer output." , "false" ) ;

	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist ) ; //Only report the SASA for those residues specified as packable for the given taskoperations. If not specified, compute over all residues.

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Computes the overall sasa of the pose.", attlist );
}

std::string TotalSasaFilterCreator::keyname() const {
	return TotalSasaFilter::class_name();
}

protocols::filters::FilterOP
TotalSasaFilterCreator::create_filter() const {
	return utility::pointer::make_shared< TotalSasaFilter >();
}

void TotalSasaFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TotalSasaFilter::provide_xml_schema( xsd );
}


}
}
