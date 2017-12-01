// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/MembAccesResidueLipophilicityFilter.cc
/// @brief return the matching (%) between span file determined topology and of the actual pose
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)
// Project Headers

#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <map>
#include <numeric/random/random.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/MembAccesResidueLipophilicityFilter.hh>
#include <protocols/simple_filters/MembAccesResidueLipophilicityFilterCreator.hh>
#include <protocols/simple_moves/ddG.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <iostream>
#include <fstream>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <protocols/membrane/util.hh>
#include <core/conformation/membrane/Exceptions.hh>
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/scoring/sasa/util.hh>
#include <core/scoring/membrane/MPSpanInsertionEnergy.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.MembAccesResidueLipophilicityFilter" );

protocols::filters::FilterOP
MembAccesResidueLipophilicityFilterCreator::create_filter() const { return protocols::filters::FilterOP( new MembAccesResidueLipophilicityFilter ); }

std::string
MembAccesResidueLipophilicityFilterCreator::keyname() const { return "MembAccesResidueLipophilicity"; }

MembAccesResidueLipophilicityFilter::~MembAccesResidueLipophilicityFilter(){}

void
MembAccesResidueLipophilicityFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data );
	threshold_ = tag->getOption< core::Real >( "threshold", 0.0 );
	verbose_ = tag->getOption< bool >( "verbose", false );
	ignore_burial_ = tag->getOption< bool >( "ignore_burial", false );
}

bool
MembAccesResidueLipophilicityFilter::apply( core::pose::Pose const & pose ) const {
	core::Real total_dg = compute( pose );
	bool const status = total_dg <= threshold_;
	if ( status ) {
		TR << "total ddG is " << total_dg << " which passes the threshold " << threshold_ << std::endl;
	} else {
		TR << "total ddG is " << total_dg << " and fails with threshold " << threshold_ << std::endl;
	}
	return status;
}

void
MembAccesResidueLipophilicityFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real total_dg = compute( pose );
	out << "total ddG is " << total_dg << std::endl;
}

core::Real
MembAccesResidueLipophilicityFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real total_dg = compute( pose );
	TR << "total ddG is " << total_dg << std::endl;
	return( total_dg );
}

core::Real
MembAccesResidueLipophilicityFilter::compute( core::pose::Pose const & pose ) const {
	runtime_assert( task_factory_ != nullptr  );
	core::pack::task::PackerTaskCOP packer_task( task_factory_->create_task_and_apply_taskoperations( pose ) );
	TR << "calculating dG of insertion in kcal/mol per resiude " << std::endl;
	if ( ignore_burial_ ) TR << "IGNORING BURIAL" << std::endl;
	utility::vector1< core::Real > sasas = core::scoring::sasa::rel_per_res_sc_sasa( pose  );
	core::scoring::membrane::MPSpanInsertionEnergy mpsie = core::scoring::membrane::MPSpanInsertionEnergy();
	core::Real dG( 0.0 );
	core::Real z( 0.0 );

	core::Real total_dg( 0.0 );
	if ( verbose_ ) TR << "res\tz\tdG\tsasa\tdG*sasa" << std::endl;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		core::conformation::Residue rsd = pose.residue( i  );
		if ( ! rsd.is_protein() ) continue;
		if ( ! packer_task->being_packed( i ) ) {
			TR << "skipping " << i << rsd.aa() << std::endl;
			continue;
		} // skip res by task operations
		if ( verbose_ ) TR << "TO " << i << rsd.aa() << " " << packer_task->being_packed( i ) << " " << packer_task->being_designed(i) << std::endl;
		z = pose.conformation().membrane_info()->residue_z_position( pose.conformation(),i );
		dG = mpsie.spline_by_z( static_cast< char > (oneletter_code_from_aa(rsd.aa())),  z);
		if ( verbose_ ) TR << i << static_cast< char > (oneletter_code_from_aa(rsd.aa())) << "\t" << z << "\t" << std::setprecision (2) << std::fixed << dG << "\t" << sasas[ i ] << "\t" << dG*sasas[ i ] << std::endl;
		total_dg += (( ignore_burial_ ) ? dG : dG * sasas[ i ]);
	}
	return( total_dg );
}

void MembAccesResidueLipophilicityFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "if score lower than threshold, pass" , "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "verbose" , xsct_rosetta_bool , "whehter to print a table of dG and SASA" , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "ignore_burial" , xsct_rosetta_bool , "disregard the SASA in the calculation" , "false" );
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist ) ; //Only report the SASA for those residues specified as packable for the given taskoperations. If not specified, compute over all residues.
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Computes the dG from dsTbL profiles of every residue in every span, and multiply it by its relative SASA", attlist );
}

void MembAccesResidueLipophilicityFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MembAccesResidueLipophilicityFilter::provide_xml_schema( xsd );
}

std::string MembAccesResidueLipophilicityFilter::class_name() {
	return "MembAccesResidueLipophilicity";
}

std::string MembAccesResidueLipophilicityFilter::name() const {
	return class_name();
}

}
}
