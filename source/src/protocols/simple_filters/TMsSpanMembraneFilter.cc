// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/TMsSpanMembraneFilter.cc
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
#include <protocols/simple_filters/TMsSpanMembraneFilter.hh>
#include <protocols/simple_filters/TMsSpanMembraneFilterCreator.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

#include <iostream>
#include <fstream>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <protocols/membrane/util.hh>
#include <core/conformation/membrane/Exceptions.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.TMsSpanMembraneFilter" );

protocols::filters::FilterOP
TMsSpanMembraneFilterCreator::create_filter() const { return protocols::filters::FilterOP( new TMsSpanMembraneFilter ); }

std::string
TMsSpanMembraneFilterCreator::keyname() const { return "TMsSpanMembrane"; }

TMsSpanMembraneFilter::~TMsSpanMembraneFilter(){}

void
TMsSpanMembraneFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	threshold_ = tag->getOption< core::Real >( "threshold", 0.5 );
	output_ = tag->getOption< std::string >( "output_file", "" );
	min_distance_ = tag->getOption< core::Real >( "min_distance", 20.0 );
	flank_ = tag->getOption< core::Size >( "flank", 1 );
	required_dist_ = tag->getOption< core::Size >( "required_distance", 10 );
}

bool
TMsSpanMembraneFilter::apply( core::pose::Pose const & pose ) const {
	core::Real tms_agree = compute( pose );
	bool const status = tms_agree >= threshold_;
	return status;
}

void
TMsSpanMembraneFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real tms_agree = compute( pose );
	out << "matching pose topology to span results in " << tms_agree << std::endl;
}

core::Real
TMsSpanMembraneFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real tms_agree = compute( pose );
	return( tms_agree );
}

core::Real
TMsSpanMembraneFilter::compute( core::pose::Pose const & pose ) const {

	utility::vector1< utility::pointer::shared_ptr< core::conformation::membrane::Span > > const spans = pose.conformation().membrane_info()->spanning_topology()->get_spans();

	// iterate over all spans in the span file, and check if they uphold two things:
	// 1. is the absolute distance (on the Z axis) bigger than a minimal threshold
	// 2. is the span's orientation upheld
	for ( core::Size span_i=1; span_i <= spans.size(); span_i++ ) {

		// get the Z of the span's starting residue
		core::conformation::Residue const & start_residue = pose.residue(spans[ span_i ]->start());
		core::Real start_z = pose.conformation().membrane_info()->residue_z_position( pose.conformation(), start_residue.seqpos() );

		// get the Z of the span's ending residue
		core::conformation::Residue const & end_residue = pose.residue(spans[ span_i ]->end());
		core::Real end_z = pose.conformation().membrane_info()->residue_z_position( pose.conformation(), end_residue.seqpos() );

		// if the Z coordinate distance between the TM edges is smaller than min_distance_, fail the pose
		core::Real z_distance = std::abs( start_z - end_z );
		if ( z_distance < min_distance_ ) {
			TR << "distance between residue " << spans[ span_i  ]->start() << " and residue " << spans[ span_i  ]->end() << " is " << z_distance << std::endl;
			return( 0 );
		}

		// if span edges are not close enough to membrane edges, fail the structure
		if ( std::abs(start_z) < required_dist_ || std::abs(end_z) < required_dist_ ) {
			TR << "either start or end (or both) is not close enough to membrane edge, required distance is " << required_dist_ << std::endl;
			TR << "start of span is at " << start_z << " end of span is at " << end_z << std::endl;
			return( 0 );
		}

		// orientation=1 means from in to out
		// check if the span's orientation is correct
		if ( spans[ span_i ]->orientation() == 1 ) { // means it's in -> out
			if ( start_z > 0 || end_z < 0 ) {
				TR << "although span num " << span_i << " is in->out " << " either start is out or end is in, start(z)=" << start_z << ", end(z)=" << end_z << std::endl;
				return( 0 );
			}
		} else { // means out -> in
			if ( start_z < 0 || end_z > 0 ) {
				TR << "although span num " << span_i << " is out->in " << " either start is in or end is out, start(z)=" << start_z << ", end(z)=" << end_z << std::endl;
				return( 0 );
			}
		}
	}

	// test if residues that shouldn't be in the membrnae are
	for ( core::Size res_i=1; res_i <= pose.size(); res_i++ ) {
		core::conformation::Residue const & residue_i = pose.residue(res_i);
		if ( ! residue_i.is_protein() || residue_i.is_virtual_residue() ) continue;
		core::Real z = pose.conformation().membrane_info()->residue_z_position( pose.conformation(), residue_i.seqpos() );
		if ( -15 <= z && z <= 15 ) {
			// residue is in the membrane
			bool found_in_span = false;
			for ( core::Size span_i=1; span_i <= spans.size(); span_i++ ) {
				if ( spans[ span_i ]->start()-flank_ <= res_i && spans[ span_i ]->end()+flank_ ) {
					found_in_span = true;
					break;
				}
			}
			if ( !found_in_span ) {
				TR << "residue num " << res_i << " is in the membrane, but is not in any span" << std::endl;
				return( 0 );
			}
		}
	}
	return( 1 );
}

std::string TMsSpanMembraneFilter::name() const {
	return class_name();
}

void TMsSpanMembraneFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold" , xsct_real , "leave 0.5 for boolean behaviour" , "0.5" )
		+ XMLSchemaAttribute::attribute_w_default( "output_file" , xs_string , "file to which print output" , "" )
		+ XMLSchemaAttribute::attribute_w_default( "min_distance" , xsct_real , "minimal distance between TMH edges over the Z axis (membrane depth)" , "20" )
		+ XMLSchemaAttribute::attribute_w_default( "flank" , xsct_non_negative_integer , "by how many residues 'enlarge' the TMH to allow a more fuzzy test" , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "required_distance" , xsct_non_negative_integer , "the distance the edges should be from the membrane edge" , "10" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Computes the overall sasa of the pose.", attlist );
}

void TMsSpanMembraneFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	TMsSpanMembraneFilter::provide_xml_schema( xsd );
}

std::string TMsSpanMembraneFilter::class_name() {
	return "TMsSpanMembrane";
}

}
}
