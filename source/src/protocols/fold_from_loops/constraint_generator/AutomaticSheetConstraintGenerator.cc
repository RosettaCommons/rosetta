// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/AutomaticSheetConstraintGenerator.cc
/// @brief Finds angle constraints between the paired residues in betas
/// @author Jaume Bonet (jaume.bonet@gmail.com)

// Unit headers
#include <protocols/fold_from_loops/constraint_generator/AutomaticSheetConstraintGenerator.hh>
#include <protocols/fold_from_loops/constraint_generator/AutomaticSheetConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/AtomPairConstraintGenerator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>
#include <protocols/constraint_generator/util.hh>

// Core headers
#include <core/id/SequenceMapping.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/func/SOGFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/SecondaryStructureSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/func/Func.fwd.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/constants.hh>

static basic::Tracer TR( "protocols.fold_from_loops.AutomaticSheetConstraintGenerator" );

namespace protocols {
namespace fold_from_loops {
namespace constraint_generator {

protocols::constraint_generator::ConstraintGeneratorOP
AutomaticSheetConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new AutomaticSheetConstraintGenerator );
}

std::string
AutomaticSheetConstraintGeneratorCreator::keyname() const
{
	return AutomaticSheetConstraintGenerator::class_name();
}

AutomaticSheetConstraintGenerator::AutomaticSheetConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( AutomaticSheetConstraintGenerator::class_name() ),
	sd_( default_sd() ),
	distance_( default_distance() ),
	angle_tolerance_( 0.35 ),
	bb_dihedral_tolerance_( 0.52 ),
	weight_( 1.0 )
{}

AutomaticSheetConstraintGenerator::~AutomaticSheetConstraintGenerator() = default;

protocols::constraint_generator::ConstraintGeneratorOP
AutomaticSheetConstraintGenerator::clone() const
{
	return AutomaticSheetConstraintGeneratorOP( new AutomaticSheetConstraintGenerator( *this ) );
}

void
AutomaticSheetConstraintGenerator::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{

	sd( tag->getOption< core::Real >( "sd", default_sd() ) );
	distance( tag->getOption< core::Real >( "distance", default_distance() ) );

}

core::scoring::constraints::ConstraintCOPs
AutomaticSheetConstraintGenerator::apply( core::pose::Pose const & pose ) const
{
	using namespace protocols::constraint_generator;
	using namespace core::select::residue_selector;
	using namespace numeric::conversions;
	using namespace numeric;
	using namespace core::scoring::constraints;

	core::scoring::constraints::ConstraintCOPs constraints;
	SecondaryStructureSelector structure("E");
	structure.set_minE(2);
	ResidueRanges ranges;
	ResidueRanges selected;
	ranges.from_subset( structure.apply( pose ) );
	TR << represent_residue_selector( structure.apply( pose ) ) << std::endl;
	for ( core::Size i=1; i <= ranges.size(); ++i ) {
		for ( core::Size ii=ranges[i].start(); ii <= ranges[i].stop(); ++ii ) {
			for ( core::Size j=i+1; j <= ranges.size(); ++j ) {
				core::Real min  = 1000;
				core::Size sele = 0;
				for ( core::Size jj=ranges[j].start(); jj <= ranges[j].stop(); ++jj ) {
					core::Real const dist( pose.residue( ii ).xyz( "CA" ).distance( pose.residue( jj ).xyz( "CA" ) ) );
					if ( dist < distance_ and dist < min ) {
						sele = jj;
					}
				}
				if ( sele != 0 ) {
					ResidueRange selepair(ii, sele);
					selected.push_back( selepair );
				}
			}
		}
	}
	for ( auto i : selected ) {
		core::Real const dist( pose.residue( i.start() ).xyz( "CA" ).distance( pose.residue( i.stop() ).xyz( "CA" ) ) );
		core::Real angle1( angle_of(
			pose.residue( i.start() ).xyz( "N" ),
			pose.residue( i.start() ).xyz( "C" ),
			pose.residue( i.stop() ).xyz( "N" )
			) );
		core::Real angle2( angle_of(
			pose.residue( i.stop() ).xyz( "N" ),
			pose.residue( i.stop() ).xyz( "C" ),
			pose.residue( i.start() ).xyz( "N" )
			) );
		core::Real angle3( angle_of(
			pose.residue( i.start() ).xyz( "N" ),
			pose.residue( i.start() ).xyz( "C" ),
			pose.residue( i.stop() ).xyz( "C" )
			) );
		core::Real angle4( angle_of(
			pose.residue( i.stop() ).xyz( "N" ),
			pose.residue( i.stop() ).xyz( "C" ),
			pose.residue( i.start() ).xyz( "C" )
			) );
		TR << i.start() << " " << i.stop() << ": " << dist << "<->" << angle1 << "<->" << angle2 << std::endl;
		core::id::AtomID const resi_n( pose.residue_type( i.start() ).atom_index( "N" ), i.start() );
		core::id::AtomID const resi_c( pose.residue_type( i.start() ).atom_index( "C" ), i.start() );
		core::id::AtomID const resi_o( pose.residue_type( i.start() ).atom_index( "O" ), i.start() );
		core::id::AtomID const resj_n( pose.residue_type( i.stop()  ).atom_index( "N" ), i.stop()  );
		core::id::AtomID const resj_c( pose.residue_type( i.stop()  ).atom_index( "C" ), i.stop()  );
		core::id::AtomID const resj_o( pose.residue_type( i.stop()  ).atom_index( "O" ), i.stop()  );

		constraints.push_back( ConstraintCOP( core::scoring::constraints::ConstraintOP(
			new AngleConstraint(resi_n, resi_c, resj_n, create_bb_angle_func( angle1 ))
			)));
		constraints.push_back( ConstraintCOP( core::scoring::constraints::ConstraintOP(
			new AngleConstraint(resj_n, resj_c, resi_n, create_bb_angle_func( angle2 ))
			)));
		constraints.push_back( ConstraintCOP( core::scoring::constraints::ConstraintOP(
			new AngleConstraint(resi_n, resi_c, resj_c, create_bb_angle_func( angle3 ))
			)));
		constraints.push_back( ConstraintCOP( core::scoring::constraints::ConstraintOP(
			new AngleConstraint(resj_n, resj_c, resi_c, create_bb_angle_func( angle4 ))
			)));
	}

	return constraints;
}

core::scoring::func::FuncOP
AutomaticSheetConstraintGenerator::create_bb_angle_func( core::Real const ideal_angle ) const
{
	using core::scoring::constraints::BoundFunc;
	using namespace core::scoring::func;

	return weighted_func( FuncOP( new BoundFunc(
		ideal_angle-angle_tolerance_,
		ideal_angle+angle_tolerance_, sqrt(1.0/42.0), "angle_bb") ) );
}

core::scoring::func::FuncOP
AutomaticSheetConstraintGenerator::create_bb_dihedral_func( core::Real const ideal_dihedral ) const
{
	using core::scoring::constraints::OffsetPeriodicBoundFunc;
	using namespace core::scoring::func;
	core::Real const periodicity = numeric::constants::f::pi;
	return weighted_func( FuncOP( new OffsetPeriodicBoundFunc(
		ideal_dihedral-bb_dihedral_tolerance_,
		ideal_dihedral+bb_dihedral_tolerance_,
		std::sqrt(1.0/42.0), "dihed_bb", periodicity, 0.0 ) ) );

}

core::scoring::func::FuncOP
AutomaticSheetConstraintGenerator::weighted_func( core::scoring::func::FuncOP func ) const
{
	return protocols::constraint_generator::scalar_weighted( func, weight_ );
}

void
AutomaticSheetConstraintGenerator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	using namespace protocols::constraint_generator;
	AttributeList attlist, innerlist, outerlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "sd", xsct_real, "standar deviation", std::to_string(default_sd()) )
		+ XMLSchemaAttribute::attribute_w_default( "distance", xsct_real, "Max distance to consider two residues paired",std::to_string(default_distance()) );

	ConstraintGeneratorFactory::xsd_constraint_generator_type_definition_w_attributes(
		xsd,
		class_name(),
		"Generates constraints between sheet residues putatively paired. Runs DSSP, so it might change a previously assignated secondary structure.",
		attlist );
}

void
AutomaticSheetConstraintGeneratorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	AutomaticSheetConstraintGenerator::provide_xml_schema( xsd );
}

}
} //protocols
} //fold_from_loops
