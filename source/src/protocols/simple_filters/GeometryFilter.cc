// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/GeometryFilter.cc
/// @brief
/// @author Possu Huang (possu@u.washington.edu)
/// @author Lei Shi (shilei@u.washington.edu)
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
#include <core/select/residue_selector/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <map>
#include <cmath>
#include <numeric/random/random.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/GeometryFilter.hh>
#include <protocols/simple_filters/GeometryFilterCreator.hh>
#include <protocols/simple_moves/ddG.hh>
//#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.GeometryFilter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP GeometryFilterCreator::create_filter() const { return protocols::filters::FilterOP( new GeometryFilter ); }

// XRW TEMP std::string
// XRW TEMP GeometryFilterCreator::keyname() const { return "Geometry"; }

GeometryFilter::GeometryFilter() :
	filters::Filter( "GeometryFilter" ),
	omega_cutoff_( 165.0 ),
	cart_bonded_cutoff_( 20.0 ),
	filename_( "none" ),
	cst_cutoff_( 10000.0 ),
	start_( 1 ),
	end_( 100000 ),
	count_bad_residues_( false ),
	selector_()
{}

GeometryFilter::~GeometryFilter() = default;

void
GeometryFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	omega_cutoff_ = tag->getOption<core::Real>( "omega", omega_cutoff_ );
	cart_bonded_cutoff_ = tag->getOption<core::Real>( "cart_bonded", cart_bonded_cutoff_ );
	filename_ = tag->getOption< std::string >( "cstfile", filename_ );
	cst_cutoff_ = tag->getOption< core::Real >( "cst_cutoff", cst_cutoff_ );
	start_ = tag->getOption< core::Size>( "start", start_ );
	end_ = tag->getOption< core::Size >( "end", end_ );
	count_bad_residues_ = tag->getOption< bool >( "count_bad_residues", count_bad_residues_ );
	selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data );
}

bool
GeometryFilter::apply( core::pose::Pose const & pose ) const {
	core::Size const bad_residues = compute( pose );
	bool status = false;
	if ( bad_residues == 0 ) {
		TR << "passing." << std::endl;
		status = true;
	} else TR << "failing." << std::endl;
	return status;
}

void
GeometryFilter::report( std::ostream & /*out*/, core::pose::Pose const & pose ) const {
	/*core::Real const dist =*/ compute( pose );
}

core::Real
GeometryFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Size const bad_residues = compute( pose );
	if ( count_bad_residues_ ) {
		return core::Real( bad_residues );
	} else {
		if ( bad_residues == 0 ) {
			return core::Real( 1.0 );
		} else {
			return core::Real( 0.0 );
		}
	}
}


// take per-residue cart_bonded score for evaluating out-liers
core::Size
GeometryFilter::compute( core::pose::Pose const & pose ) const {
	using namespace core::scoring;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using namespace protocols::simple_filters;
	using namespace protocols::simple_moves;
	using core::scoring::ScoreType;

	core::pose::Pose copy_pose;
	if ( is_symmetric( pose ) ) {
		extract_asymmetric_unit( pose, copy_pose, false);
	} else {
		copy_pose = pose;
	}

	copy_pose.update_residue_neighbors();

	core::Size bad_residues = 0;

	// scoring is necessary for Interface to work reliably
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
	scorefxn->set_weight( cart_bonded, 1.0);
	(*scorefxn)(copy_pose);

	// find subset of residues to scan
	core::select::residue_selector::ResidueSubset scan_me;
	if ( selector_ ) {
		scan_me = selector_->apply( copy_pose );
	} else {
		scan_me = core::select::residue_selector::ResidueSubset( copy_pose.size(), true );
	}

	core::Size const start = std::max( start_, Size( 1 ) );
	core::Size const stop = std::min( copy_pose.size(), end_ );
	TR << "Scan residues between " << start << " and " << stop << std::endl;
	for ( Size resnum = start; resnum < stop; resnum++ ) {
		if ( !scan_me[ resnum ] ) continue;

		if ( copy_pose.fold_tree().is_cutpoint( resnum+1 ) || copy_pose.fold_tree().is_jump_point( resnum+1 ) ) continue;
		// TL: skip checking omega values at the end of chains
		//     they are zero and meaningless, but will cause this filter to fail
		if ( copy_pose.chain( resnum ) != copy_pose.chain( resnum+1 ) ) continue;

		// TL: skip if the current residue is not a protein residue, prevents filtering on bogus numbers
		if ( !copy_pose.residue( resnum ).is_protein() ) continue;

		core::Real const weight( (*scorefxn)[ core::scoring::ScoreType( cart_bonded ) ] );
		core::Real const score( copy_pose.energies().residue_total_energies( resnum )[ core::scoring::ScoreType( cart_bonded ) ]);
		core::Real weighted_score = weight * score ;

		// also gets omega values:
		core::Real omega = copy_pose.omega(resnum);

		TR << "residue " << resnum << " name " << copy_pose.residue( resnum ).name3() << " cart_bonded term: " << weighted_score ;
		TR << " omega angle: " << omega << std::endl;

		if ( (std::abs(omega) > 180-omega_cutoff_ && std::abs(omega) < omega_cutoff_ ) && copy_pose.residue( resnum+1 ).name3() == "PRO" )  {
			TR << "omega " << resnum <<" "<< copy_pose.residue( resnum+1 ).name3() << " fail " << std::endl;
			++bad_residues;
		}

		if ( std::abs(omega) < omega_cutoff_ && copy_pose.residue( resnum+1 ).name3() != "PRO" )  {
			TR << "omega " << resnum <<" "<< copy_pose.residue( resnum+1 ).name3() << " fail " << std::endl;
			++bad_residues;
		}

		if ( weighted_score >= cart_bonded_cutoff_ ) {
			TR << "cart_bond " << resnum <<" "<< copy_pose.residue( resnum+1 ).name3() << " fail " << std::endl;
			++bad_residues;
		}
	}

	//check the last residues
	Size resnum = std::min( copy_pose.size(), end_ );
	core::Real const weight( (*scorefxn)[ core::scoring::ScoreType( cart_bonded ) ] );
	core::Real const score( copy_pose.energies().residue_total_energies( resnum )[ core::scoring::ScoreType( cart_bonded ) ]);
	core::Real weighted_score = weight * score ;
	TR << "residue " << resnum << " name " << copy_pose.residue( resnum ).name3() << " cart_bonded term: " << weighted_score ;
	TR << weighted_score  << std::endl;

	if ( filename_ != "none" ) {
		TR << "Evaluate constraint energy?" << std::endl;
		ConstraintSetMoverOP cst_set_mover( new ConstraintSetMover() );
		cst_set_mover->constraint_file( filename_ );
		cst_set_mover->apply( copy_pose );

		// only evaluate atom_pair constraints for now
		scorefxn->set_weight( core::scoring::atom_pair_constraint, 1.0 );
		(*scorefxn)( copy_pose );
		ScoreTypeFilter const constraint_filter( scorefxn , atom_pair_constraint, cst_cutoff_ );
		bool CScore(constraint_filter.apply( copy_pose ));
		if ( !CScore ) {
			++bad_residues;
		}
	}

	return bad_residues;
}

std::string GeometryFilter::name() const {
	return class_name();
}

std::string GeometryFilter::class_name() {
	return "Geometry";
}

void GeometryFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("omega", xsct_real, "cutoff for omega angle of peptide plane, Cis-proline is also considered. Works for multiple chains", "165.0")
		+ XMLSchemaAttribute::attribute_w_default("cart_bonded", xsct_real, "bond angle and length penalty score", "20.0")
		+ XMLSchemaAttribute::attribute_w_default("cstfile", xs_string, "if specified, the given constraint file will be used to introduce constraints into the pose. Only atom pair constraints will be used. The constraint scores will be checked to see if they are lower than cst_cutoff. If the constraint score is greater than cst_cutoff, the pose will fail.", "filename")
		+ XMLSchemaAttribute::attribute_w_default("cst_cutoff", xsct_real, "cutoff for use with cstfile option", "10000.0")
		+ XMLSchemaAttribute::attribute_w_default("start", xsct_non_negative_integer, "starting residue number to scan", "1")
		+ XMLSchemaAttribute::attribute_w_default("end", xsct_non_negative_integer, "ending residue number", "100000")
		+ XMLSchemaAttribute::attribute_w_default("count_bad_residues", xsct_rosetta_bool, "If true, the number of residues failing the filter will be computed and returned as the filter score by report_sm(). If false, the filter score will be either 1.0 (i.e. all residues pass) or 0.0 (i.e. a residue failed the filter). Default: false", "false");

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "If specified, only residues selected by the user-specified residue selector will be scanned. By default, all residues are selected. If start and/or end are also set, only residues selected by the residue_selector AND within the range [start, end] will be scanned.");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "Filters poses based upon their geometry (omega angle, for example(", attlist );
}

std::string GeometryFilterCreator::keyname() const {
	return GeometryFilter::class_name();
}

protocols::filters::FilterOP
GeometryFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new GeometryFilter );
}

void GeometryFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GeometryFilter::provide_xml_schema( xsd );
}



}
}
