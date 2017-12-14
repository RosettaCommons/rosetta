// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @brief
/// @author Neil King ( neilking@uw.edu )
/// @author Javier Castellanos ( javiercv@uw.edu )

// Unit headers
#include <devel/matdes/Symmetrizer.hh>
#include <devel/matdes/SymmetrizerMoverCreator.hh>

// Package headers
#include <devel/matdes/SymmetrizerSampler.hh>

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "devel.matdes.Symmetrizer" );


namespace devel {
namespace matdes {

using namespace core;
using namespace utility;

// -------------  Mover Creator -------------
// XRW TEMP std::string
// XRW TEMP SymmetrizerMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return Symmetrizer::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SymmetrizerMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new Symmetrizer );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP Symmetrizer::mover_name()
// XRW TEMP {
// XRW TEMP  return "Symmetrizer";
// XRW TEMP }
// -------------  Mover Creator -------------

Symmetrizer::Symmetrizer() :
	symm_file_(""),
	radial_disp_(0.0),
	radial_disp_min_(0.0),
	radial_disp_max_(0.0),
	angle_(0.0),
	angle_min_(0.0),
	angle_max_(0.0),
	radial_disp_delta_(0.0),
	angle_delta_(0.0),
	symmetry_axis_('z'),
	explore_grid_(false),
	sampling_mode_("single_dock")
{ }

protocols::moves::MoverOP
Symmetrizer::clone() const {
	return protocols::moves::MoverOP( new Symmetrizer( *this ) );
}

protocols::moves::MoverOP
Symmetrizer::fresh_instance() const {
	return protocols::moves::MoverOP( new Symmetrizer() );
}


Real
Symmetrizer::get_angle() {
	if ( explore_grid_ || sampling_mode_ == "grid" ) {
		return SymmetrizerSampler::get_instance()->get_angle();
	} else if ( sampling_mode_ == "uniform" ) {
		return angle_min_ + ( angle_max_ - angle_min_) * numeric::random::rg().uniform();
	} else if ( sampling_mode_ == "gaussian" ) {
		return angle_ + angle_delta_ * numeric::random::rg().gaussian();
	} else {
		return angle_;
	}
}


Real
Symmetrizer::get_radial_disp() {
	if ( explore_grid_ || sampling_mode_ == "grid" ) {
		return SymmetrizerSampler::get_instance()->get_radial_disp();
	} else if ( sampling_mode_ == "uniform" ) {
		return radial_disp_min_ + ( radial_disp_max_ - radial_disp_min_) * numeric::random::rg().uniform();
	} else if ( sampling_mode_ == "gaussian" ) {
		return radial_disp_ + radial_disp_delta_ * numeric::random::rg().gaussian();
	} else {
		return radial_disp_;
	}
}

void
Symmetrizer::apply(Pose & pose) {
	using core::conformation::symmetry::SymmetryInfoCOP;
	using core::conformation::symmetry::SymDof;
	using core::pose::Pose;
	using Vec = numeric::xyzVector<Real>;
	using Mat = numeric::xyzMatrix<Real>;
	TR << "mode: " << sampling_mode_ << std::endl;
	core::pose::symmetry::make_symmetric_pose(pose, symm_file_);
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	std::map<Size,SymDof> dofs = sym_info->get_dofs();
	int sym_jump = 0;
	for ( auto & dof : dofs ) {
		Size jump_num = dof.first;
		if ( sym_jump == 0 ) {
			sym_jump = jump_num;
		} else {
			utility_exit_with_message("Can only handle one subunit!");
		}
	}
	if ( sym_jump == 0 ) {
		utility_exit_with_message("No jump defined!");
	}

	core::kinematics::Jump j = pose.jump(sym_jump);

	const Vec init_trans = pose.jump(sym_jump).get_translation();
	const Mat init_rot = pose.jump(sym_jump).get_rotation();

	Real radial_disp = get_radial_disp();
	Real angle = get_angle();
	TR << "radial_disp = " << radial_disp << " angle = " << angle << std::endl;
	core::pose::setPoseExtraScore(pose, "radial_disp", radial_disp);
	core::pose::setPoseExtraScore(pose, "angle", angle);
	Vec translation;
	Mat rotation;
	switch(symmetry_axis_) {
	case 'x' :
		translation = Vec(get_radial_disp(),0,0) + init_trans;
		rotation = Mat(numeric::x_rotation_matrix_degrees( angle ) * init_rot);
		break;

	case 'y' :
		translation = Vec(0, get_radial_disp(), 0) + init_trans;
		rotation = Mat(numeric::y_rotation_matrix_degrees( angle )* init_rot);
		break;

	case 'z' :
		translation = Vec(0,0, get_radial_disp()) + init_trans;
		rotation = Mat(numeric::z_rotation_matrix_degrees( angle ) * init_rot);
		break;

	default :
		utility_exit_with_message(std::string(1, symmetry_axis_) +
			" is not a valid axis (x,y or z). Use lower case");
	}
	j.set_translation( translation );
	j.set_rotation( rotation );
	pose.set_jump(sym_jump,j);

	if ( explore_grid_ ) {
		SymmetrizerSampler::get_instance()->step();
	}
}

void
Symmetrizer::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const & ) {

	// Turn symmetry hacks on
	if ( !basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].user() ) {
		basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].value( "dummy" );
	}

	using std::string;
	symm_file_ = tag->getOption<string>( "symm_file" );
	symmetry_axis_ = tag->getOption<char>("axis", 'z');

	angle_min_ = tag->getOption<Real>("angle_min", 0.0);
	angle_max_ = tag->getOption<Real>("angle_max", 0.0);
	radial_disp_min_ = tag->getOption<Real>("radial_disp_min", 0.0);
	radial_disp_max_ = tag->getOption<Real>("radial_disp_max", 0.0);
	sampling_mode_= tag->getOption<std::string>("sampling_mode", "uniform");
	TR << "Setting sampling mode to " << sampling_mode_ << std::endl;
	explore_grid_ = tag->getOption<bool>("grid", false);
	if ( sampling_mode_ == "grid" || explore_grid_ ) {
		auto angle_step = tag->getOption<Real>("angle_step");
		auto radial_disp_step = tag->getOption<Real>("radial_disp_step");
		TR << "Setting the exploration grid." << std::endl;
		SymmetrizerSampler::get_instance()->set_angle_range(angle_min_, angle_max_, angle_step);
		SymmetrizerSampler::get_instance()->set_radial_disp_range(radial_disp_min_, radial_disp_max_, radial_disp_step);
	} else if ( sampling_mode_ == "gaussian" ) {
		angle_delta_ = tag->getOption<Real>("angle_delta", 0.0);
		radial_disp_delta_ = tag->getOption<Real>("radial_disp_delta", 0.0);
	}
	radial_disp_ = tag->getOption<Real>( "radial_disp" ,0.0);
	angle_ = tag->getOption<Real>( "angle",0.0 );

}

std::string Symmetrizer::get_name() const {
	return mover_name();
}

std::string Symmetrizer::mover_name() {
	return "Symmetrizer";
}

void Symmetrizer::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	// Ultimately this should go in utility because it's used in multiple places
	// To avoid a laborious recompile, we will instead 'namespace' it
	XMLSchemaRestriction sampling_mode_enumeration;
	sampling_mode_enumeration.name( "Symmetrizer_sampling_mode_choices" );
	sampling_mode_enumeration.base_type( xs_string );
	sampling_mode_enumeration.add_restriction( xsr_enumeration, "grid" );
	sampling_mode_enumeration.add_restriction( xsr_enumeration, "uniform" );
	sampling_mode_enumeration.add_restriction( xsr_enumeration, "gaussian" );
	xsd.add_top_level_element( sampling_mode_enumeration );

	XMLSchemaRestriction axis_pattern;
	axis_pattern.name( "axis_char" );
	axis_pattern.base_type( xsct_char );
	axis_pattern.add_restriction( xsr_pattern,  "x|y|z" );

	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("symm_file", xs_string, "Symmetry definition file to create fully symmetric pose from input asymmetric unit.")
		+ XMLSchemaAttribute::attribute_w_default("axis", "axis_char" , "Axis around which to symmetrize.", "z")
		+ XMLSchemaAttribute::attribute_w_default("angle_min", xsct_real, "For use with uniform sampling. ", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("angle_max", xsct_real, "For use with uniform sampling.", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("sampling_mode", "Symmetrizer_sampling_mode_choices", "The default is 'uniform' , but this can also be 'grid' or 'gaussian'", "uniform")
		+ XMLSchemaAttribute::attribute_w_default("grid", xsct_rosetta_bool, "XRW TODO (what does it DO?) The default is false; but if it is set to true, you must also set angle_step and radial_disp_step.", "false")
		+ XMLSchemaAttribute("angle_step", xsct_real, "xRW TODO (what does it DO?) You must use this option when sampling_mode=grid OR when grid=true")
		+ XMLSchemaAttribute("radial_disp_step", xsct_real, "You must use this option when sampling_mode=grid OR when grid=1")
		+ XMLSchemaAttribute::attribute_w_default("angle_delta", xsct_real, "The default is '0.0', but when using sampling_mode=gaussian, but you change this value.", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("radial_disp_min", xsct_real, "XSD XRW TO DO", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("radial_disp_max", xsct_real, "XSD XRW TO DO", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("radial_disp_delta", xsct_real, "The default is '0.0', but when using sampling_mode=gaussian, but you change this value.", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("radial_disp", xsct_real, "Translation of asymmetric unit if necessary", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("angle", xsct_real, "Rotation angle of asymmetric unit if necessary.", "0.0");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Functionally an optimization mover; will take a pose with sufficiently small deviations from symmetry and resolve them.", attlist );
}

std::string SymmetrizerMoverCreator::keyname() const {
	return Symmetrizer::mover_name();
}

protocols::moves::MoverOP
SymmetrizerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new Symmetrizer );
}

void SymmetrizerMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	Symmetrizer::provide_xml_schema( xsd );
}


} // matdes
} // devel
