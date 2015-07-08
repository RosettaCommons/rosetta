// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
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

static thread_local basic::Tracer TR( "devel.matdes.Symmetrizer" );


namespace devel {
namespace matdes {

using namespace core;
using namespace utility;

// -------------  Mover Creator -------------
std::string
SymmetrizerMoverCreator::keyname() const
{
	return SymmetrizerMoverCreator::mover_name();
}

protocols::moves::MoverOP
SymmetrizerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new Symmetrizer );
}

std::string
SymmetrizerMoverCreator::mover_name()
{
	return "Symmetrizer";
}
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
	if(explore_grid_ || sampling_mode_ == "grid" )
		return SymmetrizerSampler::get_instance().get_angle();
	else if(sampling_mode_ == "uniform")
		return angle_min_ + ( angle_max_ - angle_min_) * numeric::random::rg().uniform();
	else if(sampling_mode_ == "gaussian")
		return angle_ + angle_delta_ * numeric::random::rg().gaussian();
	else
		return angle_;
}


Real
Symmetrizer::get_radial_disp() { 
	if(explore_grid_ || sampling_mode_ == "grid" )
		return SymmetrizerSampler::get_instance().get_radial_disp();
	else if(sampling_mode_ == "uniform") 
		return radial_disp_min_ + ( radial_disp_max_ - radial_disp_min_) * numeric::random::rg().uniform();
	else if(sampling_mode_ == "gaussian")
		return radial_disp_ + radial_disp_delta_ * numeric::random::rg().gaussian();
	else
		return radial_disp_;
}

void
Symmetrizer::apply(Pose & pose) {
	using core::conformation::symmetry::SymmetryInfoCOP;
	using core::conformation::symmetry::SymDof;
	using core::pose::Pose;
	typedef numeric::xyzVector<Real> Vec;
	typedef numeric::xyzMatrix<Real> Mat;
	TR << "mode: " << sampling_mode_ << std::endl;
	core::pose::symmetry::make_symmetric_pose(pose, symm_file_);
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	std::map<Size,SymDof> dofs = sym_info->get_dofs();
 	int sym_jump = 0;
 	for(std::map<Size,SymDof>::iterator i = dofs.begin(), end = dofs.end(); i != end; ++i ) {
   	Size jump_num = i->first;
   	if (sym_jump == 0) {
	 		sym_jump = jump_num;
   	} else {
   		utility_exit_with_message("Can only handle one subunit!");
   	}
 	}
 	if (sym_jump == 0)
   	utility_exit_with_message("No jump defined!");

	core::kinematics::Jump j = pose.jump(sym_jump);

	const	Vec init_trans = pose.jump(sym_jump).get_translation();
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

				default:
					utility_exit_with_message(std::string(1, symmetry_axis_) +
							" is not a valid axis (x,y or z). Use lower case");
	}
	j.set_translation( translation );
	j.set_rotation( rotation );
	pose.set_jump(sym_jump,j); 
	
	if(explore_grid_) {
		SymmetrizerSampler::get_instance().step();
	}
}

void 
Symmetrizer::parse_my_tag( TagCOP const tag,
										 basic::datacache::DataMap &,
										 Filters_map const &,
										 Movers_map const &,
										 Pose const & ) {

	// Turn symmetry hacks on
	if( !basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].user() )
		basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].value( "dummy" );

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
	if( sampling_mode_ == "grid" || explore_grid_) {
		Real angle_step = tag->getOption<Real>("angle_step");
		Real radial_disp_step = tag->getOption<Real>("radial_disp_step");
		TR << "Setting the exploration grid." << std::endl;
		SymmetrizerSampler::get_instance().set_angle_range(angle_min_, angle_max_, angle_step);
		SymmetrizerSampler::get_instance().set_radial_disp_range(radial_disp_min_, radial_disp_max_, radial_disp_step);
	} else if (sampling_mode_ == "gaussian") {
		angle_delta_ = tag->getOption<Real>("angle_delta", 0.0);
	radial_disp_delta_ = tag->getOption<Real>("radial_disp_delta", 0.0);
	}
	radial_disp_ = tag->getOption<Real>( "radial_disp" ,0.0);
	angle_ = tag->getOption<Real>( "angle",0.0 );

}
void Symmetrizer::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & ,
				utility::lua::LuaObject const & ,
				protocols::moves::MoverCacheSP ) {

	// Turn symmetry hacks on
	basic::options::option[basic::options::OptionKeys::symmetry::symmetry_definition].value( "dummy" );

	symm_file_ = def["symm_file"].to<std::string>();
	symmetry_axis_ = def["axis"] ? def["axis"].to<std::string>()[0] : 'z';

	explore_grid_ = def["grid"] ? def["grid"].to<bool>() : false;
	if(explore_grid_) {
		Real angle_min = def["angle_min"].to<Real>();
		Real angle_max = def["angle_max"].to<Real>();
		Real angle_step = def["angle_step"].to<Real>();
		Real radial_disp_min = def["radial_disp_min"].to<Real>();
		Real radial_disp_max = def["radial_disp_max"].to<Real>();
		Real radial_disp_step = def["radial_disp_step"].to<Real>();
		TR << "Setting the exploration grid." << std::endl;
		SymmetrizerSampler::get_instance().set_angle_range(angle_min, angle_max, angle_step);
		SymmetrizerSampler::get_instance().set_radial_disp_range(radial_disp_min, radial_disp_max, radial_disp_step);
	} else {
		radial_disp_ = def["radial_disp"] ? def["radial_disp"].to<Real>() : 0.0;
		angle_ = def["angle"] ? def["angle"].to<Real>() : 0.0;
	}
}

} // matdes
} // devel
