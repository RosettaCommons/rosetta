
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
#include <devel/matdes/SymmetricMultimerDesign.hh>
#include <devel/matdes/SymmetricMultimerDesignMoverCreator.hh>

// Package headers

// project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/jd2/parser/BluePrint.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

namespace devel {
namespace matdes {

using namespace core;
using namespace utility;

// -------------  Mover Creator -------------
std::string
SymmetrizeCreator::keyname() const
{
	return SymmetrizeCreator::mover_name();
}

MoverOP
SymmetrizeCreator::create_mover() const {
	return new ConstrainedDesign;
}

std::string
SymmetrizeCreator::mover_name()
{
	return "Symmetrize";
}
// -------------  Mover Creator -------------

void
Symmetrize::apply(Pose & pose) {
	core::pose::symmetry::make_symmetric_pose(pose, sym_file_);
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	std::map<Size,SymDof> dofs = sym_info->get_dofs();
 	int sym_jump = 0;
 	for(std::map<Size,SymDof>::iterator i = dofs.begin(); i != dofs.end(); i++) {
   	Size jump_num = i->first;
   	if (sym_jump == 0) {
	 		sym_jump = jump_num;
   	} else {
   		utility_exit_with_message("Can only handle one subunit!");
   	}
 	}
 	if (sym_jump == 0)
   	utility_exit_with_message("No jump defined!");

	core::kinematics::Jump j;
	//fpd if we have multiple symm DOFs or we have radial_disp==angle==0 then use input conformation
	j = pose.jump(sym_jump);

	Mat init_rot = pose.jump(sym_jump).get_rotation();

	
	Vec translation;
	Mat rotation;
	switch(symmetry_axis) {
				case 'x' : 
					translation = Vec(get_radial_disp(),0,0);
					rotation = Mat(numeric::x_rotation_matrix_degrees(get_angle()* init_rot));
				break;

				case 'y' : 
					translation = Vec(0, get_radial_disp(), 0);
					rotation = Mat(numeric::y_rotation_matrix_degrees(get_angle() * init_rot));
				break;

				case 'z' : 
					translation = Vec(0,0, get_radial_disp());
					rotation = Mat(numeric::z_rotation_matrix_degrees(get_angle()* init_rot));
				break;
	}
	j.set_translation( translation );
	j.set_rotation( rotation );
	pose.set_jump(sym_jump,j); 
}

}
}
