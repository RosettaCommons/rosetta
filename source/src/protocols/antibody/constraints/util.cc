// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/util.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/constraints/util.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>

#include <protocols/antibody/AntibodyInfo.hh>

#include <numeric/conversions.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>

#include <iostream>
#include <fstream>
#include <cctype>

#include <boost/algorithm/string.hpp>

static basic::Tracer TR("antibody.constraints.util");

namespace protocols {
namespace antibody {
namespace constraints {
	using namespace protocols::antibody;
	using namespace protocols::antibody::clusters;
	using utility::vector1;

// AMW: do not pass constraint_type by reference, it'll kill a unit test
bool
cdr_has_res_constraints(AntibodyInfoCOP ab_info, core::pose::Pose & pose, CDRNameEnum cdr, std::string const constraint_type){
	using namespace core::scoring::constraints;

	core::Size start_res = ab_info->get_CDR_start(cdr, pose);
	core::Size end_res = ab_info->get_CDR_end(cdr, pose);

	std::map< core::Size, bool> cst_found;

	//Initialize our map of whether the CDR residue has constraints.
	for (core::Size i = start_res; i <= end_res; ++i){
		cst_found[i] = false;
	}
	utility::vector1< ConstraintCOP > csts = pose.constraint_set()->get_all_constraints();
	for (core::Size i = 1; i <= csts.size(); ++i){
		if (csts[i]->type() != constraint_type){ continue; }

		utility::vector1< core::Size > residues = csts[i]->residues();
		for (core::Size x = 1; x <= residues.size(); ++x){
			cst_found[residues[x]] = true;
		}
	}

	//Check that all residues have constraints of the particular type:
	for (core::Size i = start_res ; i <= end_res; ++i){
		if (! cst_found[i]){ return false; }
	}

	return true;
}

void
add_harmonic_dihedral_cst_general(
	AntibodyInfoCOP ab_info, core::pose::Pose & pose,
	CDRNameEnum const cdr,
	core::Real phi_sd_deg /* 20.0 */, core::Real psi_sd_deg /* 30.0 */)
{
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	if ( !ab_info->has_CDR( cdr ) ){
		return;
	}


	using core::id::AtomID;
	core::Real phi_sd_rad = numeric::conversions::radians(phi_sd_deg);
	core::Real psi_sd_rad = numeric::conversions::radians(psi_sd_deg);

	for (core::Size i = ab_info->get_CDR_start(cdr, pose); i <= ab_info->get_CDR_end(cdr, pose); ++i){

		//Current residue designated as 1
		AtomID C_0 =   AtomID(pose.residue(i-1).atom_index("C"), i-1);
		AtomID N_1 =   AtomID(pose.residue(i).atom_index("N"), i );
		AtomID CA_1 = AtomID(pose.residue(i).atom_index("CA"), i );
		AtomID C_1 =   AtomID(pose.residue(i).atom_index("C"), i );
		AtomID N_2 =   AtomID(pose.residue(i+1).atom_index("N"), i+1 );

		core::Real phi = numeric::conversions::radians(pose.phi(i));
		core::Real psi = numeric::conversions::radians(pose.psi(i));

		CircularHarmonicFuncOP phi_func( new core::scoring::func::CircularHarmonicFunc(phi, phi_sd_rad) );
		CircularHarmonicFuncOP psi_func( new core::scoring::func::CircularHarmonicFunc(psi, psi_sd_rad) );

		DihedralConstraintOP phi_cst( new DihedralConstraint(C_0, N_1, CA_1, C_1, phi_func) );
		DihedralConstraintOP psi_cst( new DihedralConstraint(N_1, CA_1, C_1, N_2, psi_func) );

		pose.add_constraint(phi_cst);
		pose.add_constraint(psi_cst);

	}
}

} //constraints
} //antibody
} //protocols
