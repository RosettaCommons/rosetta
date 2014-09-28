// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

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

#include <iostream>
#include <fstream>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <numeric/conversions.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

static basic::Tracer TR("antibody.constraints.util");

namespace protocols {
namespace antibody {
namespace constraints {
	using namespace protocols::antibody;
	using namespace protocols::antibody::clusters;
	using utility::vector1;
	
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
add_harmonic_cluster_cst_or_coordinate_cst(AntibodyInfoOP ab_info, core::pose::Pose & pose, CDRNameEnum const cdr, core::Real coordinate_cst_sd){
	bool constraint_add_successful = add_harmonic_cluster_constraint(ab_info, pose, ab_info->get_CDR_cluster(cdr)->cluster());
	
	if (! constraint_add_successful ){
		core::Size start_res = ab_info->get_CDR_start(cdr, pose);
		core::Size end_res = ab_info->get_CDR_end(cdr, pose);
			
		TR << "Adding coordinate constraints for " << ab_info->get_CDR_name(cdr) << std::endl;
		core::scoring::constraints::add_coordinate_constraints(pose, start_res, end_res, coordinate_cst_sd, false /* include_sc */);
	}
	
}

void
add_harmonic_cluster_cst_or_coordinate_cst(AntibodyInfoOP ab_info, core::pose::Pose & pose, CDRNameEnum const cdr, CDRClusterEnum const cluster, core::Real coordinate_cst_sd){
	bool constraint_add_successful = add_harmonic_cluster_constraint(ab_info, pose, cluster);
	
	if (! constraint_add_successful ){
		core::Size start_res = ab_info->get_CDR_start(cdr, pose);
		core::Size end_res = ab_info->get_CDR_end(cdr, pose);
			
		TR << "Adding coordinate constraints for " << ab_info->get_CDR_name(cdr) << std::endl;
		core::scoring::constraints::add_coordinate_constraints(pose, start_res, end_res, coordinate_cst_sd, false /* include_sc */);
	}
	
}

void
add_harmonic_cluster_cst_or_dihedral_cst(
	AntibodyInfoOP ab_info, core::pose::Pose & pose, 
	CDRNameEnum const cdr,
	core::Real phi_sd_deg /* 23.0 */, core::Real psi_sd_deg /* 42.0 */)
{
	bool constraint_add_successful = add_harmonic_cluster_constraint(ab_info, pose, ab_info->get_CDR_cluster(cdr)->cluster());
	if (! constraint_add_successful){
		add_harmonic_dihedral_cst_general(ab_info, pose, cdr, phi_sd_deg, psi_sd_deg);
	}
}

void
add_harmonic_cluster_cst_or_dihedral_cst(
	AntibodyInfoOP ab_info, core::pose::Pose & pose, 
	CDRNameEnum const cdr, CDRClusterEnum const cluster,
	core::Real phi_sd_deg /* 23.0 */, core::Real psi_sd_deg /* 42.0 */)
{
	bool constraint_add_successful = add_harmonic_cluster_constraint(ab_info, pose, cluster);
	if (! constraint_add_successful){
		TR << "Adding general dihedral constraints" << std::endl;
		add_harmonic_dihedral_cst_general(ab_info, pose, cdr, phi_sd_deg, psi_sd_deg);
	}
}

void
add_harmonic_dihedral_cst_general(
	AntibodyInfoOP ab_info, core::pose::Pose & pose,
	CDRNameEnum const cdr,
	core::Real phi_sd_deg /* 23.0 */, core::Real psi_sd_deg /* 42.0 */) 
{
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	
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


std::map< CDRNameEnum, bool>
add_harmonic_cluster_constraints(AntibodyInfoOP ab_info, core::pose::Pose & pose){
	std::map< CDRNameEnum, bool> result;
	for (core::Size i = 1; i <= core::Size(ab_info->get_total_num_CDRs()); ++i){
		CDRNameEnum cdr_name = static_cast<CDRNameEnum>(i);
		result[cdr_name] = add_harmonic_cluster_constraint(ab_info, pose, ab_info->get_CDR_cluster(cdr_name)->cluster());
	}
	return result;
}

std::map< CDRNameEnum, bool>
add_harmonic_cluster_constraints(AntibodyInfoOP ab_info, core::pose::Pose & pose, utility::vector1< core::scoring::constraints::ConstraintCOP > constraints){
	
	std::map< CDRNameEnum, bool> result;
	for (core::Size i = 1; i <= core::Size(ab_info->get_total_num_CDRs()); ++i){
		CDRNameEnum cdr_name = static_cast<CDRNameEnum>(i);
		result[cdr_name] = add_harmonic_cluster_constraint(ab_info, pose, ab_info->get_CDR_cluster(cdr_name)->cluster(), constraints);
	}
	return result;
}

bool
add_harmonic_cluster_constraint(AntibodyInfoCOP ab_info, core::pose::Pose & pose, CDRClusterEnum const cluster){

	using namespace core::scoring::constraints;

	
	std::string fname = get_harmonic_cluster_constraint_filename(ab_info, cluster);
	if (fname=="NA"){return false;}
	try {
		ConstraintSetOP cst = ConstraintIO::get_instance()->read_constraints(fname, ConstraintSetOP( new ConstraintSet ), pose);

		pose.add_constraints(cst->get_all_constraints());
		return true;	
	}
	catch(utility::excn::EXCN_Exception &excn){
		TR<< "Problem adding dihedral constraints for CDR cluster." <<std::endl;
		return false;
	}

}

bool
add_harmonic_cluster_constraint(AntibodyInfoCOP ab_info, core::pose::Pose & pose, CDRClusterEnum const cluster, utility::vector1< core::scoring::constraints::ConstraintCOP > constraints){

	using namespace core::scoring::constraints;
	
	if (cluster==NA) return false;
	
	std::string fname = get_harmonic_cluster_constraint_filename(ab_info, cluster);
	if (fname=="NA"){return false;}
	
	try {
		ConstraintSetOP cst = ConstraintIO::get_instance()->read_constraints(fname, ConstraintSetOP( new ConstraintSet ), pose);

		vector1< ConstraintCOP > local_csts = cst->get_all_constraints();
		pose.add_constraints(local_csts);

		constraints.insert(constraints.end(), local_csts.begin(), local_csts.end());
		return true;
	} catch (utility::excn::EXCN_Exception &excn){
		TR<< "Problem adding dihedral constraints for CDR cluster." <<std::endl;
		return false;
	}
	
}


std::string
get_harmonic_cluster_constraint_filename(AntibodyInfoCOP ab_info, CDRClusterEnum const cluster ){
	
	using namespace basic::options;
	
	std::string fname;
	std::string cluster_type = ab_info->get_cluster_name(cluster);
	if (cluster_type=="NA") {
		TR<< "Cannot add cluster dihedral constraint to cdr cluster of type NA.  Skipping."<<std::endl;
		return "NA";
	}
	std::string path = "sampling/antibodies/cluster_based_constraints/CircularHarmonic/";
	std::string extension = ".txt";
	std::string specific_path = path + cluster_type + extension;
	fname = option[ OptionKeys::in::path::database ](1).name() + specific_path;
	if( !utility::file::file_exists(fname)) {
		TR<< "Fname "<<fname<<" Does not exist.  No constraint will be added."<<std::endl;
		return "NA";
	}
	return fname;
}

} //constraints
} //antibody
} //protocols
