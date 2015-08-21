// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/public/scenarios/chemically_conjugated_docking/Gp_quantification_metrics.hh
/// @brief contains helper quantification metrics for the original publication of the UBQ_Gp code
/// @author Steven Lewis

#ifndef INCLUDED_apps_public_scenarios_chemically_conjugated_docking_Gp_quantification_metrics_HH
#define INCLUDED_apps_public_scenarios_chemically_conjugated_docking_Gp_quantification_metrics_HH

// Unit Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/chemical/ResidueType.hh>

#include <core/id/AtomID.hh>

#include <protocols/rigid/RB_geometry.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>

#include <utility/FixedSizeLexicographicalIterator.hh>
#include <utility/FixedSizeLexicographicalIterator.tmpl.hh>
#include <utility/fixedsizearray1.hh>

#include <sstream> //ostringstream
//#include <string>

namespace apps {
namespace public1 { //public is a reserved keyword
namespace scenarios {
namespace chemically_conjugated_docking {


/// @details do ubiquitin-ras pair distance measurement and reporting
void ubq_ras_distance(core::pose::Pose const & pose, basic::Tracer & /*TR*/, protocols::jd2::JobOP job_me){

	//Generate distance metrics for desired residue pairs

	//ras residues: Y32 (yes) and W35 (no) - same as resid

	//UBQ residues / resid:
	//L8   174
	//L43  209
	//I44  210
	//V70  236

	//hard-code (oops) which residues to pair distances
	core::Size const ras_array[] = {32, 35};
	core::Size const ras_residues_num(sizeof(ras_array)/sizeof(ras_array[0]));
	utility::vector1<core::Size> const ras_residues(ras_array, ras_array+ras_residues_num);

	core::Size const ubq_array[] = {174, 209, 210, 236};
	core::Size const ubq_residues_num(sizeof(ubq_array)/sizeof(ubq_array[0]));
	utility::vector1<core::Size> const ubq_residues(ubq_array, ubq_array+ubq_residues_num);

	//set up iterator to enumerate pairs
	utility::fixedsizearray1< core::Size, 2 > num_to_be_paired(ras_residues_num);
	num_to_be_paired[2] = ubq_residues_num;
	utility::FixedSizeLexicographicalIterator< 2 > lex( num_to_be_paired );

	std::ostringstream score_name;
	//apparently standard lexicographical iterator for loop - stolen from ZnHash
	for ( lex.begin(); !lex.at_end(); ++lex ) {
		core::Size const ras(ras_residues[lex[1]]);
		core::Size const ubq(ubq_residues[lex[2]]);

		//TR << "residue pair " << lex[1] << " " << ras << " " << lex[2] << " " << ubq << " " << std::endl;

		//calculate distance
		core::id::AtomID const ras_ca( pose.residue_type(ras).atom_index("CA"), ras);
		core::id::AtomID const ubq_ca( pose.residue_type(ubq).atom_index("CA"), ubq);

		core::Real const distance(pose.xyz(ras_ca).distance(pose.xyz(ubq_ca)));

		score_name << "distance_CA_ras" << ras << "_ubq" << ubq;
		job_me->add_string_real_pair(score_name.str(), distance);
		//TR << score_name.str() << " " << distance << std::endl;

		//clear stringstream
		score_name.str("");
		score_name.clear();
	}

	return;
}


void ubq_ras_rotation_angle(
	core::pose::Pose const & pose,
	basic::Tracer & TR,
	protocols::jd2::JobOP job_me,
	core::Size const GTPase_target,
	bool const ubiquitin){

	//Generate "rotation angle" metric

	//determine center of mass of ras & ubiquitin not kept in scope - dangerous
	//because if UBQ backbone moves, center might move; that will cause problems -
	//but this is how it was originally calculated
	// UBQ:ras settings
	// core::Size const ras_begin(1);
	// core::Size const ras_end(166);
	// core::Size const ubq_begin(167);
	// core::Size const ubq_end(242);

	//PDZ:ras settings
	// core::Size const ras_begin(1);
	// core::Size const ras_end(166);
	// core::Size const ubq_begin(167);
	// core::Size const ubq_end(260);//262 or 263, depending on length

	// core::Size const ubq_calculated_center(protocols::geometry::core::pose::residue_center_of_mass(pose, ubq_begin, ubq_end));
	// core::Size const ras_calculated_center(core::pose::residue_center_of_mass(pose, ras_begin, ras_end));

	// TR << "Calculated centers of mass: ubq/pdz: " << ubq_calculated_center << " ras: " << ras_calculated_center << std::endl;

	core::Size const ubq_center(ubiquitin ? 209 : 224); //expected values for ubiquitin and PDZ, respectively
	core::Size const ras_center(79); //expected value
	core::Size const ras_linker(GTPase_target);
	core::Size const ras_reference(112);
	//This is close in space to the ras center and should suffice for being a fixed point-of-reference for dihedral calculation

	//report
	TR << "ras " << ras_center << " ubq " << ubq_center << " ras_reference " << ras_linker << " ras_reference " << ras_reference << std::endl;
	job_me->add_string_real_pair("ras_center", ras_center);
	job_me->add_string_real_pair("ubq_center", ubq_center);
	job_me->add_string_real_pair("ras_linker", ras_linker);
	job_me->add_string_real_pair("ras_reference", ras_reference);

	//get positions
	core::id::AtomID const ras_reference_ca_ID( pose.residue_type(ras_reference).atom_index("CA"), ras_reference );
	core::id::AtomID const ras_center_ca_ID( pose.residue_type(ras_center).atom_index("CA"), ras_center );
	core::id::AtomID const ras_linker_ca_ID( pose.residue_type(ras_linker).atom_index("CA"), ras_linker );
	core::id::AtomID const ubq_center_ca_ID( pose.residue_type(ubq_center).atom_index("CA"), ubq_center );

	numeric::xyzVector<core::Real> const & ras_reference_xyz(pose.xyz(ras_reference_ca_ID));
	numeric::xyzVector<core::Real> const & ras_center_xyz(pose.xyz(ras_center_ca_ID));
	numeric::xyzVector<core::Real> const & ras_linker_xyz(pose.xyz(ras_linker_ca_ID));
	numeric::xyzVector<core::Real> const & ubq_center_xyz(pose.xyz(ubq_center_ca_ID));

	//calculate ras_linker-ubq_center distance
	core::Real const the_distance(ubq_center_xyz.distance(ras_linker_xyz));
	//TR << the_distance << std::endl;
	job_me->add_string_real_pair("ubq-center_ras-linker_distance", the_distance);

	//angle_of gives the angle about B for the triangle ABC (with arguments in that order)
	core::Real const the_angle(numeric::conversions::degrees(angle_of(ubq_center_xyz, ras_linker_xyz, ras_center_xyz)));
	//TR << the_angle << std::endl;
	job_me->add_string_real_pair("ubq-center_ras-linker_ras-center_angle", the_angle);

	//calculate a rotation angle for ubiquitin
	core::Real const the_dihedral(numeric::dihedral_degrees(ubq_center_xyz, ras_linker_xyz, ras_center_xyz, ras_reference_xyz));
	//TR << the_dihedral << std::endl;
	job_me->add_string_real_pair("ubq-center_ras-linker_ras-center_ras-reference_dihedral", the_dihedral);
	return;
}

/// @details This function is specific to the original system for which this code
///was written - if you are not trying to duplicate the initial results you
///should remove it!
void create_extra_output(
	core::pose::Pose const & pose,
	basic::Tracer & TR,
	bool const ubiquitin /*do pair-matching or not*/,
	core::Size const GTPase_target){

	//get a handle to the Job object for reporting
	using protocols::jd2::JobDistributor;
	protocols::jd2::JobOP job_me( JobDistributor::get_instance()->current_job() );

	if ( ubiquitin ) ubq_ras_distance(pose, TR, job_me);

	ubq_ras_rotation_angle(pose, TR, job_me, GTPase_target, ubiquitin);

	return;
}

}//chemically_conjugated_docking
}//scenarios
}//public1
}//apps

#endif //INCLUDED_apps_public_scenarios_chemically_conjugated_docking_Gp_quantification_metrics_HH
