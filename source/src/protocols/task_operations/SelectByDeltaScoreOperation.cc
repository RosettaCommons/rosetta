// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Restrict design to residues passing a user-specified threshold on a given score type.
/// @author Jacob Bale (balej@uw.edu)

// Unit Headers
#include <protocols/task_operations/SelectByDeltaScoreOperation.hh>
#include <protocols/task_operations/SelectByDeltaScoreOperationCreator.hh>

// Project Headers
#include <protocols/rosetta_scripts/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>


// C++ Headers

static basic::Tracer TR( "protocols.task_operations.SelectByDeltaScoreOperation" );

namespace protocols {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

core::pack::task::operation::TaskOperationOP
SelectByDeltaScoreOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new SelectByDeltaScoreOperation );
}

void SelectByDeltaScoreOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SelectByDeltaScoreOperation::provide_xml_schema( xsd );
}

std::string SelectByDeltaScoreOperationCreator::keyname() const
{
	return SelectByDeltaScoreOperation::keyname();
}

/// @brief default constructor
SelectByDeltaScoreOperation::SelectByDeltaScoreOperation() :
	TaskOperation(),
	scorefxn_(/* NULL */),
	threshold_( 100 ),
	reference_pose_(/* NULL */)
{}

/// @brief destructor
SelectByDeltaScoreOperation::~SelectByDeltaScoreOperation()= default;

/// @brief make clone
core::pack::task::operation::TaskOperationOP
SelectByDeltaScoreOperation::clone() const{
	return core::pack::task::operation::TaskOperationOP( new SelectByDeltaScoreOperation( *this ) );
}

core::scoring::ScoreFunctionOP SelectByDeltaScoreOperation::scorefxn() const{ return scorefxn_; }
void SelectByDeltaScoreOperation::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){ scorefxn_ = scorefxn; }

core::scoring::ScoreType SelectByDeltaScoreOperation::score_type() const{ return score_type_; }
void SelectByDeltaScoreOperation::score_type( core::scoring::ScoreType const score_type ){ score_type_ = score_type; }

std::string SelectByDeltaScoreOperation::score_type_name() const{ return score_type_name_; }
void SelectByDeltaScoreOperation::score_type_name( std::string const & score_type_name ){ score_type_name_ = score_type_name; }

core::Real SelectByDeltaScoreOperation::threshold() const{ return threshold_; }
void SelectByDeltaScoreOperation::threshold( core::Real threshold ) { threshold_ = threshold; }

bool SelectByDeltaScoreOperation::lower() const{ return lower_; }
void SelectByDeltaScoreOperation::lower( bool lower ){ lower_ = lower; }

core::pose::PoseOP SelectByDeltaScoreOperation::reference_pose() const{ return reference_pose_; }
void SelectByDeltaScoreOperation::reference_pose( core::pose::PoseOP reference_pose ){ reference_pose_ = reference_pose; }

bool SelectByDeltaScoreOperation::individual_hbonds() const{ return individual_hbonds_; }
void SelectByDeltaScoreOperation::individual_hbonds( bool individual_hbonds ){ individual_hbonds_ = individual_hbonds; }

void
SelectByDeltaScoreOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	core::pose::Pose refp( *reference_pose_ );
	//refp.dump_pdb("reference_pose.pdb");
	//Loop over the independent residues, assess score_type for each residue, and restrict residues to repacking based on whether or not they pass the user-provided threshold.
	core::Size nres_asymmetric_unit, ref_nres_asymmetric_unit;
	if ( core::conformation::symmetry::is_symmetric( refp.conformation() ) ) {
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(refp);
		ref_nres_asymmetric_unit = symm_info->num_independent_residues();
	} else {
		ref_nres_asymmetric_unit = refp.size();
	}
	if ( core::conformation::symmetry::is_symmetric( pose.conformation() ) ) {
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
		nres_asymmetric_unit = symm_info->num_independent_residues();
	} else {
		nres_asymmetric_unit = pose.size();
	}

	if ( nres_asymmetric_unit != ref_nres_asymmetric_unit ) {
		utility_exit_with_message( "Reference pose and current pose have a different number of residues" );
	}

	core::Real delta_score;
	core::pose::Pose p = pose;
	std::string pymol_selection("select score_type_selection, resi ");
	core::scoring::hbonds::HBondSet hbset, hbset_ref;
	core::scoring::EnergyMap em, ref_em;
	if ( individual_hbonds_ ) {
		core::scoring::methods::EnergyMethodOptions myopt = scorefxn()->energy_method_options();
		myopt.hbond_options().decompose_bb_hb_into_pair_energies(true);
		scorefxn()->set_energy_method_options(myopt);
		scorefxn()->score( refp );
		scorefxn()->score( p );
		// get the HBondSet
		core::scoring::hbonds::fill_hbond_set(p,false,hbset);
		core::scoring::hbonds::fill_hbond_set(refp,false,hbset_ref);
	} else {
		scorefxn()->score( refp );
		scorefxn()->score( p );
	}

	for ( core::Size i=1; i<=nres_asymmetric_unit; i++ ) {
		if ( individual_hbonds_ ) {
			core::Real resi_hbond_bb_sc_energy = 0;
			core::Real ref_resi_hbond_bb_sc_energy = 0;
			for ( core::Size j=1; j<=nres_asymmetric_unit; j++ ) {
				for ( Size ihb = 1; ihb <= hbset_ref.nhbonds(); ++ihb ) {
					core::scoring::hbonds::HBond const & hb(hbset_ref.hbond(ihb));
					if ( hb.don_res()==i && hb.acc_res()==j ) {
						if ( !hb.don_hatm_is_protein_backbone() && hb.acc_atm_is_protein_backbone() ) {
							ref_resi_hbond_bb_sc_energy += hb.energy();
							TR << "Refpose donor sc resi: " << i << refp.residue(i).name3() << ". Refpose acceptor bb resi: " << j << refp.residue(j).name3() << std::endl;
						}
					}
					if ( hb.don_res()==j && hb.acc_res()==i ) {
						if ( hb.don_hatm_is_protein_backbone() && !hb.acc_atm_is_protein_backbone() ) {
							ref_resi_hbond_bb_sc_energy += hb.energy();
							TR << "Refpose acceptor sc resi: " << i << refp.residue(i).name3() << ". Refpose donor bb resi: " << j << refp.residue(j).name3() << std::endl;
						}
					}
				}
				for ( Size ihb = 1; ihb <= hbset.nhbonds(); ++ihb ) {
					core::scoring::hbonds::HBond const & hb(hbset.hbond(ihb));
					if ( hb.don_res()==i && hb.acc_res()==j ) {
						if ( !hb.don_hatm_is_protein_backbone() && hb.acc_atm_is_protein_backbone() ) {
							resi_hbond_bb_sc_energy += hb.energy();
							TR << "Pose donor sc resi: " << i << refp.residue(i).name3() << ". Pose acceptor bb resi: " << j << refp.residue(j).name3() << std::endl;
						}
					}
					if ( hb.don_res()==j && hb.acc_res()==i ) {
						if ( hb.don_hatm_is_protein_backbone() && !hb.acc_atm_is_protein_backbone() ) {
							resi_hbond_bb_sc_energy += hb.energy();
							TR << "Pose acceptor sc resi: " << i << refp.residue(i).name3() << ". Pose donor bb resi: " << j << refp.residue(j).name3() << std::endl;
						}
					}
				}
			}
			delta_score = resi_hbond_bb_sc_energy - ref_resi_hbond_bb_sc_energy;
			TR << "Residue: CurrentPose - RefPose " << score_type_name_ << " : " << i << " " << resi_hbond_bb_sc_energy << " - " << ref_resi_hbond_bb_sc_energy << " = " << delta_score << std::endl;
		} else {
			ref_em = refp.energies().residue_total_energies( i );
			ref_em *= scorefxn()->weights();
			em = p.energies().residue_total_energies( i );
			em *= scorefxn()->weights();
			delta_score = em[score_type_] - ref_em[score_type_];
			TR << "Residue: CurrentPose - RefPose " << score_type_name_ << " : " << i << " " << em[score_type_] << " - " << ref_em[score_type_] << " = " << delta_score << std::endl;
		}
		if ( (delta_score <= threshold_ && lower_ == false) || ( delta_score >= threshold_ && lower_ == true) ) {
			task.nonconst_residue_task(i).prevent_repacking();
		} else {
			core::Size output_resi = i;
			if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
				output_resi = pose.pdb_info()->number( i );
			}
			pymol_selection.append(ObjexxFCL::string_of(output_resi) + "+");
		}
	}
	TR << pymol_selection << std::endl;
}

void
SelectByDeltaScoreOperation::parse_tag( TagCOP tag , DataMap & data)
{
	scorefxn(protocols::rosetta_scripts::parse_score_function(tag, data));
	score_type(core::scoring::score_type_from_name( tag->getOption<std::string>( "score_type", "total_score" ) ) );
	score_type_name( tag->getOption<std::string>( "score_type", "total_score" ) );
	threshold(tag->getOption< core::Real >("threshold", 100 ));
	lower(tag->getOption< bool >("lower", false ));
	individual_hbonds(tag->getOption< bool >("individual_hbonds", false ));
	if ( tag->hasOption("reference_name") ) {
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag,data );
	} else if ( tag->hasOption("reference_pdb") ) {
		std::string reference_pdb_filename( tag->getOption< std::string >( "reference_pdb", "" ) );
		reference_pose_ = core::import_pose::pose_from_file( reference_pdb_filename , core::import_pose::PDB_file);
	}
}

void SelectByDeltaScoreOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	rosetta_scripts::attributes_for_parse_score_function( attributes );

	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "score_type", xs_string, "Loop over the independent residues, assess score_type for each residue.",  "total_score"  )
		+ XMLSchemaAttribute::attribute_w_default(  "threshold", xsct_real, "Restrict residues to repack based on whether they pass this threshold.",  "100"  )
		+ XMLSchemaAttribute::attribute_w_default(  "lower", xsct_rosetta_bool, "Higher or lower than threshold.",  "false"  )
		+ XMLSchemaAttribute::attribute_w_default(  "individual_hbonds", xsct_rosetta_bool, "If true, get Hbond set.",  "false"  )
		+ XMLSchemaAttribute( "reference_name", xs_string , "Name of reference pose." )
		+ XMLSchemaAttribute( "reference_pdb", xs_string , "Reference pdb filename." );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "Restrict design to residues passing a user-specified threshold on a given score type." );
}


} //namespace task_operations
} //namespace protocols
