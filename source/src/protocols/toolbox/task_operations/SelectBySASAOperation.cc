// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief  Restrict design to residues matching user-specified SASA criteria in the monomeric, bound, or unbound state.
/// @author Jacob Bale (balej@uw.edu)

// Unit Headers
#include <protocols/toolbox/task_operations/SelectBySASAOperation.hh>
#include <protocols/toolbox/task_operations/SelectBySASAOperationCreator.hh>

// Project Headers
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/sasa.hh>
#include <core/types.hh>
#include <core/id/AtomID_Map.hh>
#include <protocols/rigid/RigidBodyMover.hh>
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

static THREAD_LOCAL basic::Tracer TR( "protocols.toolbox.task_operations.SelectBySASAOperation" );

namespace protocols {
namespace toolbox {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace utility::tag;

core::pack::task::operation::TaskOperationOP
SelectBySASAOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new SelectBySASAOperation );
}

void SelectBySASAOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SelectBySASAOperation::provide_xml_schema( xsd );
}

std::string SelectBySASAOperationCreator::keyname() const
{
	return SelectBySASAOperation::keyname();
}


SelectBySASAOperation::SelectBySASAOperation( std::string mode, std::string state, core::Real probe_radius, core::Real core_asa, core::Real surface_asa, std::string jump_nums, std::string sym_dof_names, bool core, bool boundary, bool surface, bool verbose ):
	mode_(mode),
	state_(state),
	probe_radius_(probe_radius),
	core_asa_(core_asa),
	surface_asa_(surface_asa),
	jump_nums_(jump_nums),
	sym_dof_names_(sym_dof_names),
	core_(core),
	boundary_(boundary),
	surface_(surface),
	verbose_(verbose)
{}

SelectBySASAOperation::~SelectBySASAOperation() {}

core::pack::task::operation::TaskOperationOP SelectBySASAOperation::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new SelectBySASAOperation( *this ) );
}

void
SelectBySASAOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using core::id::AtomID;
	std::string layername = "";
	bool prev = 0;
	if ( core_ ) {
		layername.append("core");
		prev = 1;
	}
	if ( boundary_ ) {
		if ( prev ) {
			layername.append("_boundary");
		} else {
			layername.append("boundary");
		}
		prev = 1;
	}
	if ( surface_ ) {
		if ( prev ) {
			layername.append("_surface");
		} else {
			layername.append("surface");
		}
		prev = 1;
	}

	if ( !prev ) {
		utility_exit_with_message("The layers are not set properly. Core, boundary, and surface are all set to false. At least one layer needs to be selected.");
	}

	std::string selected_pos("select " + layername + ", resi ");
	std::string core_pos("select core, resi ");
	std::string boundary_pos("select boundary, resi ");
	std::string surface_pos("select surface, resi ");
	core::pose::Pose mono, sasa_pose;

	core::Size nsubposes = 1;
	if ( state_ == "monomer" ) {
		if ( core::pose::symmetry::is_symmetric(pose) ) {
			core::pose::symmetry::extract_asymmetric_unit( pose, mono , false );
		} else {
			mono = pose;
		}
		nsubposes = mono.conformation().num_chains();
	}

	core::Size res_count = 0;
	for ( core::Size i = 1; i <= nsubposes; i++ ) {
		utility::vector1<bool> indy_resi;

		if ( state_ == "monomer" ) {
			sasa_pose = *mono.split_by_chain(i);
		} else {
			sasa_pose = pose;
			if ( core::pose::symmetry::is_symmetric(sasa_pose) ) {
				indy_resi = core::pose::symmetry::symmetry_info(sasa_pose)->independent_residues();
			}
		}

		if ( state_ == "unbound" ) {
			int sym_aware_jump_id = 0;
			if ( sym_dof_names_ != "" ) {
				utility::vector1<std::string> sym_dof_name_list = utility::string_split( sym_dof_names_ , ',' );
				for ( Size i = 1; i <= sym_dof_name_list.size(); i++ ) {
					sym_aware_jump_id = core::pose::symmetry::sym_dof_jump_num( sasa_pose, sym_dof_name_list[i] );
					protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( sasa_pose, sym_aware_jump_id ) );
					translate->step_size( 1000.0 );
					translate->apply( sasa_pose );
				}
			} else {
				utility::vector1<core::Size> jump_list = utility::string_split( jump_nums_ , ',',core::Size());
				for ( Size i = 1; i <= jump_list.size(); i++ ) {
					sym_aware_jump_id = core::pose::symmetry::get_sym_aware_jump_num( sasa_pose, jump_list[i] );
					protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( sasa_pose, sym_aware_jump_id ) );
					translate->step_size( 1000.0 );
					translate->apply( sasa_pose );
				}
			}
		}

		//Loop over the independent residues or the residues in the monomer(s), calculate SASA for each, and assign each to core, boundary, or surface.
		// Define atom_map for main-chain and CB, or for sidechains depending on mode.
		core::id::AtomID_Map< bool > atom_mask;
		core::id::AtomID_Map< core::Real > atom_sasa;
		utility::vector1< core::Real > sasas, final_sasas;
		core::pose::initialize_atomid_map( atom_mask, sasa_pose, false );
		core::pose::initialize_atomid_map( atom_sasa, sasa_pose, 0.0);
		core::Size itype;

		for ( core::Size i = 1; i <= sasa_pose.size(); ++i ) {
			if ( sasa_pose.residue(i).name3()=="GLY" ) { // Don't try CBs with GLYs
				itype = 4;
			} else {
				itype = 5;
			}
			for ( core::Size j = 1; j<=sasa_pose.residue(i).nheavyatoms(); ++j ) {
				if ( (mode_ == "mc") && (j > itype) ) {
					continue;
				} else {
					core::id::AtomID atom( j, i );
					atom_mask.set( atom, true );
				}
			}
		}

		// calc sasa
		core::scoring::calc_per_atom_sasa( sasa_pose, atom_sasa, sasas, probe_radius_, false, atom_mask );

		if ( mode_ == "sc" ) {
			utility::vector1<core::Real> sc_sasas(sasa_pose.size(),0.0);
			for ( Size i = 1; i <= sasa_pose.size(); i++ ) {
				// Use CA as the side chain for Glys
				if ( sasa_pose.residue(i).name3()=="GLY" ) sc_sasas[i] += atom_sasa[AtomID(2,i)];
				for ( Size j = 5; j <= sasa_pose.residue(i).nheavyatoms(); j++ ) {
					sc_sasas[i] += atom_sasa[AtomID(j,i)];
				}
			}
			final_sasas = sc_sasas;
		} else {
			final_sasas = sasas;
		}

		// Prevent repacking at resis that do match the user-specified parameters.
		//bool prevent_repacking;
		for ( Size iaa=1; iaa<=sasa_pose.size(); iaa++ ) {
			if ( core::pose::symmetry::is_symmetric(sasa_pose) ) {
				if ( !indy_resi[iaa] ) {
					continue;
				}
			}
			if ( sasa_pose.residue( iaa ).is_protein() ) {
				bool prevent_repacking = 1;
				res_count++;
				core::Size output_resi = res_count;
				if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
					if ( pose.pdb_info() ) {
						output_resi = pose.pdb_info()->number( res_count );
					}
				}
				TR.Debug << iaa << " res_count = " << res_count << " sasa = " << final_sasas[iaa] << std::endl;
				if ( final_sasas[ iaa ] <= core_asa_ ) {
					if ( core_ ) {
						selected_pos.append(ObjexxFCL::string_of(output_resi) + "+");
						prevent_repacking = 0;
					}
					core_pos.append(ObjexxFCL::string_of(output_resi) + "+");
				} else if ( final_sasas[iaa] >= surface_asa_ ) {
					if ( surface_ ) {
						selected_pos.append(ObjexxFCL::string_of(output_resi) + "+");
						prevent_repacking = 0;
					}
					surface_pos.append(ObjexxFCL::string_of(res_count) + "+");
				} else if ( (final_sasas[iaa] >= core_asa_) && (final_sasas[iaa] <= surface_asa_) ) {
					if ( boundary_ ) {
						selected_pos.append(ObjexxFCL::string_of(output_resi) + "+");
						prevent_repacking = 0;
					}
					boundary_pos.append(ObjexxFCL::string_of(output_resi) + "+");
				}
				if ( prevent_repacking ) {
					task.nonconst_residue_task(res_count).prevent_repacking();
				}
			}
		}
	}
	if ( verbose_ ) {
		TR << selected_pos << std::endl;
		TR << core_pos << std::endl;
		TR << boundary_pos << std::endl;
		TR << surface_pos << std::endl;
	}
}

void
SelectBySASAOperation::parse_tag( TagCOP tag , DataMap & )
{
	mode_ = tag->getOption< std::string >("mode", "sc" );
	state_ = tag->getOption< std::string >("state", "monomer" );
	probe_radius_ = tag->getOption<core::Real>("probe_radius", 2.2);
	core_asa_ = tag->getOption<core::Real>("core_asa", 0);
	surface_asa_ = tag->getOption<core::Real>("surface_asa", 30);
	jump_nums_ = tag->getOption<std::string>("jumps", "1");
	sym_dof_names_ = tag->getOption<std::string>("sym_dof_names","");
	core_ = tag->getOption< bool >("core", 0 );
	boundary_ = tag->getOption< bool >("boundary", 0 );
	surface_ = tag->getOption< bool >("surface", 0 );
	verbose_ = tag->getOption< bool >("verbose", 0 );
}

void SelectBySASAOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "mode", xs_string, "sc" )
		+ XMLSchemaAttribute::attribute_w_default(  "state", xs_string, "monomer" )
		+ XMLSchemaAttribute::attribute_w_default(  "probe_radius", xs_decimal, "2.2" )
		+ XMLSchemaAttribute::attribute_w_default(  "core_asa", xs_decimal, "0" )
		+ XMLSchemaAttribute::attribute_w_default(  "surface_asa", xs_decimal, "30" )
		+ XMLSchemaAttribute::attribute_w_default(  "jumps", xsct_int_cslist, "1" )
		+ XMLSchemaAttribute::attribute_w_default(  "sym_dof_names", xs_string, "" )
		+ XMLSchemaAttribute::attribute_w_default(  "core", xs_boolean, "false" )
		+ XMLSchemaAttribute::attribute_w_default(  "boundary", xs_boolean, "false" )
		+ XMLSchemaAttribute::attribute_w_default(  "surface", xs_boolean, "false" )
		+ XMLSchemaAttribute::attribute_w_default(  "verbose", xs_boolean, "false" );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}


} //namespace task_operations
} //namespace toolbox
} //namespace protocols
