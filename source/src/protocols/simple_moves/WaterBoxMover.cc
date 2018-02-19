// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file WaterBoxMover.cc

// Unit headers
#include <protocols/simple_moves/WaterBoxMover.hh>
#include <protocols/simple_moves/WaterBoxMoverCreator.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator_Private.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/conformation/symmetry/SymmetricConformation.fwd.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/packer_neighbors.hh>

#include <basic/Tracer.hh>
#include <fstream>
#include <iostream>

using basic::Error;
using basic::Warning;
static basic::Tracer TR("protocols.moves.WaterBoxMover");

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <core/pose/chains_util.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

namespace protocols {
namespace simple_moves {

std::string
WaterBoxMoverCreator::keyname() const {
	return WaterBoxMoverCreator::mover_name();
}

protocols::moves::MoverOP
WaterBoxMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new WaterBoxMover );
}

std::string
WaterBoxMover::mover_name() {
	return WaterBoxMoverCreator::mover_name();
}

std::string
WaterBoxMoverCreator::mover_name() {
	return "WaterBoxMover";
}

void WaterBoxMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	WaterBoxMover::provide_xml_schema( xsd );
}

void WaterBoxMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	auto ct_gen = define_water_box_mover_schema();
	ct_gen->write_complex_type_to_schema(xsd);
}

WaterBoxMover::WaterBoxMover()
: Mover("WaterBoxMover") {
	point_waters_ = true;
	make_rotatable_ = false;
	make_point_ = false;
	remove_all_ = false;
	bb_sol_ = false;
	sc_sol_ = false;
	lig_sol_ = false;
}

std::string
WaterBoxMover::get_name() const {
	return WaterBoxMoverCreator::mover_name();
}

protocols::moves::MoverOP
WaterBoxMover::clone() const {
	return protocols::moves::MoverOP( new WaterBoxMover( *this ) );
}

protocols::moves::MoverOP
WaterBoxMover::fresh_instance() const {
	return protocols::moves::MoverOP( new WaterBoxMover );
}


// helper function to delete virtualized waters from a pose
void
WaterBoxMover::delete_virtual_waters( core::pose::Pose  & pose ) {
	core::Size nres = pose.total_residue();

	for ( core::Size i=nres; i>=1; --i ) {   // don't worry about renumbering
		core::conformation::Residue const &res = pose.residue(i);
		if ( res.name() == "HOH_V" || res.name() == "PWAT_V" ) {
			pose.conformation().delete_residue_slow( i );
		} else if ( remove_all_ && (res.name() == "HOH" || res.name() == "PWAT") ) {
			pose.conformation().delete_residue_slow( i );
		}
	}
	renumber_pdbinfo_based_on_conf_chains(pose);
}


// build initial point waters of PWAT exactly on top of
// backbone N or O atoms -- use the proximity of the point
// waters to backbone atoms to determine where to build
// the rotamers for each backbone polar gropu in RotamerSet_.cc
void
WaterBoxMover::bb_sol( Pose & pose, core::pack::task::PackerTaskCOP task ) {
	if ( task_ == 0 ) {
		task_ = task;
		TR << "setting packer task from argument for bb_sol" << std::endl;
	}
	bb_sol( pose );
}

void
WaterBoxMover::bb_sol( Pose & pose ) {
	core::chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );
	core::Size nres = pose.total_residue();
	for ( int i=1; i<=(int)nres; ++i ) {
		core::conformation::Residue const &res = pose.residue(i);
		if ( res.is_water() ) continue;
		if ( task_->pack_residue(i) || task_->design_residue(i) ) {
			core::conformation::ResidueOP vrt_wat = core::conformation::ResidueFactory::create_residue( rsd_set->name_map("PWAT") );
			core::Vector const OH1( vrt_wat->xyz("V1") - vrt_wat->xyz("O") );
			core::Vector const OH2( vrt_wat->xyz("V2") - vrt_wat->xyz("O") );
			core::conformation::ResidueOP new_res = core::conformation::ResidueOP( new core::conformation::Residue( *vrt_wat ) );

			if ( res.is_protein() ) {

				// build point water on top of O atom
				Size aatm = res.atom_index("O");
				new_res->set_xyz("O",  res.xyz(aatm));
				new_res->set_xyz("V1", res.xyz(aatm)+OH1);
				new_res->set_xyz("V2", res.xyz(aatm)+OH2);
				pose.append_residue_by_jump( *new_res, i );
				TR << "solvating backbone of residue: " << i << " " << res.name() << std::endl;

				// special cases for termini including chain breaks
				if ( i>1 ) {
					core::conformation::Residue const &prev_res = pose.residue(i-1);
					// use prev_res check for :CtermTruncation cases that don't have OXT atom
					if ( res.is_upper_terminus() && ( prev_res.chain() != res.chain() ) ) {
						Size aatm = res.atom_index("OXT");
						new_res->set_xyz("O",  res.xyz(aatm));
						new_res->set_xyz("V1", res.xyz(aatm)+OH1);
						new_res->set_xyz("V2", res.xyz(aatm)+OH2);
						pose.append_residue_by_jump( *new_res, i );
					}
				}

				if ( i<(int)nres ) {
					core::conformation::Residue const &next_res = pose.residue(i+1);
					// use next_res chain check for :NtermTruncation cases that don't have 1H, 2H and 3H atoms
					if ( res.is_lower_terminus() && ( next_res.chain() == res.chain() ) ) {
						Size aatm = res.atom_index("1H");
						new_res->set_xyz("O",  res.xyz(aatm));
						new_res->set_xyz("V1", res.xyz(aatm)+OH1);
						new_res->set_xyz("V2", res.xyz(aatm)+OH2);
						pose.append_residue_by_jump( *new_res, i );
						aatm = res.atom_index("2H");
						new_res->set_xyz("O",  res.xyz(aatm));
						new_res->set_xyz("V1", res.xyz(aatm)+OH1);
						new_res->set_xyz("V2", res.xyz(aatm)+OH2);
						pose.append_residue_by_jump( *new_res, i );
						if ( res.name3() == "PRO" ) continue; // 3H is not found in N-terminal PRO
						aatm = res.atom_index("3H");
						new_res->set_xyz("O",  res.xyz(aatm));
						new_res->set_xyz("V1", res.xyz(aatm)+OH1);
						new_res->set_xyz("V2", res.xyz(aatm)+OH2);
						pose.append_residue_by_jump( *new_res, i );
					} else if ( res.name3() != "PRO" && ( !res.is_lower_terminus() ||
							( res.is_lower_terminus() && ( next_res.chain() != res.chain() ) ) ) ) {
						// build point water on top of H atom
						Size hatm = res.atom_index("H");
						new_res->set_xyz("O",  res.xyz(hatm));
						new_res->set_xyz("V1", res.xyz(hatm)+OH1);
						new_res->set_xyz("V2", res.xyz(hatm)+OH2);
						pose.append_residue_by_jump( *new_res, i );
					}
				}
			}
		}
	}
}

void
WaterBoxMover::sc_sol( Pose & pose, core::pack::task::PackerTaskCOP task ) {
	if ( task_ == 0 ) {
		task_ = task;
		TR << "setting packer task from argument for sc_sol" << std::endl;
	}
	sc_sol( pose );
}

void
WaterBoxMover::sc_sol( Pose & pose ) {
	core::chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );
	core::Size nres = pose.total_residue();
	for ( int i=1; i<=(int)nres; ++i ) {
		core::conformation::Residue const &res = pose.residue(i);
		if ( res.is_water() ) continue;
		if ( task_->pack_residue(i) || task_->design_residue(i) ) {
			core::conformation::ResidueOP vrt_wat = core::conformation::ResidueFactory::create_residue( rsd_set->name_map("PWAT") );
			core::Vector const OH1( vrt_wat->xyz("V1") - vrt_wat->xyz("O") );
			core::Vector const OH2( vrt_wat->xyz("V2") - vrt_wat->xyz("O") );
			core::conformation::ResidueOP new_res = core::conformation::ResidueOP( new core::conformation::Residue( *vrt_wat ) );

			if ( res.is_protein() ) {
				// build point water on top of CB atom
				if ( res.name3() == "GLY" ) continue;
				Size aatm = res.atom_index("CB");
				new_res->set_xyz("O",  res.xyz(aatm));
				new_res->set_xyz("V1", res.xyz(aatm)+OH1);
				new_res->set_xyz("V2", res.xyz(aatm)+OH2);
				pose.append_residue_by_jump( *new_res, i );
			}
		}
	}
}

void
WaterBoxMover::lig_sol( Pose & pose ) {
	core::chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );
	core::Size nres = pose.total_residue();
	for ( int i=1; i<=(int)nres; ++i ) {
		core::conformation::Residue const &res = pose.residue(i);
		if ( task_->pack_residue(i) || task_->design_residue(i) ) {
			core::conformation::ResidueOP vrt_wat = core::conformation::ResidueFactory::create_residue( rsd_set->name_map("PWAT") );
			core::Vector const OH1( vrt_wat->xyz("V1") - vrt_wat->xyz("O") );
			core::Vector const OH2( vrt_wat->xyz("V2") - vrt_wat->xyz("O") );
			core::conformation::ResidueOP new_res = core::conformation::ResidueOP( new core::conformation::Residue( *vrt_wat ) );

			// ligands: all donors and acceptors
			if ( !res.is_protein() && res.aa() != core::chemical::aa_vrt ) {
				// acceptor atoms
				for ( core::chemical::AtomIndices::const_iterator
						anum  = res.accpt_pos().begin(),
						anume = res.accpt_pos().end(); anum != anume; ++anum ) {
					core::Size const aatm( *anum );
					new_res->set_xyz("O",  res.xyz(aatm));
					new_res->set_xyz("V1", res.xyz(aatm)+OH1);
					new_res->set_xyz("V2", res.xyz(aatm)+OH2);
					pose.append_residue_by_jump( *new_res, i );
				}
				// donor atoms
				for ( core::chemical::AtomIndices::const_iterator
						hnum  = res.Hpos_polar().begin(),
						hnume = res.Hpos_polar().end(); hnum != hnume; ++hnum ) {
					core::Size hatm(*hnum);
					new_res->set_xyz("O",  res.xyz(hatm));
					new_res->set_xyz("V1", res.xyz(hatm)+OH1);
					new_res->set_xyz("V2", res.xyz(hatm)+OH2);
					pose.append_residue_by_jump( *new_res, i );
				}
			}
		}
	}
}

void
WaterBoxMover::apply( Pose & pose ) {
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::pack::task;

	core::chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );

	// score pose to make sure energy graph is correct
	core::scoring::ScoreFunctionOP test_score = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
	test_score->set_weight(core::scoring::fa_rep, 0.1);
	(*test_score)(pose);

	// if present, task_factory_ always overrides/regenerates task_
	if ( task_factory_ != 0 ) {
		task_ = task_factory_->create_task_and_apply_taskoperations( pose );
	} else if ( task_ == 0 ) {
		task_ = TaskFactory::create_packer_task( pose );
	}

	// clean up
	if ( del_virt_ ) {
		delete_virtual_waters(pose);
	}

	// solvate backbone with point water residues
	if ( bb_sol_ ) {
		bb_sol(pose);
	}

	if ( sc_sol_ ) {
		sc_sol(pose);
	}

	if ( lig_sol_ ) {
		lig_sol(pose);
	}

	if ( del_virt_ || remove_all_ || bb_sol_ || sc_sol_ || lig_sol_ ) return;

	core::Size nres = pose.total_residue();
	// special case -- if this flag is given, convert all point waters to rotatable waters
	if ( make_rotatable_ || make_point_ ) {
		for ( core::Size i=1; i<=nres; ++i ) {
			core::conformation::Residue const &res = pose.residue(i);
			if ( res.name() == "PWAT" || res.name() == "HOH" ) {
				ResidueOP rsd_new = core::conformation::ResidueFactory::create_residue( rsd_set->name_map(make_point_?"PWAT":"HOH") );
				core::Vector OH1( rsd_new->xyz(make_point_?"V1":"H1") - rsd_new->xyz("O") );
				core::Vector OH2( rsd_new->xyz(make_point_?"V2":"H2") - rsd_new->xyz("O") );
				core::Vector O = res.xyz("O");
				rsd_new->set_xyz("O", O);
				rsd_new->set_xyz(make_point_?"V1":"H1", O+OH1);
				rsd_new->set_xyz(make_point_?"V2":"H2", O+OH2);
				pose.replace_residue( i, *rsd_new, false );
			}
		}
		return;
	}
	renumber_pdbinfo_based_on_conf_chains(pose); // rep: update pdb_info to account for added wate  s
}

void
WaterBoxMover::add_water( Pose & pose, core::Vector O, core::Size resnum ) {
	using namespace core::conformation;
	using namespace core::pack::task;

	core::chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );

	// canonical water
	ResidueOP vrt_wat = point_waters_ ?
		ResidueFactory::create_residue( rsd_set->name_map("PWAT_V") ) :
		ResidueFactory::create_residue( rsd_set->name_map("HOH_V") );
	std::string name1 = point_waters_ ? "V1" : "H1";
	std::string name2 = point_waters_ ? "V2" : "H2";
	core::Vector const OH1( vrt_wat->xyz(name1) - vrt_wat->xyz("O") );
	core::Vector const OH2( vrt_wat->xyz(name2) - vrt_wat->xyz("O") );

	ResidueOP new_res = ResidueOP( new Residue( *vrt_wat ) );
	new_res->set_xyz("O", O);
	new_res->set_xyz(name1, O+OH1);
	new_res->set_xyz(name2, O+OH2);
	pose.append_residue_by_jump( *new_res, resnum );
	renumber_pdbinfo_based_on_conf_chains(pose);
}


void
WaterBoxMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &  ) {

	using namespace core::conformation;
	using namespace core::pack::task;

	task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, datamap );
	point_waters_ = tag->getOption<bool>("point_waters",true);
	make_rotatable_ = tag->getOption<bool>("make_rotatable",false);
	make_point_ = tag->getOption<bool>("make_point",false);
	del_virt_ = tag->getOption<bool>("delete_virtual",true);
	remove_all_ = tag->getOption<bool>("remove_all",false);
	bb_sol_ = tag->getOption<bool>("bb_sol",false);
	sc_sol_ = tag->getOption<bool>("sc_sol",false);
	lig_sol_ = tag->getOption<bool>("lig_sol",false);
}

utility::tag::XMLSchemaComplexTypeGeneratorOP
WaterBoxMover::define_water_box_mover_schema() {
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default(
		"remove_all", xsct_rosetta_bool,
		"XSD XRW TO DO",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"bb_sol", xsct_rosetta_bool,
		"XSD XRW TO DO",
		"false");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"sc_sol", xsct_rosetta_bool,
		"XSD XRW TO DO",
		"false");

	rosetta_scripts::attributes_for_parse_task_operations( attlist );
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );

	ct_gen->complex_type_naming_func(&moves::complex_type_name_for_mover)
		.element_name( mover_name() )
		.description(
		"This mover cad solvate a pose in one of two ways: "
		"1.) a uniform grid within a specified distance (\"shell\") of selected residues; "
		"2.) based on statistics of where waters are found in high resolution crystal "
		"structures about polar side chains and backbone groups.")
		.add_attributes( attlist )
		.add_optional_name_attribute();
	return ct_gen;
}

} // moves
} // protocols
