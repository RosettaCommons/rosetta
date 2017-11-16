// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/InsertResMover.fwd.hh
/// @details inserts ideal residues into pose. Useful for extending helices
///
/// @author TJ Brunette tjbrunette@gmail.com
///
// Unit headers
#include <protocols/pose_length_moves/InsertResMover.hh>
#include <protocols/pose_length_moves/InsertResMoverCreator.hh>
#include <protocols/moves/Mover.hh>


// Core Headers
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <boost/assign/list_of.hpp>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.pose_length_moves.InsertResMover" );

namespace protocols {
namespace pose_length_moves {
using namespace core;
using namespace std;
using utility::vector1;



InsertResMover::InsertResMover():moves::Mover("InsertResMover"){}

// XRW TEMP std::string InsertResMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return InsertResMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string InsertResMover::mover_name(){
// XRW TEMP  return "InsertResMover";
// XRW TEMP }


// XRW TEMP protocols::moves::MoverOP
// XRW TEMP InsertResMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new InsertResMover );
// XRW TEMP }

numeric::xyzVector<core::Real>
InsertResMover::center_of_mass(core::pose::Pose const & pose) {
	int nAtms = 0;
	numeric::xyzVector<core::Real> massSum(0.,0.,0.), CoM;
	for ( core::Size ii =1; ii <= pose.size(); ++ii ) {
		if ( pose.residue_type(ii).aa() == core::chemical::aa_vrt ) continue;
		for ( core::Size iatom = 1; iatom <= pose.residue_type(ii).nheavyatoms(); ++iatom ) {
			core::conformation::Atom const & atom( pose.residue(ii).atom(iatom) );
			massSum += atom.xyz();
			nAtms++;
		}
	}
	CoM = massSum / (core::Real)nAtms;
	return CoM;
}


void InsertResMover::extendRegion(core::pose::PoseOP poseOP, Size chain_id, Size length){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace chemical;
	core::conformation::ResidueCOPs residues(core::pose::get_chain_residues(*poseOP, chain_id));
	Size res_position = residue_;
	if ( residue_>poseOP->size() ) {
		res_position=poseOP->size();
	}
	Size inPoseResidue = residues[res_position]->seqpos();
	Real tmpPhi =  poseOP->phi(inPoseResidue);
	Real tmpPsi =  poseOP->psi(inPoseResidue);
	Real tmpOmega = poseOP->omega(inPoseResidue);
	if ( ideal_ ) {
		if ( resType_ == "E" || resType_ == "L" ) {
			tmpPhi = -135;
			tmpPsi = 135;
			tmpOmega = 180;
		}
		if ( resType_ == "H" ) {
			tmpPhi = -57.8;
			tmpPsi = -47.0;
			tmpOmega = 180.0;
		}
	}
	if ( useInputAngles_ ) {
		tmpPhi = phi_;
		tmpPsi = psi_;
		tmpOmega = omega_;
	}
	if ( steal_angles_from_res_!= 0 ) {
		tmpPhi =  poseOP->phi(steal_angles_from_res_);
		tmpPsi =  poseOP->psi(steal_angles_from_res_);
		tmpOmega = poseOP->omega(steal_angles_from_res_);
	}
	kinematics::FoldTree backupFt = poseOP->fold_tree();
	core::conformation::ResidueOP new_rsd( nullptr );
	string build_aa_type_one_letter =option[OptionKeys::remodel::generic_aa];
	string build_aa_type = name_from_aa(aa_from_oneletter_code(build_aa_type_one_letter[0]));
	debug_assert( poseOP != nullptr );
	core::chemical::ResidueTypeSetCOP rs( poseOP->residue_type_set_for_pose() );
	kinematics::FoldTree ft;
	//Code below is a bunch of junk. I need a uniform way to deal with broken fold trees
	if ( inPoseResidue != 1 && inPoseResidue!=poseOP->total_residue() ) {
		numeric::xyzVector<core::Real> CoM;
		CoM = center_of_mass(*poseOP);
		new_rsd = core::conformation::ResidueFactory::create_residue( rs->name_map("VRT") );
		new_rsd->atom(1).xyz(CoM);
		//Set up fold tree for this one case --------------------------------------------------
		poseOP->append_residue_by_jump( *new_rsd,poseOP->total_residue());
		if ( grow_toward_Nterm_ ) {
			ft.add_edge(1,inPoseResidue-1,core::kinematics::Edge::PEPTIDE);
			ft.add_edge(1,poseOP->size(),1);
			ft.add_edge(poseOP->size(),inPoseResidue,core::kinematics::Edge::PEPTIDE);
		} else {
			ft.add_edge(1,inPoseResidue,core::kinematics::Edge::PEPTIDE);
			ft.add_edge(1,poseOP->size(),1);
			ft.add_edge(poseOP->size(),inPoseResidue+1,core::kinematics::Edge::PEPTIDE);
		}
		poseOP->fold_tree(ft);
	} else {
		ft.add_edge(1,poseOP->size(),core::kinematics::Edge::PEPTIDE);
		poseOP->fold_tree(ft);
	}
	poseOP->fold_tree(ft);
	if ( grow_toward_Nterm_ ) {
		for ( Size ii=0; ii<length; ++ii ) {
			new_rsd = core::conformation::ResidueFactory::create_residue( rs->name_map(build_aa_type) );
			poseOP->conformation().safely_prepend_polymer_residue_before_seqpos( *new_rsd,inPoseResidue, true);
		}
	} else {
		for ( Size ii=0; ii<length; ++ii ) {
			new_rsd = core::conformation::ResidueFactory::create_residue( rs->name_map(build_aa_type) );
			poseOP->conformation().safely_append_polymer_residue_after_seqpos( *new_rsd,inPoseResidue+ii, true);
		}
	}
	renumber_pdbinfo_based_on_conf_chains(*poseOP,true,false,false,false);
	for ( Size ii=0; ii<=length; ++ii ) {
		poseOP->set_phi(inPoseResidue+ii, tmpPhi );
		poseOP->set_psi(inPoseResidue+ii, tmpPsi );
		poseOP->set_omega(inPoseResidue+ii, tmpOmega );
	}
}


void InsertResMover::apply(core::pose::Pose & pose) {
	Size chain_id;
	if ( chain_=="999" ) {
		utility::vector1< Size > chains = get_chains(pose );
		chain_id = chains[1];
	} else {
		if ( !has_chain(chain_,pose) ) {
			utility_exit_with_message( "Chain does not exit:"+chain_);
		}
		chain_id = get_chain_id_from_chain(chain_,pose);
	}
	for ( Size length=lowAddRes_; length<=highAddRes_; ++length ) {
		core::pose::PoseOP tmpPoseOP = pose.clone();
		extendRegion(tmpPoseOP,chain_id,length);
		posesToOutput_.push_back(tmpPoseOP);
		posesOutputed_.push_back(false);
	}
	core::pose::PoseOP poseOP = get_additional_output();
	if ( poseOP==nullptr ) {
		TR << "no succeful extensions" << std::endl;
		set_last_move_status(protocols::moves::FAIL_RETRY);
	} else {
		pose = *poseOP;
	}
}


// XRW TEMP std::string InsertResMover::get_name() const {
// XRW TEMP  return "InsertResMover";
// XRW TEMP }

void
InsertResMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ){
	chain_ = ( tag->getOption< std::string >( "chain", "999") );
	resType_ = ( tag->getOption< std::string >( "resType", "H") );
	steal_angles_from_res_= (tag->getOption<Size>("steal_angles_from_res",0));
	residue_ = ( tag->getOption< Size >( "residue", 1) );
	grow_toward_Nterm_ = ( tag->getOption< bool >( "grow_toward_Nterm", false) );
	ideal_ = ( tag->getOption< bool >( "ideal", true) );
	useInputAngles_ = ( tag->getOption< bool >( "ideal", false) );
	phi_ = ( tag->getOption< Real >( "phi", -57.8 ) );
	psi_ = ( tag->getOption< Real >( "psi", -47.0 ) );
	omega_ = ( tag->getOption< Real >( "omega", 180.0 ) );
	if ( ideal_==false ) {
		TR << "stealing phi,psi,omega from residue location" << std::endl;
	}
	std::string additionResString( tag->getOption< std::string >( "additionalResidue", "1") ); //would be valid to do 1-4
	utility::vector1< std::string > additionRes_split( utility::string_split( additionResString , ',' ) );
	if ( additionRes_split.size()==2 ) {
		lowAddRes_ = atoi(additionRes_split[1].c_str());
		highAddRes_ = atoi(additionRes_split[2].c_str());
	}
	if ( additionRes_split.size()==1 ) {
		lowAddRes_ = atoi(additionRes_split[1].c_str());
		highAddRes_ = atoi(additionRes_split[1].c_str());
	}
	TR << "additionalRes" << lowAddRes_  <<"," << highAddRes_ << std::endl;
}

core::pose::PoseOP InsertResMover::get_additional_output(){
	TR << "NOT YET IMPLEMENTED CORRECTLY WILL THROW OFF SELECTOR****************************" << std::endl;
	for ( Size ii=1; ii<=posesToOutput_.size(); ii++ ) {
		if ( !posesOutputed_[ii] ) {
			posesOutputed_[ii]=true;
			return(posesToOutput_[ii]);
		}
	}
	set_last_move_status(protocols::moves::FAIL_RETRY);
	return nullptr;
}

std::string InsertResMover::get_name() const {
	return mover_name();
}

std::string InsertResMover::mover_name() {
	return "InsertResMover";
}

void InsertResMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default(
		"chain", xs_string,
		"XSD_XRW: TO DO",
		"999");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"resType", xs_string,
		"XSD_XRW: TO DO",
		"H");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"steal_angles_from_res", xsct_non_negative_integer,
		"XSD_XRW: TO DO",
		"0");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"residue", xsct_non_negative_integer,
		"XSD_XRW: TO DO",
		"1");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"grow_toward_Nterm", xsct_rosetta_bool,
		"XSD_XRW: TO DO",
		"false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"ideal", xsct_rosetta_bool,
		"XSD_XRW: TO DO",
		"true");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"phi", xsct_real,
		"XSD_XRW: TO DO",
		"-57.8");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"psi", xsct_real,
		"XSD_XRW: TO DO",
		"-47.0");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"omega", xsct_real,
		"XSD_XRW: TO DO",
		"180.0");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"additionalResidue", xs_string,
		"XSD_XRW: TO DO",
		"1");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Inserts ideal residues into pose. Useful for extending helices",
		attlist );
}

std::string InsertResMoverCreator::keyname() const {
	return InsertResMover::mover_name();
}

protocols::moves::MoverOP
InsertResMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new InsertResMover );
}

void InsertResMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InsertResMover::provide_xml_schema( xsd );
}


}//pose_length_moves
}//protocols
