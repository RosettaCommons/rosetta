// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/ParatopeSiteConstraintMover.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.hh>
#include <protocols/antibody/constraints/ParatopeSiteConstraintMoverCreator.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>

#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/SiteConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR("protocols.antibody.constraints.ParatopeSiteConstraintMover");

namespace protocols {
namespace antibody {
namespace constraints {
using utility::vector1;
ParatopeSiteConstraintMover::ParatopeSiteConstraintMover() :
	protocols::moves::Mover(),
	ab_info_(/* NULL */)
{
	set_defaults();
}

ParatopeSiteConstraintMover::ParatopeSiteConstraintMover(AntibodyInfoCOP ab_info) :
	protocols::moves::Mover(),
	current_func_(/* NULL */)
{
	ab_info_ = ab_info;
	set_defaults();
}

ParatopeSiteConstraintMover::~ParatopeSiteConstraintMover(){}

ParatopeSiteConstraintMover::ParatopeSiteConstraintMover( ParatopeSiteConstraintMover const & src ):
	protocols::moves::Mover( src ),
	cdrs_to_apply_( src.cdrs_to_apply_),
	paratope_residues_( src.paratope_residues_ ),
	antigen_chains_( src.antigen_chains_ ),
	interface_distance_( src.interface_distance_ )

{
	if ( src.ab_info_ ) ab_info_ =  AntibodyInfoOP( new AntibodyInfo( *src.ab_info_));
	if ( src.current_func_ ) current_func_ = current_func_->clone();
}


void
ParatopeSiteConstraintMover::set_defaults() {
	cdrs_to_apply_.clear();
	cdrs_to_apply_.resize(6, true);

	antigen_chains_.clear();
	paratope_residues_.clear();
	interface_distance_ = 10;
	//cst_set_ = new core::scoring::constraints::ConstraintSet();
}

void
ParatopeSiteConstraintMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap &,
	Filters_map const & ,
	moves::Movers_map const & ,
	Pose const & pose
){
	//Paratope Constraint options
	if ( tag->hasOption("paratope_cdrs") ) {
		cdrs_to_apply_ = get_cdr_bool_from_tag(tag, "paratope_cdrs");
	} else {
		cdrs_to_apply_.clear();
		cdrs_to_apply_.resize(6, true);
	}

	interface_distance_ = tag->getOption< core::Real >("interface_dis", interface_distance_);

	if ( tag->hasOption("antigen_chains") ) {
		utility::vector1<std::string> chain_strings = utility::string_split_multi_delim(tag->getOption<std::string>("antigen_chains"), ":,'`~+*&|;.");
		antigen_chains_ = core::pose::get_chain_ids_from_chains(chain_strings, pose);
	}

	if ( tag->hasOption("paratope_residues_pdb") && tag->hasOption("paratope_residues") ) {
		utility_exit_with_message("Cannot specify both paratope_residues_pdb and paratope_residues.");
	}

	//Rosetta Numberings
	if ( tag->hasOption("paratope_residues") ) {
		TR << "Using paratope as user set residues." << std::endl;
		paratope_residues_.clear();
		paratope_residues_.resize(pose.size(), false);
		utility::vector1<std::string> residues = utility::string_split_multi_delim(tag->getOption<std::string>("paratope_residues"), ",'`~+*&|;. ");
		for ( core::Size i = 1; i <= residues.size(); ++i ) {
			paratope_residues_[ utility::string2Size( residues[ i ])] = true;
		}

	}

	//PDB Numbering
	if ( tag->hasOption("paratope_residues_pdb") ) {
		TR << "Using paratope as user set residues." << std::endl;

		paratope_residues_ = design::get_resnums_from_strings_with_ranges(pose, utility::string_split_multi_delim(tag->getOption<std::string>("paratope_residues_pdb"), ",; ") );

	}

}

void
ParatopeSiteConstraintMover::constrain_to_paratope_cdrs(const vector1<CDRNameEnum>& paratope_cdrs){
	cdrs_to_apply_.clear();
	cdrs_to_apply_.resize(6, false);
	for ( core::Size i = 1; i <= paratope_cdrs.size(); ++i ) {
		cdrs_to_apply_[core::Size(paratope_cdrs[i])] = true;
	}

}

void
ParatopeSiteConstraintMover::constrain_to_paratope_cdrs(const vector1<bool>& paratope_cdrs) {
	cdrs_to_apply_ = paratope_cdrs;
}

void
ParatopeSiteConstraintMover::constrain_to_antigen_chains(const vector1<core::Size>& antigen_chains){
	antigen_chains_ = antigen_chains;
}

void
ParatopeSiteConstraintMover::constrain_to_paratope_residues(const vector1<bool>& paratope_residues) {
	paratope_residues_ = paratope_residues;
}

void
ParatopeSiteConstraintMover::set_constraint_func(core::scoring::func::FuncOP constraint_func){
	current_func_ = constraint_func;
}

void
ParatopeSiteConstraintMover::set_interface_distance(core::Real interface_distance){
	interface_distance_ = interface_distance;
}

void
ParatopeSiteConstraintMover::remove(core::pose::Pose& pose, bool reset_paratope_residues){
	using namespace core::scoring::constraints;

	if ( ! reset_paratope_residues ) {
		debug_assert(paratope_residues_.size() == pose.size());
	} else {
		this->setup_paratope_residues_from_cdrs(pose);
	}

	if ( antigen_chains_.size() == 0 ) antigen_chains_ = ab_info_->get_antigen_chain_ids(pose);

	vector1<ConstraintOP> csts_to_be_removed;
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ! paratope_residues_[i] ) continue;

		for ( core::Size x = 1; x <= antigen_chains_.size(); ++x ) {
			SiteConstraintOP res_constraint = setup_constraints(pose, i, utility::to_string(core::pose::get_chain_from_chain_id(x, pose)));
			csts_to_be_removed.push_back(res_constraint);
		}
	}
	pose.remove_constraints(csts_to_be_removed, true);
}

void
ParatopeSiteConstraintMover::setup_paratope_residues_from_cdrs(core::pose::Pose const & pose){

	paratope_residues_.clear();
	paratope_residues_.resize(pose.size(), false);
	for ( core::Size i =1; i <= core::Size(CDRNameEnum_total); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);

		if ( cdrs_to_apply_[cdr] ) {
			for ( core::Size x = ab_info_->get_CDR_start(cdr, pose); x <= ab_info_->get_CDR_end(cdr, pose); ++x ) {
				paratope_residues_[x] = true;
			}
		}
	}
}

void
ParatopeSiteConstraintMover::apply(core::pose::Pose& pose){

	using namespace core::scoring::constraints;

	if ( ! ab_info_ ) {
		ab_info_ = AntibodyInfoCOP( AntibodyInfoOP( new AntibodyInfo(pose) ) );
	}
	//Check if antigen is present
	if ( ! ab_info_->antigen_present() ) {
		TR <<"Antigen not present!  Could not apply constraints" << std::endl;
		set_last_move_status(protocols::moves::FAIL_BAD_INPUT);
		return;
	}
	//If the antibody is camelid, remove trying to setup constraints to light chain
	if ( ab_info_->is_camelid() ) {
		cdrs_to_apply_[l1] = false;
		cdrs_to_apply_[l2] = false;
		cdrs_to_apply_[l3] = false;
	}

	//Check any settings, set defaults from our antibody info.
	if ( antigen_chains_.size() == 0 ) antigen_chains_ = ab_info_->get_antigen_chain_ids(pose);
	if ( paratope_residues_.size() == 0 ) setup_paratope_residues_from_cdrs(pose);

	if ( paratope_residues_.size() != pose.size() ) {
		TR << "Paratope residues does not match total residues. Using settings from paratope CDRs."<<std::endl;
		setup_paratope_residues_from_cdrs(pose);
	}
	//Check any set function
	if ( ! current_func_ ) {
		current_func_ = core::scoring::func::FuncOP( new core::scoring::func::FlatHarmonicFunc(0, 1, interface_distance_) );
	}

	debug_assert(paratope_residues_.size() == pose.size());

	//Ready to go!
	ConstraintCOPs current_csts = pose.constraint_set()->get_all_constraints();

	//pose.constraint_set()->show(TR);
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ! paratope_residues_[i] ) continue;

		for ( core::Size x = 1; x <= antigen_chains_.size(); ++x ) {
			SiteConstraintOP res_constraint = setup_constraints(pose, i, utility::to_string(core::pose::get_chain_from_chain_id(antigen_chains_[x], pose)));

			//Use find  - this may take some time.  Also, the func comparisons seem to be by value - I really don't think they would ever be the same since we clone.
			//Looks like clone of atom pair constraint DOES NOT clone the func, so we should be ok.
			if ( std::find(current_csts.begin(), current_csts.end(), res_constraint) == current_csts.end() ) {
				pose.add_constraint(res_constraint);
				//TR<< "Added constraint: " << i << " to chain "<<utility::to_string(core::pose::get_chain_from_chain_id(antigen_chains_[x], pose)) << std::endl;
			}
		}
	}
}

core::scoring::constraints::SiteConstraintOP
ParatopeSiteConstraintMover::setup_constraints(core::pose::Pose const & pose, core::Size resnum, const std::string chain){
	using namespace core::scoring::constraints;

	//core::scoring::constraints::ConstraintSetOP atom_constraints = new ConstraintSet();
	//AmbiguousConstraintOP res_constraint = new AmbiguousConstraint();


	SiteConstraintOP atom_constraint( new SiteConstraint() );
	atom_constraint->setup_csts(
		resnum,
		"CA",
		chain,
		pose,
		current_func_);

	return atom_constraint;

}

std::string ParatopeSiteConstraintMover::get_name() const {
	return mover_name();
}

std::string ParatopeSiteConstraintMover::mover_name() {
	return "ParatopeSiteConstraintMover";
}

void ParatopeSiteConstraintMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attributes_for_get_cdr_bool_from_tag(
		attlist, "paratope_cdrs",
		"Specifically set the paratope as these CDR.");

	attlist + XMLSchemaAttribute::attribute_w_default(
		"interface_distance", xsct_real,
		"Distance in Angstroms for the interface, which effects when the SiteConstraint penalty begins.",
		"10.0");

	attlist + XMLSchemaAttribute(
		"antigen_chains", xs_string,
		"Specify the particular antigen to create the SiteConstraint to");

	attlist + XMLSchemaAttribute(
		"paratope_residues_pdb", xs_string,
		"Set specific residues as the paratope. (Ex: 24L,26L-32L, 44H-44H:A). "
		"Can specify ranges or individual residues as well as insertion codes "
		"(Ex: 44H:A with A being insertion code).");

	attlist + XMLSchemaAttribute(
		"paratope_residues", xs_string,
		"Set paratope_residues instead of paratope_residues_pdb as the internal "
		"rosetta residue numbers (Ex: 14,25,26). Internal rosetta numbering "
		"parsing does not currently support ranges.");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Adds SiteConstraints from the Antibody Paratope to the antigen, "
		"defined for simplicity as the CDRs. Individual residues of the paratope can be set, "
		"or specific CDRs of the paratope can be set as well. These help to keep only the "
		"paratope in contact with the antigen (as apposed to the framework) during rigid-body "
		"movement. See the Constraint File Overview for more information on manually adding "
		"SiteConstraints. Do not forget to add the atom_pair_constraint term to your scorefunction. "
		"A weight of .01 for the SiteConstraints seems optimum. "
		"Default paratope is defined as all 6 CDRs (or 3 if working with a camelid antibody).",
		attlist );
}

std::string ParatopeSiteConstraintMoverCreator::keyname() const {
	return ParatopeSiteConstraintMover::mover_name();
}

protocols::moves::MoverOP
ParatopeSiteConstraintMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ParatopeSiteConstraintMover );
}

void ParatopeSiteConstraintMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ParatopeSiteConstraintMover::provide_xml_schema( xsd );
}


} //constraints
} //antibody
} //protocols
