// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


/// @file protocols/protein_interface_design/movers/SetupHotspotConstraintsLoopsMover.cc
/// @brief
/// @author Jacob Corn (jecorn@u.washington.edu), Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/movers/SetupHotspotConstraintsLoopsMover.hh>
#include <protocols/protein_interface_design/movers/SetupHotspotConstraintsLoopsMoverCreator.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <complex>

#include <protocols/moves/MoverStatus.hh>
#include <core/types.hh>

#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/util.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/Loops.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <core/kinematics/FoldTree.hh>


#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/protein_interface_design/movers/SetupHotspotConstraintsMover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/BackboneStubConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <core/fragment/FragSet.hh>
#include <utility/exit.hh>

#include <core/scoring/electron_density/util.hh>

#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/util.hh>
#include <protocols/moves/MoverContainer.hh>
// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

//silly using/typedef
#include <core/id/AtomID_Map.hh>

#include <core/chemical/AA.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/pack_rotamers.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <protocols/comparative_modeling/LoopRelaxMover.hh>

#include <core/import_pose/import_pose.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
// Utility headers
#include <basic/options/option_macros.hh>

#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;

static basic::Tracer tr( "protocols.protein_interface_design.movers.SetupHotspotConstraintsLoopsMover" );

// XRW TEMP std::string
// XRW TEMP SetupHotspotConstraintsLoopsMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SetupHotspotConstraintsLoopsMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SetupHotspotConstraintsLoopsMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SetupHotspotConstraintsLoopsMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SetupHotspotConstraintsLoopsMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "SetupHotspotConstraintsLoops";
// XRW TEMP }

SetupHotspotConstraintsLoopsMover::SetupHotspotConstraintsLoopsMover() :
	protocols::moves::Mover( "SetupHotspotConstraintsLoopsMover" ),
	chain_to_design_( 2 ),
	CB_force_constant_( 0.5 ),
	worst_allowed_stub_bonus_( 0.0 ),
	apply_self_energies_( true ),
	bump_cutoff_( 4.0 ),
	apply_ambiguous_constraints_( true ),
	colonyE_( false ),
	resfile_( "NONE" ),
	loop_start_( 0 ),
	loop_stop_( 0 )
{}

protocols::moves::MoverOP
SetupHotspotConstraintsLoopsMover::clone() const{
	return protocols::moves::MoverOP( new SetupHotspotConstraintsLoopsMover( *this ) );
}

protocols::moves::MoverOP
SetupHotspotConstraintsLoopsMover::fresh_instance() const{
	return protocols::moves::MoverOP( new SetupHotspotConstraintsLoopsMover() );
}

SetupHotspotConstraintsLoopsMover::SetupHotspotConstraintsLoopsMover(
	protocols::hotspot_hashing::HotspotStubSetCOP hotspot_stub_set
) :
	// core::Size const chain_to_design,
	// core::Real const & CB_force_constant( 0.5 ),
	//core::Real const & worst_allowed_stub_bonus,
	//bool const apply_self_energies,
	//core::Real const & bump_cutoff,
	// bool const apply_ambiguous_constraints( true ),
	//bool const colonyE
	protocols::moves::Mover( "SetupHotspotConstraintMover" ),
	chain_to_design_( 2 ),
	CB_force_constant_( 0.5 ),
	worst_allowed_stub_bonus_( 0.0 ),
	apply_self_energies_( true ),
	bump_cutoff_( 4.0 ),
	apply_ambiguous_constraints_( true ),
	colonyE_( false ),
	resfile_( "NONE" ),
	loop_start_( 0 ),
	loop_stop_( 0 )

	//  chain_to_design_(chain_to_design),
	//  CB_force_constant_(CB_force_constant),
	//  worst_allowed_stub_bonus_(worst_allowed_stub_bonus),
	//  apply_self_energies_(apply_self_energies),
	//  bump_cutoff_(bump_cutoff),
	//  apply_ambiguous_constraints_(apply_ambiguous_constraints),
	//  colonyE_( colonyE)
{
	//  packer_task_ = packer_task->clone();
	hotspot_stub_set_ = protocols::hotspot_hashing::HotspotStubSetOP( new protocols::hotspot_hashing::HotspotStubSet( *hotspot_stub_set ) );
}

SetupHotspotConstraintsLoopsMover::SetupHotspotConstraintsLoopsMover( SetupHotspotConstraintsLoopsMover const & init ) :
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( init ),
	chain_to_design_( init.chain_to_design_),
	CB_force_constant_(init.CB_force_constant_),
	worst_allowed_stub_bonus_(init.worst_allowed_stub_bonus_),
	apply_self_energies_(init.apply_self_energies_),
	bump_cutoff_(init.bump_cutoff_),
	apply_ambiguous_constraints_(init.apply_ambiguous_constraints_),
	colonyE_( init.colonyE_ ),
	resfile_( init.resfile_ ),
	loop_start_( init.loop_start_ ),
	loop_stop_( init.loop_stop_ )
{
	hotspot_stub_set_ = protocols::hotspot_hashing::HotspotStubSetOP( new protocols::hotspot_hashing::HotspotStubSet( *init.hotspot_stub_set_ ) );
}

core::Size
SetupHotspotConstraintsLoopsMover::generate_csts(
	core::pose::Pose const& pose,
	core::scoring::constraints::ConstraintCOPs& constraints
) {
	if ( colonyE_ ) {
		protocols::hotspot_hashing::HotspotStubSetOP colonyE_set = hotspot_stub_set_->colonyE();
		hotspot_stub_set_ = colonyE_set;
	}
	core::id::AtomID fixed_atom(1, 1);
	core::Real worst_allowed_stub_bonus = 0;
	// bool apply_self_energies = false; // Unused variable causes warning.

	// Take in a scaffold pose (with PackerTask, for its DesignMap), and a set of stubs.
	// Each repacked residue will get one "AmbiguousConstraint".
	// This AmbiguousConstraint will contain a series of BackboneStubConstraints (one for each valid stub)
	runtime_assert( CB_force_constant_ > -1E-6 ); // these can't be negative
	runtime_assert( worst_allowed_stub_bonus < 1E-6 ); // these can't be positive

	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
	if ( resfile_ =="NONE" && basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		core::pack::task::parse_resfile(pose, *task);
	} else if ( resfile_ != "NONE" ) {
		core::pack::task::parse_resfile(pose, *task, resfile_ );
	}
	tr.Debug << "finished reading resfile " << std::endl;
	task->show( tr.Trace );
	// *****associate the stub set with the unbound pose
	// *****associate the stub set with the unbound pose
	// *****hotspot_stub_set_->pair_with_scaffold( pose, partner );

	protocols::filters::FilterCOP true_filter( protocols::filters::FilterOP( new protocols::filters::TrueFilter ) );
	for ( core::Size resnum=loop_start_; resnum <= loop_stop_; ++resnum ) {
		//  if ( task->pack_residue(resnum) ) {
		hotspot_stub_set_->pair_with_scaffold( pose, pose.chain( resnum ), true_filter );
		break;
		//}
	}

	tr.Info << "Making hotspot constraints..." << std::endl;
	//Size scaffold_seqpos(0);
	Size ct_cst( 0 );
	for ( core::Size resnum=loop_start_; resnum <= loop_stop_; ++resnum ) {
		tr.Debug << "make constraints for position " << resnum << " ..." << std::endl;
		// Check that this position is allowed to be used for stub constraints
		//  if ( ! task->pack_residue(resnum) ) continue;

		// sets the index used by the hotspot for its associated scaffold
		//scaffold_seqpos = resnum - pose.conformation().chain_begin( pose.chain( resnum ) );  // set but never used

		// Start the vector which will become a single AmbiguousConstraint, if apply_ambiguous_constraints is true
		utility::vector1< core::scoring::constraints::ConstraintCOP > ambig_csts;
		// Loop over all allowed AAs at this position
		std::list< core::chemical::ResidueTypeCOP > allowed_aas = task->residue_task( resnum ).allowed_residue_types();
		for ( std::list< core::chemical::ResidueTypeCOP >::const_iterator restype = allowed_aas.begin();
				restype != allowed_aas.end(); ++restype ) {
			tr.Debug << "looking for suitable hotspots for " << (*restype)->name3() << " at pos " << resnum << std::endl;
			// Loop over all stubs with this restype
			hotspot_hashing::HotspotStubSet::Hotspots res_stub_set( hotspot_stub_set_->retrieve( (*restype )->name3() ) );
			tr.Debug << "found " << res_stub_set.size() << " suitable hotspots" << std::endl;
			for ( auto & hs_stub : res_stub_set ) {

				// prevent Gly/Pro constraints
				if ( (hs_stub.second->residue()->aa() == core::chemical::aa_gly) || (hs_stub.second->residue()->aa() == core::chemical::aa_pro && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline] ) ) {
					tr.Info << "WARNING - Gly/Pro stubs cannot be used for constraints." << std::endl;
					continue;
				}

				// prevent Gly/Pro constraints
				if ( (pose.residue(resnum).aa() == core::chemical::aa_gly) || (pose.residue(resnum).aa() == core::chemical::aa_pro && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline]) ) {
					tr.Debug << "WARNING - Position " << resnum << " is currently Gly/Pro and cannot be used for stub constraints." << std::endl;
					continue;
				}

				core::Real stub_bonus_value = hs_stub.second->bonus_value();
				tr.Trace << "stub_bonus: " << stub_bonus_value << std::endl;
				if ( stub_bonus_value < worst_allowed_stub_bonus ) {
					hs_stub.second->set_scaffold_status( resnum, protocols::hotspot_hashing::accept );
					//tr.Info << " SuccSelfEnergy=" << stub_bonus_value << std::endl;
					// ****** accept the pairing -- do we really want this? better to just reject, since bb fit doesn't necessarily mean good pair
					// ****** hs_stub->scaffold_status( resnum, accept );

					// Build a BackboneStubConstraint from this stub
					if ( apply_ambiguous_constraints_ ) {
						// Push it onto ambig_csts for this residue
						ct_cst++;
						ambig_csts.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::BackboneStubConstraint( pose, resnum, fixed_atom, *(hs_stub.second->residue()), stub_bonus_value, CB_force_constant_ ) ) );
					} else {
						ct_cst++;
						// Apply it directly
						constraints.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::BackboneStubConstraint( pose, resnum, fixed_atom, *(hs_stub.second->residue()), stub_bonus_value, CB_force_constant_ ) ) );
					}
				} else hs_stub.second->set_scaffold_status( resnum, protocols::hotspot_hashing::reject );
				//else tr.Info << " FailSelfEnergy=" << stub_bonus_value << std::endl;
				// ****** reject the pairing
				// ******else hs_stub->scaffold_status( resnum, reject );
			}
		}

		// Finally, add the constraint corresponding to this resnum to the main set
		if ( ( apply_ambiguous_constraints_ ) && ( ambig_csts.size() > 0 ) ) {
			constraints.push_back( core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AmbiguousConstraint(ambig_csts) ) );
		}
	}
	return ct_cst;
}


void
SetupHotspotConstraintsLoopsMover::apply( core::pose::Pose & pose ) {
	if ( std::abs(CB_force_constant_) > 1E-9 ) {
		core::scoring::constraints::ConstraintCOPs constraints;
		Size ct_cst = generate_csts( pose, constraints );
		constraints = pose.add_constraints( constraints );
		tr.Info << "Applied " << ct_cst << " hotspots in " << constraints.size() << " constraints to the pose." << std::endl;
	} else { //CB_force is <1e-19
		core::scoring::constraints::ConstraintSetOP empty_constraint_set( new core::scoring::constraints::ConstraintSet );
		pose.constraint_set( empty_constraint_set );
	}
}

// XRW TEMP std::string
// XRW TEMP SetupHotspotConstraintsLoopsMover::get_name() const {
// XRW TEMP  return "SetupHotspotConstraintsLoopsMover";
// XRW TEMP }

/// This needs to be parsed before all other movers b/c it changes scorefxns
void
SetupHotspotConstraintsLoopsMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	using core::Real;
	chain_to_design_ = tag->getOption<Size>( "redesign_chain", 2 );
	resfile_ = tag->getOption<std::string>("resfile","NONE");

	CB_force_constant_ = tag->getOption<Real>( "cb_force", 0.5 );

	worst_allowed_stub_bonus_ = tag->getOption<Real>( "worst_allowed_stub_bonus", 0 );
	apply_self_energies_ = tag->getOption<bool>( "apply_stub_self_energies", false );
	bump_cutoff_ = tag->getOption<Real>( "apply_stub_bump_cutoff", 10. );
	apply_ambiguous_constraints_ = tag->getOption<bool>( "pick_best_energy_constraint", true );
	// core::Real const bb_stub_cst_weight( tag->getOption< core::Real >( "backbone_stub_constraint_weight", 1.0 ) );  // Unused variable causes warning.
	loop_start_ = tag->getOption<Size>("start", 0 );
	loop_stop_ = tag->getOption<Size>("stop", 0 );
	if ( loop_start_ <= 0 ) {
		utility_exit_with_message( "please provide loop-start with 'start' option" );
	}
	if ( loop_stop_ <= loop_start_ ) {
		utility_exit_with_message( "loop-stop has to be larger than loop-start. assign with 'stop'" );
	}
	tr.Debug << "setup hotspots for loop " << loop_start_ << " to " << loop_stop_ << std::endl;
	colonyE_ = tag->getOption<bool>( "colonyE", false );

	hotspot_stub_set_ = protocols::hotspot_hashing::HotspotStubSetOP( new hotspot_hashing::HotspotStubSet );
	if ( tag->hasOption( "stubfile" ) ) {
		std::string const hotspot_fname( tag->getOption<std::string>( "stubfile", "stubs.pdb" ) );
		hotspot_stub_set_->read_data( hotspot_fname );
	}
	utility::vector1< TagCOP > const branch_tags( tag->getTags() );
	for ( TagCOP const curr_tag : branch_tags ) {
		if ( curr_tag->getName() == "HotspotFiles" ) {
			utility::vector1< TagCOP > const branch_tags2( curr_tag->getTags() );
			for ( TagCOP const curr_tag2 : branch_tags2 ) {
				std::string const file_name( curr_tag2->getOption< std::string >( "file_name" ) );
				std::string const nickname( curr_tag2->getOption< std::string >( "nickname" ) );
				auto const stub_num( curr_tag2->getOption< core::Size >( "stub_num", 100000 ) );
				hotspot_hashing::HotspotStubSetOP temp_stubset( new hotspot_hashing::HotspotStubSet );
				temp_stubset->read_data( file_name );
				temp_stubset->remove_random_stubs_from_set( temp_stubset->size() - stub_num );
				hotspot_stub_set_->add_stub_set( *temp_stubset );
				tr.Info<<"Read stubset from file "<<file_name<<" and associating it with name "<<nickname<<'\n';
				tr.Info<<stub_num<<" stubs kept in memory\n";
				data.add( "hotspot_library", nickname, temp_stubset );
			}
		} else {
			utility_exit_with_message( curr_tag->getName() + " not recognized by SetupHotspotConstraints, did you mean HotspotFiles?" );
		}
	}

	tr.Info<<"applying hotspot hashing constraints to pose with " << " cb_force weight of "<<CB_force_constant_<<", apply ambiguous constraints set to "<<apply_ambiguous_constraints_<< " and colonyE set to " << colonyE_ << "\n";
	data.add( "constraints" , "hotspot_stubset", hotspot_stub_set_ );

	tr.Info.flush();
}

SetupHotspotConstraintsLoopsMover::~SetupHotspotConstraintsLoopsMover() = default;

std::string SetupHotspotConstraintsLoopsMover::get_name() const {
	return mover_name();
}

std::string SetupHotspotConstraintsLoopsMover::mover_name() {
	return "SetupHotspotConstraintsLoops";
}

std::string subtag_for_SetupHotspotConstraintsLoopsMover( std::string const & foo ) {
	return "subtag_for_SHCLM_" + foo + "_type";
}

void SetupHotspotConstraintsLoopsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "redesign_chain", xsct_non_negative_integer, "The chain on which to redesign, numbered sequentially from 1", "2" )
		+ XMLSchemaAttribute::attribute_w_default( "resfile", xs_string, "Resfile applied", "NONE" )
		+ XMLSchemaAttribute::attribute_w_default( "cb_force", xsct_real, "Force applied to CB atoms in hotspot constraints", "0.5" )
		+ XMLSchemaAttribute::attribute_w_default( "worst_allowed_stub_bonus", xsct_real, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "apply_stub_self_energies", xsct_rosetta_bool, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "apply_stub_bump_cutoff", xsct_real, "XRW TO DO", "10" )
		+ XMLSchemaAttribute::attribute_w_default( "pick_best_energy_constraint", xsct_rosetta_bool, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::required_attribute( "start", xsct_positive_integer, "First loop residue" )
		+ XMLSchemaAttribute::required_attribute( "stop", xsct_positive_integer, "Final loop residue" )
		+ XMLSchemaAttribute::attribute_w_default( "colonyE", xsct_rosetta_bool, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "stubfile", xs_string, "File containing hotspot stubs", "stubs.pdb" );

	AttributeList subtag_attributes, subsubtag_attributes;
	subsubtag_attributes + XMLSchemaAttribute::required_attribute( "file_name", xs_string, "File name of the hot spot stub file" )
		+ XMLSchemaAttribute::required_attribute( "nickname", xs_string, "Nickname for the hot spot stub file" )
		+ XMLSchemaAttribute::attribute_w_default( "stub_num", xsct_positive_integer, "XRW TO DO", "100000" );


	utility::tag::XMLSchemaSimpleSubelementList subsubelements;
	subsubelements.add_simple_subelement( "HotspotFile"/*name never read from branch_tags2*/, subsubtag_attributes, "Tags for each file"/*, 0minoccurs*/ );

	XMLSchemaComplexTypeGenerator subtag_gen;
	subtag_gen.complex_type_naming_func( & subtag_for_SetupHotspotConstraintsLoopsMover )
		.element_name( "HotspotFiles" )
		.description( "Wrapper for tags describing each of several hotspot files to read" )
		.add_attributes( subtag_attributes )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subsubelements )
		.write_complex_type_to_schema( xsd );

	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_already_defined_subelement( "HotspotFiles", & subtag_for_SetupHotspotConstraintsLoopsMover/*, 0*/ );
	//.complex_type_naming_func( & subtag_for_SetupHotspotConstraintsLoopsMover );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, ssl );
}

std::string SetupHotspotConstraintsLoopsMoverCreator::keyname() const {
	return SetupHotspotConstraintsLoopsMover::mover_name();
}

protocols::moves::MoverOP
SetupHotspotConstraintsLoopsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetupHotspotConstraintsLoopsMover );
}

void SetupHotspotConstraintsLoopsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetupHotspotConstraintsLoopsMover::provide_xml_schema( xsd );
}


} //movers
} //protein_interface_design
} //protocols
