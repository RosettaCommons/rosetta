// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/glycopeptide_docking/GlycopeptideDockingProtocol.cc
/// @brief A protocol to predict glycopeptide_docking and elongation of glycans on peptide/glycopeptide, lipid and protein substrates by glycosyltransferase. Current implementation only includes glycosylation of peptides and glycopeptides. It has been tested on GalNAcT2-peptide and GalNAcT12-glypeptide systems to predict peptide substrate specificity(T2) and di-glycopeptide binding sites(T12). Also, limited testing on N-linked glycosylation by enzyme ApNGT.
/// @author Yashes Srinivasan (yashess@gmail.com), Sai Pooja Mahajan (saipooja@gmail.com)
///
// Unit headers
#include <protocols/glycopeptide_docking/GlycopeptideDockingProtocol.hh>
#include <protocols/glycopeptide_docking/GlycopeptideDockingProtocolCreator.hh>
#include <protocols/glycopeptide_docking/GlycopeptideDockingFlags.hh>
#include <protocols/glycopeptide_docking/GlycopeptideDockingHighResRefinement.hh>
#include <protocols/glycopeptide_docking/GlycopeptideDockingLowResRefinement.hh>
#include <protocols/glycopeptide_docking/utils.hh>
#include <protocols/carbohydrates/SimpleGlycosylateMover.hh>
#include <protocols/carbohydrates/GlycanSampler.hh>

// Protocol headers
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/constraint_movers/ConstraintSetMover.hh>
#include <protocols/docking/util.hh>
#include <protocols/docking/metrics.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh> //needed for rmsd_no_super
#include <core/scoring/ScoreType.hh>
#include <core/import_pose/import_pose.hh> //pose_from_file
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/pose/extra_pose_info_util.hh>
// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <string>
#include <list>
#include <map>

static basic::Tracer TR( "protocols.glycopeptide_docking.GlycopeptideDockingProtocol" );

namespace protocols {
namespace glycopeptide_docking {

using namespace std;
using namespace utility;
using namespace core;
using namespace protocols;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
GlycopeptideDockingProtocol::GlycopeptideDockingProtocol():
	protocols::moves::Mover( GlycopeptideDockingProtocol::mover_name() )
{
	using namespace glycopeptide_docking;
	using namespace constraint_movers;
	using namespace scoring;
	using namespace pose;

	type( "GlycopeptideDockingProtocol" );

	flags_ = GlycopeptideDockingFlagsOP( utility::pointer::make_shared< GlycopeptideDockingFlags>() );
	if ( flags_->get_constraints_file().size() > 0 ) constraint_setter_ = ConstraintSetMoverOP( utility::pointer::make_shared< ConstraintSetMover >());
	std::cout<<"Preparing score function"<<std::endl;
	prepare_score_function();
	std::cout<<"score function done"<<std::endl;
	//pymol_mover_ = protocols::moves::PyMOLMoverOP(utility::pointer::make_shared< protocols::moves::PyMOLMover>());
	//pymol_mover_->keep_history(true);

}


/// @brief Constructor with flags object as argument
GlycopeptideDockingProtocol::GlycopeptideDockingProtocol(GlycopeptideDockingFlags &flags):
	protocols::moves::Mover( GlycopeptideDockingProtocol::mover_name() )
{
	using namespace glycopeptide_docking;
	using namespace constraint_movers;
	using namespace scoring;
	using namespace pose;

	type( "GlycopeptideDockingProtocol" );

	flags_ = GlycopeptideDockingFlagsOP( utility::pointer::make_shared< GlycopeptideDockingFlags>(flags) );
	if ( flags_->get_constraints_file().size() > 0 ) constraint_setter_ = ConstraintSetMoverOP( utility::pointer::make_shared< ConstraintSetMover >());

	prepare_score_function();
	//pymol_mover_ = protocols::moves::PyMOLMoverOP(utility::pointer::make_shared< protocols::moves::PyMOLMover>());
	//pymol_mover_->keep_history(true);

}


////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
GlycopeptideDockingProtocol::GlycopeptideDockingProtocol( GlycopeptideDockingProtocol const & object_to_copy ):
	protocols::moves::Mover( object_to_copy )
{
	copy_data( *this, object_to_copy );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Assignment operator
GlycopeptideDockingProtocol &
GlycopeptideDockingProtocol::operator=( GlycopeptideDockingProtocol const & object_to_copy )
{
	// Abort self-assignment.
	if ( this != &object_to_copy ) {
		moves::Mover::operator=( object_to_copy );
		copy_data( *this, object_to_copy );
	}
	return *this;
}


////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
GlycopeptideDockingProtocol::~GlycopeptideDockingProtocol(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @details Apply 1. Setups sampling options that require pose, including setting up
/// glycosylation fold tree, substrate and enzyme chains. If user specifies residues to be
/// preglycosylated before sampling, apply does it. 2. Randomizes peptide backbone torsions if requested.
/// 3. Runs low-resolution sampling of peptide/glycopeptide-enzyme if enabled. 4. Runs high-resolution
/// of peptide/glycopeptide-enzyme of enabled. 5. Records metrics.
void
GlycopeptideDockingProtocol::apply( core::pose::Pose& pose )
{
	using namespace scoring;
	using namespace glycopeptide_docking;

	/* TODO: include <core/io/raw_data/ScoreMap.hh>
	// Use to add intermediate data to the scorefile */
	if ( !ref_pose_ ) ref_pose_ = pose::PoseOP( utility::pointer::make_shared< pose::Pose>( pose ) );

	//Debug information
	TR.Debug << "Seed " << numeric::random::rg().get_seed() <<std::endl;
	// TR.Debug << pose.residue( flags_->get_sugar_donor() ).name() << endl;

	TR<<"setting up substrate chain"<<std::endl;
	// Update chains with pose information
	flags_->set_substrate_chain( pose );
	flags_->set_enzyme_chain( pose );
	flags_->show();

	/* TODO:
	The foldtree needs to be setup correctly for this code to
	work - fix foldtree setup after glycosylation before enabling this block
	utility::vector1<core::Size> const positions =  flags_->get_residues_to_preglycosylate();
	utility::vector1<std::string> sugar_names =flags_->get_sugars_to_preglycosylate();
	if ( !positions.empty() ) {
	if (positions.size() == sugar_names.size()){
	glycosylate_residues(pose,positions,sugar_names);
	* Modifying the ref pose here to match for rms util
	* Makes sense since the input pose does not have sugars
	* or generate a residue selector to pass to rms util without
	* new sugars TODO
	glycosylate_residues(*ref_pose_,positions,sugar_names);
	} else {
	TR << "number of elements in glycosylation positions and sugar names does not match. Will not preglycosylate." << endl;
	TR << "Positions: " << positions << endl;
	TR << "Names: " << sugar_names << endl;
	}
	}*/

	// Setup foldtree
	//ft_substrate is currently unused but could be used to equilibrate the glycopeptide
	TR << endl << "Setting up fold tree..." << endl;
	ft_substrate_ = core::kinematics::FoldTreeOP( utility::pointer::make_shared< core::kinematics::FoldTree>( pose.fold_tree()) );
	setup_glycosylation_foldtree( pose , flags_, ft_docking_);
	TR << endl << " Fold tree: " << endl;
	TR << " " << pose.fold_tree() << endl;

	// prepare_score_function();
	// Should be ok now - the scorefunction needs be cleaned up if exiting early
	//from low-res or high-res sampling - this can cause the sf_ to have modified temperatures
	// from simulated annealing or modified weight from ramping of attractive and repulsive
	// potentials since the constructor that sets the score function is only called once.

	// Set up constraints
	if ( constraint_setter_ ) {
		TR << endl << " Setting up constraints..." << endl;
		constraint_setter_->constraint_file( flags_->get_constraints_file() );
		constraint_setter_->apply( pose );
	}


	/* TODO: Add mover to place objects given a set of constraints */

	// Print information about the starting structure ( for debugging purpsoses )

	TR << endl << " Starting score: " << endl;
	sf_->show( std::cout, pose );
	jd2::JobOP job( jd2::JobDistributor::get_instance()->current_job() );

	/* If score_only is enabled, sampling is skipped.
	// This allows generation of score files with
	// glycosylation specific terms. Great for rescoring with
	application.*/

	if ( !flags_->score_only() ) {
		if ( flags_->randomize_substrate_torsions() ) {
			TR << endl
				<< " Randomizing substrate conformation..." << endl;
			core::uint residue_start(flags_->first_residue_substrate());
			core::uint residue_stop(flags_->last_residue_substrate());

			for ( core::uint residue(residue_start); residue <= residue_stop; residue++ ) {
				if ( residue != flags_->anchor_residue_substrate() ) {
					pose.set_psi(residue, numeric::random::rg().uniform() * 360);
					pose.set_phi(residue, numeric::random::rg().uniform() * 360);
				}
			}
		}
		/* if ( flags_->randomize_residue()) {
		TR << endl << " Randomizing residue conformation..." << endl;
		pose.set_psi( flags_->get_randomize_residue(), numeric::random::rg().uniform() * 360 );
		pose.set_phi( flags_->get_randomize_residue(), numeric::random::rg().uniform() * 360 );
		}*/

		protocols::carbohydrates::GlycanSamplerOP sugar_sampler(utility::pointer::make_shared<protocols::carbohydrates::GlycanSampler>()); //Most cases currently are mono and di saccharides so N=20 should be fine
		// TODO: Sugar cycles can be a flag or the number can be empirically tested.
		// Using small number of cycles. Most cases currently are mono and di saccharides.
		sugar_sampler->set_rounds(20);
		sugar_sampler->set_scorefunction(sf_);

		if ( flags_->get_allow_glycan_torsion_moves() ) {
			sugar_sampler->apply(pose);
		}

		/////  LOW-RES REFINEMENT  /////
		if ( flags_->low_res_refinement() ) {
			TR << endl
				<< " Low-res refinement..." << endl;
			GlycopeptideDockingLowResRefinement lowres_protocol(flags_, sf_, ft_docking_, ft_substrate_);
			lowres_protocol.apply(pose);
			if ( flags_->get_allow_glycan_torsion_moves() ) {
				sugar_sampler->apply(pose);
			}

			if ( flags_->debug_pdbs() ) {
				std::string pdbname("final_lowres");
				write_debug_pdb(pose, job->nstruct_max(), job->nstruct_index(), pdbname);
			}

		}

		///// HIGH_RES REFINEMENT /////

		if ( flags_->high_res_refinement() ) {
			TR << endl << " High-res refinement..." << endl;
			GlycopeptideDockingHighResRefinement highres_protocol( flags_, sf_ , ft_docking_, ft_substrate_);
			highres_protocol.show(std::cout);
			highres_protocol.apply( pose );
			//This should be included within the highres protocol eventually
			//Needs testing. Experimental feature.
			if ( flags_->get_allow_glycan_torsion_moves() ) sugar_sampler->apply(pose);
		}
	}

	///// POSE METRICS /////

	// Print information about the final structure
	if ( flags_->debug_pdbs() ) {
		std::string pdbname("final");
		write_debug_pdb(pose, job->nstruct_max(), job->nstruct_index(),pdbname);
	}
	vector1<int> movable_jumps(1, flags_->jump_num_substrate());
	record_pose_metrics(pose, flags_, movable_jumps, ref_pose_);
}

void
GlycopeptideDockingProtocol::set_ref_pose_from_filename( std::string const & filename )
{
	ref_pose_ = import_pose::pose_from_file( filename, core::import_pose::PDB_file );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
GlycopeptideDockingProtocol::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
GlycopeptideDockingProtocol::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap&
)
{
	//overwrite options if set at command line
	//TODO: need to make other options available in rosetta scripts
	//untested
	flags_->set_glycosylation_residue(tag->getOption<core::Size>( "donor_residue" )) ;
	flags_->set_sugar_donor(tag->getOption<core::Size>( "sugar_residue")) ;
	flags_->set_glycosylation_refinement(tag->getOption<bool>("lowres",true),tag->getOption<bool>("high_res",true));
	flags_->set_backbone_moves_substrate(tag->getOption<bool>("backbone_moves_highres",true));
	flags_->set_score_only(tag->getOption<bool>("score_only",false));
}

void GlycopeptideDockingProtocol::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the provided AttributeList.
	attlist
		+ XMLSchemaAttribute::required_attribute("glycosylation_residue",xsct_residue_number,"Indicate the site of glycosylation in internal numbering. Alternatively, this is indicate a glycosylated residue. The peptide foldtree emnates outwards from this site.")
		+ XMLSchemaAttribute::required_attribute("donor_residue",xsct_residue_number,"Indicate the sugar donor in internal numbering.")
		//TODO: add options to indicate enzyme chain and peptide chain
		+ XMLSchemaAttribute::attribute_w_default("highres",xsct_rosetta_bool ,"Enable/disable high resolution sampling.","true")
		+ XMLSchemaAttribute::attribute_w_default("lowres",xsct_rosetta_bool,"Enable/disable low resolution sampling.","false")
		+ XMLSchemaAttribute::attribute_w_default("backbone_moves_highres",xsct_rosetta_bool,"Enable/disable backbone moves in high resolution.","true")
		+ XMLSchemaAttribute::attribute_w_default("score_only",xsct_rosetta_bool,"Do not run protocol. Just setup foldtree and output scores.","false");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This is a protocol to predicts the glycosyaltion peptides and glycosylated peptides by glycosyltransferases. The protocol evaluated the peptide-enzyme complex in context of the glycosylation.", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycopeptideDockingProtocol::fresh_instance() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< GlycopeptideDockingProtocol>() );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
GlycopeptideDockingProtocol::clone() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< GlycopeptideDockingProtocol>( *this ) );
}

std::string GlycopeptideDockingProtocol::get_name() const {
	return mover_name();
}

std::string GlycopeptideDockingProtocol::mover_name() {
	return "GlycopeptideDockingProtocol";
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
GlycopeptideDockingProtocolCreator::create_mover() const
{
	return protocols::moves::MoverOP( utility::pointer::make_shared< GlycopeptideDockingProtocol>() );
}

std::string
GlycopeptideDockingProtocolCreator::keyname() const
{
	return GlycopeptideDockingProtocol::mover_name();
}

void GlycopeptideDockingProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GlycopeptideDockingProtocol::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

/// @brief Copy private data, called by copy constructor
void
GlycopeptideDockingProtocol::copy_data( GlycopeptideDockingProtocol & object_to_copy_to, GlycopeptideDockingProtocol const & object_to_copy_from )
{
	object_to_copy_to.flags_ = object_to_copy_from.flags_;
	object_to_copy_to.sf_ = object_to_copy_from.sf_;
	object_to_copy_to.mc_ = object_to_copy_from.mc_;
	object_to_copy_to.constraint_setter_ = object_to_copy_from.constraint_setter_;
	object_to_copy_to.ref_pose_ = object_to_copy_from.ref_pose_;
}

/// @brief Prepares the score function used by the refinement and abinitio protocols
void
GlycopeptideDockingProtocol::prepare_score_function()
{
	///  TODO: test with beta scorefunction and allow
	/// use of specified scorefunction.
	using namespace scoring;

	core::Real Hbond_mult = 1.0;
	core::Real elec_mult = 1.0;
	core::Real sol_mult = 1.0;

	sf_ = get_score_function();

	sf_->set_weight( sugar_bb, 1.0 );
	sf_->set_weight( hbond_sr_bb, sf_->get_weight( hbond_sr_bb ) * Hbond_mult );
	sf_->set_weight( hbond_lr_bb, sf_->get_weight( hbond_lr_bb ) * Hbond_mult );
	sf_->set_weight( hbond_bb_sc, sf_->get_weight( hbond_bb_sc ) * Hbond_mult );
	sf_->set_weight( hbond_sc, sf_->get_weight( hbond_sc ) * Hbond_mult );
	sf_->set_weight( fa_elec, sf_->get_weight( fa_elec ) * elec_mult );
	sf_->set_weight( fa_sol, sf_->get_weight( fa_sol ) * sol_mult );

	// default - when constraints are provided to keep ligands/peptide/motifs in place
	sf_->set_weight( atom_pair_constraint, 10.0 );  /* TODO: Test out weights? Probably make this an option */
}



std::ostream &
operator<<( std::ostream & os, GlycopeptideDockingProtocol const & mover )
{
	mover.show(os);
	return os;
}

} //glycopeptide_docking
} //protocols

