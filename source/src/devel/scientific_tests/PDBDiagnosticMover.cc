// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/scientific_tests/PDBDiagnosticMover.cc
/// @brief does some very lightweight modeling on a PDB.  Meant to be run against the whole PDB as a scientific test
/// @author Steven Lewis (smlewi@gmail.com)

// Unit headers
#include <devel/scientific_tests/PDBDiagnosticMover.hh>
#include <devel/scientific_tests/PDBDiagnosticMoverCreator.hh>

// Core headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>

#include <core/chemical/ResidueType.hh>

#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/MinMover.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/PDB_diagnostic.OptionKeys.gen.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/options/option.hh>
//#include <utility/excn/Exceptions.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "devel.scientific_tests.PDBDiagnosticMover" );

namespace devel {
namespace scientific_tests {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PDBDiagnosticMover::PDBDiagnosticMover():
	protocols::moves::Mover( PDBDiagnosticMover::mover_name() )
{
	TR << "PDBDiagnosticMover ctor" << std::endl;

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PDBDiagnosticMover::PDBDiagnosticMover( PDBDiagnosticMover const & src ):
	protocols::moves::Mover( src )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PDBDiagnosticMover::~PDBDiagnosticMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
PDBDiagnosticMover::apply( core::pose::Pose & pose ){
	//Diagnostic plan:
	//1) try reading PDB in (done by JD, already complete at this step)
	//2) calculate statistics based on the constellation of ResidueTypes in the Pose
	//3) try scoring
	//4) try packing
	//5) try minimizing

	using protocols::jd2::JobDistributor;
	protocols::jd2::JobOP job_me( JobDistributor::get_instance()->current_job() );
	std::string this_pdb_name(JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );
	core::Size const nres(pose.size());
	TR << this_pdb_name << " nres " << nres << std::endl;

	//step 2: accumulate residue type stats
	residue_type_statistics(pose, job_me, nres);

	//if option, cut and run! Useful for not-mega-clusters where packing a virus/ribosome will crush memory.
	if ( basic::options::option[ basic::options::OptionKeys::PDB_diagnostic::reading_only ].value() ) {
		return;
	}

	//create a ScoreFunction from commandline options
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

	//step 3: scoring
	using protocols::simple_moves::ScoreMoverOP;
	using protocols::simple_moves::ScoreMover;
	ScoreMoverOP score_mover(new ScoreMover(score_fxn));
	score_mover->set_verbose(false);
	score_mover->apply(pose);

	if ( !basic::options::option[ basic::options::OptionKeys::PDB_diagnostic::skip_pack_and_min ].value() ) {

		//packing crashes on zero nres
		if ( nres > 0 ) {

			//step 4: packing
			//note the main function forces -repack_only
			//set up SF and TF defaults (manually for some reason)
			//create a task factory: this will create a new PackerTask for each input pose
			using core::pack::task::operation::TaskOperationCOP;
			core::pack::task::TaskFactoryOP main_task_factory( new core::pack::task::TaskFactory );
			main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
			using protocols::minimization_packing::PackRotamersMoverOP;
			using protocols::minimization_packing::PackRotamersMover;
			PackRotamersMoverOP pack_rotamers(new protocols::minimization_packing::PackRotamersMover());
			pack_rotamers->task_factory( main_task_factory );
			pack_rotamers->score_function( score_fxn );
			pack_rotamers->apply(pose);
		} //skip all this if no residues!

		//step 5: minimizing
		using protocols::minimization_packing::MinMoverOP;
		using protocols::minimization_packing::MinMover;
		MinMoverOP min_mover(new protocols::minimization_packing::MinMover());
		min_mover->score_function( score_fxn );
		min_mover->apply(pose);
	}//skip_pack_and_min

	return;
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
PDBDiagnosticMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
PDBDiagnosticMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

}
void PDBDiagnosticMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	/* Option_Group( 'PDB_diagnostic',
	Option( 'reading_only', 'Boolean', desc="if true, become a no-op (only structure reading is tested). Useful for not-mega-clusters where packing a virus/ribosome will crush memory.", default="false" ),
	Option( 'skip_pack_and_min', 'Boolean', desc="Skip the packing and minimization steps (leaving in the scoring and ResidueType analyses). Good for speed, or for shallower tests.", default="false"),*/

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "DOCUMENTATION STRING", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PDBDiagnosticMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new PDBDiagnosticMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PDBDiagnosticMover::clone() const
{
	return protocols::moves::MoverOP( new PDBDiagnosticMover( *this ) );
}

std::string PDBDiagnosticMover::get_name() const {
	return mover_name();
}

std::string PDBDiagnosticMover::mover_name() {
	return "PDBDiagnosticMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
PDBDiagnosticMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new PDBDiagnosticMover );
}

std::string
PDBDiagnosticMoverCreator::keyname() const
{
	return PDBDiagnosticMover::mover_name();
}

void PDBDiagnosticMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PDBDiagnosticMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

///@details this function sums statistics: for each is_?? function in ResidueType, it accumulates how many residues of the Pose have that quality, and dumps that data to the scorefile via the Job's interface.
void PDBDiagnosticMover::residue_type_statistics( core::pose::Pose const & pose, protocols::jd2::JobOP job_me, core::Size const nres ){
	//I wonder if it's faster to pass nres in the function signature or just recalculate it?  I wonder if I've already spent more time thinking about it than the computer would spend, over all of human history, on the slower option?  I wonder if the wasted disk space for this comment outweighs having thought about it?

	//counter variable for each thing to track
	//core::Size is_??(0);
	core::Size is_polymer(0);
	core::Size is_sidechain_thiol(0);
	core::Size is_disulfide_bonded(0);
	core::Size is_sidechain_amine(0);
	core::Size is_protein(0);
	core::Size is_alpha_aa(0);
	core::Size is_beta_aa(0);
	core::Size is_gamma_aa(0);
	core::Size is_sri(0);
	core::Size is_triazolemer(0);
	core::Size is_d_aa(0);
	core::Size is_l_aa(0);
	core::Size is_achiral_backbone(0);
	core::Size is_DNA(0);
	core::Size is_RNA(0);
	core::Size is_coarse(0);
	core::Size is_NA(0);
	core::Size is_peptoid(0);
	core::Size is_carbohydrate(0);
	core::Size is_ligand(0);
	core::Size is_lipid(0);
	core::Size is_metal(0);
	core::Size is_metalbinding(0);
	core::Size is_membrane(0);
	core::Size is_surface(0);
	core::Size has_sc_orbitals(0);
	core::Size is_polar(0);
	core::Size is_charged(0);
	core::Size is_aromatic(0);
	core::Size is_cyclic(0);
	core::Size is_terminus(0);
	core::Size is_lower_terminus(0);
	core::Size is_upper_terminus(0);
	core::Size is_branch_point(0);
	core::Size is_acetylated_nterminus(0);
	core::Size is_methylated_cterminus(0);
	core::Size is_virtual_residue(0);
	core::Size is_adduct(0);

	//loop over pose residues, accumulating properties
	for ( core::Size i(1); i<=nres; ++i ) {
		//if(pose.residue_type(i).??) ++??;
		if ( pose.residue_type(i).is_polymer() ) ++is_polymer;
		if ( pose.residue_type(i).is_sidechain_thiol() ) ++is_sidechain_thiol;
		if ( pose.residue_type(i).is_disulfide_bonded() ) ++is_disulfide_bonded;
		if ( pose.residue_type(i).is_sidechain_amine() ) ++is_sidechain_amine;
		if ( pose.residue_type(i).is_protein() ) ++is_protein;
		if ( pose.residue_type(i).is_alpha_aa() ) ++is_alpha_aa;
		if ( pose.residue_type(i).is_beta_aa() ) ++is_beta_aa;
		if ( pose.residue_type(i).is_gamma_aa() ) ++is_gamma_aa;
		if ( pose.residue_type(i).is_sri() ) ++is_sri;
		if ( pose.residue_type(i).is_triazolemer() ) ++is_triazolemer;
		if ( pose.residue_type(i).is_d_aa() ) ++is_d_aa;
		if ( pose.residue_type(i).is_l_aa() ) ++is_l_aa;
		if ( pose.residue_type(i).is_achiral_backbone() ) ++is_achiral_backbone;
		if ( pose.residue_type(i).is_DNA() ) ++is_DNA;
		if ( pose.residue_type(i).is_RNA() ) ++is_RNA;
		if ( pose.residue_type(i).is_coarse() ) ++is_coarse;
		if ( pose.residue_type(i).is_NA() ) ++is_NA;
		if ( pose.residue_type(i).is_peptoid() ) ++is_peptoid;
		if ( pose.residue_type(i).is_carbohydrate() ) ++is_carbohydrate;
		if ( pose.residue_type(i).is_ligand() ) ++is_ligand;
		if ( pose.residue_type(i).is_lipid() ) ++is_lipid;
		if ( pose.residue_type(i).is_metal() ) ++is_metal;
		if ( pose.residue_type(i).is_metalbinding() ) ++is_metalbinding;
		if ( pose.residue_type(i).is_membrane() ) ++is_membrane;
		if ( pose.residue_type(i).is_surface() ) ++is_surface;
		if ( pose.residue_type(i).has_sc_orbitals() ) ++has_sc_orbitals;
		if ( pose.residue_type(i).is_polar() ) ++is_polar;
		if ( pose.residue_type(i).is_charged() ) ++is_charged;
		if ( pose.residue_type(i).is_aromatic() ) ++is_aromatic;
		if ( pose.residue_type(i).is_cyclic() ) ++is_cyclic;
		if ( pose.residue_type(i).is_terminus() ) ++is_terminus;
		if ( pose.residue_type(i).is_lower_terminus() ) ++is_lower_terminus;
		if ( pose.residue_type(i).is_upper_terminus() ) ++is_upper_terminus;
		if ( pose.residue_type(i).is_branch_point() ) ++is_branch_point;
		if ( pose.residue_type(i).is_acetylated_nterminus() ) ++is_acetylated_nterminus;
		if ( pose.residue_type(i).is_methylated_cterminus() ) ++is_methylated_cterminus;
		if ( pose.residue_type(i).is_virtual_residue() ) ++is_virtual_residue;
		if ( pose.residue_type(i).is_adduct() ) ++is_adduct;

	}

	//dump all these counts into scorefile
	//job_me->add_string_real_pair("??", ??);

	job_me->add_string_real_pair("nres", nres);

	job_me->add_string_real_pair("is_polymer", is_polymer);
	job_me->add_string_real_pair("is_sidechain_thiol", is_sidechain_thiol);
	job_me->add_string_real_pair("is_disulfide_bonded", is_disulfide_bonded);
	job_me->add_string_real_pair("is_sidechain_amine", is_sidechain_amine);
	job_me->add_string_real_pair("is_protein", is_protein);
	job_me->add_string_real_pair("is_alpha_aa", is_alpha_aa);
	job_me->add_string_real_pair("is_beta_aa", is_beta_aa);
	job_me->add_string_real_pair("is_gamma_aa", is_gamma_aa);
	job_me->add_string_real_pair("is_sri", is_sri);
	job_me->add_string_real_pair("is_triazolemer", is_triazolemer);
	job_me->add_string_real_pair("is_d_aa", is_d_aa);
	job_me->add_string_real_pair("is_l_aa", is_l_aa);
	job_me->add_string_real_pair("is_achiral_backbone", is_achiral_backbone);
	job_me->add_string_real_pair("is_DNA", is_DNA);
	job_me->add_string_real_pair("is_RNA", is_RNA);
	job_me->add_string_real_pair("is_coarse", is_coarse);
	job_me->add_string_real_pair("is_NA", is_NA);
	job_me->add_string_real_pair("is_peptoid", is_peptoid);
	job_me->add_string_real_pair("is_carbohydrate", is_carbohydrate);
	job_me->add_string_real_pair("is_ligand", is_ligand);
	job_me->add_string_real_pair("is_lipid", is_lipid);
	job_me->add_string_real_pair("is_metal", is_metal);
	job_me->add_string_real_pair("is_metalbinding", is_metalbinding);
	job_me->add_string_real_pair("is_membrane", is_membrane);
	job_me->add_string_real_pair("is_surface", is_surface);
	job_me->add_string_real_pair("has_sc_orbitals", has_sc_orbitals);
	job_me->add_string_real_pair("is_polar", is_polar);
	job_me->add_string_real_pair("is_charged", is_charged);
	job_me->add_string_real_pair("is_aromatic", is_aromatic);
	job_me->add_string_real_pair("is_cyclic", is_cyclic);
	job_me->add_string_real_pair("is_terminus", is_terminus);
	job_me->add_string_real_pair("is_lower_terminus", is_lower_terminus);
	job_me->add_string_real_pair("is_upper_terminus", is_upper_terminus);
	job_me->add_string_real_pair("is_branch_point", is_branch_point);
	job_me->add_string_real_pair("is_acetylated_nterminus", is_acetylated_nterminus);
	job_me->add_string_real_pair("is_methylated_cterminus", is_methylated_cterminus);
	job_me->add_string_real_pair("is_virtual_residue", is_virtual_residue);
	job_me->add_string_real_pair("is_adduct", is_adduct);

	return;
}

std::ostream &
operator<<( std::ostream & os, PDBDiagnosticMover const & mover )
{
	mover.show(os);
	return os;
}

} //devel
} //scientific_tests
