// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/AntibodyDesignProtocol.cc
/// @brief Handles the Antibody Design Protocol.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/AntibodyDesignProtocol.hh>
#include <protocols/antibody/design/AntibodyDesignProtocolCreator.hh>
#include <protocols/antibody/design/util.hh>

#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyDesignMover.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/clusters/util.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyNumberingConverterMover.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/snugdock/SnugDock.hh>
#include <protocols/antibody/snugdock/SnugDockProtocol.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <protocols/antibody/constraints/util.hh>
#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.hh>
#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.hh>
#include <protocols/antibody/constraints/CDRDihedralConstraintMover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/jd2/util.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/DeleteChainsMover.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>

//Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/BasicDataCache.hh>


#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.design.AntibodyDesignProtocol" );
namespace protocols {
namespace antibody {
namespace design {

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::scoring;
using namespace protocols::simple_moves;
using namespace core::import_pose;
using namespace protocols::antibody::clusters;

using core::pose::Pose;
using std::string;
using utility::vector1;

AntibodyDesignProtocol::AntibodyDesignProtocol() : protocols::moves::Mover(),
	graft_designer_(/* NULL */),
	cdr_dihedral_cst_mover_(/* NULL */),
	scorefxn_(/* NULL */),
	scorefxn_min_(/* NULL */),
	ab_info_(/* NULL */)

{
	protocols::moves::Mover::type( "AntibodyDesign" );
	read_cmd_line_options();

}


AntibodyDesignProtocol::~AntibodyDesignProtocol(){}

AntibodyDesignProtocol::AntibodyDesignProtocol( AntibodyDesignProtocol const & src ):
	protocols::moves::Mover( src ),
	run_graft_designer_( src.run_graft_designer_ ),
	run_snugdock_( src.run_snugdock_ ),
	run_relax_(src.run_relax_),
	remove_antigen_(src.remove_antigen_),
	instruction_file_( src.instruction_file_ )

{
	using namespace constraints;

	if ( src.ab_info_ ) ab_info_ = AntibodyInfoOP( new AntibodyInfo( *src.ab_info_ ));
	if ( src.scorefxn_ ) scorefxn_ = src.scorefxn_->clone();
	if ( src.scorefxn_min_ ) scorefxn_min_ = src.scorefxn_min_->clone();

	if ( src.graft_designer_ ) graft_designer_ = AntibodyDesignMoverOP( new AntibodyDesignMover( *src.graft_designer_ ));
	if ( src.cdr_dihedral_cst_mover_ ) cdr_dihedral_cst_mover_ = CDRDihedralConstraintMoverOP( new CDRDihedralConstraintMover( *src.cdr_dihedral_cst_mover_ ));


}


// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AntibodyDesignProtocolCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AntibodyDesignProtocol );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AntibodyDesignProtocolCreator::keyname() const {
// XRW TEMP  return AntibodyDesignProtocol::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AntibodyDesignProtocol::mover_name(){
// XRW TEMP  return "AntibodyDesignProtocol";
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP AntibodyDesignProtocol::get_name() const {
// XRW TEMP  return "AntibodyDesignProtocol";
// XRW TEMP }

//protocols::moves::MoverOP
//AntibodyDesignProtocol::clone() const {
// return protocols::moves::MoverOP( new AntibodyDesignProtocol(*this) );
//}


void
AntibodyDesignProtocol::read_cmd_line_options(){

	run_graft_designer_= true;
	set_run_snugdock(option [OptionKeys::antibody::design::run_snugdock]());
	set_run_relax(option [OptionKeys::antibody::design::run_relax]());

	remove_antigen_ = option [OptionKeys::antibody::design::remove_antigen]();

	if ( basic::options::option [basic::options::OptionKeys::antibody::design::cdr_instructions].user() ) {
		instruction_file_ = basic::options::option [basic::options::OptionKeys::antibody::design::cdr_instructions]();
		TR << "Instructions file: " << instruction_file_ << std::endl;
	}

}

protocols::moves::MoverOP
AntibodyDesignProtocol::clone() const {
	return protocols::moves::MoverOP( new AntibodyDesignProtocol(*this) );
}

protocols::moves::MoverOP
AntibodyDesignProtocol::fresh_instance() const {
	return protocols::moves::MoverOP( new AntibodyDesignProtocol );
}

void
AntibodyDesignProtocol::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data,
	Filters_map const & ,
	moves::Movers_map const & ,
	Pose const &
){
	using namespace core::scoring;
	AntibodyEnumManager manager = AntibodyEnumManager();

	scorefxn_ = get_ab_design_global_scorefxn(tag, data);
	scorefxn_min_ = get_ab_design_min_scorefxn(tag, data);


	//A little redundancy
	if ( tag->hasOption("instruction_file") ) {
		instruction_file_ = tag->getOption< std::string >("instruction_file");
	} else if ( tag->hasOption("instructions_file") ) {
		instruction_file_ = tag->getOption< std::string >("instructions_file");
	} else if ( tag->hasOption("cdr_instructions_file") ) {
		instruction_file_ = tag->getOption< std::string >("cdr_instructions_file");
	}

	if ( tag->hasOption("seq_design_cdrs") ) {
		utility::vector1<std::string> cdr_strings = utility::string_split_multi_delim(tag->getOption< std::string>("seq_design_cdrs"), ":,'`~+*&|;.");
		for ( core::Size i = 1; i <= cdr_strings.size(); ++i ) {
			CDRNameEnum cdr_enum =manager.cdr_name_string_to_enum( cdr_strings[ i ] );
			seq_design_override_.push_back(cdr_enum);
		}
	}

	if ( tag->hasOption("graft_design_cdrs") ) {
		utility::vector1<std::string> cdr_strings = utility::string_split_multi_delim(tag->getOption< std::string>("graft_design_cdrs"), ":,'`~+*&|;.");
		for ( core::Size i = 1; i <= cdr_strings.size(); ++i ) {
			CDRNameEnum cdr_enum =manager.cdr_name_string_to_enum( cdr_strings[ i ] );
			graft_design_override_.push_back(cdr_enum);
		}
	}
	if ( tag->hasOption("primary_cdrs") ) {
		utility::vector1<std::string> cdr_strings = utility::string_split_multi_delim(tag->getOption< std::string>("primary_cdrs"), ":,'`~+*&|;.");
		for ( core::Size i = 1; i <= cdr_strings.size(); ++i ) {
			CDRNameEnum cdr_enum =manager.cdr_name_string_to_enum( cdr_strings[ i ] );
			primary_cdrs_.push_back(cdr_enum);
		}
	}

	run_snugdock_ = tag->getOption<bool>("run_snugdock", run_snugdock_);
	run_relax_ = tag->getOption<bool>("run_relax", run_relax_);
	remove_antigen_ = tag->getOption<bool>("remove_antigen", remove_antigen_);

}

void
AntibodyDesignProtocol::setup_scorefxns(){
	using namespace basic::options;

	if ( ! scorefxn_ ) {
		scorefxn_ = get_ab_design_global_scorefxn();
	}

	if ( ! scorefxn_min_ ) {
		scorefxn_min_ = get_ab_design_min_scorefxn();
	}

}

void
AntibodyDesignProtocol::set_run_snugdock(bool setting){
	run_snugdock_ = setting;
}

void
AntibodyDesignProtocol::set_run_relax(bool setting){
	run_relax_ = setting;
}

void
AntibodyDesignProtocol::set_scorefxn(ScoreFunctionOP scorefxn){
	scorefxn_ = scorefxn;
}

void
AntibodyDesignProtocol::set_scorefxn_min(core::scoring::ScoreFunctionOP scorefxn){
	scorefxn_min_ = scorefxn;
}

void
AntibodyDesignProtocol::set_instruction_file_path(std::string instruction_file){
	instruction_file_ = instruction_file;
}

void
AntibodyDesignProtocol::model_post_design(core::pose::Pose & pose){

	//If no constraints are set - add constraints to pose CDRs.
	core::Real start = scorefxn_->score(pose);

	if ( run_snugdock_ ) {
		snugdock::SnugDockOP snug( new snugdock::SnugDock() );
		AntibodyInfoOP updated_ab_info( new AntibodyInfo(pose, AHO_Scheme, North) ); //Should be a reset method in AbInfo..

		snug->set_scorefxn(scorefxn_);
		snug->set_antibody_info(updated_ab_info);//Updated info for movemaps, foldtrees, etc.
		//snug->number_of_high_resolution_cycles(25); //Half as many cycles to save some time as we do a full dualspace relax next.
		//snug->debug();
		snug->apply(pose);
		core::Real post_snugdock = scorefxn_->score(pose);
		TR <<"start:            " << start <<std::endl;
		TR <<"postSD:           " << post_snugdock << std::endl;
	}

	if ( run_relax_ ) {
		protocols::relax::FastRelaxOP rel( new protocols::relax::FastRelax() );

		ScoreFunctionOP dualspace_scorefxn = scorefxn_->clone(); //May need to be strictly Talaris2013_cart;
		dualspace_scorefxn->set_weight_if_zero(cart_bonded, .5);
		dualspace_scorefxn->set_weight(pro_close, 0);

		rel->set_scorefxn(dualspace_scorefxn);
		rel->dualspace(true);
		rel->max_iter(200);
		rel->minimize_bond_angles(true);
		rel->apply(pose);

		core::Real post_modeling = scorefxn_->score(pose);

		TR <<"start:            " << start <<std::endl;
		TR <<"postDsREL:        " << post_modeling << std::endl;
	}
}

void
AntibodyDesignProtocol::setup_design_mover() {

	graft_designer_ = AntibodyDesignMoverOP( new AntibodyDesignMover( ab_info_) );
	graft_designer_->set_scorefunction(scorefxn_);
	graft_designer_->set_scorefunction_min(scorefxn_min_);
	graft_designer_->set_instruction_file(instruction_file_);

	if ( primary_cdrs_.size() > 0 ) {
		graft_designer_->set_primary_cdrs( primary_cdrs_ );
	}
	if ( seq_design_override_.size() >0 ) {
		graft_designer_->set_seq_design_cdrs( seq_design_override_ );
	}
	if ( graft_design_override_.size() > 0 ) {
		graft_designer_->set_graft_design_cdrs( graft_design_override_ );
	}

}

void
AntibodyDesignProtocol::reorder_poses(utility::vector1<core::pose::PoseOP>& poses){
	//Can be refactored to use utility::TopScoreSelector
	//From mc algorithm, you can have multiple poses that are equivalent...

	if ( poses.size() <= 1 ) return;
	utility::vector1<core::pose::PoseOP> sorted_poses;
	sorted_poses.push_back(poses[1]);

	for ( core::Size i = 1; i <= poses.size(); ++i ) {
		core::Real scU = (*scorefxn_)(*poses[i]);
		vector1<core::pose::PoseOP>::iterator pose_it = sorted_poses.begin();

		for ( core::Size x = 1; x <= sorted_poses.size(); ++x ) {
			core::Real scS = (*scorefxn_)(*sorted_poses[x]);
			if ( scU < scS ) {
				sorted_poses.insert(pose_it+x-1, poses[i]);
				break;
			}
			//If the unsorted pose is not lower in energy than any of the sorted poses, add it to the end
			if ( x == sorted_poses.size() ) {
				sorted_poses.push_back(poses[i]);
			}
		}
	}
	debug_assert(sorted_poses.size() == poses.size());
	poses = sorted_poses;
}

void
AntibodyDesignProtocol::output_ensemble(vector1<core::pose::PoseOP> ensemble, core::Size range_start, core::Size range_end, std::string prefix){

	for ( core::Size i = range_start; i <= range_end; ++i ) {
		//Filter here
		std::string tag = prefix+"_"+utility::to_string(i);
		TR << "Outputting ensemble " << i << std::endl;
		add_cluster_comments_to_pose( *(ensemble[ i ]), ab_info_ );
		check_fix_aho_cdr_numbering(ab_info_, *(ensemble[ i ]) );

		//Note that this does NOT write to a scorefile by default.  To pass it to the scorefile,
		// you need -other_pose_to_scorefile.  This is now working in MPI.

		if ( basic::options::option[ basic::options::OptionKeys::antibody::output_ab_scheme].user() ) {
			AntibodyNumberingConverterMover converter = AntibodyNumberingConverterMover();
			converter.apply(*(ensemble[i]));
		}

		protocols::jd2::output_intermediate_pose( *(ensemble[i]), tag );
	}
}

void
AntibodyDesignProtocol::apply(core::pose::Pose& pose){

	///Setup Objects///


	//if (! protocols::antibody::clusters::check_if_pose_renumbered_for_clusters(pose)){
	// utility_exit_with_message("PDB must be numbered correctly to identify North CDR clusters.  Please see Antibody Design documentation.");
	//}

	using namespace protocols::analysis;

	ab_info_ = AntibodyInfoOP( new AntibodyInfo(pose, AHO_Scheme, North) );
	ab_info_->show(std::cout);

	if ( remove_antigen_ && ab_info_->antigen_present() ) {
		DeleteChainsMover remove_chains_mover = DeleteChainsMover();
		remove_chains_mover.set_chains(ab_info_->get_antigen_chain_string(), pose);
		remove_chains_mover.apply( pose );

		//Reinit AbInfo.  Will call reinit function once we actually have that.
		ab_info_ = AntibodyInfoOP(new AntibodyInfo(pose, AHO_Scheme, North) );

		TR << "Antigen chains deleted from pose.." << std::endl;
	}

	setup_scorefxns();
	setup_design_mover();

	core::Real native_score = (*scorefxn_)(pose);
	TR <<"Score without constraints:" << std::endl;
	scorefxn_->show(TR, pose);

	//Will go in SHOW.
	TR << std::endl;
	TR << "/// Run Antibody Designer: "<< std::boolalpha << run_graft_designer_ <<std::endl;
	TR << "//      ->Post SnugDock:  "<< std::boolalpha << run_snugdock_ << std::endl;
	TR << "//      ->Post DS Relax:  "<< std::boolalpha << run_relax_ << std::endl;
	TR << std::endl;

	vector1<core::pose::PoseOP> pose_ensemble;
	vector1<core::pose::PoseOP> final_pose_ensemble;
	if ( run_graft_designer_ ) {
		TR <<"Running Graft Design Protocol" <<std::endl;
		graft_designer_->apply(pose);
		pose_ensemble = graft_designer_->get_top_designs();

		for ( core::Size i = 1; i <= pose_ensemble.size(); ++i ) {

			//Use a filter on ensembles generated?  Could have ~15 ensembles. What about the main pose?
			//If it passes filter, add to final_pose_ensemble.
			final_pose_ensemble.push_back(pose_ensemble[i]);
		}
	} else {
		//ab_info_->setup_CDR_clusters(pose);
		//current_constraint_result = protocols::antibody::add_harmonic_cluster_constraints(ab_info_, pose);
		final_pose_ensemble.push_back( core::pose::PoseOP( new Pose() ));
		*(final_pose_ensemble[1]) = pose;
	}

	//Any post-design modeling? Any filters?
	//Reorder pose ensemble before output
	//reorder_poses(final_pose_ensemble); Need debugging - no time to run debugger

	protocols::analysis::InterfaceAnalyzerMoverOP analyzer( new protocols::analysis::InterfaceAnalyzerMover(get_dock_chains_from_ab_dock_chains(ab_info_, "A_LH"), false, scorefxn_, false , true, false) );

	if ( run_snugdock_ || run_relax_ ) {
		for ( core::Size i = 1; i <= final_pose_ensemble.size(); ++i ) {
			ab_info_->setup_CDR_clusters( *final_pose_ensemble[i], false );
			add_cluster_comments_to_pose( *final_pose_ensemble[i], ab_info_ );
			check_fix_aho_cdr_numbering( ab_info_, *final_pose_ensemble[i] );

			if ( option [ OptionKeys::antibody::design::run_interface_analyzer ]() ) {
				analyzer->init_on_new_input(*final_pose_ensemble[i]);
				analyzer->apply(*final_pose_ensemble[i]);
				analyzer->add_score_info_to_pose(*final_pose_ensemble[i]);
			}

			scorefxn_->score( *final_pose_ensemble[i] );
			output_ensemble( final_pose_ensemble, 1, final_pose_ensemble.size(), "pre_model" );
		}
	}


	TR << "Running Post-Graft modeling" << std::endl;
	for ( core::Size i = 1; i <= final_pose_ensemble.size(); ++i ) {
		this->model_post_design(*final_pose_ensemble[ i ]);
	}

	TR <<"Native: "<< native_score << std::endl;
	for ( core::Size i = 1; i <= final_pose_ensemble.size(); ++i ) {
		TR << "Pose " << i << std::endl;
		scorefxn_->show( TR, *final_pose_ensemble[ i ] );
		ab_info_->setup_CDR_clusters( *final_pose_ensemble[ i ], false );
		add_cluster_comments_to_pose( *final_pose_ensemble[ i ], ab_info_ );
		check_fix_aho_cdr_numbering( ab_info_, *final_pose_ensemble[ i ] );

		if ( option [ OptionKeys::antibody::design::run_interface_analyzer ]() ) {
			analyzer->init_on_new_input(*final_pose_ensemble[i]);
			analyzer->apply(*final_pose_ensemble[i]);
			analyzer->add_score_info_to_pose(*final_pose_ensemble[i]);
		}

	}

	//Output final ensembles.
	pose = *final_pose_ensemble[1];
	if ( final_pose_ensemble.size() > 1 ) {
		output_ensemble(final_pose_ensemble, 2, final_pose_ensemble.size(), "ensemble");
	}

	if ( ! option [OptionKeys::out::file::pdb_comments]() ) {
		TR << "Added pose comments for cluster info.  Use -pdb_comments option to have it output to pdb file." << std::endl;
	}

	if ( basic::options::option[ basic::options::OptionKeys::antibody::output_ab_scheme].user() ) {
		AntibodyNumberingConverterMover converter = AntibodyNumberingConverterMover();
		converter.apply(pose);
	}

	TR << "Complete." << std::endl;



}

std::string AntibodyDesignProtocol::get_name() const {
	return mover_name();
}

std::string AntibodyDesignProtocol::mover_name() {
	return "AntibodyDesignProtocol";
}

void AntibodyDesignProtocol::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	AttributeList attlist;

	attributes_for_get_ab_design_global_scorefxn(attlist);
	attributes_for_get_ab_design_min_scorefxn(attlist);

	XMLSchemaRestriction ABcdr_enum;
	ABcdr_enum.name("ABcdr_definitions");
	ABcdr_enum.base_type(xs_string);
	AntibodyEnumManager ABenum_manager;
	for ( auto& cdr_def : ABenum_manager.get_recognized_cdr_definitions() ) {
		ABcdr_enum.add_restriction(xsr_enumeration, cdr_def);
	}
	xsd.add_top_level_element(ABcdr_enum);

	attlist + XMLSchemaAttribute(
		"seq_design_cdrs", xs_string,
		"CDR regions to be Sequence-Designed");

	attlist + XMLSchemaAttribute(
		"graft_design_cdrs", xs_string,
		"CDR regions to be Graft-Designed");

	attlist + XMLSchemaAttribute(
		"primary_cdrs", xs_string,
		"Manually set the CDRs which can be chosen in the outer cycle. \n"
		"These should be on for either Sequence-Design or Graft-Design. \n"
		"Normally, the outer cycles are whatever CDRs we are designing, including CDRs which are sequence-design only.  \n"
		"Use this if you are primarily interested in specific CDRs (such as graft-designing H3 and allowing H1 and L3 to sequence design during the inner cycle.)");

	attlist + XMLSchemaAttribute(
		"instruction_file", xs_string,
		"XSD XRW: TO DO");
	attlist + XMLSchemaAttribute(
		"instructions_file", xs_string,
		"used if not instruction_file attribute specified");
	attlist + XMLSchemaAttribute(
		"cdr_instructions_file", xs_string,
		"used if not instructions_file attribute specififed");

	attlist + XMLSchemaAttribute(
		"run_snugdock", xsct_rosetta_bool,
		"Run snugdock at the end");
	attlist + XMLSchemaAttribute(
		"run_relax", xsct_rosetta_bool,
		"Run relax at the end and after snugdock (optional)");
	attlist + XMLSchemaAttribute(
		"remove_antigen", xsct_rosetta_bool,
		"Removes antigen at the very beginning of the protocol");

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Main AntibodyDesignProtocol, tieing together individual movers and classes. "
		"Main mover for application.",
		attlist );
}

std::string AntibodyDesignProtocolCreator::keyname() const {
	return AntibodyDesignProtocol::mover_name();
}

protocols::moves::MoverOP
AntibodyDesignProtocolCreator::create_mover() const {
	return protocols::moves::MoverOP( new AntibodyDesignProtocol );
}

void AntibodyDesignProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AntibodyDesignProtocol::provide_xml_schema( xsd );
}



} //Design
} //Antibody
} //Protocols
