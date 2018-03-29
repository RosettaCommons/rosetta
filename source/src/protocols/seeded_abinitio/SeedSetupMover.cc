// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief mover to set up packer task and movemap for seededabinitio
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

// Unit headers
#include <protocols/seeded_abinitio/SeedSetupMover.hh>
#include <protocols/seeded_abinitio/SeedSetupMoverCreator.hh>
#include <protocols/seeded_abinitio/SeededAbinitio_util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <basic/datacache/DataMap.hh>


#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/task_operations/RestrictChainToRepackingOperation.hh>
#include <protocols/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/task_operations/PreventResiduesFromRepackingOperation.hh>
#include <protocols/task_operations/RestrictResiduesToRepackingOperation.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/types.hh>

//scoring
#include <core/scoring/ScoreFunction.hh>


#include <basic/Tracer.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <protocols/scoring/Interface.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

#include <core/chemical/AtomType.hh>
#include <core/conformation/Conformation.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace seeded_abinitio {

using namespace core;

static basic::Tracer TR( "protocols.seeded_abinitio.SeedSetupMover" );

// XRW TEMP std::string
// XRW TEMP SeedSetupMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SeedSetupMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SeedSetupMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SeedSetupMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SeedSetupMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "SeedSetupMover";
// XRW TEMP }

SeedSetupMover::~SeedSetupMover() = default;

SeedSetupMover::SeedSetupMover() :
	protocols::moves::Mover( SeedSetupMover::mover_name() ){
	//movemap options
	chi_chain1_ = false;
	chi_chain2_ = false;

	//packer tasks options
	packtask_ = false;//until it is more stable
	repack_target_ = true;
	repack_foldpose_ = true;
	design_target_ = false;
	design_foldpose_ = true;
	allow_all_aas_ = false;
	norepack_res_ = "";
	design_res_ = "" ;
	//clearing containers
	//task_factory_->clear();
}

protocols::moves::MoverOP
SeedSetupMover::clone() const {
	return( protocols::moves::MoverOP( new SeedSetupMover( *this ) ) );
}

protocols::moves::MoverOP
SeedSetupMover::fresh_instance() const {
	return protocols::moves::MoverOP( new SeedSetupMover );
}

void
SeedSetupMover::task_factory( core::pack::task::TaskFactoryOP tf ) {
	task_factory_ = tf;
}


void SeedSetupMover::clear_task_factory(){
	if ( task_factory_ ) {
		task_factory_->clear();
	}
}

/*
void
SeedSetupMover::clear_task(){
task_ = NULL;
}
*/

/// getters
core::pack::task::TaskFactoryOP &
SeedSetupMover::task_factory() {
	return task_factory_;
}


bool is_part ( utility::vector1< core::Size > vec, core::Size pos ){
	bool is_part_of_this = false;
	for ( core::Size i = 1; i <= vec.size(); i++ ) {
		if ( i == pos ) is_part_of_this = true;
	}
	return is_part_of_this;
}


void
SeedSetupMover::set_packerTasks_target_and_seeds (  core::pose::Pose & pose ,
	protocols::loops::Loops & seeds,
	utility::vector1< core::Size > & designable_residues, // ){
	utility::vector1< core::Size > & norepack_res, // ){//,
	core::pack::task::TaskFactoryOP & tf ){

	//primitive starting version of the packer tasks for seeded_abinitio
	//eventually, this shoudl be more sophisticated

	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace protocols::task_operations;

	Size num_chains = pose.conformation().num_chains();


	///////////////////////////////
	///1. set repacking behavior:
	///////////////////////////////

	TR<<"disallowing target chain(s) to repack" << std::endl;
	if ( num_chains > 1  && repack_target_ ) {
		PreventChainFromRepackingOperationOP pcfr( new PreventChainFromRepackingOperation ) ;
		Size chains_norepack = num_chains;
		if ( repack_foldpose_ ) {
			chains_norepack -= 1;
			TR << "allowing repacking of fold pose" <<std::endl;
		}
		for ( Size chain = 1; chain <= chains_norepack ; chain++ ) {
			pcfr->chain( chain );
			//pcfr.apply( pose, *tasks );
		}
		tf->push_back( pcfr );
	}

	/// disallow specified residues to repack:
	if ( norepack_res.size() != 0 ) {
		TR.Debug<<"disallow " << norepack_res.size() << " to repack " <<std::endl;
		using namespace protocols::task_operations;
		PreventResiduesFromRepackingOperationOP prfr( new PreventResiduesFromRepackingOperation ) ;
		prfr->set_residues( norepack_res );
		tf->push_back( prfr );
	}

	/*
	/// disallow specified residues to repack:
	for( const Size res , norepack_res ){
	using namespace core::pack::task::operation;
	PreventRepackingOP pr = new PreventRepacking;
	pr->include_residue( res );
	tf->push_back( pr );
	}
	*/
	/// do not repack disulfides
	NoRepackDisulfidesOP nrd( new NoRepackDisulfides );
	tf->push_back( nrd );


	//////////////////////////////////
	/// 2. restrict design:
	//////////////////////////////////

	if ( !design_target_ ) {
		RestrictChainToRepackingOperationOP rctr( new RestrictChainToRepackingOperation );
		for ( Size chain = 1; chain <= num_chains - 1 ; chain++ ) {
			rctr->chain( chain );
		}
		tf->push_back( rctr );
	}
	if ( design_target_ ) {
		TR<<"WARNING, are you sure you want to design the target chains? designing target chain"  << std::endl;
	}


	/// setting up design restrictions for new fold pose with seeds:
	if ( design_foldpose_ ) {
		utility::vector1< core::Size > residues;
		residues.clear();

		for ( core::Size i = pose.conformation().chain_begin( num_chains ); i <= pose.size(); i++ ) {
			bool allow_design( false );
			if ( !seeds.is_loop_residue( i ) || is_part( designable_residues, i ) ) {
				TR.Debug<<"is not seed: "<< i << " allow design" << std::endl;
				allow_design = true;
			}
			if ( !allow_design ) {
				residues.push_back( i );
			}
		}

		if ( residues.size() ) {
			TR.Debug<<"The following residues will be repacked only: ";
			for ( core::Size const res : residues ) {
				TR.Debug<<res<<", ";
			}
			TR.Debug<<std::endl;
			RestrictResiduesToRepackingOperationOP rrtr( new RestrictResiduesToRepackingOperation );
			rrtr->set_residues( residues );
			tf->push_back( rrtr );
		}
	}
}


void define_movemap_chains(
	core::pose::Pose & pose,
	core::kinematics::MoveMapOP & movemap,
	protocols::loops::Loops & seeds,
	bool chi_chain1,
	bool chi_chain2,
	bool interface_chi1,
	bool,
	core::Real interface_distance_cutoff
){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	Size num_chains = pose.conformation().num_chains();

	movemap->set_bb( false );
	movemap->set_chi( false );
	movemap->set_jump( false );

	//setting chi to true if specified
	if ( num_chains > 1 ) {

		//setting up interface object
		//multy chain possible, make jump a variable
		Size rb_jump = 1;
		//kinematics::FoldTree ori_ft = pose.fold_tree();
		pose.update_residue_neighbors(); // o/w fails assertion `graph_state_ == GOOD`
		protocols::scoring::Interface interface_obj(rb_jump );
		interface_obj.distance( interface_distance_cutoff );
		interface_obj.calculate( pose );
		//protocols::scoring::Interface interface_obj(rb_jump);

		TR<<"disallowing all chains to be movable but the last one" <<std::endl;
		for ( Size i = 1; i <= pose.conformation().chain_end( num_chains - 1 ); ++i ) {
			//to be safe
			movemap->set_bb( i , false );
			movemap->set_chi( i, false );
			if ( chi_chain1 ) {
				movemap->set_chi( i, true );
			}
			if ( interface_chi1 && interface_obj.is_interface(i) ) {
				movemap->set_chi( i, true );
			}
			//TR.Debug<<"i: "<<i << ", bb: "<< movemap->get_bb(i)<<", chi: "<<movemap->get_chi(i)<<  std::endl;
		}
	}
	//iterating through the last chain (this is the one that gets folded, and allowing the different degrees of freedom
	for ( Size pos = pose.conformation().chain_begin( num_chains ); pos <= pose.conformation().chain_end( num_chains ) ; pos++ ) {
		if ( !seeds.is_loop_residue( pos ) ) {
			movemap->set_bb( pos, true );
			movemap->set_chi(pos, true );
		}
		if ( seeds.is_loop_residue( pos ) && chi_chain2 ) {
			movemap->set_chi( true );
		}
	}

	TR<<"setting movemap method:"<<std::endl;
	for ( Size i = 1 ; i <=pose.size(); ++i ) {
		TR.Debug<<"i: "<<i<<", bb: "<< movemap->get_bb(i)<< ", chi: "<<movemap->get_chi(i)<<std::endl;
	}

	//option for movable residues within the seeds?
	//coudl be something like so:
	//<movableSeed_residues N=1 C=1/>
} //end define movemap

///adjustment since parse time specified residues are different numbered than run time residues
utility::vector1< core::Size >
adjust_des_residues( pose::Pose & pose,
	std::string design_residues ){

	utility::vector1< std::string > const design_keys( utility::string_split( design_residues, ',' ) );
	utility::vector1< core::Size > design_res;

	for ( std::string const & key : design_keys ) {
		core::Size const resnum( core::pose::parse_resnum( key, pose ));
		TR.Debug<<"design within seed, residue: "<< key <<", parsed: "<< resnum <<std::endl;
		design_res.push_back( resnum);
		TR<<"allowing design for "<<key<<std::endl;
	}
	//TR.Debug<<"runtime designable: " << design_res <<std::endl;
	return design_res;
}//end parsing design residues

void
SeedSetupMover::apply( core::pose::Pose & pose ){

	using namespace core::pack;
	using namespace core::pack::task;

	/// re-parsing seed residues at runtime since they might have changed within the trajectory due to length changes

	/// 1. re-parsing the input elements
	utility::vector1 <Size > designable_residues_motif;
	all_seeds_ = parse_seeds( pose , seed_vector_ );
	designable_residues_motif = adjust_des_residues( pose , design_res_);
	utility::vector1 <Size> norepack_residues = adjust_des_residues( pose, norepack_res_ );

	/// 2. set movemap
	debug_assert( movemap_factory_ );
	core::kinematics::MoveMapOP movemap;
	if ( movemap_factory_ ) {
		movemap = movemap_factory_->create_movemap_from_pose( pose );
	} else {
		movemap = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
	}
	define_movemap_chains( pose , movemap , all_seeds_, chi_chain1_, chi_chain2_, interface_chi1_, interface_chi2_, interface_distance_cutoff_ );

	/// 3. compute new task operations for seeds and target

	//reset taskfactory
	clear_task_factory();

	// todo: write setters adn getters for teh private variables....
	// very hacky and inconsistent use of private variables
	if ( packtask_ ) {
		set_packerTasks_target_and_seeds( pose, all_seeds_ , designable_residues_motif , norepack_residues , task_factory_ ); //, repack_target_, design_foldpose_, design_target_ )

		//TR<<"at the end of apply movemap:"<<std::endl;
		for ( Size i = 1 ; i <= pose.size(); ++i ) {
			TR.Debug <<"position: "<< i <<", bb: "<< movemap->get_bb(i)<< ", chi: "<<movemap->get_chi(i)<<std::endl;
		}

		PackerTaskOP ptask = task_factory_->create_task_and_apply_taskoperations( pose );
		if ( design_ ) {
			pack_rotamers( pose, *scorefxn_repack_ , ptask );
		}
	}

	TR.flush();
}//end apply

// XRW TEMP std::string
// XRW TEMP SeedSetupMover::get_name() const {
// XRW TEMP  return SeedSetupMover::mover_name();
// XRW TEMP }

void
SeedSetupMover::parse_my_tag( TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const &,
	core::pose:: Pose const & ){

	TR<<"SeedSetupMover has been invoked"<<std::endl;

	//adding the movemap to the datamap
	movemap_factory_ = protocols::rosetta_scripts::parse_movemap_factory_legacy( tag, data );

	//adding the taskfactory to the datamap
	//task_factory_ = new core::pack::task::TaskFactory;
	//temporarily inactivating this option.....
	//task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data, task_factory_ ) );

	chi_chain2_ = tag->getOption< bool >("chi_chain2", false );
	chi_chain1_ = tag->getOption< bool >("chi_chain1", false );

	interface_chi1_ = tag->getOption< bool >("interface_chi1", false );
	interface_chi2_ = tag->getOption< bool >("interface_chi2", false );
	interface_distance_cutoff_ = tag->getOption< core::Real >("interface_distance_cutoff" , 8 );

	//repacking for packer tasks options
	repack_target_ = tag->getOption< bool >("repack_target", true );
	repack_foldpose_ = tag->getOption< bool >("repack_foldpose", true );

	//for design option for packertasks
	design_target_ = tag->getOption< bool >("design_target", false );
	design_foldpose_ = tag->getOption< bool >("design_foldpose", true );
	allow_all_aas_ = tag->getOption< bool >("allow_all_aas", false );

	/// this mover can perform repacking and design mover -- which is not recommaned
	scorefxn_repack_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_repack", data )->clone();
	scorefxn_minimize_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn_minimize", data )->clone();
	design_ = tag->getOption< bool >( "design" , false );

	/// read input seeds
	utility::vector0< TagCOP > const & branch_tags( tag->getTags() );
	for ( TagCOP const btag : branch_tags ) {

		if ( btag->getName() == "Seeds" ) { //need an assertion for the presence of these or at least for the option file

			std::string const beginS( btag->getOption<std::string>( "begin" ) );
			std::string const endS( btag->getOption<std::string>( "end" ) );
			std::pair <std::string,std::string> seedpair;
			seedpair.first  = beginS;
			TR.Debug <<"parsing seeds: " << beginS << " " <<endS <<std::endl;
			seedpair.second = endS;
			seed_vector_.push_back( seedpair );
		}//end seeds
	}//end b-tags

	//allow design within the seed/motif
	if ( tag->hasOption( "allow_design" ) ) {
		design_res_ = tag->getOption< std::string >( "allow_design"  );
	}

	//dont allow these residues to repack
	if ( tag->hasOption( "norepack_res" ) ) {
		norepack_res_ = tag->getOption< std::string >( "norepack_res" );
	}

	if ( tag->hasOption("packtask") ) {
		packtask_ = tag->getOption< bool >( "packtask" );
		TR<<"packer task factory set to: " << std::endl;
	}

}//end parse my tag

std::string SeedSetupMover::get_name() const {
	return mover_name();
}

std::string SeedSetupMover::mover_name() {
	return "SeedSetupMover";
}

void SeedSetupMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default("chi_chain2", xsct_rosetta_bool, "XRW TO DO","0")
		+ XMLSchemaAttribute::attribute_w_default("chi_chain1", xsct_rosetta_bool, "XRW TO DO","0")
		+ XMLSchemaAttribute::attribute_w_default("interface_chain1", xsct_rosetta_bool, "XRW TO DO","0")
		+ XMLSchemaAttribute::attribute_w_default("interface_chain2", xsct_rosetta_bool, "XRW TO DO","0")
		+ XMLSchemaAttribute::attribute_w_default("interface_distance_cutoff", xsct_real, "XRW TO DO","0")
		+ XMLSchemaAttribute::attribute_w_default("repack_target", xsct_rosetta_bool, "XRW TO DO","1")
		+ XMLSchemaAttribute::attribute_w_default("repack_foldpose", xsct_rosetta_bool, "XRW TO DO","1")
		+ XMLSchemaAttribute::attribute_w_default("design_target", xsct_rosetta_bool, "XRW TO DO","0")
		+ XMLSchemaAttribute::attribute_w_default("design_foldpose", xsct_rosetta_bool, "XRW TO DO","1")
		+ XMLSchemaAttribute::attribute_w_default("allow_all_aas", xsct_rosetta_bool, "XRW TO DO","0");

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist, "scorefxn_repack" );
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist, "scorefxn_minimize" );
	attlist
		+ XMLSchemaAttribute::attribute_w_default("design", xsct_rosetta_bool, "XRW TO DO","0");

	//Subelements
	XMLSchemaSimpleSubelementList subelement_list;
	AttributeList subelement_attributes;
	subelement_attributes
		+ XMLSchemaAttribute::required_attribute("begin", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute::required_attribute("end", xs_string, "XRW TO DO");
	subelement_list.add_simple_subelement("Seeds", subelement_attributes, "XRW TO DO");

	attlist
		+ XMLSchemaAttribute("allow_design", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute("norepack_res", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute("packtask", xsct_rosetta_bool, "XRW TO DO");

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, subelement_list );
}

std::string SeedSetupMoverCreator::keyname() const {
	return SeedSetupMover::mover_name();
}

protocols::moves::MoverOP
SeedSetupMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SeedSetupMover );
}

void SeedSetupMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SeedSetupMover::provide_xml_schema( xsd );
}


}
}//protocol


