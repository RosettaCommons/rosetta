// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file
/// @brief
/// @author Jorge Fallas ( jaf18@uw.edu )

// Unit headers
#include <protocols/matdes/SchemePlaceMotifsMover.hh>
#include <protocols/matdes/SchemePlaceMotifsMoverCreator.hh>
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <numeric/random/random_xyz.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <protocols/matdes/MotifHitsRotamersOperation.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/symmetry/util.hh>

static basic::Tracer TR("protocols.matdes.SchemePlaceMotifsMover");

//static numeric::random::RandomGenerator::RandomGenerator() ;

namespace protocols {
namespace matdes {

using namespace core;
using namespace utility;
using core::pose::Pose;


// -------------  Mover Creator -------------
// XRW TEMP std::string
// XRW TEMP SchemePlaceMotifsMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return SchemePlaceMotifsMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SchemePlaceMotifsMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SchemePlaceMotifsMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SchemePlaceMotifsMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "SchemePlaceMotifs";
// XRW TEMP }

// -------------  Mover Creator -------------

SchemePlaceMotifsMover::SchemePlaceMotifsMover() : halt_on_error_(false)
{
	core::scoring::motif::MotifHashManager::get_instance();
	scorefxn_ = core::scoring::get_score_function()     ;
	task_factory_ = NULL  ;
}

protocols::moves::MoverOP
SchemePlaceMotifsMover::clone() const {
	return protocols::moves::MoverOP(new SchemePlaceMotifsMover( *this ));
}

protocols::moves::MoverOP
SchemePlaceMotifsMover::fresh_instance() const {
	return protocols::moves::MoverOP(new SchemePlaceMotifsMover());
}


void
SchemePlaceMotifsMover::apply(Pose & pose) {

	////////////////////////////////////////////////////////////////////////////
	///////////////////////// get matching motifs //////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	core::scoring::motif::ResPairMotifQuery q(pose);
	core::scoring::motif::MotifHits hits;
	if ( motif_sets_.size() ) {
		TR << "using motif set names from XML" << std::endl;
		for ( std::string const & msetname : motif_sets_ ) {
			core::scoring::motif::MotifHash const & mh( *core::scoring::motif::MotifHashManager::get_instance()->get_motif_hash_by_fname(msetname) );
			int n = mh.get_matching_motifs(q,hits);
			TR << n << " new hits from mset: " << msetname << ", total: " << hits.size() << std::endl;
		}
	} else {
		TR << "using all motif sets from CLI" << std::endl;
		core::scoring::motif::MotifHashManager::get_instance()->get_matching_motifs(q,hits);

	}
	if ( halt_on_error_ && hits.size()==0 ) {
		utility_exit_with_message("no matching motifs found!");
	}
	std::cout << "found " << hits.size() << " matching motifs" << std::endl;
	if ( dumpfile_!="" && hits.size() ) {
		for ( Size i = 0; i < 1; ++i ) std::cout << "!!!!!!!!!!!!!!! DUMPING MOTIFS " << dumpfile_ << "!!!!!!!!!!!!!!!!" << std::endl;
		static int dump_count = 0;
		hits.dump_motifs_pdb(dumpfile_+"_"+ObjexxFCL::string_of(++dump_count)+".pdb");
	}

	////////////////////////////////////////////////////////////////////////////
	///////////////////////// remove interface sidechains //////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	utility::vector1<bool> has_motif(pose.size(),false) ;
	for ( core::scoring::motif::MotifHit const & hit : hits ) {
		if ( hit.motif.type1() == core::scoring::motif::RM_BB ) has_motif[hit.residue1] = true ;
		if ( hit.motif.type2() == core::scoring::motif::RM_BB )  has_motif[hit.residue2] = true ;
	}

	core::pack::task::PackerTaskOP ptask = task_factory_->create_task_and_apply_taskoperations(pose);
	if ( mode_ == "basic" ) {
		for ( core::Size ir = 1; ir <= pose.size(); ++ir ) {
			if ( pose.residue(ir).is_protein() && has_motif[ir] ) {
				core::chemical::ResidueTypeCOP new_restype( core::pose::get_restype_for_pose(pose, "ALA", pose.residue_type(ir).mode() ) );
				core::pose::replace_pose_residue_copying_existing_coordinates( pose, ir, *new_restype );
			}
		}
		//core::pack::task::PackerTaskOP ptask = task_factory_->create_task_and_apply_taskoperations(pose);
		for ( core::Size ir =1 ; ir <= pose.size(); ++ir )  ptask->nonconst_residue_task(ir).restrict_to_repacking()      ;
	} else {
		utility::vector1<bool>  dummy(20,false)    ;
		utility::vector1<utility::vector1<bool> >  aas_for_design(pose.size(),dummy)    ;
		TR << "about to restrict aas to motifs" << std::endl ;
		if ( allowed_aas_ == "motifs" ) {
			for ( core::scoring::motif::MotifHit const & hit : hits ) {
				aas_for_design.at(hit.residue2).at(chemical::aa_from_oneletter_code(hit.motif.aa2())) = true ;
				aas_for_design.at(hit.residue1).at(chemical::aa_from_oneletter_code(hit.motif.aa1())) = true ;
			}
			for ( size_t i = 1; i <= aas_for_design.size(); i++ ) TR << aas_for_design[i] << std::endl ;
			//core::pack::task::PackerTaskOP ptask = task_factory_->create_task_and_apply_taskoperations(pose);
			for ( core::Size ir = 1; ir <= pose.size(); ++ir ) {
				ptask->nonconst_residue_task(ir).restrict_absent_canonical_aas(aas_for_design[ir]) ;
			}
		}
	}
	core::pack::rotamer_set::RotamerSetOperationOP    mot_rot = core::pack::rotamer_set::RotamerSetOperationOP(new MotifHitsRotamersOperation(hits))  ;
	//core::pack::task::PackerTaskOP ptask = task_factory_->create_task_and_apply_taskoperations(pose);
	ptask->append_rotamerset_operation(mot_rot)           ;
	TR << *ptask << std::endl;
	TR << "about to pack" << std::endl ;
	if ( core::pose::symmetry::is_symmetric(pose) )  protocols::simple_moves::symmetry::SymPackRotamersMover(scorefxn_, ptask).apply(pose)   ;
	else protocols::simple_moves::PackRotamersMover(scorefxn_, ptask).apply(pose)       ;
	////////////////////////////////////////////////////////////////////////////
	//////////////////////// add cloud constraints /////////////////////////////
	////////////////////////////////////////////////////////////////////////////
	/*
	core::scoring::methods::FloatingPointCloudRestraintEnergyFuncType type = core::scoring::methods::FloatingPointCloudRestraintEnergyFuncType_HARMINIC;

	core::scoring::methods::FloatingPointCloudRestraintEnergyInfoOP fpcinfo;
	fpcinfo = new core::scoring::methods::FloatingPointCloudRestraintEnergyInfo;
	pose.data().set( core::pose::datacache::CacheableDataType::FLOATING_POINT_CLOUD_INFO, fpcinfo );
	utility::vector1< core::scoring::methods::FloatingPoints > & clouds = fpcinfo->point_clouds();
	for ( core::scoring::motif::MotifHit const & mh : hits ){
	// numeric::xyzTransform<core::Real> x = numeric::random::gaussian_random_xform(360.0,20.0);
	core::scoring::methods::FloatingPoints cloud;
	for (std::string const & atomname : utility::string_split("N CA C O CB",' ')){
	if( pose.residue(mh.residue1).has(atomname) && mh.mpose().residue(1).has(atomname) ){
	cloud.push_back( core::scoring::methods::FloatingPoint(
	mh.mpose().residue(1).xyz(atomname),
	core::id::AtomID( pose.residue(mh.residue1).atom_index(atomname), mh.residue1 ),
	1.0,
	type
	) );
	}
	if( pose.residue(mh.residue2).has(atomname) && mh.mpose().residue(2).has(atomname) ){
	cloud.push_back( core::scoring::methods::FloatingPoint(
	mh.mpose().residue(2).xyz(atomname),
	core::id::AtomID( pose.residue(mh.residue2).atom_index(atomname), mh.residue2 ),
	1.0,
	type
	) );
	}
	}
	clouds.push_back(cloud);
	// mh.mpose().dump_pdb("mtest.pdb");
	// utility_exit_with_message("foo");
	}
	fpcinfo->reinitialize_point_map();
	pose.dump_pdb("test.pdb");

	// fpcinfo->dump_debug_pdb("fcp_test0.pdb");
	// fpcinfo->align_clouds_to_pose(pose);
	// fpcinfo->dump_debug_pdb("fcp_test1.pdb");
	// utility_exit_with_message("foo");
	*/
}

void
SchemePlaceMotifsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & /*filtermap*/,
	protocols::moves::Movers_map const & /*movermap*/,
	core::pose::Pose const & /*pose*/
){
	tag_name_ = tag->getOption< std::string >("name","");
	dumpfile_ = tag->getOption<std::string>("dumpfile","");
	mode_ = tag->getOption<std::string>("mode","basic");
	allowed_aas_ = tag->getOption<std::string>("allowed_aas","motifs");
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data ) ;
	task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data ) ;
	halt_on_error_ = tag->getOption<bool>("halt_on_error",false);
	std::string const motif_sets = tag->getOption< std::string >("motif_sets","");
	TR << "checking motif sets: " << motif_sets << std::endl;
	for ( std::string const & mset : utility::string_split(motif_sets,',') ) {
		if ( mset=="" ) continue;
		if ( !core::scoring::motif::MotifHashManager::get_instance()->have_motif_set_named(mset) ) {
			for ( std::string const & s : core::scoring::motif::MotifHashManager::get_instance()->motif_set_names() ) {
				std::cout << "candidate: " << s << std::endl;
			}
			utility_exit_with_message("MotifResiduesOperation: unknown motif set: "+mset+"\nspecify -mh:path:motifs");
		}
		motif_sets_.push_back(mset);
	}
	if ( motif_sets_.size()==0 ) {
		TR.Warning << "no motif sets in XML " << tag_name_ << ", will use all from cli" << std::endl;
		// utility_exit_with_message("MotifResiduesOperation: no motif sets specified, use motif_sets='foo,bar,baz'!");
	}


}

std::string SchemePlaceMotifsMover::get_name() const {
	return mover_name();
}

std::string SchemePlaceMotifsMover::mover_name() {
	return "SchemePlaceMotifs";
}

void SchemePlaceMotifsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "name", xs_string, "Tag name.", "XRW TO DO") //should this be default?
		+ XMLSchemaAttribute::attribute_w_default( "dumpfile", xs_string, "This is an option to generate a extra pdb that dumps the motifs used for packing, it's meant for troubleshooting.", "XRW TO DO") //should this be default or just optional?
		+ XMLSchemaAttribute::attribute_w_default( "mode", xs_string, "The option 'standard' will just skip the mode setting. 'standard' is the recommended use, it requires taskops and then looks for motifs in the designable residues. 'basic' uses Will's scheme logic to decide where the motifs are being placed (and therefore which residues are designable).", "basic" )
		+ XMLSchemaAttribute::attribute_w_default( "allowed_aas", xs_string, "This means that it will use the motifs in the motif files. Do not change this, there is currently no other option.", "motifs" )
		+ XMLSchemaAttribute::attribute_w_default( "halt_on_error", xsct_rosetta_bool, "XRW TO DO: unsure?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "motif_sets", xs_string, "Path to motif sets file to use.", "XRW TO DO" ) ; //XRW TO DO: does this look correct?

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist ) ;
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist ) ;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This mover will use only motif pairs found in input motif_sets or in the default motif sets to design interface residues.", attlist );
}

std::string SchemePlaceMotifsMoverCreator::keyname() const {
	return SchemePlaceMotifsMover::mover_name();
}

protocols::moves::MoverOP
SchemePlaceMotifsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SchemePlaceMotifsMover );
}

void SchemePlaceMotifsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SchemePlaceMotifsMover::provide_xml_schema( xsd );
}


} // matdes
} // protocols
