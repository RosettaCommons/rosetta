// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file
/// @brief
/// @author Jacob Bale ( balej@uw.edu )

// Unit headers
#include <protocols/matdes/SchemePlaceMotifsMover.hh>
#include <protocols/matdes/SchemePlaceMotifsMoverCreator.hh>
//#include <core/scoring/methods/FloatingPointCloudRestraintEnergyInfo.hh>
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>


#include <boost/foreach.hpp>

static basic::Tracer TR("protocols.matdes.SchemePlaceMotifsMover");


namespace protocols {
namespace matdes {

using namespace core;
using namespace utility;
using core::pose::Pose;


// -------------  Mover Creator -------------
std::string
SchemePlaceMotifsMoverCreator::keyname() const
{
	return SchemePlaceMotifsMoverCreator::mover_name();
}

protocols::moves::MoverOP
SchemePlaceMotifsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SchemePlaceMotifsMover );
}

std::string
SchemePlaceMotifsMoverCreator::mover_name()
{
	return "SchemePlaceMotifs";
}
// -------------  Mover Creator -------------

SchemePlaceMotifsMover::SchemePlaceMotifsMover() : halt_on_error_(false)
{
	core::scoring::motif::MotifHashManager::get_instance(); 
}

protocols::moves::MoverOP
SchemePlaceMotifsMover::clone() const {
	return protocols::moves::MoverOP( new SchemePlaceMotifsMover( *this ) );
}

protocols::moves::MoverOP
SchemePlaceMotifsMover::fresh_instance() const {
	return protocols::moves::MoverOP( new SchemePlaceMotifsMover() );
}


void
SchemePlaceMotifsMover::apply(Pose & pose) {

	////////////////////////////////////////////////////////////////////////////
	///////////////////////// get matching motifs //////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	core::scoring::motif::ResPairMotifQuery q(pose);
	core::scoring::motif::MotifHits hits;

	if( motif_sets_.size() ){
		TR << "using motif set names from XML" << std::endl;
		BOOST_FOREACH(std::string const & msetname, motif_sets_ ){
			core::scoring::motif::MotifHash const & mh( *core::scoring::motif::MotifHashManager::get_instance()->get_motif_hash_by_fname(msetname) );
			int n = mh.get_matching_motifs(q,hits);
			TR << n << " new hits from mset: " << msetname << ", total: " << hits.size() << std::endl;
		} 
	} else {
		TR << "using all motif sets from CLI" << std::endl;
		core::scoring::motif::MotifHashManager::get_instance()->get_matching_motifs(q,hits);
	}
	if( halt_on_error_ && hits.size()==0 ){
		utility_exit_with_message("no matching motifs found!");
	}

	std::cout << "found " << hits.size() << " matching motifs" << std::endl;
	if( dumpfile_!="" && hits.size() ){
		for(Size i = 0; i < 1; ++i) std::cout << "!!!!!!!!!!!!!!! DUMPING MOTIFS " << dumpfile_ << "!!!!!!!!!!!!!!!!" << std::endl;
		static int dump_count = 0;
		hits.dump_motifs_pdb(dumpfile_+"_"+ObjexxFCL::string_of(++dump_count)+".pdb");
	}

	////////////////////////////////////////////////////////////////////////////
	///////////////////////// remove interface sidechains //////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	for(core::Size ir = 1; ir <= pose.n_residue(); ++ir){
		if(pose.residue(ir).is_protein()){
			core::pose::replace_pose_residue_copying_existing_coordinates( pose, ir, pose.residue(ir).residue_type_set().name_map("ALA") );
		}
	}

	////////////////////////////////////////////////////////////////////////////
	//////////////////////// add cloud constraints /////////////////////////////
	////////////////////////////////////////////////////////////////////////////

	/*
	core::scoring::methods::FloatingPointCloudRestraintEnergyFuncType type = core::scoring::methods::FloatingPointCloudRestraintEnergyFuncType_HARMINIC;
	
	core::scoring::methods::FloatingPointCloudRestraintEnergyInfoOP fpcinfo( new core::scoring::methods::FloatingPointCloudRestraintEnergyInfo );
	pose.data().set( core::pose::datacache::CacheableDataType::FLOATING_POINT_CLOUD_INFO, fpcinfo );
	utility::vector1< core::scoring::methods::FloatingPoints > & clouds = fpcinfo->point_clouds();
	BOOST_FOREACH( core::scoring::motif::MotifHit const & mh, hits ){
		// numeric::xyzTransform<core::Real> x = numeric::random::gaussian_random_xform(360.0,20.0);
		core::scoring::methods::FloatingPoints cloud;
		BOOST_FOREACH(std::string atomname, utility::string_split("N CA C O CB",' ')){
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
	basic::datacache::DataMap & /*datamap*/,
	protocols::filters::Filters_map const & /*filtermap*/,
	protocols::moves::Movers_map const & /*movermap*/,
	core::pose::Pose const & /*pose*/
){
	tag_name_ = tag->getOption< std::string >("name","");
	dumpfile_ = tag->getOption<std::string>("dumpfile","");
	halt_on_error_ = tag->getOption<bool>("halt_on_error",false);
	std::string const motif_sets = tag->getOption< std::string >("motif_sets","");
	TR << "checking motif sets: " << motif_sets << std::endl;
	BOOST_FOREACH(std::string const & mset, utility::string_split(motif_sets,',')){
		if( mset=="" ) continue;
		if(!core::scoring::motif::MotifHashManager::get_instance()->have_motif_set_named(mset)){
			BOOST_FOREACH(std::string s,core::scoring::motif::MotifHashManager::get_instance()->motif_set_names()){
				std::cout << "candidate: " << s << std::endl;
			}
			utility_exit_with_message("MotifResiduesOperation: unknown motif set: "+mset+"\nspecify -mh:path:motifs");
		}
		motif_sets_.push_back(mset);
	}
	if(motif_sets_.size()==0){
		TR << "WARNING: no motif sets in XML " << tag_name_ << ", will use all from cli" << std::endl;
		// utility_exit_with_message("MotifResiduesOperation: no motif sets specified, use motif_sets='foo,bar,baz'!");
	}


}

} // matdes
} // protocols
