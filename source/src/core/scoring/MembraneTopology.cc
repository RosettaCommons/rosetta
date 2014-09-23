// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/MembraneTopology.cc
/// @brief  MembraneTopology
/// @author Bjorn Wallner
///


// Unit headers
#include <core/scoring/MembraneTopology.hh>


// Project headers
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>



// just for debugging
//#include <ObjexxFCL/format.hh>

// C++


namespace core {
namespace scoring {

static thread_local basic::Tracer TR( "core.scoring.MembraneTopology" );

MembraneTopology::MembraneTopology( MembraneTopology const & src ) :
	CacheableData()
{
	helix_id_=src.helix_id_;
	span_=src.span_;
	full_span_=src.full_span_;
	relative_tmh_ori_=src.relative_tmh_ori_;
	total_tmhelix_=src.total_tmhelix_;
	N_term_inside_=src.N_term_inside_;
	beta_barrel_=src.beta_barrel_;
	depth_=src.depth_;
	LipidExposure_=src.LipidExposure_;
	LipidBurial_=src.LipidBurial_;
	LipoDefined_=src.LipoDefined_;
	tmregion_=src.tmregion_;
	allow_scoring_=src.allow_scoring_;
	allow_tmh_scoring_=src.allow_tmh_scoring_;
	tmh_inserted_=src.tmh_inserted_;
	init_=src.init_;
  initialized_=src.initialized_;
  	  spanfile_=src.spanfile_;
}

std::string
MembraneTopology::read_in_spanfile()
{
    using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if(option[in::file::spanfile].user())
	{
		//At this point, assert if don't have a spanfile
		spanfile_ = option[in::file::spanfile].value();
		TR.Debug << "spanfile used by TMHTopologySamplerClaimer:  " << spanfile_ << std::endl;
		TR.Debug << "spanfile used by TMHTopologySamplerClaimer:  " << option[in::file::spanfile].value() << std::endl;
	}else{
		utility_exit_with_message( "[ERROR] Error opening spanfile '" + spanfile_ + "'" );
	}
	return spanfile_;
}

void
MembraneTopology::initialize(std::string const & spanfile)
{
    using namespace basic::options;
	using namespace basic::options::OptionKeys;
		//	std::string spanfile("BRD7.span");
	TR << "Initialize Membrane spanning regions with " << spanfile << std::endl;
	std::string line;
	utility::io::izstream stream (spanfile);
	getline(stream,line);//header for file. Usually says TM region prediction for...
	TR << line << std::endl;
	getline(stream,line);
	Size total_residue(0);
	{
		std::istringstream l(line);
		TR << line << std::endl;
		l>> total_tmhelix_ >> total_residue;
	}

	Size const max_tmhelix(total_tmhelix_);
	span_.dimension(max_tmhelix,2);
	full_span_.dimension(max_tmhelix,2);
	relative_tmh_ori_.dimension(max_tmhelix,max_tmhelix);
	helix_id_.dimension(max_tmhelix);


	getline(stream,line); //n2c
	TR << line << std::endl;
	getline(stream,line); //parallel
	TR << line << std::endl;
	for(Size i=1;i<=total_tmhelix_;++i)
	{
		getline(stream,line);
		{
			std::istringstream l(line);
			l >> span_(i,1) >> span_(i,2);
		}

		TR << line << std::endl ;
		helix_id_(i)=i;
	}
	if(total_tmhelix_==0)
	{
		utility_exit_with_message("bad format for spanfile total_tmhelix=0");
	}
	stream.close();
	stream.clear();



	for ( Size reg1 = 1; reg1 <= total_tmhelix_; ++reg1 ) {
		for ( Size reg2 = 1; reg2 <= total_tmhelix_; ++reg2 ) {
				relative_tmh_ori_(reg1,reg2)=1;
		}
	}

	full_span_(1,1)=1;
	full_span_(total_tmhelix_,2)=total_residue;

	for ( Size reg1 = 2; reg1 <= total_tmhelix_; ++reg1 ) {
		full_span_(reg1,1)=span_(reg1,1);
		full_span_(reg1-1,2)=span_(reg1,1);
	}
	depth_.resize(total_residue,0.0);
	tmregion_.resize(total_residue,false);
	allow_scoring_.resize(total_residue,true);
	allow_tmh_scoring_.resize(total_tmhelix_,true);
	tmh_inserted_=total_tmhelix_;

	for(Size i=1;i<=total_residue;++i) {
		for ( Size reg1 = 1; reg1 <= total_tmhelix_; ++reg1 ) {
			if(i>=span_(reg1,1) && i<=span_(reg1,2)) {
				tmregion_[i]=true;
				continue;
			}
		}
		TR << "tmregion " << i << " " << tmregion_[i] << std::endl;
	}

	//init Lipo
	if(option[in::file::lipofile].user())
	{
		Size NumberOfConstraints(0);
		Size resnum;
		Real exposure;

		std::string lipofile(option[OptionKeys::in::file::lipofile]());
		TR << "init lipo using " << lipofile << std::endl;
		LipidExposure_.resize(total_residue,0.0);
		LipidBurial_.resize(total_residue,0.0);
		stream.open(lipofile);

		if(stream) {
			getline(stream,line);
			getline(stream,line);
			while(!stream.eof()) {
				std::istringstream l(line);

				l >> resnum;
				l >> exposure;
				if(exposure>0) {
					LipidExposure_[resnum]+=exposure;
					NumberOfConstraints++;
					} else {
					LipidBurial_[resnum]+=std::abs(exposure);
					NumberOfConstraints++;
				}
				TR <<  resnum << " " << exposure << " " << LipidExposure_[resnum] << " " << LipidBurial_[resnum] << std::endl;
				getline(stream,line);
			}
			stream.close();
			stream.clear();
			TR << NumberOfConstraints << " exposure constraints read!" << std::endl;
			if(NumberOfConstraints>0)
			{
				LipoDefined_=true;
			}
		}else{
			TR << "unable to open " << lipofile << std::endl;
		}



	}



	init_=true;
  //pba
  initialized_=true;
	return;
}

void
MembraneTopology::shift_span( Size shift )
{

	for(Size i=1;i<=total_tmhelix_;++i)
		{
			span_(i,1)+=shift;
			span_(i,2)+=shift;
			full_span_(i,1)+=shift;
			full_span_(i,2)+=shift;
			if(full_span_(i,1)<=0)
			{
				utility_exit_with_message("Span out of bounds");
			}
		}
}

//void
//MembraneTopology::attach_to_pose(pose::Pose & pose)
//{
//	if ( pose.data().has( basic::MEMBRANE_TOPOLOGY ) ) {
//
//		core::scoring::MembraneTopologyOP topology = *this; //new core::scoring::MembraneTopology;
//		return *( static_cast< core::scoring::MembraneTopology * >( pose.data().get_ptr( basic::MEMBRANE_TOPOLOGY )() ));
//	}
	// else
//	core::scoring::MembraneTopologyOP topology = *this;
//	pose.data().set( basic::MEMBRANE_TOPOLOGY, *this); //topology );

//}

void
MembraneTopology::print() const
{
	TR << "TOTAL TMH " << total_tmhelix_ << std::endl;
	for(Size i=1;i<=total_tmhelix_;++i)
		{
			TR << "SPAN " << i << " " << span_(i,1) << " " << span_(i,2) << std::endl;
		}
}


void
MembraneTopology::get_subset( utility::vector1< Size > & TMH_list, MembraneTopology & src)
{

	// Assume list is sorted for now...
	// Will add a sorter here later....

	if(TMH_list.size()>src.total_tmhelix_)
		{
			utility_exit_with_message("Too long TMH_list");
		}
	Size const len(TMH_list.size());
	span_.dimension(len,2);
	full_span_.dimension(len,2);
	relative_tmh_ori_.dimension(len,len);


	for(Size i=1;i<=len;++i)
		{
			span_(i,1)=src.span_(TMH_list[i],1);
			span_(i,2)=src.span_(TMH_list[i],2);
			full_span_(i,1)=src.full_span_(TMH_list[i],1);
			full_span_(i,2)=src.full_span_(TMH_list[i],2);
			helix_id_(i)=TMH_list[i];
			for(Size j=1;j<=TMH_list.size();++j)
				{
					relative_tmh_ori_(i,j)=src.relative_tmh_ori_(TMH_list[i],TMH_list[j]);
				}
		}



	total_tmhelix_=TMH_list.size();




}
//pbadebug
MembraneTopology const &
MembraneTopology_from_pose( pose::Pose const & pose )
{
  // ////using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;

  assert( pose.data().has( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) );
  return *( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology const > ( pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
}

/// @details Either returns a non-const reference to the cenlist object already stored
/// in the pose, or creates a new cenlist object, places it in the pose, and returns
/// a non-const reference to it.
MembraneTopology &
nonconst_MembraneTopology_from_pose( pose::Pose & pose )
{
  // ////using core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY;

  if ( pose.data().has( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ) {
    return *( utility::pointer::static_pointer_cast< core::scoring::MembraneTopology > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) ));
  }
  // else
  MembraneTopologyOP membrane_topology( new MembraneTopology );
  pose.data().set( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY, membrane_topology );
  return *membrane_topology;
}
//pbadebug
/*
MembraneTopology &
MembraneTopology_from_pose( pose::Pose const & pose ) const
{
	return *( static_cast< MembraneTopology const * >( pose.data().get_const_ptr( basic::MEMBRANE_TOPOLOGY )() ));
}
*/
}
}
