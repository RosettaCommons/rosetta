// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueIEFilter.cc
/// @brief
/// @author Sagar Khare (khares@u.washington.edu)

#include <protocols/simple_filters/ResidueIEFilter.hh>
#include <protocols/simple_filters/ResidueIEFilterCreator.hh>

#include <protocols/filters/Filter.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>
#include <core/pose/selection.hh>
#include <protocols/simple_filters/EnergyPerResidueFilter.hh>
#define foreach BOOST_FOREACH

namespace protocols {
namespace simple_filters {

static basic::Tracer tr( "protocols.simple_filters.ResidueIEFilter" );

protocols::filters::FilterOP
ResidueIEFilterCreator::create_filter() const { return new ResidueIEFilter; }

std::string
ResidueIEFilterCreator::keyname() const { return "ResidueIE"; }

ResidueIEFilter::ResidueIEFilter( 
	utility::vector1< core::Size > const resnums,
	std::string const restype,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::scoring::ScoreType const score_type,
	core::Real const threshold,
	bool const whole_pose,
	bool const whole_interface,
	core::Size const rb_jump,
	core::Real const interface_distance_cutoff,
	core::Real max_penalty,
	core::Real penalty_factor
	) : 
	filters::Filter( "ResidueIE" ),
	resnums_( resnums ),
	restype_ (restype),
	score_type_( score_type ),
	threshold_( threshold ),
	whole_pose_ ( whole_pose ),
	whole_interface_ ( whole_interface ),
	rb_jump_ ( rb_jump ),
	interface_distance_cutoff_ ( interface_distance_cutoff ),
	max_penalty_ (max_penalty),
	penalty_factor_ (penalty_factor)
	{
		using namespace core::scoring;

		if( scorefxn ) scorefxn_ = new core::scoring::ScoreFunction( *scorefxn );
		if( score_type_ != total_score ) {
			core::Real const old_weight( scorefxn_->get_weight( score_type_ ) );
			scorefxn_->reset();
			scorefxn_->set_weight( score_type_, old_weight );

		}
		use_resE_ = false;
	}

ResidueIEFilter::ResidueIEFilter( ResidueIEFilter const &init ) :
	Filter( init ), resnums_( init.resnums_ ),
	restype_( init.restype_ ),
	score_type_( init.score_type_ ),
	threshold_( init.threshold_ ),
	whole_pose_ (init.whole_pose_),
	whole_interface_ (init.whole_interface_),
	rb_jump_ (init.rb_jump_),
	interface_distance_cutoff_ ( init.interface_distance_cutoff_),
	max_penalty_ ( init.max_penalty_),
	penalty_factor_ ( init.penalty_factor_),
	use_resE_ ( init.use_resE_ )
{
	using namespace core::scoring;
	if( init.scorefxn_ ) scorefxn_ = new core::scoring::ScoreFunction( *init.scorefxn_ );
}

core::scoring::ScoreFunctionOP
ResidueIEFilter::scorefxn() const{
	return scorefxn_;
}

void
ResidueIEFilter::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

core::scoring::ScoreType
ResidueIEFilter::score_type() const{
	return score_type_;
}

void
ResidueIEFilter::score_type( core::scoring::ScoreType score_type ){
	score_type_ = score_type;
}

core::Real
ResidueIEFilter::threshold() const{
	return threshold_;
}

void
ResidueIEFilter::threshold( core::Real const th ){
	threshold_ = th;
}

utility::vector1 <core::Size>
ResidueIEFilter::resnums() const{
	return resnums_;
}

void
ResidueIEFilter::resnums( utility::vector1<core::Size> const & rn ){
	resnums_ = rn;
}

ResidueIEFilter::~ResidueIEFilter() {}

void
ResidueIEFilter::parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap & data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	using namespace core::scoring;
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	score_type_ = core::scoring::score_type_from_name( tag->getOption<std::string>( "score_type", "total_score" ) );
	threshold_ = tag->getOption<core::Real>( "energy_cutoff", 0.0 );
	whole_pose_ = tag->getOption<bool>( "whole_pose" , 0 );
	whole_interface_ = tag->getOption<bool>( "interface" , 0 );
	rb_jump_ = tag->getOption<core::Size>( "jump_number", pose.num_jump() );
	interface_distance_cutoff_ = tag->getOption<core::Real>( "interface_distance_cutoff" , 8.0 );
	restype_ =  tag->getOption<std::string>( "restype3", "TRP" );
	penalty_factor_ =  tag->getOption<core::Real>( "penalty_factor", 1.0 );
	max_penalty_ =  tag->getOption<core::Real>( "max_penalty", 1000.0 );

	if(tag->hasOption("residues")) {
		tr << " the tag residues is seen by the program" << std::endl;
		resnums_ = core::pose::get_resnum_list(tag, "residues", pose);

		if(resnums_.empty()){
			whole_pose_=1;
      whole_interface_=0;
			tr<< "Failed to parse residues: " << tag->getOption<std::string> ("residues") << ". Using whole pose." << std::endl;
		}
    else
    {
      whole_pose_=0;
      whole_interface_=0;
    }
	}

	runtime_assert(tag->hasOption("residues") || whole_pose_ || whole_interface_);

	use_resE_ = tag->getOption<bool>( "use_resE" , 0 );
}

bool
ResidueIEFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real const penalty( compute( pose ) );
	if (penalty>max_penalty_) return false;
	return true;
}	



void
ResidueIEFilter::report( std::ostream & out, core::pose::Pose const & pose ) const 
{
		core::Real const penalty( compute( pose ) );
		out << "Total penalty for restype "<< restype_ << "is "<< penalty << std::endl;
}


core::Real
ResidueIEFilter::report_sm( core::pose::Pose const & pose ) const
{

	core::Real const penalty( compute( pose ) );
	return( penalty );
}



core::Real
ResidueIEFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;
	using namespace core::graph;


  if ( whole_interface_)
  {
		tr << "Detecting target resnums from interface." << std::endl;
		resnums_.clear();

    if ( pose.conformation().num_chains() < 2 ) {
      tr << "pose must contain at least two chains!" << std::endl;
      return false;
    }

    core::pose::Pose in_pose = pose;

    (*scorefxn_)(in_pose);
    protocols::scoring::Interface interface_obj(rb_jump_);

    in_pose.update_residue_neighbors();

    interface_obj.distance( interface_distance_cutoff_ );
    interface_obj.calculate( in_pose );
    for ( core::Size resnum = 1; resnum <= pose.total_residue(); ++resnum)
    {
      if ( in_pose.residue(resnum).is_protein()  &&  interface_obj.is_interface( resnum ) && (in_pose.residue_type(resnum).name3() == restype_) )
      {
        resnums_.push_back( resnum );
      }
    }
  }// whole_interface_

  else if ( whole_pose_ )
  {
		tr << "Detecting target resnums from whole pose." << std::endl;
		resnums_.clear();

    core::pose::Pose in_pose = pose;
    for ( core::Size resnum = 1; resnum <= in_pose.total_residue(); ++resnum) 
    {
      if ( in_pose.residue(resnum).is_protein()  && (in_pose.residue_type(resnum).name3() == restype_) ) resnums_.push_back( resnum );
    }
  }//whole_pose_

  std::unique( resnums_.begin(), resnums_.end() );
  tr << "The following residues will be considered for interaction energy calculation:"<< std::endl;
  foreach (core::Size const res, resnums_){
    tr << pose.residue_type( res ).name3() << res <<" + ";
    if (!(pose.residue_type( res ).name3() == restype_)) {
      tr << "Residue " << res << " in pose is of type "<< pose.residue_type( res ).name3() << ". Requested restype3 is "<< restype_<<". Skipping!"<<std::endl;
      return false;
    }
  }

	core::pose::Pose in_pose = pose;
	(*scorefxn_)(in_pose);
	core::Real penalty (0.0);
  EnergyMap const curr_weights = in_pose.energies().weights();
  
	if (resnums_.size() == 0) {
    tr << "No residues found. Skipping calculation."<< std::endl;

    return (0.0);
  }
	
  foreach (core::Size const res, resnums_)
  {


		core::Real res_intE (0.0); 

		if (use_resE_)
    {
			protocols::simple_filters::EnergyPerResidueFilter const eprf(res, scorefxn_, score_type_, 100000.0/*dummy threshold*/);
			res_intE = eprf.compute( pose );
		}
		else
    {     
      core::Real res_intE_samechain (0.0);
      core::Real res_intE_differentchain (0.0);

      int res_chain = in_pose.chain(res);

      //Fill residue energies by traversing energy graph
			for ( EdgeListConstIterator egraph_it = in_pose.energies().energy_graph().get_node( res )->const_edge_list_begin(); egraph_it != in_pose.energies().energy_graph().get_node( res )->const_edge_list_end(); ++egraph_it)
      {
				EnergyEdge const * Eedge = static_cast< EnergyEdge const * > (*egraph_it);
				res_intE += Eedge->dot( curr_weights );



        if (in_pose.chain(Eedge->get_other_ind(res)) == res_chain)
        {
          res_intE_samechain += Eedge->dot( curr_weights );
        }
        else
        {
          res_intE_differentchain += Eedge->dot( curr_weights );
        }
			}//for each egraph_it
      tr << "Residue "<< pose.residue_type( res ).name3()<<res<< " has a intra-chain interaction energy "<< res_intE_samechain << std::endl;
      tr << "Residue "<< pose.residue_type( res ).name3()<<res<< " has a inter-chain interaction energy "<< res_intE_differentchain << std::endl;
		}

		tr << "Residue "<< pose.residue_type( res ).name3()<<res<< " has an (interaction) energy "<< res_intE <<", threshold is "<< threshold_<<" and penalty is ";

		if (res_intE > threshold_) 
    {
			penalty+= (res_intE - threshold_);
			tr<< (res_intE - threshold_) << std::endl;
		}
		else 
    {
      tr << " 0"<<std::endl;
    }
	} //foreach res
	
	return( penalty * penalty_factor_ );
}

}
}
