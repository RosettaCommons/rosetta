// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/EnergyerResidueFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#include <protocols/simple_filters/EnergyPerResidueFilter.hh>
#include <protocols/simple_filters/EnergyPerResidueFilterCreator.hh>

#include <protocols/filters/Filter.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <protocols/scoring/Interface.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/Energies.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <set>

namespace protocols {
namespace simple_filters {

static basic::Tracer energy_per_residue_filter_tracer( "protocols.simple_filters.EnergyPerResidueFilter" );

protocols::filters::FilterOP
EnergyPerResidueFilterCreator::create_filter() const { return new EnergyPerResidueFilter; }

std::string
EnergyPerResidueFilterCreator::keyname() const { return "EnergyPerResidue"; }

EnergyPerResidueFilter::EnergyPerResidueFilter(
	core::Size const resnum,
	core::scoring::ScoreFunctionCOP scorefxn,
	core::scoring::ScoreType const score_type,
	core::Real const threshold,
	bool const whole_interface,
	bool const whole_protein,
	bool const select_resnums,
	bool const select_around_resnums,
	core::Real const around_shell,
  std::string const string_resnums,
  std::string const string_around_resnums,
	core::Size const rb_jump,
	core::Real const interface_distance_cutoff,
	bool const bb_bb
	) :
	filters::Filter( "EnergyPerResidue" ),
	resnum_( resnum ),
	score_type_( score_type ),
	threshold_( threshold ),
	whole_interface_ ( whole_interface ),
	whole_protein_ ( whole_protein ),
	select_resnums_ ( select_resnums ),
	select_around_resnums_ ( select_around_resnums ),
	string_resnums_ (string_resnums ),
	string_around_resnums_ (string_around_resnums ),
	around_shell_ ( around_shell ),
	rb_jump_ ( rb_jump ),
	interface_distance_cutoff_ ( interface_distance_cutoff ),
	bb_bb_ ( bb_bb )
	{
		using namespace core::scoring;

		if( scorefxn ) scorefxn_ = scorefxn->clone();
		if( score_type_ != total_score ) {
			core::Real const old_weight( scorefxn_->get_weight( score_type_ ) );
			scorefxn_->reset();
			scorefxn_->set_weight( score_type_, old_weight );
	}
}


EnergyPerResidueFilter::EnergyPerResidueFilter( EnergyPerResidueFilter const &init ) :
	Filter( init ), resnum_( init.resnum_ ),
	score_type_( init.score_type_ ),
	threshold_( init.threshold_ ),
	whole_interface_ (init.whole_interface_),
	whole_protein_ (init.whole_protein_),
	select_resnums_ (init.select_resnums_),
	select_around_resnums_ (init.select_around_resnums_),
	string_resnums_ (init.string_resnums_ ),
	string_around_resnums_ (init.string_around_resnums_ ),
	around_shell_(init.around_shell_),
	rb_jump_ (init.rb_jump_),
	interface_distance_cutoff_ ( init.interface_distance_cutoff_),
	bb_bb_ ( init.bb_bb_ )
{
	using namespace core::scoring;
	if( init.scorefxn_ ) scorefxn_ = init.scorefxn_->clone();
}

core::scoring::ScoreFunctionOP
EnergyPerResidueFilter::scorefxn() const{
	return scorefxn_;
}

void
EnergyPerResidueFilter::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
	scorefxn_ = scorefxn;
}

core::scoring::ScoreType
EnergyPerResidueFilter::score_type() const{
	return score_type_;
}

void
EnergyPerResidueFilter::score_type( core::scoring::ScoreType score_type ){
	score_type_ = score_type;
}

core::Real
EnergyPerResidueFilter::threshold() const{
	return threshold_;
}

void
EnergyPerResidueFilter::threshold( core::Real const th ){
	threshold_ = th;
}

bool
EnergyPerResidueFilter::bb_bb() const{
	return bb_bb_;
}

void
EnergyPerResidueFilter::bb_bb( bool const b_b ){
	bb_bb_ = b_b;
}

core::Size
EnergyPerResidueFilter::resnum() const{
	return resnum_;
}

void
EnergyPerResidueFilter::resnum( core::Size const rn ){
	resnum_ = rn;
}

EnergyPerResidueFilter::~EnergyPerResidueFilter() {}

void
EnergyPerResidueFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
	using namespace core::scoring;

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	score_type_ = core::scoring::score_type_from_name( tag->getOption<std::string>( "score_type", "total_score" ) );
	threshold_ = tag->getOption<core::Real>( "energy_cutoff", 0.0 );
	whole_interface_ = tag->getOption<bool>( "whole_interface" , 0 );
	whole_protein_ = tag->getOption<bool>( "whole_protein" , 0 );
	rb_jump_ = tag->getOption<core::Size>( "jump_number", 1 );
	interface_distance_cutoff_ = tag->getOption<core::Real>( "interface_distance_cutoff" , 8.0 );
	bb_bb_ = tag->getOption< bool >("bb_bb", false );
	resnum_=0;

	if (!whole_protein_) {
		select_resnums_=false;
		string_resnums_ = "";
		select_around_resnums_=false;
		string_around_resnums_ ="";
		if (whole_interface_) {
			if (pose.conformation().num_chains() < 2 ) {
				throw utility::excn::EXCN_RosettaScriptsOption( "ERROR: it doesn't make sense for interface since you have only one chain");
		    } else {
				energy_per_residue_filter_tracer<<"energies for all interface residues with a distance cutoff of "
				<< interface_distance_cutoff_ << " A will be calculated "
				<< "jump_number is set to "<< rb_jump_
				<< " and scorefxn " << rosetta_scripts::get_score_function_name(tag) <<" will be used" <<std::endl;
			  }
		} else if ( tag->hasOption( "resnums" ) ) {
			select_resnums_=true;
			string_resnums_ = tag->getOption< std::string >( "resnums" );
			energy_per_residue_filter_tracer<<"filter on energy per residue for" << string_resnums_ << std::endl;
		} else if ( tag->hasOption( "around_resnums" ) ) {
			select_around_resnums_=true;
			string_around_resnums_ = tag->getOption< std::string >( "around_resnums" );
		  around_shell_ = tag->getOption<core::Real>( "around_shell", 8.0 );
			energy_per_residue_filter_tracer<<"energies for residues around specified residues " << string_around_resnums_  <<" using a cutoff distance (around_shell) of " << around_shell_ << std::endl;
		} else {
     resnum_ = core::pose::get_resnum( tag, pose );
     energy_per_residue_filter_tracer<<"EnergyPerResidueFilter for residue "<<resnum_<<" of score_type "<<score_type_<<" with cutoff "<<threshold_<<std::endl;
		}
	} else {
		 whole_protein_=1;
     energy_per_residue_filter_tracer<<"EnergyPerResidueFilter default for whole_protein"<<std::endl;
  }
}

bool
EnergyPerResidueFilter::apply( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;
	using ObjexxFCL::FArray1D_bool;

  utility::vector1<bool> use_all_residues;
  use_all_residues.resize(pose.total_residue(),false);

	if ( whole_interface_ ) {
	energy_per_residue_filter_tracer << "apply reassign interface residues" << std::endl;
  core::pose::Pose in_pose = pose;
  FArray1D_bool partner1_( in_pose.total_residue(), false );
  in_pose.fold_tree().partition_by_jump( rb_jump_, partner1_);
  protocols::scoring::Interface interface_obj(rb_jump_);
  in_pose.update_residue_neighbors();
  interface_obj.distance( interface_distance_cutoff_ );
  interface_obj.calculate( in_pose );
  (*scorefxn_)( in_pose );

    for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
 	 	if( interface_obj.is_interface( resid ) ) {
 	 		use_all_residues[resid]=true;
 	   }	
 	 }
	} else if ( select_resnums_ ) {
    std::set< core::Size > const res_vec( core::pose::get_resnum_list( string_resnums_, pose ) );
    foreach( core::Size const res, res_vec)
					use_all_residues[ res ]=true;
  } else if ( select_around_resnums_ ) {
  if (string_around_resnums_.size()!=0) {
    std::set< core::Size > const res_vec( core::pose::get_resnum_list( string_around_resnums_, pose ) );
    foreach( core::Size const res, res_vec) {
					use_all_residues[ res ]=true;

		for ( core::Size i = 1; i <= pose.total_residue(); ++i) {
        if ( i == res ) continue; //dont need to check if is nbr of self
				if ( use_all_residues[ i ] ) continue; //dont need to check one that is marked already
        core::Real const distance( pose.residue( i ).xyz( pose.residue( i ).nbr_atom() ).distance( pose.residue( res ).xyz( pose.residue( res ).nbr_atom() )) );
				if( distance <= around_shell_ )
						use_all_residues[ i ]=true;
   }
	}
  } else {
     utility_exit_with_message( "around_resnums is true but no residues were specified" );
  }
	} else if ( resnum_ !=0 )  {
						use_all_residues[ resnum_ ]=true;
  } else {
		energy_per_residue_filter_tracer<< "apply reassign everything to be true" << std::endl;
		std::fill (use_all_residues.begin(),use_all_residues.end(),true);
  }

  //tracer output
  energy_per_residue_filter_tracer<<"Apply Filter Scoretype "<<name_from_score_type( score_type_ )<<" total of " << use_all_residues.size()<<": ";
  for ( core::Size i = 1; i <= use_all_residues.size()-1; ++i) {
      if (use_all_residues[i]) {
            energy_per_residue_filter_tracer << i << "+";
      }
  }

  if (use_all_residues[use_all_residues.size()])
            energy_per_residue_filter_tracer << use_all_residues.size()<<std::endl;
  else
            energy_per_residue_filter_tracer << std::endl;

	bool pass = false;
	core::Real energy=0;

  for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
					if ( !pose.residue(resid).is_protein() ) continue;
					if (use_all_residues[resid]) {
						energy= compute( pose, resid ) ;
						energy_per_residue_filter_tracer<<"Scoretype "<<name_from_score_type( score_type_ )<<" for residue: " << pose.pdb_info()->number(resid)<<" " <<pose.residue( resid ).name3() <<" is "<<energy<<std::endl;

						if ( energy > threshold_ ) {
							pass=false;
							break;
						} else
						  pass=true;
  			}
	}
  return pass;

}


void
EnergyPerResidueFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	using namespace core::scoring;
	using ObjexxFCL::FArray1D_bool;

  utility::vector1<bool> use_all_residues;
  use_all_residues.resize(pose.total_residue(),false);

  if ( whole_interface_ ) {
	energy_per_residue_filter_tracer << "report reassign interface residues" << std::endl;
  core::pose::Pose in_pose = pose;
  FArray1D_bool partner1_( in_pose.total_residue(), false );
  in_pose.fold_tree().partition_by_jump( rb_jump_, partner1_);
  protocols::scoring::Interface interface_obj(rb_jump_);
  in_pose.update_residue_neighbors();
  interface_obj.distance( interface_distance_cutoff_ );
  interface_obj.calculate( in_pose );
  (*scorefxn_)( in_pose );

    for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
    if( interface_obj.is_interface( resid ) ) {
      use_all_residues[resid]=true;
     }
   }
  } else if ( select_resnums_ ) {
    std::set< core::Size > const res_vec( core::pose::get_resnum_list( string_resnums_, pose ) );
    foreach( core::Size const res, res_vec)
          use_all_residues[ res ]=true;
  } else if ( select_around_resnums_ ) {
    if (string_around_resnums_.size()!=0) {
    std::set< core::Size > const res_vec( core::pose::get_resnum_list( string_around_resnums_, pose ) );
    foreach( core::Size const res, res_vec) {
          use_all_residues[ res ]=true;

    for ( core::Size i = 1; i <= pose.total_residue(); ++i) {
        if ( i == res ) continue; //dont need to check if is nbr of self
        if ( use_all_residues[ i ] ) continue; //dont need to check one that is marked already
        core::Real const distance( pose.residue( i ).xyz( pose.residue( i ).nbr_atom() ).distance( pose.residue( res ).xyz( pose.residue( res ).nbr_atom() )) );
        if( distance <= around_shell_ ) {
            use_all_residues[ i ]=true;
				}
   }
  }
  } else {
     utility_exit_with_message( "around_resnums is true but no residues were specified" );
  }
  } else if ( resnum_ !=0 )  {
            use_all_residues[ resnum_ ]=true;
  } else {
		energy_per_residue_filter_tracer<< "report reassign everything to be true" << std::endl;
		std::fill (use_all_residues.begin(),use_all_residues.end(),true);
  }

  //tracer output
  energy_per_residue_filter_tracer<<"Apply Filter Scoretype "<<name_from_score_type( score_type_ )<<" total of " << use_all_residues.size()<<": ";
  for ( core::Size i = 1; i <= use_all_residues.size()-1; ++i) {
      if (use_all_residues[i]) {
						energy_per_residue_filter_tracer << i << "+";
			}
  }
  if (use_all_residues[use_all_residues.size()])
            energy_per_residue_filter_tracer << use_all_residues.size()<<std::endl;
  else
            energy_per_residue_filter_tracer << std::endl;

  core::Real energy=0;
  for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
          if ( !pose.residue(resid).is_protein() ) continue;
          if (use_all_residues[resid])
            energy= compute( pose, resid ) ;
            out<<"Scoretype "<<name_from_score_type( score_type_ )<<" for residue: " << pose.pdb_info()->number(resid)<<" " <<pose.residue( resid ).name3() <<" is "<<energy<<std::endl;
  }

}

core::Real
EnergyPerResidueFilter::report_sm( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;
  using ObjexxFCL::FArray1D_bool;

  utility::vector1<bool> use_all_residues;
  use_all_residues.resize(pose.total_residue(),false);

  if ( whole_interface_ ) {
	energy_per_residue_filter_tracer << "report_sm reassign interface residues" << std::endl;
  core::pose::Pose in_pose = pose;
  FArray1D_bool partner1_( in_pose.total_residue(), false );
  in_pose.fold_tree().partition_by_jump( rb_jump_, partner1_);
  protocols::scoring::Interface interface_obj(rb_jump_);
  in_pose.update_residue_neighbors();
  interface_obj.distance( interface_distance_cutoff_ );
  interface_obj.calculate( in_pose );
  (*scorefxn_)( in_pose );

    for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
    if( interface_obj.is_interface( resid ) ) {
      use_all_residues[resid]=true;
     }
   }
  } else if ( select_resnums_ ) {
    std::set< core::Size > const res_vec( core::pose::get_resnum_list( string_resnums_, pose ) );
    foreach( core::Size const res, res_vec)
          use_all_residues[ res ]=true;
  } else if ( select_around_resnums_ ) {
		if (string_around_resnums_.size()!=0) {
        std::set< core::Size > const res_vec( core::pose::get_resnum_list( string_around_resnums_, pose ) );
        foreach( core::Size const res, res_vec) {
              use_all_residues[ res ]=true;
  
        for ( core::Size i = 1; i <= pose.total_residue(); ++i) {
            if ( i == res ) continue; //dont need to check if is nbr of self
            if ( use_all_residues[ i ] ) continue; //dont need to check one that is marked already
            core::Real const distance( pose.residue( i ).xyz( pose.residue( i ).nbr_atom() ).distance( pose.residue( res ).xyz( pose.residue( res ).nbr_atom() )) );
            if( distance <= around_shell_ ){
                use_all_residues[ i ]=true;
					}
       }
      }
 		 } else {
 		 		utility_exit_with_message( "around_resnums is true but no residues were specified" );
 		 }
  } else if ( resnum_ !=0 )  {
            use_all_residues[ resnum_ ]=true;
  } else {
		energy_per_residue_filter_tracer<< "report_sm reassign everything to be true" << std::endl;
		std::fill (use_all_residues.begin(),use_all_residues.end(),true);
  }

  //tracer output
  energy_per_residue_filter_tracer<<"Apply Filter Scoretype "<<name_from_score_type( score_type_ )<<" total of " << use_all_residues.size()<<": ";
  for ( core::Size i = 1; i <= use_all_residues.size()-1; ++i) {
      if (use_all_residues[i]) {
            energy_per_residue_filter_tracer << i << "+";
      }
  }
  if (use_all_residues[use_all_residues.size()])
            energy_per_residue_filter_tracer << use_all_residues.size()<<std::endl;
  else
            energy_per_residue_filter_tracer << std::endl;

  core::Real energy=999999;
  core::Real ind_energy=0;
  core::Size num_res=0;
  for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
          if ( !pose.residue(resid).is_protein() ) continue;
          if (use_all_residues[resid]) {
        		ind_energy+= compute( pose, resid ) ;
        		num_res++;
           }
  }
  energy=ind_energy/num_res;

	return( energy );
}

core::Real
EnergyPerResidueFilter::compute( core::pose::Pose const & pose , core::Size const resid ) const
{
	using namespace core::scoring;

	core::pose::Pose in_pose = pose;
	if( ( (*scorefxn_)[ interchain_env ] > 0.0 )
		&& ( (*scorefxn_)[ interchain_vdw ] > 0.0 )
		&& ( (*scorefxn_)[fa_rep] == 0.0 )
		&& ( (*scorefxn_)[fa_atr] == 0.0 ) )
		{
			if( in_pose.is_fullatom() ) {
			core::util::switch_to_residue_type_set( in_pose, core::chemical::CENTROID );
		}
	}
	else {
		if( in_pose.is_centroid() ) {
			core::util::switch_to_residue_type_set( in_pose, core::chemical::FA_STANDARD );
		}
	}

	in_pose.update_residue_neighbors();
	(*scorefxn_)( in_pose );
	core::Real weighted_score;
	if( score_type_ == total_score ) weighted_score = in_pose.energies().residue_total_energies( resid )[ ScoreType( score_type_ )];
	else {

		if( bb_bb_ ){
			energy_per_residue_filter_tracer << "decomposing bb hydrogen bond terms" << std::endl;
    	core::scoring::methods::EnergyMethodOptionsOP energy_options(new core::scoring::methods::EnergyMethodOptions(scorefxn_->energy_method_options()));
    	energy_options->hbond_options().decompose_bb_hb_into_pair_energies(true);
    	scorefxn_->set_energy_method_options(*energy_options);
	    (*scorefxn_)( in_pose );
		}

		core::Real const weight( (*scorefxn_)[ ScoreType( score_type_ ) ] );
		core::Real const score( in_pose.energies().residue_total_energies( resid )[ ScoreType( score_type_ ) ]);
		weighted_score = weight * score ;
	}
	return( weighted_score );
}

core::Real
EnergyPerResidueFilter::compute( core::pose::Pose const & pose ) const
{
  using namespace core::scoring;

  core::pose::Pose in_pose = pose;
  if( ( (*scorefxn_)[ interchain_env ] > 0.0 )
    && ( (*scorefxn_)[ interchain_vdw ] > 0.0 )
    && ( (*scorefxn_)[fa_rep] == 0.0 )
    && ( (*scorefxn_)[fa_atr] == 0.0 ) )
    {
      if( in_pose.is_fullatom() ) {
      core::util::switch_to_residue_type_set( in_pose, core::chemical::CENTROID );
    }
  }
  else {
    if( in_pose.is_centroid() ) {
      core::util::switch_to_residue_type_set( in_pose, core::chemical::FA_STANDARD );
    }
  }

  in_pose.update_residue_neighbors();
  (*scorefxn_)( in_pose );
  core::Real weighted_score;
  if( score_type_ == total_score ) weighted_score = in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( score_type_ )];
  else {

    if( bb_bb_ ){
      energy_per_residue_filter_tracer << "decomposing bb hydrogen bond terms" << std::endl;
      core::scoring::methods::EnergyMethodOptionsOP energy_options(new core::scoring::methods::EnergyMethodOptions(scorefxn_->energy_method_options()));
      energy_options->hbond_options().decompose_bb_hb_into_pair_energies(true);
      scorefxn_->set_energy_method_options(*energy_options);
    }

    core::Real const weight( (*scorefxn_)[ ScoreType( score_type_ ) ] );
    core::Real const score( in_pose.energies().residue_total_energies( resnum_ )[ ScoreType( score_type_ ) ]);
    weighted_score = weight * score ;
  }
  return( weighted_score );
}

}
}
