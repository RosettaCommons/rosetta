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
#include <protocols/moves/DataMap.hh>
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
	rb_jump_ ( rb_jump ),
	interface_distance_cutoff_ ( interface_distance_cutoff ),
	bb_bb_ ( bb_bb )
	{
		using namespace core::scoring;

		if( scorefxn ) scorefxn_ = new core::scoring::ScoreFunction( *scorefxn );
		if( score_type_ != total_score ) {
			core::Real const old_weight( scorefxn_->get_weight( score_type_ ) );
			scorefxn_->reset();
			scorefxn_->set_weight( score_type_, old_weight );

		}
	}
EnergyPerResidueFilter::EnergyPerResidueFilter( EnergyPerResidueFilter const &init ) :
	//utility::pointer::ReferenceCount(),
	Filter( init ), resnum_( init.resnum_ ),
	score_type_( init.score_type_ ),
	threshold_( init.threshold_ ),
	whole_interface_ (init.whole_interface_),
	rb_jump_ (init.rb_jump_),
	interface_distance_cutoff_ ( init.interface_distance_cutoff_),
	bb_bb_ ( init.bb_bb_ )
{
	using namespace core::scoring;
	if( init.scorefxn_ ) scorefxn_ = new core::scoring::ScoreFunction( *init.scorefxn_ );
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
EnergyPerResidueFilter::parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap & data, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
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
	
	if (whole_interface_) {
		resnum_ = 1;
		energy_per_residue_filter_tracer<<"energies for all interface residues with a distance cutoff of "
		<< interface_distance_cutoff_ << " A will be calculated \n"
		<< "jump_number is set to "<< rb_jump_
		<< "\n and scorefxn " << rosetta_scripts::get_score_function_name(tag) <<" will be used" <<std::endl;
	}

  if ( whole_protein_) {
		resnum_ = 1;
		energy_per_residue_filter_tracer<<"energies for all residues in protein" << std::endl;
	} 

	if (!whole_interface_ && !whole_protein_) {
		resnum_ = core::pose::get_resnum( tag, pose );
		energy_per_residue_filter_tracer<<"EnergyPerResidueFilter for residue "<<resnum_<<" of score_type "<<score_type_<<" with cutoff "<<threshold_<<std::endl;
	}
}

bool
EnergyPerResidueFilter::apply( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;
  using ObjexxFCL::FArray1D_bool;

	if ( whole_interface_)	{
		if ( pose.conformation().num_chains() < 2 ) {
			energy_per_residue_filter_tracer << "pose must contain at least two chains!" << std::endl;
			return false;

		}
		else {
			energy_per_residue_filter_tracer<<" \n \t --------- computing  ---------- \n \t-------- interface energies   --------\n \t \t------------------- \n" << std::endl;
		  core::pose::Pose in_pose = pose;
      FArray1D_bool partner1_( in_pose.total_residue(), false );
      in_pose.fold_tree().partition_by_jump( rb_jump_, partner1_);
      protocols::scoring::Interface interface_obj(rb_jump_);
      in_pose.update_residue_neighbors();
      interface_obj.distance( interface_distance_cutoff_ );
      interface_obj.calculate( in_pose );
      (*scorefxn_)( in_pose );
			
			bool pass;
			core::Real energy;
			for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
					if ( !in_pose.residue(resid).is_protein() ) continue;
					if( interface_obj.is_interface( resid ) ) {
						energy= compute( pose, resid ) ;
						if ( energy > threshold_ ) {
							pass=false;
							break;
						} else
						  pass=true;
					}//interface residues
				}
			return pass;
			}

	} 
   if ( whole_protein_ ) {
      energy_per_residue_filter_tracer<<" \n \t --------- computing  ---------- \n \t-------- all energies   --------\n \t \t------------------- \n" << std::endl;
			bool pass;
			core::Real energy;
		  core::pose::Pose in_pose = pose;
			for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
					if ( !in_pose.residue(resid).is_protein() ) continue;
					energy= compute( pose ,resid ) ;
					if ( energy > threshold_ ) {
							pass=false;
							break;
					} else
						  pass=true;
			}

      return pass;
	} 

	if (!whole_interface_ && !whole_protein_) {
		core::Real const energy( compute( pose, resnum_ ) );
		energy_per_residue_filter_tracer<<"Scoretype "<<name_from_score_type( score_type_ )<<" for residue: " << pose.pdb_info()->number(resnum_)<<" " <<pose.residue( resnum_).name3() <<" is "<<energy<<". ";
		bool const pass( energy <= threshold_ );
		if( pass ) energy_per_residue_filter_tracer<<"passing."<<std::endl;
		else energy_per_residue_filter_tracer<<"failing."<<std::endl;

		return pass;
  }
}


void
EnergyPerResidueFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	using namespace core::scoring;
	using ObjexxFCL::FArray1D_bool;

	if( whole_interface_ ) {
		core::pose::Pose in_pose = pose;
		FArray1D_bool partner1_( in_pose.total_residue(), false );
		in_pose.fold_tree().partition_by_jump( rb_jump_, partner1_);
		protocols::scoring::Interface interface_obj(rb_jump_);
		in_pose.update_residue_neighbors();
		interface_obj.distance( interface_distance_cutoff_ );
		interface_obj.calculate( in_pose );
		(*scorefxn_)( in_pose );

		out<<ObjexxFCL::fmt::A(9, "chain")<<
			ObjexxFCL::fmt::A( 9, "res")<<
			ObjexxFCL::fmt::A( 9, "AA")<<
			ObjexxFCL::fmt::A( 9, "total")<<
			ObjexxFCL::fmt::A( 9, "contact")<<
			ObjexxFCL::fmt::A( 9, "fa_atr")<<
			ObjexxFCL::fmt::A( 9, "fa_rep")<<
			ObjexxFCL::fmt::A( 9, "hb_bb_sc")<<
			ObjexxFCL::fmt::A( 9, "hb_sc")<<
			ObjexxFCL::fmt::A( 9, "fa_sol")<<
			ObjexxFCL::fmt::A( 9, "fa_dun")<<
			ObjexxFCL::fmt::A( 9, "fa_pair")<<"\n";
		for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
			if ( !in_pose.residue(resid).is_protein() ) continue;
			if( interface_obj.is_interface( resid ) ) { // in interface

				core::Real total=in_pose.energies().residue_total_energies( resid )[ ScoreType( total_score ) ];
				core::Real weighted_fa_atr=( (*scorefxn_)[ ScoreType( fa_atr) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( fa_atr ) ]);
				core::Real weighted_fa_rep=( (*scorefxn_)[ ScoreType( fa_rep) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( fa_rep) ]);
				core::Real weighted_hbond_bb_sc=( (*scorefxn_)[ ScoreType( hbond_bb_sc) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( hbond_bb_sc ) ]);
				core::Real weighted_hbond_sc=( (*scorefxn_)[ ScoreType( hbond_sc) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( hbond_sc ) ]);
				core::Real weighted_fa_sol=( (*scorefxn_)[ ScoreType( fa_sol) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( fa_sol ) ]);
				core::Real weighted_contact_score = weighted_fa_atr + weighted_fa_rep + weighted_hbond_bb_sc + weighted_hbond_sc + weighted_fa_sol;

				core::Real weighted_fa_dun=( (*scorefxn_)[ ScoreType( fa_dun) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( fa_dun ) ]);
				core::Real weighted_fa_pair=( (*scorefxn_)[ ScoreType( fa_pair ) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( fa_pair ) ]);

				out<<
					ObjexxFCL::fmt::A (9 , in_pose.pdb_info()->chain( resid) )<<
					ObjexxFCL::fmt::I(9,0, in_pose.pdb_info()->number(resid))<<
					ObjexxFCL::fmt::A (9,in_pose.residue( resid).name3())
					<<ObjexxFCL::fmt::F (9 , 3, total) <<" "
					<<ObjexxFCL::fmt::F (9 , 3, weighted_contact_score)<<" "
					<<ObjexxFCL::fmt::F (9 , 3, weighted_fa_atr) <<" "
					<<ObjexxFCL::fmt::F (9 , 3, weighted_fa_rep) <<" "
					<<ObjexxFCL::fmt::F (9 , 3, weighted_hbond_bb_sc )<<" "
					<<ObjexxFCL::fmt::F (9 , 3, weighted_hbond_sc )<<" "
					<<ObjexxFCL::fmt::F (9 , 3, weighted_fa_sol )<<" "
					<<ObjexxFCL::fmt::F (9 , 3, weighted_fa_dun )<<" "
					<<ObjexxFCL::fmt::F (9 , 3, weighted_fa_pair )<<"\n";

			}
		}
	} 
  if ( whole_protein_ ) {
    core::pose::Pose in_pose = pose;
    (*scorefxn_)( in_pose );
    out<<ObjexxFCL::fmt::A(9, "chain")<<
      ObjexxFCL::fmt::A( 9, "res")<<
      ObjexxFCL::fmt::A( 9, "AA")<<
      ObjexxFCL::fmt::A( 9, "total")<<
      ObjexxFCL::fmt::A( 9, "contact")<<
      ObjexxFCL::fmt::A( 9, "fa_atr")<<
      ObjexxFCL::fmt::A( 9, "fa_rep")<<
      ObjexxFCL::fmt::A( 9, "hb_bb_sc")<<
      ObjexxFCL::fmt::A( 9, "hb_sc")<<
      ObjexxFCL::fmt::A( 9, "fa_sol")<<
      ObjexxFCL::fmt::A( 9, "fa_dun")<<
      ObjexxFCL::fmt::A( 9, "fa_pair")<<"\n";
    for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
      if ( !in_pose.residue(resid).is_protein() ) continue;

        core::Real total=in_pose.energies().residue_total_energies( resid )[ ScoreType( total_score ) ];
        core::Real weighted_fa_atr=( (*scorefxn_)[ ScoreType( fa_atr) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( fa_atr ) ]);
        core::Real weighted_fa_rep=( (*scorefxn_)[ ScoreType( fa_rep) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( fa_rep) ]);
        core::Real weighted_hbond_bb_sc=( (*scorefxn_)[ ScoreType( hbond_bb_sc) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( hbond_bb_sc ) ]);
        core::Real weighted_hbond_sc=( (*scorefxn_)[ ScoreType( hbond_sc) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( hbond_sc ) ]);
        core::Real weighted_fa_sol=( (*scorefxn_)[ ScoreType( fa_sol) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( fa_sol ) ]);
        core::Real weighted_contact_score = weighted_fa_atr + weighted_fa_rep + weighted_hbond_bb_sc + weighted_hbond_sc + weighted_fa_sol;

        core::Real weighted_fa_dun=( (*scorefxn_)[ ScoreType( fa_dun) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( fa_dun ) ]);
        core::Real weighted_fa_pair=( (*scorefxn_)[ ScoreType( fa_pair ) ] ) * ( in_pose.energies().residue_total_energies( resid )[ ScoreType( fa_pair ) ]);

        out<<
          ObjexxFCL::fmt::A (9 , in_pose.pdb_info()->chain( resid) )<<
          ObjexxFCL::fmt::I(9,0, in_pose.pdb_info()->number(resid))<<
          ObjexxFCL::fmt::A (9,in_pose.residue( resid).name3())
          <<ObjexxFCL::fmt::F (9 , 3, total) <<" "
          <<ObjexxFCL::fmt::F (9 , 3, weighted_contact_score)<<" "
          <<ObjexxFCL::fmt::F (9 , 3, weighted_fa_atr) <<" "
          <<ObjexxFCL::fmt::F (9 , 3, weighted_fa_rep) <<" "
          <<ObjexxFCL::fmt::F (9 , 3, weighted_hbond_bb_sc )<<" "
          <<ObjexxFCL::fmt::F (9 , 3, weighted_hbond_sc )<<" "
          <<ObjexxFCL::fmt::F (9 , 3, weighted_fa_sol )<<" "
          <<ObjexxFCL::fmt::F (9 , 3, weighted_fa_dun )<<" "
          <<ObjexxFCL::fmt::F (9 , 3, weighted_fa_pair )<<"\n";
    }
	} 

	if (!whole_interface_ && !whole_protein_) {
		core::Real const energy( compute( pose, resnum_ ) );
		out<<"Scoretype "<<name_from_score_type( score_type_ )<<" for residue " << resnum_ <<" is "<<energy<<". ";
		bool const pass( energy <= threshold_ );
		if( pass ) out<<"passing."<<'\n';
		else out<<"failing."<<'\n';
  }
}


core::Real
EnergyPerResidueFilter::report_sm( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;
  using ObjexxFCL::FArray1D_bool;
  core::Real energy;

  if( whole_interface_ ) {

    core::pose::Pose in_pose = pose;
    FArray1D_bool partner1_( in_pose.total_residue(), false );
    in_pose.fold_tree().partition_by_jump( rb_jump_, partner1_);
    protocols::scoring::Interface interface_obj(rb_jump_);
    in_pose.update_residue_neighbors();
    interface_obj.distance( interface_distance_cutoff_ );
    interface_obj.calculate( in_pose );
    (*scorefxn_)( in_pose );

    core::Real ind_energy=0;
    core::Size num_res=0;
    for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
        if ( !in_pose.residue(resid).is_protein() ) continue;
					if( interface_obj.is_interface( resid ) ) {
        		ind_energy+= compute( pose, resid ) ;
        		num_res++;
					}
    }
    energy=ind_energy/num_res;


	} 
	if ( whole_protein_ ) {
      core::Real ind_energy=0;
			core::Size num_res=0;
      core::pose::Pose in_pose = pose;
      for ( core::Size resid = 1; resid <= pose.total_residue(); ++resid) {
          if ( !in_pose.residue(resid).is_protein() ) continue;
          ind_energy+= compute( pose, resid ) ;
				  num_res++;
      }
			energy=ind_energy/num_res;
	}
   
	if (!whole_interface_ && !whole_protein_) {
		energy=compute( pose, resnum_ );
  }

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
