// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/parse_filters.cc 
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Project Headers
#include <core/types.hh>
#include <utility/exit.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
//#include <protocols/moves/ResidueMover.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility Headers
#include <basic/Tracer.hh>

// Unit Headers

// C++ headers
#include <map>
#include <string>

#include <core/chemical/AA.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

using namespace core;
using namespace core::scoring;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.parse_filters" );

namespace protocols {
namespace protein_interface_design {
using namespace protocols::filters;

using namespace protocols::moves;

using namespace utility::tag;
using namespace std;

void
ResiduesInInterfaceFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	residues_in_interface_threshold_ = tag->getOption<core::Size>( "residues", 20 );
	rb_jump_ = tag->getOption<core::Size>( "jump_number", 1 );

	TR<<"residues in interface filter over jump number " << rb_jump_ << " with threshold "<<residues_in_interface_threshold_<<std::endl;
}

void
ScoreTypeFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	using namespace core::scoring;

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	score_type_ = core::scoring::score_type_from_name( tag->getOption<string>( "score_type", "total_score" ) );
	if( ! tag->hasOption( "threshold" ) ) throw utility::excn::EXCN_RosettaScriptsOption("Must specify 'threshold' for ScoreTypeFilter.");
	score_type_threshold_ = tag->getOption<core::Real>( "threshold" );

	TR<<"ScoreType filter for score_type "<<score_type_<<" with threshold "<<score_type_threshold_<<std::endl;
}

void
InterfaceSasaFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	lower_threshold_ = tag->getOption<core::Real>( "threshold", 800 );
	jump( tag->getOption< core::Size >( "jump", 1 ));
	hydrophobic_ = tag->getOption<bool>( "hydrophobic", false );
	polar_ = tag->getOption<bool>( "polar", false );
	runtime_assert( !hydrophobic_ || !polar_ );
	if( jump() != 1 && ( polar_ || hydrophobic_ ) )
		throw utility::excn::EXCN_RosettaScriptsOption( "ERROR: presently, only total sasa is supported across a jump other than 1. Remove polar and hydrophobic flags and try again." );

	TR<<"SasaFilter with lower threshold of "<<lower_threshold_<<" Ang^2 and jump "<<jump()<<'\n';
	if( hydrophobic_ )
		TR<<"Only reporting hydrophobic sasa\n";
	if( polar_ )
		TR<<"Only reporting polar sasa\n";
	TR.flush();
}

void
AlaScan::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	distance_threshold_ = tag->getOption<core::Real>( "interface_distance_cutoff", 8.0 );
	chain1_ = tag->getOption< bool >( "partner1", 0 );
	chain2_ = tag->getOption< bool >( "partner2", 1 );
	jump_ = tag->getOption< Size >( "jump", 1 );
	runtime_assert( chain1_ || chain2_ );
	repeats_ = tag->getOption< core::Size >( "repeats", 1 );
	symmetry_ = tag->getOption< bool >( "symmetry", 0 );
	repack( tag->getOption< bool >( "repack", 1 ) );

	if ( symmetry_ ) {
		using namespace core::scoring::symmetry;
		scorefxn_ = new SymmetricScoreFunction( * protocols::rosetta_scripts::parse_score_function( tag, data ) );
		TR<<"Symmetric AlaScan with distance threshold of "<<distance_threshold_<<" Ang "<<". jump="<<jump_<<" partner1="<<chain1_<<", partner2="<<chain2_<<" using "<<repeats_<<" repeats."<<std::endl;
		return;
	}
	using namespace core::scoring;
	scorefxn_ = new ScoreFunction( * protocols::rosetta_scripts::parse_score_function( tag, data ) );
	TR<<"AlaScan with distance threshold of "<<distance_threshold_<<" Ang "<<". jump="<<jump_<<" partner1="<<chain1_<<", partner2="<<chain2_<<" using "<<repeats_<<" repeats repack "<<repack()<<std::endl;
}

void
NeighborTypeFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
  residue_types_.assign( chemical::num_canonical_aas, false );
	utility::vector0< TagCOP > const & neighbor_type_tags( tag->getTags() );
	for( utility::vector0< TagCOP >::const_iterator nt_it=neighbor_type_tags.begin(); nt_it!=neighbor_type_tags.end(); ++nt_it ) {
    TagCOP const nt_tag_ptr = *nt_it;
    if( nt_tag_ptr->getName() == "Neighbor" ) {
			std::string const type( nt_tag_ptr->getOption<string>( "type" ) );
			residue_types_[ chemical::aa_from_name( type ) ] = true;
    }
	}
	target_residue_ = protocols::rosetta_scripts::get_resnum( tag, pose );
	distance_threshold_ = tag->getOption<core::Real>( "distance", 8.0 );

	TR<<"NeighborTypeFilter with distance threshold of "<<distance_threshold_<<" around residue "<<target_residue_<<std::endl;
}

void
ResidueBurialFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	target_residue_ = protocols::rosetta_scripts::get_resnum( tag, pose );
	distance_threshold_ = tag->getOption<core::Real>( "distance", 8.0 );
	neighbors_ = tag->getOption<core::Size>( "neighbors", 1 );

	TR<<"ResidueBurialFilter with distance threshold of "<<distance_threshold_<<" around residue "<<target_residue_<<" with "<<neighbors_<<" neighbors."<<std::endl;
}

void
ResidueDistanceFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	res1_ = protocols::rosetta_scripts::get_resnum( tag, pose, "res1_" );
	res2_ = protocols::rosetta_scripts::get_resnum( tag, pose, "res2_" );
	distance_threshold_ = tag->getOption<core::Real>( "distance", 8.0 );

	TR<<"ResidueDistanceFilter with distance threshold of "<<distance_threshold_<<" between residues "<<res1_<<" and "<<res2_<<std::endl;
}

void
DdgFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, Filters_map const & , Movers_map const & , core::pose::Pose const & )
{
	using namespace core::scoring;

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	ddg_threshold_ = tag->getOption<core::Real>( "threshold", -15 );
	rb_jump_ = tag->getOption< core::Size >( "jump", 1 );
	repeats( tag->getOption< core::Size >( "repeats", 1 ) );
	repack( tag->getOption< bool >( "repack", 1 ) );
	symmetry_ = tag->getOption<bool>( "symmetry", 0 );

	if( repeats() > 1 && !repack() )
		throw utility::excn::EXCN_RosettaScriptsOption( "ERROR: it doesn't make sense to have repeats if repack is false, since the values converge very well." );

	if ( symmetry_ )
		TR<<"ddg filter with threshold "<< ddg_threshold_<<" repeats="<<repeats()<<" and scorefxn "<<scorefxn_name<<" with symmetry " <<std::endl;
	else
		TR<<"ddg filter with threshold "<< ddg_threshold_<<" repeats="<<repeats()<<" and scorefxn "<<scorefxn_name<<" over jump "<<rb_jump_<<" and repack "<<repack()<<std::endl;
}

void
HbondsToResidueFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	partners_ = tag->getOption<core::Size>( "partners" );
	energy_cutoff_ = tag->getOption<core::Real>( "energy_cutoff", -0.5 );
	backbone_ = tag->getOption<bool>( "backbone", 0 );
	sidechain_ = tag->getOption<bool>( "sidechain", 1 );
	resnum_ = protocols::rosetta_scripts::get_resnum( tag, pose );

	TR<<"Hbonds to residue filter for resnum "<<resnum_<<" with "<<partners_<<" hbonding partners"<<std::endl;
}

void
HbondsToAtomFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & pose )
{
  partners_ = tag->getOption<core::Size>( "partners" );
  energy_cutoff_ = tag->getOption<core::Real>( "energy_cutoff", -0.5 );
  bb_bb_ = tag->getOption<bool>( "bb_bb", 0 );
  backbone_ = tag->getOption<bool>( "backbone", 0 );
  sidechain_ = tag->getOption<bool>( "sidechain", 1 );
  resnum_ = core::pose::get_resnum( tag, pose );

  if ( tag->hasOption( "atomname" ) ) {
    atomdesg_ = tag->getOption< std::string >( "atomname" );
  } else {
      throw utility::excn::EXCN_RosettaScriptsOption("Need to set atomname");
  }

  TR<<"Hbonds to atom filter for resnum "<<resnum_<<" and name " << atomdesg_ <<" with "<<partners_<<" hbonding partners"<<std::endl;
}

void
EnergyPerResidueFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	using namespace core::scoring;

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	score_type_ = core::scoring::score_type_from_name( tag->getOption<string>( "score_type", "total_score" ) );
	threshold_ = tag->getOption<core::Real>( "energy_cutoff", 0.0 );
	whole_interface_ = tag->getOption<bool>( "whole_interface" , 0 );
	rb_jump_ = tag->getOption<core::Size>( "jump_number", 1 );
	interface_distance_cutoff_ = tag->getOption<core::Real>( "interface_distance_cutoff" , 8.0 );

	if (whole_interface_==1 ) {
		resnum_ = 1;
		TR<<"energies for all interface residues with a distance cutoff of "
		<< interface_distance_cutoff_ << " A will be calculated \n"
		<< "jump_number is set to "<< rb_jump_
		<< "\n and scorefxn " <<scorefxn_name <<" will be used" <<std::endl;
	}
	else {
		resnum_ = protocols::rosetta_scripts::get_resnum( tag, pose );
		TR<<"EnergyPerResidueFilter for residue "<<resnum_<<" of score_type "<<score_type_<<" with cutoff "<<threshold_<<std::endl;
		}
}

void
BuriedUnsatHbondFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	jump_num_ = tag->getOption<core::Size>( "jump_number", 1 );
	upper_threshold_ = tag->getOption<core::Size>( "cutoff", 20 );

	TR<<"Buried Unsatisfied Hbond filter over jump number " << jump_num_ << " with cutoff " << upper_threshold_ << std::endl;
}

void
TerminusDistanceFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	jump_num_ = tag->getOption<core::Size>( "jump_number", 1 );
	distance_ = tag->getOption<core::Size>( "distance", 5 );

	TR<<"Distance From Terminus filter over jump number " << jump_num_ << " with cutoff " << distance_ << std::endl;
}


} // protein_interface_design
} // devel
