// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief app for matching motifs to a protein DNA interface
/// @author xingjiepan

#include <devel/init.hh>

// Project Headers (protocols)
#include <protocols/dna/DnaDesignDef.hh>
#include <protocols/dna/DnaInterfaceFinder.hh>
#include <protocols/dna/util.hh>
#include <protocols/motifs/motif_utils.hh>
#include <protocols/motifs/MotifSearch.hh>

// Project Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.dna_motif_matcher" );

// C++ Headers
#include <string>
#include <iostream>

// Option Key Includes
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/vector0.hh>





////////////////////////////////////////////////////////////////////////////////

void dna_motif_matcher()
{
  // 
  utility::vector1< std::string > pdb_files( basic::options::start_files() );
  // Set up scoring
  core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
  
  for ( utility::vector1< std::string >::const_iterator pdb_file( pdb_files.begin() );
        pdb_file != pdb_files.end(); ++pdb_file ) 
  {
    std::string pdb_name( *pdb_file );
    TR << "Working on file: " << pdb_name << std::endl;
    std::string pdb_prefix( utility::string_split( utility::string_split( pdb_name, '/' ).back(), '.' ).front() );
    if ( basic::options::option[ basic::options::OptionKeys::out::prefix ].user() ) 
    {
      pdb_prefix = basic::options::option[ basic::options::OptionKeys::out::prefix ]();
    }

    // Create the pose object from pdb file
    core::pose::PoseOP pose( new core::pose::Pose );
    core::import_pose::pose_from_file( *pose, pdb_name , core::import_pose::PDB_file);
    core::scoring::dna::set_base_partner( *pose );//Copied from sthyme's code. She said that she added this for bug fixing. This function adds to the pose a BasePartner object which is a vector of length #residues that records which base is paired to which.
    
    // Get protein positions
    utility::vector1< core::Size > protein_positions; // Protein positions
    core::Size nres( pose->size() ); // Number of residues
    for ( core::Size p_index(1); p_index<=nres; ++p_index ) 
    {
      if ( pose->residue_type( p_index ).is_protein() ) 
      {
        protein_positions.push_back( p_index );
      }
    }
    
    // Get target DNAs
    protocols::dna::DnaDesignDefOPs targeted_dna;
    if ( basic::options::option[ basic::options::OptionKeys::dna::design::dna_defs ].user() ) 
    {
      protocols::dna::load_dna_design_defs_from_options( targeted_dna );
      TR << "Designed DNAs are " << targeted_dna << std::endl;
    }
    if ( targeted_dna.empty() )
    {
      throw "Targeted DNA is not defined!";
    }
    
    // Save target DNA pairs into vectors
    utility::vector1< core::Size > dna_design_positions; // Target DNA positions
    utility::vector1< core::Size > base_partners; // Target DNA pairs
    dna_design_positions = protocols::motifs::defs2vector( *pose, targeted_dna );
    core::scoring::dna::BasePartner const & partner( core::scoring::dna::retrieve_base_partner_from_pose( *pose ) ); // Partners of target DNAs
    for ( core::Size i(1); i<=dna_design_positions.size(); ++i ) 
    {
      base_partners.push_back( dna_design_positions[i] );
      base_partners.push_back( partner[dna_design_positions[i]] );
    }
    
    // For each base type
    utility::vector1< std::string > base_types;
    base_types.push_back("GUA");
    base_types.push_back("ADE");
    base_types.push_back("CYT");
    base_types.push_back("THY");
    for( utility::vector1< std::string >::const_iterator t_base( base_types.begin() ),
          t_base_end( base_types.end() ); t_base != t_base_end; ++t_base)
    {
      // Make corresponding DNA mutations
      for( utility::vector1< core::Size >::const_iterator p_dna( dna_design_positions.begin() ) , 
            p_dna_end( dna_design_positions.end() ); p_dna != p_dna_end; ++p_dna)
      {
        protocols::motifs::make_base_pair_mutation( *pose, *p_dna, core::chemical::aa_from_name(*t_base) );
      }
      //pose->dump_pdb(std::cout, dna_design_positions);
      
      // Get DNA-protein interface
      protocols::dna::DnaInterfaceFinderOP interface( new protocols::dna::DnaInterfaceFinder( 10 * 10, 3.9 * 3.9, 6., true ) );
      interface->determine_protein_interface( *pose, protein_positions, base_partners );
      
      // Get design positions
      utility::vector1< core::Size > design_positions; //Design positions, aka, protein residues in the interface
      TR << "The protein design positions (not pdb indices) are: ";
      for(utility::vector1< core::Size >::const_iterator p_index( protein_positions.begin() ),
          end( protein_positions.end() ); p_index != end; ++p_index)
      {
        if( interface->protein_neighbors().at(*p_index).contact() )
        {
            design_positions.push_back(*p_index);
            TR << *p_index <<" ";
        }
      }
      TR << std::endl;
      //pose->dump_pdb(std::cout, design_positions);
    
      // Motif search on protein-DNA interface
      protocols::motifs::MotifSearchOP motif_search( new protocols::motifs::MotifSearch );
      motif_search->run( *pose, design_positions );
      
    }// loop over 4 base types
  }// loop over all pdb files

  TR << "SUCCESSFUL COMPLETION" << std::endl;
}

////////////////////////////////////////////////////////////////////////////////

int main( int argc, char * argv [] )
{

  try 
  {
    devel::init( argc, argv );
    dna_motif_matcher();

  } catch ( utility::excn::EXCN_Base const & e ) {
    std::cout << "caught exception " << e.msg() << std::endl;
    return -1;
  }

}
