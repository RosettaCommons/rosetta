// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/mpi_refinement/Clusterer.cc
/// @brief
/// @author Hahnbeom Park

#include <protocols/mpi_refinement/util.hh>
#include <protocols/mpi_refinement/Clusterer.hh>
#include <protocols/wum/SilentStructStore.hh>
#include <core/import_pose/import_pose.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/pose/Pose.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace mpi_refinement {

static basic::Tracer TR("MPI.LHR.C");

Clusterer::Clusterer()
{
  set_defaults();
}

Clusterer::~Clusterer(){}

void
Clusterer::set_defaults()
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  //simlimit_ = option[ lh::rms_limit ]();
  similarity_method_ = option[ lh::similarity_method ](); 
  similarity_measure_ = option[ lh::similarity_measure ](); 
  sim_replace_obj_ = option[ lh::sim_replace_obj ]();
  simtol_ = option[ lh::similarity_tolerance ]();
  method_ = "energy_sort";
  
}

protocols::wum::SilentStructStore
Clusterer::apply( protocols::wum::SilentStructStore structs,
		core::Size const ncluster,
		core::Real const dist_cut ) const
{
  protocols::wum::SilentStructStore clustered;
  if( method_.compare( "energy_sort" ) == 0 ){
    clustered = energy_sort_cluster( structs, ncluster, dist_cut );
  } else {
    TR << "Unknown clustering method: " << method_ << "!" << std::endl;
  }
  
  return clustered;
}

bool
Clusterer::get_distance( core::io::silent::SilentStructOP ss1, 
			 core::io::silent::SilentStructOP ss2,
			 core::Real &distance,
			 core::Real const dist_cut )  const
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  core::Real dumm;
  bool is_close( false );

  if( similarity_measure_.compare( "Sscore" ) == 0 ){
    distance = CA_Sscore( ss1, ss2, dumm, 2.0 );
    if( distance > dist_cut ){
      is_close = true;
    } else { 
      is_close = false;
    }
	
  } else if (similarity_measure_.compare( "rmsd" ) == 0 ){
    dumm = CA_Sscore( ss1, ss2, distance, 2.0 );
    if( distance < dist_cut ){
      is_close = true;
    } else { 
      is_close = false;
    }
	
  } else if (similarity_measure_.compare( "looprmsd" ) == 0 ){
    std::string loopstr = option[ lh::loop_string ]();
    utility::vector1< core::Size > loopres = loopstring_to_loopvector( loopstr );
    dumm = CA_Sscore( ss1, ss2, distance, loopres, 2.0 );
    if( distance < dist_cut ){
      is_close = true;
    } else { 
      is_close = false;
    }
	} else {
		utility_exit_with_message( "Unknown similarity measure( -lh::similarity_measure )!" );
	}

  return is_close;
}

protocols::wum::SilentStructStore
Clusterer::energy_sort_cluster( protocols::wum::SilentStructStore structs,
				core::Size const ncluster,
				core::Real const dist_cut ) const
{
  core::Size const n( structs.size() );
  runtime_assert( n >= ncluster );

  // First sort based on score
  structs.sort_by( sim_replace_obj_ );

  protocols::wum::SilentStructStore selected;
  utility::vector1< core::Size > selected_id;

  for( core::Size i = 0; i < n-1; ++i ){
    core::io::silent::SilentStructOP ss1 = structs.get_struct( i );

    bool is_close( false );
    for( core::Size j = 0; j < selected.size(); ++j ){
      core::io::silent::SilentStructOP ss2 = selected.get_struct( j );

      core::Real distance;
      is_close = get_distance( ss1, ss2, distance, dist_cut );
      if( is_close ) break;
    }

    if( !is_close ){
      selected_id.push_back( i );
      selected.add( ss1 );
    }

    if( selected.size() >= ncluster ) return selected;
  }

  // If remain, fill up from lowest energy
  core::Size iadd( 0 );
  while( selected.size() < ncluster ){
    if( !selected_id.contains( iadd  ) )
      selected.add( structs.get_struct( iadd ) );
    iadd++;
  }

  return selected;
}

}
}
