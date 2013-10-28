// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AlignEndsMover.cc
/// @brief

// Unit headers
#include <devel/splice/AlignEndsMover.hh>
#include <devel/splice/AlignEndsMoverCreator.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static basic::Tracer TR("devel.splice.AlignEndsMover");
#include <utility/tag/Tag.hh>

// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.hh>
#include <string>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <core/import_pose/import_pose.hh>
#include <devel/splice/FindEndpointsOperation.hh>
#include <protocols/toolbox/superimpose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <numeric/xyzVector.hh>
#include <algorithm>

namespace devel {
namespace splice {

using namespace::protocols;

std::string
AlignEndsMoverCreator::keyname() const
{
	return AlignEndsMoverCreator::mover_name();
}

protocols::moves::MoverOP
AlignEndsMoverCreator::create_mover() const {
	return new AlignEndsMover;
}

std::string
AlignEndsMoverCreator::mover_name()
{
	return "AlignEnds";
}

AlignEndsMover::AlignEndsMover(): moves::Mover("AlignEnds"),
	distance_threshold_( 18.0 ),
	neighbors_( 6 ),
	N_terminal_count_( 3 ),
	odd_( true ),
	even_( true ),
	template_pose_( NULL ),
	stagger_( 0 )
{
}

AlignEndsMover::~AlignEndsMover()
{
}

void
AlignEndsMover::template_pose( core::pose::PoseOP p ){ template_pose_ = p; }

core::pose::PoseOP
AlignEndsMover::template_pose() const{ return template_pose_; }

utility::vector1< core::Size >
AlignEndsMover::reference_positions( core::pose::Pose const & pose ) const{
	using namespace core::scoring::dssp;

	Dssp dssp( pose );
	dssp.dssp_reduced();

	utility::vector1< core::Size > strand_positions;
	strand_positions.clear();
	for( core::Size i = 1; i <= pose.total_residue() - 2; ++i ){ /// 3-strand positions in a row
      if( dssp.get_dssp_secstruct( i ) == 'E' && dssp.get_dssp_secstruct( i + 1) == 'E' && dssp.get_dssp_secstruct( i + 2 ) == 'E' )
          strand_positions.push_back( i );
  }
  utility::vector1< core::Size > strand_ntermini;
  strand_ntermini.clear();
  core::Size prev_resnum( 99999 );
  foreach( core::Size const resi, strand_positions ){ /// find strand ntermini
      if( resi != prev_resnum + 1 )
          strand_ntermini.push_back( resi );
      prev_resnum = resi;
  }
  utility::vector1< core::Size > ntermini_w_neighbors;
  ntermini_w_neighbors.clear();
	foreach( core::Size const res, strand_ntermini ){ /// strand_ntermini with enough neighbors
      core::Size const resi_neighbors( neighbors_in_vector( pose, res, strand_ntermini, distance_threshold(), dssp ));
      if( resi_neighbors >= neighbors() ){
          ntermini_w_neighbors.push_back( res );
          TR<<res<<" has "<<resi_neighbors<<" neighbors"<<std::endl;
      }//fi resi_neighbors
  }//foreach res

  utility::vector1< core::Size > ntermini_w_close_neighbors; // ntermini with 2 close neighbors
  ntermini_w_close_neighbors.clear();
	foreach( core::Size const res, ntermini_w_neighbors ){
      core::Size const close_neighbors( neighbors_in_vector( pose, res, strand_ntermini, 10.0, dssp ));
      if( close_neighbors >= 2 )
          ntermini_w_close_neighbors.push_back( res );
      else{
          core::Size const close_neighbors_one_down( neighbors_in_vector( pose, res+1, strand_ntermini, 10.0, dssp ));
          if( close_neighbors_one_down >= 2 )
              ntermini_w_close_neighbors.push_back( res + 1);
      }
      TR<<res<<" has "<<close_neighbors<<" close neighbors"<<std::endl;
  }

	utility::vector1< core::Size > positions_for_alignment;
	positions_for_alignment.clear();
	bool curr_odd( true );
	foreach( core::Size const res, ntermini_w_close_neighbors ){
		if( ( curr_odd && odd() ) || ( !curr_odd && even() ) ){
			for( core::Size count = 0; count < N_terminal_count(); ++count )
				positions_for_alignment.push_back( res + count );
		}//fi
	}//foreach
	return positions_for_alignment;
}

utility::vector1< numeric::xyzVector< core::Real > >
Ca_coords( core::pose::Pose const & pose, utility::vector1< core::Size > const positions ){
	utility::vector1< numeric::xyzVector< core::Real > > coords;

	coords.clear();
	foreach( core::Size const pos, positions ){
		coords.push_back( pose.residue( pos ).xyz( "CA" ) );
	}
	return coords;
}

void
AlignEndsMover::apply( Pose & pose ){
	using namespace protocols::toolbox;
	using namespace core::scoring::dssp;

	utility::vector1< core::Size > const template_positions( reference_positions( *template_pose() ) ), pose_positions( reference_positions( pose ) );
	TR<<"Aligning positions on template: ";
	foreach( core::Size const p, template_positions ){
		TR<<p<<'+';}
	TR<<std::endl;
	TR<<"To reference positions: ";
	foreach( core::Size const p, pose_positions ){
		TR<<p<<'+';}
	TR<<std::endl;

	utility::vector1< numeric::xyzVector< core::Real > > init_coords( Ca_coords( pose, pose_positions ) ), ref_coords( Ca_coords( *template_pose(), template_positions ) );

	numeric::xyzMatrix< core::Real > rotation;
	numeric::xyzVector< core::Real > to_init_center, to_fit_center;

	std::rotate( init_coords.begin(), init_coords.begin() + ( N_terminal_count() * stagger()), init_coords.end());

	superposition_transform( init_coords, ref_coords, rotation, to_init_center, to_fit_center );

	apply_superposition_transform( pose, rotation, to_init_center, to_fit_center );
}

std::string
AlignEndsMover::get_name() const {
	return AlignEndsMoverCreator::mover_name();
}

moves::MoverOP
AlignEndsMover::clone() const
{
	return new AlignEndsMover( *this );
}

moves::MoverOP
AlignEndsMover::fresh_instance() const
{
	return new AlignEndsMover;
}

void
AlignEndsMover::parse_my_tag(
	utility::tag::TagPtr const tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	distance_threshold( tag->getOption< core::Real >( "distance_threshold", 18.0 ));
	neighbors( tag->getOption< core::Size >( "neighbors", 6 ) );
	N_terminal_count( tag->getOption< core::Size >( "N_terminal_count", 3 ) );
	odd( tag->getOption< bool >( "odd", true ) );
	even( tag->getOption< bool >( "even", true ) );
	std::string const template_fname( tag->getOption< std::string >( "template_pose" ) );
	template_pose( core::import_pose::pose_from_pdb( template_fname ) );
	stagger( tag->getOption< core::Size >( "stagger", 0 ) );

	TR<<"distance_threshold: "<<distance_threshold()<<" neighbors: "<<neighbors()<<" N_terminal_count: "<<N_terminal_count()<<" odd: "<<odd()<<" even: "<<even()<<" template_pose: "<<template_fname<<" stagger: "<<stagger()<<std::endl;
}
} // simple_moves
} // protocols

