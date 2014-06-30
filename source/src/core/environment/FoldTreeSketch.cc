// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/FoldTreeSketch.cc
/// @author Justin Porter

// Unit Headers
#include <core/environment/FoldTreeSketch.hh>

// Package headers
#include <core/environment/EnvCore.hh>

// Project headers
#include <core/kinematics/FoldTree.hh>

#include <utility/string_util.hh>
#include <numeric/random/random.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers
#include <algorithm>
#include <iterator>
#include <list>
#include <stack>
#include <numeric>

#include <boost/foreach.hpp>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

static numeric::random::RandomGenerator RG(624581); // <- Magic number, do not change it!!!
static basic::Tracer tr("core.environment.FoldTreeSketch", basic::t_info);

namespace core {
namespace environment {

EXCN_FTSketchGraph::EXCN_FTSketchGraph( Size by,
                                        Size on,
                                        std::string const& action,
                                        std::string const& reason ):
Parent( "Error in FoldTreeSketch: Unsuccessful "+
        action+" by Node (seqid) "+utility::to_string(by)+
        " on Node (seqid)"+utility::to_string(on)+"."+reason )
{}

EXCN_FTSketchGraph::EXCN_FTSketchGraph( std::string const& message ):
  Parent( "Error in FoldTreeSketch: "+message )
{}

FoldTreeSketch::FoldTreeSketch():
  ReferenceCount(),
  nodes_(),
  n_jumps_( 0 ),
  n_cuts_( 0 )
{}

FoldTreeSketch::FoldTreeSketch( Size const length ):
  ReferenceCount(),
  nodes_(),
  n_jumps_( 0 ),
  n_cuts_( -1 )
{
  append_peptide( length );
}

FoldTreeSketch::FoldTreeSketch( core::kinematics::FoldTree const& ft ):
  ReferenceCount(),
  nodes_(),
  n_jumps_( 0 ),
  n_cuts_( -1 )
{
  append_peptide( ft.nres() );

  BOOST_FOREACH( core::Size cutpoint, ft.cutpoints() ){
    insert_cut( cutpoint );
  }

  for( int i = 1; i <= (int) ft.num_jump(); ++i ){
    insert_jump( (Size) ft.upstream_jump_residue(i),
                 (Size) ft.downstream_jump_residue(i) );
  }
}

FoldTreeSketch::FoldTreeSketch( FoldTreeSketch const& rhs ):
  ReferenceCount(),
  nodes_(utility::vector1< NodeOP >()),
  n_jumps_( 0 ),
  n_cuts_( rhs.n_cuts_ )
{

  for ( core::Size i = 1; i <= rhs.nodes_.size(); ++i ) {
    nodes_.push_back( new Node( *rhs.nodes_[ i ] ) );
  }
  n_jumps_ = rhs.n_jumps_;
}

bool FoldTreeSketch::has_cut( Size const p ) const {
  range_check( p );
  if( p == nres() ){
    return false;
  }

  bool has_neighbor = nodes_[p]->has_peptide_neighbor( nodes_[p+1].get() );

  return !has_neighbor;
}

bool FoldTreeSketch::has_jump( Size const p1, Size const p2 ) const {
  range_check( p1 );
  range_check( p2 );
  return nodes_[p1]->has_jump_neighbor( nodes_[p2].get() );
}

void FoldTreeSketch::insert_cut( Size const seqid ) {
  range_check( seqid );
  range_check( seqid+1 );

  if( has_cut( seqid ) ){
    throw EXCN_FTSketchGraph( "Rejected cut insertion at "
                              + utility::to_string( seqid )
                              + " because no edge to cut exists." );
  }

  nodes_[seqid]->rm_peptide_neighbor( nodes_[ seqid+1 ].get() );
  n_cuts_ += 1;
}

void FoldTreeSketch::insert_jump( Size const p1, Size const p2 ) {
  range_check( p1 );
  range_check( p2 );
  if( p1 == p2 ){
    throw EXCN_FTSketchGraph( "Jumps cannot connect a node to itself." );
  }

  if( has_jump( p1, p2 ) ){
     throw EXCN_FTSketchGraph( "Rejected jump insertion at "
                               + utility::to_string( p1 ) + ", "
                               + utility::to_string( p2 )
                               + " because the jump already exists." );
  }

  nodes_[p1]->add_jump_neighbor( nodes_[p2].get() );
  n_jumps_ += 1;
}

void FoldTreeSketch::append_peptide( Size length ){
  if( length < 1 ){
    throw utility::excn::EXCN_RangeError( "New FoldTreeSketch polymer stretches must be length >= 1. Obviously." );
  }

  Size old_size = nodes_.size();

  nodes_.push_back( new Node( 1 + old_size ) );

  // If a peptide with more than one residue is being appended, connect them with a single peptide edge
  for( Size i = old_size + 2; i <= old_size + length; ++i ){
    nodes_.push_back( new Node( i ) );
    nodes_[ i ]->add_peptide_neighbor( nodes_[ i - 1 ].get() );
  }

  n_cuts_ += 1;
}

void FoldTreeSketch::render( core::kinematics::FoldTree& ft ) const{
  std::set< std::pair< Size, Size > > jumps;
  std::set< Size > cuts;

  // TODO: Checks of putative fold tree integrity go here.

  // Run through the graph and pull out jumps and cuts.
  for( Size i = 1; i <= nodes_.size() - 1; ++i ){
    if( !nodes_[i]->has_peptide_neighbor( nodes_[i+1].get() ) ){
      cuts.insert( i );
    }
    nodes_[i]->collect_jumps( jumps );
  }

  if( cuts.size() != jumps.size() ){
    tr.Debug << "Cuts were: ";
    BOOST_FOREACH( Size cut, cuts ){ tr.Debug << cut << ", "; }
    tr.Debug << std::endl << "Jumps were: ";
    for( std::set< std::pair< Size, Size > >::iterator it = jumps.begin(); it != jumps.end(); ++it ) {
      tr.Debug << "(" << it->first << "," << it->second << "), ";
    }
    tr.Debug << std::endl;

    throw EXCN_FTSketchGraph( "Number of jumps ("+utility::to_string( jumps.size() )+
                              ") and number of cuts ("+utility::to_string(cuts.size())+
                              ") must be equal for rendering. Try randomized cut deletion first." );
  }

  // Populate the stupid FCL arrays for use in the fold tree constructor thing.
  ObjexxFCL::FArray1D_int cut_array( (int) cuts.size() );
  ObjexxFCL::FArray2D_int jump_array( 2, (int) jumps.size() );

  std::set< std::pair< core::Size, core::Size > >::iterator jump_it = jumps.begin();
  std::set< core::Size >::iterator cut_it = cuts.begin();
  for( int i = 1; i <= (int) jumps.size(); ++i ){
    cut_array(i) = (int) *cut_it;
    jump_array( 1, i ) = (int) jump_it->first;
    jump_array( 2, i ) = (int) jump_it->second;

    ++cut_it;
    ++jump_it;
  }

  // With everything set up, now the fold tree method can be used.
  bool success = ft.tree_from_jumps_and_cuts( (int) nodes_.size(), (int) jumps.size(),
                                              jump_array, cut_array);
  if( !success ){
    throw EXCN_FTSketchGraph( "FoldTree generation was not successful." );
  }
}

core::kinematics::FoldTreeOP FoldTreeSketch::render() const{
  kinematics::FoldTreeOP ft = new kinematics::FoldTree();
  render( *ft );
  return ft;
}

Size FoldTreeSketch::nres() const {
  return nodes_.size();
}

Size FoldTreeSketch::num_jumps() const {
  return n_jumps_;
}

Size FoldTreeSketch::num_cuts() const {
  return n_cuts_;
}

void FoldTreeSketch::range_check( Size const seqpos ) const {
  if( seqpos > nres() || seqpos < 1 ){
    throw EXCN_FTSketchGraph( "Sequence position "+utility::to_string(seqpos)+
                              " is not a valid sequence position (nres="+
                              utility::to_string( nres() )+")");
  }
}

utility::vector1< Size > const FoldTreeSketch::cycle( core::Size const start_resid ) const{
  std::stack< NodeCAP > cycle_path;

  range_check( start_resid );
  NodeCAP start = nodes_[ start_resid ].get();

  start->has_cycle( cycle_path, start );

  // If it's a disconnected net, we might not've found the cycle. Check for
  // disconnected nodes and recurse on the first disconnected node we find.
  if( cycle_path.empty() ){
    for( utility::vector1< NodeOP >::const_iterator n_it = nodes_.begin();
        n_it != nodes_.end(); ++n_it ){
      if( !(*n_it)->visited() ){
        return cycle( (*n_it)->seqid() );
      }
    }
  }

  utility::vector1< Size > out_vect;
  while( !cycle_path.empty() ){
    out_vect.push_back( cycle_path.top()->seqid() );
    cycle_path.pop();
  }

  for( utility::vector1< NodeOP >::const_iterator n_it = nodes_.begin();
       n_it != nodes_.end(); ++n_it ){
    (*n_it)->unvisit();
  }

  return out_vect;
}

core::Size FoldTreeSketch::insert_cut( utility::vector1< Real > const& bias ){
  if( bias.size() != this->nres() ){
    throw EXCN_FTSketchGraph( "Cut bias array does not match FoldTreeSize." );
  }

  core::Real const sum = std::accumulate( bias.begin(), bias.end(), 0 );

  if( sum <= 0.0 ){
    throw EXCN_FTSketchGraph( "Cut bias array sum <= 0.0" );
  }

  core::Real rand = RG.uniform() * sum;

  for( Size seqpos = 1; seqpos <= nres(); ++seqpos ){
    rand -= bias[seqpos];
    if( rand <= 0 ){
      tr.Debug << "FoldTreeSketch inserting random cut at " << seqpos << std::endl;
      insert_cut( seqpos );
      return seqpos;
    }
  }

  std::ostringstream ss;
  ss << "Random cut insertion failed. Check your bias array for bad values. Cut bias array: " << bias;
  throw EXCN_FTSketchGraph( ss.str() );
}

std::set< core::Size > FoldTreeSketch::remove_cycles( utility::vector1< Real > const& bias ){
  if( bias.size() != this->nres() ){
    throw EXCN_FTSketchGraph( "Cut bias array does not match FoldTreeSize." );
  }

  utility::vector1< Size > cycle = this->cycle();
  if( cycle.size() != 0 && !cuttable( cycle ) ){
    throw EXCN_FTSketchGraph( "All-jump cycles are not automatically resolvable: "+
                              utility::to_string( cycle )+"." );
  }

  std::set< Size > new_cuts;

  while( cycle.size() != 0 ){
    core::Real sum = 0.0;

    if( tr.Trace.visible() )
      tr.Trace << "found cycle: " << utility::to_string( cycle ) << std::endl;

    for( utility::vector1< Size >::iterator resid_it = cycle.begin();
        resid_it != cycle.end(); ++resid_it ){
      sum += bias[ *resid_it ];
    }
    if( sum == 0 ){
      throw EXCN_FTSketchGraph( "All-zero cut biased cycles cannot be automatically resolved. Cycle: "
                                +utility::to_string( cycle )+", Biases: "+utility::to_string( bias )+"." );
    }

    // use a randomly generated number [0, sum] to choose which peptide edge to cut
    core::Real rand = RG.uniform() * sum;
    for( Size cycle_id = 1; cycle_id <= cycle.size(); ++cycle_id ){
      Size seqpos = cycle[cycle_id];
      rand -= bias[ seqpos ];
      if( rand <= 0 ){ //this seqpos was selected by the random number
        // Verify that this seqpos and the next seqpos is in the cycle, and that
        // there is a peptide bond between them.
        if( std::find( cycle.begin(), cycle.end(), seqpos+1 ) != cycle.end() &&
            nodes_[seqpos]->has_peptide_neighbor( nodes_[ seqpos + 1 ].get() ) ){
          nodes_[seqpos]->rm_peptide_neighbor( nodes_[ seqpos + 1 ].get() );
          new_cuts.insert( seqpos );
          tr.Debug << "removing cycle by cutting at " << seqpos << std::endl;
        } else {
          tr.Debug << "seqpos " << seqpos << " was selected, but no peptide edge exists to cut. Trying again." << std::endl;
        }
        break;
      }
    }

    cycle = this->cycle();
  }

  return new_cuts;
}

bool FoldTreeSketch::cuttable( utility::vector1< Size > const& cycle ) const {
  for( Size cycle_id = 1; cycle_id <= cycle.size(); ++cycle_id ){
    Size seqpos = cycle[cycle_id];
    Size next_seqpos = cycle[ ( cycle_id % cycle.size() )+1 ];
    if( nodes_[ seqpos ]->has_peptide_neighbor( nodes_[next_seqpos].get() ) ){
      return true;
    }
  }

  return false;
}

///////////////////////////////////////
/// NODE IMPLEMENTATION METHODS
///////////////////////////////////////

FoldTreeSketch::Node::Node( Size i ):
  seqid_( i ),
  pep_prev_( NULL ),
  pep_next_( NULL ),
  parent_( NULL )
{}

void FoldTreeSketch::Node::add_peptide_neighbor( NodeAP n ){
  assert( n );
  assert( n != this );

  if ( ( this->seqid() - n->seqid() ) == 1 ) {
    if ( pep_prev_ &&
         pep_prev_ != this ) {
      throw EXCN_FTSketchGraph( this->seqid(), n->seqid(), "add_pep_neighbor",
                                "Node already has a previous peptide member." );
    }
    pep_prev_ = n;
    assert( !n->pep_next_ || n->pep_next_ == this );
    n->pep_next_ = this;
  } else if ( ( (int) this->seqid() - (int) n->seqid() ) == -1 ){
    if( pep_next_ &&
        pep_next_ != this ){
      throw EXCN_FTSketchGraph( this->seqid(), n->seqid(), "add_pep_neighbor",
                               "Node already has a next peptide member." );
    }
    pep_next_ = n;
    assert( !n->pep_prev_ || n->pep_prev_ == this );
    n->pep_prev_ = this;
  } else {
    throw EXCN_FTSketchGraph( this->seqid(), n->seqid(), "add_pep_neighbor",
                              "Peptide connections exist only between sequence-adacjent nodes." );
  }
}

void FoldTreeSketch::Node::add_jump_neighbor( NodeAP n ){
  assert( n );
  assert( n != this );

  jump_neighbors_.insert( n.get() );
  n->jump_neighbors_.insert( this );
}

bool FoldTreeSketch::Node::has_jump_neighbor( NodeCAP n ) const {
  assert( n );
  return jump_neighbors_.find( const_cast< Node* >( n.get() ) ) != jump_neighbors_.end();
}

bool FoldTreeSketch::Node::has_peptide_neighbor( NodeCAP n ) const {
  assert( n );
  // For some reason, raw pointer comparison doesn't work here. Use id instead.

  bool is_next = ( pep_next_ && n->seqid() == pep_next_->seqid() );
  bool is_prev = ( pep_prev_ && n->seqid() == pep_prev_->seqid() );

  return ( is_prev || is_next );
}

bool FoldTreeSketch::Node::has_neighbor( NodeCAP n ) const {
  return has_peptide_neighbor( n ) || has_jump_neighbor( n );
}

void FoldTreeSketch::Node::rm_jump_neighbor( NodeAP n ){
  assert( n );
  assert( n != this );

  jump_neighbors_.erase( n.get() );
  n->jump_neighbors_.erase( this );
}

void FoldTreeSketch::Node::rm_peptide_neighbor( NodeAP n ){
  assert( n );
  assert( n != this );

  if( pep_next_ == n ){
    pep_next_->pep_prev_ = NULL;
    pep_next_ = NULL;
  } else if ( pep_prev_ == n ){
    pep_prev_->pep_next_ = NULL;
    pep_next_ = NULL;
  } else {
    throw EXCN_FTSketchGraph( this->seqid(), n->seqid(), "rm_peptide_neighbor",
                              "Neighbor does not exist to be removed." );
  }
}

Size FoldTreeSketch::Node::seqid() const {
  return seqid_;
}

bool FoldTreeSketch::Node::has_cycle( std::stack< NodeCAP >& path, NodeCAP caller ) const {
  assert( caller );

  // Base case: if we've been visited before, there's a cycle.
  if( parent_ != 0 ){
    NodeCAP p = caller;
    while( p != parent_ ){
      path.push( p );
      p = p->parent_;
    }
    return true;
  }

  // If we haven't been visited, recurse on neighbors and note that we've
  // now been visited by caller.
  parent_ = caller;

  if( pep_next_ != caller &&
      pep_next_ != 0 ) {
    if( pep_next_->has_cycle( path, this ) )
      return true;
  }

  if( pep_prev_ != caller &&
      pep_prev_ != 0 ) {
    if( pep_prev_->has_cycle( path, this ) ){
      return true;
    }
  }

  for( EdgeList::const_iterator it = jump_neighbors_.begin(); it != jump_neighbors_.end(); ++it ){
    if( it->get() != caller.get() ){
      if( (*it)->has_cycle(path, this) ){
        return true;
      }
    }
  }

  return false;
}

void FoldTreeSketch::Node::collect_jumps( std::set< std::pair< Size, Size > >& jumps ) const {
  for( EdgeList::const_iterator it = jump_neighbors_.begin(); it != jump_neighbors_.end(); ++it ){
    Size seqid1 = seqid();
    Size seqid2 = (*it)->seqid();

    if( seqid1 > seqid2 ){
      jumps.insert( std::make_pair( seqid2, seqid1 ) );
    } else {
      jumps.insert( std::make_pair( seqid1, seqid2 ) );
    }
  }
}

} // environment
} // core
