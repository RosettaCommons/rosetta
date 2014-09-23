// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/HotspotDisjointedFoldTreeMover
/// @author Sarel Fleishman (sarelf@uw.edu)

//unit header
#include <protocols/protein_interface_design/movers/HotspotDisjointedFoldTreeMover.hh>
#include <protocols/protein_interface_design/movers/HotspotDisjointedFoldTreeMoverCreator.hh>
//project header
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/protein_interface_design/movers/PlaceUtils.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/random/random.hh>
#include <string>
#include <protocols/loops/loops_main.hh>
#include <boost/foreach.hpp>
#include <protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.hh>
//#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/simple_filters/RotamerBoltzmannWeight.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;
using namespace core::chemical;
using namespace core::kinematics;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.HotspotDisjointedFoldTreeMover" );

std::string HotspotDisjointedFoldTreeMoverCreator::keyname() const
{
	return HotspotDisjointedFoldTreeMoverCreator::mover_name();
}

protocols::moves::MoverOP
HotspotDisjointedFoldTreeMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new HotspotDisjointedFoldTreeMover );
}

std::string
HotspotDisjointedFoldTreeMoverCreator::mover_name() {
	return "HotspotDisjointedFoldTree";
}

HotspotDisjointedFoldTreeMover::HotspotDisjointedFoldTreeMover( ) :
	protocols::moves::Mover( HotspotDisjointedFoldTreeMoverCreator::mover_name()  ),
	ddG_threshold_( 1.0 ),
	chain_( 2 ),
	interface_radius_( 8.0 ),
	scorefxn_( /* NULL */ )
{
	residues_.clear();
}

HotspotDisjointedFoldTreeMover::~HotspotDisjointedFoldTreeMover( ) {}

std::string HotspotDisjointedFoldTreeMover::get_name() const {
	return HotspotDisjointedFoldTreeMoverCreator::mover_name();
}

///@details generates a foldtree that links the nearest residues on chain 1 to key residues on chain 2 and breaks the chain around the key residues on chain 2
core::kinematics::FoldTreeOP
HotspotDisjointedFoldTreeMover::make_disjointed_foldtree( core::pose::Pose const & pose ) const
{
	using namespace core::kinematics;
	FoldTreeOP ft( new FoldTree );

	TR<<"Fold tree before disjointed foldtree:\n"<<pose.fold_tree()<<std::endl;
	ft->clear();
	runtime_assert( chain() == 2 );
/// THIS will only work with chain==2, though reworking it should not be too difficult
	core::Size head( *get_residues().begin()-1 );
	std::set< core::Size > residues_on_target;
	residues_on_target.clear();
	core::Size jump( 1 );
	BOOST_FOREACH( core::Size const r, get_residues() ){
/// Connect the nearest residues on chain1 to key residues on chain2 leaving breaks in chain2 around each key residue
		if( head < r - 1 ){
///connect segments intervening between hot spots
			ft->add_edge( head, r-1, Edge::PEPTIDE );
			ft->add_edge( r-1, r, jump );
			jump++;
		}
/// connect anchor residue on target with hotspot residue
		core::Size const target_res( find_nearest_residue_to_coord( pose, pose.residue( r ).xyz( nearest_atom_for_constraint( pose.residue( r ) )), chain() == 2 ? 1 : 2 ) );
		ft->add_edge( target_res, r, jump );
		residues_on_target.insert( target_res );
		jump++;
		head = r + 1;
	}
	core::Size target_head( pose.conformation().chain_begin( 1 ) );
	BOOST_FOREACH( core::Size const r, residues_on_target ){
/// connect chain1 with no breaks
		ft->add_edge( target_head, r, Edge::PEPTIDE );
		target_head = r;
	}
/// connect the last anchor residue on the target chain with the last residue on the target chain
	core::Size const target_end( *residues_on_target.rbegin() );
	ft->add_edge( target_end, pose.conformation().chain_end( 1 ), Edge::PEPTIDE );

/// refold chain2 from the first key residue to the beginning of the chain and from the last key residue to its end
	core::Size const begin( *get_residues().begin() );
	core::Size const end( *get_residues().rbegin() );
	if( begin - 1 >= pose.conformation().chain_begin( chain() ) ){
		ft->add_edge( begin, begin - 1, jump );
		ft->add_edge( begin - 1, pose.conformation().chain_begin( chain() ), Edge::PEPTIDE );
		jump++;
	}
	if( end + 1 <= pose.conformation().chain_end( chain() ) ){
		ft->add_edge( end, end + 1, jump );
		ft->add_edge( end + 1, pose.conformation().chain_end( chain() ), Edge::PEPTIDE );
	}
	ft->delete_self_edges();
	ft->reorder( 1 );
	TR<<"Fold tree after disjointed foldtree:\n"<<*ft<<std::endl;
	return( ft );
}

protocols::moves::MoverOP
HotspotDisjointedFoldTreeMover::clone() const {
	return( protocols::moves::MoverOP( new HotspotDisjointedFoldTreeMover( *this ) ) );
}

void
HotspotDisjointedFoldTreeMover::apply( core::pose::Pose & pose )
{
	using namespace protocols::toolbox::task_operations;
	using namespace core::pack::task::operation;

	ProteinInterfaceDesignOperationOP pido( new ProteinInterfaceDesignOperation );
	pido->repack_chain1( false );
	pido->repack_chain2( true );
	pido->design_chain2( true );
	pido->interface_distance_cutoff( interface_radius() );
	pido->jump( 1 );

	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory );
	tf->push_back( pido );

	if( ddG_threshold() <= 100 ){
		protocols::simple_filters::RotamerBoltzmannWeight rbw;
		rbw.ddG_threshold( ddG_threshold() );
		rbw.scorefxn( scorefxn() );
		rbw.repacking_radius( interface_radius() );
		rbw.task_factory( tf );
		utility::vector1< core::Size > const ala_scan_res( rbw.first_pass_ala_scan( pose ) );
		BOOST_FOREACH( core::Size const r, ala_scan_res ){ add_residue( r ); }
	}// fi ddG_threshold
	TR<<"Making a disjointed fold tree for residues: ";
	BOOST_FOREACH( core::Size const r, get_residues() )
		TR<<r<<" ";
	TR<<std::endl;
	FoldTreeOP ft( make_disjointed_foldtree( pose ) );
	pose.fold_tree( *ft );
	protocols::loops::add_cutpoint_variants( pose );
}

void
HotspotDisjointedFoldTreeMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & pose )
{
	ddG_threshold( tag->getOption< core::Real >( "ddG_threshold", 1.0 ) );
	if( ddG_threshold() >= 100.0 )
		TR<<"Ala scan calculation will not be carried out. Only residues specifically chosen in the script will be selected"<<std::endl;
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	chain( tag->getOption< core::Size >( "chain", 2 ) );
	interface_radius( tag->getOption< core::Real >( "radius", 8.0 ) );
	utility::vector1< core::Size > v1 = core::pose::get_resnum_list( tag, "resnums", pose );
	BOOST_FOREACH( core::Size const r, v1 ){ add_residue( r ); }

	runtime_assert( ddG_threshold() <= 100.0 || get_residues().size() );
	TR<<"HotspotDisjointedFoldTreeMover with: chain "<<chain()<<" ddG_threshold "<<ddG_threshold()<<" and residues ";
	BOOST_FOREACH( core::Size const r, get_residues() ){ TR<<r<<" "; }
	TR<<std::endl;
}

void
HotspotDisjointedFoldTreeMover::add_residue( core::Size const r ){
	residues_.insert( r );
}

std::set< core::Size >
HotspotDisjointedFoldTreeMover::get_residues() const{
	return( residues_ );
}

void
HotspotDisjointedFoldTreeMover::chain( core::Size const c ){
	chain_ = c;
}

core::Size
HotspotDisjointedFoldTreeMover::chain() const{
	return( chain_ );
}

core::Real
HotspotDisjointedFoldTreeMover::ddG_threshold() const{
	return( ddG_threshold_ );
}

void
HotspotDisjointedFoldTreeMover::ddG_threshold( core::Real const d ){
	ddG_threshold_ = d;
}

void
HotspotDisjointedFoldTreeMover::scorefxn( core::scoring::ScoreFunctionOP scorefxn ){
  scorefxn_ = scorefxn;
}

core::scoring::ScoreFunctionOP
HotspotDisjointedFoldTreeMover::scorefxn() const{
  return scorefxn_;
}

void
HotspotDisjointedFoldTreeMover::interface_radius( core::Real const rad )
{
  interface_radius_ = rad;
}

core::Real
HotspotDisjointedFoldTreeMover::interface_radius() const
{
  return interface_radius_;
}
}//movers
}//protein_interface_design
}//protocols
