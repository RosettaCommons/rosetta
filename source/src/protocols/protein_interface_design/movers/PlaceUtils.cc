// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/PlaceUtils.cc
/// @brief a collection of utility functions for placement classes
/// @author Sarel Fleishman (sarelf@u.washington.edu)
#include <protocols/protein_interface_design/movers/PlaceSimultaneouslyMover.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/AA.hh>
#include <numeric/xyzVector.hh>
#include <core/scoring/func/HarmonicFunc.hh>


#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/pack_rotamers.hh>
//#include <core/pack/rotamer_trials.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>

//#include <protocols/docking/DockingProtocol.hh>
#include <protocols/moves/Mover.hh>
#include <core/chemical/ResidueType.hh>
//#include <protocols/moves/ResidueMover.hh>
#include <protocols/hotspot_hashing/HotspotStub.hh>
#include <protocols/hotspot_hashing/HotspotStub.fwd.hh>  // REQUIRED FOR WINDOWS
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <protocols/moves/MoverStatus.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <boost/foreach.hpp>
// Utility Headers
#include <utility/exit.hh>

// Unit Headers
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/hotspot.OptionKeys.gen.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>

// C++ headers
#include <map>
#include <algorithm>

#include <utility/vector0.hh>

//Auto Headers
#include <protocols/simple_moves/DesignRepackMover.hh>



using namespace core::scoring;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.PlaceUtils" );

namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;
using namespace core;

/// @details a utility function that tests whether two residues are aligned properly.
/// Useful for testing whether a single stub should be placed at a given scaffold position
bool
test_res_res_aln( core::conformation::Residue const & res1, core::conformation::Residue const & res2, core::Real & C_N_angle, core::Real & CB_CA_angle  )
{
	using namespace numeric;

	core::Real const threshold( 0.5 );//what angular threshold to use in accepting a stub. 0 implies 90o, 0.5 implies 60o

	Vector const res1_CA = res1.xyz( "CA" );
	Vector const res1_CB = res1.xyz( "CB" );
	Vector const res1_N  = res1.xyz(  "N" );
	Vector const res1_C  = res1.xyz(  "C" );

	Vector const res2_CA = res2.xyz( "CA" );
	Vector const res2_CB = res2.xyz( "CB" );
	Vector const res2_N  = res2.xyz(  "N" );
	Vector const res2_C  = res2.xyz(  "C" );

	Vector const res1_CB_CA = (res1_CB - res1_CA);
	Vector const res1_C_N   = (res1_C  - res1_N);

	Vector const res2_CB_CA = (res2_CB - res2_CA);
	Vector const res2_C_N   = (res2_C  - res2_N);

	C_N_angle = angle_of( res2_C_N, res1_C_N );
	CB_CA_angle = angle_of( res2_CB_CA, res1_CB_CA );

	if( cos_of( res1_C_N, res2_C_N ) <= threshold ||
			cos_of( res1_CB_CA, res2_CB_CA ) <= threshold )
		return( false );
	return( true );
}

core::scoring::constraints::ConstraintCOPs
add_coordinate_constraints( pose::Pose & pose, core::Size const resnum, core::conformation::Residue const rsd_i, core::scoring::func::HarmonicFuncOP coord_cst_func, core::id::AtomID const anchor_atom )
{
	using namespace core::scoring::constraints;

	ConstraintCOPs cst;


	{ // defining constraints scope
		using namespace core::conformation;
		using namespace core::chemical;

		ResidueType const & rsd_type( rsd_i.type() );
		AA const aa = rsd_i.aa();
		switch ( aa ){
			using namespace core::id;
			case( aa_phe ) :
			case( aa_trp ) :
			case( aa_tyr ) :
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CG" ), resnum ), anchor_atom, rsd_i.xyz( "CG" ),coord_cst_func));
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CD1" ), resnum ), anchor_atom, rsd_i.xyz( "CD1" ),coord_cst_func ));
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CD2" ), resnum ), anchor_atom, rsd_i.xyz( "CD2" ),coord_cst_func ));
				break;
			case( aa_gln ) :
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CD" ), resnum ), anchor_atom, rsd_i.xyz( "CD" ),coord_cst_func ));
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OE1" ), resnum ), anchor_atom, rsd_i.xyz( "OE1" ),coord_cst_func ));
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NE2" ), resnum ), anchor_atom, rsd_i.xyz( "NE2" ),coord_cst_func ));
				break;
			case( aa_glu ) :
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CD" ), resnum ), anchor_atom, rsd_i.xyz( "CD" ),coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OE1" ), resnum ), anchor_atom, rsd_i.xyz( "OE1" ),coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OE2" ), resnum ), anchor_atom, rsd_i.xyz( "OE2" ),coord_cst_func ) );
				break;
			case( aa_arg ) :
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NH1" ), resnum ), anchor_atom, rsd_i.xyz( "NH1" ),coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NE" ), resnum ), anchor_atom, rsd_i.xyz( "NE" ),coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NH2" ), resnum ) , anchor_atom, rsd_i.xyz( "NH2" ),coord_cst_func  ) );
				break;
			case( aa_lys ) :
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NZ" ), resnum ), anchor_atom, rsd_i.xyz( "NZ" ),coord_cst_func ) );
				break;
			case( aa_asn ) :
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CG" ), resnum ), anchor_atom, rsd_i.xyz( "CG" ),coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OD1" ), resnum ), anchor_atom, rsd_i.xyz( "OD1" ),coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "ND2" ), resnum ), anchor_atom, rsd_i.xyz( "ND2" ),coord_cst_func ) );
				break;
			case( aa_his ) :
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CG" ), resnum ), anchor_atom, rsd_i.xyz( "CG" ),coord_cst_func ));
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "ND1" ), resnum ) , anchor_atom, rsd_i.xyz( "ND1" ),coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "NE2" ), resnum ), anchor_atom, rsd_i.xyz( "NE2" ),coord_cst_func ) );
				break;
			case( aa_asp ) :
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CG" ), resnum ), anchor_atom, rsd_i.xyz( "CG" ),coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OD1" ), resnum ), anchor_atom, rsd_i.xyz( "OD1" ),coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OD2" ), resnum ), anchor_atom, rsd_i.xyz( "OD2" ),coord_cst_func ) );
				break;
			case( aa_met ) :
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "SD" ), resnum ), anchor_atom, rsd_i.xyz( "SD" ),coord_cst_func ) );
				break;
			case( aa_leu ):
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CG" ), resnum ), anchor_atom, rsd_i.xyz( "CG" ), coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CD1" ), resnum ), anchor_atom, rsd_i.xyz( "CD1" ), coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CD2" ), resnum ), anchor_atom, rsd_i.xyz( "CD2" ), coord_cst_func ) );
				break;
			case( aa_ile ):
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CG1" ), resnum ), anchor_atom, rsd_i.xyz( "CG1" ), coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CD1" ), resnum ), anchor_atom, rsd_i.xyz( "CD1" ), coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CG2" ), resnum ), anchor_atom, rsd_i.xyz( "CG2" ), coord_cst_func ) );
				break;
			case( aa_cys ):
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "SG" ), resnum ), anchor_atom, rsd_i.xyz( "SG" ), coord_cst_func ) );
				// Preserve the CB, since the chi angle is important in disulfides
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CB" ), resnum ), anchor_atom, rsd_i.xyz( "CB" ), coord_cst_func ) );
				break;
			case( aa_thr ):
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CB" ), resnum ), anchor_atom, rsd_i.xyz( "CB" ), coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OG1" ), resnum ), anchor_atom, rsd_i.xyz( "OG1" ), coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CG2" ), resnum ), anchor_atom, rsd_i.xyz( "CG2" ), coord_cst_func ) );
				break;
			case( aa_ala ):
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CB" ), resnum ), anchor_atom, rsd_i.xyz( "CB" ), coord_cst_func ) );
				break;
			case( aa_val ):
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CB" ), resnum ), anchor_atom, rsd_i.xyz( "CB" ), coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CG1" ), resnum ), anchor_atom, rsd_i.xyz( "CG1" ), coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "CG2" ), resnum ), anchor_atom, rsd_i.xyz( "CG2" ), coord_cst_func ) );
				break;
			case( aa_ser ):
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "HG" ), resnum ), anchor_atom, rsd_i.xyz( "HG" ), coord_cst_func ) );
				cst.push_back( new CoordinateConstraint( AtomID( rsd_type.atom_index( "OG" ), resnum ), anchor_atom, rsd_i.xyz( "OG" ), coord_cst_func ) );
				break;
			default :
				utility_exit_with_message( "ERROR: Residue not supported by Placement coordinate constraint machinery" );
		}
	}//defining constraints scope
	TR<<"Constraining residue "<<pose.residue( resnum ).name()<<resnum<<std::endl;
	cst = pose.add_constraints( cst );
	return( cst );
}

/// @details utility function for finding the nearest residue on a chain to
/// a given coordinate. Useful in finding the nearest residue on the host chain to
/// a coordinate constraint
core::Size
find_nearest_residue_to_coord( pose::Pose const & pose, numeric::xyzVector< core::Real > coord, core::Size const host_chain )
{
	core::Real min_dist( 1000000.0 );
	core::Size nearest_res( 0 );
	for( core::Size i( pose.conformation().chain_begin( host_chain ) ),
								end( pose.conformation().chain_end  ( host_chain ) );
								i<=end; ++i ){
		core::Real distance( pose.residue(i).xyz( pose.residue(i).nbr_atom() ).distance( coord ) );
		if( distance <= min_dist ){
			nearest_res = i;
			min_dist = distance;
		}
	}
	runtime_assert( nearest_res );
	return( nearest_res );
}

std::string
nearest_atom_for_constraint( core::conformation::Residue const residue )
{
	using namespace core::chemical;
	AA const aa = residue.aa();
	switch ( aa ){
		using namespace core::id;
		case( aa_phe ) :
		case( aa_trp ) :
		case( aa_asn ) :
		case( aa_his ) :
		case( aa_asp ) :
		case( aa_leu ) :
		case( aa_pro ) :
		case( aa_tyr ) : return( "CG" );
		case( aa_glu ) :
		case( aa_gln ) : return( "CD" );
		case( aa_arg ) : return( "NE" );
		case( aa_lys ) : return( "NZ" );
		case( aa_met ) : return( "SD" );
		case( aa_ile ) : return( "CG1" );
		case( aa_cys ) : return( "SG" );
		case( aa_val ) :
		case( aa_ala ) :
		case( aa_thr ) : return( "CB" );
		case( aa_ser ) : return( "HG" );
		case( aa_gly ) : return( "CA" );
		default :
			utility_exit_with_message( "ERROR: Residue not supported by Placement coordinate constraint machinery" );
	}
	return("");
}

core::scoring::constraints::ConstraintCOPs
add_coordinate_constraints( pose::Pose & pose, core::conformation::Residue const source, core::Size const host_chain, core::Size const resnum, core::Real const coord_sdev, core::scoring::func::HarmonicFuncOP & coord_cst_func )
{
	using namespace core::scoring::constraints;

	ConstraintCOPs cst;

	core::Size const fixed_res( find_nearest_residue_to_coord( pose, source.xyz( nearest_atom_for_constraint( source ) ) ,host_chain == 2 ? 1 : 2 ));
//	core::Size const fixed_res( host_chain == 1 ? pose.total_residue() : 1 );
	TR<<"Anchor residue for the coordinate constraint is "<<fixed_res<<std::endl;
	std::string atom_id( "CB" );
	if( pose.residue( fixed_res ).aa() == core::chemical::aa_gly )
		atom_id = "CA";
	core::id::AtomID const anchor_atom( core::id::AtomID( pose.residue( fixed_res ).atom_index( atom_id ), fixed_res ) );

	if( !coord_cst_func ) coord_cst_func = new core::scoring::func::HarmonicFunc( 0.0, 0.0 );
	coord_cst_func->sd( coord_sdev );
	cst = add_coordinate_constraints( pose, resnum, source, coord_cst_func, anchor_atom );
	return( cst );
}

/// @details common function to both grafting protocols for parsing out the movers that need to be task aware. Returns an empty task factory. If the tag is labeled TaskAware, all of the DesignRepackMovers within it will be assigned
/// this task factory.
void
generate_taskfactory_and_add_task_awareness( utility::tag::TagCOP tag, Movers_map const & movers, basic::datacache::DataMap & data, core::pack::task::TaskFactoryOP & task_factory ){
	using namespace utility::tag;
	using namespace core::pack::task;
	if( !data.has( "TaskFactory", "placement" ) ){
		task_factory = new core::pack::task::TaskFactory;
		data.add( "TaskFactory", "placement", task_factory );
	}
	else
		task_factory = data.get< TaskFactory * >( "TaskFactory", "placement" );
	if( tag->getName() != "NotifyMovers" ) return;
	utility::vector0< TagCOP > const & ta_tags( tag->getTags() );
	BOOST_FOREACH( TagCOP const ta_tag, ta_tags ){
		std::string const mover_name( ta_tag->getOption< std::string >( "mover_name" ) );
		std::map< std::string const, MoverOP >::const_iterator find_mover( movers.find( mover_name ));
		bool const mover_found( find_mover != movers.end() );
		if( mover_found ){
			simple_moves::DesignRepackMoverOP drOP = dynamic_cast< simple_moves::DesignRepackMover * >( find_mover->second.get() );
			if( drOP ){// don't do anything with non-DesignRepackMovers
				TR<<"Setting the task factory of mover "<<find_mover->first<<" to be aware of PlaceSimultaneously's rotamer and sidechain choices.\n";
				drOP->task_factory( task_factory );
			}//fi
		}//fi mover-found
	}// foreach ta_tag
	TR.flush();
}

/// @details utility function for setting up a scorefxn dominated by the stub constraints.
core::scoring::ScoreFunctionOP
make_stub_scorefxn(){
	using namespace core::scoring;
	ScoreFunctionOP stub_scorefxn = new ScoreFunction;
	stub_scorefxn->reset();
	stub_scorefxn->set_weight( backbone_stub_constraint, 10.0 );
	stub_scorefxn->set_weight( fa_rep, 0.44 );
	stub_scorefxn->set_weight( fa_dun, 0.56 );
	stub_scorefxn->set_weight( coordinate_constraint, 1.0 );

	return( stub_scorefxn );
}

utility::vector1< std::pair< protocols::hotspot_hashing::HotspotStubSetOP, std::pair< protocols::hotspot_hashing::HotspotStubOP, core::Size > > >
parse_stub_sets( utility::tag::TagCOP tag, core::pose::Pose const & pose, core::Size const host_chain, basic::datacache::DataMap data ){
	using namespace utility::tag;
	using namespace protocols::hotspot_hashing;
	utility::vector1< PlaceSimultaneouslyMover::StubSetStubPos > stub_sets;

	stub_sets.clear();
	bool contain_StubSet( false );
	utility::vector0< utility::tag::TagCOP > const btags( tag->getTags() );
	utility::vector0< utility::tag::TagCOP >::const_iterator ss_tag = btags.begin();
	for( ; ss_tag != btags.end(); ++ss_tag ){
		if( (*ss_tag)->getName() == "StubSets" ){
			contain_StubSet = true;
			break;
		}
	}//for stubset_tag
	if( !contain_StubSet )
		return stub_sets;
	utility::vector0< utility::tag::TagCOP > const stubset_tags( (*ss_tag)->getTags() );
	BOOST_FOREACH( utility::tag::TagCOP const stubset_tag, stubset_tags ){
		std::string const stub_fname = stubset_tag->getOption< std::string >( "stubfile" );
		HotspotStubSetOP stubset = new HotspotStubSet;
		if( data.has( "hotspot_library", stub_fname ) ){
			stubset = data.get< protocols::hotspot_hashing::HotspotStubSet * >( "hotspot_library", stub_fname );
			TR<<"Associated mover with an already read stubset named "<<stub_fname<<std::endl;
		}
		else{
			TR<<"Associating mover with stub file "<<stub_fname<<std::endl;
			stubset->read_data( stub_fname );
		}
        stub_sets.push_back(std::make_pair(stubset, std::make_pair(new HotspotStub(), 0)));  // REQUIRED FOR WINDOWS
		//stub_sets.push_back( PlaceSimultaneouslyMover::StubSetStubPos( stubset, std::pair< HotspotStubOP, core::Size >( 0, 0 ) ) );

		core::pose::PoseOP ala_pose = new core::pose::Pose( pose );
		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( *ala_pose ));
		task->initialize_from_command_line().or_include_current( true );

		utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, false );
		allowed_aas[ chemical::aa_ala ] = true;

		core::Size const chain_begin( ala_pose->conformation().chain_begin( host_chain ) );
		core::Size const chain_end( ala_pose->conformation().chain_end( host_chain ) );

		for ( core::Size i = 1; i <= pose.total_residue(); i++) {
			if ( !pose.residue(i).is_protein() ) continue;
			if( i >= chain_begin && i <=chain_end ) {
				core::Size const restype( ala_pose->residue(i).aa() );
				if( ( restype == chemical::aa_pro && !basic::options::option[basic::options::OptionKeys::hotspot::allow_proline] )|| restype == chemical::aa_gly )
					task->nonconst_residue_task(i).prevent_repacking();
				else
					task->nonconst_residue_task(i).restrict_absent_canonical_aas( allowed_aas );
			}//fi
			else {
				task->nonconst_residue_task( i ).prevent_repacking();
			}
		}//for i
		if( basic::options::option[basic::options::OptionKeys::packing::resfile].user() )
			core::pack::task::parse_resfile(*ala_pose, *task);

		core::scoring::ScoreFunctionOP scorefxn( get_score_function() );
		pack::pack_rotamers( *ala_pose, *scorefxn, task);
		(*scorefxn)( *ala_pose );
		stubset->pair_with_scaffold( *ala_pose, host_chain, new protocols::filters::TrueFilter );
	}//foreach stubset_tag
	return( stub_sets );
}


} //movers
} //protein_interface_design
} //protocols
