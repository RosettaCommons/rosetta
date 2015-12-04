// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/connection/BridgeChains.cc
/// @brief The BridgeChains
/// @detailed
/// @author Tom Linsky

//Unit Headers
#include <protocols/denovo_design/connection/BridgeChains.hh>
#include <protocols/denovo_design/connection/BridgeChainsCreator.hh>

//Project Headers
#include <protocols/denovo_design/components/Picker.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>

//Protocol Headers
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/SegmentRebuild.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/methods/fold_tree_functions.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

//Core Headers
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constants.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/ABEGOManager.hh>

//Basic Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/string_util.hh>
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

//C++ Headers
#include <stack>
#include <math.h>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.connection.BridgeChains" );

///////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace connection {
///////////////////////////////////////////////////////////////////////////////

// coordinate cst rcg class
class CoordinateCstRCG : public protocols::forge::remodel::RemodelConstraintGenerator {
public:
	CoordinateCstRCG() { residues_.clear(); }
	virtual ~CoordinateCstRCG() {}
	virtual std::string get_name() const { return "CoordinateCstRCG"; }
	virtual void generate_remodel_constraints( core::pose::Pose const & pose ) {
		for ( core::Size idx=1; idx<=residues_.size(); ++idx ) {
			add_constraint( create_coordinate_cst( pose, residues_[idx] ) );
		}
	}

	void add_residue( core::Size const res ) { residues_.push_back(res); }

	/// @brief creates a Ca coordinate constraint for residue resi
	core::scoring::constraints::ConstraintOP
	create_coordinate_cst( core::pose::Pose const & pose,
		core::Size const resi ) const
	{
		core::Size atom( pose.residue_type(resi).nbr_atom() );
		if ( pose.residue_type(resi).has("CA") ) {
			atom = pose.residue_type(resi).atom_index("CA");
		}

		return core::scoring::constraints::ConstraintOP(
			new core::scoring::constraints::CoordinateConstraint(
			core::id::AtomID(atom,resi),
			core::id::AtomID(pose.residue(1).nbr_atom(),1),
			pose.residue(resi).xyz(atom),
			core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc(0.0, 0.5) ) ) );
	}

private:
	utility::vector1< core::Size > residues_;
};
typedef utility::pointer::shared_ptr< CoordinateCstRCG > CoordinateCstRCGOP;

////////////////////////////////////////////////////////////////////////////////

std::string
BridgeChainsCreator::keyname() const
{
	return BridgeChainsCreator::mover_name();
}

protocols::moves::MoverOP
BridgeChainsCreator::create_mover() const
{
	return protocols::moves::MoverOP( new BridgeChains() );
}

std::string
BridgeChainsCreator::mover_name()
{
	return "BridgeChains";
}

///  ---------------------------------------------------------------------------------
///  BridgeChains main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
BridgeChains::BridgeChains() :
	Connection(),
	scorefxn_(),
	frag_picker_( components::PickerOP( new components::Picker() ) )
{
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
BridgeChains::~BridgeChains() {}


/// Return a copy of ourselves
protocols::moves::MoverOP
BridgeChains::clone() const
{
	return protocols::moves::MoverOP( new BridgeChains( *this ) );
}

/// @brief return a fresh instance of ourselves
protocols::moves::MoverOP
BridgeChains::fresh_instance() const
{
	return protocols::moves::MoverOP( new BridgeChains() );
}

/// @brief return a fresh instance of ourselves
std::string
BridgeChains::get_name() const
{
	return BridgeChainsCreator::mover_name();
}

void
BridgeChains::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	Connection::parse_my_tag( tag, data, filters, movers, pose );

	if ( tag->hasOption( "scorefxn" ) ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
		TR << "score function " << tag->getOption< std::string >( "scorefxn" ) << " is used. " << std::endl;
	}
}

/// @brief configures based on a permutation
/// @throw EXCN_Setup if no valid connection endpoints are found
void
BridgeChains::setup_permutation( components::StructureData & perm ) const
{
	Connection::setup_permutation( perm );

	if ( !segments_fixed( perm ) ) {
		std::stringstream ss;
		ss << id() << ": two components were given as input thare are not fixed relative to one another." << perm << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( ss.str() );
	}
}

/// @brief Does the work of remodeling the connection
void
BridgeChains::apply_connection( components::StructureData & perm ) const
{
	core::scoring::ScoreFunctionCOP my_scorefxn = scorefxn();
	if ( !my_scorefxn ) {
		std::stringstream err;
		err << "BridgeChains: You must set a valid scorefunction to "
			<< id() << " before connecting" << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}

	// save phi/psi/omega
	core::Real const c1e_omega = perm.pose()->omega( perm.segment( loop_lower(perm) ).stop() );

	// strip terminal residues
	perm.delete_trailing_residues( loop_lower(perm) );
	perm.delete_leading_residues( loop_upper(perm) );

	if ( do_remodel() ) {
		// store original in case there is a failure
		components::StructureDataOP orig = perm.clone();

		try {
			apply_constraints( perm );
			// remodel the loop
			build_loop( perm );
		} catch ( utility::excn::EXCN_Base const & e ) {
			perm = *orig;
			throw EXCN_ConnectionFailed( id() );
		}
		remove_constraints( perm );
	}

	perm.declare_covalent_bond(
		loop_lower(perm), perm.segment( loop_lower(perm) ).length(), "C",
		loop_upper(perm), 1, "N" );
	perm.set_omega( perm.segment( loop_lower(perm) ).stop(), c1e_omega );

	// remove jump/cutpoint
	perm.delete_jump_and_intervening_cutpoint( loop_lower(perm), loop_upper(perm) );
}

/// @brief builds the loop
void
BridgeChains::build_loop( components::StructureData & perm ) const
{
	debug_assert( perm.pose() );
	static bool const const_fold_tree = true;

	protocols::moves::DsspMover dssp;
	perm.apply_mover( dssp );

	core::pose::Pose const & pose = *(perm.pose());

	// information about actual new loop residues
	core::Size const loopstart = perm.segment( lower_segment_id( perm ) ).cterm_resi() + 1;
	core::Size const loopend = perm.segment( upper_segment_id( perm ) ).nterm_resi() - 1;
	core::Size const insert_length = build_len( perm );

	// setup left and right overlap regions
	core::Size const left = build_left( perm );
	core::Size const right = build_right( perm );
	debug_assert( left <= right );

	// connection ss and abego
	std::string const conn_ss = build_ss(perm);
	std::string const conn_abego = build_abego(perm);
	TR << "Going to build ss = " << conn_ss << " abego = " << conn_abego << std::endl;

	utility::vector1< std::string > complete_abego =
		core::sequence::ABEGOManager().get_symbols( pose, 1, pose.total_residue(), 1 );

	utility::vector1< std::string > const abego_ins = abego_insert( complete_abego, conn_abego, left, right, loopstart, loopend );
	std::string const ss = ss_insert( pose, conn_ss, left, right, loopstart, loopend );
	std::string const aa = aa_insert( pose, insert_length, left, right, loopstart, loopend );

	std::string complete_ss;
	complete_ss += perm.ss().substr(0,left-1);
	complete_ss += ss;
	complete_ss += perm.ss().substr(right,std::string::npos);
	debug_assert( complete_ss.size() == pose.total_residue() );

	for ( core::Size i=left, c=1; i<=right; ++i ) {
		complete_abego[i] = abego_ins[c];
		++c;
	}

	TR << "Original size is " << pose.total_residue() << "; Going to build the secondary structure: " << ss << " between residues " << left << " and " << right << " to cover the loop starting at residue " << loopstart << " and ending at residue " << loopend << " with abego " << abego_ins << std::endl;

	// setup cutpoints
	core::Size const cutres = perm.segment( loop_lower(perm) ).cterm_resi();

	// setup coordinate csts
	CoordinateCstRCGOP coord_cst( new CoordinateCstRCG() );
	for ( core::Size i=left; i<loopstart; ++i ) {
		TR << "Creating cst for residue " << i << std::endl;
		coord_cst->add_residue( i );
	}
	for ( core::Size i=loopend+1; i<=right; ++i ) {
		TR << "Creating cst for residue " << i << std::endl;
		coord_cst->add_residue( i );
	}

	protocols::loops::LoopsOP loops( new protocols::loops::Loops() );
	debug_assert( cutres >= left );
	debug_assert( cutres <= right );
	loops->add_loop( left, right, cutres );

	protocols::moves::MoverOP remodel =
		create_remodel_mover(
		pose,
		loops,
		const_fold_tree,
		complete_ss,
		complete_abego,
		left,
		right );

	// switch to centroid if necessary
	bool const input_centroid = perm.pose()->is_centroid();
	core::pose::Pose stored;

	if ( !input_centroid ) {
		stored = *(perm.pose());
		perm.switch_residue_type_set( "centroid" );
	}

	perm.apply_mover( coord_cst );
	perm.apply_mover( remodel );

	if ( remodel->get_last_move_status() != protocols::moves::MS_SUCCESS ) {
		throw EXCN_ConnectionFailed( "during remodel in " + id() );
	}

	// Rebuild fold tree
	utility::vector1< std::string > roots;
	roots.push_back( lower_segment_id(perm) );
	roots.push_back( upper_segment_id(perm) );
	perm.consolidate_movable_groups( roots );
	TR << "After remodel: " << perm.pose()->fold_tree() << std::endl;
	TR << "CLOSED THE LOOP" << std::endl;

	// switch back to FA if necessary
	if ( !input_centroid ) {
		perm.switch_residue_type_set( "fa_standard" );
		copy_rotamers( perm, stored );
	}
}

/// @brief creates mover that does remodeling of loop residues
protocols::moves::MoverOP
BridgeChains::create_remodel_mover(
	core::pose::Pose const & pose,
	protocols::loops::LoopsOP loops,
	bool const const_fold_tree,
	std::string const & complete_ss,
	StringVec const & complete_abego,
	core::Size const left,
	core::Size const right ) const
{
	debug_assert( frag_picker_ );

	protocols::forge::remodel::RemodelLoopMoverOP remodel(
		new protocols::forge::remodel::RemodelLoopMover( loops )
	);
	debug_assert( remodel );
	remodel->set_keep_input_foldtree( const_fold_tree );

	core::fragment::FragSetOP frag3;
	core::fragment::FragSetOP frag9;
	if ( do_remodel() ) {
		frag9 = frag_picker_->pick_and_cache_fragments( complete_ss, complete_abego, left, right, 9 );
		frag3 = frag_picker_->pick_and_cache_fragments( complete_ss, complete_abego, left, right, 3 );
		remodel->add_fragments( frag9 );
		remodel->add_fragments( frag3 );
	}

	// turn off interchain score terms if we only have one chain in the final structure
	core::scoring::ScoreFunctionOP sfx = scorefxn()->clone();
	if ( pose.conformation().num_chains() == 2 ) {
		sfx->set_weight( core::scoring::interchain_env, 0.0 );
		sfx->set_weight( core::scoring::interchain_pair, 0.0 );
		sfx->set_weight( core::scoring::interchain_contact, 0.0 );
		sfx->set_weight( core::scoring::interchain_vdw, 0.0 );
	}
	remodel->scorefunction( *sfx );

	// setup movemap for remodel so that only the target loops can move
	core::kinematics::MoveMap mm;
	// initialize explicitly to all false
	for ( core::Size res=1; res<=pose.total_residue(); ++res ) {
		mm.set_bb( res, false );
		mm.set_chi( res, false );
	}
	// then set loop residues (plus overlap) to true
	for ( core::Size l=left; l<=right; ++l ) {
		mm.set_bb( l, true );
		mm.set_chi( l, true );
	}
	remodel->false_movemap( mm );

	return remodel;
}

/// @brief checks to ensure that both pieces being connected are fixed relative to one another
/// @details If segments are fixed, this returns true.  If the segments
/// are allowed to move relative to one another, this returns false.
/// This function is used primarily for checking user input
bool
BridgeChains::segments_fixed( components::StructureData const & perm ) const
{
	if ( !do_remodel() ) {
		return true;
	}

	std::set< core::Size > mj1, mj2;
	std::set< std::string > visited;
	std::stack< std::string > ids;
	ids.push( lower_segment_id(perm) );
	while ( ids.size() ) {
		std::string const test_id = ids.top();
		ids.pop();
		TR.Debug << "looking at " << test_id << " group " << perm.segment(test_id).movable_group << std::endl;
		visited.insert( test_id );
		mj1.insert( perm.segment(test_id).movable_group );
		if ( !perm.segment(test_id).has_free_lower_terminus() ) {
			if ( visited.find( perm.segment(test_id).lower_segment() ) == visited.end() ) {
				ids.push( perm.segment(test_id).lower_segment() );
			}
		}
		if ( !perm.segment(test_id).has_free_upper_terminus() ) {
			if ( visited.find( perm.segment(test_id).upper_segment() ) == visited.end() ) {
				ids.push( perm.segment(test_id).upper_segment() );
			}
		}
	}
	visited.clear();
	ids.push( upper_segment_id(perm) );
	while ( ids.size() ) {
		std::string const test_id = ids.top();
		ids.pop();
		TR.Debug << "looking2 at " << test_id << " group " << perm.segment(test_id).movable_group << std::endl;
		visited.insert( test_id );
		mj2.insert( perm.segment(test_id).movable_group );
		if ( !perm.segment(test_id).has_free_lower_terminus() ) {
			if ( visited.find(perm.segment(test_id).lower_segment()) == visited.end() ) {
				ids.push( perm.segment(test_id).lower_segment() );
			}
		}
		if ( !perm.segment(test_id).has_free_upper_terminus() ) {
			if ( visited.find(perm.segment(test_id).upper_segment()) == visited.end() ) {
				ids.push( perm.segment(test_id).upper_segment() );
			}
		}
	}

	bool fixed = false;
	for ( std::set< core::Size >::const_iterator it1 = mj1.begin(); it1 != mj1.end(); ++it1 ) {
		if ( mj2.find(*it1) != mj2.end() ) {
			fixed = true;
			break;
		}
	}
	return fixed;
}

/// @brief using the motif list, find the desired abego for each position in the connection
utility::vector1< std::string >
BridgeChains::abego_insert(
	StringVec const & complete_abego,
	std::string const & connection_abego,
	core::Size const left,
	core::Size const right,
	core::Size const end1,
	core::Size const start2 ) const
{
	debug_assert( ! complete_abego.empty() );
	debug_assert( complete_abego.size() >= right );
	utility::vector1< std::string > retval;
	// residues to be connected are obliterated by the vlb, so we need to add them to the ss/abego
	for ( core::Size i=left; i<end1; ++i ) {
		retval.push_back( complete_abego[ i ] );
	}
	for ( core::Size i=1, s=connection_abego.size(); i<=s; ++i ) {
		std::string a = "";
		a += connection_abego[i-1];
		retval.push_back( a );
	}
	for ( core::Size i=start2+1; i<=right; ++i ) {
		// residues on the right must also be added
		retval.push_back( complete_abego[ i ] );
	}
	return retval;
}

/// @brief using the motif list and input pose, find the desired aa sequence for each position in the connection. default="V"
std::string
BridgeChains::aa_insert(
	core::pose::Pose const & pose,
	core::Size const connection_len,
	core::Size const left,
	core::Size const right,
	core::Size const end1,
	core::Size const start2 ) const
{
	std::string retval;
	for ( core::Size i=left; i<end1; ++i ) {
		retval += pose.residue(i).name1();
	}
	for ( core::Size i=1; i<=connection_len; ++i ) {
		retval += 'V';
	}
	for ( core::Size i=start2+1; i<=right; ++i ) {
		retval += pose.residue(i).name1();
	}
	return retval;
}

/// @brief using the motif list and input pose, find the desired secondary structure for each position in the connection
std::string
BridgeChains::ss_insert(
	core::pose::Pose const & pose,
	std::string const & connection_ss,
	core::Size const left,
	core::Size const right,
	core::Size const end1,
	core::Size const start2 ) const
{
	assert( right <= pose.total_residue() );
	assert( left >= 1 );
	TR << "end1=" << end1 << " start2=" << start2 << " left=" << left << " right=" << right << std::endl;
	std::string retval;
	for ( core::Size i=left; i<end1; ++i ) {
		retval += pose.secstruct(i);
	}
	retval += connection_ss;
	for ( core::Size i=start2+1; i<=right; ++i ) {
		retval += pose.secstruct(i);
	}
	return retval;
}

/// @brief returns the scorefunction, creating it if necessary
/// @details Default scorefunction is fldsgn_cen
core::scoring::ScoreFunctionCOP
BridgeChains::scorefxn() const
{
	if ( !scorefxn_ ) {
		TR.Info << "Using default scorefunction fldsgn_cen.wts" << std::endl;
		return core::scoring::ScoreFunctionFactory::create_score_function( "fldsgn_cen.wts" );
	}
	return scorefxn_;
}

} // namespace connection
} // namespace denovo_design
} // namespace protocols
