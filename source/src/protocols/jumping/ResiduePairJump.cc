// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/jumping/ResiduePairJump.cc
/// @brief a class to create jump transform from a pair of residues
/// @details
///  This class is to create possible jump transforms between a pair of residues
///  when their sidechains are locked in certain geometry constraints. It starts
///  from the predefined constraints and takes backbone-independent rotamer conformation
///  into account and reversely generate their backbone positions. Then a jump transform is
///  measured.
/// @author Chu Wang

// Unit Headers
#include <protocols/jumping/ResiduePairJump.hh>

// Package Headers

// Project Headers
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/JumpSRFD.hh>
#include <core/fragment/JumpingFrame.hh>
#include <core/graph/Graph.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// ObjexxFCL Headers
// Utility headers
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <utility/exit.hh>

//// C++ headers
#include <string>

#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jumping {


ResiduePairJumpSingle::ResiduePairJumpSingle() {}

ResiduePairJumpSingle::~ResiduePairJumpSingle() = default;


/// @brief constructed by residue_type
ResiduePairJumpSingle::ResiduePairJumpSingle(
	core::chemical::ResidueType const & residue_type
) :
	residueType_( residue_type.clone() ),
	fixResidue_( ! residue_type.is_protein() )
{}
/////////////////////////////////////////////////////////////////////
//ResiduePairJumpSingle
void
ResiduePairJumpSingle::set_jumpAtoms(
	utility::vector1< std::string > const & jump_atoms
)
{
	jumpAtoms_ = jump_atoms;
}

/////////////////////////////////////////////////////////////////////
void
ResiduePairJumpSingle::set_jumpAtoms(
	core::Size i,
	std::string const & atom_name
)
{
	jumpAtoms_.resize(i,"");
	jumpAtoms_[i] = atom_name;
}
/////////////////////////////////////////////////////////////////////
void
ResiduePairJumpSingle::set_cstAtoms(
	utility::vector1< std::string > const & cst_atoms
)
{
	cstAtoms_ = cst_atoms;
}
/////////////////////////////////////////////////////////////////////
void
ResiduePairJumpSingle::set_cstAtoms(
	core::Size i,
	std::string const & atom_name
)
{
	cstAtoms_.resize(i,"");
	cstAtoms_[i] = atom_name;
}

core::chemical::ResidueTypeOP
ResiduePairJumpSingle::residueType() const
{
	return residueType_;
}


/////////////////////////////////////////////////////////////////////
//ResiduePairJump


// empty constructor
ResiduePairJump::ResiduePairJump() {}
ResiduePairJump::~ResiduePairJump() = default;

// constructed by two input ResidueTypes
ResiduePairJump::ResiduePairJump(
	core::chemical::ResidueType const & residue1,
	core::chemical::ResidueType const & residue2
)
{
	residues_.push_back( ResiduePairJumpSingleOP( new ResiduePairJumpSingle( residue1 ) ) );
	residues_.push_back( ResiduePairJumpSingleOP( new ResiduePairJumpSingle( residue2 ) ) );
}

void
ResiduePairJump::add_residue_pair(
	core::chemical::ResidueType const & residue1,
	core::chemical::ResidueType const & residue2
)
{
	residues_.clear();
	cstInfoMap_.clear();
	residues_.push_back( ResiduePairJumpSingleOP( new ResiduePairJumpSingle( residue1 ) ) );
	residues_.push_back( ResiduePairJumpSingleOP( new ResiduePairJumpSingle( residue2 ) ) );
}
/////////////////////////////////////////////////////////////////////
void
ResiduePairJump::add_residue_single(
	core::chemical::ResidueType const & residue
)
{
	if ( residues_.size() < 2 ) {
		cstInfoMap_.clear();
		residues_.push_back( ResiduePairJumpSingleOP( new ResiduePairJumpSingle(residue) ) );
	} else {
		utility_exit_with_message("ResiduePairJump can only take two residues\n");
	}
}


/////////////////////////////////////////////////////////////////////
/// @details add a value to cst type
void
ResiduePairJump::set_cstInfo(
	cstType type,
	core::Real value
)
{
	auto it = cstInfoMap_.find( type );
	if ( it == cstInfoMap_.end() ) {
		// not in map yet, add it
		utility::vector1< core::Real > value_vector;
		value_vector.push_back( value );
		cstInfoMap_.insert( std::make_pair( type, value_vector ) );
	} else {
		utility::vector1< core::Real > & value_vector( it->second );
		bool redundant = false;
		for ( core::Size i = 1; i <= value_vector.size(); ++i ) {
			if ( value == value_vector[i] ) {
				redundant = true;
				break;
			}
		}
		if ( ! redundant ) value_vector.push_back( value );
	}
}
/////////////////////////////////////////////////////////////////////
/// @details add a vector of values to cstInfo type
void
ResiduePairJump::set_cstInfo(
	cstType type,
	utility::vector1< core::Real > const & values
)
{
	for ( core::Size i = 1; i <= values.size(); ++i ) {
		set_cstInfo( type, values[i] );
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
core::fragment::FrameOP
ResiduePairJump::generate_frame()
{
	using namespace core;
	using namespace core::fragment;

	// create a frame with proper data structures Frame->FragData->(UpJumpSRFD and DownJumpSRFD)
	core::kinematics::RT jump_rt;
	UpJumpSRFDOP up_jump( new UpJumpSRFD(residues_[1]->residueType()->name1()) );
	DownJumpSRFDOP down_jump( new DownJumpSRFD( jump_rt, jumpAtoms(2), jumpAtoms(1), residues_[2]->residueType()->name1() ) );
	FragDataOP jump_frag( new FragData() );
	jump_frag->add_residue(up_jump);
	jump_frag->add_residue(down_jump);
	JumpingFrameOP ResiduePairJumpFrame( new JumpingFrame( 1, 2, 2 ) ); // from residue 1 to 2
	ResiduePairJumpFrame->set_pos( 1, 1);
	ResiduePairJumpFrame->set_pos( 2, 2);
	ResiduePairJumpFrame->add_fragment( jump_frag );

	diversify_dof_conformers();

	for ( Size i = 1; i <= dof_conformers_.size(); ++i ) {
		apply_dof_conformer( dof_conformers_[i] );
		//miniPose_->dump_pdb("chutmp_div"+right_string_of(i,4, '0') + ".pdb");
		ResiduePairJumpFrame->steal( *miniPose_ );
	}

	return ResiduePairJumpFrame;
}
/////////////////////////////////////////////////////////////////////
void
ResiduePairJump::init_mini_pose()
{
	using namespace core;

	runtime_assert( cstAtoms_defined() );

	miniPose_ = core::pose::PoseOP( new core::pose::Pose() );

	// initial configuration for how these two residues are connected
	//  Real init_phi = numeric::conversions::radians(180.0);
	//  Real init_theta = numeric::conversions::radians(60.0);
	//  Real init_d = 2.0;
	//  // connection numbers for new connections added to these residues
	//  Size anchor_connection, root_connection;
	// operation for each of the two residues created and added
	for ( Size i = 1; i <= 2; ++i ) {
		ResiduePairJumpSingleOP rsd = residues_[i];
		//  // create atomIcoord
		//   chemical::AtomICoor icoor( init_phi, init_theta, init_d,
		//    rsd->cstAtoms(1), rsd->cstAtoms(2), rsd->cstAtoms(3), *(rsd->residueType()) );
		//   // add a residue connection
		//   Size connection = rsd->residueType()->add_residue_connection( rsd->cstAtoms(1) );
		//   rsd->residueType()->residue_connection(connection).icoor( icoor );
		conformation::ResidueOP new_residue( conformation::ResidueFactory::create_residue( *( rsd->residueType() ) ) );

		if ( i == 1 ) {
			//anchor_connection = connection;
			miniPose_->append_residue_by_jump( *new_residue, 1 ); // the first residue attached by jump
		} else {
			//root_connection = connection;
			// the second residue is attched to the first residue by bond with defined connection
			//miniPose_->append_residue_by_bond( *new_residue, true, root_connection, 1, anchor_connection, false );
			miniPose_->append_residue_by_jump( *new_residue, 1, residues_[1]->cstAtoms(1), rsd->cstAtoms(1), false );
		}
	}
	// build rotamers for each of the sidechains if applicable and add
	build_sidechain_rotamers();

	build_cst_conformer_jumps();

	//setup_cstTypeToDofMap();

	//miniPose_->dump_pdb("chutmp_miniPose.pdb");

	return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// void ResiduePairJump::setup_cstTypeToDofMap()
// {
//  // erase any existing data
//  cstTypeToDofMap_.clear();

//  // six cst types will be mapped to DOFs of the three atoms in residue 2
//  ResiduePairJumpSingleOP rsd = residues_[2];
//  core::chemical::ResidueTypeOP rsd_type = rsd->residueType();
//  core::id::AtomID atom_id1( rsd_type->atom_index(rsd->cstAtoms(1)), 2 );
//  core::id::AtomID atom_id2( rsd_type->atom_index(rsd->cstAtoms(2)), 2 );
//  core::id::AtomID atom_id3( rsd_type->atom_index(rsd->cstAtoms(3)), 2 );
//  core::id::AtomID bogus_id( 0, 0 );

//  core::id::DOF_ID dof_id12( miniPose_->atom_tree().bond_length_dof_id( atom_id1, atom_id2 ) );
//  core::id::DOF_ID dof_id23( miniPose_->atom_tree().bond_length_dof_id( atom_id2, atom_id3 ) );
//  core::id::DOF_ID dof_id13( miniPose_->atom_tree().bond_length_dof_id( atom_id1, atom_id3 ) );
//  // sanity check: downstream cstAtoms 1 2 and 3 have to be connected in atom_tree, either 1->2->3 or 2<-1->3
//  if ( ! ( dof_id12.valid() && ( dof_id23.valid() || dof_id13.valid() ) ) ) {
//   utility_exit_with_message("setup_cstTypeToDofMap() downstream cstAtoms are not connected in atom tree!\n");
//  }

//  for ( cstInfoMapIterator it = cstInfoMap_.begin(), it_end = cstInfoMap_.end(); it != it_end; ++it ) {
//   cstType cst_type = it->first;
//   core::id::AtomID atom_id;
//   core::id::DOF_Type dof_type;
//   if ( cst_type == disAB )  {
//    atom_id = atom_id1;
//    dof_type = core::id::D;
//   } else if ( cst_type == angleA ) {
//    atom_id = atom_id1;
//    dof_type = core::id::THETA;
//   } else if ( cst_type == dihedralA ) {
//    atom_id = atom_id1;
//    dof_type = core::id::PHI;
//   } else if ( cst_type == angleB ) {
//    atom_id = atom_id2;
//    dof_type = core::id::THETA;
//   } else if ( cst_type == dihedralAB ) {
//    atom_id = atom_id2;
//    dof_type = core::id::PHI;
//   } else if ( cst_type == dihedralB ) {
//    runtime_assert( dof_id23.valid() ); // dihedralB only meaningful when atom3 is connected atom2
//    atom_id = atom_id3;
//    dof_type = core::id::PHI;
//   } else if ( cst_type == rot1 || cst_type == rot2 ) {
//    // for sidechain rot1 and rot2, they do not correspond to one single DOF
//    // so set it as BOGUS_ID to prevent accidental lookup
//    atom_id = bogus_id;
//    dof_type = core::id::PHI;
//   } else {
//    utility_exit_with_message("unknown cstType defined in ResiduePairJump\n" );
//   }
//   core::id::DOF_ID dof_id( atom_id, dof_type);
//   cstTypeToDofMap_.insert( std::make_pair( cst_type, dof_id ) );
//  }

//  return;
// }
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResiduePairJump::diversify_dof_conformers(
	dofType type,
	core::Size max_index
)
{
	if ( max_index == 0 ) return;

	if ( dof_conformers_.size() == 0 ) {
		// is empty, generate initial vector
		for ( core::Size i = 1; i <= max_index; ++i ) {
			std::map< dofType, core::Size > map;
			map.insert( std::make_pair( type, i ) );
			dof_conformers_.push_back( map );
		}
	} else {
		// conformers exist, diversify them if not redundant
		utility::vector1< std::map< dofType, core::Size > > tmp_conformers= dof_conformers_;
		dof_conformers_.clear();
		for ( core::Size j = 1; j <= tmp_conformers.size(); ++j ) {
			std::map< dofType, core::Size > map( tmp_conformers[j] );
			if ( map.find(type) == map.end() ) {
				map.insert( std::make_pair( type, 0 ) );
				for ( core::Size i = 1; i <= max_index; ++i ) {
					map[type] = i;
					dof_conformers_.push_back( map );
				}
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResiduePairJump::diversify_dof_conformers()
{

	dof_conformers_.clear();
	diversify_dof_conformers( rot1, rotsets_[1]->num_rotamers() );
	diversify_dof_conformers( rot2, rotsets_[2]->num_rotamers() );
	diversify_dof_conformers( cstJump, cst_jumps_.size() );

	return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
core::pose::PoseOP
ResiduePairJump::apply_dof_conformer(
	std::map< dofType, core::Size > const & conformer_map
)
{
	 for ( auto const & it : conformer_map ) {
		dofType type = it.first;
		core::Size index = it.second;
		if ( ( type == rot1 ) || (type == rot2) ) { // sidechain rotamer dofs
			core::Size seqpos = (type==rot1) ? 1 : 2;
			core::conformation::ResidueCOP rotamer = rotsets_[seqpos]->rotamer(index);
			for ( core::Size i = 1; i <= rotamer->nchi(); ++i ) {
				miniPose_->set_chi( i, seqpos, rotamer->chi()[i] );
			}
		} else { // cst_jump dofs
			ResiduePairJumpSingleOP rsd = residues_[2];
			core::id::AtomID jump_atom_id( rsd->residueType()->atom_index(rsd->cstAtoms(1)), 2 );
			miniPose_->set_jump( jump_atom_id, cst_jumps_[index] );
		}
	}
	return miniPose_;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResiduePairJump::build_sidechain_rotamers()
{
	using namespace core::scoring;

	runtime_assert( miniPose_->total_residue() == 2 );

	// set up a scorefxn and packer_neighbor_graph for build_rotamers function
	ScoreFunctionOP scorefxn( get_score_function() );
	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( *miniPose_ );
	task->initialize_from_command_line().restrict_to_repacking();
	(*scorefxn)(*miniPose_);
	utility::vector1< bool > task_mask( miniPose_->total_residue(), true );
	for ( core::Size which = 1; which <= residues_.size(); ++which ) {
		ResiduePairJumpSingleOP rsd = residues_[which];

		// create one residue pose and generate an empty rotset
		core::pack::rotamer_set::RotamerSetFactory rsf;
		core::pack::rotamer_set::RotamerSetOP rotset = rsf.create_rotamer_set( miniPose_->residue(which) );
		rotset->set_resid(which);
		// unless this residue needs to be fixed, build rotamers into this RotamerSet
		if ( ! rsd->fixResidue() ) {
			//put phi/psi at popular beta-strand conformation for most unconstrained rotamer conforamtions
			miniPose_->set_phi(which, -90.0);
			miniPose_->set_psi(which, 120.0);
			// set up a simple packer task with this single residue enabled
			task_mask[which] = true;
			task->restrict_to_residues( task_mask );
			core::graph::GraphOP packer_neighbor_graph( core::pack::create_packer_graph( *miniPose_, *scorefxn, task ) );
			// build possible rotamers for this residue
			rotset->build_rotamers( *miniPose_, *scorefxn, *task, packer_neighbor_graph );
			task_mask[which] = false;
		}
		// add it into rotsets_ container
		rotsets_.push_back( rotset );
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResiduePairJump::build_cst_conformer_jumps()
{
	typedef  numeric::xyzVector< core::Real > Vector;

	// generate all combinations of cst_conformers
	diversify_cst_conformers();

	// set up atoms to define the jump transform, from cstAtoms
	utility::vector1< Vector > atoms_xyz;
	for ( Size i = 1; i <= 2; ++i ) {
		for ( Size j = 1; j <= 3 ; ++j ) {
			atoms_xyz.push_back( miniPose_->xyz( core::id::NamedAtomID( jumpAtoms(i)[j], i ) ) );
		}
	}
	// for each cst conformer conmbination, generate a jump transform
	for ( Size i = 1; i <= cst_conformers_.size(); ++i ) {
		std::map< cstType, core::Size > const & conformer( cst_conformers_[i] );
		utility::vector1< core::Real > cst_values(6);
		 for ( auto const & it : conformer ) {
			cstType type = it.first;
			Size index = it.second;
			core::Real value = cstInfoMap_.find(type)->second[index];
			if ( type != disAB ) value = numeric::conversions::radians(value);
			cst_values[Size(type)] = value;
		}
		core::kinematics::Jump jump;
		jump.from_bond_cst( atoms_xyz, cst_values );
		cst_jumps_.push_back( jump );
	}
	return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResiduePairJump::diversify_cst_conformers(
	cstType type,
	core::Size max_index
)
{
	if ( cst_conformers_.size() == 0 ) {
		// is empty, generate initial vector
		for ( core::Size i = 1; i <= max_index; ++i ) {
			std::map< cstType, core::Size > map;
			map.insert( std::make_pair( type, i ) );
			cst_conformers_.push_back( map );
		}
	} else {
		// conformers exist, diversify them if not redundant
		utility::vector1< std::map< cstType, core::Size > > tmp_conformers= cst_conformers_;
		cst_conformers_.clear();
		for ( core::Size j = 1; j <= tmp_conformers.size(); ++j ) {
			std::map< cstType, core::Size > map( tmp_conformers[j] );
			if ( map.find(type) == map.end() ) {
				map.insert( std::make_pair( type, 0 ) );
				for ( core::Size i = 1; i <= max_index; ++i ) {
					map[type] = i;
					cst_conformers_.push_back( map );
				}
			}
		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
ResiduePairJump::diversify_cst_conformers()
{
	for ( auto & it : cstInfoMap_ ) {
		diversify_cst_conformers( it.first, it.second.size() );
	}
	return;
}
///////////////////////////////////////////////////////////////////////////////////////////////
} //protocols
} //jumping


