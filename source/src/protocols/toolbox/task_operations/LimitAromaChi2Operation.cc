// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/LimitAromaChi2Operation.cc
/// @brief  eliminate aromatic rotamers, of which chi2 are around 0, 180 degree.
/// @detail Chi2=0, 180 rotamers of aromatic residues ( PHE, TYR, HIS ) are not observed in nature very much,
/// however Rosetta really like them. This is really pathology. For design purpose, we don't need them actually.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


// Unit Headers
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2OperationCreator.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <utility/tag/Tag.hh>


#include <core/pack/rotamer_set/RotamerSet.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#ifdef WIN32
	#include <core/graph/Graph.hh>
#endif

namespace protocols {
namespace toolbox {
namespace task_operations {


LimitAromaChi2_RotamerSetOperation::LimitAromaChi2_RotamerSetOperation() :
	RotamerSetOperation(),
	chi2max_( 110 ),
	chi2min_(  70 ),
	include_trp_( false )
{}

LimitAromaChi2_RotamerSetOperation::LimitAromaChi2_RotamerSetOperation( Real const chi2max, Real const chi2min ) :
	RotamerSetOperation(),
	chi2max_( chi2max ),
	chi2min_( chi2min ),
	include_trp_( false )
{}

LimitAromaChi2_RotamerSetOperation::~LimitAromaChi2_RotamerSetOperation() {}

core::pack::rotamer_set::RotamerSetOperationOP
LimitAromaChi2_RotamerSetOperation::clone() const
{
	return core::pack::rotamer_set::RotamerSetOperationOP( new LimitAromaChi2_RotamerSetOperation( *this ) );
}

void
LimitAromaChi2_RotamerSetOperation::alter_rotamer_set(
		Pose const &,
		ScoreFunction const &,
		//mjo commenting out 'ptask' because it is unused and causes a warning
		PackerTask const & /*ptask*/,
		GraphCOP,
		RotamerSet & rotamer_set
)
{

	utility::vector1< bool > rotamers_to_delete;
	rotamers_to_delete.resize( rotamer_set.num_rotamers() );
	Size irot( 0 ), num_to_delete(0);
	for( Rotamers::const_iterator it=rotamer_set.begin(), ite=rotamer_set.end(); it != ite; ++it ) {

		irot ++;
		rotamers_to_delete[ irot ] = false;
		core::conformation::ResidueOP rop ( *it );
		if( rop->aa() == core::chemical::aa_tyr ||
				rop->aa() == core::chemical::aa_phe ||
				rop->aa() == core::chemical::aa_his ||
			 (rop->aa() == core::chemical::aa_trp && include_trp_) )
		{

			runtime_assert( rop->nchi() >= 2 );
			utility::vector1< Real > chi( rop->chi() );

			if( chi[ 2 ] < 0 ) {
				chi[ 2 ] += 180.0;
			}

			if( chi[ 2 ] > chi2max_ || chi[ 2 ] < chi2min_ ) {
				rotamers_to_delete[ irot ] = true;
				num_to_delete++;
			}
		}

	}
	//flo nov 2010: if all the rotamers in the rotamer set are to be deleted,
	//and there is no input rotamer at that position, this will cause a program exit
	//if this is the case (rare), it's probably best to not delete any rotamers
	if( (num_to_delete == rotamer_set.num_rotamers() ) && (rotamer_set.id_for_current_rotamer() == 0 ) ){
		//std::cerr << "shit condition at position " << rotamer_set.resid() << ", not deleting any of the " << rotamer_set.num_rotamers() << " rotamers." << std::endl;
		return;
	}
	rotamer_set.drop_rotamers( rotamers_to_delete );

} // alter_rotamer_set

///////////////////////////////////////////////////////////////////////////////////////////
core::pack::task::operation::TaskOperationOP
LimitAromaChi2OperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new LimitAromaChi2Operation );
}

/// @brief defauot constructor
LimitAromaChi2Operation::LimitAromaChi2Operation():
	TaskOperation(),
	chi2max_( 110 ),
	chi2min_(  70 ),
	include_trp_( false )
{}

/// @brief destructor
LimitAromaChi2Operation::~LimitAromaChi2Operation(){}

/// @brief clone
core::pack::task::operation::TaskOperationOP
LimitAromaChi2Operation::clone() const {
	return core::pack::task::operation::TaskOperationOP( new LimitAromaChi2Operation( *this ) );
}

//mjo commenting out 'pose' because it is unused and causes a warning
/// @brief
void
LimitAromaChi2Operation::apply( Pose const & /*pose*/, PackerTask & task ) const
{
	LimitAromaChi2_RotamerSetOperationOP rso( new LimitAromaChi2_RotamerSetOperation( chi2max_, chi2min_ ) );
	rso->include_trp( include_trp_ );
	task.append_rotamerset_operation( rso );
}

void
LimitAromaChi2Operation::parse_tag( TagCOP tag , DataMap & )
{
	chi2max( tag->getOption< Real >( "chi2max", 110.0 ) );
	chi2min( tag->getOption< Real >( "chi2min", 70.0 ) );
	include_trp( tag->getOption< bool >( "include_trp", 0 ) );
}

void
LimitAromaChi2Operation::parse_def( utility::lua::LuaObject const & def)
{
	chi2max( def["chi2max"] ? def["chi2max"].to<Real>() : 110.0 );
	chi2min( def["chi2min"] ? def["chi2min"].to<Real>() : 70.0 );
}



} // TaskOperations
} // toolbox
} // protocols

