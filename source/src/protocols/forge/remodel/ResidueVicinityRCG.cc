// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/ResidueVicinityRCG.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, april 2009

// unit headers
#include <protocols/forge/remodel/ResidueVicinityRCG.hh>


//package headers
#include <protocols/forge/remodel/ResidueVicinityCstGeneratorCreator.hh>

//project headers
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/id/SequenceMapping.hh>
#include <basic/Tracer.hh>
#include <protocols/constraints_additional/AmbiguousMultiConstraint.hh>

//utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

static basic::Tracer tr( "protocols.forge.remodel.ResidueVicinityRCG" );

namespace protocols{
namespace forge{
namespace remodel{

ResidueVicinityInfo::ResidueVicinityInfo(
	core::Size old_seqpos,
	utility::vector1< core::Size > const & residue_atoms,
	utility::vector1< core::Size > const & loopres_atoms,
	core::Size desired_remodelres_in_vicinity
) :
	old_seqpos_( old_seqpos ),
	residue_atoms_( residue_atoms ),
	loopres_atoms_(loopres_atoms),
	dis_(NULL),
	loop_ang_(NULL),
	targ_ang_(NULL),
	loop_dih_(NULL),
	targ_dih_(NULL),
	lt_dih_(NULL),
	desired_remodelres_in_vicinity_(desired_remodelres_in_vicinity)
{}

ResidueVicinityInfo::~ResidueVicinityInfo() {}

core::scoring::constraints::FuncOP
ResidueVicinityInfo::dis() const {
	return dis_; }

core::scoring::constraints::FuncOP
ResidueVicinityInfo::loop_ang() const{
	return loop_ang_; }

	core::scoring::constraints::FuncOP
ResidueVicinityInfo::targ_ang() const{
	return targ_ang_; }

	core::scoring::constraints::FuncOP
ResidueVicinityInfo::loop_dih() const{
	return loop_dih_; }

	core::scoring::constraints::FuncOP
ResidueVicinityInfo::targ_dih() const{
	return targ_dih_; }

core::scoring::constraints::FuncOP
ResidueVicinityInfo::lt_dih() const{
	return lt_dih_; }

void
ResidueVicinityInfo::set_dis( core::scoring::constraints::FuncOP dis ){
	dis_ = dis; }

	void
ResidueVicinityInfo::set_loop_ang( core::scoring::constraints::FuncOP loop_ang ){
	loop_ang_ = loop_ang; }

	void
ResidueVicinityInfo::set_targ_ang( core::scoring::constraints::FuncOP targ_ang ){
	targ_ang_ = targ_ang; }

void
ResidueVicinityInfo::set_loop_dih( core::scoring::constraints::FuncOP loop_dih ){
	loop_dih_ = loop_dih; }

void
ResidueVicinityInfo::set_targ_dih( core::scoring::constraints::FuncOP targ_dih ){
	targ_dih_ = targ_dih; }

void
ResidueVicinityInfo::set_lt_dih( core::scoring::constraints::FuncOP lt_dih ){
	lt_dih_ = lt_dih; }


//ResidueVicinityConstraintsCreator Functions
protocols::moves::MoverOP
ResidueVicinityCstGeneratorCreator::create_mover() const
{
	return new ResidueVicinityRCG();
}

std::string
ResidueVicinityCstGeneratorCreator::keyname() const
{
	return ResidueVicinityCstGeneratorCreator::mover_name();
}

std::string
ResidueVicinityCstGeneratorCreator::mover_name()
{
	return "ResidueVicinityCstCreator";
}

ResidueVicinityRCG::ResidueVicinityRCG()
	: RemodelConstraintGenerator(),
		lstart_( 0 ),
		lstop_( 0 )
{}

ResidueVicinityRCG::ResidueVicinityRCG( ResidueVicinityRCG const & rval )
	: RemodelConstraintGenerator( rval ),
		lstart_( rval.lstart_ ),
		lstop_( rval.lstop_ )
{}

ResidueVicinityRCG::ResidueVicinityRCG(
	core::Size lstart,
	core::Size lstop,
	utility::vector1< ResidueVicinityInfoOP > const & rv_infos
) : RemodelConstraintGenerator(),
		lstart_(lstart), lstop_(lstop),
		rv_infos_(rv_infos)
{}

ResidueVicinityRCG::~ResidueVicinityRCG() {}

void
ResidueVicinityRCG::parse_my_tag( TagCOP const tag,
																	basic::datacache::DataMap & data,
																	protocols::filters::Filters_map const & filters,
																	protocols::moves::Movers_map const & movers,
																	core::pose::Pose const & pose )
{
	RemodelConstraintGenerator::parse_my_tag( tag, data, filters, movers, pose );
	lstart( tag->getOption< core::Size >( "lstart", lstart_ ) );
	if ( lstart_ == 0 ) {
		throw utility::excn::EXCN_RosettaScriptsOption("lstart must be specified in ResidueVicinityCstGenerator mover");
	}
	lstop( tag->getOption< core::Size >( "lstart", lstop_ ) );
	if ( lstop_ == 0 ) {
		throw utility::excn::EXCN_RosettaScriptsOption("lstop must be specified in ResidueVicinityCstGenerator mover");
	}
}

std::string
ResidueVicinityRCG::get_name() const
{
	return ResidueVicinityCstGeneratorCreator::mover_name();
}


protocols::moves::MoverOP
ResidueVicinityRCG::fresh_instance() const
{
	return new ResidueVicinityRCG();
}

protocols::moves::MoverOP
ResidueVicinityRCG::clone() const
{
	return new ResidueVicinityRCG( *this );
}

/// @brief for every ResidueVicinityInfo (RVI) that this generator has, an AmbiguousConstraint between the
/// @brief neighbor atoms of all residues in the remodel region and the RVI are generated.
void
ResidueVicinityRCG::generate_remodel_constraints(
	core::pose::Pose const & pose  )
{

	this->clear_constraints();

	//pose.dump_pdb("remodel_start.pdb"); //debug

	//std::set< Interval > regions = vlb()->manager().intervals_containing_undefined_positions();

	core::Size remstart( lstart_ );
	core::Size remend( lstop_ );

	if( this->seqmap() ){
		remstart = (*(this->seqmap() ))[ remstart ];
		remend = (*(this->seqmap() ))[ remend ];
	}


	tr << "setting up constraints for res " << remstart << " to res " << remend << std::endl;

	for( utility::vector1< ResidueVicinityInfoOP >::const_iterator rv_it = rv_infos_.begin(), rv_end = rv_infos_.end();
			 rv_it != rv_end; ++rv_it ){


		//now generate constraints between all nbr atoms of the remodel region positions and all specified atoms
		//in the rvi. Note: if more than one atom is specified in the rvi, an ambiguous constraint will be
		//created between all specified atoms and the residue nb atom
		utility::vector1< core::scoring::constraints::ConstraintCOP > rv_csts;

		for( core::Size i = remstart; i <= remend; ++i){

			utility::vector1< core::scoring::constraints::ConstraintCOP > respair_csts;

			generate_remodel_constraints_for_respair(pose, i, *(*rv_it), respair_csts );

			if( respair_csts.size() == 1 ){
				rv_csts.push_back( respair_csts[1] );
				//debug
				//tr << "respair_cst pushing back single constraint for loopres " << i << " and rvinfo with old seqpos " << (*rv_it)->old_seqpos() << std::endl;
			}

			else{
				rv_csts.push_back( new core::scoring::constraints::AmbiguousConstraint( respair_csts ) );
				//tr << "respair_cst pushing back ambig constraint for loopres " << i << " and rvinfo with old seqpos " << (*rv_it)->old_seqpos() << " containing " << respair_csts.size() << " members" << std::endl;
			}

		} //loop over all residues of the remodel region


			//and finally add the constraints to the RGC
		if( (*rv_it)->desired_remodelres_in_vicinity() == 1 ){
			this->add_constraint( new core::scoring::constraints::AmbiguousConstraint( rv_csts ) );
			//debug
			//tr << "creating ambig constraint of size " << rv_csts.size() << std::endl;
		}

		else{
			this->add_constraint( new protocols::constraints_additional::AmbiguousMultiConstraint( (*rv_it)->desired_remodelres_in_vicinity(), rv_csts ) );

			//debug
			//tr << "creating ambig multi constraint of size " << rv_csts.size() << " and n " << (*rv_it)->desired_remodelres_in_vicinity() <<  std::endl;
		}

	} //loop over residue vicinity infos

} //generate constraints


void
ResidueVicinityRCG::generate_remodel_constraints_for_respair(
	core::pose::Pose const & pose,
	core::Size const loopres,
	ResidueVicinityInfo const & rv_info,
	utility::vector1< core::scoring::constraints::ConstraintCOP > & csts
){

	using namespace core::scoring::constraints;
	csts.clear();

	//core::Size target_pos = ( *(vlb()->manager().sequence_mapping()) )[ rv_info->old_seqpos() ];

	core::Size target_pos( rv_info.old_seqpos() );

	if( this->seqmap() ){
		target_pos = (*(this->seqmap() ))[ target_pos ];
	}

	//generate the penalty function for all constraints of this rvi

	//FuncOP cst_func = new BoundFunc(
	//rv_info.distance() - rv_info.tolerance(), rv_info.distance() + rv_info.tolerance(), 0.1, "rem_cstfunc" );

	for( core::Size i =1; i <= rv_info.residue_atoms().size(); ++i ){

		for( core::Size j =1; j <= rv_info.loopres_atoms().size(); ++j ){

			//csts.push_back( new AtomPairConstraint( core::id::AtomID( rv_info.loopres_atoms()[j], loopres ),core::id::AtomID( rv_info.residue_atoms()[i], target_pos ), cst_func ) );

			utility::vector1< core::scoring::constraints::ConstraintCOP > csts_this_pair;

			if( rv_info.dis() ) csts_this_pair.push_back(	new AtomPairConstraint( core::id::AtomID( rv_info.loopres_atoms()[j], loopres ),core::id::AtomID( rv_info.residue_atoms()[i], target_pos ), rv_info.dis() )	);
			//tr << "creating dist constraints between looppos " << loopres << " atom " <<  pose.residue_type( loopres ).atom_name( rv_info->loopres_atoms()[j] ) << " and targ pos " << target_pos << " atom " << pose.residue_type( target_pos ).atom_name( rv_info->residue_atoms()[i] ) << std::endl;

			if( rv_info.loop_ang() ){
				csts_this_pair.push_back( new AngleConstraint( core::id::AtomID( rv_info.loopres_base_atoms()[j], loopres ), core::id::AtomID( rv_info.loopres_atoms()[j], loopres ), core::id::AtomID( rv_info.residue_atoms()[i], target_pos ), rv_info.loop_ang() ) );

				//tr << "creating loop ang constraints between looppos " << loopres << " base atom " << pose.residue_type( loopres ).atom_name( rv_info.loopres_base_atoms()[j] ) << ", atom " <<  pose.residue_type( loopres ).atom_name( rv_info.loopres_atoms()[j] ) << " and targ pos " << target_pos << " atom " << pose.residue_type( target_pos ).atom_name( rv_info.residue_atoms()[i] ) << std::endl;
			}

			if( rv_info.targ_ang() ){
				csts_this_pair.push_back( new AngleConstraint( core::id::AtomID( rv_info.loopres_atoms()[j], loopres ), core::id::AtomID( rv_info.residue_atoms()[i], target_pos ), core::id::AtomID( rv_info.residue_base_atoms()[i], target_pos ), rv_info.targ_ang() ) );

				//tr << "creating targ ang constraints between looppos " << loopres << " atom " << pose.residue_type( loopres ).atom_name( rv_info.loopres_atoms()[j] ) << ",targ atom " <<  pose.residue_type( target_pos ).atom_name( rv_info.residue_atoms()[i] ) << " and targ pos " << target_pos << "base atom " << pose.residue_type( target_pos ).atom_name( rv_info.residue_base_atoms()[i] ) << std::endl;
			}

			if( rv_info.loop_dih() ){
				csts_this_pair.push_back( new DihedralConstraint( core::id::AtomID( rv_info.loopres_base2_atoms()[j], loopres ), core::id::AtomID( rv_info.loopres_base_atoms()[j], loopres ), core::id::AtomID( rv_info.loopres_atoms()[j], loopres ), core::id::AtomID( rv_info.residue_atoms()[i], target_pos ), rv_info.loop_dih() ) );

				//tr << "creating loop dih constraints between looppos " << loopres << " base2 atom " << pose.residue_type( loopres ).atom_name( rv_info.loopres_base2_atoms()[j] )<<  ", base atom " << pose.residue_type( loopres ).atom_name( rv_info.loopres_base_atoms()[j] ) << ", atom " <<  pose.residue_type( loopres ).atom_name( rv_info.loopres_atoms()[j] ) << " and targ pos " << target_pos << " atom " << pose.residue_type( target_pos ).atom_name( rv_info.residue_atoms()[i] ) << std::endl;
			}

			if( rv_info.targ_dih() ){
				csts_this_pair.push_back( new DihedralConstraint( core::id::AtomID( rv_info.loopres_atoms()[j], loopres ), core::id::AtomID( rv_info.residue_atoms()[i], target_pos ), core::id::AtomID( rv_info.residue_base_atoms()[i], target_pos ), core::id::AtomID( rv_info.residue_base2_atoms()[i], target_pos ), rv_info.targ_dih() ) );
				//tr << "creating targ dih constraints between looppos " << loopres << " atom " << pose.residue_type( loopres ).atom_name( rv_info.loopres_atoms()[j] ) << ",targ atom " <<  pose.residue_type( target_pos ).atom_name( rv_info.residue_atoms()[i] ) << " and targ pos " << target_pos << " base atom " << pose.residue_type( target_pos ).atom_name( rv_info.residue_base_atoms()[i] ) << ", base 2 atom " << pose.residue_type( target_pos ).atom_name( rv_info.residue_base2_atoms()[i] )<< std::endl;
			}

			if( rv_info.lt_dih() ){
				csts_this_pair.push_back( new DihedralConstraint( core::id::AtomID( rv_info.loopres_base_atoms()[j], loopres ), core::id::AtomID( rv_info.loopres_atoms()[j], loopres ), core::id::AtomID( rv_info.residue_atoms()[i], target_pos ), core::id::AtomID( rv_info.residue_base_atoms()[i], target_pos ), rv_info.lt_dih() ) );
			}

			if( csts_this_pair.size() == 1 ) csts.push_back( csts_this_pair[1] );
			else if( csts_this_pair.size() > 1 ) csts.push_back( new MultiConstraint( csts_this_pair ) );

			//debug
			//lil stupid: we only need the pose in this function for the debug info
			//but i don't feel like changing headers for debug at the moment,
			//so let's do some call to a pose function to prevent the compiler warning
			pose.total_residue();

		} //loop over loopres residue atoms

	} //loop over targ residue atoms

}

void
ResidueVicinityRCG::lstart( core::Size const lstart )
{
	lstart_ = lstart;
}

void
ResidueVicinityRCG::lstop( core::Size const lstop )
{
	lstop_ = lstop;
}

void
ResidueVicinityRCG::set_rv_infos( utility::vector1< ResidueVicinityInfoOP > const & rv_infos )
{
	rv_infos_ = rv_infos;
}


} //namespace remodel
} //namespace forge
} //namespace protocols
