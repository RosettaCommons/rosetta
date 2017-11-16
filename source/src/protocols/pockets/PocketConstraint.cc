// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pockets/PocketConstraint.cc
///
/// @brief
/// @author David Johnson


#include <protocols/pockets/PocketConstraint.hh>
#include <protocols/pockets/PocketConstraintCreator.hh>
#include <protocols/pockets/PocketGrid.hh>


#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>
#include <string>
#include <ObjexxFCL/string.functions.hh>
#include <iostream>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#ifndef _WIN32

#include <utility/vector1.hh>

#endif

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/vector.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace pockets {

PocketConstraintCreator::PocketConstraintCreator() {}
PocketConstraintCreator::~PocketConstraintCreator() = default;

core::scoring::constraints::ConstraintOP
PocketConstraintCreator::create_constraint() const {
	return core::scoring::constraints::ConstraintOP( new PocketConstraint );
}

std::string PocketConstraintCreator::keyname() const
{
	return "Pocket";
}


static basic::Tracer TR( "core.scoring.constraints.PocketConstraint" );

void PocketConstraint::init(core::pose::Pose const & pose){
	using namespace basic::options;
	seqpos_ = 0;
	weight_ = option[ OptionKeys::constraints::pocket_constraint_weight ]();
	dumppdb_=option[ OptionKeys::pocket_grid::pocket_dump_pdbs ]();
	totalres_=pose.size();
	if ( option[ OptionKeys::pocket_grid::pocket_num_angles ] <1 ) {
		std::cout<<"PocketConstraint: invalid number of angles specified.  Exiting."<<std::endl;
		exit(999);
	} else {
		angles_ = option[ OptionKeys::pocket_grid::pocket_num_angles ];
	}

}

void PocketConstraint::read_def(
	std::istream & line_stream,
	core::pose::Pose const & pose,
	core::scoring::func::FuncFactory const & /* func_factory */)
{
	init(pose);
	//std::string tmp;
	std::string resid("");
	if ( (line_stream >> weight_>>resid ) ) {
		residues_ = protocols::pockets::PocketGrid::getRelaxResidues(pose, resid);
		if ( residues_.size() == 0 ) {
			std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
			exit(1);
		}
	} else {
		std::cout << "ERROR!! Invalid PocketConstraint specification" << std::endl;
		exit(1);
	}
	pocketgrid_ = protocols::pockets::PocketGridOP( new protocols::pockets::PocketGrid(residues_) );

}

void PocketConstraint::show_def( std::ostream&  out , core::pose::Pose const& /* pose */ ) const {
	out << "PocketConstraint::show_def() " << std::endl;
}


PocketConstraint::PocketConstraint():Constraint( core::scoring::pocket_constraint ){}

PocketConstraint::PocketConstraint(
	core::pose::Pose const & pose
):
	Constraint( core::scoring::pocket_constraint )
{
	using namespace basic::options;
	// for now, set the constraint to depend on ALL atom positions, ie. if ANYTHING moves we have to update the constraint
	// later, we could pre-define residues near the selected residue and make the constraint depend only on these
	// This is the residue we'll backrub around!!
	//int const central_relax_pdb_number = option[ OptionKeys::pocket_grid::central_relax_pdb_num ];
	init(pose);

	std::string resid(option[ OptionKeys::pocket_grid::central_relax_pdb_num ]);
	int  central_relax_pdb_number;
	char chain = ' ';
	std::size_t fpos( resid.find(':') );
	if ( fpos != std::string::npos ) {
		central_relax_pdb_number = ObjexxFCL::int_of( resid.substr(0,fpos) );
		if ( fpos != resid.size()-1 ) {
			chain = resid[ fpos+1 ];
		}
	} else {
		central_relax_pdb_number = ObjexxFCL::int_of( resid );
	}

	for ( int j = 1, resnum = pose.size(); j <= resnum; ++j ) {
		if ( pose.pdb_info()->number(j) == central_relax_pdb_number ) {
			//seqpos_ = j;
			if ( chain != ' ' ) {
				if ( pose.pdb_info()->chain(j) == chain ) {
					seqpos_ = j;
				}
			} else {
				seqpos_ = j;
			}
		}
	}

	//  Do not crash yet;
	//      if ( seqpos_ == 0 ) {
	//      std::cout << "ERROR!! Could not find residue to backrub around" << std::endl;
	//      exit(1);
	//      }

	if ( seqpos_ != 0 ) {
		pocketgrid_ = protocols::pockets::PocketGridOP( new protocols::pockets::PocketGrid( pose.conformation().residue(seqpos_) ) );
	}

	// JK NOTE: WE'RE NOT USING THE "FUNC" SYSTEM, THIS COULD BE ADDED LATER....

}

PocketConstraint::PocketConstraint( const PocketConstraint& old ):
	Constraint( core::scoring::pocket_constraint )
{
	seqpos_ = old.seqpos_;
	totalres_ = old.totalres_;
	pocketgrid_ = old.pocketgrid_;
	dumppdb_ = old.dumppdb_;
}


PocketConstraint::~PocketConstraint() = default;


void PocketConstraint::set_target_res( core::pose::Pose const & pose, Size new_seqpos ){
	if ( new_seqpos>pose.size() ) {
		std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
		exit(1);
	}
	if ( seqpos_ != 0 ) {
		seqpos_=new_seqpos;
		pocketgrid_ = protocols::pockets::PocketGridOP( new protocols::pockets::PocketGrid( pose.conformation().residue(seqpos_) ) );
	} else {
		std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
		exit(1);
	}
}


void PocketConstraint::set_target_res_pdb( core::pose::Pose const & pose, std::string resid ){
	//std::cout<<size_x_<<" "<<size_y_<<" "<<size_z_<<"\n";
	int  central_relax_pdb_number;
	char chain = ' ';
	std::size_t fpos( resid.find(':') );
	if ( fpos != std::string::npos ) {
		central_relax_pdb_number = ObjexxFCL::int_of( resid.substr(0,fpos) );
		if ( fpos != resid.size()-1 ) {
			chain = resid[ fpos+1 ];
		}
	} else {
		central_relax_pdb_number = ObjexxFCL::int_of( resid );
	}

	seqpos_ = 0;
	for ( int j = 1, resnum = pose.size(); j <= resnum; ++j ) {
		if ( pose.pdb_info()->number(j) == central_relax_pdb_number ) {
			//seqpos_ = j;
			if ( chain != ' ' ) {
				if ( pose.pdb_info()->chain(j) == chain ) {
					seqpos_ = j;
				}
			} else {
				seqpos_ = j;
			}
		}
	}

	if ( seqpos_ != 0 ) {
		pocketgrid_ = protocols::pockets::PocketGridOP( new protocols::pockets::PocketGrid( pose.conformation().residue(seqpos_) ) );
	} else {
		std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
		exit(1);
	}

}


// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
PocketConstraint::score( core::scoring::func::XYZ_Func const & xyz_func, core::scoring::EnergyMap const & weights, core::scoring::EnergyMap & emap ) const
{
	using namespace basic::options;
	bool debug = option[ OptionKeys::pocket_grid::pocket_dump_pdbs ]();
	//std::cout<< "hi\n";
	//TR<<"hi\n";
	if ( weights[ this->score_type() ] == 0 ) return;
	if ( seqpos_==0 && residues_.size()==0 ) {
		std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
		exit(1);
	}


	core::Real cst_avg = 0;
	core::Real largestPocketVol;

	if ( debug )  TR<<"Pocket Volumes: ";

	for ( core::Size angleCount=0; angleCount < (angles_ - 1); ++angleCount ) {
		pocketgrid_ -> randomAngle();
		if ( seqpos_ != 0 ) {
			core::conformation::Residue const & curr_rsd ( xyz_func.residue(seqpos_) );
			pocketgrid_->autoexpanding_pocket_eval( curr_rsd, xyz_func, totalres_ );
		} else {
			pocketgrid_->autoexpanding_pocket_eval( residues_, xyz_func, totalres_ );
		}
		core::Real largestPocketVol=pocketgrid_->netTargetPocketVolume();
		cst_avg += largestPocketVol;
		if ( debug ) TR<<largestPocketVol<<" ";
	}

	pocketgrid_ -> zeroAngle();
	if ( seqpos_ != 0 ) {
		core::conformation::Residue const & curr_rsd ( xyz_func.residue(seqpos_) );
		pocketgrid_->autoexpanding_pocket_eval( curr_rsd, xyz_func, totalres_ );
	} else {
		pocketgrid_->autoexpanding_pocket_eval( residues_, xyz_func, totalres_ );
	}
	core::Real cst_val = -1.;

	//  core::Real largestPocketVol=pocketgrid_->netTargetPocketVolume();

	if ( dumppdb_ ) pocketgrid_->dumpGridToFile();

	// core::Real vol=pocketgrid_->targetPocketVolume(surf_score, bur_score);
	// core::Real sa=pocketgrid_->targetPocketSolventSurface();
	// core::Real psa=pocketgrid_->targetPocketProteinSurface();
	// core::Real hpsa=pocketgrid_->targetPocketHydrophobicProteinSurface();
	// core::Real ppsa=pocketgrid_->targetPocketPolarProteinSurface();
	// core::Real nps=pocketgrid_->targetPocketHeuristicScore();

	//core::Real largestPocketVol=pocketgrid_->largestTargetPocketVolume();
	largestPocketVol=pocketgrid_->netTargetPocketVolume();
	//core::Real largestPocketVol=0;

	cst_avg += largestPocketVol;
	if ( debug ) TR<<largestPocketVol<<" ";
	cst_avg /= angles_;
	if ( debug ) TR<<"Average: "<<cst_avg<<std::endl;
	cst_val *= (cst_avg);
	//cst_val *= (vol);
	//cst_val *= (vol*vol/sa);
	//  std::cout<<"Vol: "<<vol<<" Solvent surface: "<<sa<<" Prot surface: "<<psa<<" Hprot surface: "<<hpsa<<" Pprot surface: "<<ppsa<<" New Score: "<<nps<<std::endl;
	//std::cout<<"done4\n";
	//  cst_val *= distance(CB_curr,CA_curr);
	emap[ this->score_type() ] += cst_val*weight_;
	//std::cout<<cst_val<<" done5\n";

}


void
PocketConstraint::fill_f1_f2(
	core::id::AtomID const & ,
	core::scoring::func::XYZ_Func const & ,
	core::Vector & ,
	core::Vector & ,
	core::scoring::EnergyMap const & weights
) const
{

	if ( weights[ this->score_type() ] == 0 ) return;

	using namespace basic::options;
	if ( !option[ OptionKeys::constraints::pocket_zero_derivatives ]() ) {
		TR << "ERROR - derivatives not yet implemented for PocketConstraints." << std::endl;
		std::exit(1);
	}

	return;

}


core::scoring::constraints::ConstraintOP PocketConstraint::clone() const {
	PocketConstraintOP myclone( new PocketConstraint( *this ) );
	myclone->pocketgrid_ = PocketGridOP( new PocketGrid( *pocketgrid_ ));
	return myclone;
}

bool PocketConstraint::operator == ( core::scoring::constraints::Constraint const & other ) const
{
	if ( ! same_type_as_me( other ) || ! other.same_type_as_me( *this ) ) return false;
	PocketConstraint const & other_downcast( static_cast< PocketConstraint const & >( other ) );

	if ( seqpos_ != other_downcast.seqpos_ ) return false;
	if ( totalres_ != other_downcast.totalres_ ) return false;
	if ( angles_ != other_downcast.angles_ ) return false;
	if ( weight_     != other_downcast.weight_ ) return false;
	if ( ! ( pocketgrid_ == other_downcast.pocketgrid_ || ( pocketgrid_ && other_downcast.pocketgrid_ && *pocketgrid_ == *other_downcast.pocketgrid_ )) ) return false;
	if ( atom_ids_   != other_downcast.atom_ids_ ) return false;
	if ( dumppdb_    != other_downcast.dumppdb_ ) return false;
	if ( residues_   != other_downcast.residues_ ) return false;

	return true;
}
bool PocketConstraint::same_type_as_me( core::scoring::constraints::Constraint const & other ) const
{
	return dynamic_cast< PocketConstraint const *  > (&other);
}


} // namespace constraints_additional
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::pockets::PocketConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::constraints::Constraint >( this ) );
	arc( CEREAL_NVP( seqpos_ ) ); // core::Size
	arc( CEREAL_NVP( totalres_ ) ); // core::Size
	arc( CEREAL_NVP( angles_ ) ); // core::Size
	arc( CEREAL_NVP( weight_ ) ); // core::Real
	arc( CEREAL_NVP( pocketgrid_ ) ); // protocols::pockets::PocketGridOP
	arc( CEREAL_NVP( atom_ids_ ) ); // utility::vector1<AtomID>
	arc( CEREAL_NVP( dumppdb_ ) ); // _Bool
	arc( CEREAL_NVP( residues_ ) ); // std::vector<core::conformation::ResidueCOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::pockets::PocketConstraint::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::constraints::Constraint >( this ) );
	arc( seqpos_ ); // core::Size
	arc( totalres_ ); // core::Size
	arc( angles_ ); // core::Size
	arc( weight_ ); // core::Real
	arc( pocketgrid_ ); // protocols::pockets::PocketGridOP
	arc( atom_ids_ ); // utility::vector1<AtomID>
	arc( dumppdb_ ); // _Bool
	std::vector< core::conformation::ResidueOP > local_residues;
	arc( local_residues ); // std::vector<core::conformation::ResidueCOP>
	// booo! too bad this doesn't work residues_ = local_residues; // copy the non-const pointer(s) into the const pointer(s)
	for ( std::vector< core::conformation::ResidueOP >::const_iterator iter = local_residues.begin();
			iter != local_residues.end(); ++iter ) {
		residues_.push_back( *iter );
	}
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::pockets::PocketConstraint );
CEREAL_REGISTER_TYPE( protocols::pockets::PocketConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_pockets_PocketConstraint )
#endif // SERIALIZATION
