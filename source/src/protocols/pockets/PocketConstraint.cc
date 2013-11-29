// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

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
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <basic/Tracer.hh>
#include <string>
#include <ObjexxFCL/string.functions.hh>
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/out.OptionKeys.gen.hh>
#ifndef _WIN32
// AUTO-REMOVED #include <sys/time.h>

#include <utility/vector1.hh>

#endif

namespace protocols {
namespace pockets {

PocketConstraintCreator::PocketConstraintCreator() {}
PocketConstraintCreator::~PocketConstraintCreator() {}

core::scoring::constraints::ConstraintOP
PocketConstraintCreator::create_constraint() const {
	return new PocketConstraint;
}

std::string PocketConstraintCreator::keyname() const
{
	return "Pocket";
}


static basic::Tracer TR("core.scoring.constraints.PocketConstraint");

void PocketConstraint::init(core::pose::Pose const & pose){
	using namespace basic::options;
	seqpos_ = 0;
	weight_ = option[ OptionKeys::constraints::pocket_constraint_weight ]();
	dumppdb_=option[ OptionKeys::pocket_grid::pocket_dump_pdbs ]();
	totalres_=pose.total_residue();
	if (option[ OptionKeys::pocket_grid::pocket_num_angles ] <1){
		std::cout<<"PocketConstraint: invalid number of angles specified.  Exiting."<<std::endl;
		exit(999);
	} else {
		angles_ = option[ OptionKeys::pocket_grid::pocket_num_angles ];
	}

}

void PocketConstraint::read_def(
  std::istream & line_stream,
  core::pose::Pose const & pose,
  core::scoring::constraints::FuncFactory const & /* func_factory */)
{
	init(pose);
	std::string tmp;
	std::string resid("");
	if ((line_stream >> weight_>>resid )){
		residues_ = protocols::pockets::PocketGrid::getRelaxResidues(pose, resid);
		if ( residues_.size() == 0 ) {
			std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
			exit(1);
		}
	}
	else{
		std::cout << "ERROR!! Invalid PocketConstraint specification" << std::endl;
		exit(1);
	}
	pocketgrid_ = new protocols::pockets::PocketGrid(residues_);

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
		if (fpos != resid.size()-1 ) {
			chain = resid[ fpos+1 ];
		}
	} else {
		central_relax_pdb_number = ObjexxFCL::int_of( resid );
	}

	for ( int j = 1, resnum = pose.total_residue(); j <= resnum; ++j ) {
		if ( pose.pdb_info()->number(j) == central_relax_pdb_number ) {
			//seqpos_ = j;
			if (chain != ' '){
				if ( pose.pdb_info()->chain(j) == chain ) {
					seqpos_ = j;
				}
			}else{
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
		pocketgrid_ = new protocols::pockets::PocketGrid( pose.conformation().residue(seqpos_) );
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


PocketConstraint::~PocketConstraint() {}


void PocketConstraint::set_target_res( core::pose::Pose const & pose, Size new_seqpos ){
	if (new_seqpos>pose.total_residue()){
		std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
		exit(1);
	}
	if ( seqpos_ != 0 ) {
		seqpos_=new_seqpos;
		pocketgrid_ = new protocols::pockets::PocketGrid( pose.conformation().residue(seqpos_) );
	}else{
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
		if (fpos != resid.size()-1 ) {
			chain = resid[ fpos+1 ];
		}
	} else {
		central_relax_pdb_number = ObjexxFCL::int_of( resid );
	}

	seqpos_ = 0;
	for ( int j = 1, resnum = pose.total_residue(); j <= resnum; ++j ) {
		if ( pose.pdb_info()->number(j) == central_relax_pdb_number ) {
			//seqpos_ = j;
			if (chain != ' '){
				if ( pose.pdb_info()->chain(j) == chain ) {
					seqpos_ = j;
				}
			}else{
				seqpos_ = j;
			}
		}
	}

	if ( seqpos_ != 0 ) {
		pocketgrid_ = new protocols::pockets::PocketGrid( pose.conformation().residue(seqpos_) );
	} else {
		std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
		exit(1);
	}

}




// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
PocketConstraint::score( core::scoring::constraints::XYZ_Func const & xyz_func, core::scoring::EnergyMap const & weights, core::scoring::EnergyMap & emap ) const
{
  using namespace basic::options;
  bool debug = option[ OptionKeys::pocket_grid::pocket_dump_pdbs ]();
	//std::cout<< "hi\n";
	//TR<<"hi\n";
	if ( weights[ this->score_type() ] == 0 ) return;
	if (seqpos_==0 && residues_.size()==0){
		std::cout << "ERROR!! Invalid residue to backrub around" << std::endl;
		exit(1);
	}


	core::Real cst_avg = 0;
	core::Real largestPocketVol;

  if (debug)  TR<<"Pocket Volumes: ";

	for (core::Size angleCount=0; angleCount < (angles_ - 1); ++angleCount){
		pocketgrid_ -> randomAngle();
		if (seqpos_ != 0){
			core::conformation::Residue const & curr_rsd ( xyz_func.residue(seqpos_) );
			pocketgrid_->autoexpanding_pocket_eval( curr_rsd, xyz_func, totalres_ );
		}else{
			pocketgrid_->autoexpanding_pocket_eval( residues_, xyz_func, totalres_ );
		}
		core::Real largestPocketVol=pocketgrid_->netTargetPocketVolume();
		cst_avg += largestPocketVol;
    if (debug) TR<<largestPocketVol<<" ";
	}

	pocketgrid_ -> zeroAngle();
	if (seqpos_ != 0){
		core::conformation::Residue const & curr_rsd ( xyz_func.residue(seqpos_) );
		pocketgrid_->autoexpanding_pocket_eval( curr_rsd, xyz_func, totalres_ );
	}else{
		pocketgrid_->autoexpanding_pocket_eval( residues_, xyz_func, totalres_ );
	}
	core::Real cst_val = -1.;

	//  core::Real largestPocketVol=pocketgrid_->netTargetPocketVolume();

	if (dumppdb_) pocketgrid_->dumpGridToFile();

	//	core::Real vol=pocketgrid_->targetPocketVolume(surf_score, bur_score);
	//	core::Real sa=pocketgrid_->targetPocketSolventSurface();
	//	core::Real psa=pocketgrid_->targetPocketProteinSurface();
	//	core::Real hpsa=pocketgrid_->targetPocketHydrophobicProteinSurface();
	//	core::Real ppsa=pocketgrid_->targetPocketPolarProteinSurface();
	//	core::Real nps=pocketgrid_->targetPocketHeuristicScore();

	//core::Real largestPocketVol=pocketgrid_->largestTargetPocketVolume();
	largestPocketVol=pocketgrid_->netTargetPocketVolume();
	//core::Real largestPocketVol=0;

	cst_avg += largestPocketVol;
	if (debug) TR<<largestPocketVol<<" ";
  cst_avg /= angles_;
  if (debug) TR<<"Average: "<<cst_avg<<std::endl;
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
	core::scoring::constraints::XYZ_Func const & ,
	core::Vector & ,
	core::Vector & ,
	core::scoring::EnergyMap const & weights
) const
{

	if ( weights[ this->score_type() ] == 0 ) return;

	using namespace basic::options;
	if (!option[ OptionKeys::constraints::pocket_zero_derivatives ]()){
		TR << "ERROR - derivatives not yet implemented for PocketConstraints." << std::endl;
		std::exit(1);
	}

	return;

}


core::scoring::constraints::ConstraintOP PocketConstraint::clone() const {
	return core::scoring::constraints::ConstraintOP( new PocketConstraint( *this ) );
}


} // namespace constraints_additional
} // namespace protocols
