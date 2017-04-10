// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Brandon Frenz
/// @author Frank DiMaio

#include <protocols/loop_grower/SheetSampler.hh>

#include <iostream>
#include <fstream>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <map>

#include <core/scoring/ScoreFunction.hh>
#include <core/sequence/Sequence.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>
#include <numeric/random/WeightedSampler.hh>

#include <basic/Tracer.hh>

#include <boost/unordered/unordered_map.hpp>
#include <core/id/SequenceMapping.hh>

//possibily duplicate includes here
#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/fragment/JumpingFrame.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/dssp/StrandPairing.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>

#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/RandomSheetBuilder.hh>
#include <protocols/jumping/SheetBuilder.hh>



namespace protocols {
namespace loop_grower {

using namespace core;
using namespace protocols::moves;
using namespace protocols::loops;

void
SheetSampler::apply(core::pose::Pose & pose){
		runtime_assert( start_-end_ >=3);

		//get starting score
		Real bestscore = (*sf_)(pose);
		core::pose::Pose bestpose = pose;
		Size total_residues = pose.total_residue();

		// add residues
		core::Size insert_baseL = start_ + (end_ - start_)/2;
		core::Size insert_baseR = insert_baseL-1;

		std::string chemicaltype;
		if(pose.is_fullatom()){
			chemicaltype = core::chemical::FA_STANDARD;
		}else if(pose.is_centroid()){
			chemicaltype = core::chemical::CENTROID;
		}else{
			chemicaltype = core::chemical::CENTROID_ROT;
		}
		core::chemical::ResidueType const &ala_type = core::chemical::ChemicalManager::get_instance()->residue_type_set( chemicaltype )->name_map("ALA");
		core::conformation::ResidueOP newres( core::conformation::ResidueFactory::create_residue(ala_type) );

		//disable sheet growing of strands that have a proline anchor
		bool build_left = true;
		if( pose.residue(insert_baseL).name3() == "PRO" || pose.residue(insert_baseL+2).name3() == "PRO"){
				build_left = false;
		}
		bool build_right = true;
		if( pose.residue(insert_baseR).name3() == "PRO" || pose.residue(insert_baseR+2).name3() == "PRO" ) {
				build_right = false;
		}

		// strand(left)
		if( build_left ){
				pose.conformation().append_residue_by_jump(*newres, insert_baseL );
				Size strand_start=pose.total_residue();
				Size strand_ref = strand_start;

				//this code assumes append residue by jump always appends the residue at the very end. This should be verified.
				for(Size i=end_; i>insert_baseL; i--){
						pose.conformation().safely_prepend_polymer_residue_before_seqpos(*newres, strand_start, true);
						strand_ref++;
				}
				for(Size i=insert_baseL-1; i>=start_; i--){
						pose.conformation().safely_append_polymer_residue_after_seqpos(*newres, pose.total_residue(), true);
				}
				Size strand_stop=pose.total_residue();

				for(Size j=strand_start; j<=strand_stop; j++){
						pose.set_phi( j, -140.0 );
						pose.set_psi( j, 135.0 );
						pose.set_omega( j, 180.0 );
				}

				// align this to insert-base
				if(ideal_sheets_){
						Size sheet_size = strand_stop-strand_start+1;
						alignStrand( pose, insert_baseL, strand_ref, strand_start, sheet_size);
				}
		
		}
		Real afterfirst = (*sf_)(pose);

		// strand(right)
		if( build_right ) {
				pose.conformation().append_residue_by_jump(*newres, insert_baseR );
				Size strand_start=pose.total_residue();
				Size strand_ref = strand_start;

				// align this to insert-base
				alignPerfectCA( pose, pose.total_residue(), insert_baseR);

				for(Size i=end_; i>insert_baseR; i--){
						pose.conformation().safely_prepend_polymer_residue_before_seqpos(*newres, strand_start, true);
						strand_ref++;
				}
				for(Size i=insert_baseR-1; i>=start_; i--){
						pose.conformation().safely_append_polymer_residue_after_seqpos(*newres, pose.total_residue(), true);
				}
				Size strand_stop=pose.total_residue();

				for(Size j=strand_start; j<=strand_stop; j++){
						pose.set_phi( j, -140.0 );
						pose.set_psi( j, 135.0 );
						pose.set_omega( j, 180.0 );
				}
				// align this to insert-base
				if(ideal_sheets_){
						Size sheet_size = strand_stop-strand_start+1;
						alignStrand( pose, insert_baseR, strand_ref, strand_start, sheet_size);
				}
		}
		Real aftersecond = (*sf_)(pose);
		//Delete strands that clash
		Size totalnew = end_-start_+1;
		bool keep_first = true;
		bool keep_second = true;
		if (bestscore + clashtolerance_ < afterfirst ) keep_first = false;
		if( afterfirst + clashtolerance_ < aftersecond ) keep_second = false;
		
		TR << " start " << bestscore << " after first " << afterfirst << " " << aftersecond << " ct " << clashtolerance_ << std::endl;
		if( !keep_first and !keep_second ){
				pose.conformation().delete_residue_range_slow(total_residues+1, pose.total_residue());
		}else if ( !keep_first || ( !keep_second && !build_left )  ){
				pose.conformation().delete_residue_range_slow(total_residues+1, total_residues+totalnew);
		}else if ( !keep_second && build_left ){
				pose.conformation().delete_residue_range_slow(total_residues+totalnew+1, pose.total_residue());
		}

}

void
SheetSampler::alignPerfectCA( core::pose::Pose & pose, Size moving, Size ref ) {
	numeric::xyzMatrix< core::Real > R;
	numeric::xyzVector< core::Real > preT(0,0,0), postT(0,0,0);

	ObjexxFCL::FArray2D< core::Real > final_coords( 3, 8 );
	ObjexxFCL::FArray2D< core::Real > init_coords( 3, 8 );

	numeric::xyzVector< core::Real > w_1 = pose.residue(moving).atom(" C  ").xyz();
	numeric::xyzVector< core::Real > w_2 = pose.residue(moving).atom(" O  ").xyz();
	numeric::xyzVector< core::Real > w_3 = pose.residue(moving).atom(" N  ").xyz();
	numeric::xyzVector< core::Real > w_4 = pose.residue(moving).atom(" H  ").xyz();
	preT = 0.25*( w_1+w_2+w_3+w_4 );

	numeric::xyzVector< core::Real > x_1 = pose.residue(ref).atom(" C  ").xyz();
	numeric::xyzVector< core::Real > x_2 = pose.residue(ref).atom(" O  ").xyz();
	numeric::xyzVector< core::Real > x_3 = pose.residue(ref).atom(" N  ").xyz();
	numeric::xyzVector< core::Real > x_4 = pose.residue(ref).atom(" H  ").xyz();
	numeric::xyzVector< core::Real > x_12 = (x_2-x_1).normalize();
	numeric::xyzVector< core::Real > x_34 = (x_4-x_3).normalize();
	numeric::xyzVector< core::Real > y_1 = x_4 + 2.9*x_34;
	numeric::xyzVector< core::Real > y_2 = x_4 + 1.9*x_34;
	numeric::xyzVector< core::Real > y_3 = x_2 + 2.9*x_12;
	numeric::xyzVector< core::Real > y_4 = x_2 + 1.9*x_12;
	postT = 0.25*( y_1+y_2+y_3+y_4 );

	for ( int j=0; j<3; ++j ) {
		init_coords(j+1,1) = w_1[j];
		init_coords(j+1,2) = w_2[j];
		init_coords(j+1,3) = w_3[j];
		init_coords(j+1,4) = w_4[j];
		final_coords(j+1,1) = y_1[j];
		final_coords(j+1,2) = y_2[j];
		final_coords(j+1,3) = y_3[j];
		final_coords(j+1,4) = y_4[j];
	}

	for ( int i=1; i<=4; ++i ) {
		for ( int j=0; j<3; ++j ) {
			init_coords(j+1,i) -= preT[j];
			final_coords(j+1,i) -= postT[j];
		}
	}

	// get optimal superposition
	// rotate >init< to >final<
	ObjexxFCL::FArray1D< numeric::Real > ww( 8, 1.0 );
	ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::Real ctx;

	numeric::model_quality::findUU( init_coords, final_coords, ww, 4, uu, ctx );
	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );

	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector<core::Real> > coords;
	for ( int j=1; j<=(int)pose.residue(moving).natoms(); ++j ) {
		core::id::AtomID src(j,moving);
		ids.push_back(src);
		coords.push_back( postT + (R*(pose.xyz( src )-preT)) );
	}
	pose.batch_set_xyz(ids, coords);
}

void 
printxyz(numeric::xyzVector< core::Real > printvec ){
		TR << printvec[0] << " " << printvec[1] << " " << printvec[2] << std::endl;;
}

void
SheetSampler::alignStrand( core::pose::Pose & pose, Size ref, Size moving, Size strand_start, Size strand_size ) {

	Size second_moving = moving-2;
	Size second_ref = ref+2;
	//TR << " ref second ref " << ref << " " << second_ref << std::endl << " moving second moving " << moving << " " << second_moving << std::endl;
	numeric::xyzMatrix< core::Real > R;
	numeric::xyzVector< core::Real > preT(0,0,0), postT(0,0,0);

	ObjexxFCL::FArray2D< core::Real > final_coords( 3, 8 );
	ObjexxFCL::FArray2D< core::Real > init_coords( 3, 8 );

	numeric::xyzVector< core::Real > w_1 = pose.residue(moving).atom(" C  ").xyz();
	numeric::xyzVector< core::Real > w_2 = pose.residue(moving).atom(" O  ").xyz();
	numeric::xyzVector< core::Real > w_3 = pose.residue(moving).atom(" N  ").xyz();
	numeric::xyzVector< core::Real > w_4 = pose.residue(moving).atom(" H  ").xyz();
	numeric::xyzVector< core::Real > w2_1 = pose.residue(second_moving).atom(" C  ").xyz();
	numeric::xyzVector< core::Real > w2_2 = pose.residue(second_moving).atom(" O  ").xyz();
	numeric::xyzVector< core::Real > w2_3 = pose.residue(second_moving).atom(" N  ").xyz();
	numeric::xyzVector< core::Real > w2_4 = pose.residue(second_moving).atom(" H  ").xyz();
	preT = 0.125*( w_1+w_2+w_3+w_4+w2_1+w2_2+w2_3+w2_4 );

	//first ref
	numeric::xyzVector< core::Real > x_1 = pose.residue(ref).atom(" C  ").xyz();
	numeric::xyzVector< core::Real > x_2 = pose.residue(ref).atom(" O  ").xyz();
	numeric::xyzVector< core::Real > x_3 = pose.residue(ref).atom(" N  ").xyz();
	numeric::xyzVector< core::Real > x_4 = pose.residue(ref).atom(" H  ").xyz();
	numeric::xyzVector< core::Real > x_12 = (x_2-x_1).normalize();
	numeric::xyzVector< core::Real > x_34 = (x_4-x_3).normalize();
	numeric::xyzVector< core::Real > y_1 = x_4 + 2.9*x_34;
	numeric::xyzVector< core::Real > y_2 = x_4 + 1.9*x_34;
	numeric::xyzVector< core::Real > y_3 = x_2 + 2.9*x_12;
	numeric::xyzVector< core::Real > y_4 = x_2 + 1.9*x_12;
	//secondref
	numeric::xyzVector< core::Real > x2_1 = pose.residue(second_ref).atom(" C  ").xyz();
	numeric::xyzVector< core::Real > x2_2 = pose.residue(second_ref).atom(" O  ").xyz();
	numeric::xyzVector< core::Real > x2_3 = pose.residue(second_ref).atom(" N  ").xyz();
	numeric::xyzVector< core::Real > x2_4 = pose.residue(second_ref).atom(" H  ").xyz();
	numeric::xyzVector< core::Real > x2_12 = (x2_2-x2_1).normalize();
	numeric::xyzVector< core::Real > x2_34 = (x2_4-x2_3).normalize();
	numeric::xyzVector< core::Real > y2_1 = x2_4 + 2.9*x2_34;
	numeric::xyzVector< core::Real > y2_2 = x2_4 + 1.9*x2_34;
	numeric::xyzVector< core::Real > y2_3 = x2_2 + 2.9*x2_12;
	numeric::xyzVector< core::Real > y2_4 = x2_2 + 1.9*x2_12;
	/*TR << "second C ";
	printxyz(y2_1);
	TR << "second O ";
	printxyz(y2_2);
	TR << "second N ";
	printxyz(y2_3);
	TR << " second H ";
	printxyz(y2_4);*/
	
	//center of mass
	postT = 0.125*( y_1+y_2+y_3+y_4+y2_1+y2_2+y2_3+y2_4 );

	//write coords to disk
	/*std::ofstream outbeam;
	outbeam.open("alignmentpoints.txt", std::ofstream::app);
	outbeam << "y1 1 " << y_1[0] << " " << y_1[1] << " " << y_1[2] << std::endl;
	outbeam << "y1 2 " << y_2[0] << " " << y_2[1] << " " << y_2[2] << std::endl;
	outbeam << "y1 3 " << y_3[0] << " " << y_3[1] << " " << y_3[2] << std::endl;
	outbeam << "y1 4 " << y_4[0] << " " << y_4[1] << " " << y_4[2] << std::endl;
	outbeam << "y2 1 " << y2_1[0] << " " << y2_1[1] << " " << y2_1[2] << std::endl;
	outbeam << "y2 2 " << y2_2[0] << " " << y2_2[1] << " " << y2_2[2] << std::endl;
	outbeam << "y2 3 " << y2_3[0] << " " << y2_3[1] << " " << y2_3[2] << std::endl;
	outbeam << "y2 4 " << y2_4[0] << " " << y2_4[1] << " " << y2_4[2] << std::endl;
	outbeam.close();*/

	//fill the array
	for ( int j=0; j<3; ++j ) {
		init_coords(j+1,1) = w_1[j];
		init_coords(j+1,2) = w_2[j];
		init_coords(j+1,3) = w_3[j];
		init_coords(j+1,4) = w_4[j];
		init_coords(j+1,5) = w2_1[j];
		init_coords(j+1,6) = w2_2[j];
		init_coords(j+1,7) = w2_3[j];
		init_coords(j+1,8) = w2_4[j];
		
		final_coords(j+1,1) = y_1[j];
		final_coords(j+1,2) = y_2[j];
		final_coords(j+1,3) = y_3[j];
		final_coords(j+1,4) = y_4[j];
		final_coords(j+1,5) = y2_1[j];
		final_coords(j+1,6) = y2_2[j];
		final_coords(j+1,7) = y2_3[j];
		final_coords(j+1,8) = y2_4[j];
	}

	for ( int i=1; i<=8; ++i ) {
		for ( int j=0; j<3; ++j ) {
			init_coords(j+1,i) -= preT[j];
			final_coords(j+1,i) -= postT[j];
		}
	}

	// get optimal superposition
	// rotate >init< to >final<
	ObjexxFCL::FArray1D< numeric::Real > ww( 8, 1.0 );
	ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
	numeric::Real ctx;

	numeric::model_quality::findUU( init_coords, final_coords, ww, 8, uu, ctx );
	R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
	R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
	R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );

	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector<core::Real> > coords;
	for( Size i=strand_start; i<=strand_start+strand_size-1; i++){
		for ( int j=1; j<=(int)pose.residue(i).natoms(); ++j ) {
			core::id::AtomID src(j,i);
			ids.push_back(src);
			coords.push_back( postT + (R*(pose.xyz( src )-preT)) );
		}
	}
	pose.batch_set_xyz(ids, coords);
}

Real 
SheetSampler::sheethbonds(core::pose::Pose& pose, Size lower, Size upper){
	
		Real hbond_energies = 0;	
		//populate hbond set:
		core::scoring::hbonds::HBondSet set1;
		pose.update_residue_neighbors();
		set1.setup_for_residue_pair_energies( pose, false, false );
		
		
		//query hbond set, get energies:
		for (Size i = 1; i<= set1.nhbonds(); i++) {
		
		  core::scoring::hbonds::HBond bond = set1.hbond( i );
		  Size accResNum = bond.acc_res();
		  Size donResNum = bond.don_res();
		  //get acc and donor residues from sequence numbers
		  core::conformation::Residue accRes = pose.residue( accResNum );
		  core::conformation::Residue donRes = pose.residue( donResNum );
			bool accinrange = false;
			bool doninrange = false;
			if( accResNum > lower && accResNum <= upper ) accinrange = true;
			if( donResNum > lower && donResNum <= upper ) doninrange = true;
			if( doninrange && accinrange) continue;
			if( accinrange || doninrange ){
		  		Real energy = bond.energy();
		  		Real weight = bond.weight();
		  		hbond_energies +=  weight*energy;
			}
		}
		return hbond_energies;

}

} //loop_grower		
} //protocols
