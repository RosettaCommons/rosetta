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
/// @author Frank DiMaio
/// @author Srivatsan Raman
#include <protocols/rbsegment_relax/RBSegmentMover.hh>
#include <protocols/rbsegment_relax/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/dssp/Dssp.hh>

// Rosetta Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/conformation/symmetry/util.hh>
#include <core/pose/symmetry/util.hh>

#include <core/chemical/VariantType.hh>
#include <core/kinematics/MoveMap.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Random number generator
#include <numeric/xyzVector.io.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray1D.hh>


#include <string>

//Auto Headers
#include <core/id/NamedAtomID.hh>
#include <core/pose/util.hh>
#include <core/id/SequenceMapping.hh>


namespace protocols {
namespace rbsegment_relax {

using namespace core;

static thread_local basic::Tracer TR( "protocols::moves::RBSegmentMover" );


/// @brief set up constraints from RB segs
///       currently uses bounded constraints on each CA ... make CST type flag-selectable????
void set_rb_constraints(
	core::pose::Pose & pose,
	core::pose::Pose const &cst_pose,
	utility::vector1< protocols::rbsegment_relax::RBSegment > const & rbsegs ,
	core::id::SequenceMapping const & resmap,
	core::Real cst_width,
	core::Real cst_stdev,
	core::Size cst_seqwidth
                   ) {

	using core::scoring::constraints::ConstraintOP;
	using core::scoring::constraints::ConstraintCOP;

	core::scoring::constraints::ConstraintSetOP new_csts( new core::scoring::constraints::ConstraintSet );

	for (int i =  1; i <= (int)rbsegs.size(); ++i) {
		for (core::Size j=1; j<=rbsegs[i].nContinuousSegments(); ++j) {
			int rbstart = rbsegs[i][j].start(), rbend = rbsegs[i][j].end();

			for (int cst_i=rbstart; cst_i<=rbend; ++cst_i) {

				// make ambiguous cst from ( cst_i-cst_seqwidth , cst_i+cst_seqwidth )
				core::scoring::constraints::ConstraintCOPs CSTs_i;


				for (int cst_ij = std::max(1,cst_i-(int)cst_seqwidth), cst_ij_end=std::min(cst_pose.total_residue(),cst_i+cst_seqwidth);
				     cst_ij<=cst_ij_end; ++cst_ij) {
					// only constrain protein residues
					if (!cst_pose.residue(cst_ij).is_protein()) continue;

					core::scoring::func::FuncOP fx( new core::scoring::constraints::BoundFunc(0,cst_width,cst_stdev,"xyz") );
					core::scoring::constraints::ConstraintOP newcst_ij( new core::scoring::constraints::CoordinateConstraint(
											 core::id::AtomID( pose.residue(resmap[cst_i]).atom_index("CA"), resmap[cst_i]),
											 core::id::AtomID( pose.residue(pose.total_residue()).atom_index("ORIG"), pose.total_residue()),
											 cst_pose.residue(cst_ij).xyz( "CA" ),
											  fx ) );


					// make sure not randomized <<<< kind of hacky
					core::Real len = (cst_pose.residue(cst_ij).xyz( "CA" )).length();
					if ( len < 500 )
						CSTs_i.push_back( newcst_ij );
				}


				if ( CSTs_i.size() > 0 ) {
					TR << "Adding " << CSTs_i.size() << " ambiguous constraints for res " << resmap[cst_i]
					   << " (= input " << cst_i << ")" << std::endl;
					new_csts->add_constraint( ConstraintCOP( ConstraintOP( new core::scoring::constraints::AmbiguousConstraint( CSTs_i ) ) ) );
				}

			}
		}
	}
	pose.constraint_set( new_csts );
}


/// @brief set up constraints on complete pose (not just RB segments)
///       currently uses bounded constraints ... make CST type flag-selectable????
void set_constraints(
	core::pose::Pose & pose,
	core::pose::Pose const &cst_pose,
	core::Real cst_width,
	core::Real cst_stdev,
	core::Size cst_seqwidth
                   ) {

	using core::scoring::constraints::ConstraintOP;
	using core::scoring::constraints::ConstraintCOP;

	if (pose.total_residue() != cst_pose.total_residue()) {
		TR.Warning << "set_constraints() error: #res in cst_pose (" << cst_pose.total_residue()
		            << ") != #res in pose (" << pose.total_residue() << ").  Continuing..." << std::endl;
	}
	int nres = (int)std::min( pose.total_residue() , cst_pose.total_residue() );

	core::scoring::constraints::ConstraintSetOP new_csts( new core::scoring::constraints::ConstraintSet );
	for (int cst_i =  1; cst_i < nres; ++cst_i) {
		// make ambiguous cst from ( cst_i-cst_seqwidth , cst_i+cst_seqwidth )
		core::scoring::constraints::ConstraintCOPs CSTs_i;

		for (int cst_ij=std::max(1,cst_i-(int)cst_seqwidth), cst_ij_end=std::min(nres,cst_i+(int)cst_seqwidth);
				 cst_ij<=cst_ij_end; ++cst_ij) {
			// only constrain protein residues
			if (!cst_pose.residue(cst_ij).is_protein()) continue;

			core::scoring::func::FuncOP fx( new core::scoring::constraints::BoundFunc(0,cst_width,cst_stdev,"xyz") );
			core::scoring::constraints::ConstraintOP newcst_ij( new core::scoring::constraints::CoordinateConstraint(
				core::id::AtomID( pose.residue(cst_i               ).atom_index("CA"  ), cst_i),
				core::id::AtomID( pose.residue(pose.total_residue()).atom_index("ORIG"), pose.total_residue()),
				cst_pose.residue(cst_ij).xyz( "CA" ),
				fx
			) );

			// make sure not randomized <<<< kind of hacky
			core::Real len = (cst_pose.residue(cst_ij).xyz( "CA" )).length();
			if ( len < 500 )
				CSTs_i.push_back( newcst_ij );
		}


		if ( CSTs_i.size() > 0 ) {
			TR << "Adding " << CSTs_i.size() << " ambiguous constraints for res " << cst_i << std::endl;
			new_csts->add_constraint( ConstraintCOP( ConstraintOP( new core::scoring::constraints::AmbiguousConstraint( CSTs_i ) ) ) );
		}
	}
	pose.constraint_set( new_csts );
}


/// @brief Helper function to set up a pose; unlike alt version keep loops (use cutpoint variants)
///   unlike version in loops_main, this uses RBSegment structure to build multi-level topology
/// returns jump residues
utility::vector1< core::Size > setup_pose_rbsegs_keep_loops(
              core::pose::Pose &pose,
              utility::vector1< protocols::rbsegment_relax::RBSegment > const &rbsegs,
              loops::Loops const &loops,
              core::kinematics::MoveMapOP mm) {
 	using namespace core::kinematics;

	//fpd pose should be rooted on vrt!
	runtime_assert( pose.residue_type( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt );

	FoldTree const &f_in = core::conformation::symmetry::get_asymm_unit_fold_tree( pose.conformation() );
	// nres points to last protein residue;
	Size totres = f_in.nres();
	Size nres = totres - 1;
 	core::Size vrtid = nres+1;
 	core::Size nrbsegs( rbsegs.size() );

 	// "star" topology fold tree
	utility::vector1< core::Size > cuts = f_in.cutpoints(), jump_res;
	utility::vector1< std::pair<core::Size,core::Size> > jumps;
	for( loops::Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it ) {
		//if (pose.residue(it->stop()).is_upper_terminus() && it->stop() != nres)
		//	cuts.push_back( it->stop() );

		//fpd segment with no rbsegs
		if (pose.residue(it->start()).is_lower_terminus() && pose.residue(it->stop()).is_upper_terminus()) {
			core::Size midpt = (it->start()+it->stop()) / 2;
			jumps.push_back (std::pair<core::Size,core::Size>( vrtid, midpt ) );
		}

		if (pose.residue(it->start()).is_lower_terminus() || pose.residue(it->stop()).is_upper_terminus()) continue;
		cuts.push_back( it->cut() );
	}
	//cuts.push_back( nres );

	for (int i=1; i <= (int)nrbsegs; ++i ) {
		// jump from vrt to midpt of 1st seg
		core::Size midpt = (rbsegs[i][1].start()+rbsegs[i][1].end()) / 2;
		jumps.push_back (std::pair<core::Size,core::Size>( vrtid, midpt ) );
		jump_res.push_back ( midpt );
	}
	for (int i=1; i <= (int)nrbsegs; ++i ) {
		for (int j=2; j<= (int)rbsegs[i].nContinuousSegments() ; ++j ) {
			// jump from midpt of 1st seg here
			core::Size midpt1 = (rbsegs[i][1].start()+rbsegs[i][1].end()) / 2;
			core::Size midpt2 = (rbsegs[i][j].start()+rbsegs[i][j].end()) / 2;
			jumps.push_back (std::pair<core::Size,core::Size>( midpt1, midpt2 ) );
			jump_res.push_back ( midpt2 );
		}
	}

	ObjexxFCL::FArray2D_int fjumps( 2, jumps.size() );
	ObjexxFCL::FArray1D_int fcuts ( cuts.size() );
	for ( Size i=1; i<=jumps.size(); ++i ) {
		fjumps(1,i) = std::min( jumps[i].first , jumps[i].second );
		fjumps(2,i) = std::max( jumps[i].first , jumps[i].second );
		// DEBUG -- PRINT JUMPS AND CUTS
		TR.Error << " jump " << i << " : " << fjumps(1,i) << " , " << fjumps(2,i) << std::endl;
	}
	for ( Size i = 1; i<=cuts.size(); ++i ) {
		fcuts(i) = cuts[i];
		TR.Error << " cut " << i << " : " << fcuts(i) << std::endl;
	}

	kinematics::FoldTree f;
	bool valid_tree = f.tree_from_jumps_and_cuts( nres+1, jumps.size(), fjumps, fcuts );
	runtime_assert( valid_tree );
	f.reorder( vrtid );


	TR << "New (asu) fold tree: " << f << std::endl;
	core::pose::symmetry::set_asymm_unit_fold_tree( pose , f );
	TR << "New fold tree: " << pose.fold_tree() << std::endl;

	// movemap
	mm->clear();
	for( loops::Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it )
		for ( int j=(int)(it->start()) ; j<=(int)(it->stop()); ++j )
			mm->set_bb(j,true);
	for (int i=1; i <= (int)nrbsegs; ++i )
		mm->set_jump(i,true);

	core::pose::PDBInfoOP pdbinfo_old;
	if (pose.pdb_info())
		pdbinfo_old = core::pose::PDBInfoOP( new core::pose::PDBInfo( *(pose.pdb_info()) ) );

	// cb variants
	for( protocols::loops::Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it ) {
	 	Size const loop_cut(it->cut());
	 	if ( !pose.residue(loop_cut).is_lower_terminus() && !pose.residue(loop_cut).is_upper_terminus() ) {
	 		if ( ! pose.residue(loop_cut).has_variant_type(core::chemical::CUTPOINT_LOWER) )
	 			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, loop_cut );
	 		if ( ! pose.residue(loop_cut+1).has_variant_type(core::chemical::CUTPOINT_UPPER) )
	 			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, loop_cut+1 );
	 	}
	}

	// copy back B factors
	if (pdbinfo_old) {
		for (Size resid = 1; resid <= nres; ++resid) {
			core::conformation::Residue const &rsd_i = pose.residue(resid);
			Size natoms_old = pdbinfo_old->natoms( resid ), natoms_new = rsd_i.natoms();
			for (Size atmid = 1; atmid <= natoms_new; ++atmid) {
				Real B = pdbinfo_old->temperature( resid, std::min(atmid,natoms_old) );
				pose.pdb_info()->temperature( resid, atmid, B );
			}
		}
		pose.pdb_info()->obsolete( false );
	}

	return jump_res;
}

////
//// set up star topology foldtree
void
setup_star_topology( core::pose::Pose & pose ) {
	utility::vector1< RBSegment > rigid_segs, rb_chunks;
	utility::vector1< core::Size > jumps;
	protocols::loops::Loops loops;
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );

	core::pose::Pose pose_asu;
	core::pose::symmetry::extract_asymmetric_unit(pose, pose_asu);

	guess_rbsegs_from_pose( pose_asu, rigid_segs, rb_chunks, loops );  // do this on asu pose
	jumps = setup_pose_rbsegs_keep_loops( pose,  rigid_segs , loops,  movemap );
}

////
//// set up disconnected
void
setup_disconnected( core::pose::Pose & pose ) {
 	using namespace core::kinematics;

	//fpd pose should be rooted on vrt!
	runtime_assert( pose.residue_type( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt );

	FoldTree const &f_in = core::conformation::symmetry::get_asymm_unit_fold_tree( pose.conformation() );
	// nres points to last protein residue;
	Size totres = f_in.nres();
	Size nres = totres - 1;
 	core::Size vrtid = nres+1;

 	// "star" topology fold tree
	utility::vector1< core::Size > cuts;
	utility::vector1< std::pair<core::Size,core::Size> > jumps;
	for( core::Size i=1; i<=nres; ++i) {
		// fpd no terminal CB variants; connect termini w/ peptide bonds
	 	if ( pose.residue(i).is_lower_terminus() ) {
			// nothing
		} else if ( i<nres && pose.residue(i+1).is_upper_terminus() ) {
			// jump only
			jumps.push_back (std::pair<core::Size,core::Size>( vrtid, i ) );
	 	} else if ( pose.residue(i).is_upper_terminus() ) {
			// cut only
			cuts.push_back( i );
		} else {
			cuts.push_back( i );
			jumps.push_back (std::pair<core::Size,core::Size>( vrtid, i ) );
		}
	}

	ObjexxFCL::FArray2D_int fjumps( 2, jumps.size() );
	ObjexxFCL::FArray1D_int fcuts ( cuts.size() );
	for ( Size i=1; i<=jumps.size(); ++i ) {
		fjumps(1,i) = std::min( jumps[i].first , jumps[i].second );
		fjumps(2,i) = std::max( jumps[i].first , jumps[i].second );
	}
	for ( Size i = 1; i<=cuts.size(); ++i ) {
		fcuts(i) = cuts[i];
	}

	kinematics::FoldTree f;
	bool valid_tree = f.tree_from_jumps_and_cuts( nres+1, jumps.size(), fjumps, fcuts );
	runtime_assert( valid_tree );
	f.reorder( vrtid );

	TR << "New (asu) fold tree: " << f << std::endl;
	core::pose::symmetry::set_asymm_unit_fold_tree( pose , f );

	core::pose::PDBInfoOP pdbinfo_old;
	if (pose.pdb_info())
		pdbinfo_old = core::pose::PDBInfoOP( new core::pose::PDBInfo( *(pose.pdb_info()) ) );

	// cb variants
	for( core::Size i=1; i<=nres; ++i) {
	 	if ( !pose.residue(i).is_lower_terminus() && !pose.residue(i).is_upper_terminus() ) {
	 		if ( ! pose.residue(i).has_variant_type(core::chemical::CUTPOINT_LOWER) )
	 			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, i );
	 		if ( ! pose.residue(i).has_variant_type(core::chemical::CUTPOINT_UPPER) )
	 			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, i );
	 	}
	}

	// copy back B factors
	if (pdbinfo_old) {
		for (Size resid = 1; resid <= nres; ++resid) {
			core::conformation::Residue const &rsd_i = pose.residue(resid);
			Size natoms_old = pdbinfo_old->natoms( resid ), natoms_new = rsd_i.natoms();
			for (Size atmid = 1; atmid <= natoms_new; ++atmid) {
				Real B = pdbinfo_old->temperature( resid, std::min(atmid,natoms_old) );
				pose.pdb_info()->temperature( resid, atmid, B );
			}
		}
		pose.pdb_info()->obsolete( false );
	}

}


/// @brief Helper function to set up a loop-removed pose
void setup_pose_from_rbsegs(
             utility::vector1< protocols::rbsegment_relax::RBSegment > const &rbsegs ,
             core::pose::Pose const &pose_in ,
             core::pose::Pose &pose_out ,
             core::id::SequenceMapping &resmap,
             core::kinematics::MoveMap &mm,
             bool fix_ligand  ) {
	using namespace core::kinematics;

	core::Size nres( pose_in.total_residue() );
	core::Size nres_rb = 0;

	// each ligand is its own rb-segment (for the purposes of fold-tree generation)
	utility::vector1< protocols::rbsegment_relax::RBSegment > rbsegs_with_ligands = rbsegs;
	core::Size last_peptide_res = nres;
	while ( !pose_in.residue( last_peptide_res ).is_protein() )
		last_peptide_res--;
	for (core::Size i=last_peptide_res+1; i<=nres; ++i) {
		if ( pose_in.residue( i ).aa() != core::chemical::aa_vrt ) {
			rbsegs_with_ligands.push_back( protocols::rbsegment_relax::RBSegment( i, i, 'X' ) );
			TR << "setup_pose_from_rbsegs: Ligand at " << i << std::endl;
		}
	}

	// count rb reses
	for ( core::Size i=1; i <= rbsegs_with_ligands.size(); ++i )
		for (core::Size j=1; j<=rbsegs_with_ligands[i].nContinuousSegments(); ++j)
			nres_rb += rbsegs_with_ligands[i][j].end() - rbsegs_with_ligands[i][j].start() + 1;

	resmap.resize( nres, nres_rb );

	pose_out.clear();
	int rb_ctr = 0;
	for ( core::Size i=1; i <= rbsegs_with_ligands.size(); ++i ) {

		for (core::Size j=1; j<=rbsegs_with_ligands[i].nContinuousSegments(); ++j) {
			core::Size rb_start  = rbsegs_with_ligands[i][j].start() ,
			           rb_end    = rbsegs_with_ligands[i][j].end()   ;
			int nsegment_reses = rb_end - rb_start + 1;
			rb_ctr++;

			char secstruct = rbsegs_with_ligands[i][j].char_type();

			if (nsegment_reses > 0) {
				core::conformation::Residue new_res_0 = pose_in.residue( rb_start );
				pose_out.append_residue_by_jump( new_res_0 , pose_out.total_residue(), "" , "" , true );
				resmap[rb_start] = pose_out.total_residue();

				for (int k=1; k<nsegment_reses; ++k) {
					core::conformation::Residue new_res_k = pose_in.residue( rb_start+k );
					pose_out.append_residue_by_bond( new_res_k );
					resmap[rb_start + k] = pose_out.total_residue();

					// set secstruct
					if ( (secstruct == 'H' || secstruct == 'E') && k!=nsegment_reses-1 )
						pose_out.set_secstruct( pose_out.total_residue(), secstruct );
				}
			}
		}
	}

	// virtual res as root
	addVirtualResAsRoot( pose_out );

	// "star" topology fold tree
	core::kinematics::FoldTree newF;
	core::Size nstart = 1, njump = 1;
	for ( core::Size i=1; i <= rbsegs_with_ligands.size(); ++i ) {
		core::Size in_start  = rbsegs_with_ligands[i][1].start() ,
				       in_end    = rbsegs_with_ligands[i][1].end();
		core::Size out_start  = nstart,
				       out_end    = nstart + (in_end - in_start);
		core::Size out_midpt  = (out_end + out_start)/2;

		newF.add_edge( nres_rb+1, out_midpt, njump );
		newF.add_edge( out_midpt, out_start,    Edge::PEPTIDE );
		newF.add_edge( out_midpt, out_end  , Edge::PEPTIDE );
		nstart = out_end+1;

		// let all jumps except those to ligands move
		if ( !fix_ligand || i <= rbsegs.size() )
			mm.set_jump( njump , true );

		njump++;

		// in strands let bb move
		if (rbsegs_with_ligands[i].isSheet()) {
			for ( core::Size j=out_start; j<=out_end; ++j) {
				mm.set_bb( j, true );
			}
		}

		core::Size in_midpt = out_midpt;
		for (core::Size j=2; j<=rbsegs_with_ligands[i].nContinuousSegments(); ++j) {
			in_start  = rbsegs_with_ligands[i][j].start();
			in_end    = rbsegs_with_ligands[i][j].end();
			out_start  = nstart;
			out_end    = nstart + (in_end - in_start);
			out_midpt  = (out_end + out_start)/2;

			newF.add_edge( in_midpt, out_midpt, njump );
			newF.add_edge( out_midpt, out_start, Edge::PEPTIDE );
			newF.add_edge( out_midpt, out_end  , Edge::PEPTIDE );
			nstart = out_end+1;
			mm.set_jump( njump , false );
			njump++;
		}
	}

	newF.reorder( nres_rb+1 );
	pose_out.fold_tree( newF );

	TR << "New fold tree: " << newF << std::endl;
}


/// @brief Helper function to restore a fully-connected pose
void restore_pose_from_rbsegs(
             utility::vector1< protocols::rbsegment_relax::RBSegment > const &rbsegs ,
             core::pose::Pose const &pose_in ,
             core::pose::Pose &pose_out /* input/output */ )
{
	using namespace core::kinematics;

	// each ligand is its own rb-segment (for the purposes of fold-tree generation)
	utility::vector1< protocols::rbsegment_relax::RBSegment > rbsegs_with_ligands = rbsegs;

	int res_rb_counter = 1;
	for ( core::Size i=1; i <= rbsegs_with_ligands.size(); ++i ) {

		for (core::Size j=1; j<=rbsegs_with_ligands[i].nContinuousSegments(); ++j) {

			core::Size rb_start  = rbsegs_with_ligands[i][j].start() ,
								 rb_end    = rbsegs_with_ligands[i][j].end()   ;
			int nsegment_reses = rb_end - rb_start + 1;

			// xyz copy??
			if (nsegment_reses > 0) {
					pose_out.copy_segment( nsegment_reses, pose_in,  rb_start, res_rb_counter  );
			}
			// copy secstruct
			for (int k=(int)rb_start; k<=(int)rb_end; ++k)
				pose_out.set_secstruct( k, pose_in.secstruct( res_rb_counter+k-rb_start ) );
			res_rb_counter += nsegment_reses;
		}

	}

	// virtual res as root
	addVirtualResAsRoot( pose_out );

	// make the virt atom the root
	core::kinematics::FoldTree newF(pose_out.fold_tree());
	newF.reorder( pose_out.total_residue() );
	pose_out.fold_tree( newF );

	TR << "New fold tree: " << newF << std::endl;
}

////
//// set up foldtree, variants, movemap, etc.
void guess_rbsegs_from_pose(
	core::pose::Pose const & pose,
	utility::vector1< RBSegment > & rigid_segs,
	utility::vector1< RBSegment > & rb_chunks,
	protocols::loops::Loops & loops
) {
	core::Size nres = pose.total_residue()-1; // terminal VRT

	rigid_segs.clear();
	rb_chunks.clear();
	loops.clear();

	// dssp parse
	core::scoring::dssp::Dssp secstruct( pose );
	ObjexxFCL::FArray1D< char > dssp_pose( nres );
	secstruct.dssp_reduced (dssp_pose);
	//secstruct.insert_ss_into_pose( pose );

	// find all helices > 5 residues
	//          strands > 3 residues
	utility::vector1< RBSegment > simple_segments;
	bool in_helix=false, in_strand=false;
	int ss_start = -1;
	for (int i=1; i<=(int)nres; ++i) {
		bool is_lower_term = pose.residue(i).is_lower_terminus();

		// strand end
		if ((is_lower_term || dssp_pose(i) != 'E') && in_strand) {
			in_strand = false;
			if (i-ss_start >= 3)
				simple_segments.push_back( RBSegment( ss_start, i-1, 'E' ) );
		}
		// helix end
		if ((is_lower_term || dssp_pose(i) != 'H') && in_helix) {
			in_helix = false;
			if (i-ss_start >= 5)
				simple_segments.push_back( RBSegment( ss_start, i-1, 'H' ) );
		}
		// strand start
		if (dssp_pose(i) == 'E' && !in_strand)  {
			in_strand = true;
			ss_start = i;
		}
		// helix start
		if (dssp_pose(i) == 'H' && !in_helix)  {
			in_helix = true;
			ss_start = i;
		}
	}

	// put at least 2 "loop residues" inbetween each RB segment as long as there is no intervening cutpoint
	// always eat into the helix (even if this leaves a helix < 3 reses)
	for (int i=1; i< (int)simple_segments.size(); ++i) {
		if (simple_segments[i+1][1].start() - simple_segments[i][1].end() <= 2) {
			core::chemical::ResidueType const &restype_i = pose.residue_type( simple_segments[i][1].end()+1 );
			if (restype_i.is_upper_terminus() || restype_i.is_lower_terminus())
				continue;

			if (simple_segments[i][1].char_type() == 'H') {
				simple_segments[i][1].set_end( simple_segments[i+1][1].start()-3 );
			} else if (simple_segments[i+1][1].char_type() == 'H') {
				simple_segments[i+1][1].set_start( simple_segments[i][1].end()+3 );
			} else {
				// eat into longer strand
				if (simple_segments[i+1][1].length() > simple_segments[i][1].length() )
					simple_segments[i+1][1].set_start( simple_segments[i][1].end()+3 );
				else
					simple_segments[i][1].set_end( simple_segments[i+1][1].start()-3 );
			}
		}
	}

	// check for "1-residue" near termini; extend RB segments if necessary
	for (int i=1; i<=(int)simple_segments.size(); ++i) {
		core::Size start_i = simple_segments[i][1].start(), stop_i = simple_segments[i][1].end();
		if (start_i>1 && pose.residue(start_i-1).is_lower_terminus()) {
			TR << start_i << " --> " << start_i-1 << std::endl;
			simple_segments[i][1].set_start(start_i-1);
		}
		if (stop_i<nres  && (pose.residue(stop_i+1).is_upper_terminus() || (stop_i+1 == nres) ) ) {
			TR << stop_i << " --> " << stop_i+1 << std::endl;
			simple_segments[i][1].set_end(stop_i+1);
		}
	}

	// auto-gen loops
	if ( simple_segments.size() > 1 ) {
		std::sort( simple_segments.begin(), simple_segments.end(), RB_lt());
		int start_res=1, end_res=simple_segments[1][1].start()-1;
		int nsegs = simple_segments.size();

		if (end_res > start_res)
			loops.push_back( protocols::loops::Loop(start_res, end_res, 0, 0.0, false) );

		for (int i=1; i<nsegs; ++i) {
			start_res = simple_segments[i][1].end()+1;
			end_res   = simple_segments[i+1][1].start()-1;
			if (end_res > start_res)
				loops.push_back( protocols::loops::Loop(start_res, end_res, 0, 0.0, false) );
		}
		start_res = simple_segments[nsegs][1].end()+1;
		end_res   = nres;
		if (end_res > start_res)
			loops.push_back( protocols::loops::Loop(start_res, end_res, 0, 0.0, false) );
	}

	// split loops on cutpoints from original pose
	utility::vector1< int > cuts_in = pose.fold_tree().cutpoints();
	std::sort( cuts_in.begin(), cuts_in.end() );
	for (Size i=1; i<=loops.size(); ++i) {
		for (Size j=1; j<=cuts_in.size(); ++j) {
			if (loops[i].start() <= Size(cuts_in[j]) && loops[i].stop() > Size(cuts_in[j])) {
				TR << "splitting [" << loops[i].start() << " , " << loops[i].stop() << "] at " << cuts_in[j] << std::endl;
				core::Size new_start = cuts_in[j]+1, new_stop = loops[i].stop();
				loops[i].set_stop(cuts_in[j]);
				if (new_stop > new_start)
					loops.push_back( protocols::loops::Loop(new_start, new_stop, 0, 0.0, false) );
				else // 1 residue segment in original foldtree ...
					simple_segments.push_back( RBSegment( new_start, new_stop, 'X' ) );
			}
		}
	}

	// now combine paired strands into a compound segment
	//   look for NH--O distance < 2.6A (?)
	utility::vector1< utility::vector1< int > > compound;
	utility::vector1< core::Size > parent_seg( simple_segments.size(), 0 );
	for (int i=1; i< (int)simple_segments.size(); ++i) {
		if (simple_segments[i][1].char_type() != 'E') continue;

		utility::vector1< int > this_seg( 1, i );
		for (int j=i+1; j<=(int)simple_segments.size(); ++j) {
			if (simple_segments[j][1].char_type() != 'E') continue;

			// foreach res in i,j
			bool found = false;
			for (int ii=(int)simple_segments[i][1].start(); ii<=(int)simple_segments[i][1].end() && !found; ++ii)
			for (int jj=(int)simple_segments[j][1].start(); jj<=(int)simple_segments[j][1].end() && !found; ++jj) {
				core::Real d2=10;

				if (pose.residue(ii).aa() != core::chemical::aa_pro && pose.residue(jj).aa() != core::chemical::aa_pro)
					d2 = std::min(
						(pose.residue(ii).atom("H").xyz() - pose.residue(jj).atom("O").xyz()).length_squared() ,
						(pose.residue(ii).atom("O").xyz() - pose.residue(jj).atom("H").xyz()).length_squared() );
				else if (pose.residue(jj).aa() != core::chemical::aa_pro)
					d2 = (pose.residue(ii).atom("O").xyz() - pose.residue(jj).atom("H").xyz()).length_squared();
				else if (pose.residue(ii).aa() != core::chemical::aa_pro)
					d2 = (pose.residue(jj).atom("O").xyz() - pose.residue(ii).atom("H").xyz()).length_squared();

				if (d2 < 2.6*2.6) {
					this_seg.push_back(j);
					if (parent_seg[i] == 0)
						parent_seg[i] = i;
					if (parent_seg[j] == 0)
						parent_seg[j] = parent_seg[i];
					else {
						// tricky case ... j is already mapped
						// in this case map everything mapped to i to parent_seg[j]
						for (int k=1; k<j; ++k)
							if ((int)parent_seg[k] == i) parent_seg[k] = parent_seg[j];
					}

					TR << "Merging " << j << " (" << simple_segments[j][1].start() << "," << simple_segments[j][1].end() << ") ";
					TR << " to " << i << " (" <<  simple_segments[i][1].start() << "," << simple_segments[i][1].end() << ") " << std::endl;
					found = true;
				}
			}
		}
	}

	// make the compound segments
	for (int i=1; i< (int)simple_segments.size(); ++i) {
		if (((int)parent_seg[i]) != i) continue; // not compound or already added
		utility::vector1< RBSegment > thisLockSeg;
		for (core::Size j=i; j<=simple_segments.size(); ++j) {
			if ( ((int)parent_seg[j]) == i )
				thisLockSeg.push_back( simple_segments[ j ] );
		}
		rigid_segs.push_back( RBSegment( thisLockSeg ) );
	}
	// add in all other simple segs
	for (int i=1; i<=(int)simple_segments.size(); ++i)
		if (parent_seg[i] == 0)
			rigid_segs.push_back( simple_segments[i] );

	// sort loops & rbsegs, choose cutpoints
	std::sort( rigid_segs.begin(), rigid_segs.end(), RB_lt());
	std::sort( loops.v_begin(), loops.v_end(), protocols::loops::Loop_lt());
	loops.auto_choose_cutpoints( pose );

	// define chunks
	rb_chunks = rigid_segs;
	for (int i=1; i<=(int)rigid_segs.size(); ++i) {
		for (int j=1; j<=(int)rigid_segs[i].nContinuousSegments() ; ++j ) {
			core::Size c_low=1, c_high=nres;
			for (int k=1; k<=(int)loops.size(); ++k) {
				if (loops[k].cut() < rigid_segs[i][j].start() && loops[k].cut() > c_low )
					c_low = loops[k].cut();
				if (loops[k].cut()+1 > rigid_segs[i][j].end() && loops[k].cut()+1 < c_high )
					c_high = loops[k].cut()+1;
			}
			rb_chunks[i][j].set_start(c_low);
			rb_chunks[i][j].set_end(c_high);
		}
	}
}

//// remap all segs
void remap_rb_segments(
                 utility::vector1< RBSegment > const &rbsegs,
                 utility::vector1< RBSegment > &rbsegs_remap,
                 core::id::SequenceMapping const &resmap ) {
	for ( RBConsIt it_seg = rbsegs.begin(), it_seg_end = rbsegs.end(); it_seg != it_seg_end; ++it_seg ) {
		RBSegment it_remap = it_seg->remap( resmap );
		rbsegs_remap.push_back( it_remap );
	}
}


}
}
