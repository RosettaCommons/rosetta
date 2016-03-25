// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief protocols for folding into density
/// @details
/// @author Frank DiMaio

#include <protocols/cryst/TwoDLattice.hh>
#include <protocols/cryst/refinable_lattice_creator.hh>
#include <protocols/cryst/util.hh>

#include <core/types.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>

#include <core/scoring/electron_density/util.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymRotamerTrialsMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/relax/FastRelax.hh>



#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <core/scoring/constraints/util.hh>

#include <utility/excn/Exceptions.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>


#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace cryst {

static THREAD_LOCAL basic::Tracer TR("protocols.cryst.refinable_lattice");

using namespace protocols;
using namespace core;
using namespace kinematics;
using namespace scoring;
using namespace scoring::symmetry;
using namespace conformation;
using namespace conformation::symmetry;


////////////////////////////////////////////////////

// creators
std::string
MakeLayerMoverCreator::keyname() const { return MakeLayerMoverCreator::mover_name(); }

protocols::moves::MoverOP
MakeLayerMoverCreator::create_mover() const { return protocols::moves::MoverOP(new MakeLayerMover); }

std::string
MakeLayerMoverCreator::mover_name() { return "MakeLayerMover"; }

////////////////////////////////////////////////////

void
MakeLayerMover::apply( core::pose::Pose & pose ) {
	// initialize wg_ from pose CRYST1 line
	io::CrystInfo ci = pose.pdb_info()->crystinfo();
	runtime_assert(ci.A()*ci.B()*ci.C() != 0);  // TODO: allow these to be randomized

	wg_.set_wallpaper_group(ci.spacegroup());

	core::Real angle = ci.beta();
	if ( wg_.setting() == wgHEXAGONAL ) angle=ci.gamma();

	wg_.set_parameters(ci.A(),ci.B(), angle);

	Size rootres = place_near_origin( pose );

	core::pose::Pose posebase;
	utility::vector1<Size> Ajumps, Bjumps, monomer_jumps, monomer_anchors;
	Size base_monomer;

	build_layer_of_virtuals( posebase, Ajumps, Bjumps, monomer_anchors, base_monomer);

	Size nvrt = posebase.total_residue();
	Size nres_monomer = pose.total_residue();

	detect_connecting_subunits( pose, posebase, monomer_anchors, base_monomer );
	add_monomers_to_layer( pose, posebase, monomer_anchors, monomer_jumps, rootres );

	Size nsubunits = monomer_anchors.size();

	core::pose::PDBInfoOP pdbinfo_old=pose.pdb_info(), pdbinfo_new;
	pose = posebase;

	conformation::symmetry::SymmetryInfo syminfo;
	setup_xtal_symminfo( pose, nsubunits, nvrt, base_monomer, nres_monomer, Ajumps, Bjumps, monomer_jumps, syminfo );

	//fpd hack for disulfides
	bool symmdetectdiulf = basic::options::option[ basic::options::OptionKeys::symmetry::detect_bonds ]();
	basic::options::option[ basic::options::OptionKeys::symmetry::detect_bonds ].value(false);
	pose::symmetry::make_symmetric_pose( pose, syminfo );
	basic::options::option[ basic::options::OptionKeys::symmetry::detect_bonds ].value(symmdetectdiulf);

	pdbinfo_new = pose::PDBInfoOP( new pose::PDBInfo( pose, true ) );
	core::pose::symmetry::make_symmetric_pdb_info( pose, pdbinfo_old, pdbinfo_new );
	pdbinfo_new->set_crystinfo(ci);
	pose.pdb_info( pdbinfo_new );

	//fpd foldtree above assumes no monomer jumps
	//fpd here, we reinitialize the foldtree to let the SymmetryInfo machinery take care of
	//    monomer jumps
	//fpd this also has the advantage of renumbering jumps in a consistent way so that
	//    all symmetric foldtree manipulations work well with symm poses made from this function
	core::kinematics::FoldTree ft = core::conformation::symmetry::get_asymm_unit_fold_tree(pose.conformation());
	core::conformation::symmetry::symmetrize_fold_tree(pose.conformation(), ft);
	pose.fold_tree( ft );

	// force symmetrization
	core::conformation::symmetry::SymmetryInfoCOP symminfo_new =
		dynamic_cast<core::conformation::symmetry::SymmetricConformation const & >( pose.conformation()).Symmetry_Info();
	for ( core::Size i=pose.fold_tree().num_jump(); i>=1; --i ) {
		if ( symminfo_new->jump_is_independent(i) && symminfo_new->jump_clones(i).size() > 0 ) {
			core::kinematics::Jump j_i = pose.jump( i );
			pose.set_jump( i, j_i );
		}
	}

	// update disulf info
	pose.conformation().detect_disulfides();
}


Size
MakeLayerMover::place_near_origin (
	Pose & pose
) {
	Size rootpos=0;
	Size nres = pose.total_residue();

	Vector com(0,0,0);
	for ( Size i=1; i<= nres; ++i ) {
		com += pose.residue(i).xyz("CA");
		if ( pose.residue(i).is_upper_terminus() ) break;
	}
	com /= nres;

	Real mindis2(1e6);
	for ( Size i=1; i<= nres; ++i ) {
		Real const dis2( com.distance_squared(  pose.residue(i).xyz("CA") ) );
		if ( dis2 < mindis2 ) {
			mindis2 = dis2;
			rootpos = i;
		}
		if ( pose.residue(i).is_upper_terminus() ) break;
	}

	Size nsymm = wg_.nsymmops();
	Size bestxform=0;
	Vector bestoffset(0,0,0);
	mindis2=1e6;
	com = wg_.c2f()*com;
	for ( Size i=1; i<=nsymm; ++i ) {
		Vector foffset = wg_.symmop(i).get_rotation()*com + wg_.symmop(i).get_translation(), rfoffset;
		rfoffset[0] = pos_mod( foffset[0], 1.0 );
		rfoffset[1] = pos_mod( foffset[1], 1.0 );
		rfoffset[2] = 0.0;
		Real dist = (wg_.f2c()*rfoffset).length_squared();
		if ( dist<mindis2 ) {
			mindis2=dist;
			bestxform=i;
			bestoffset = foffset - rfoffset;
		}
	}

	numeric::xyzMatrix<Real> R = wg_.f2c()*wg_.symmop(bestxform).get_rotation()*wg_.c2f();
	numeric::xyzVector<Real> T = wg_.f2c()*(wg_.symmop(bestxform).get_translation() - bestoffset);
	pose.apply_transform_Rx_plus_v( R,T );

	return rootpos;
}

void
MakeLayerMover::detect_connecting_subunits(
	Pose const & monomer_pose,
	Pose const & pose,
	utility::vector1<Size> & monomer_anchors,
	Size &basesubunit
) {
	utility::vector1<Size> new_monomer_anchors;
	Size new_basesubunit=0;

	// get pose radius
	Size const num_monomers( monomer_anchors.size() ), nres_monomer( monomer_pose.total_residue () );

	Vector com(0,0,0);
	Real radius = 0;
	utility::vector1<Vector> monomer_cas(nres_monomer);

	Vector T0 = pose.residue(monomer_anchors[basesubunit]).xyz("ORIG");
	runtime_assert( T0.length() < 1e-6);

	for ( Size i=1; i<= nres_monomer; ++i ) {
		Vector ca_i = monomer_pose.residue(i).xyz("CA");
		monomer_cas[i] = ca_i;
		radius = std::max( (ca_i).length_squared() , radius );
	}
	radius = sqrt(radius);

	// make master first
	new_monomer_anchors.push_back(monomer_anchors[basesubunit]);
	new_basesubunit = 1;

	for ( Size i=1; i<=num_monomers; ++i ) {
		if ( i==basesubunit ) continue;

		// pass 1 check vrt-vrt dist to throw out very distant things
		Vector T = pose.residue(monomer_anchors[i]).xyz("ORIG");
		Real disVRT = T.length();
		if ( disVRT>contact_dist_+2*radius ) continue;

		// pass 2 check ca-ca dists
		Vector X = pose.residue(monomer_anchors[i]).xyz("X") - T;
		Vector Y = pose.residue(monomer_anchors[i]).xyz("Y") - T;
		Vector Z = X.cross( Y );

		numeric::xyzMatrix<core::Real> R = numeric::xyzMatrix<core::Real>::cols(X,Y,Z);

		bool contact=false;
		for ( Size j=1; j<= nres_monomer && !contact; ++j ) {
			Vector Ri = R*monomer_cas[j] + T;
			for ( Size k=1; k<= nres_monomer && !contact; ++k ) {
				contact = ((Ri-monomer_cas[k]).length_squared() < contact_dist_*contact_dist_);
			}
		}
		if ( contact ) {
			new_monomer_anchors.push_back(monomer_anchors[i]);
		}
	}
	basesubunit = new_basesubunit;
	monomer_anchors = new_monomer_anchors;
}


void
MakeLayerMover::add_monomers_to_layer(
	Pose const & monomer_pose,
	Pose & pose,
	utility::vector1<Size> const & monomer_anchors,
	utility::vector1<Size> & monomer_jumps,
	Size rootpos
) {
	Size const num_monomers( monomer_anchors.size() ), nres_monomer( monomer_pose.total_residue () );

	monomer_jumps.clear();
	Size n_framework_jumps = pose.fold_tree().num_jump();

	Size nres_protein(0);
	for ( Size i=1; i<= num_monomers; ++i ) {
		//std::cerr << "inserting " << i << " of " << num_monomers << std::endl;
		Size const anchor( monomer_anchors[i] + nres_protein ); // since we've already done some insertions
		Size const old_nres_protein( nres_protein );
		pose.insert_residue_by_jump( monomer_pose.residue(rootpos), old_nres_protein+1, anchor ); ++nres_protein;
		for ( Size j=rootpos-1; j>=1; --j ) {
			pose.prepend_polymer_residue_before_seqpos( monomer_pose.residue(j), old_nres_protein+1, false ); ++nres_protein;
		}
		for ( Size j=rootpos+1; j<= monomer_pose.total_residue(); ++j ) {
			if ( monomer_pose.residue(j).is_lower_terminus() ) {
				pose.insert_residue_by_jump( monomer_pose.residue(j), nres_protein+1, nres_protein ); ++nres_protein;
			} else {
				pose.append_polymer_residue_after_seqpos( monomer_pose.residue(j), nres_protein, false ); ++nres_protein;
			}
		}
		monomer_jumps.push_back( n_framework_jumps+i );
	}
	for ( Size i=1; i<= num_monomers; ++i ) {
		pose.conformation().insert_chain_ending( nres_monomer*i );
	}

	core::kinematics::FoldTree f( pose.fold_tree() );
	f.reorder( nres_protein+1 );
	pose.fold_tree(f);
}


void
MakeLayerMover::build_layer_of_virtuals(
	core::pose::Pose & posebase,
	utility::vector1<Size> &Ajumps,
	utility::vector1<Size> &Bjumps,
	utility::vector1<Size> &subunit_anchors,
	Size &basesubunit
) {
	numeric::xyzVector<Size> EXTEND((Size)2,(Size)2,(Size)0);
	numeric::xyzVector<int> nvrts_dim(2*EXTEND[0]+1, 2*EXTEND[1]+1, 1);
	posebase.clear();

	// get gridspace in each dimension
	numeric::xyzVector< Size > grid = wg_.get_nsubdivisions();
	numeric::xyzVector< Size > trans_dofs = wg_.get_trans_dofs();

	Vector O, Ax, Bx, Cx, Ay, By, Cy;
	if ( wg_.setting() == wgHEXAGONAL ) {
		Ax = wg_.f2c()*Vector(1,0,0); Ax.normalize();
		Bx = wg_.f2c()*Vector(0,1,0); Bx.normalize();
		Cx = Vector(0,0,1); Cx.normalize();
	} else {
		Ax = Vector(1,0,0); Ax.normalize();
		Bx = Vector(0,1,0); Bx.normalize();
		Cx = Vector(0,0,1); Cx.normalize();
	}

	// we don't care where y is as long as it is perpendicular
	Ay = Bx - Ax.dot(Bx)*Ax; Ay.normalize();
	By = Cx - Bx.dot(Cx)*Bx; By.normalize();
	Cy = Ax - Cx.dot(Ax)*Cx; Cy.normalize();

	ObjexxFCL::FArray2D<int> vrtX, vrtY;
	vrtX.dimension(nvrts_dim[0]*grid[0]+1,1,1); vrtX=0;
	vrtY.dimension(nvrts_dim[0]*grid[0]+1,nvrts_dim[1]*grid[1]+1,1); vrtY=0;

	// 1 expand A (j==k==1)
	for ( int i=1; i<=(int)(nvrts_dim[0]*grid[0]+1); ++i ) {
		Vector fX( (Real)(i-1-(Real)(EXTEND[0]*grid[0]))/(Real)grid[0], -1.0*EXTEND[1], 0 );
		O = wg_.f2c()*fX;

		// now add 3 virtuals to the pose, with X pointing toward A,B,C, respectively
		ResidueOP vrt_x = make_vrt(O,Ax,Ay);
		ResidueOP vrt_y = make_vrt(O,Bx,By);

		if ( i==1 ) {
			posebase.append_residue_by_bond( *vrt_x );    vrtX(1,1) = 1;
			posebase.append_residue_by_jump( *vrt_y, 1);  vrtY(1,1) = 2;
		} else {
			posebase.append_residue_by_jump( *vrt_x, vrtX(i-1,1));  vrtX(i,1) = posebase.total_residue();
			Ajumps.push_back(posebase.fold_tree().num_jump());
			posebase.append_residue_by_jump( *vrt_y, vrtX(i  ,1));  vrtY(i,1) = posebase.total_residue();
		}
	}

	// expand B (k==1)
	for ( int i=1; i<=(int)(nvrts_dim[0]*grid[0]+1); ++i ) {
		for ( int j=2; j<=(int)(nvrts_dim[1]*grid[1]+1); ++j ) {
			Vector fX( (i-1-(Real)(EXTEND[0]*grid[0]))/(Real)grid[0], (j-1-(Real)(EXTEND[1]*grid[1]))/(Real)grid[1], 0 );
			O = wg_.f2c()*fX;

			// now add 3 virtuals to the pose, with X pointing toward A,B,C, respectively
			ResidueOP vrt_y = make_vrt(O,Bx,By);

			posebase.append_residue_by_jump( *vrt_y, vrtY(i,j-1)); vrtY(i,j) = posebase.total_residue();
			Bjumps.push_back(posebase.fold_tree().num_jump());
		}
	}

	// add "hanging" virtuals
	for ( int s=1; s<=(int)wg_.nsymmops(); ++s ) {
		numeric::xyzMatrix<Real> R_i = wg_.symmop(s).get_rotation();
		numeric::xyzVector<Real> T_i = wg_.symmop(s).get_translation();

		// T_i -> indices
		for ( int i=-(int)EXTEND[0]; i<=(int)EXTEND[0]; ++i ) {
			for ( int j=-(int)EXTEND[1]; j<=(int)EXTEND[1]; ++j ) {
				// find lattice anchor
				int x_i = (int)std::floor( (i+EXTEND[0]+T_i[0])*grid[0] + 1.5 );
				int y_i = (int)std::floor( (j+EXTEND[1]+T_i[1])*grid[1] + 1.5 );

				O = posebase.residue(vrtY(x_i,y_i)).xyz("ORIG");

				if ( wg_.setting() == wgHEXAGONAL ) {
					Ax = wg_.f2c()*R_i*wg_.c2f()*Vector(1,0,0); Ax.normalize();
					Ay = wg_.f2c()*R_i*wg_.c2f()*Vector(0,1,0); Ay.normalize();
				} else {
					Ax = R_i*Vector(1,0,0); Ax.normalize();
					Ay = R_i*Vector(0,1,0); Ay.normalize();
				}
				ResidueOP vrt_sub = make_vrt(O,Ax,Ay);

				posebase.append_residue_by_jump( *vrt_sub, vrtY(x_i,y_i));

				subunit_anchors.push_back(posebase.total_residue());
				if ( s==1 && i==0 && j==0 ) {
					basesubunit = subunit_anchors.size();
				}
			}
		}
	}
}


void
MakeLayerMover::setup_xtal_symminfo(
	Pose const & pose,
	Size const num_monomers,
	Size const num_virtuals,
	Size const base_monomer,
	Size const nres_monomer,
	utility::vector1<Size> const &Ajumps,
	utility::vector1<Size> const &Bjumps,
	utility::vector1<Size> const &monomer_jumps,
	conformation::symmetry::SymmetryInfo & symminfo
) {

	for ( Size i=1; i<= num_monomers; ++i ) {
		if ( i != base_monomer ) {
			Size const offset( (i-1)*nres_monomer ), base_offset( (base_monomer-1)*nres_monomer);
			for ( Size j=1; j<= nres_monomer; ++j ) {
				symminfo.add_bb_clone ( base_offset+j, offset+j );
				symminfo.add_chi_clone( base_offset+j, offset+j );
			}
		}
	}

	// subunit jump clones
	Size const base_monomer_jump( monomer_jumps[ base_monomer ] );
	for ( Size i=1; i<=monomer_jumps.size(); ++i ) {
		if ( monomer_jumps[i]!= base_monomer_jump ) {
			symminfo.add_jump_clone( base_monomer_jump, monomer_jumps[i], 0.0 );
		}
	}


	// unit cell clones
	using core::conformation::symmetry::SymDof;
	std::map< Size, SymDof > symdofs;

	Size Amaster=Ajumps[1], Bmaster=Bjumps[1];
	numeric::xyzVector<core::Size> linked_dofs=wg_.get_trans_dofs();
	if ( linked_dofs[1]==1 ) { Bmaster=Ajumps[1]; /*std::cerr << "clone B->A" << std::endl;*/ }

	for ( Size i=2; i<=Ajumps.size(); ++i ) {
		symminfo.add_jump_clone( Amaster, Ajumps[i], 0.0 );
	}
	for ( Size i=1; i<=Bjumps.size(); ++i ) {
		if ( Bmaster!=Bjumps[i] ) {
			symminfo.add_jump_clone( Bmaster, Bjumps[i], 0.0 );
		}
	}

	SymDof symdof_a;
	SymDof symdof_b;

	core::Size nrot_dofs=wg_.get_nrot_dofs();
	if ( nrot_dofs == 1 ) {
		symdof_b.read( "x y" );
	} else {
		symdof_b.read( "x" );
	}
	symdof_a.read( "x" );

	if ( moving_lattice_ ) {
		symdofs[ Ajumps[1] ] = symdof_a;
		if ( Bmaster==Bjumps[1] ) {
			symdofs[ Bjumps[1] ] = symdof_b;
		}
	}

	// jump names
	for ( Size v=1; v<=pose.num_jump(); ++v ) symminfo.set_jump_name(v, "v_"+utility::to_string(v));

	symminfo.set_jump_name(Ajumps[1], "A");
	symminfo.set_jump_name(Bjumps[1], "B");

	SymDof symdof_m;
	symdof_m.read( wg_.get_moveable_dofs()+" angle_x angle_y angle_z"); // even z!!
	symdofs[ base_monomer_jump ] = symdof_m;
	symminfo.set_jump_name(base_monomer_jump, "SUB");

	symminfo.set_dofs( symdofs );

	symminfo.num_virtuals( num_virtuals );
	symminfo.set_use_symmetry( true );
	symminfo.set_flat_score_multiply( pose.total_residue(), 0 );
	symminfo.set_nres_subunit( nres_monomer );

	Size const nres_protein( num_monomers * nres_monomer );
	for ( Size i=1; i<= nres_protein; ++i ) {
		if ( symminfo.bb_is_independent( i ) ) symminfo.set_score_multiply( i, 2 );
		else symminfo.set_score_multiply( i, 1 );
	}
	symminfo.update_score_multiply_factor();
}

// parse_my_tag
void
MakeLayerMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & ,
	moves::Movers_map const & ,
	core::pose::Pose const & /*pose*/ )
{
	if ( tag->hasOption( "contact_dist" ) ) {
		contact_dist_ = tag->getOption<core::Real>( "contact_dist" );
	}
}

}
}
