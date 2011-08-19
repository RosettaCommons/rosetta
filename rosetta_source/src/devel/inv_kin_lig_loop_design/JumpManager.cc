// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/JumpManager.cc
///
/// @brief
/// @author

#include <cmath>
#include <devel/inv_kin_lig_loop_design/JumpManager.hh>
#include <devel/inv_kin_lig_loop_design/ResID.hh>
#include <devel/inv_kin_lig_loop_design/std_extra.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

#include <numeric/random/random.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/types.hh>
#include <core/coarse/Translator.fwd.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/io/pdb/file_data.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomPointer.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/AtomWithDOFChange.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/DomainMap.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
//#include <core/optimization/MinimizerMap.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <basic/MetricValue.fwd.hh>
// AUTO-REMOVED #include <basic/OStream.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <devel/inv_kin_lig_loop_design/Loop.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/all.fwd.hh>
#include <utility/io/icstream.fwd.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/ocstream.fwd.hh>
#include <utility/io/orstream.fwd.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/BodyPosition.fwd.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/Quaternion.fwd.hh>
#include <numeric/all.fwd.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzTriple.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <numeric/random/random.fwd.hh>
#include <numeric/random/uniform.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>


namespace devel {

  namespace inv_kin_lig_loop_design {

    using namespace std;

    bool JumpManager::is_allowable_jump_atom( core::pose::Pose& pose, int seqpos, int atom ) {

      if( atom == 0 ) return false;

      cout << "cout is allowable " << seqpos << " " << atom << endl;
      if( pose.residue(seqpos).is_protein() ) {


				if( pose.residue(seqpos).atom_is_backbone(atom) ) {
					// protein backbone atom => all heavy atoms can be jump atoms
					return !pose.residue(seqpos).atom_is_hydrogen(atom);
				}
				else {
					// protein sidechain atom => only atoms defining chi angles can be jump atoms
					for( Size k = 1; k <= pose.residue(seqpos).nchi(); ++k ) {
						core::chemical::AtomIndices const& atom_indices = pose.residue(seqpos).chi_atoms(k);
						for( Size l = 1; l <= atom_indices.size(); ++l ) {
							if( atom == static_cast<int>(atom_indices[l]) ) {
								return true;
							}
						}
					}
					return false;
				}
      }
      else if( pose.residue(seqpos).is_ligand() ) {
				return true;
      }
      else {
				cout << "WARNING - JumpManager::is_allowable_jump_atom - no rules for non-protein, non-ligand allowable jump atoms" << endl;
				assert( false );
				return false;
      }

    }

    string JumpManager::get_jump_atom_for_hbond( core::pose::Pose& pose, int seqpos, string hbond_atom_name ) {

      int hbond_atom = pose.residue(seqpos).atom_index(hbond_atom_name);

      int abase1 = pose.residue(seqpos).atom_base(hbond_atom);
      int abase2 = pose.residue(seqpos).atom_base(abase1);

      int rval = -1;

      if( is_allowable_jump_atom( pose, seqpos, hbond_atom ) ) {
				cout << "get_jump_atom_for_hbond = " << seqpos << " " << hbond_atom << "=>" << hbond_atom << endl;
				rval = hbond_atom;
      }
      else if( is_allowable_jump_atom(pose,seqpos,abase1 ) ) {
				cout << "get_jump_atom_for_hbond = " << seqpos << " " << hbond_atom << "=>" << abase1 << endl;
				assert( abase1 != 0 );
				rval = abase1;
      }
      else if( is_allowable_jump_atom(pose,seqpos,abase2 ) ) {
				cout << "get_jump_atom_for_hbond = " << seqpos << " " << hbond_atom << "=>" << abase2 << endl;
				assert( abase2 != 0 );
				rval = abase2;
      }
      else {
				cout << "WARNING none of these are allowable jump atoms: " << seqpos << "." << hbond_atom << "," << abase1 << "," << abase2 << endl;
				assert( false );
      }

      assert( rval != -1 );
      return pose.residue(seqpos).atom_name(rval);

    }

		//     std::pair<int,int> JumpManager::get_jump_atoms_for_hbond(core::pose::Pose& pose,
		//							     int seqpos_from, int hbond_atom_from,
		//							     int seqpos_to,   int hbond_atom_to ) {

		//       int jump_atom_from = get_jump_atom_for_hbond( pose, seqpos_from, hbond_atom_from );
		//       int jump_atom_to   = get_jump_atom_for_hbond( pose, seqpos_to,   hbond_atom_to   );



		//       return std::pair<int,int>( jump_atom_from, jump_atom_to );

		//     }


    namespace {

      core::kinematics::RT operator*(core::kinematics::RT const& a, core::kinematics::RT const& b ) {
				core::kinematics::RT rval;
				rval.set_rotation( a.get_rotation() * b.get_rotation() );
				rval.set_translation( a.get_translation() + a.get_rotation() * b.get_translation() );
				return rval;
      }

    }


    namespace {

      core::kinematics::Stub get_stub( core::pose::Pose& pose, int seqpos, int atom ) {

				assert( atom != 0 );

				core::conformation::Residue const& r = pose.residue(seqpos);

				int abase1 = r.atom_base(atom);
				int abase2 = r.atom_base(abase1);

				typedef core::kinematics::Stub::Vector vector_type;

				if( abase1 == 0 ) {
					cout << "WARNING abase1 == 0, resetting to atom-1" << endl;
					abase1 = atom-1;
				}

				if( abase2 == 0 ) {
					cout << "WARNING abase2 == 0, resetting to abase1-1" << endl;
					abase2 = abase1-1;
				}

				cout << "getting stub for " << seqpos << " " << atom << " " << abase1 << " " << abase2 << " " << r.atom_name(atom) << " " << r.atom_name(abase1) << " " << r.atom_name(abase2) << endl;

				assert( abase1 != 0 );
				assert( abase2 != 0 );

				return core::kinematics::Stub( r.atom(atom).xyz(),
																			 r.atom(atom).xyz(),
																			 r.atom( abase1 ).xyz(),
																			 r.atom( abase2 ).xyz() );

      }

      //core::kinematics::RT get_hbond_rt( double dist, double xBAH, double xAHD, double xBAHD, bool h2a ) {

      namespace {

				int random_sign() {
					return 2*numeric::random::random_range(0,1) - 1;
				}

      }

	}

      core::kinematics::RT JumpManager::get_hbond_rt( TagPtr tag, bool h2a ) {

				double const dist = tag->getOption<double>("dist",1.8)  + numeric::random::gaussian() * tag->getOption<double>("sigma_dist",0.1);
				double const xBAH = (tag->getOption<double>("xBAH",120) + numeric::random::gaussian() * tag->getOption<double>("sigma_xBAH",5) ) * random_sign();
				double const xAHD = (tag->getOption<double>("xAHD",180) + numeric::random::gaussian() * tag->getOption<double>("sigma_xAHD",5) );

				//double const xAHD = (tag->getOption<double>("xBAHD",180) + numeric::random::gaussian() * tag->getOption<double>("sigma_xAHD",5) );
				double const xBAHD = numeric::random::uniform() * 360;
				double const xBA   = numeric::random::uniform() * 360;
				double const xHD   = numeric::random::uniform() * 360;

				const double DEG2RAD = M_PI/180.0;

				// gives the transformation between hbond acceptor/donor

				static const numeric::xyzVector<double> x(1,0,0);
				static const numeric::xyzVector<double> y(0,1,0);
				static const numeric::xyzVector<double> z(0,0,1);

				core::kinematics::RT qBAH, qAHD, q_flip, qHA, qBAHD, qBA, qHD;

				qBAH.set_rotation( numeric::rotation_matrix<>( y, DEG2RAD*xBAH));
				q_flip.set_rotation( numeric::rotation_matrix<>( y, DEG2RAD*180.0));
				qHA.set_translation( numeric::xyzVector<double>(dist,0,0) );
				qAHD.set_rotation( numeric::rotation_matrix<>( y, DEG2RAD*xAHD ) );
				qBAHD.set_rotation( numeric::rotation_matrix<>(x, DEG2RAD*xBAHD ) );
				qBA.set_rotation( numeric::rotation_matrix<>(x, DEG2RAD*xBA ) );
				qHD.set_rotation( numeric::rotation_matrix<>(x, DEG2RAD*xHD ) );

				if( h2a ) {
					return qHD * qAHD * qBAHD * q_flip * qHA * qBAH * qBA;
				}
				else {
					return qBA * qBAH * qBAHD * q_flip * qHA * qAHD * qHD;
				}

      }


    void JumpManager::set_random_hbond_jump(core::pose::Pose& pose, Loop const& loop ) {

      //     core::kinematics::Jump JumpManager::get_random_hbond_jump(core::pose::Pose& pose, int jumpno,
      //							      int seqpos_from, int hbond_atom_from, int jump_atom_from,
      //							      int seqpos_to,   int hbond_atom_to ,  int jump_atom_to ) {
      //cout << "get random hbond jump, from=" << pose.residue(seqpos_from).atom_name(hbond_atom_from) << "/" << pose.residue(seqpos_from).atom_name(jump_atom_from) << "=>" << pose.residue(seqpos_to).atom_name(hbond_atom_to) << "/" << pose.residue(seqpos_to).atom_name(jump_atom_to) << endl;

      core::kinematics::Edge const & edge = pose.fold_tree().jump_edge( loop.jumpno );
      int from = edge.start();
      int to = edge.stop();
      int hbond_from_atom = pose.residue(from).atom_index( loop.tag->getTag("from")->getOption<std::string>("atom") );
      int hbond_to_atom   = pose.residue(to).atom_index( loop.tag->getTag("to")->getOption<std::string>("atom") );

      bool h2a = pose.residue(from).atom_is_hydrogen(hbond_from_atom);
      cout << "h2a=" << h2a << endl;

      core::kinematics::Stub
				sf = pose.conformation().upstream_jump_stub( loop.jumpno ),
				sa = get_stub(pose,from,hbond_from_atom),
				sh = get_stub(pose,to,hbond_to_atom),
				st = pose.conformation().downstream_jump_stub( loop.jumpno );

      assert( loop.tag->hasTag("hbond") );

      core::kinematics::RT qfa(sf,sa);
      //qfa.from_stubs(sf,sa);
      core::kinematics::RT qah = get_hbond_rt( loop.tag->getTag("hbond"), h2a );
      core::kinematics::RT qht(sh,st);
      //qht.from_stubs(sh,st);

      core::kinematics::Jump jump( qfa * qah * qht );

      pose.set_jump(loop.jumpno, jump );

    }

		/*
    core::kinematics::Jump JumpManager::get_jump_from_template( TagPtr tag_segment ) {

			assert( tag_segment->getName() == "anchored_loop" );
			assert( tag_segment->hasTag("template") );

			core::pose::PoseOP pose( new core::pose::Pose );

			string const filename = tag_segment->getTag("template")->getOption<string>("filename");
			ResID const res_id_from = tag_segment->getTag("template")->getOption<ResID>("from");
			ResID const res_id_to   = tag_segment->getTag("template")->getOption<ResID>("to");

			cout << "reading filename " << filename << endl;

			core::io::pdb::pose_from_pdb( *pose, filename );

			cout << "getting resid map" << endl;

			devel::inv_kin_lig_loop_design::resids_type resids = ResID::get_resids( *pose );

			// need to set the same foldtree as we have in this

			Residue const* r_from = find_or_throw(resids,res_id_from);
			Residue const* r_to   = find_or_throw(resids,res_id_to);

			string atom_name_from = tag_segment->getTag("from")->getOption<string>("atom");
			string atom_name_to   = tag_segment->getTag("to")->getOption<string>("atom");

			cout << "adding jump between " << r_from->seqpos() << "." << atom_name_from << " and " << r_to->seqpos() << "." << atom_name_to << endl;


			// build a fold tree with the same jump as we have in the primary pose

			int from = r_from->seqpos();
			int to = r_to->seqpos();
			int jumpno = 1;

			core::kinematics::FoldTree ft;
			core::kinematics::Edge edge( from, to, jumpno );
			edge.start_atom() = r_from->atom_index( atom_name_from );
			edge.stop_atom() = r_to->atom_index( atom_name_to );
			ft.add_edge(edge);
			ft.add_edge( to, 1, core::kinematics::Edge::PEPTIDE );
			ft.add_edge( to, from-1, core::kinematics::Edge::PEPTIDE ); // NB: this is a big huge assumption
			ft.reorder( from );

			cout << "get_jump: setting new fold tree" << endl;

			pose->fold_tree(ft);
			pose->energies();

			core::id::AtomID id_from( r_from->atom_index(atom_name_from), r_from->seqpos() );
			core::id::AtomID id_to  ( r_to  ->atom_index(atom_name_to), r_to  ->seqpos() );

			//        core::kinematics::Atom const& atom_from = pose->atom_tree().atom( id_from );
			//        core::kinematics::Atom const& atom_to   = pose->atom_tree().atom( id_to   );

			cout << "to:" << endl;
			cout << r_to->seqpos() << endl;
			//print_stubs(pose,atom_to);

			core::kinematics::Jump rval( pose->jump(jumpno) );

			return rval;
    }

    void JumpManager::set_template_jump( core::pose::Pose& pose, Loop const& loop ) {
      core::kinematics::Jump jump = get_jump_from_template( loop.tag );
      pose.set_jump( loop.jumpno, jump );
    }

*/
    core::kinematics::RT JumpManager::get_template_rt( TagPtr tag_segment ) {

			assert( tag_segment->getName() == "anchored_loop" );
			assert( tag_segment->hasTag("template") );

			TagPtr tag_template = tag_segment->getTag("template");

			core::pose::PoseOP pose( new core::pose::Pose );

			string const filename = tag_template->getOption<string>("filename");
			ResID const res_id_from = tag_template->getOption<ResID>("from");
			ResID const res_id_to   = tag_template->getOption<ResID>("to");

			cout << "reading filename " << filename << endl;

			core::import_pose::pose_from_pdb( *pose, filename );

			cout << "getting residues" << endl;
			devel::inv_kin_lig_loop_design::resids_type resids = ResID::get_resids( *pose );

			Residue const* r_from = find_or_throw(resids,res_id_from);
			Residue const* r_to   = find_or_throw(resids,res_id_to);

			string const tag_start_atom = tag_template->hasOption("from_atom") ? tag_template->getOption<string>("from_atom") : tag_segment->getTag("from")->getOption<string>("atom");
			string const tag_stop_atom  = tag_template->hasOption("to_atom")   ? tag_template->getOption<string>("to_atom")   : tag_segment->getTag("to")->getOption<string>("atom");

			int from = r_from->seqpos();
			int to = r_to->seqpos();

			int const from_atom = pose->residue(from).atom_index(tag_start_atom);
			int const to_atom   = pose->residue(to).atom_index(tag_stop_atom);

      core::kinematics::Stub
				//sf = pose->conformation().upstream_jump_stub( jumpno ),
				sa = get_stub(*pose,from,from_atom),
				sh = get_stub(*pose,to,to_atom);
				//st = pose->conformation().downstream_jump_stub( jumpno );

			core::kinematics::RT qah(sa,sh);

			return qah;

    }

		void JumpManager::set_template_jump( core::pose::Pose& pose, Loop const& loop ) {

			assert( loop.tag->hasTag("template") );

      //     core::kinematics::Jump JumpManager::get_random_hbond_jump(core::pose::Pose& pose, int jumpno,
      //							      int seqpos_from, int hbond_atom_from, int jump_atom_from,
      //							      int seqpos_to,   int hbond_atom_to ,  int jump_atom_to ) {
      //cout << "get random hbond jump, from=" << pose.residue(seqpos_from).atom_name(hbond_atom_from) << "/" << pose.residue(seqpos_from).atom_name(jump_atom_from) << "=>" << pose.residue(seqpos_to).atom_name(hbond_atom_to) << "/" << pose.residue(seqpos_to).atom_name(jump_atom_to) << endl;

      core::kinematics::Edge const & edge = pose.fold_tree().jump_edge( loop.jumpno );
      int from = edge.start();
      int to = edge.stop();
      int from_atom = pose.residue(from).atom_index( loop.tag->getTag("from")->getOption<std::string>("atom") );
      int to_atom   = pose.residue(to).atom_index( loop.tag->getTag("to")->getOption<std::string>("atom") );

      //bool h2a = pose.residue(from).atom_is_hydrogen(hbond_from_atom);
      //cout << "h2a=" << h2a << endl;

      core::kinematics::Stub
				sf = pose.conformation().upstream_jump_stub( loop.jumpno ),
				sa = get_stub(pose,from,from_atom),
				sh = get_stub(pose,to,to_atom),
				st = pose.conformation().downstream_jump_stub( loop.jumpno );

      core::kinematics::RT qfa(sf,sa);
      //qfa.from_stubs(sf,sa);
      core::kinematics::RT qah = get_template_rt( loop.tag );
      core::kinematics::RT qht(sh,st);
      //qht.from_stubs(sh,st);

      //core::kinematics::Jump jump( qfa * qah * qht );
			core::kinematics::Jump jump( qfa * qah * qht );

      pose.set_jump(loop.jumpno, jump );

    }

  } // namespace LoopDesign

} // namespace Devel
