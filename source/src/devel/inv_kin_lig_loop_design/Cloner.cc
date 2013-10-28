// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/Cloner.cc
///
/// @brief
/// @author
#include <devel/inv_kin_lig_loop_design/Cloner.hh>
#include <devel/inv_kin_lig_loop_design/Ints.hh>
#include <devel/inv_kin_lig_loop_design/std_extra.hh>
#include <devel/inv_kin_lig_loop_design/JumpManager.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>


#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>

// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/util.hh>

#include <core/chemical/VariantType.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>



#define FORVC(Iter,Type,Vec)  for( utility::vector0<Type>::const_iterator Iter  = (Vec).begin(); Iter != (Vec).end(); ++Iter)
#define FORTAGS(Iter,Name,ATagCOP) const utility::vector0<TagCOP>& TagCOP ## Name ## Iter ## __LINE__ = (ATagCOP)->getTags(#Name); for( utility::vector0<TagCOP>::const_iterator Iter = (TagCOP ## Name ## Iter ## __LINE__ ).begin(); Iter != (TagCOP ## Name ## Iter ## __LINE__ ).end(); ++Iter )

using namespace std;

//using core::conformation::Residue;

namespace devel {

  namespace inv_kin_lig_loop_design {

    namespace {

      // =================================================================
      // ==================== collapse_random_options ====================
      // =================================================================

      void set_nres(TagCOP const& tag, string which ) {
				using devel::inv_kin_lig_loop_design::Ints;

				if( tag->hasOption(which) ) {
					const Ints ints = tag->getOption<Ints>(which);
					const int  nres = ints.getRandomInt();
					cout << "InvKinLig:LoopDesign::set_nres - setting " << tag->getName() << "-" << tag->getOption<string>("name","<unnamed>") << " " << ints << " => " << nres << endl;
					const_cast< utility::tag::Tag * >(tag())->setOption(which,nres);
				}
      }  // set_nres

      TagCOP collapse_random_options( TagCOP tag ) {

				TagCOP const rval = tag->clone();

				FORTAGS(i,loop,rval) {
					set_nres(*i,"nres_pre");
					set_nres(*i,"nres_post");
				} // i

				FORTAGS(i,anchored_loop,rval) {
					set_nres(*i,"nres_pre");
					set_nres(*i,"nres_post");
				} // i

				FORTAGS(i,deletion,rval) {
				} // i

				utility::vector0<TagCOP> vVary;
				/// apl -- replacing std::vector specific += operator from std_extra.hh with insert() calls.
				utility::vector0<TagCOP> looptags = rval->getTags("loop");
				vVary.insert(vVary.end(),looptags.begin(),looptags.end());

				utility::vector0<TagCOP> anchortags = rval->getTags("anchored_loop");
				vVary.insert(vVary.end(),anchortags.begin(),anchortags.end());

				utility::vector0<TagCOP> deletiontags = rval->getTags("deletion");
				vVary.insert(vVary.end(),deletiontags.begin(),deletiontags.end());

				FORVC(i,TagCOP,vVary) {
					set_nres(*i,"vary_pre");
					set_nres(*i,"vary_post");
				} // i

				return rval;

      }

      // ====================================================
      // ==================== get_indels ====================
      // ====================================================

      segments_type get_indels(TagCOP const& tag0, resids_type const& resids ) {

				//cout << "get_indels: " << endl;

				segments_type rval;

				FORTAGS(i,loop,tag0) {
					TagCOP tag = *i;
					//      <loop nres_pre=3 nres_post=3 weight_cb=1 ss=EEEELL>
					//           <begin res_id=A:178>
					//           <end res_id=A:184>
					//      </loop>
					//cout << *tag << endl;

					Segment indel;
					indel.tag = tag;
					indel.lo_res = find_or_throw( resids, tag->getTag("begin")->getOption<ResID>("res_id") );
					indel.hi_res = find_or_throw( resids, tag->getTag("end")->getOption<ResID>("res_id") );
					indel.type  = Segment::LOOP;
					indel.nres_pre  = tag->getOption<int>("nres_pre");
					indel.nres_post = tag->getOption<int>("nres_post");
					indel.aas += vector<core::chemical::AA>( indel.nres_pre + indel.nres_post, core::chemical::aa_ala );

					rval.push_back(indel);

				} // i

				FORTAGS(i,anchored_loop,tag0) {
					TagCOP tag = *i;
					//cout << *tag << endl;
					//   <anchored_loop nres_pre=2 nres_post=5 ss=EEHLLLLH weight_cb_begin=2 weight_cb_end=1>
					//           <begin res_id=A:150>
					//           <end   res_id=A:157>
					//           <from  res_id=L:1   type=ligand>
					//           <to    res_id=A:153 type=amide res_type="GLN">
					//           <template pdbfile=ctl.pdb mdlfile=ctl_lig.mdl mdl_chain_id=L mdl_res_type="CTL" from=L:1 to=A:153>
					//      </anchored_loop>

					Segment indel;
					indel.tag = tag;
					indel.lo_res = find_or_throw(resids, tag->getTag("begin")->getOption<ResID>("res_id" ) );
					indel.hi_res = find_or_throw(resids, tag->getTag("end")->getOption<ResID>("res_id" ) );
					indel.from_res = find_or_throw( resids, tag->getTag("from")->getOption<ResID>("res_id" ) );
					indel.type  = Segment::ANCHORED_LOOP;

					indel.nres_pre  = tag->getOption<int>("nres_pre");
					indel.nres_post = tag->getOption<int>("nres_post");

					indel.aas += vector<core::chemical::AA>( indel.nres_pre, core::chemical::aa_ala );
					indel.aas.push_back( core::chemical::aa_from_name( tag->getTag("to")->getOption<string>("res_type") ) );
					indel.aas += vector<core::chemical::AA>( indel.nres_post, core::chemical::aa_ala );

					rval.push_back(indel);

				} // i

				FORTAGS(i,deletion,tag0) {
					TagCOP tag = *i;
					//cout << *tag << endl;

					Segment indel;
					indel.lo_res = find_or_throw(resids, tag->getTag("begin")->getOption<ResID>("res_id") );
					indel.type = Segment::DELETION;
					indel.nres_pre  = 0;
					indel.nres_post = 0;

					if( tag->hasTag("end") ) {
						indel.hi_res = find_or_throw( resids, tag->getTag("end")->getOption<ResID>("res_id") );
					}
					else {
						indel.hi_res = indel.lo_res; // !!! 1 residue deletion
					}

					rval.push_back(indel);

				} // i

				//cout << "done with get_indels" << endl;

				return rval;

      } // get_indels

    } // namespace

    // ================================================
    // ==================== Cloner ====================
    // ================================================

    Cloner::Cloner( TagCOP const& tag0, core::pose::PoseOP pose0) : tag0( collapse_random_options(tag0) ), pose0(pose0) {
      cout << "after random options collapsed: " << endl << tag0 << endl;
    }

    core::kinematics::FoldTree Cloner::getFoldTree() {
			//     core::kinematics::FoldTree get_fold_tree(TagCOP tag0,
			//					     core::pose::PoseOP pose0,
			//					     core::pose::PoseOP pose1,
			//					     resids_type& resids,
			//					     clones_type& clones,
			//					     segments_type& segments ) {


      core::kinematics::FoldTree ft;

      //       ft.defaultRoot(tube);
      //       ft.defaultJumps(tube);

      int n_jump = 1;

      for( Size k = 0; k < segments.size(); ++k ) {

				Residue* lo_res = segments[k].lo_res; // apl -- this has got to go
				Residue* hi_res = segments[k].hi_res; // apl -- this has got to go

				if( segments[k].type == Segment::ORIGINAL ) {
					ft.add_edge( find_or_throw(clones,lo_res)->seqpos(),
											 find_or_throw(clones,hi_res)->seqpos(),
											 core::kinematics::Edge::PEPTIDE );
				}
				else if( segments[k].type == Segment::LOOP ) {

					int lo = find_or_throw(clones,lo_res)->seqpos();
					int hi = find_or_throw(clones,hi_res)->seqpos();

					cout << "InvKinLigLoopDesign::Cloner::adding loop segment " << segments[k].nres_pre << " " << segments[k].nres_post << " " << lo << " " << hi << endl;

					ft.add_edge( lo - 1, lo - 1 + segments[k].nres_pre,  core::kinematics::Edge::PEPTIDE );
					ft.add_edge( hi + 1, hi + 1 - segments[k].nres_post, core::kinematics::Edge::PEPTIDE );
					ft.add_edge( lo-1, hi+1, n_jump++ );

					//	  ft.add_edge( lo_res->seqpos()-1, lo_res->seqpos() + segments[k].nres_pre, core::kinetic::Edge::PEPTIDE );
					//	  ft.add_edge( hi_res->seqpos() - segments[k].nres_post, hi_res->seqpos()+1, core::kinetic::Edge::PEPTIDE );

					//assert( lo_res->seqpos() + segments[k].nres_pre + 1 == hi_res->seqpos() - segments[k].nres_post );

					cout << "InvKinLigLoopDesign::Cloner::break is between: " << lo - 1 + segments[k].nres_pre << " " << hi + 1 - segments[k].nres_post << endl;
					cout << "InvKinLigLoopDesign::Cloner::jump is between " << lo - 1 << " " << hi + 1 << endl;

				}
				else if ( segments[k].type == Segment::ANCHORED_LOOP ) {
					//assert( false );

					// will want the atom info at this point....

					int lo = find_or_throw(clones,lo_res)->seqpos();
					int hi = find_or_throw(clones,hi_res)->seqpos();

					int from = find_or_throw( clones, segments[k].from_res )->seqpos();
					int to   = lo + segments[k].nres_pre;

					segments[k].jumpno = n_jump;

					//ft.add_edge( from, to, n_jump++ );

					cout << "InvKinLigLoopDesign::Cloner::creating anchored loop jump edge" << endl;

					core::kinematics::Edge edge(from,to,n_jump++); // NB: we want this to behave as a jump but also have atom info, ie, not a chemical bond
					//	  edge.start_atom() = get_edge_start_atom( pose1, from, segments[k].tag );
					//	  edge.stop_atom() = get_edge_stop_atom( pose1, to, segments[k].tag );

					string const& tag_start_atom = segments[k].tag->getTag("from")->getOption<string>("atom");
					string const& tag_stop_atom = segments[k].tag->getTag("to")->getOption<string>("atom");

					int tag_start_atom_index = pose1->residue(from).atom_index( tag_start_atom );
					int tag_stop_atom_index  = pose1->residue(to).atom_index( tag_stop_atom );

					devel::inv_kin_lig_loop_design::JumpManager jm;
					edge.start_atom() = jm.get_jump_atom_for_hbond( *pose1, from, tag_start_atom );
					edge.stop_atom()  = jm.get_jump_atom_for_hbond( *pose1, to, tag_stop_atom );

					cout << "InvKinLigLoopDesign::Cloner::jump edge: " << tag_start_atom << "=" << tag_start_atom_index << "," << edge.start_atom() << " " << tag_stop_atom << "=" << tag_stop_atom_index << "," << edge.stop_atom() << endl;
					cout << "InvKinLigLoopDesign::Cloner::edge.has_atom_info() " << edge.has_atom_info() << " start=" << edge.start_atom() << " stop=" << edge.stop_atom() << endl;

					ft.add_edge( edge );

					ft.add_edge( to, lo, core::kinematics::Edge::PEPTIDE );
					ft.add_edge( to, hi, core::kinematics::Edge::PEPTIDE );
					ft.add_edge( lo-1, hi+1, n_jump++ );

					//	  cout << "adding loop segment " << segments[k].nres_pre << " " << segments[k].nres_post << " " << lo << " " << hi << endl;
					//	  ft.add_edge( lo - 1, lo - 1 + segments[k].nres_pre,  core::kinematics::Edge::PEPTIDE );
					//	  ft.add_edge( hi + 1, hi + 1 - segments[k].nres_post, core::kinematics::Edge::PEPTIDE );
					//	  ft.add_edge( lo-1, hi+1, n_jump++ );
					//	  //	  ft.add_edge( lo_res->seqpos()-1, lo_res->seqpos() + segments[k].nres_pre, core::kinetic::Edge::PEPTIDE );
					//	  //	  ft.add_edge( hi_res->seqpos() - segments[k].nres_post, hi_res->seqpos()+1, core::kinetic::Edge::PEPTIDE );
					//	  //assert( lo_res->seqpos() + segments[k].nres_pre + 1 == hi_res->seqpos() - segments[k].nres_post );

					cout << "InvKinLigLoopDesign::Cloner::jump between " << from << " " << to << endl;
					cout << "InvKinLigLoopDesign::Cloner::edge between " << lo << " " << from << endl;
					cout << "InvKinLigLoopDesign::Cloner::edge between " << from << " " << hi << endl;
					cout << "InvKinLigLoopDesign::Cloner::jump between " << lo-1 << " " << hi + 1 << endl;

				}
				else if ( segments[k].type == Segment::LIGAND ) {
					assert( segments[k].lo_res == segments[k].hi_res );
					ft.add_edge( 1, find_or_throw(clones,lo_res)->seqpos(), n_jump++ );
				}

      }


      cout << "InvKinLigLoopDesign::Cloner::before reorder: " << ft << endl;
      cout << "InvKinLigLoopDesign::Cloner::reordering fold tree... " << endl;
      ft.reorder(1);
      cout << "InvKinLigLoopDesign::Cloner::after reorder: " << ft << endl;

      return ft;

    } // get_fold_tree

		//     core::kinematics::FoldTree Cloner::getFoldTree() {
		//       //return core::kinematics::FoldTree();
		//       return get_fold_tree(tag0,pose0,pose1,resids,clones,segments);
		//     } // Cloner::getFoldTree

    void append_seq( core::pose::Pose& pose0, int sp_lo, int sp_hi, utility::vector1< core::chemical::ResidueTypeCOP >& seq ) {

      //seq.clear();

      string s;

      for( int i = sp_lo; i <= sp_hi; ++i ) {
				seq.push_back( core::chemical::ResidueTypeCOP( pose0.residue_type(i) ) );
				s.push_back( pose0.residue(i).name1() );
      }

      cout << "InvKinLigLoopDesign: creating a pose of sequence: " << s << endl;

    }

    utility::vector1< core::chemical::ResidueTypeCOP > get_seq_from_aas( vector< core::chemical::AA > const& aas ) {

      core::chemical::ResidueTypeSetCAP residue_set
				( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

      utility::vector1< core::chemical::ResidueTypeCOP > rval;
      for( Size i = 0; i < aas.size(); ++i ) {
				core::chemical::ResidueTypeCOPs const& res_types = residue_set->aa_map( aas[i] );
				assert( res_types.size() != 0 );
				rval.push_back( res_types[1] );
      }
      return rval;
    }

    namespace {

      void make_pose_from_sequence(core::pose::Pose & pose, utility::vector1< core::chemical::ResidueTypeCOP > const& sequence ) {

				using namespace core::chemical;
				// clear all of the old data in the pose
				pose.clear();

				// setup the pose by appending the appropriate residues residues
				for ( Size seqpos = 1; seqpos <= sequence.size(); ++seqpos ) {

					bool const is_lower_terminus( seqpos == 1 );
					//bool const is_upper_terminus( seqpos == sequence.size() );
					//bool const is_terminus( is_lower_terminus || is_upper_terminus ); // redundant, but for convenience


					core::chemical::ResidueType const & rsd_type = *(sequence[seqpos]);
					core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( rsd_type ) );

					if ( is_lower_terminus || !new_rsd->is_polymer() ) {
						pose.append_residue_by_jump( *new_rsd, 1 /*pose.total_residue()*/ );
					} else {
						pose.append_residue_by_bond( *new_rsd, true );
					}
				} // for seqpos
				// pose.conformation().insert_chain_ending( pose.total_residue() - 1 );
      } // make_pose_from_sequence

    }

    Segment get_segment_from_indel( Segment const& indel ) {
      Segment rval;

      assert( indel.aas.size() != 0 );

      rval.pose = core::pose::PoseOP( new core::pose::Pose );

      make_pose_from_sequence( *rval.pose, get_seq_from_aas(indel.aas) );

      rval.lo_res = const_cast< core::conformation::Residue * > ( & rval.pose->residue( 1 ) ); /// apl -- note -- you should never do this; I'm only doing it to keep legacy code compiling.
      rval.hi_res = const_cast< core::conformation::Residue * > ( & rval.pose->residue( rval.pose->total_residue() ) ); /// apl -- note -- you should never do this; I'm only doing it to keep legacy code compiling.
      //for( core::conformation::ResidueOPs::iterator iter = rval.pose->res_begin(); iter != rval.pose->res_end(); ++iter ) {
		//		rval.hi_res = iter->get(); // rbegin()
      //}

      rval.type = indel.type;
      rval.nres_pre  = indel.nres_pre;
      rval.nres_post = indel.nres_post;
      rval.from_res = indel.from_res;
      rval.tag = indel.tag;

      return rval;
    }

    // ==============================================================
    // ==================== get_pose_with_indels ====================
    // ==============================================================

    core::pose::PoseOP get_pose_with_indels( core::pose::PoseOP pose0, segments_type const& indels, segments_type& segments, clones_type& clones ) {

      //cout << "get_pose_with_indels" << endl;

      segments.clear();

      //core::conformation::ResidueOPs::iterator iter = pose0->res_begin();

      Residue* r_prev = 0; // apl -- this has got to go
      Residue* r = const_cast< core::conformation::Residue * > ( & pose0->residue( 1 ) );// apl -- this has got to go

      Segment segment;

      segment.pose  = pose0;
      segment.lo_res = r;
      segment.type  = Segment::ORIGINAL;

		Size ii( 1 );
		while( ii <= pose0->total_residue() ) {
      ///while( iter != pose0->res_end() ) {

				for( Size k = 0; k < indels.size(); ++k ) {

					if( r == indels[k].lo_res ) {

						segment.hi_res = r_prev;
						segments.push_back( segment );

						do {
							//cout << "skipping " << r->seqpos() << endl;
							//++iter;
							//assert( iter != pose0->res_end() );
							++ii;
							assert( ii <= pose0->total_residue() );

							r_prev = r;
							r = const_cast< core::conformation::Residue * > ( & pose0->residue( ii ) );//iter->get();

						} while( r != indels[k].hi_res );

						segment.lo_res = r;


						// now need to instantiate a new pose for the segment
						if( indels[k].aas.size() != 0 ) {
							//cout << "creating new pose from indel" << endl;
							Segment indel = get_segment_from_indel( indels[k] );
							segments.push_back( indel );
						}

					}

				}


				if( !r->is_polymer() ) {
					//cout << "found a ligand" << endl;

					if( r_prev && r_prev->is_polymer() ) {
						segment.hi_res = r_prev;
						segments.push_back( segment );
					}

					Segment ligand;
					ligand.lo_res  = r;
					ligand.hi_res = r;
					ligand.type = Segment::LIGAND;
					ligand.pose = pose0;

					segments.push_back(ligand);

				}
				else {
					//cout << "including " << r->seqpos() << endl;
				}

				//++iter;
				++ii;
				r_prev = r;
				//if( iter != pose0->res_end() ) {
				if ( ii <= pose0->total_residue() ) {
					r = const_cast< core::conformation::Residue * > ( & pose0->residue( ii ) );//iter->get();
				}

				if( r_prev && !r_prev->is_polymer() ) {
					segment.lo_res = r;
				}

      }

      if( r_prev && r_prev->is_polymer() ) {
				segment.hi_res = r_prev;
				segments.push_back( segment );
      }


      // Step 2 - determine the total size of the new pose

      int total_size = 0;

      utility::vector1< core::chemical::ResidueTypeCOP > seq;

      for( Size k = 0; k < segments.size(); ++k ) {
				int size = segments[k].hi_res->seqpos() - segments[k].lo_res->seqpos() + 1;
				//cout << segments[k].lo_res->seqpos() << ":" << segments[k].hi_res->seqpos() << " " << size << endl;
				total_size += size;

				append_seq(*segments[k].pose,segments[k].lo_res->seqpos(),segments[k].hi_res->seqpos(),seq);
      }

      //cout << "creating a pose" << endl;

      core::pose::PoseOP rval( new core::pose::Pose );

			//       cout << "getting a residue set" << endl;
			//       core::chemical::ResidueTypeSetCAP residue_set
			//	( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );


      //string seq(total_size,'A');
      //cout << "creating a pose from sequence: " << endl; // << seq << endl;
      make_pose_from_sequence( *rval, seq );

      //resize( *rval, total_size );


      clones.clear();

      Size k, begin,size;
      for( k = 0, begin = 1; k < segments.size(); ++k, begin += size ) {
				size = segments[k].size();
				//cout << "copying segment " << k << endl;
				rval->copy_segment(size, *segments[k].pose, begin, segments[k].lo_res->seqpos() );
      }




      for( k = 0, begin = 1; k < segments.size(); ++k, begin += size ) {
				size = segments[k].size();

				if( segments[k].type == Segment::ANCHORED_LOOP ) {
					int lo = begin;
					int hi = begin + size - 1;
					assert( lo != 1 && hi != static_cast<int>(rval->n_residue()) );

					//cout << "set_cutpoints: Adding cutpoints between " << lo-1 << "," << lo << " and " << hi << "," << hi+1 << endl;

					core::pose::add_variant_type_to_pose_residue( *rval, core::chemical::CUTPOINT_LOWER, lo - 1 );
					core::pose::add_variant_type_to_pose_residue( *rval, core::chemical::CUTPOINT_UPPER, lo );
					core::pose::add_variant_type_to_pose_residue( *rval, core::chemical::CUTPOINT_LOWER, hi );
					core::pose::add_variant_type_to_pose_residue( *rval, core::chemical::CUTPOINT_UPPER, hi + 1 );

				}
				else if( segments[k].type == Segment::LOOP ) {
					int lo = begin;
					//int hi = begin + size - 1;
					int mid = lo + segments[k].nres_pre;
					//cout << "adding variant: " << lo << "," << mid << "," << mid+1 << endl;

					core::pose::add_variant_type_to_pose_residue( *rval, core::chemical::CUTPOINT_LOWER, mid - 1);
					core::pose::add_variant_type_to_pose_residue( *rval, core::chemical::CUTPOINT_UPPER, mid );
				}

      }

      for( k = 0, begin = 1; k < segments.size(); ++k, begin += size ) {

				size = segments[k].size();

				for( int i = 0; i < static_cast<int>(size); ++i ) {

					Residue* r0 = const_cast<Residue*>(&(segments[k].pose->residue( segments[k].lo_res->seqpos() + i ) ) ); // this has got to go
					Residue* r1 = const_cast<Residue*>(& (rval->residue(begin + i)) ); // this has got to go
					//cout << "clones: " << r0 << "\t" << r1 << endl;
					clones[ r0 ] = r1;
				}

      }

      //cout << "done with cloning, clones.size() = " << clones.size() << endl;

      return rval;

    }

		//     Residue const* get_cloned_residue( Cloner::ResIDMap_type const& residss,
		//				       Cloner::CloneMap_type const& clones,
		//				       ResID const& res_id ) {

		//       return find_or_default( clones, find_or_throw( residss, res_id ), static_cast<Residue const*>(0) );
		//     }



    core::pose::PoseOP Cloner::clone() {

      resids = ResID::get_resids( *pose0 );

      //vector< Indel > vIndels = get_indels(tag0,resids);
      indels = get_indels(tag0,resids);

      cout << "InvKinLigLoopDesign::Cloner:: number of indels = " << indels.size() << endl;

			//       for( Size k = 0; k < indels.size(); ++k ) {
			//	cout << "hi " << k << endl;
			//	cout << indels[k].begin << " " << indels[k].end << endl;
			//	cout << "indel: " << k << ": " << flush << ResID(*indels[k].begin) << "-" << ResID(*indels[k].end)<< endl;
			//       }

      //cout << "calling get_pose_with_indels" << endl;

      pose1 = get_pose_with_indels( pose0, indels, segments, clones );

      //cout << "returning from clone" << endl;

      return pose1;
    }

		//     core::kinematics::Stub get_stub(Residue const* r, string const& a, string const& b, string const& c ) {
		//       return core::kinematics::Stub( r->atom( a ).xyz(),
		//				     r->atom( b ).xyz(),
		//				     r->atom( c ).xyz() );
		//     }

		//     core::kinematics::Stub get_stub(Residue const* r, TagCOP tag ) {

		// //       typedef numeric::xyzVector< Real > Vector;
		// //       if( tag->hasOption("moiety") ) {
		// //	string const moiety = tag->getOption<string>("moiety");
		// //       }
		// //       else if( ) {
		// //       }

		//     }

		//     void print_stubs( core::pose::PoseOP pose, core::kinematics::Atom const& atom ) {

		//       Residue const& r = pose->residue( atom.id().rsd() );

		//       cout << r.seqpos() << endl;

		//        cout << "0: " << atom.id().rsd() << endl;
		//        cout << "1: " << atom.stub_atom1_id().rsd() << endl;
		//        cout << "2: " << atom.stub_atom2_id().rsd() << endl;
		//        cout << "3: " << atom.stub_atom3_id().rsd() << endl;

		//        cout << "0: " << r.atom_name( atom.id().atomno() ) << endl;
		//        cout << "1: " << r.atom_name( atom.stub_atom1_id().atomno() ) << endl;
		//        cout << "2: " << r.atom_name( atom.stub_atom2_id().atomno() ) << endl;
		//        cout << "3: " << r.atom_name( atom.stub_atom3_id().atomno() ) << endl;

		//        cout << "i1: " << atom.input_stub_atom1_id().rsd() << endl;
		//        cout << "i2: " << atom.input_stub_atom2_id().rsd() << endl;
		//        cout << "i3: " << atom.input_stub_atom3_id().rsd() << endl;

		//        Residue const& r_in = pose->residue( atom.input_stub_atom1_id().rsd() );

		//        cout << "i1: " << r_in.atom_name( atom.input_stub_atom1_id().atomno() ) << endl;
		//        cout << "i2: " << r_in.atom_name( atom.input_stub_atom2_id().atomno() ) << endl;
		//        cout << "i3: " << r_in.atom_name( atom.input_stub_atom3_id().atomno() ) << endl;

		//     }

    void set_secstruct( core::pose::Pose& pose, int const a, int const b, const string& ss ) {

      cout << "InvKinLigLoopDesign::Cloner::setting secondary structure of " << a << "," << b << " to " << ss << endl;

      const int n = ss.size();

      const int diff = b-a+1;

      if( n == diff ) {
				for( int i = 0; i < n; ++i ) {
					pose.set_secstruct(a+i,ss[i]);
				}
      }
      else {
				// handling the case where the secondary structure string doesn't exactly match up with the length of residues

				int n_pre = 0, n_post = 0;
				if( n > diff ) {
					n_pre  = diff / 2;
					n_post = diff - n_pre;
				}
				else if( n < diff ) {
					n_pre  = n / 2;
					n_post = n - n_pre;
				}
				else {
					assert( false );
				}

				//Residue* r;
				int i;

				for( i = 0; i < n_pre; ++i ) {
					pose.set_secstruct(a+i,ss[i]);
				} // i

				for( i = 0; i < n_post; ++i ) {
					pose.set_secstruct(b-i,ss[n - 1 - i]);
				} // i

      }

    } // State::set_secstruct


    void Cloner::setInitialConfig() {

      core::pose::set_ss_from_phipsi(*pose1);

      //cout << "indels.size() = " << indels.size() << endl;
      //cout << "segments.size() = " << segments.size() << endl;

      for( Size k = 0; k < segments.size(); ++k ) {

				if( segments[k].type == Segment::LOOP ) {
					int lo = find_or_throw( clones, segments[k].lo_res  )->seqpos();
					int hi = find_or_throw( clones, segments[k].hi_res )->seqpos();

					// repair icoor for lo
					//	  repair_icoor(pose1,lo-1);
					//	  repair_icoor(pose1,lo);
					//	  repair_icoor(pose1,hi);
					//	  repair_icoor(pose1,hi+1);

					cout << "InvKinLigLoopDesign::Cloner::inserting ideal geometry at polymer bond " << lo-1 << "," << lo << endl;

					pose1->conformation().insert_ideal_geometry_at_polymer_bond( lo-1 );

					cout << "InvKinLigLoopDesign::Cloner::inserting ideal geometry at polymer bond " << hi << "," << hi+1 << endl;

					pose1->conformation().insert_ideal_geometry_at_polymer_bond( hi );

					// NB: need to be a little more careful to preserve the right angles...

				}
				else if( segments[k].type == Segment::ANCHORED_LOOP ) {

					// don't really need this
					cout << "InvKinLigLoopDesign::Cloner::setting jumpno " << segments[k].jumpno << " to initial value" << endl;

					JumpManager jm;
					if( segments[k].tag->hasTag("template") ) {
						jm.set_template_jump( *pose1, Loop(segments[k], clones) );
					}
					else if( segments[k].tag->hasTag("hbond") ) {
						jm.set_random_hbond_jump( *pose1, Loop(segments[k], clones) );
					}
					else {
						assert( false );
					}

				}

				if( segments[k].type == Segment::LOOP ||
						segments[k].type == Segment::ANCHORED_LOOP ) {

					int lo = find_or_throw( clones, segments[k].lo_res  )->seqpos();
					int hi = find_or_throw( clones, segments[k].hi_res )->seqpos();

					for( int i = lo; i <= hi ; ++i ) {
						pose1->set_phi( i, -135.0 );
						pose1->set_psi( i,  135.0 );
						pose1->set_omega( i, 180.0 );
						pose1->set_secstruct(i,'L');
					}

					if( segments[k].tag.get() && segments[k].tag->hasOption("ss") ) {
						set_secstruct(*pose1,lo,hi,segments[k].tag->getOption<string>("ss") );
					}

				}

      }

      cout << "InvKinLigLoopDesign::Cloner::secstruct= " << endl;
      for( int i = 1; i <= static_cast<int>(pose1->n_residue()); ++i ) {
				//cout << i << " " << pose1->secstruct(i);
				cout << pose1->secstruct(i);
      }
      cout << endl;

    }

    // ==========================================================
    // ==================== Cloner::getLoops ====================
    // ==========================================================

    vector< Loop > Cloner::getLoops() {
      vector< Loop > rval;
      for( Size k = 0; k < segments.size(); ++k ) {
				if( segments[k].type == Segment::LOOP ||
						segments[k].type == Segment::ANCHORED_LOOP ) {
					rval.push_back( Loop( segments[k], clones ) );
				}
      }
      return rval;
    }

  } // namespace LoopDesign
} // namespace devel
