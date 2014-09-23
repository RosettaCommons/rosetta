// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/copy_dofs/CopyDofs.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/copydofs/CopyDofs.hh>
#include <core/pose/copydofs/util.hh>
#include <core/conformation/Residue.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>

#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "core.pose.copydofs.CopyDofs" );

using namespace core;

namespace core {
namespace pose {
namespace copydofs {

	//Constructor
	CopyDofs::CopyDofs( pose::MiniPose const & template_pose,
											std::map < id::AtomID , id::AtomID > const & atom_id_map,
											std::map< id::AtomID, Size > const & atom_id_domain_map ):
		scratch_pose_( template_pose ),
		atom_id_map_( atom_id_map ),
		atom_id_domain_map_( atom_id_domain_map ),
		atom_id_domain_map_inputted_( true )
	{}

	CopyDofs::CopyDofs( pose::MiniPose const & template_pose,
											std::map < id::AtomID , id::AtomID > const & atom_id_map ):
		scratch_pose_( template_pose ),
		atom_id_map_( atom_id_map ),
		atom_id_domain_map_inputted_( false )
	{}

	//Destructor
	CopyDofs::~CopyDofs()
	{}

	////////////////////////////////////////////////////////
	void
	CopyDofs::apply( pose::Pose & pose ){

		figure_out_atom_id_domain_map( pose );
		figure_out_dofs( pose );
		apply_dofs( pose, copy_dofs_info_ );

	}

	////////////////////////////////////////////////////////
	void
	CopyDofs::figure_out_atom_id_domain_map( pose::Pose & pose ){
		if ( atom_id_domain_map_inputted_ ) return;
		atom_id_domain_map_ = blank_atom_id_domain_map( pose );
	}

	////////////////////////////////////////////////////////
	void
	CopyDofs::figure_out_dofs( pose::Pose & pose ){
		using namespace core::id;
		using namespace core::kinematics;
		using namespace core::pose;

		copy_dofs_info_.clear();

		//Useful ID's.
		AtomID current_atom_scratch_atom_id( 1, 1);//( current_atom_scratch_atomno, current_atom_scratch_rsd );
		AtomID input_stub_atom1_scratch_atom_id( 1, 1);
		AtomID input_stub_atom2_scratch_atom_id( 1, 1);
		AtomID input_stub_atom3_scratch_atom_id( 1, 1);

		AtomID stub_atom1_scratch_atom_id( 1, 1);
		AtomID stub_atom2_scratch_atom_id( 1, 1);
		AtomID stub_atom3_scratch_atom_id( 1, 1);
		AtomID reference_scratch_atom_id( 1, 1);
		AtomID dummy_atom_id( 1, 1);

		pose::Pose const & reference_pose( pose );

		for ( std::map < id::AtomID , id::AtomID >::const_iterator
						it=atom_id_map_.begin(), it_end = atom_id_map_.end(); it != it_end; ++it ) {

			bool verbose( false );

			Size const i = (it->first).rsd(); //Residue index in big pose.
			Size const j = (it->first).atomno(); //Atom-number index in big pose.

 			core::kinematics::tree::AtomCOP current_atom = reference_pose.atom_tree().atom_dont_do_update( AtomID(j,i) ).get_self_ptr();

			if ( !get_scratch_atom_id( current_atom_scratch_atom_id, atom_id_map_, current_atom ) ) {
				if ( verbose ){ std::cout << "No current atom id? " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;}
				continue;
			}

			//////////////////////////////////
			// JUMP
			if ( current_atom->is_jump() ) {
				//Special case.

				// Root?
				if ( !current_atom->parent() ) {
					if ( verbose ){ std::cout << "No parent " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;}
					continue;
				}

				core::kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );
				core::kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
				core::kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() );
				if ( !get_scratch_atom_id( input_stub_atom1_scratch_atom_id, atom_id_map_, input_stub_atom1 ) ) {
					if ( verbose ){ std::cout << "No JUMP input_stub_atom1? " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;}
					continue;
				}
				if ( !get_scratch_atom_id( input_stub_atom2_scratch_atom_id, atom_id_map_, input_stub_atom2 ) ) {
					if ( verbose ){ std::cout << "No JUMP input_stub_atom2? " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;}
					continue;
				}
				if ( !get_scratch_atom_id( input_stub_atom3_scratch_atom_id, atom_id_map_, input_stub_atom3 ) ) {
					if ( verbose ){ std::cout << "No JUMP input_stub_atom3? " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;}
					continue;
				}
				Stub const input_stub( scratch_pose_.xyz( input_stub_atom1_scratch_atom_id ),
															 scratch_pose_.xyz( input_stub_atom2_scratch_atom_id ),
															 scratch_pose_.xyz( input_stub_atom3_scratch_atom_id ) );

				core::kinematics::tree::AtomCOP stub_atom1( current_atom->stub_atom1() );
				core::kinematics::tree::AtomCOP stub_atom2( current_atom->stub_atom2() );
				core::kinematics::tree::AtomCOP stub_atom3( current_atom->stub_atom3() );
				if ( !get_scratch_atom_id( stub_atom1_scratch_atom_id, atom_id_map_, stub_atom1 ) ) {
					if ( verbose ){ std::cout << "No JUMP stub_atom1? " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;}
					continue;
				}
				if ( !get_scratch_atom_id( stub_atom2_scratch_atom_id, atom_id_map_, stub_atom2 ) ) {
					if ( verbose ){ std::cout << "No JUMP stub_atom2? " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl; }
					continue;
				}
				if ( !get_scratch_atom_id( stub_atom3_scratch_atom_id, atom_id_map_, stub_atom3 ) ) {
					Size const stub3_rsd = stub_atom3->id().rsd();
					Size const stub3_atmno = stub_atom3->id().atomno();
					if ( verbose ){  std::cout << "No JUMP stub_atom3? " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) <<
							"   stub_atom3 " << pose.residue( stub3_rsd ).name() << stub3_rsd << " " <<	 pose.residue( stub3_rsd ).atom_name( stub3_atmno ) <<
							std::endl;
					}
					continue;
				}
				Stub const stub( scratch_pose_.xyz( stub_atom1_scratch_atom_id ),
												 scratch_pose_.xyz( stub_atom2_scratch_atom_id ),
												 scratch_pose_.xyz( stub_atom3_scratch_atom_id ) );

				Jump const jump( input_stub, stub );

				if ( verbose ) {
					std::cout << "copy_dofs set jump --> " << pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << " " << jump << std::endl;

					std::cout << "input_stub defined by: " << std::endl;
					std::cout << "  " << input_stub_atom1->id().rsd() << " " << pose.residue( input_stub_atom1->id().rsd() ).atom_name( input_stub_atom1->id().atomno() ) << std::endl;
					std::cout << "  " << input_stub_atom2->id().rsd() << " " << pose.residue( input_stub_atom2->id().rsd() ).atom_name( input_stub_atom2->id().atomno() ) << std::endl;
					std::cout << "  " << input_stub_atom3->id().rsd() << " " << pose.residue( input_stub_atom3->id().rsd() ).atom_name( input_stub_atom3->id().atomno() ) << std::endl;
					std::cout << " should match --> " << std::endl;
					// 				std::cout << "  " << input_stub_atom1_scratch_atom_id.rsd() << " " << scratch_pose_.residue( input_stub_atom1_scratch_atom_id.rsd() ).atom_name( input_stub_atom1_scratch_atom_id.atomno() ) << std::endl;
					// 				std::cout << "  " << input_stub_atom2_scratch_atom_id.rsd() << " " << scratch_pose_.residue( input_stub_atom2_scratch_atom_id.rsd() ).atom_name( input_stub_atom2_scratch_atom_id.atomno() ) << std::endl;
					// 				std::cout << "  " << input_stub_atom3_scratch_atom_id.rsd() << " " << scratch_pose_.residue( input_stub_atom3_scratch_atom_id.rsd() ).atom_name( input_stub_atom3_scratch_atom_id.atomno() ) << std::endl;

					std::cout << "stub defined by: " << std::endl;
					std::cout << "  " << stub_atom1->id().rsd() << " " << pose.residue( stub_atom1->id().rsd() ).atom_name( stub_atom1->id().atomno() ) << std::endl;
					std::cout << "  " << stub_atom2->id().rsd() << " " << pose.residue( stub_atom2->id().rsd() ).atom_name( stub_atom2->id().atomno() ) << std::endl;
					std::cout << "  " << stub_atom3->id().rsd() << " " << pose.residue( stub_atom3->id().rsd() ).atom_name( stub_atom3->id().atomno() ) << std::endl;
					std::cout << " should match --> " << std::endl;
					// 				std::cout << "  " << stub_atom1_scratch_atom_id.rsd() << " " << scratch_pose_.residue( stub_atom1_scratch_atom_id.rsd() ).atom_name( stub_atom1_scratch_atom_id.atomno() ) << std::endl;
					// 				std::cout << "  " << stub_atom2_scratch_atom_id.rsd() << " " << scratch_pose_.residue( stub_atom2_scratch_atom_id.rsd() ).atom_name( stub_atom2_scratch_atom_id.atomno() ) << std::endl;
					// 				std::cout << "  " << stub_atom3_scratch_atom_id.rsd() << " " << scratch_pose_.residue( stub_atom3_scratch_atom_id.rsd() ).atom_name( stub_atom3_scratch_atom_id.atomno() ) << std::endl;
					std::cout << "OLD " << pose.jump( AtomID( j, i ) ) << std::endl;
				}

				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////
				// WARNING: NOT COMPLETE YET.
				// Need to check atom_id_domain_map_ to make sure this is jump
				//  is OK
				//////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////
				pose.set_jump( AtomID( j, i ), jump );

				if ( verbose ) std::cout << "NEW " << pose.jump( AtomID( j, i ) ) << std::endl;

				//				pose.dump_pdb( "after_jump_change.pdb" );
				continue;
			}

			//////////////////////////////////
			// D
			core::kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );
			if ( !get_scratch_atom_id( input_stub_atom1_scratch_atom_id, atom_id_map_, input_stub_atom1 ) ) {
				if ( verbose ) { std::cout << "No D input_stub_atom1? " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl; }
				continue;
			}

			if ( check_domain_map( atom_id_domain_map_, current_atom->id(), input_stub_atom1->id() ) ){

				Real const d = ( scratch_pose_.xyz( current_atom_scratch_atom_id ) -
												 scratch_pose_.xyz( input_stub_atom1_scratch_atom_id ) ).length();

				copy_dofs_info_.push_back( std::make_pair( DOF_ID( AtomID( j, i), D), d ) );

			} else {
				if ( verbose ) { std::cout << "D not OK to change? " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl; }
			}


			//////////////////////////////////
			// THETA

			// Following does not generally work..
			//		if (  input_stub_atom1->is_jump() && input_stub_atom1->parent() /* root is special case */ ) {
			//			if ( verbose ) std::cout << "input_stub_atom1 is jump... skipping THETA,PHI --> " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;
			//			continue;
			//		}

			// There's one more case ... if this atom is at a junction between the copied pose and the rest of the structure, that could be bad news.
			// happens with op2, op1 at junctions.
			if ( false ){

				bool problem_with_sister( false );
				for ( Size n = 0; n < input_stub_atom1->n_children(); n++ ) {
					core::kinematics::tree::AtomCOP child_atom( input_stub_atom1->child( n ) );
					if (!get_scratch_atom_id( dummy_atom_id, atom_id_map_, child_atom ) ) {
						if ( verbose ) { std::cout << "No THETA ... Sister atom outside region of interest? " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl; }
						problem_with_sister = true;
						break;
					}
				}

				////////////////////////////////////////////////////////////////////////////////////////////////
				if (problem_with_sister) continue;
			}

			core::kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
			if ( !get_scratch_atom_id( input_stub_atom2_scratch_atom_id, atom_id_map_, input_stub_atom2 ) ) {
				if ( verbose ) std::cout << "No THETA input_stub_atom2? " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;
				continue;
			}

			if ( input_stub_atom2_scratch_atom_id == current_atom_scratch_atom_id /* part of a jump triumvirate*/ ) {
				if ( verbose ) { std::cout << "Part of jump triumvirate, No THETA " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl; }
				continue;
			}

			// To check domain movement, need to compare to grandparent *and* all sisters!
			utility::vector1< AtomID > upstream_atom_ids, downstream_atom_ids;
			downstream_atom_ids.push_back( current_atom->id() );
			upstream_atom_ids.push_back( input_stub_atom2->id() );
			for ( Size n = 0; n < input_stub_atom1->n_children(); n++ ) {
				if( input_stub_atom1->child( n ) == current_atom ) continue;
				upstream_atom_ids.push_back( input_stub_atom1->child( n )->id() );
			}


			//if ( check_domain_map( atom_id_domain_map_, current_atom->id(), input_stub_atom2->id() ) ){
			if ( check_domain_map( atom_id_domain_map_, upstream_atom_ids, downstream_atom_ids ) ){

				Real const theta = angle_radians(
																				 scratch_pose_.xyz( current_atom_scratch_atom_id ),
																				 scratch_pose_.xyz( input_stub_atom1_scratch_atom_id ),
																				 scratch_pose_.xyz( input_stub_atom2_scratch_atom_id ) );

				copy_dofs_info_.push_back( std::make_pair( DOF_ID( AtomID( j, i), THETA), numeric::constants::d::pi - theta ) );


				if ( verbose ) {
					AtomID const & granny = input_stub_atom2->id();
					std::cout << "Good job, change THETA! " << pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << " " << ( atom_id_domain_map_.find( current_atom->id() ) )->second << ";   stub_atom2: " << pose.residue( granny.rsd() ).name1() << granny.rsd()<< " " << pose.residue( granny.rsd() ).atom_name( granny.atomno() )  << " " <<  ( atom_id_domain_map_.find( input_stub_atom2->id() ) )->second << std::endl;
				}

			} else {
				if ( verbose ) { std::cout << "THETA not OK to change? " << pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << ( atom_id_domain_map_.find( current_atom->id() ) )->second  << std::endl; }
			}

			//////////////////////////////////
			// PHI
			//		if ( input_stub_atom2->is_jump() && input_stub_atom2->parent() /* root is special case*/ ) {
			//			if ( verbose ) std::cout << "input_stub_atom2 is jump... skipping PHI --> " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;
			//			continue;
			//		}
			/////////////////////////////////////////////////////////////////

			bool stub3_is_external_sister( false );
			bool stub3_is_external_granny( false );

			core::kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() ); // the real great-grandparent, not some phi_offset garbage.
			core::kinematics::tree::AtomCOP reference_atom( current_atom->input_stub_atom3() );
			core::kinematics::tree::AtomCOP atom_to_move( current_atom );

			if ( input_stub_atom3 == current_atom ) {
				if ( verbose ) { std::cout << "Part of jump triumvirate, No PHI " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl; }
				continue;
			}

			if ( !get_scratch_atom_id( input_stub_atom3_scratch_atom_id, atom_id_map_, input_stub_atom3 ) ||
					 ( input_stub_atom3_scratch_atom_id == current_atom_scratch_atom_id ) /* part of a jump triumvirate*/
					 ) {

				// There are still a couple ways to save this atom.
				// one way: perhaps its stub3 atom is a "sister" lying outside the region of interest,
				//  whose PHI still depends on what's happening *inside* the region of interest
				// An example is the backbone H after protein chainbreaks -- the new pose thinks its PHI is
				// an offset of the C of the previous residue, but it can also be positioned based on
				// the location of C of the next residue.

				if ( !input_stub_atom1->is_jump() &&
						 !input_stub_atom3->is_jump() &&
						 input_stub_atom3->input_stub_atom3() == input_stub_atom1->input_stub_atom2()  &&
						 get_scratch_atom_id( reference_scratch_atom_id, atom_id_map_,
																	input_stub_atom1->input_stub_atom2() ) ){

					reference_atom = input_stub_atom1->input_stub_atom2();

					stub3_is_external_sister = true;
					atom_to_move = input_stub_atom3; // move the current atom by moving the atom from which there's an offset.

					if ( verbose ) { std::cout << "SPECIAL CASE input_stub_atom3 is external sister! " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl; }

				} else { // check for an "aunt" that lives inside the pose

					// Second way. This happened to me for coarse-grained RNA...
					for ( Size n = 0; n < input_stub_atom2->n_children(); n++ ) {
						reference_atom = input_stub_atom2->child( n );
						if( reference_atom == input_stub_atom1 ) continue;
						if ( get_scratch_atom_id( reference_scratch_atom_id, atom_id_map_, reference_atom ) ){
							stub3_is_external_granny = true;
							atom_to_move = current_atom;
							break;
						}
					}

					if ( !stub3_is_external_granny ) {
						// OK, could not salvage this atom's dihedral.
						if ( verbose ) {std::cout << "No PHI input_stub_atom3? " << " skipping " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl; }
						continue;
						//				}
					} else {
						if ( verbose ) { std::cout << "SPECIAL CASE input_stub_atom3 is external granny! " <<  pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl; }
					}

				}

			}

			// What are that atoms whose relative positions would be changed by changing PHI?
			// Downstream: current atom and sisters.
			upstream_atom_ids.clear();
			downstream_atom_ids.clear();
			for ( Size n = 0; n < input_stub_atom1->n_children(); n++ ) {
				if( input_stub_atom1->child( n ) == input_stub_atom3 ) continue;
				downstream_atom_ids.push_back( input_stub_atom1->child( n )->id() );
			}
			if ( stub3_is_external_sister ) downstream_atom_ids.push_back( input_stub_atom3->id() );

			// Upstream: great-grandparent and aunts.
			if ( input_stub_atom2->parent() ) upstream_atom_ids.push_back( input_stub_atom2->parent()->id() );
			//		upstream_atom_ids.push_back( input_stub_atom3->id() ); // this might actually be a sister, but it will move. right?
			for ( Size n = 0; n < input_stub_atom2->n_children(); n++ ) {
				if( input_stub_atom2->child( n ) == input_stub_atom1 ) continue;
				upstream_atom_ids.push_back( input_stub_atom2->child( n )->id() );
			}


			if ( check_domain_map( atom_id_domain_map_, upstream_atom_ids, downstream_atom_ids ) ) {

				if ( stub3_is_external_sister || stub3_is_external_granny ) {

					Real const phi = dihedral_radians(
																						scratch_pose_.xyz( current_atom_scratch_atom_id ),
																						scratch_pose_.xyz( input_stub_atom1_scratch_atom_id ),
																						scratch_pose_.xyz( input_stub_atom2_scratch_atom_id ),
																						scratch_pose_.xyz( reference_scratch_atom_id ) );

					Real const reference_phi = dihedral_radians(
																											reference_pose.xyz( current_atom->id() ),
																											reference_pose.xyz( input_stub_atom1->id() ),
																											reference_pose.xyz( input_stub_atom2->id() ),
																											reference_pose.xyz( reference_atom->id() ) );

					copy_dofs_info_.push_back( std::make_pair( DOF_ID( atom_to_move->id(), PHI),  phi - reference_phi + reference_pose.dof( DOF_ID( atom_to_move->id(), PHI ) ) ) );

					if ( verbose ) std::cout << "Good JOB, CHANGED PHI THROUGH OFFSET!  " << pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;

				} else if (	check_domain_map( atom_id_domain_map_, current_atom->id(), input_stub_atom3->id() ) ){

					Real const phi = dihedral_radians(
																						scratch_pose_.xyz( current_atom_scratch_atom_id ),
																						scratch_pose_.xyz( input_stub_atom1_scratch_atom_id ),
																						scratch_pose_.xyz( input_stub_atom2_scratch_atom_id ),
																						scratch_pose_.xyz( input_stub_atom3_scratch_atom_id ) );

					copy_dofs_info_.push_back( std::make_pair( DOF_ID( AtomID( j, i ), PHI),  phi ) );

					if ( verbose ) std::cout << "Good JOB, CHANGED PHI!  " << pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;

				} else {
					if ( verbose ) {
						AtomID const & great_granny = input_stub_atom3->id();
						std::cout << "special case but current and input_stub_atom3 in same domain? "
											<< pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j )
											<< " " << ( atom_id_domain_map_.find( current_atom->id() ) )->second
											<< " in map? " << ( atom_id_domain_map_.find( current_atom->id() ) != atom_id_domain_map_.end() )
											<< ";   stub_atom3: " << pose.residue( great_granny.rsd() ).name1()
											<< great_granny.rsd()<< " " << pose.residue( great_granny.rsd() ).atom_name( great_granny.atomno() )
											<< " " <<  ( atom_id_domain_map_.find( great_granny ) )->second
											<< " in map? " << ( atom_id_domain_map_.find( great_granny ) != atom_id_domain_map_.end() )
											<< std::endl;
					}
				}

			} else {
				if ( verbose ) {
					std::cout << "PHI not OK to change? " << pose.residue( i ).name1() << i << " " << pose.residue( i ).atom_name( j ) << std::endl;
				}
			}

		} // loop over AtomIDs

		//		pose.dump_pdb( "after_res"+ ObjexxFCL::string_of( i )+".pdb" );

	}

	///////////////////////////////////////////////////////////////////
	bool
	CopyDofs::get_scratch_atom_id( id::AtomID & other_scratch_atom_id,
																 std::map< core::id::AtomID, core::id::AtomID> const & atom_id_map,
																 core::kinematics::tree::AtomCOP other_atom )
	{

		using namespace core::id;

		if ( !other_atom ) return false;

		std::map< AtomID, AtomID >::const_iterator iter( atom_id_map.find( other_atom->id() ) );
		if ( iter == atom_id_map.end() ) return false; // atom not present in scratch pose!

		other_scratch_atom_id = iter->second;

	return true;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	CopyDofs::check_domain_map( std::map< id::AtomID, Size > const & atom_id_domain_map,
										id::AtomID const & atom_id1,
										id::AtomID const & atom_id2 ){

		std::map< id::AtomID, Size >::const_iterator it1 = atom_id_domain_map.find( atom_id1 );
		std::map< id::AtomID, Size >::const_iterator it2 = atom_id_domain_map.find( atom_id2 );

		Size domain1( 999 );
		Size domain2( 999 );

		if ( it1 != atom_id_domain_map.end() ) domain1 = it1->second;
		if ( it2 != atom_id_domain_map.end() ) domain2 = it2->second;

		if ( domain1 == 0)  return true; // domain "0" means OK to change.

		if ( domain2 == 0)  return true; // domain "0" means OK to change.

		//in different domains is OK.    [The only exception is the evil 999 --> code for a totally fixed atom.]
		if ( domain1 < 999 && domain2 < 999 &&  domain1 != domain2 ) return true;

		return false;

	}


	////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	CopyDofs::check_domain_map( std::map< id::AtomID, Size > const & atom_id_domain_map,
										utility::vector1< id::AtomID > const & atom_ids1,
										utility::vector1< id::AtomID > const & atom_ids2 ){

		for ( Size i = 1; i <= atom_ids1.size(); i++ ) {
			for ( Size j = 1; j <= atom_ids2.size(); j++ ) {
				if ( !check_domain_map( atom_id_domain_map, atom_ids1[ i ], atom_ids2[ j ] ) ) return false;
			}
		}
		return true;
	}

} //copy_dofs
} //pose
} //core
