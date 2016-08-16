// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   SetupNCSMover.cc
/// @brief  Sets up NCS restraints
/// @author Frank DiMaio

// Unit headers
#include <protocols/simple_moves/symmetry/SetupNCSMoverCreator.hh>
#include <protocols/simple_moves/symmetry/SetupNCSMover.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>

#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/constraints/DihedralPairConstraint.hh>
#include <core/scoring/constraints/DistancePairConstraint.hh>
#include <core/scoring/func/TopOutFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>

#include <basic/Tracer.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace simple_moves {
namespace symmetry {

static THREAD_LOCAL basic::Tracer TZ( "protocols.simple_moves.symmetry.SetupNCSMover" );

// creators
std::string
SetupNCSMoverCreator::keyname() const {
	return SetupNCSMoverCreator::mover_name();
}

protocols::moves::MoverOP
SetupNCSMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetupNCSMover );
}

std::string
SetupNCSMoverCreator::mover_name() {
	return "SetupNCS";
}

////////////////////
////////////////////

NCSResMapping::NCSResMapping( core::pose::Pose &pose ) {
	ngroups_ = 0;
	nres_ = pose.total_residue();
}

utility::vector1< core::Size >
NCSResMapping::get_equiv( core::Size resid ) {
	utility::vector1< core::Size > retval(1,resid);
	if ( resid <= nres_ ) { // tolerate residues that may have been added after NCS was applied
		for ( core::Size i = 1; i <= ngroups_; ++i ) {
			if ( mapping_[i][resid] != 0 ) retval.push_back( mapping_[i][resid] );
		}
	}
	return retval;
}

core::Size
NCSResMapping::get_equiv( core::Size groupID, core::Size resid ) {
	if ( groupID > ngroups_ ) return 0;
	return mapping_[groupID][resid];
}

void
NCSResMapping::set_equiv( core::Size groupID, core::Size res1, core::Size res2 ) {
	runtime_assert( res1 <= nres_ && res2 <= nres_ );

	if ( groupID > ngroups_ ) {
		mapping_.resize(groupID,utility::vector1< core::Size >(nres_,0));
		ngroups_ = groupID;
	}

	mapping_[groupID][res1] = res2;
	mapping_[groupID][res2] = res1;
}


////////////////////
////////////////////

SetupNCSMover::SetupNCSMover() : protocols::moves::Mover("SetupNCSMover") {
	set_defaults();
}

SetupNCSMover::SetupNCSMover( std::string src, std::string tgt ) : protocols::moves::Mover("SetupNCSMover") {
	add_group( src, tgt );
	set_defaults();
}

SetupNCSMover::SetupNCSMover( std::string src, utility::vector1<std::string> tgt ): protocols::moves::Mover("SetupNCSMover") {
	for ( Size i=1; i<=tgt.size(); ++i ) {
		add_group( src, tgt[i] );
	}
	set_defaults();
}

SetupNCSMover::SetupNCSMover( utility::vector1<std::string> src, utility::vector1<std::string> tgt ) : protocols::moves::Mover("SetupNCSMover") {
	runtime_assert( src.size() == tgt.size() );
	for ( Size i=1; i<=tgt.size(); ++i ) {
		add_group( src[i], tgt[i] );
	}
	set_defaults();
}

SetupNCSMover::~SetupNCSMover(){}

void SetupNCSMover::set_defaults() {
	wt_ = 0.01;  //fpd this seems reasonable in fullatom
	limit_ = 10.0;  // in degrees
	bb_ = chi_ = true;
	symmetric_sequence_ = false;
	sd_ = 1000; // for Harmonic, distance. Value determined empirically after test with relaxed structures
	distance_pair_ = false;
}

void SetupNCSMover::add_group( std::string src, std::string tgt ) {
	src_.push_back( src );
	tgt_.push_back( tgt );
}

void SetupNCSMover::add_groupE( std::string src, std::string tgt ) {
	srcE_.push_back( src );
	tgtE_.push_back( tgt );
}

void SetupNCSMover::add_groupD( std::string src, std::string tgt ) {
	srcD_.push_back( src );
	tgtD_.push_back( tgt );
}


void SetupNCSMover::apply( core::pose::Pose & pose ) {
	using namespace std;
	using namespace utility;
	using namespace core;
	using namespace core::id;
	using namespace core::scoring::constraints;

	assert( src_.size() == tgt_.size() );

	// store NCS data in pose datacache
	NCSResMappingOP ncs;
	if ( pose.data().has( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ) {
		ncs = ( utility::pointer::static_pointer_cast< protocols::simple_moves::symmetry::NCSResMapping > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING ) ));
	}
	ncs = NCSResMappingOP( new NCSResMapping(pose) );
	pose.data().set( core::pose::datacache::CacheableDataType::NCS_RESIDUE_MAPPING, ncs );

	// map residue ranges -> resid pairs (using PDBinfo if necessary)
	for ( Size i=1; i<=src_.size(); ++i ) {
		utility::vector1<Size> src_i = core::pose::get_resnum_list_ordered( src_[i], pose );
		utility::vector1<Size> tgt_i = core::pose::get_resnum_list_ordered( tgt_[i], pose );

		runtime_assert( src_i.size() == tgt_i.size() );
		if ( src_i.size() == 0 ) utility_exit_with_message("Error creating NCS constraints: no groups found");

		for ( Size j=1; j<=src_i.size(); ++j ) {
			core::Size resnum_src = src_i[j], resnum_tgt = tgt_i[j];
			id::AtomID id_a1,id_a2,id_a3,id_a4, id_b1,id_b2,id_b3,id_b4;

			if ( TZ.Debug.visible() ) {
				if ( j<=10 ) {
					TZ.Debug << "Add constraint " << resnum_src << " : " << resnum_tgt << std::endl;
				} else if ( j==11 ) {
					TZ.Debug << " ... " << src_i.size()-10 << "  omitted" << std::endl;
				}
			}

			ncs->set_equiv( i, resnum_src, resnum_tgt );

			//replace residue identity to match reference positions
			if ( symmetric_sequence_ ) {
				pose.replace_residue( resnum_tgt, pose.residue(resnum_src), true );
			}

			if ( bb_ ) {
				// fpd better safe than sorry
				runtime_assert (pose.residue(resnum_src).mainchain_torsions().size() == pose.residue(resnum_tgt).mainchain_torsions().size());
				for ( Size k=1, k_end = pose.residue(resnum_src).mainchain_torsions().size(); k<= k_end; ++k ) {
					id::TorsionID tors_src(resnum_src,BB,k);
					id::TorsionID tors_tgt(resnum_tgt,BB,k);

					pose.conformation().get_torsion_angle_atom_ids(tors_src,  id_a1, id_a2, id_a3, id_a4 );
					pose.conformation().get_torsion_angle_atom_ids(tors_tgt,  id_b1, id_b2, id_b3, id_b4 );

					if ( id_a1.atomno()==0 || id_a2.atomno()==0 || id_a3.atomno()==0 || id_a4.atomno()==0 ) continue;
					if ( id_b1.atomno()==0 || id_b2.atomno()==0 || id_b3.atomno()==0 || id_b4.atomno()==0 ) continue;

					// make the cst
					pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new DihedralPairConstraint( id_a1, id_a2, id_a3, id_a4, id_b1, id_b2, id_b3, id_b4, core::scoring::func::FuncOP( new core::scoring::func::TopOutFunc( wt_, 0.0, limit_ ) ) ) ) ) );
				}
			}

			if ( chi_ ) {
				if ( pose.residue(resnum_src).aa() != pose.residue(resnum_tgt).aa() ) {
					TZ.Error << "Trying to constrain sidechain torsions of different residue types!" << std::endl;
					TZ.Error << " >>> " << resnum_src << " vs " << resnum_tgt << std::endl;
					continue;
				}

				// cys/cyd messes this up I think
				if ( pose.residue(resnum_src).nchi() != pose.residue(resnum_tgt).nchi() ) {
					TZ.Error << "Trying to constrain sidechains with different number of chi angles!" << std::endl;
					continue;
				}

				for ( Size k=1, k_end = pose.residue(resnum_src).nchi(); k<= k_end; ++k ) {
					id::TorsionID tors_src(resnum_src,CHI,k);
					id::TorsionID tors_tgt(resnum_tgt,CHI,k);

					pose.conformation().get_torsion_angle_atom_ids(tors_src,  id_a1, id_a2, id_a3, id_a4 );
					pose.conformation().get_torsion_angle_atom_ids(tors_tgt,  id_b1, id_b2, id_b3, id_b4 );

					if ( id_a1.atomno()==0 || id_a2.atomno()==0 || id_a3.atomno()==0 || id_a4.atomno()==0 ) continue;
					if ( id_b1.atomno()==0 || id_b2.atomno()==0 || id_b3.atomno()==0 || id_b4.atomno()==0 ) continue;

					// make the cst
					pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new DihedralPairConstraint( id_a1, id_a2, id_a3, id_a4, id_b1, id_b2, id_b3, id_b4, core::scoring::func::FuncOP( new core::scoring::func::TopOutFunc( wt_, 0.0, limit_ ) ) ) ) ) );
				}
			}
		}
	}


	//symmetrize ONLY sequence (eg. N- and C-term)
	if ( symmetric_sequence_ ) {
		assert( srcE_.size() == tgtE_.size() );

		// map residue ranges -> resid pairs (using PDBinfo if necessary)
		for ( Size i=1; i<=srcE_.size(); ++i ) {
			utility::vector1<Size> srcE_i = core::pose::get_resnum_list_ordered( srcE_[i], pose );
			utility::vector1<Size> tgtE_i = core::pose::get_resnum_list_ordered( tgtE_[i], pose );
			runtime_assert( srcE_i.size() == tgtE_i.size() );

			if ( srcE_i.size() == 0 ) {
				utility_exit_with_message("Error creating NCS constraints: " + srcE_[i] +" : "+tgtE_[i]);
			}

			for ( Size j=1; j<=srcE_i.size(); ++j ) {
				core::Size resnum_srcE = srcE_i[j], resnum_tgtE = tgtE_i[j];

				TZ.Debug << "Symmetrizing residues " << resnum_srcE << " <--> " << resnum_tgtE << std::endl;

				//replace residue identity to match reference positions
				pose.replace_residue( resnum_tgtE, pose.residue(resnum_srcE), true );
			}
		}
	}


	//calculate distance pairing constraints
	if ( distance_pair_ ) {
		assert( srcD_.size() == tgtD_.size() );

		core::Real test1=srcD_.size();
		core::Real test2=tgtD_.size();
		TZ.Debug << "Size src " << test1 << " Size tgt " << test2 << std::endl;


		// map residue ranges -> resid pairs (using PDBinfo if necessary)
		for ( Size i=1; i<=srcD_.size(); ++i ) {
			utility::vector1<Size> srcD_i = core::pose::get_resnum_list_ordered( srcD_[i], pose );

			core::Real temp_d=srcD_i.size();
			TZ.Debug << "src 1 size " << temp_d << std::endl;

			utility::vector1<Size> tgtD_i = core::pose::get_resnum_list_ordered( tgtD_[i], pose );
			runtime_assert( srcD_i.size() == tgtD_i.size() );
			runtime_assert( srcD_i.size() % 2 == 0 );


			if ( srcD_i.size() == 0 ) {
				utility_exit_with_message("Error creating NCS distance pair constraints: " + srcD_[i] + " : " + tgtD_[i] );
			}

			for ( Size j=1; j<=srcD_i.size()/2; ++j ) {
				core::Size resnum_src1 = srcD_i[j], resnum_tgt1 = tgtD_i[j], resnum_src2 = srcD_i[j+srcD_i.size()/2], resnum_tgt2 = tgtD_i[j+srcD_i.size()/2];
				id::AtomID id_ad1,id_ad2, id_bd1,id_bd2;

				TZ.Debug << "Add distance pair constraint " << resnum_src1 << " - " << resnum_src2 << " <--> " << resnum_tgt1 << " - " << resnum_tgt2 << std::endl;

				//get the ca atoms
				id_ad1 = id::AtomID(pose.residue(resnum_src1).atom_index("CA") , resnum_src1);
				id_ad2 = id::AtomID(pose.residue(resnum_src2).atom_index("CA") , resnum_src2);
				id_bd1 = id::AtomID(pose.residue(resnum_tgt1).atom_index("CA") , resnum_tgt1);
				id_bd2 = id::AtomID(pose.residue(resnum_tgt2).atom_index("CA") , resnum_tgt2);

				// make the cst
				pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new DistancePairConstraint( id_ad1, id_ad2, id_bd1, id_bd2, core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc( 0.0, sd_ ) ) ) ) ) ); // for Harmonic

			}
		}
	}
}

void SetupNCSMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	filters::Filters_map const & /*filters*/,
	moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/ ) {
	if ( tag->hasOption( "bb" ) ) bb_ = tag->getOption< bool >( "bb" );
	if ( tag->hasOption( "chi" ) ) chi_ = tag->getOption< bool >( "chi" );
	if ( tag->hasOption( "wt" ) ) wt_ = tag->getOption< core::Real >( "wt" );
	if ( tag->hasOption( "limit" ) ) limit_ = tag->getOption< core::Real >( "limit" );
	if ( tag->hasOption( "symmetric_sequence" ) ) symmetric_sequence_ = tag->getOption< bool >( "symmetric_sequence" );
	if ( tag->hasOption( "sd" ) ) sd_ = tag->getOption< core::Real >( "sd" ); // for Harmonic, distance
	if ( tag->hasOption( "distance_pair" ) ) distance_pair_ = tag->getOption< bool >( "distance_pair" );

	// now parse ncs groups <<< subtags
	utility::vector1< TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< TagCOP >::const_iterator tag_it;
	for ( tag_it = branch_tags.begin(); tag_it!=branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() == "NCSgroup" || (*tag_it)->getName() == "ncsgroup" ) {
			std::string src_i = (*tag_it)->getOption<std::string>( "source" );
			std::string tgt_i = (*tag_it)->getOption<std::string>( "target" );
			if ( symmetric_sequence_ ) {
				TZ.Debug << "Symmetrizing sequences " << src_i << " -> " << tgt_i << std::endl;
			}
			TZ.Debug << "Adding NCS cst " << src_i << " -> " << tgt_i << std::endl;
			add_group( src_i, tgt_i );
		}

		if ( (*tag_it)->getName() == "NCSend" || (*tag_it)->getName() == "ncsend" ) {
			std::string srcE_i = (*tag_it)->getOption<std::string>( "source" );
			std::string tgtE_i = (*tag_it)->getOption<std::string>( "target" );
			if ( symmetric_sequence_ ) {
				TZ.Debug << "Symmetrizing sequences at ends" << srcE_i << " -> " << tgtE_i << std::endl;
			}
			add_groupE( srcE_i, tgtE_i );
		}

		//group in the format <NCSdistance source="2A,9A" target="28A,35A"/> or
		//          <NCSdistance source="2A-20A,42A-60A" target="128A-156A,168A-186A"/>
		if ( (*tag_it)->getName() == "NCSdistance" || (*tag_it)->getName() == "ncsdistance" ) {
			std::string srcD_i = (*tag_it)->getOption<std::string>( "source" );
			std::string tgtD_i = (*tag_it)->getOption<std::string>( "target" );
			add_groupD( srcD_i, tgtD_i );
		}

	}
}

std::string
SetupNCSMover::get_name() const {
	return SetupNCSMoverCreator::mover_name();
}


} // symmetry
} // moves
} // protocols

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::simple_moves::symmetry::NCSResMapping::NCSResMapping() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::simple_moves::symmetry::NCSResMapping::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( mapping_ ) ); // utility::vector1<utility::vector1<core::Size> >
	arc( CEREAL_NVP( ngroups_ ) ); // core::Size
	arc( CEREAL_NVP( nres_ ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::simple_moves::symmetry::NCSResMapping::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( mapping_ ); // utility::vector1<utility::vector1<core::Size> >
	arc( ngroups_ ); // core::Size
	arc( nres_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::simple_moves::symmetry::NCSResMapping );
CEREAL_REGISTER_TYPE( protocols::simple_moves::symmetry::NCSResMapping )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_simple_moves_symmetry_SetupNCSMover )
#endif // SERIALIZATION
