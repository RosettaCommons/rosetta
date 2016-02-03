// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/flexpack/FlexbbRotamerSets.cc
/// @brief
/// @author Florian Richter (floric@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), oct 08


#include <protocols/flexpack/rotamer_set/FlexbbRotamerSets.hh>

#include <protocols/flexpack/interaction_graph/FlexbbInteractionGraph.hh>
#include <protocols/flexpack/interaction_graph/OTFFlexbbInteractionGraph.hh>
#include <protocols/flexpack/interaction_graph/PrecomputedFlexbbInteractionGraph.hh>

/// Project headers
#include <core/conformation/Residue.hh>
#include <core/fragment/Frame.hh>
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <basic/Tracer.hh>

// needed for windows build
#ifdef WIN32
#include <platform/types.hh> // ssize_t
#endif

/// Package headers

#include <core/io/pdb/pdb_writer.hh>
#include <utility/string_util.hh>


/// ObjexxFCL headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray4D.hh>

// C++
#include <fstream>

#include <core/conformation/AbstractRotamerTrie.hh>
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSet.hh>
#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace flexpack {
namespace rotamer_set {

using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer TR( "protocols.flexpack.rotamer_set.FlexbbRotamerSets" );

FlexbbRotamerSets::FlexbbRotamerSets( core::pack::task::PackerTaskCOP task ) :
	task_( task )
{
}

FlexbbRotamerSets::~FlexbbRotamerSets() {}


FlexbbRotamerSetCOP
FlexbbRotamerSets::rotset_for_moltenres( Size molten_resid, Size bbconf /*= 1*/ ) const
{ return rotamers_[ molten_resid ][ bbconf ]; }

FlexbbRotamerSetCOP
FlexbbRotamerSets::rotset_for_residue( Size resid, Size bbconf /*= 1*/ ) const
{ return rotamers_[ resid_2_moltenres_[ resid] ][ bbconf ]; }

core::conformation::ResidueCOP
FlexbbRotamerSets::rotamer_for_residue( Size resid, Size rotindex_on_residue ) const
{ return rotamer_for_moltenres( resid_2_moltenres_[ resid ], rotindex_on_residue ); }


/// @brief function to set up the internal mapping data structures
void
FlexbbRotamerSets::set_frames(
	core::pose::Pose const & pose,
	utility::vector1< core::fragment::FrameCOP > const & frames
)
{

	//to accurately build the different backbone conformations, we need a modifiable copy of the pose
	core::pose::Pose helper_pose = pose;


	core::Size num_frames = frames.size();
	nbbconfs_ = 0;

	flexsegment_span_.resize( num_frames );
	nbbconfs_for_flexseg_.resize( num_frames );


	// Question: when are we setting up the other data structures?
	// either this function or in a separate one (just like set_task in the fixbb RotamerSet)
	// or in a different one? for now, the following segment of code is a duplication of
	// the code in core::pack::rotamer_set::RotamerSets::set_task
	nmoltenres_ = task_->num_to_be_packed();
	total_residue_ = pose.total_residue();

	resid_2_moltenres_.resize( total_residue_ );
	moltenres_2_resid_.resize( nmoltenres_ );
	rotamers_.resize( nmoltenres_ );
	nrotamers_for_moltenres_.resize( nmoltenres_ );
	nrots_for_moltenres_bbconf_.resize( nmoltenres_ );
	nrotoffset_for_moltenres_.resize( nmoltenres_ );
	nrotoffset_for_moltenres_bbconf_.resize( nmoltenres_ );
	moltenres_2_flexseg_.resize( nmoltenres_ );
	bbconf_for_rotamer_of_moltenres_.resize( nmoltenres_ );


	conformations_for_flexible_segments_.resize( nmoltenres_ );

	core::Size count_moltenres(0);
	for ( core::Size ii = 1; ii<= total_residue_; ++ii ) {

		if ( task_->pack_residue( ii ) ) {
			count_moltenres++;
			resid_2_moltenres_[ ii ] = count_moltenres;
			moltenres_2_resid_[ count_moltenres ] = ii;
			moltenres_2_flexseg_[ count_moltenres ] = 0; //will be filled with the correct values later
		} else resid_2_moltenres_[ ii ] = 0;
	}
	assert( count_moltenres == nmoltenres_ );

	//for( utility::vector1< core::fragment::FrameOP >::const_iterator frame_it = frames.begin(); frame_it != frames.end(); ++frame_it ){
	for ( core::Size frame_count = 1; frame_count <= num_frames; ++frame_count ) {

		core::Size cur_frame_start( frames[ frame_count ]->start() ), cur_frame_end( frames[ frame_count]->end() );

		flexsegment_span_[ frame_count ] =  std::pair< core::Size, core::Size>( cur_frame_start, cur_frame_end ) ;
		nbbconfs_for_flexseg_[ frame_count ] = frames[ frame_count ]->nr_frags();
		nbbconfs_ += nbbconfs_for_flexseg_[ frame_count ];

		for ( core::Size frameres_count = cur_frame_start; frameres_count <= cur_frame_end; ++frameres_count ) {

			moltenres_2_flexseg_[ resid_2_moltenres_[ frameres_count ] ] = frame_count;

		}

		//now we need to translate the information in the fragments of this frame  into actual  3D residue coordinates
		for ( core::Size frag_count = 1; frag_count <= frames[ frame_count ]->nr_frags(); ++frag_count ) {
			utility::vector1< core::conformation::ResidueOP > fragment_res;
			build_residue_vector_from_fragment( helper_pose, frames[ frame_count ], frag_count, fragment_res );

			for ( core::Size rescount = cur_frame_start; rescount <= cur_frame_end; ++rescount ) {
				conformations_for_flexible_segments_[ resid_2_moltenres_[ rescount ] ].push_back( fragment_res[ rescount - cur_frame_start + 1] );
			}

		} //iterator over all fragments

	} //iterator over all frames

	//finally, we have to fill the conformations_for_flexseg array for all moltenres that are not part of a frame
	for ( core::Size ii = 1; ii <= nmoltenres_; ++ii ) {

		if ( conformations_for_flexible_segments_[ ii ].size() == 0 ) { //means this was not part of any fragment
			conformations_for_flexible_segments_[ ii ].push_back( core::conformation::ResidueCOP( core::conformation::ResidueOP( new core::conformation::Residue( pose.residue( moltenres_2_resid_[ ii ] ) ) ) ) );
		}
	} //loop over moltenres

} //set_frames

FlexbbRotamerSets::Size
FlexbbRotamerSets::nbackbone_conformations() const
{
	return nbbconfs_;
}


/// @brief function to figure out the flexpack neighbor graph, see core::pack::create_packer_graph for fixbb version
/// @brief problem: Ca/Cb not fixed, so we need to add a certain distance to the interaction energies. this distance
/// @brief will depend on the average deviation of the CBs in the flexible segments for each residue
core::graph::GraphOP
FlexbbRotamerSets::flexpack_neighbor_graph(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::pack::task::PackerTaskCOP task
) const
{

	utility::vector1< core::Distance > residue_radii = core::pack::find_residue_max_radii( pose, task);

	for ( core::Size i = 1; i <= total_residue_; ++i ) {
		if ( resid_2_moltenres_[ i ] != 0 ) residue_radii[ i ] += determine_res_cb_deviation( pose, resid_2_moltenres_[ i ] );
	}

	return core::pack::create_packer_graph( pose, sfxn, task, total_residue_, residue_radii );

} //flexpack_neighbor_graph


/// @brief function to determine the maximum deviation in the position of Cbs for the bbconfs of one residue
/// This should use the "neighbor atom" instead of CB... CB  exist for neither GLY nor non-protein residues.
core::Distance
FlexbbRotamerSets::determine_res_cb_deviation(
	core::pose::Pose const & pose,
	core::Size const moltenres
) const
{

	core::Distance max_sq_dist( 0.0 );
	core::Size const resid = moltenres_2_resid_[ moltenres ];

	core::PointPosition cur_xyz = pose.residue( resid ).xyz( pose.residue( resid ).nbr_atom() );

	//utility::vector1< core::conformation::ResidueCOP > const & confs_this_position = conformations_for_flexible_segments[ resid ];

	for ( utility::vector1< core::conformation::ResidueCOP >::const_iterator res_it = conformations_for_flexible_segments_[ moltenres ].begin();
			res_it != conformations_for_flexible_segments_[ moltenres ].end(); ++res_it ) {

		core::Distance sq_dist = cur_xyz.distance_squared( (*res_it)->xyz( (*res_it)->nbr_atom() ) );

		if ( sq_dist > max_sq_dist ) max_sq_dist = sq_dist;
	}

	return sqrt( max_sq_dist);

} //determine_res_cb_deviation


void
FlexbbRotamerSets::build_rotamers(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::graph::Graph const & flexpack_neighbor_graph
)
{

	for ( core::Size ii = 1; ii <= nmoltenres_; ++ii ) {

		core::Size cur_resid = moltenres_2_resid_[ ii ];
		//for( core::Size jj = 1; jj<= conformations_for_flexible_segments_[ cur_resid ].size(); ++jj ){
		for ( utility::vector1< core::conformation::ResidueCOP >::const_iterator res_it = conformations_for_flexible_segments_[ ii ].begin();
				res_it != conformations_for_flexible_segments_[ ii ].end(); ++res_it ) {

			FlexbbRotamerSetOP rotset( new FlexbbRotamerSet() );
			rotset->set_owner( get_self_weak_ptr() );
			rotset->set_resid( cur_resid );
			rotset->set_existing_residue( *res_it );

			rotset->build_rotamers(pose, sfxn, *task_, flexpack_neighbor_graph.get_self_ptr() );

			rotamers_[ ii ].push_back( rotset );
			//TR << "Built: " << rotset->num_rotamers() << " for moltenres " << ii
			// << " (resid: " << cur_resid << ") ";
			//if ( conformations_for_flexible_segments_[ ii ].size() != 1 ) {
			// TR << " backone # " << rotamers_[ii].size();
			//}
			//TR << std::endl;;
		} //loop over bb conformations for moltenres

	} //loop over moltenres

	update_offset_data() ;


	if ( basic::options::option[basic::options::OptionKeys::packing::dump_rotamer_sets].user() ) {

		this->dump_pdbs( pose, "flex_rotamers");
	}

} //build_rotamers


/// @brief function to dump all the rotamer sets out to pdb files, more or less only used for debugging
/// @brief currently modelled after the core/pack/rotamer_set/RotamerSets dump_pdb function
void
FlexbbRotamerSets::dump_pdbs( core::pose::Pose const & pose, std::string const & filename_base ) const
{

	// model 0 -- just the non-moving residues
	// model N -- Nth rotamer from each set
	using ObjexxFCL::format::I;

	std::string fix_filename = filename_base + "_fixbb.pdb";

	// open file
	std::ofstream base_out( fix_filename.c_str() );

	// write model 0
	Size model(0), atom_counter(0);

	base_out << "MODEL" << I(9,model) << '\n';
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( task_->pack_residue(i) ) continue;
		core::io::pdb::dump_pdb_residue( pose.residue(i), atom_counter, base_out );
	}
	base_out << "ENDMDL\n";


	while ( true ) {
		bool found_a_rotamer( false );
		++model;

		for ( Size ii=1; ii<= nmoltenres_; ++ii ) {
			//Size const resid( moltenres_2_resid[ ii ] );
			FlexbbRotamerSetOP rotset( rotamers_[ ii ][1] );
			if ( rotset->num_rotamers() >= model ) {
				if ( !found_a_rotamer ) {
					found_a_rotamer = true;
					base_out << "MODEL" << I(9,model) << '\n';
				}
				core::io::pdb::dump_pdb_residue( *(rotset->rotamer( model )), atom_counter, base_out );
			}
		} //loop over moltenres
		if ( found_a_rotamer ) {
			base_out << "ENDMDL\n";
		} else {
			break;
		}
	}
	base_out.close();

	//now the stuff for flexible segments
	core::Size curconf(2);

	bool all_confs_covered(false);

	while ( ! all_confs_covered ) {

		std::string cur_filename = filename_base + "_flexconf" + utility::to_string( curconf - 1 ) + ".pdb";
		std::ofstream cur_out( cur_filename.c_str() );

		all_confs_covered = true;

		model = 0;
		atom_counter = 0;

		while ( true ) {

			bool found_a_rotamer( false );
			++model;

			for ( core::Size ii =1; ii <= nmoltenres_; ++ii ) {

				if ( rotamers_[ ii ].size() < curconf ) continue;
				else {

					all_confs_covered = false;
					FlexbbRotamerSetOP rotset( rotamers_[ ii ][ curconf ] );

					if ( rotset->num_rotamers() >= model ) {
						if ( !found_a_rotamer ) {
							found_a_rotamer = true;
							cur_out << "MODEL" << I(9,model) << '\n';
						}
						core::io::pdb::dump_pdb_residue( *(rotset->rotamer( model )), atom_counter, cur_out );
					}
				}

			} //loop over moltenres
			if ( found_a_rotamer ) cur_out << "ENDMDL\n";

			else break;
		} // loop over models

		cur_out.close();
		curconf++;

	} //loop over all confs


} //dump_pdbs

core::uint FlexbbRotamerSets::nrotamers() const { return nrotamers_; }
core::uint FlexbbRotamerSets::nrotamers_for_moltenres( core::uint moltres ) const
{
	return nrotamers_for_moltenres_[ moltres ];
}

core::uint FlexbbRotamerSets::nmoltenres() const
{
	return nmoltenres_;
}

core::uint FlexbbRotamerSets::total_residue() const
{
	return total_residue_;
}

core::uint
FlexbbRotamerSets::moltenres_2_resid( core::uint moltres ) const
{
	return moltenres_2_resid_[ moltres ];
}


core::uint
FlexbbRotamerSets::resid_2_moltenres( core::uint resid ) const
{
	return resid_2_moltenres_[ resid ];
}


core::uint
FlexbbRotamerSets::moltenres_for_rotamer( core::uint rotid ) const
{
	return moltenres_for_rotamer_[ rotid ];
}


core::uint
FlexbbRotamerSets::res_for_rotamer( core::uint rotid ) const
{
	return moltenres_2_resid_[ moltenres_for_rotamer_[ rotid ]];
}


core::conformation::ResidueCOP
FlexbbRotamerSets::rotamer( core::uint rotid ) const
{
	Size moltenres = moltenres_for_rotamer_[ rotid ];
	Size rotid_on_node = rotid - nrotoffset_for_moltenres_[ moltenres ];
	Size bb = bbconf_for_rotamer_of_moltenres_[ moltenres ][ rotid_on_node ];
	Size rotid_on_bb = rotid - nrotoffset_for_moltenres_bbconf_[ moltenres ][ bb ];
	return rotamers_[ moltenres ][ bb ]->rotamer( rotid_on_bb );

}

core::conformation::ResidueCOP
FlexbbRotamerSets::rotamer_for_moltenres( core::uint moltenres_id, core::uint rotamerid ) const
{
	Size bb = bbconf_for_rotamer_of_moltenres_[ moltenres_id ][ rotamerid ];
	Size rotid_on_bb = rotamerid + nrotoffset_for_moltenres_[ moltenres_id ] - nrotoffset_for_moltenres_bbconf_[ moltenres_id ][ bb ];
	return rotamers_[ moltenres_id ][ bb ]->rotamer( rotid_on_bb );
}


core::uint
FlexbbRotamerSets::nrotamer_offset_for_moltenres( core::uint moltres ) const
{
	return nrotoffset_for_moltenres_[ moltres ];
}


/// @details -- unimplemented!
core::uint
FlexbbRotamerSets::rotid_on_moltenresidue( core::uint /*rotid*/ ) const
{
	utility_exit_with_message("UNIMPLEMENTED");
	return 0;
}

/// @brief convert moltenres rotid to id in full rotamer enumeration

core::uint
FlexbbRotamerSets::moltenres_rotid_2_rotid( core::uint /*moltenres*/, core::uint /*moltenresrotid */) const
{
	utility_exit_with_message("UNIMPLEMENTED");
	return 0;
}

void
FlexbbRotamerSets::build_residue_vector_from_fragment(
	core::pose::Pose & pose,
	core::fragment::FrameCOP frame,
	core::Size frag_num,
	utility::vector1< core::conformation::ResidueOP > & fragment_res
)
{

	if ( frame->apply( frag_num, pose ) != frame->length() ) utility_exit_with_message("unknown error when trying to apply a fragment to a pose in setting up FlexbbRotamerSets.");

	for ( core::Size rescount = frame->start(); rescount <= frame->end(); ++rescount ) {

		fragment_res.push_back( core::conformation::ResidueOP( new core::conformation::Residue ( pose.residue( rescount ) ) ) );
	}
} //build residue vector from fragment


core::Size
FlexbbRotamerSets::nbbconfs_for_moltenres( core::Size moltenres ) const
{
	if ( moltenres_2_flexseg_[ moltenres ] == 0 ) return 1;

	else return nbbconfs_for_flexseg_[ moltenres_2_flexseg_[ moltenres ] ];
}

core::Size
FlexbbRotamerSets::nbbconfs_for_res( core::Size resid ) const
{
	if ( resid_2_moltenres_[ resid ] == 0 ) return 1;

	else return nbbconfs_for_moltenres( resid_2_moltenres_[ resid ] );
}

//core::conformation::ResidueCOP
//FlexbbRotamerSets::rotamer_for_moltenres( core::Size moltenres, core::Size rotindex_on_residue ) const
//{
// core::Size bb_conf( bbconf_for_rotamer_of_moltenres_[ moltenres ][ rotindex_on_residue + nrotoffset_for_moltenres_[ moltenres ] ] );
// return rotamers_[ moltenres ][ bb_conf ]->rotamer( rotindex_on_residue - nrotoffset_for_moltenres_bbconf_[ moltenres ][ bb_conf ] );
//}

//core::conformation::ResidueCOP
//FlexbbRotamerSets::rotamer( core::Size rotindex ) const
//{
// core::Size moltenres( moltenres_for_rotamer_[ rotindex ] );
// core::Size bb_conf( bbconf_for_rotamer_of_moltenres_[ moltenres ][ rotindex ] );
// return rotamers_[ moltenres ][ bb_conf ]->rotamer( rotindex - nrotoffset_for_moltenres_bbconf_[ moltenres ][ bb_conf ] );
//}

void
FlexbbRotamerSets::update_offset_data()
{
	nrotamers_ = 0;
	moltenres_for_rotamer_.clear();

	for ( core::Size ii = 1; ii <= nmoltenres_; ++ii ) {
		nrotamers_for_moltenres_[ ii ] = 0;

		nrots_for_moltenres_bbconf_[ ii ].resize( conformations_for_flexible_segments_[ ii ].size() );
		nrotoffset_for_moltenres_bbconf_[ ii ].resize( conformations_for_flexible_segments_[ ii ].size() );
		bbconf_for_rotamer_of_moltenres_[ ii ].clear();

		assert( nrots_for_moltenres_bbconf_[ ii ].size() == rotamers_[ ii ].size() );

		for ( core::Size jj = 1; jj <= rotamers_[ ii ].size(); ++jj ) {

			core::Size cur_numrots = rotamers_[ ii ][ jj ]->num_rotamers();

			nrots_for_moltenres_bbconf_[ ii ][ jj ] = cur_numrots;
			//TR << "nrots_for_moltenres_bbconf_[ "<<ii<<" ][ " <<jj <<" ]  is " << nrots_for_moltenres_bbconf_[ ii ][ jj ] <<  ",   ";
			nrotamers_for_moltenres_[ ii ] += cur_numrots;
			//TR << "nrotamers_for_moltenres_[ " << ii << " ] is " << nrotamers_for_moltenres_[ ii ] << " ;   ";

			for ( core::Size rot = 1; rot <= cur_numrots; ++rot ) {
				moltenres_for_rotamer_.push_back( ii );
				bbconf_for_rotamer_of_moltenres_[ ii ].push_back( jj );
			}

			if ( jj > 1 ) {
				nrotoffset_for_moltenres_bbconf_[ ii ][ jj ] = nrotoffset_for_moltenres_bbconf_[ ii ][ jj-1 ] + nrots_for_moltenres_bbconf_[ ii ][ jj-1 ];
				//TR << "nrotoffset_for_moltenres_bbconf_[ " << ii << " ][ " << jj << " ] is " << nrotoffset_for_moltenres_bbconf_[ ii ][ jj ] << std::endl;
			} else {
				if ( ii > 1 ) {
					nrotoffset_for_moltenres_bbconf_[ ii ][ jj ] = nrotoffset_for_moltenres_[ ii - 1 ] + nrotamers_for_moltenres_[ ii-1 ];
					//TR << "nrotoffset_for_moltenres_bbconf_[ " << ii << " ][ " << jj << " ] is " << nrotoffset_for_moltenres_bbconf_[ ii ][ jj ] << std::endl;
				} else {
					nrotoffset_for_moltenres_bbconf_[ ii ][ jj ] = 0;
					//TR << "nrotoffset_for_moltenres_bbconf_[ " << ii << " ][ " << jj << " ] is " << nrotoffset_for_moltenres_bbconf_[ ii ][ jj ] << std::endl;
				}
			}


		} //iterator over rotamer sets for moltenres

		if ( ii > 1 ) nrotoffset_for_moltenres_[ ii ] = nrotoffset_for_moltenres_[ ii - 1 ] + nrotamers_for_moltenres_[ ii - 1 ] ;
		else nrotoffset_for_moltenres_[ ii ] = 0;
		//TR << "nrotoffset_for_moltenres_[ " << ii << " ] is " << nrotoffset_for_moltenres_[ ii ] << std::endl;

		nrotamers_ += nrotamers_for_moltenres_[ ii ];

	} //loop over moltenres


} //update_offset_data


void
FlexbbRotamerSets::precompute_energies(
	Pose const & pose,
	ScoreFunction const & sfxn,
	core::graph::GraphCOP flexpack_neighbor_graph,
	interaction_graph::FlexbbInteractionGraph & flexbb_ig
) const
{
	using namespace interaction_graph;

	/// Dispatch based on the downcast.
	if ( dynamic_cast< OTFFlexbbInteractionGraph * > ( &flexbb_ig ) ) {
		OTFFlexbbInteractionGraph & otfig =
			static_cast< OTFFlexbbInteractionGraph & > ( flexbb_ig );

		compute_one_body_energies_for_otf_ig( pose, sfxn, flexpack_neighbor_graph, otfig );
	} else {
		assert( dynamic_cast< PrecomputedFlexbbInteractionGraph * > ( &flexbb_ig ) );
		PrecomputedFlexbbInteractionGraph & precomp_ig =
			static_cast< PrecomputedFlexbbInteractionGraph & > ( flexbb_ig );

		precompute_all_energies( pose, sfxn, flexpack_neighbor_graph, precomp_ig );
	}
}

void
FlexbbRotamerSets::precompute_all_energies(
	Pose const & /*pose*/,
	ScoreFunction const & /*sfxn*/,
	core::graph::GraphCOP /*flexpack_neighbor_graph*/,
	interaction_graph::PrecomputedFlexbbInteractionGraph & /*flexbb_ig*/
) const
{
	std::cout << "Yo Mama so fat!" << std::endl;

} // precompute all energies


void
FlexbbRotamerSets::compute_one_body_energies_for_otf_ig(
	Pose const & pose,
	ScoreFunction const & sfxn,
	core::graph::GraphCOP flexpack_neighbor_graph,
	interaction_graph::OTFFlexbbInteractionGraph & flexbb_ig
) const
{
	using namespace utility;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::methods;

	/// 1. Compute energies with the static portions of the structure.
	for ( Size ii = 1; ii <= nmoltenres_; ++ii ) {
		compute_onebody_interactions_with_background( ii, pose, sfxn, flexpack_neighbor_graph, flexbb_ig );
	}

	vector1< vector1< Size > > regular_representatives( nmoltenres_ );
	vector1< vector1< Size > > proline_representatives( nmoltenres_ );
	vector1< vector1< Size > > glycine_representatives( nmoltenres_ );

	/// 2. Collect a set of representative rotamers for each backbone conformation for each moltenresidue.
	/// These representatives will be used for computing bb/bb and bb/sc energies.  They should span the
	/// distinct kinds of backbones; for proteins, these are a. proline backbones (with Npro atom type for
	/// the backbone N) b. glycine backbones (with an extra 2HA where a CB would usually be) and c. regular
	/// backbones (for the other 18 amino acids).  If you wanted to design RNA/DNA hybrid molecules
	///  and the backbone were gaining and loosing the 2' hydroxyl, between RNA and DNA substitutions
	/// during then you would have to modify the code below as well as the code
	/// inside the OTFFlexbbInteracionGraph
	for ( Size ii = 1; ii <= nmoltenres_; ++ii ) {
		Size iinbb = nbbconfs_for_moltenres( ii );
		regular_representatives[ ii ].resize( iinbb );
		proline_representatives[ ii ].resize( iinbb );
		glycine_representatives[ ii ].resize( iinbb );
		std::fill( regular_representatives[ ii ].begin(), regular_representatives[ ii ].end(), 0 );
		std::fill( proline_representatives[ ii ].begin(), proline_representatives[ ii ].end(), 0 );
		std::fill( glycine_representatives[ ii ].begin(), glycine_representatives[ ii ].end(), 0 );
		for ( Size jj = 1; jj <= iinbb; ++jj ) {
			FlexbbRotamerSetCOP jjrotset = rotamers_[ ii ][ jj ];
			Size jjntypes = jjrotset->get_n_residue_groups();
			bool jjregfound( false ), jjprofound( false ), jjglyfound( false );
			for ( Size kk = 1; kk <= jjntypes; ++kk ) {
				Size kkrep = jjrotset->get_residue_type_begin( kk );
				AA kkaa = jjrotset->rotamer( kkrep )->aa();
				if ( ! jjprofound && kkaa == aa_pro ) {
					proline_representatives[ ii ][ jj ] = kkrep;
					jjprofound = true;
				} else if ( ! jjglyfound && kkaa == aa_gly ) {
					glycine_representatives[ ii ][ jj ] = kkrep;
					jjglyfound = true;
				} else if ( ! jjregfound && kkaa != aa_gly && kkaa != aa_pro ) {
					regular_representatives[ ii ][ jj ] = kkrep;
					jjregfound = true;
				}
				if ( jjprofound && jjglyfound && jjregfound ) break;
			}
		}
	}

	for ( Size ii = 1; ii <= nmoltenres_; ++ii ) {
		Size ii_resid = moltenres_2_resid_[ ii ];
		for ( core::graph::Graph::EdgeListConstIter
				li    = flexpack_neighbor_graph->get_node( ii_resid )->const_upper_edge_list_begin(),
				liend = flexpack_neighbor_graph->get_node( ii_resid )->const_upper_edge_list_end();
				li != liend; ++li ) {
			Size const jj_resid = (*li)->get_second_node_ind();
			Size const jj = resid_2_moltenres_[ jj_resid ];

			if ( jj == 0 ) continue;

			compute_sr_one_body_energies_for_flexsets(
				ii, jj,
				regular_representatives,
				proline_representatives,
				glycine_representatives,
				pose, sfxn, flexbb_ig );
		}
	}

	/// Long range interactions
	// Iterate across the long range energy functions and use the iterators generated
	// by the LRnergy container object
	for ( ScoreFunction::LR_2B_MethodIterator
			lr_iter = sfxn.long_range_energies_begin(),
			lr_end  = sfxn.long_range_energies_end();
			lr_iter != lr_end; ++lr_iter ) {
		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*lr_iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue; // only score non-empty energies.
		// Potentially O(N^2) operation...

#ifdef WIN32
		using namespace platform;
#endif

		for ( core::Size ii = 1; ii <= nmoltenres_; ++ ii ) {
			core::Size const ii_resid = moltenres_2_resid_[ ii ];

			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( ii_resid ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii_resid );
					(*rni) != (*rniend); ++(*rni) ) {
				Size const jj_resid = rni->upper_neighbor_id();

				core::Size const jj = resid_2_moltenres_[ jj_resid ];
				if ( ii > jj ) continue; // compute against upper, moltenres neighbors only
				assert( jj != 0 );

				if ( ! flexbb_ig.get_edge_exists( ii, jj ) ) {
					flexbb_ig.add_edge( ii, jj );
				}
				//std::cout << "flexbb_ig.note_long_range_interactions_exist_for_edge( " << ii << ", " << jj << ");" << std::endl;
				flexbb_ig.note_long_range_interactions_exist_for_edge( ii, jj );
			}
		}
	}


}

/// @details Computes the backbone/backbone and backbone/sidechain energies for a
/// pair of molten residues.  The interactions are separated into three kinds of backbone
/// types: regular backbones, proline backbones and glycine backbones.  The different
/// solvation for proline's backbone nitrogen makes the proline correction necessary;
/// the new carbon-hbond term makes the glycine correction necessary.
void
FlexbbRotamerSets::compute_sr_one_body_energies_for_flexsets(
	Size lowermoltenres,
	Size uppermoltenres,
	utility::vector1< utility::vector1< Size > > const & regular_representatives,
	utility::vector1< utility::vector1< Size > > const & proline_representatives,
	utility::vector1< utility::vector1< Size > > const & glycine_representatives,
	Pose const & pose,
	ScoreFunction const & sfxn,
	interaction_graph::OTFFlexbbInteractionGraph & flexbb_ig
) const
{
	//std::cout << "ONE BODY ENERGIES: " << lowermoltenres  << " " << uppermoltenres << std::endl;
	///bool sought_pair( lowermoltenres == 27 && uppermoltenres == 28 );
	bool sought_pair = false;

	flexbb_ig.add_edge( lowermoltenres, uppermoltenres );

	utility::vector1< Size > const & lregrep( regular_representatives[ lowermoltenres ] );
	utility::vector1< Size > const & uregrep( regular_representatives[ uppermoltenres ] );
	utility::vector1< Size > const & lprorep( proline_representatives[ lowermoltenres ] );
	utility::vector1< Size > const & uprorep( proline_representatives[ uppermoltenres ] );
	utility::vector1< Size > const & lglyrep( glycine_representatives[ lowermoltenres ] );
	utility::vector1< Size > const & uglyrep( glycine_representatives[ uppermoltenres ] );

	bool const sameflexseg = flexsegid_for_moltenres( lowermoltenres ) == flexsegid_for_moltenres( uppermoltenres );

	Size const lnbb = nbbconfs_for_moltenres( lowermoltenres );
	Size const unbb = nbbconfs_for_moltenres( uppermoltenres );

	Size const lnrots = nrotamers_for_moltenres_[ lowermoltenres ];
	Size const unrots = nrotamers_for_moltenres_[ uppermoltenres ];

	Size const BBREG = 1;
	Size const BBPRO = 2;
	Size const BBGLY = 3;
	Size const NBBCLASSES = 3;

	/// 1. backbone backbone energies.
	FArray4D< PackerEnergy > bbbb_energies( NBBCLASSES, NBBCLASSES, lnbb, unbb, 0.0 );

	utility::vector1< Size > lower_bbrepresentatives( NBBCLASSES );
	utility::vector1< Size > upper_bbrepresentatives( NBBCLASSES );
	for ( Size ii = 1; ii <= lnbb; ++ii ) {
		FlexbbRotamerSetCOP iirotset( rotamers_[ lowermoltenres ][ ii ] );
		//std::fill( lower_bbrepresentatives.begin(), lower_bbrepresentatives.end(), 0 );
		lower_bbrepresentatives[ BBREG ] = lregrep[ ii ];
		lower_bbrepresentatives[ BBPRO ] = lprorep[ ii ];
		lower_bbrepresentatives[ BBGLY ] = lglyrep[ ii ];
		for ( Size jj = 1; jj <= unbb; ++jj ) {
			FlexbbRotamerSetCOP jjrotset( rotamers_[ uppermoltenres ][ jj ] );
			//std::fill( upper_bbrepresentatives.begin(), upper_bbrepresentatives.end(), 0 );
			upper_bbrepresentatives[ BBREG ] = uregrep[ jj ];
			upper_bbrepresentatives[ BBPRO ] = uprorep[ jj ];
			upper_bbrepresentatives[ BBGLY ] = uglyrep[ jj ];

			for ( Size kk = 1; kk <= NBBCLASSES; ++kk ) {
				if ( lower_bbrepresentatives[ kk ] == 0 ) continue;
				core::conformation::Residue const & kkres( * (iirotset->rotamer( lower_bbrepresentatives[ kk ] )) );
				for ( Size ll = 1; ll <= NBBCLASSES; ++ll ) {
					if ( upper_bbrepresentatives[ ll ] == 0 ) continue;
					core::conformation::Residue const & llres( * (jjrotset->rotamer( upper_bbrepresentatives[ ll ] )) );
					bbbb_energies( kk, ll, ii, jj ) =
						core::pack::rotamer_set::RotamerSets::get_bb_bbE( pose, sfxn, kkres, llres );
					if ( sought_pair ) {
						std::cout << "bbbbE: " << ii  << " " << jj << " " << kk << " " << ll << " " << bbbb_energies( kk, ll, ii, jj ) << std::endl;
					}
				}
			}
		}
	}

	if ( sought_pair ) {
		for ( Size ii = 1; ii <= lnbb; ++ii ) {

			FlexbbRotamerSetCOP iirotset( rotamers_[ lowermoltenres ][ ii ] );
			for ( Size jj = 1; jj <= unbb; ++jj ) {
				FlexbbRotamerSetCOP jjrotset( rotamers_[ uppermoltenres ][ jj ] );

				for ( Size kk = 1; kk <= iirotset->get_n_residue_types(); ++kk ) {
					Size kkrep = iirotset->get_residue_type_begin( kk );
					core::conformation::Residue const & kkrot( * iirotset->rotamer( kkrep ));
					for ( Size ll = 1; ll <= jjrotset->get_n_residue_types(); ++ll ) {
						Size llrep = jjrotset->get_residue_type_begin( ll );
						core::conformation::Residue const & llrot(* jjrotset->rotamer( llrep ));

						std::cout << "BBBB: " << ii << " " << jj << " " <<
							kkrot.aa() << " " << llrot.aa() << " e: " <<
							core::pack::rotamer_set::RotamerSets::get_bb_bbE( pose, sfxn, kkrot, llrot ) << std::endl;
					}
				}
			}
		}
	}

	{ // scope
		/// Sidechain/backbone energies for the lower-residue's rotamers.
		utility::vector1< PackerEnergy > scbb_energies_lower( NBBCLASSES );
		for ( Size ii = 1; ii <= lnrots; ++ii ) {
			Size iibb = bbconf_for_rotamer_of_moltenres_[ lowermoltenres ][ ii ];
			core::conformation::Residue const & iirot( *( rotamer_for_moltenres( lowermoltenres, ii )));
			Size iibbclass = ( iirot.aa() == core::chemical::aa_pro ? BBPRO : iirot.aa() == core::chemical::aa_gly ? BBGLY : BBREG );

			for ( Size jj = 1, jje = (sameflexseg ? 1 : unbb); jj <= jje; ++jj ) {
				Size jjbb = sameflexseg ? iibb : jj;
				FlexbbRotamerSetCOP jjrotset( rotamers_[ uppermoltenres ][ jjbb ] );
				upper_bbrepresentatives[ BBREG ] = uregrep[ jjbb ];
				upper_bbrepresentatives[ BBPRO ] = uprorep[ jjbb ];
				upper_bbrepresentatives[ BBGLY ] = uglyrep[ jjbb ];
				std::fill( scbb_energies_lower.begin(), scbb_energies_lower.end(), 0.0f );
				for ( Size kk = 1; kk <= NBBCLASSES; ++kk ) {
					if ( upper_bbrepresentatives[ kk ] == 0 ) continue;
					core::conformation::Residue const & kkrot( *(jjrotset->rotamer( upper_bbrepresentatives[ kk ] )) );
					scbb_energies_lower[ kk ] =
						core::pack::rotamer_set::RotamerSets::get_sc_bbE( pose, sfxn, iirot, kkrot );
					if ( sought_pair ) {
						std::cout << "lower bbscE: " << ii  << " " << jj << " " << kk << " " << scbb_energies_lower[ kk ] << std::endl;
					}

				}
				if ( upper_bbrepresentatives[ BBREG ] != 0  &&
						upper_bbrepresentatives[ BBPRO ] == 0 &&
						upper_bbrepresentatives[ BBGLY ] == 0 ) {
					// Neither a proline correction nor a glycine correction.
					// Set the bb/bb and bb/sc energies.
					flexbb_ig.set_ProCorrection_values_for_edge( lowermoltenres, uppermoltenres,
						lowermoltenres, ii, jjbb,
						bbbb_energies( iibbclass, BBREG, iibb, jjbb ), 0.0,
						scbb_energies_lower[ BBREG ], 0.0 );
				}
				if ( upper_bbrepresentatives[ BBPRO ] != 0 ) {
					flexbb_ig.set_ProCorrection_values_for_edge( lowermoltenres, uppermoltenres,
						lowermoltenres, ii, jjbb,
						bbbb_energies( iibbclass, BBREG, iibb, jjbb ), bbbb_energies( iibbclass, BBPRO, iibb, jjbb ),
						scbb_energies_lower[ BBREG ],                  scbb_energies_lower[ BBPRO ] );

				}
				if ( upper_bbrepresentatives[ BBGLY ] != 0 ) {
					flexbb_ig.set_GlyCorrection_values_for_edge( lowermoltenres, uppermoltenres,
						lowermoltenres, ii, jjbb,
						bbbb_energies( iibbclass, BBREG, iibb, jjbb ), bbbb_energies( iibbclass, BBGLY, iibb, jjbb ),
						scbb_energies_lower[ BBREG ],                  scbb_energies_lower[ BBGLY ] );
				}
			}
		}
	}

	{ // scope
		/// Sidechain/backbone energies for the upper-residue's rotamers.
		utility::vector1< PackerEnergy > scbb_energies_upper( NBBCLASSES );
		for ( Size ii = 1; ii <= unrots; ++ii ) {
			Size iibb = bbconf_for_rotamer_of_moltenres_[ uppermoltenres ][ ii ];
			core::conformation::Residue const & iirot( *( rotamer_for_moltenres( uppermoltenres, ii )));
			Size iibbclass = ( iirot.aa() == core::chemical::aa_pro ? BBPRO : iirot.aa() == core::chemical::aa_gly ? BBGLY : BBREG );

			for ( Size jj = 1, jje = sameflexseg ? 1 : lnbb; jj <= jje; ++jj ) {
				Size jjbb = sameflexseg ? iibb : jj;
				FlexbbRotamerSetCOP jjrotset( rotamers_[ lowermoltenres ][ jjbb ] );
				lower_bbrepresentatives[ BBREG ] = lregrep[ jjbb ];
				lower_bbrepresentatives[ BBPRO ] = lprorep[ jjbb ];
				lower_bbrepresentatives[ BBGLY ] = lglyrep[ jjbb ];
				std::fill( scbb_energies_upper.begin(), scbb_energies_upper.end(), 0.0 );

				for ( Size kk = 1; kk <= NBBCLASSES; ++kk ) {
					if ( lower_bbrepresentatives[ kk ] == 0 ) continue;
					core::conformation::Residue const & kkrot( *(jjrotset->rotamer( lower_bbrepresentatives[ kk ] )) );
					scbb_energies_upper[ kk ] =
						core::pack::rotamer_set::RotamerSets::get_sc_bbE( pose, sfxn, iirot, kkrot );
					if ( sought_pair ) {
						std::cout << "upper bbscE: " << ii  << " " << jj << " " << kk << " " << scbb_energies_upper[ kk ] <<
							" aa " << iirot.aa() <<
							" coord: " << iirot.xyz( iirot.nheavyatoms() ).x() <<
							" " << iirot.xyz( iirot.nheavyatoms() ).y() <<
							" " << iirot.xyz( iirot.nheavyatoms() ).z() <<  std::endl;
					}
				}
				if ( lower_bbrepresentatives[ BBREG ] != 0  &&
						lower_bbrepresentatives[ BBPRO ] == 0 &&
						lower_bbrepresentatives[ BBGLY ] == 0 ) {
					// Neither a proline correction nor a glycine correction.
					// Set the bb/bb and bb/sc energies.
					flexbb_ig.set_ProCorrection_values_for_edge( lowermoltenres, uppermoltenres,
						uppermoltenres, ii, jjbb,
						bbbb_energies( BBREG, iibbclass, jjbb, iibb ), 0.0,
						scbb_energies_upper[ BBREG ], 0.0 );
				}
				if ( lower_bbrepresentatives[ BBPRO ] != 0 ) {
					flexbb_ig.set_ProCorrection_values_for_edge( lowermoltenres, uppermoltenres,
						uppermoltenres, ii, jjbb,
						bbbb_energies( BBREG, iibbclass, jjbb, iibb ), bbbb_energies( BBPRO, iibbclass, jjbb, iibb ),
						scbb_energies_upper[ BBREG ],                  scbb_energies_upper[ BBPRO ] );

				}
				if ( lower_bbrepresentatives[ BBGLY ] != 0 ) {
					flexbb_ig.set_GlyCorrection_values_for_edge( lowermoltenres, uppermoltenres,
						uppermoltenres, ii, jjbb,
						bbbb_energies( BBREG, iibbclass, jjbb, iibb ), bbbb_energies( BBGLY, iibbclass, jjbb, iibb ),
						scbb_energies_upper[ BBREG ],                   scbb_energies_upper[ BBGLY ] );
				}
			}
		}
	}
}

void
FlexbbRotamerSets::compute_onebody_interactions_with_background(
	Size moltenres,
	Pose const & pose,
	ScoreFunction const & sfxn,
	core::graph::GraphCOP flexpack_neighbor_graph,
	interaction_graph::FlexbbInteractionGraph & flexbb_ig
) const
{


	utility::vector1< PackerEnergy > all_one_body_energies( nrotamers_for_moltenres( moltenres ), 0.0 );
	for ( Size ii = 1, iie = nbbconfs_for_moltenres( moltenres ); ii <= iie; ++ii ) {


		FlexbbRotamerSetCOP iirotset( rotamers_[ moltenres ][ ii ] );

		/// Build the rotamer tries now so that they are ready for rotamer/bg
		/// energy calculations.  Delete the tries at the end of this operation!
		/// Warning -- this must absolutely be rethought when precomputing all
		/// energies upfront.
		sfxn.prepare_rotamers_for_packing( pose, const_cast< FlexbbRotamerSet & > (*iirotset) );

		Size const iinrots = iirotset->num_rotamers();
		utility::vector1< PackerEnergy > ii_one_body_energies( iinrots );
		iirotset->compute_one_body_energies(
			pose, sfxn, *task_, flexpack_neighbor_graph, ii_one_body_energies );
		Size const ii_offset = nrotoffset_for_moltenres_bbconf_[ moltenres ][ ii ] - nrotoffset_for_moltenres_[ moltenres ];
		for ( Size jj = 1; jj <= iinrots; ++jj ) {
			//if ( moltenres == 17 && jj == 11 ) { std::cout << " OneBodyEnergy: " << jj << " " << ii_offset << " " << jj+ii_offset << " " << ii_one_body_energies[ jj ] << std::endl; }
			//if ( moltenres == 27 ) { std::cout << "OneBodyEnergy " << moltenres << " " << jj + ii_offset << " " << ii_one_body_energies[ jj ] << std::endl;}
			all_one_body_energies[ jj + ii_offset ] = ii_one_body_energies[ jj ];
		}

		/// For the sake of reducing memory load, clear the tries out as soon as
		/// the one body energies are finished being calculated.
		for ( Size jj = 1; jj <= core::scoring::methods::n_energy_methods; ++jj ) {
			(const_cast< FlexbbRotamerSet & > (*iirotset)).store_trie( jj, 0 );
		}
	}

	flexbb_ig.add_to_nodes_one_body_energy( moltenres, all_one_body_energies );

}

void
FlexbbRotamerSets::show( std::ostream & out ) const {
	out << "FlexbbRotamerSets with " << nmoltenres_ << " molten residues for " << total_residue_ << " total residues and " << nrotamers_ << " rotamers." << std::endl;
}


}
}
}



#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::flexpack::rotamer_set::FlexbbRotamerSets::FlexbbRotamerSets() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::flexpack::rotamer_set::FlexbbRotamerSets::save( Archive & arc ) const {
	arc( cereal::base_class< core::pack::rotamer_set::RotamerSetsBase >( this ) );
	arc( CEREAL_NVP( nmoltenres_ ) ); // Size
	arc( CEREAL_NVP( total_residue_ ) ); // Size
	arc( CEREAL_NVP( nbbconfs_ ) ); // Size
	arc( CEREAL_NVP( task_ ) ); // PackerTaskCOP
	arc( CEREAL_NVP( rotamers_ ) ); // utility::vector1<utility::vector1<rotamer_set::FlexbbRotamerSetOP> >
	arc( CEREAL_NVP( nrotamers_ ) ); // Size
	arc( CEREAL_NVP( nrotamers_for_moltenres_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( moltenres_2_resid_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( resid_2_moltenres_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( moltenres_for_rotamer_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( bbconf_for_rotamer_of_moltenres_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( conformations_for_flexible_segments_ ) ); // utility::vector1<utility::vector1<core::conformation::ResidueCOP> >
	arc( CEREAL_NVP( flexsegment_span_ ) ); // utility::vector1<std::pair<Size, Size> >
	arc( CEREAL_NVP( nbbconfs_for_flexseg_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( moltenres_2_flexseg_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( nrotoffset_for_moltenres_ ) ); // utility::vector1<Size>
	arc( CEREAL_NVP( nrots_for_moltenres_bbconf_ ) ); // utility::vector1<utility::vector1<Size> >
	arc( CEREAL_NVP( nrotoffset_for_moltenres_bbconf_ ) ); // utility::vector1<utility::vector1<Size> >
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::flexpack::rotamer_set::FlexbbRotamerSets::load( Archive & arc ) {
	arc( cereal::base_class< core::pack::rotamer_set::RotamerSetsBase >( this ) );
	arc( nmoltenres_ ); // Size
	arc( total_residue_ ); // Size
	arc( nbbconfs_ ); // Size
	std::shared_ptr< core::pack::task::PackerTask > local_task;
	arc( local_task ); // PackerTaskCOP
	task_ = local_task; // copy the non-const pointer(s) into the const pointer(s)
	arc( rotamers_ ); // utility::vector1<utility::vector1<rotamer_set::FlexbbRotamerSetOP> >
	arc( nrotamers_ ); // Size
	arc( nrotamers_for_moltenres_ ); // utility::vector1<Size>
	arc( moltenres_2_resid_ ); // utility::vector1<Size>
	arc( resid_2_moltenres_ ); // utility::vector1<Size>
	arc( moltenres_for_rotamer_ ); // utility::vector1<Size>
	arc( bbconf_for_rotamer_of_moltenres_ ); // utility::vector1<utility::vector1<Size> >
	utility::vector1< utility::vector1< std::shared_ptr< core::conformation::Residue > > > local_conformations_for_flexible_segments;
	arc( local_conformations_for_flexible_segments ); // utility::vector1<utility::vector1<core::conformation::ResidueCOP> >
	conformations_for_flexible_segments_ = local_conformations_for_flexible_segments; // copy the non-const pointer(s) into the const pointer(s)
	arc( flexsegment_span_ ); // utility::vector1<std::pair<Size, Size> >
	arc( nbbconfs_for_flexseg_ ); // utility::vector1<Size>
	arc( moltenres_2_flexseg_ ); // utility::vector1<Size>
	arc( nrotoffset_for_moltenres_ ); // utility::vector1<Size>
	arc( nrots_for_moltenres_bbconf_ ); // utility::vector1<utility::vector1<Size> >
	arc( nrotoffset_for_moltenres_bbconf_ ); // utility::vector1<utility::vector1<Size> >
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::flexpack::rotamer_set::FlexbbRotamerSets );
CEREAL_REGISTER_TYPE( protocols::flexpack::rotamer_set::FlexbbRotamerSets )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_flexpack_rotamer_set_FlexbbRotamerSets )
#endif // SERIALIZATION
