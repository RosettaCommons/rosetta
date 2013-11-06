// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_moves/GreenPacker.cc
/// @brief  packing mover that makes extensive reuse of rotamer pair energies class definition
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Unit headers
#include <protocols/simple_moves/GreenPacker.hh>

/// Core headers
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>

/// Numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

/// Utility headers
#include <utility/vector1.functions.hh>

/// C++ Headers
// AUTO-REMOVED #include <ctime>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

/// @details Auto-generated virtual destructor
MinimalRotamer::~MinimalRotamer() {}

basic::Tracer TR("protocols.simple_moves.GreenPacker");

MinimalRotamer::MinimalRotamer( core::conformation::Residue const & res ) :
	residue_type_( res.type() ),
	ideal_geometry_( false )
{
	if ( has_ideal_geometry( res ) ) {
		ideal_geometry_ = true;
		record_chi( res );
	} else {
		ideal_geometry_ = false;
		record_internal_geometry( res );
	}
}

bool
MinimalRotamer::has_ideal_geometry(
	core::conformation::Residue const & res
) const
{
	for ( Size ii = res.type().first_sidechain_atom(); ii <= res.type().nheavyatoms(); ++ii ) {
		if ( ! atom_is_ideal( res, ii ) ) {
			return false;
		}
		for ( Size jj = res.type().attached_H_begin( ii );
				jj <= res.type().attached_H_end( ii ); ++jj ) {
			if ( ! atom_is_ideal( res, jj ) ) {
				return false;
			}
		}
	}
	return true;
}

/// @details Assumption: sidechain ideal coordinates do not depend on other residues
/// ( as opposed to the backbone O for example, which depends on the coordinates of i+1.)
bool
MinimalRotamer::atom_is_ideal(
	core::conformation::Residue const & res,
	Size const atom_id
) const
{
	core::chemical::AtomICoor const & atom_icoor( res.type().icoor( atom_id ) );
	Size st1 = atom_icoor.stub_atom1().atomno();
	Size st2 = atom_icoor.stub_atom2().atomno();
	Size st3 = atom_icoor.stub_atom2().atomno();

	Real const ideal_d = atom_icoor.d();
	if ( std::abs( res.xyz( atom_id ).distance( res.xyz( st1 ) ) - ideal_d ) > 1e-8 ) {
		//TR << "Failed atom is ideal: " << res.name() << " " << res.atom_name( atom_id ) << " d: ";
		//TR << res.xyz( atom_id ).distance( res.xyz( st1 ) ) << " ideal: " << ideal_d << std::endl;
		return false;
	}

	/// It turns out that the CB placement on non-ideal backbones is almost uniformly non-ideal.
	/// It might be nice to special-case CB, describing the sidechain geometry as a (set of?)
	/// backbone-to-sidechain non-ideal-geometry (ies), and chi-angles for ideal geometry thereafter.
	/// The profiler doesn't show much time spent in determining the correspondence between
	/// rotamer sets... looping on GreenPacker->apply, 0.0% is spent in
	/// GreenPacker::find_current_and_original_rotamer_correspondence.  The time spent in "same"
	/// is clearly minimal, and so the only gain of switching to the above-described representation
	/// would be a memory savings...
	Real const ideal_theta = atom_icoor.theta();
	Real const actual_theta = numeric::constants::d::pi - numeric::angle_radians(
		res.xyz( atom_id ),
		res.xyz( st1 ),
		res.xyz( st2 ) );
	if ( std::abs( actual_theta - ideal_theta ) > 1e-8 ) {
		//TR << "Failed atom is ideal: " << res.name() << " " << res.atom_name( atom_id ) << " theta: ";
		//TR << actual_theta << " ideal: " << ideal_theta << std::endl;
		return false;
	}

	bool atom_is_last_for_some_chi = false;
	for ( Size ii = 1; ii <= res.type().nchi(); ++ii ) {
		if ( res.type().chi_atoms( ii )[ 4 ] == atom_id ) {
			atom_is_last_for_some_chi = true;
			break;
		}
	}
	if ( ! atom_is_last_for_some_chi ) {
		Real const ideal_phi = atom_icoor.phi();
		Real const actual_phi = numeric::dihedral_radians(
			res.xyz( atom_id ),
			res.xyz( st1 ),
			res.xyz( st2 ),
			res.xyz( st3 )
		);
		if ( std::abs( basic::periodic_range( actual_phi - ideal_phi, numeric::constants::d::pi ) ) > 1e-8 ) {

			//TR << "Failed atom is ideal: " << res.name() << " " << res.atom_name( atom_id ) << " phi: ";
			//TR << actual_phi << " ideal: " << ideal_phi << " diff: " << basic::periodic_range( actual_phi - ideal_phi, numeric::constants::d::pi ) << std::endl;

			return false;
		}
	}
	return true;
}

void
MinimalRotamer::record_chi( core::conformation::Residue const & res )
{
	chi_.resize( res.type().nchi() );
	for ( Size ii = 1; ii <= chi_.size(); ++ii ) {
		assert( chi_matches_coords( res, ii ));
		chi_[ ii ] = res.chi( ii );
	}
}

/// @details only record internal geometry for sidechain atoms
void
MinimalRotamer::record_internal_geometry( core::conformation::Residue const & res )
{
	internal_geometry_.resize( res.natoms(), Vector(0,0,0) );
	for ( Size ii = res.type().first_sidechain_atom(); ii <= res.type().nheavyatoms(); ++ii ) {
		record_internal_geometry( res, ii );
		for ( Size jj = res.type().attached_H_begin( ii );
				jj <= res.type().attached_H_end( ii ); ++jj ) {
			record_internal_geometry( res, jj );
		}
	}
}

void
MinimalRotamer::record_internal_geometry(
	core::conformation::Residue const & res,
	Size const atom_id
)
{
	core::chemical::AtomICoor const & atom_icoor( res.type().icoor( atom_id ) );
	Size st1 = atom_icoor.stub_atom1().atomno();
	Size st2 = atom_icoor.stub_atom2().atomno();
	Size st3 = atom_icoor.stub_atom2().atomno();

	Real actual_d = res.xyz( atom_id ).distance( res.xyz( st1 ) );
	internal_geometry_[ atom_id ][ d ] = actual_d;

	/// Don't bother subtracing pi -- just don't try to construct a rotamer using this "proper" theta
	/// (as opposed to the canonical improper bond angle...)
	Real const actual_theta = numeric::angle_radians(
		res.xyz( atom_id ),
		res.xyz( st1 ),
		res.xyz( st2 ) );
	internal_geometry_[ atom_id ][ theta ] = actual_theta;

	Real const actual_phi = numeric::dihedral_radians(
		res.xyz( atom_id ),
		res.xyz( st1 ),
		res.xyz( st2 ),
		res.xyz( st3 )
	);
	internal_geometry_[ atom_id ][ phi ] = actual_phi;
}

bool
MinimalRotamer::chi_matches_coords( core::conformation::Residue const & res, Size chi_index ) const
{
	core::chemical::AtomIndices const & chi_atoms( res.type().chi_atoms( chi_index ) );
	assert( chi_atoms.size() == 4 );
	Real const actual_chi = numeric::dihedral_degrees(
		res.xyz( chi_atoms[ 1 ] ),
		res.xyz( chi_atoms[ 2 ] ),
		res.xyz( chi_atoms[ 3 ] ),
		res.xyz( chi_atoms[ 4 ] )
	);
	if ( std::abs( basic::periodic_range( actual_chi - res.chi( chi_index ), 180 ) ) > 1e-8 ) {
		return false;
	}
	return true;
}

bool
MinimalRotamer::same( MinimalRotamer const & other ) const
{
	//TR << "same ? " << residue_type_.name() << " " << other.residue_type_.name() << " " << ideal_geometry_ << " " << other.ideal_geometry_ << std::endl;


	/// pointer comparison
	if ( & residue_type_ != & other.residue_type_ ) return false;

	if ( ideal_geometry_ && other.ideal_geometry_ ) {
		return same_chi( other );
	} else if ( ! ideal_geometry_ && ! other.ideal_geometry_ ) {
		return same_nonideal_geometry( other );
	}

	/// else, one is ideal and the other is not
	return false;
}

///  @detail This tolerance may need fiddling with
bool
MinimalRotamer::same_chi( MinimalRotamer const & other ) const
{
	assert( ideal_geometry_ && other.ideal_geometry_ );
	assert( chi_.size() == other.chi_.size() );
	for ( Size ii = 1; ii <= chi_.size(); ++ii ) {
		if ( std::abs( basic::periodic_range( chi_[ ii ] - other.chi_[ ii ], 180 ) ) > 1e-8 ) {
			//TR << "same chi failed: " << residue_type_.name() << " " << other.residue_type_.name() << " " << ii << " " << chi_[ ii ] << " " << other.chi_[ ii ] << " " << basic::periodic_range( chi_[ ii ] - other.chi_[ ii ], 180 ) << std::endl;
			return false;
		}
	}
	//TR << "SAME!" << std::endl;
	return true;
}

///  @detail These tolerances may need fiddling with
bool
MinimalRotamer::same_nonideal_geometry( MinimalRotamer const & other ) const
{
	assert( !ideal_geometry_ && !other.ideal_geometry_ );
	assert( internal_geometry_.size() == other.internal_geometry_.size() );
	for ( Size ii = 1; ii <= internal_geometry_.size(); ++ii ) {
		if ( std::abs(  internal_geometry_[ ii ][ d ] - other.internal_geometry_[ ii ][ d ] ) > 1e-8 ) {
			//TR << "same nonideal geometry d failed: " << residue_type_.name() << " " << other.residue_type_.name() << " " << ii << " " << internal_geometry_[ ii ][ d ] << " " << other.internal_geometry_[ ii ][ d ]  << " " << internal_geometry_[ ii ][ d ] - other.internal_geometry_[ ii ][ d ] << std::endl;
			return false;
		}
		if ( std::abs( internal_geometry_[ ii ][ theta ] - other.internal_geometry_[ ii ][ theta ] ) > 1e-8 ) {
			//TR << "same nonideal geometry theta failed: " << residue_type_.name() << " " << other.residue_type_.name() << " " << ii << " " << internal_geometry_[ ii ][ theta ] << " " << other.internal_geometry_[ ii ][ theta ]  << " " << internal_geometry_[ ii ][ theta ] - other.internal_geometry_[ ii ][ theta ] << std::endl;
			return false;
		}

		if ( std::abs( basic::periodic_range(
				internal_geometry_[ ii ][ phi ] - other.internal_geometry_[ ii ][ phi ],
				numeric::constants::d::pi ) ) > 1e-8 ) {
			//TR << "same nonideal geometry phi failed: " << residue_type_.name() << " " << other.residue_type_.name() << " " << ii << " " << internal_geometry_[ ii ][ phi ] << " " << other.internal_geometry_[ ii ][ phi ]  << " " <<
			//	basic::periodic_range( internal_geometry_[ ii ][ phi ] - other.internal_geometry_[ ii ][ phi ], numeric::constants::d::pi ) << std::endl;
			return false;
		}
	}
	//TR << "SAME!" << std::endl;
	return true;
}


bool
MinimalRotamer::same_residue_type(
	MinimalRotamer const & other
) const
{
	//TR << "same residue type ? " << residue_type_.name() << " " << other.residue_type_.name() << std::endl;
	return ( & residue_type_ == & other.residue_type_ );
}

core::chemical::AA
MinimalRotamer::aa() const
{
	return residue_type_.aa();
}

///// Group Discriminators

GroupDiscriminator::~GroupDiscriminator() {}

ChainGroupDiscriminator::~ChainGroupDiscriminator() {}

GroupDiscriminatorOP ChainGroupDiscriminator::clone() const
{
	return new ChainGroupDiscriminator;
}

core::Size
ChainGroupDiscriminator::group_id( Pose const & pose, Size seqpos ) const
{
	return pose.residue( seqpos ).chain();
}


	//////////////////////////////////////////////////////////////////////
UserDefinedGroupDiscriminator::~UserDefinedGroupDiscriminator() {}

GroupDiscriminatorOP UserDefinedGroupDiscriminator::clone() const
{
	return new UserDefinedGroupDiscriminator;
}

core::Size
UserDefinedGroupDiscriminator::group_id( Pose const & , Size seqpos ) const {
	return group_ids_[ seqpos ];
}

void
UserDefinedGroupDiscriminator::set_group_ids( utility::vector1< core::Size > const & group_ids_input ) {
	group_ids_ = group_ids_input;
}

/// Green Packer

GreenPacker::GreenPacker() : create_reference_data_( true ) {}

GreenPacker::~GreenPacker() {}


void
GreenPacker::apply( core::pose::Pose & pose )
{
	clock_t starttime = clock();
	if ( create_reference_data_ ) {
		setup_reference_data( pose );
	}
	repack( pose );
	clock_t stoptime = clock();
	TR << "Green Packer took: " << ((double) stoptime-starttime) / CLOCKS_PER_SEC << " seconds." << std::endl;
}

std::string
GreenPacker::get_name() const {
	return "GreenPacker";
}

void
GreenPacker::set_scorefunction( ScoreFunction const & sfxn )
{
	full_sfxn_ = sfxn.clone();

	/// create context independent and context dependent versions of this score function
	ci_sfxn_ = new ScoreFunction;
	cd_sfxn_ = new ScoreFunction;

	set_weights_for_sfxn( *ci_sfxn_, full_sfxn_->ci_2b_types(),    full_sfxn_->weights() );
	set_weights_for_sfxn( *ci_sfxn_, full_sfxn_->ci_1b_types(),    full_sfxn_->weights() );
	set_weights_for_sfxn( *ci_sfxn_, full_sfxn_->ci_lr_2b_types(), full_sfxn_->weights() );

	set_weights_for_sfxn( *cd_sfxn_, full_sfxn_->cd_2b_types(),    full_sfxn_->weights() );
	set_weights_for_sfxn( *cd_sfxn_, full_sfxn_->cd_1b_types(),    full_sfxn_->weights() );
	set_weights_for_sfxn( *cd_sfxn_, full_sfxn_->cd_lr_2b_types(), full_sfxn_->weights() );

	create_reference_data_ = true;
}

void
GreenPacker::set_weights_for_sfxn(
	ScoreFunction & sfxn,
	ScoreTypes const & scoretypes,
	EnergyMap const & weights
) const
{
	for ( Size ii = 1; ii <= scoretypes.size(); ++ii ) {
		sfxn.set_weight( scoretypes[ ii ], weights[ scoretypes[ ii ] ] );
	}
}

void
GreenPacker::set_group_discriminator( protocols::simple_moves::GroupDiscriminatorOP discriminator )
{
	group_discriminator_ = discriminator;
}


/// @details Pointer copy, atm, since there's no good way to clone task factories
/// and all of their contents.
void
GreenPacker::set_task_factory( TaskFactoryOP factory )
{
	task_factory_ = factory;
}

/// @details Pointer copy, atm, since there's no good way to clone task factories
/// and all of their contents.
void
GreenPacker::set_reference_round_task_factory( TaskFactoryOP factory )
{
	reference_task_factory_ = factory;
}



void
GreenPacker::setup_reference_data( core::pose::Pose & pose )
{
	split_pose_into_groups( pose );
	create_reference_packer_task( pose );
	create_reference_packer_neighbor_graph( pose );
	create_reference_rotamers( pose );
	compute_reference_intragroup_rpes( pose );
	create_reference_data_ = false; // next time, we'll reuse these energies
}

void
GreenPacker::repack( core::pose::Pose & pose )
{
	split_pose_into_groups( pose );
	create_fresh_task( pose );
	create_fresh_packer_neighbor_graph( pose );
	create_fresh_rotamers( pose );
	find_reference_and_current_rotamer_correspondence( pose );
	compute_energies( pose );
	run_sa( pose );
	cleanup();
}


void
GreenPacker::split_pose_into_groups(
	core::pose::Pose & pose
)
{
	group_ids_.resize( pose.total_residue() );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		group_ids_[ ii ] = group_discriminator_->group_id( pose, ii );
	}
	Size max_group_id = utility::max( group_ids_ );
	group_members_.resize( max_group_id + 1 );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		group_members_[ group_ids_[ ii ] ].push_back( ii );
	}
}

void
GreenPacker::create_reference_packer_task(
	core::pose::Pose & pose
)
{
	reference_task_ = reference_task_factory_->create_task_and_apply_taskoperations( pose );

	core::pack::pack_scorefxn_pose_handshake( pose, *full_sfxn_ );
	pose.update_residue_neighbors();
	ci_sfxn_->setup_for_packing( pose, reference_task_->repacking_residues(), reference_task_->designing_residues() );

}

// @details So that bump check does not throw out rotamers that collide only with residues
// in other groups, throw out edges from the packer neighbor graph that span groups.
void
GreenPacker::create_reference_packer_neighbor_graph(
	core::pose::Pose & pose
)
{
	reference_packer_neighbor_graph_ = core::pack::create_packer_graph( pose, *ci_sfxn_, reference_task_ );
	drop_inter_group_edges( pose, reference_packer_neighbor_graph_ );
}

void
GreenPacker::create_reference_rotamers(
	core::pose::Pose & pose
)
{
	using namespace core::pack::rotamer_set;

	reference_rotamer_sets_ = new RotamerSets;
	reference_rotamer_sets_->set_task( reference_task_ );
	reference_rotamer_sets_->build_rotamers( pose, *ci_sfxn_, reference_packer_neighbor_graph_ );

	reference_rotamer_sets_->prepare_sets_for_packing( pose, *ci_sfxn_ );


	/// Now create images of the internal geometry for each of the rotamers created.
	original_rotamers_.resize( pose.total_residue() );

	reference_moltenres_2_resid_.resize( reference_rotamer_sets_->nmoltenres() );
	reference_resid_2_moltenres_.resize( pose.total_residue() );
	std::fill( reference_resid_2_moltenres_.begin(), reference_resid_2_moltenres_.end(), 0 );

	for ( Size ii = 1; ii <= reference_rotamer_sets_->nmoltenres(); ++ii ) {
		Size const ii_resid = reference_rotamer_sets_->moltenres_2_resid( ii );
		Size const ii_nrots = reference_rotamer_sets_->rotamer_set_for_moltenresidue( ii )->num_rotamers();
		original_rotamers_[ ii_resid ].reserve( ii_nrots );

		reference_moltenres_2_resid_[ ii ] = ii_resid;
		reference_resid_2_moltenres_[ ii_resid ] = ii;

		for ( Size jj = 1; jj <= ii_nrots; ++jj ) {
			original_rotamers_[ ii_resid ].push_back(
				new protocols::simple_moves::MinimalRotamer(
				*reference_rotamer_sets_->rotamer_set_for_moltenresidue( ii )->rotamer( jj ) ) );
		}
	}
}

void
GreenPacker::compute_reference_intragroup_rpes(
	core::pose::Pose & pose
)
{
	using namespace core::pack::interaction_graph;
	InteractionGraphBaseOP ig = InteractionGraphFactory::create_interaction_graph(
		*reference_task_, *reference_rotamer_sets_, pose, *ci_sfxn_ );

	PrecomputedPairEnergiesInteractionGraphOP pig(
		dynamic_cast< PrecomputedPairEnergiesInteractionGraph * > (
		ig.get() ));

	/// if the dynamic cast failed, then the packer task has produced
	/// an on the fly interaction graph (or some other non-precomputed IG
	/// that hadn't been invented by 4/18/2008) which makes the GreenPacker
	/// pointless... just because it's pointless doesn't mean it should fail,
	/// though, so this code should be revised in the future.
	if ( ! pig ) {
		utility_exit_with_message( "GreenPacker asked to use on-the-fly interaction graph.  Why?" );
	}
	pig->initialize( *reference_rotamer_sets_ );

	reference_rotamer_sets_->precompute_two_body_energies(
		pose,
		*ci_sfxn_,
		reference_packer_neighbor_graph_,
		pig
	);

	ci_rpes_ = pig;

	/// get rid of unneeded data
	reference_task_ = 0;
	reference_rotamer_sets_ = 0;
	reference_packer_neighbor_graph_ = 0;

}

void
GreenPacker::create_fresh_task(
	core::pose::Pose & pose
)
{
	current_task_ = task_factory_->create_task_and_apply_taskoperations( pose );

	core::pack::pack_scorefxn_pose_handshake( pose, *full_sfxn_ );
	pose.update_residue_neighbors();
	full_sfxn_->setup_for_packing( pose, current_task_->repacking_residues(), current_task_->designing_residues() );

}

void
GreenPacker::create_fresh_packer_neighbor_graph(
	core::pose::Pose & pose
)
{
	current_packer_neighbor_graph_ = core::pack::create_packer_graph( pose, *full_sfxn_, current_task_ );
	current_inter_group_packer_neighbor_graph_ = new Graph( *current_packer_neighbor_graph_ );
	current_intra_group_packer_neighbor_graph_ = new Graph( *current_packer_neighbor_graph_ );
	drop_intra_group_edges( pose, current_inter_group_packer_neighbor_graph_ );
	drop_inter_group_edges( pose, current_intra_group_packer_neighbor_graph_ );
}

/// @details All edges to group 0 are considered inter-group edges.
void
GreenPacker::drop_inter_group_edges(
	core::pose::Pose & pose,
	core::graph::GraphOP packer_neighbor_graph
) const
{
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		Size const ii_grp = group_ids_[ ii ];
		for ( Graph::EdgeListIter
				iru  = packer_neighbor_graph->get_node(ii)->upper_edge_list_begin(),
				irue = packer_neighbor_graph->get_node(ii)->upper_edge_list_end();
				iru != irue; /* no increment */ ) {
			Graph::EdgeListIter irunext = iru;
			++irunext;
			if ( ii_grp != group_ids_[ (*iru)->get_second_node_ind() ] ||
				ii_grp == 0 /* implies other group also 0*/ ) {
				packer_neighbor_graph->delete_edge( *iru ); // invalidates iterator
			}
			iru = irunext; // increment;
		}
	}
}

/// @details All edges to group 0 are considered inter-group edges.
void
GreenPacker::drop_intra_group_edges(
	core::pose::Pose & pose,
	core::graph::GraphOP packer_neighbor_graph
) const
{
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		Size const ii_grp = group_ids_[ ii ];
		for ( Graph::EdgeListIter
				iru  = packer_neighbor_graph->get_node(ii)->upper_edge_list_begin(),
				irue = packer_neighbor_graph->get_node(ii)->upper_edge_list_end();
				iru != irue; /* no increment */ ) {
			Graph::EdgeListIter irunext = iru;
			++irunext;
			if ( ii_grp == group_ids_[ (*iru)->get_second_node_ind() ] &&
					ii_grp != 0 &&
					group_ids_[ (*iru)->get_second_node_ind() ] != 0 ) {
				packer_neighbor_graph->delete_edge( *iru ); // invalidates iterator
			}
			iru = irunext; // increment;
		}
	}
}


void
GreenPacker::create_fresh_rotamers(
	core::pose::Pose & pose
)
{
	current_rotamer_sets_ = new RotamerSets;
	current_rotamer_sets_->set_task( current_task_ );
	current_rotamer_sets_->build_rotamers( pose, *full_sfxn_, current_packer_neighbor_graph_ );
	current_rotamer_sets_->prepare_sets_for_packing( pose, *full_sfxn_ );

	/// Now create images of the internal geometry for each of the rotamers created.
	current_rotamers_.resize( pose.total_residue() );

	for ( Size ii = 1; ii <= current_rotamer_sets_->nmoltenres(); ++ii ) {
		Size const ii_resid = current_rotamer_sets_->moltenres_2_resid( ii );
		Size const ii_nrots = current_rotamer_sets_->rotamer_set_for_moltenresidue( ii )->num_rotamers();
		current_rotamers_[ ii_resid ].reserve( ii_nrots );

		for ( Size jj = 1; jj <= ii_nrots; ++jj ) {
			current_rotamers_[ ii_resid ].push_back(
				new protocols::simple_moves::MinimalRotamer(
				*current_rotamer_sets_->rotamer_set_for_moltenresidue( ii )->rotamer( jj ) ) );
		}
	}
}

void
GreenPacker::find_reference_and_current_rotamer_correspondence(
	core::pose::Pose & pose
)
{
	using namespace core::pack::rotamer_set;
	initialize_internal_correspondence_data( pose );

	/// Assumption: if rotamer i of the original rotamers and rotamer j of the current rotamers match,
	/// then if rotamers i+1 and rotamers j+1 do not match and i and i+1 have the same amino acid type,
	/// then i+1 will not match any rotamer in the current rotamer set.
	/// O(N) algorithm.  Optimal for O(NlgN) with a sort.
	Size correspondences_found( 0 );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		Size const ii_moltenresid = current_rotamer_sets_->resid_2_moltenres( ii );
		//TR << "find correspondences for residue " << ii << " " << ii_moltenresid << std::endl;
		if ( ii_moltenresid == 0 ) continue;

		RotamerSetCOP iirotset = current_rotamer_sets_->rotamer_set_for_residue( ii );

		Size const ii_orig_nrots = orig_rot_2_curr_rot_[ ii ].size();
		Size const ii_curr_nrots = curr_rot_2_orig_rot_[ ii ].size();

		//TR << "#orig rots: " << ii_orig_nrots << " #curr_rots: " << ii_curr_nrots << std::endl;

		Size orig_id( 1 ), curr_id( 1 );
		bool orig_prev_restype_match( false );
		Size ii_correspondences_found( 0 );
		while ( orig_id <= ii_orig_nrots && curr_id <= ii_curr_nrots ) {
			protocols::simple_moves::MinimalRotamerOP orig_rot = original_rotamers_[ ii ][ orig_id ];
			protocols::simple_moves::MinimalRotamerOP curr_rot = current_rotamers_[ ii ][ curr_id ];
			if ( orig_rot->same_residue_type( *curr_rot )) {
				if ( orig_rot->same( *curr_rot ) ) {
					curr_rot_2_orig_rot_[ ii ][ curr_id ] = orig_id;
					orig_rot_2_curr_rot_[ ii ][ orig_id ] = curr_id;
					//TR << "Correspondence found: increment both" << std::endl;
					++orig_id;
					++curr_id;
					++ii_correspondences_found;
				} else {
					//TR << "Increment original" << std::endl;
					++orig_id;
				}
				orig_prev_restype_match = true;
			} else {
				if ( orig_prev_restype_match ) {
					if ( original_rotamers_[ ii ][ orig_id - 1 ]->same_residue_type( *curr_rot ) ) {
						// orig has moved on to a new amino acid, and curr is on the same amino acid,
						// no more matches expected of the curr amino acid type to the rotamers in orig
						//TR << "Increment current" << std::endl;
						++curr_id;
						continue;
					}
				}

				/// We've arrived in a tricky situation where orig and curr don't match residue-types,
				/// and orig_prev was never found to have a restype match.  The correspondence from
				/// this point forward is not guaranteed optimal.  It will definately be suboptimal if
				/// ever curr_rots[ ii ]->aa() > curr_rots[ ii + 1]->aa().
				/// Such a situation will only arise if you're appending rotamers to a rotamer set
				/// in an out-of-order fashion... e.g. after adding rots for all 20 aa's, you add a few
				/// extra arginine rots, then these rots won't be picked up as matching.
				/// To avoid this situation from slowing your code, append extra rotamers inline:
				/// add extra arginine rotamers as soon as the rotamer set has built its canonical arginine
				/// rotamers but has not yet gone on to add rotamers for serine.
				if ( orig_rot->aa() >= curr_rot->aa() ) {
					//TR << "Increment current 2" << std::endl;
					++curr_id;
				} else {
					//TR << "Increment original 2" << std::endl;
					++orig_id;
					orig_prev_restype_match = false;
				}
			}
		} // end while
		correspondences_found += ii_correspondences_found;
		curr_rotamers_with_correspondence_[ ii ].reserve( ii_correspondences_found );
		curr_rotamers_without_correspondence_[ ii ].reserve( ii_curr_nrots - ii_correspondences_found );
		for ( Size jj = 1; jj <= ii_curr_nrots; ++jj ) {
			if ( curr_rot_2_orig_rot_[ ii ][ jj ] != 0 ) {
				curr_rotamers_with_correspondence_[ ii ].push_back( jj );
			} else {
				curr_rotamers_without_correspondence_[ ii ].push_back( jj );
			}
		}
	}

	TR << "Found correspondence for " << correspondences_found;
	TR << " rotamers out of " << current_rotamer_sets_->nrotamers() << std::endl;
}

void
GreenPacker::initialize_internal_correspondence_data(
	core::pose::Pose & pose
)
{
	orig_rot_2_curr_rot_.clear(); orig_rot_2_curr_rot_.resize( pose.total_residue() );
	curr_rot_2_orig_rot_.clear(); curr_rot_2_orig_rot_.resize( pose.total_residue() );
	curr_rotamers_with_correspondence_.clear();
	curr_rotamers_without_correspondence_.clear();
	curr_rotamers_with_correspondence_.resize( pose.total_residue() );
	curr_rotamers_without_correspondence_.resize( pose.total_residue() );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		Size const ii_moltenresid = current_rotamer_sets_->resid_2_moltenres( ii );
		if ( ii_moltenresid != 0 ) {
			orig_rot_2_curr_rot_[ ii ].resize( original_rotamers_[ ii ].size() );
			std::fill( orig_rot_2_curr_rot_[ ii ].begin(), orig_rot_2_curr_rot_[ ii ].end(), 0 );
			curr_rot_2_orig_rot_[ ii ].resize( current_rotamer_sets_->rotamer_set_for_moltenresidue( ii_moltenresid )->num_rotamers() );
			std::fill( curr_rot_2_orig_rot_[ ii ].begin(), curr_rot_2_orig_rot_[ ii ].end(), 0 );
		}
	}
}

void
GreenPacker::compute_energies(
	core::pose::Pose & pose
)
{
	using namespace core::pack::interaction_graph;
	/// 0. Create an interaction graph
	InteractionGraphBaseOP ig = InteractionGraphFactory::create_interaction_graph(
		*current_task_,
		*current_rotamer_sets_,
		pose,
		*full_sfxn_ );

	PrecomputedPairEnergiesInteractionGraphOP pig(
		dynamic_cast< PrecomputedPairEnergiesInteractionGraph * > (
		ig.get() ));

	/// if the dynamic cast failed, then the packer task has produced
	/// an on the fly interaction graph (or some other non-precomputed IG
	/// that hadn't been invented by 4/18/2008) which makes the GreenPacker
	/// pointless... just because it's pointless doesn't mean it should fail,
	/// though, so this code should be revised in the future.
	if ( ! pig ) {
		utility_exit_with_message( "GreenPacker asked to use on-the-fly interaction graph.  Why?" );
	}
	pig->initialize( *current_rotamer_sets_ );


	/// 1. Compute all pair energies across the interface, and declare those edges final.
	/// Let the rotamer sets class do all the hard work here...
	current_rotamer_sets_->precompute_two_body_energies(
		pose, *full_sfxn_,
		current_inter_group_packer_neighbor_graph_, pig
	);

	/// 2. Compute batch CD pair energies for intra-group edges
	/// Let the rotamer sets class do all the hard work here...
	current_rotamer_sets_->precompute_two_body_energies(
		pose, *cd_sfxn_,
		current_intra_group_packer_neighbor_graph_, pig
	);

	/// 3. Add in CI pair energies for intra-group edges where both rotamers have a correspondence
	add_precomputed_energies( pose, pig );

	/// 4. Compute CI pair energies for those rotamers lacking a correspondence.
	compute_absent_energies( pose, pig );

	/// 5. Compute one-body energies
	current_rotamer_sets_->compute_one_body_energies(
		pose, *full_sfxn_, current_packer_neighbor_graph_, pig );

	current_ig_ = pig;
}

void
GreenPacker::add_precomputed_energies(
	core::pose::Pose & pose,
	core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraphOP pig
)
{
	using namespace core::pack;
	using namespace core::pack::interaction_graph;
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		Size const ii_moltenres_curr = current_rotamer_sets_->resid_2_moltenres( ii );
		Size const ii_moltenres_orig = reference_resid_2_moltenres_[ ii ];

		if ( ii_moltenres_curr == 0 || ii_moltenres_orig == 0 ) continue;

		for ( ci_rpes_->reset_edge_list_iterator_for_node( ii_moltenres_orig );
				! ci_rpes_->edge_list_iterator_at_end();
				ci_rpes_->increment_edge_list_iterator() ) {
			EdgeBase const & edge( ci_rpes_->get_edge() );

			Size const jj_moltenres_orig = edge.get_other_ind( ii_moltenres_orig );
			if ( jj_moltenres_orig < ii_moltenres_orig ) continue; // only deal with upper edges

			Size const jj = reference_moltenres_2_resid_[ jj_moltenres_orig ];
			Size const jj_moltenres_curr = current_rotamer_sets_->resid_2_moltenres( jj );
			if ( jj_moltenres_curr == 0 ) continue;

			assert( dynamic_cast< PrecomputedPairEnergiesEdge const * > ( & edge ) );

			PrecomputedPairEnergiesEdge const & precomp_edge(
				static_cast< PrecomputedPairEnergiesEdge const & > ( edge ) );

			if ( !pig->get_edge_exists( ii_moltenres_curr, jj_moltenres_curr )) {
				pig->add_edge( ii_moltenres_curr, jj_moltenres_curr );
			}

			/// iterate across rotamer pairs with a correspondence.

			/// rename a few variables for brevity
			utility::vector1< Size > const & ii_corr_rots( curr_rotamers_with_correspondence_[ ii ]);
			utility::vector1< Size > const & jj_corr_rots( curr_rotamers_with_correspondence_[ jj ]);
			utility::vector1< Size > const & ii_curr_2_orig( curr_rot_2_orig_rot_[ ii ] );
			utility::vector1< Size > const & jj_curr_2_orig( curr_rot_2_orig_rot_[ jj ] );
			Size ii_ncorr_rots = ii_corr_rots.size();
			Size jj_ncorr_rots = jj_corr_rots.size();

			if ( ii_ncorr_rots == 0 || jj_ncorr_rots == 0 ) continue;

			for ( Size kk = 1; kk <= ii_ncorr_rots; ++kk ) {
				Size const kk_curr_rot = ii_corr_rots[ kk ];
				Size const kk_orig_rot = ii_curr_2_orig[ kk_curr_rot ];

				for ( Size ll = 1; ll <= jj_corr_rots.size(); ++ll ) {
					Size const ll_curr_rot = jj_corr_rots[ ll ];
					Size const ll_orig_rot = jj_curr_2_orig[ ll_curr_rot ];
					core::PackerEnergy const kkll_energy = precomp_edge.get_two_body_energy( kk_orig_rot, ll_orig_rot );
					pig->add_to_two_body_energies_for_edge(
						ii_moltenres_curr, jj_moltenres_curr,
						kk_curr_rot, ll_curr_rot,
						kkll_energy );
				}
			}
		}
	}
}

/// @details this time, iterate across edges in the packer neighbor graph
/// and compute pair interaction energies between rotamers that do not
/// correspond to any of the original rotamers with all other rotamers
/// in the neighboring set.  This has to be careful not to double-count
/// interactions between rotamer pairs that both lack a correspondence.
void
GreenPacker::compute_absent_energies(
	core::pose::Pose & pose,
	PrecomputedPairEnergiesInteractionGraphOP pig
)
{
	using namespace core::scoring;

	/// 1.  Absent short ranged context independent two body energies
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		Size const ii_moltenres_curr = current_rotamer_sets_->resid_2_moltenres( ii );
		if ( ii_moltenres_curr == 0 ) continue;

		Size const ii_nnew_rots = curr_rotamers_without_correspondence_[ ii ].size();
		for ( Graph::EdgeListIter
				iru  = current_intra_group_packer_neighbor_graph_->get_node(ii)->upper_edge_list_begin(),
				irue = current_intra_group_packer_neighbor_graph_->get_node(ii)->upper_edge_list_end();
				iru != irue; ++iru ) {

			Size const jj = (*iru)->get_second_node_ind();
			Size const jj_moltenres_curr = current_rotamer_sets_->resid_2_moltenres( jj );
			if ( jj_moltenres_curr == 0 ) continue;

			Size const jj_nnew_rots = curr_rotamers_without_correspondence_[ jj ].size();

			if ( ii_nnew_rots == 0 && jj_nnew_rots == 0 ) continue; // AWESOME, no work to be done.

			if ( ! pig->get_edge_exists( ii_moltenres_curr, jj_moltenres_curr ) ) {
				pig->add_edge( ii_moltenres_curr, jj_moltenres_curr );
			}
			compute_absent_srci_energies_for_residue_pair( pose, pig, ii, jj );
		}
	}

	/// 2.  Absent long ranged context independent two body energies
	// Iterate across the long range energy functions and use the iterators generated
	// by the LRnergy container object
	for ( ScoreFunction::LR_2B_MethodIterator
			lr_iter = ci_sfxn_->long_range_energies_begin(),
			lr_end  = ci_sfxn_->long_range_energies_end();
			lr_iter != lr_end; ++lr_iter ) {
		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*lr_iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.
		// Potentially O(N^2) operation...

		for ( Size ii = 1; ii <= pose.total_residue(); ++ ii ) {
			Size const ii_moltenres_curr = current_rotamer_sets_->resid_2_moltenres( ii );
			if ( ii_moltenres_curr == 0 ) continue;
			Size const ii_nnew_rots = curr_rotamers_without_correspondence_[ ii ].size();

			for ( ResidueNeighborConstIteratorOP
						rni = lrec->const_upper_neighbor_iterator_begin( ii ),
						rniend = lrec->const_upper_neighbor_iterator_end( ii );
						(*rni) != (*rniend); ++(*rni) ) {
				Size const jj = rni->upper_neighbor_id();

				Size const jj_moltenres_curr = current_rotamer_sets_->resid_2_moltenres( jj );
				if ( jj_moltenres_curr == 0 ) continue;
				Size const jj_nnew_rots = curr_rotamers_without_correspondence_[ jj ].size();

				if ( ii_nnew_rots == 0 && jj_nnew_rots == 0 ) continue; // AWESOME, no work to be done.

				if ( ! pig->get_edge_exists( ii_moltenres_curr, jj_moltenres_curr ) ) {
					pig->add_edge( ii_moltenres_curr, jj_moltenres_curr );
				}
				compute_absent_lrci_energies_for_residue_pair( pose, **lr_iter, pig, ii, jj );

			}
		}
	}

}

void
GreenPacker::compute_absent_srci_energies_for_residue_pair(
	core::pose::Pose & pose,
	PrecomputedPairEnergiesInteractionGraphOP pig,
	Size lower_res,
	Size upper_res
)
{
	using namespace core::conformation;
	using namespace core::scoring;

	Size const lower_res_moltenresid = current_rotamer_sets_->resid_2_moltenres( lower_res );
	assert( lower_res_moltenresid );

	Size const upper_res_moltenresid = current_rotamer_sets_->resid_2_moltenres( upper_res );
	assert( upper_res_moltenresid );

	Size const lower_res_nrots = current_rotamer_sets_->nrotamers_for_moltenres( lower_res_moltenresid );
	Size const lower_res_nnew_rots = curr_rotamers_without_correspondence_[ lower_res ].size();

	Size const upper_res_nrots = current_rotamer_sets_->nrotamers_for_moltenres( upper_res_moltenresid );
	Size const upper_res_nnew_rots = curr_rotamers_without_correspondence_[ upper_res ].size();

	if ( lower_res_nnew_rots > 0 ) {
		for ( Size ii = 1; ii <= lower_res_nnew_rots; ++ii ) {
			Size ii_rot_index = curr_rotamers_without_correspondence_[ lower_res ][ ii ];
			ResidueCOP ii_rot = current_rotamer_sets_->rotamer_set_for_residue( lower_res )->rotamer( ii_rot_index );
			for ( Size jj = 1; jj <= upper_res_nrots; ++jj ) {
				EnergyMap emap;
				ResidueCOP jj_rot = current_rotamer_sets_->rotamer_set_for_residue( upper_res )->rotamer( jj );
				ci_sfxn_->eval_ci_2b( *ii_rot, *jj_rot, pose, emap );
				Real const weighted_energy = ci_sfxn_->weights().dot( emap );
				pig->add_to_two_body_energies_for_edge(
					lower_res_moltenresid, upper_res_moltenresid,
					ii_rot_index, jj, weighted_energy );
			}
		}
	}

	if ( upper_res_nnew_rots > 0 ) {
		for ( Size ii = 1; ii <= lower_res_nrots; ++ii ) {

			// Avoid double counting pair interactions by skipping the new rotamers of the lower residue,
			// since their interaction energies with the upper residue's rotamers have already been calculated
			if ( curr_rot_2_orig_rot_[ lower_res ][ ii ] == 0 ) continue;

			ResidueCOP ii_rot = current_rotamer_sets_->rotamer_set_for_residue( lower_res )->rotamer( ii );
			for ( Size jj = 1; jj <= upper_res_nnew_rots; ++jj ) {
				Size const jj_rot_index = curr_rotamers_without_correspondence_[ upper_res ][ jj ];
				EnergyMap emap;
				ResidueCOP jj_rot = current_rotamer_sets_->rotamer_set_for_residue( upper_res )->rotamer( jj_rot_index );
				ci_sfxn_->eval_ci_2b( *ii_rot, *jj_rot, pose, emap );
				Real const weighted_energy = ci_sfxn_->weights().dot( emap );
				pig->add_to_two_body_energies_for_edge(
					lower_res_moltenresid, upper_res_moltenresid,
					ii, jj_rot_index, weighted_energy );
			}
		}
	}

}

void
GreenPacker::compute_absent_lrci_energies_for_residue_pair(
	core::pose::Pose & pose,
	core::scoring::methods::LongRangeTwoBodyEnergy const & lre,
	PrecomputedPairEnergiesInteractionGraphOP pig,
	Size lower_res,
	Size upper_res
)
{
	using namespace core::conformation;
	using namespace core::scoring;

	Size const lower_res_moltenresid = current_rotamer_sets_->resid_2_moltenres( lower_res );
	assert( lower_res_moltenresid );

	Size const upper_res_moltenresid = current_rotamer_sets_->resid_2_moltenres( upper_res );
	assert( upper_res_moltenresid );

	Size const lower_res_nrots = current_rotamer_sets_->nrotamers_for_moltenres( lower_res_moltenresid );
	Size const lower_res_nnew_rots = curr_rotamers_without_correspondence_[ lower_res ].size();

	Size const upper_res_nrots = current_rotamer_sets_->nrotamers_for_moltenres( upper_res_moltenresid );
	Size const upper_res_nnew_rots = curr_rotamers_without_correspondence_[ upper_res ].size();

	if ( lower_res_nnew_rots > 0 ) {
		for ( Size ii = 1; ii <= lower_res_nnew_rots; ++ii ) {
			Size ii_rot_index = curr_rotamers_without_correspondence_[ lower_res ][ ii ];
			ResidueCOP ii_rot = current_rotamer_sets_->rotamer_set_for_residue( lower_res )->rotamer( ii_rot_index );
			for ( Size jj = 1; jj <= upper_res_nrots; ++jj ) {
				EnergyMap emap;
				ResidueCOP jj_rot = current_rotamer_sets_->rotamer_set_for_residue( upper_res )->rotamer( jj );
				lre.residue_pair_energy( *ii_rot, *jj_rot, pose, *ci_sfxn_, emap );
				Real const weighted_energy = ci_sfxn_->weights().dot( emap );
				pig->add_to_two_body_energies_for_edge(
					lower_res_moltenresid, upper_res_moltenresid,
					ii_rot_index, jj, weighted_energy );
			}
		}
	}

	if ( upper_res_nnew_rots > 0 ) {
		for ( Size ii = 1; ii <= lower_res_nrots; ++ii ) {

			// Avoid double counting pair interactions by skipping the new rotamers of the lower residue,
			// since their interaction energies with the upper residue's rotamers have already been calculated
			if ( curr_rot_2_orig_rot_[ lower_res ][ ii ] == 0 ) continue;

			ResidueCOP ii_rot = current_rotamer_sets_->rotamer_set_for_residue( lower_res )->rotamer( ii );
			for ( Size jj = 1; jj <= upper_res_nnew_rots; ++jj ) {
				Size const jj_rot_index = curr_rotamers_without_correspondence_[ upper_res ][ jj ];
				EnergyMap emap;
				ResidueCOP jj_rot = current_rotamer_sets_->rotamer_set_for_residue( upper_res )->rotamer( jj_rot_index );
				lre.residue_pair_energy( *ii_rot, *jj_rot, pose, *ci_sfxn_, emap );
				Real const weighted_energy = ci_sfxn_->weights().dot( emap );
				pig->add_to_two_body_energies_for_edge(
					lower_res_moltenresid, upper_res_moltenresid,
					ii, jj_rot_index, weighted_energy );
			}
		}
	}
}


void
GreenPacker::run_sa(
	core::pose::Pose & pose
)
{
	core::pack::pack_rotamers_run( pose, current_task_, current_rotamer_sets_, current_ig_ );
}

/// @details Free memory that is no longer needed
void
GreenPacker::cleanup()
{
	current_task_ = 0;
	current_rotamer_sets_ = 0;
	current_rotamers_.clear();
	current_ig_ = 0;

	current_packer_neighbor_graph_ = 0;
	current_inter_group_packer_neighbor_graph_ = 0;
	current_intra_group_packer_neighbor_graph_ = 0;

	orig_rot_2_curr_rot_.clear();
	curr_rot_2_orig_rot_.clear();
	curr_rotamers_with_correspondence_.clear();
	curr_rotamers_without_correspondence_.clear();

}

void GreenPacker::store_reference_pose_geometry( Pose & pose )
{
	orig_bb_tors_.resize( pose.total_residue() );
	original_bb_rep_coords_.resize( pose.total_residue() );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		orig_bb_tors_[ ii ] = pose.residue( ii ).mainchain_torsions();
		original_bb_rep_coords_[ ii ] = pose.residue( ii ).xyz( 1 );
	}
}

void GreenPacker::compare_input_pose_geometry_to_reference( Pose & pose )
{
	assert( pose.total_residue() == orig_bb_tors_.size());
	for ( Size ii = 1; ii <= orig_bb_tors_.size(); ++ii ) {
		for ( Size jj = 1; jj <= orig_bb_tors_[ ii ].size(); ++jj ) {
			if ( std::abs( basic::periodic_range(
				orig_bb_tors_[ ii ][ jj ] - pose.residue( ii ).mainchain_torsions()[ jj ],
				180 ) ) > 1e-8 ) {
				std::cerr << "Critical Error in GreenPacker -- Backbone torsions have changed since original packing" << std::endl;
				std::cerr << "Residue " << ii << " torsion " << jj << " originally: " << orig_bb_tors_[ ii ][ jj ];
				std::cerr << " currently: " << pose.residue( ii ).mainchain_torsions()[ jj ] << std::endl;
				utility_exit_with_message("Bad torsion in GreenPacker" );
			}
		}
	}

}



}
}
