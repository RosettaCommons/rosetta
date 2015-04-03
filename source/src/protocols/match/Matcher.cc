// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/Matcher.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/Matcher.hh>

// Package headers
#include <protocols/match/Hit.hh>
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/MatcherTask.hh>
#include <protocols/match/MatchSet.hh>
#include <protocols/match/OccupiedSpaceHash.hh>

#include <protocols/match/downstream/ActiveSiteGrid.hh>
#include <protocols/match/downstream/ClassicMatchAlgorithm.hh>
#include <protocols/match/downstream/DownstreamAlgorithm.hh>
#include <protocols/match/downstream/DownstreamBuilder.hh>
#include <protocols/match/downstream/LigandConformerBuilder.hh>
#include <protocols/match/downstream/RigidLigandBuilder.hh>
#include <protocols/match/downstream/SecondaryMatcherToDownstreamResidue.hh>
#include <protocols/match/downstream/SecondaryMatcherToUpstreamResidue.hh>
#include <protocols/match/downstream/SecMatchEvaluatorFactory.hh>

#include <protocols/match/output/MatchProcessor.hh>
#include <protocols/match/output/UpstreamDownstreamCollisionFilter.hh>

#include <protocols/match/upstream/ProteinSCSampler.hh>
#include <protocols/match/upstream/ProteinUpstreamBuilder.hh>
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>
#include <protocols/match/upstream/UpstreamBuilder.hh>


// Project headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
//#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh>
#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>

#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/SemiRotamericSingleResidueDunbrackLibrary.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairGeneric.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Utility headers
#include <utility/LexicographicalIterator.hh>
#include <utility/OrderedTuple.hh>
#include <utility/string_util.hh>

// C++ headers
#include <map>
#include <string>

#include <utility/vector1.hh>

//Auto Headers
#include <protocols/match/downstream/SecMatchResiduePairEvaluator.hh>

namespace protocols {
namespace match {

static thread_local basic::Tracer TR( "protocols.match.Matcher" );

/// Construction and Destruction
Matcher::Matcher() :
	same_build_resids_for_all_csts_( true ),
	n_geometric_constraints_( 0 ),
	read_gridlig_file_( false ),
	use_input_sc_( false ),
	dynamic_grid_refinement_( false ),
	output_matches_as_singular_downstream_positioning_( false ),
	check_potential_dsbuilder_incompatibility_( false ),
	build_round1_hits_twice_( false )
{
	relevant_downstream_atoms_.clear();
}

Matcher::~Matcher() {}

void Matcher::set_upstream_pose( core::pose::Pose const & pose )
{
	upstream_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );
}

void Matcher::set_downstream_pose(
	core::pose::Pose const & pose,
	utility::vector1< core::id::AtomID > orientation_atoms
)
{
	runtime_assert( orientation_atoms.size() == 3 );

	downstream_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );
	downstream_orientation_atoms_ = orientation_atoms;
}

void
Matcher::set_original_scaffold_build_points( utility::vector1< Size > const & resids )
{
	same_build_resids_for_all_csts_ = true;
	pose_build_resids_ = resids;
	std::sort( pose_build_resids_.begin(), pose_build_resids_.end() ); // keep sorted
	per_cst_build_resids_.resize( 0 );
}

void Matcher::set_original_scaffold_build_points_for_constraint(
	Size cst_id,
	utility::vector1< Size > const & resids
)
{
	same_build_resids_for_all_csts_ = false;
	runtime_assert( n_geometric_constraints_ > 0 ); // n_geometric_constraints_ must be set first
	if ( per_cst_build_resids_.size() == 0 ) {
		per_cst_build_resids_.resize( n_geometric_constraints_ );
		pose_build_resids_.clear();
	}
	per_cst_build_resids_[ cst_id ] = resids;
	std::sort( per_cst_build_resids_[ cst_id ].begin(), per_cst_build_resids_[ cst_id ].end() );

	/// keep track of the individual resids that are built from by any geometric constraint
	std::list< Size > uniq_resids;
	for ( Size ii = 1; ii <= pose_build_resids_.size(); ++ii ) {
		uniq_resids.push_back( pose_build_resids_[ ii ] );
	}
	for ( Size ii = 1; ii <= per_cst_build_resids_[ cst_id ].size(); ++ii ) {
		uniq_resids.push_back( per_cst_build_resids_[ cst_id ][ ii ] );
	}

	uniq_resids.sort();
	uniq_resids.unique();

	pose_build_resids_.resize( uniq_resids.size() );
	std::copy( uniq_resids.begin(), uniq_resids.end(), pose_build_resids_.begin() );

}


void Matcher::set_n_geometric_constraints( Size n_constraints )
{
	n_geometric_constraints_ = n_constraints;

	hits_.resize( n_geometric_constraints_ );
	upstream_builders_.resize( n_geometric_constraints_ );
	build_set_id_for_restype_.resize( n_geometric_constraints_ );
	representative_downstream_algorithm_.resize( n_geometric_constraints_, 0 );
	//std::fill( representative_downstream_algorithm_.begin(), representative_downstream_algorithm_.end(), 0 );
	downstream_algorithms_.resize( n_geometric_constraints_ );
	geomcst_is_upstream_only_.resize( n_geometric_constraints_, false );
	downstream_builders_.resize( n_geometric_constraints_ );
	//std::fill( downstream_builders_.begin(), downstream_builders_.end(), 0 );
	per_constraint_build_points_.resize( n_geometric_constraints_ );
	geom_cst_has_primary_modification_.resize( n_geometric_constraints_ );
	std::fill( geom_cst_has_primary_modification_.begin(), geom_cst_has_primary_modification_.end(), false );
}

void Matcher::add_upstream_restype_for_constraint(
	Size cst_id,
	core::chemical::ResidueTypeCOP restype
)
{
	/// ASSUMPTION: matching from protein sidechain
	assert( restype->aa() <= core::chemical::num_canonical_aas );
	assert( build_set_id_for_restype_[ cst_id ].find( restype->name() ) == build_set_id_for_restype_[ cst_id ].end() );

	if ( ! upstream_builders_[ cst_id ] ) {
		upstream::ProteinUpstreamBuilderOP prot_sc_builder( new upstream::ProteinUpstreamBuilder );
		/// default to dunbrack sampler
		prot_sc_builder->set_sampler( upstream::ProteinSCSamplerCOP( upstream::ProteinSCSamplerOP( new upstream::DunbrackSCSampler ) ) );
		prot_sc_builder->set_use_input_sidechain( use_input_sc_ );
		upstream_builders_[ cst_id ] = prot_sc_builder;
	}

	upstream::BuildSet build_set;
	build_set.set_residue_type( restype, restype->aa() == core::chemical::aa_gly ); // HACK gly means backbone only

	runtime_assert( dynamic_cast< upstream::ProteinUpstreamBuilder * > ( upstream_builders_[ cst_id ].get() ) );

	upstream::ProteinUpstreamBuilderOP prot_sc_builder( utility::pointer::dynamic_pointer_cast< upstream::ProteinUpstreamBuilder > ( upstream_builders_[ cst_id ] ));
	prot_sc_builder->add_build_set( build_set );
	build_set_id_for_restype_[ cst_id ][ restype->name() ] = prot_sc_builder->n_build_sets();

}

void Matcher::desymmeterize_upstream_restype_for_constraint(
	Size cst_id
)
{
	upstream::ProteinUpstreamBuilderOP prot_sc_builder = utility::pointer::dynamic_pointer_cast< upstream::ProteinUpstreamBuilder > ( upstream_builders_[ cst_id ] );
	if ( ! prot_sc_builder ) {
		utility_exit_with_message( "Could not desymmeterize the upstream builder for geometric constraint " + utility::to_string( cst_id ) + ".  Desymmeterization is only available for protein residues." );
	}
	/// Replace the existing sampler with one that is desymmeterized.  This will over-write any other stored data in the DunbrackSCSampler if any
	/// other data gets added to this class.  If new data should be added, then the ProteinUpstreamBuilder needs to be modified to hand out non-const
	/// access to its sampler.
	upstream::DunbrackSCSamplerOP sampler( new upstream::DunbrackSCSampler );
	sampler->set_desymmeterize( true );
	prot_sc_builder->set_sampler( sampler );
}

void Matcher::set_sample_startegy_for_constraint(
	Size cst_id,
	core::chemical::ResidueTypeCOP restype,
	Size chi,
	upstream::SampleStrategyData const & strat
)
{
	runtime_assert( build_set_id_for_restype_[ cst_id ].find( restype->name() )
		!= build_set_id_for_restype_[ cst_id ].end() );

	//Size build_set_id = build_set_id_for_restype_[ cst_id ][ restype->name() ];
	runtime_assert( dynamic_cast< upstream::ProteinUpstreamBuilder * > ( upstream_builders_[ cst_id ].get() ) );
	upstream::ProteinUpstreamBuilderOP prot_sc_builder( utility::pointer::dynamic_pointer_cast< upstream::ProteinUpstreamBuilder > ( upstream_builders_[ cst_id ] ));

	upstream::BuildSet & build_set = prot_sc_builder->build_set( restype );

	build_set.set_sample_strategy_for_chi( chi, strat );
}

void
Matcher::set_fa_dun_cutoff_for_constraint(
	Size cst_id,
	core::chemical::ResidueTypeCOP restype,
	core::Real fa_dun_cutoff
)
{
	assert( build_set_id_for_restype_[ cst_id ].find( restype->name() )
		!= build_set_id_for_restype_[ cst_id ].end() );

	assert( dynamic_cast< upstream::ProteinUpstreamBuilder * > ( upstream_builders_[ cst_id ].get() ) );

	upstream::ProteinUpstreamBuilderOP prot_sc_builder( utility::pointer::dynamic_pointer_cast< upstream::ProteinUpstreamBuilder > ( upstream_builders_[ cst_id ] ));
	upstream::BuildSet & build_set = prot_sc_builder->build_set( restype );
	build_set.set_fa_dun_cutoff( fa_dun_cutoff );
}


void Matcher::add_external_geometry_samples_for_constraint(
	Size cst_id,
	core::chemical::ResidueTypeCOP restype,
	utility::vector1< std::string >  const & upstream_launch_atoms,
	utility::vector1< core::id::AtomID > const & downstream_3atoms,
	toolbox::match_enzdes_util::ExternalGeomSampler const & exgeom,
	Size const exgeom_id,
	bool enumerate_ligand_rotamers /* = false */,
	bool catalytic_bond /*= false */,
	bool build_round1_hits_twice /* = false */
)
{
	TR << "     Adding Classical Match Algorithm with geometry samples: " << std::endl;
	TR << "     tor_U3D1:";
	for ( Size ii = 1; ii <= exgeom.n_tor_U3D1_samples(); ++ii ) {
		TR << " " << exgeom.tor_U3D1_samples()[ ii ];
	}
	TR << std::endl;

	TR << "     ang_U2D1:";
	for ( Size ii = 1; ii <= exgeom.n_ang_U2D1_samples(); ++ii ) {
		TR << " " << exgeom.ang_U2D1_samples()[ ii ];
	}
	TR << std::endl;

	TR << "     dis_U1D1:";
	for ( Size ii = 1; ii <= exgeom.n_dis_U1D1_samples(); ++ii ) {
		TR << " " << exgeom.dis_U1D1_samples()[ ii ];
	}
	TR << std::endl;

	TR << "     tor_U2D2:";
	for ( Size ii = 1; ii <= exgeom.n_tor_U2D2_samples(); ++ii ) {
		TR << " " << exgeom.tor_U2D2_samples()[ ii ];
	}
	TR << std::endl;

	TR << "     ang_U1D2:";
	for ( Size ii = 1; ii <= exgeom.n_ang_U1D2_samples(); ++ii ) {
		TR << " " << exgeom.ang_U1D2_samples()[ ii ];
	}
	TR << std::endl;

	TR << "     tor_U1D3:";
	for ( Size ii = 1; ii <= exgeom.n_tor_U1D3_samples(); ++ii ) {
		TR << " " << exgeom.tor_U1D3_samples()[ ii ];
	}
	TR << std::endl;


	runtime_assert( upstream_launch_atoms.size() == 3 );
	runtime_assert( downstream_3atoms.size() == 3 );

	runtime_assert( build_set_id_for_restype_[ cst_id ].find( restype->name() )
		!= build_set_id_for_restype_[ cst_id ].end() );

	downstream::DownstreamBuilderOP ds_builder = create_ds_builder(
		cst_id, restype, upstream_launch_atoms,
		downstream_3atoms, enumerate_ligand_rotamers, catalytic_bond );

	//Size build_set_id = build_set_id_for_restype_[ cst_id ][ restype->name() ];
	runtime_assert( dynamic_cast< upstream::ProteinUpstreamBuilder * > (
		upstream_builders_[ cst_id ].get() ) );
	upstream::ProteinUpstreamBuilderOP prot_sc_builder(
		utility::pointer::static_pointer_cast< upstream::ProteinUpstreamBuilder > ( upstream_builders_[ cst_id ] ));

	upstream::BuildSet & build_set = prot_sc_builder->build_set( restype );

	if ( ! build_set.has_algorithm() ) {
		downstream::ClassicMatchAlgorithmOP match_algorithm( new downstream::ClassicMatchAlgorithm( cst_id ) );
		match_algorithm->set_residue_type( restype );
		build_set.set_downstream_algorithm( match_algorithm );

		downstream_algorithms_[ cst_id ].push_back( match_algorithm );
		representative_downstream_algorithm_[ cst_id ] = match_algorithm;
		all_downstream_algorithms_.push_back( match_algorithm );
		if ( cst_id == 1 && build_round1_hits_twice ) match_algorithm->set_build_round1_hits_twice();
	}

	runtime_assert( dynamic_cast< downstream::ClassicMatchAlgorithm * > ( & build_set.algorithm() ) );
	downstream::ClassicMatchAlgorithm & algorithm( static_cast< downstream::ClassicMatchAlgorithm & > (build_set.algorithm() ) );

	algorithm.add_external_geom_sampler(
		exgeom,
		exgeom_id,
		upstream_launch_atoms[ 1 ],
		upstream_launch_atoms[ 2 ],
		upstream_launch_atoms[ 3 ],
		ds_builder
	);

}

/// Initialize a secondary matcher object based on the
/// geometry from the upstream to the downstream residue types.
void Matcher::add_secondary_upstream_match_geometry_for_constraint(
	Size geom_cst_id,
	Size target_geom_cst_id,
	core::chemical::ResidueTypeCOP candidate_restype,
	core::chemical::ResidueTypeCOP target_restype,
	utility::vector1< Size > const & candidate_atids,
	utility::vector1< Size > const & target_atids,
	toolbox::match_enzdes_util::MatchConstraintFileInfoCOP mcfi,
	std::string sec_match_str,
 	core::pose::Pose const & upstream_pose
)
{
	using namespace downstream;

	runtime_assert( candidate_atids.size() == 3 );
	runtime_assert( target_atids.size() == 3 );

	runtime_assert( build_set_id_for_restype_[ geom_cst_id ].find( candidate_restype->name() )
		!= build_set_id_for_restype_[ geom_cst_id ].end() );

	runtime_assert( dynamic_cast< upstream::ProteinUpstreamBuilder * > (
		upstream_builders_[ geom_cst_id ].get() ) );
	upstream::ProteinUpstreamBuilderOP prot_sc_builder(
		utility::pointer::static_pointer_cast< upstream::ProteinUpstreamBuilder > ( upstream_builders_[ geom_cst_id ] ));

	upstream::BuildSet & build_set = prot_sc_builder->build_set( candidate_restype );

	if ( ! build_set.has_algorithm() ) {
		SecondaryMatcherToUpstreamResidueOP secondary_match_algorithm( new SecondaryMatcherToUpstreamResidue( geom_cst_id ) );
		build_set.set_downstream_algorithm( secondary_match_algorithm );
		secondary_match_algorithm->set_target_geomcst_id( target_geom_cst_id );
		downstream_algorithms_[ geom_cst_id ].push_back( secondary_match_algorithm );
		representative_downstream_algorithm_[ geom_cst_id ] = secondary_match_algorithm;
		all_downstream_algorithms_.push_back( secondary_match_algorithm );
	}

	runtime_assert( dynamic_cast< downstream::SecondaryMatcherToUpstreamResidue * > ( & build_set.algorithm() ) );
	downstream::SecondaryMatcherToUpstreamResidue & algorithm( static_cast< downstream::SecondaryMatcherToUpstreamResidue & > (build_set.algorithm() ) );

	algorithm.add_target_restype( target_restype );

	//Old code: We replaced with SecMatchResiduePairEvaulatorOP(SRPE) codes.
	//SRPE is parents of ScoringMatchRPE and GeometrySecMatchRPE.
	/*
	GeometrySecMatchRPEOP geom_evaluator = new GeometrySecMatchRPE( *mcfi, target_atids, candidate_atids );
	for ( Size ii = 1; ii <= geom_evaluator->atom_geom_rpes().size(); ++ii ) {
		TR << "    Upstream 2ndary Match: " << geom_evaluator->atom_geom_rpes()[ ii ]->print( candidate_restype, target_restype ) << std::endl;
	}
	algorithm.add_evaluator_for_target_restype( target_restype, geom_evaluator, mcfi->index() );
	*/

	//Author:Kui Chan 101409
	//Description: added score term evaluator and combine with geometrySecMatchRPE
	SecMatchResiduePairEvaluatorOP secMatch_evaluator
				= SecMatchEvaluatorFactory::create_SecMatchResiduePairEvaluatorOP( *mcfi, target_atids, candidate_atids,
								sec_match_str, upstream_pose );
	algorithm.add_evaluator_for_target_restype( target_restype, secMatch_evaluator, mcfi->index() );
	//END Kui
}

void
Matcher::add_secondary_downstream_match_geometry_for_constraint(
	Size geom_cst_id,
	core::chemical::ResidueTypeCOP candidate_restype,
	core::chemical::ResidueTypeCOP downstream_restype,
	utility::vector1< Size > const & candidate_atids,
	utility::vector1< Size > const & target_atids,
	toolbox::match_enzdes_util::MatchConstraintFileInfoCOP mcfi,
	std::string sec_match_str,
 	core::pose::Pose const & upstream_pose,
  bool catalytic_bond
)
{
	using namespace downstream;

	runtime_assert( candidate_atids.size() == 3 );
	runtime_assert( target_atids.size() == 3 );

	runtime_assert( build_set_id_for_restype_[ geom_cst_id ].find( candidate_restype->name() )
		!= build_set_id_for_restype_[ geom_cst_id ].end() );

	runtime_assert( dynamic_cast< upstream::ProteinUpstreamBuilder * > (
		upstream_builders_[ geom_cst_id ].get() ) );
	upstream::ProteinUpstreamBuilderOP prot_sc_builder(
		utility::pointer::static_pointer_cast< upstream::ProteinUpstreamBuilder > ( upstream_builders_[ geom_cst_id ] ));

	upstream::BuildSet & build_set = prot_sc_builder->build_set( candidate_restype );

	if ( ! build_set.has_algorithm() ) {
		SecondaryMatcherToDownstreamResidueOP secondary_match_algorithm( new SecondaryMatcherToDownstreamResidue( upstream_pose_, geom_cst_id ) );
		build_set.set_downstream_algorithm( secondary_match_algorithm );
		secondary_match_algorithm->set_downstream_restype( downstream_restype );
		downstream_algorithms_[ geom_cst_id ].push_back( secondary_match_algorithm );
		representative_downstream_algorithm_[ geom_cst_id ] = secondary_match_algorithm;
		all_downstream_algorithms_.push_back( secondary_match_algorithm );
		if ( catalytic_bond ) {
			utility::vector1< core::Size > catalytic_atoms(4,0);
				catalytic_atoms[ 1 ] = candidate_atids[ 2 ];
				catalytic_atoms[ 2 ] = candidate_atids[ 1 ];
			 	//target = downstream
				catalytic_atoms[ 3 ] = target_atids[ 1 ];
				catalytic_atoms[ 4 ] = target_atids[ 2 ];
			secondary_match_algorithm->set_catalytic_atoms( catalytic_atoms );
		}
	}

	runtime_assert( dynamic_cast< downstream::SecondaryMatcherToDownstreamResidue * > ( & build_set.algorithm() ) );
	downstream::SecondaryMatcherToDownstreamResidue & algorithm( static_cast< downstream::SecondaryMatcherToDownstreamResidue & > (build_set.algorithm() ) );

/*
	GeometrySecMatchRPEOP geom_evaluator = new GeometrySecMatchRPE( *mcfi, target_atids, candidate_atids );
	for ( Size ii = 1; ii <= geom_evaluator->atom_geom_rpes().size(); ++ii ) {
		TR << "    Downstream 2ndary Match: " << geom_evaluator->atom_geom_rpes()[ ii ]->print( candidate_restype, downstream_restype ) << std::endl;
	}

	algorithm.add_evaluator( geom_evaluator, mcfi->index() );
*/
	//Author:Kui Chan 101409
	//Description: added score term evaluator and combine with geometrySecMatchRPE
	SecMatchResiduePairEvaluatorOP secMatch_evaluator
				= SecMatchEvaluatorFactory::create_SecMatchResiduePairEvaluatorOP( *mcfi, target_atids, candidate_atids,
							sec_match_str, upstream_pose );
	algorithm.add_evaluator( secMatch_evaluator, mcfi->index() );
	//End Kui

}


void
Matcher::set_occupied_space_bounding_box( BoundingBox const & bb )
{
	occ_space_bounding_box_ = bb;
	TR << "Set occupied space bounding box: Lower (";
	TR << bb.lower().x() << ", ";
	TR << bb.lower().y() << ", ";
	TR << bb.lower().z() << ") Upper (";
	TR << bb.upper().x() << ", ";
	TR << bb.upper().y() << ", ";
	TR << bb.upper().z() << ")" << std::endl;
}

void Matcher::set_hash_euclidean_bin_width( Real width )
{
	//std::fill( euclidean_bin_widths_.begin(), euclidean_bin_widths_.end(), width );
	euclidean_bin_widths_ = width;
}

void Matcher::set_hash_euler_bin_width( Real width )
{
	//std::fill( euler_bin_widths_.begin(), euler_bin_widths_.end(), width );
	euler_bin_widths_ = width;
}

void Matcher::set_hash_euclidean_bin_widths( Vector widths )
{
	euclidean_bin_widths_ = widths;
}

void Matcher::set_hash_euler_bin_widths( Vector widths)
{
	euler_bin_widths_ = widths;
}

void Matcher::set_bump_tolerance( Real permitted_overlap )
{
	runtime_assert( upstream_pose_.get() );

	if ( ! bb_grid_ ) {
		bb_grid_ = bump_grid_to_enclose_pose( *upstream_pose_ );
	}
	bb_grid_->set_general_overlap_tolerance( permitted_overlap );
}

void
Matcher::initialize_from_task(
	MatcherTask const & mtask
)
{
	set_upstream_pose( *mtask.upstream_pose() );
	set_downstream_pose( *mtask.downstream_pose(), mtask.downstream_orientation_atoms() );
	set_bump_tolerance( mtask.permitted_overlap() );
	set_occupied_space_bounding_box( mtask.occ_space_bounding_box() );
	set_hash_euclidean_bin_widths( mtask.euclidean_bin_widths() );
	set_hash_euler_bin_widths( mtask.euler_bin_widths() );

	/// active site definition
	read_gridlig_file_ = mtask.gridlig_active_site_definition();
	if ( read_gridlig_file_ ) {
		gridlig_fname_ = mtask.gridlig_file_name();
		upstream_resids_and_radii_defining_active_site_.clear();
	} else {
		upstream_resids_and_radii_defining_active_site_ = mtask.upstream_resids_and_radii_defining_active_site();
	}

	downstream_atoms_required_inside_active_site_ = mtask.downstream_atoms_required_inside_active_site();
	if( mtask.only_enumerate_non_match_redundant_ligand_rotamers() ){
		relevant_downstream_atoms_ = mtask.relevant_downstream_atoms();
	}

	use_input_sc_ = mtask.use_input_sc();
	dynamic_grid_refinement_ = mtask.dynamic_grid_refinement();

	initialize_from_file( * mtask.enz_input_data(), mtask );

	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		set_original_scaffold_build_points_for_constraint(
			ii, mtask.upstream_pose_build_resids_for_geometric_constraint( ii ) );

		//Kui Native 110809
		runtime_assert( dynamic_cast< upstream::ProteinUpstreamBuilder * > ( upstream_builders_[ ii ].get() ) );
		upstream::ProteinUpstreamBuilderOP prot_sc_builder(
		utility::pointer::static_pointer_cast< upstream::ProteinUpstreamBuilder > ( upstream_builders_[ ii ] ));

		toolbox::match_enzdes_util::MatchConstraintFileInfoListCOP constraint_list = (mtask.enz_input_data())->mcfi_list( ii );
		utility::vector1< core::chemical::ResidueTypeCOP > const & upres( constraint_list->upstream_restypes() );

		for ( Size jj = 1; jj <= upres.size(); ++jj ) {

			utility::vector1< toolbox::match_enzdes_util::MatchConstraintFileInfoCOP > const & jj_mcfis(
			constraint_list->mcfis_for_upstream_restype( upres[ jj ] ));

			for ( Size kk = 1; kk <= jj_mcfis.size(); ++kk ) {
				//TR << "ii:" << ii << " jj:" << jj << " kk:" << kk << " native:" <<jj_mcfis[ kk ]->native() << std::endl;
				prot_sc_builder->set_native_flag(jj_mcfis[ kk ]->native());
			}
		}
		//Kui Native 110809
	}

	/// Should we use the match_dspos1 output pathway?
	/// TEMP -- No valid MatchEvaluator yet exists for the match_dspos1 struct; do not use
	/// the match-by-single-downstream-positioning code until one comes online.
	if ( ! mtask.consolidate_matches() && mtask.define_match_by_single_downstream_positioning() ) {
		output_matches_as_singular_downstream_positioning_ = true;
		output_match_dspos1_for_geomcst_.resize( n_geometric_constraints_, false );
		for ( Size ii = 1; ii <= mtask.geom_csts_downstream_output().size(); ++ii ) {
			output_match_dspos1_for_geomcst_[ mtask.geom_csts_downstream_output()[ ii ] ] = true;
		}
	}

	if ( mtask.build_round1_hits_twice() ) build_round1_hits_twice_ = true;
}

/// @details Inside the CST::BEGIN blocks, The following ALGORITHM_INFO:: match input data
/// may be provided to give additional data to the matcher.  Rotamer building
/// instructions for particular geometric constraints may be given.
/// Here is an example of the kinds of geometric data.  Each CHI_STRATEGY line is appropriate
/// on it's own, but they are not appropriate together.  Lines beginning with "#" are comments
/// describing the meaning of each of the lines.
///
///    ALGORITHM_INFO:: match
///       # If your upstream residue includes a proton chi (e.g. SER/THR/TYR) but the geometry
///       # of the downstream partner does not depend on the location of the proton, then
///       # use the following line to avoid enumerating rotamers that differ only in their proton coordinate.
///       IGNORE_UPSTREAM_PROTON_CHI
///
///       # Secondary Matching:
///       # You can activate secondary matching, a continuous version of the classic discrete algorithm,
///       # by adding the line
///       SECONDARY_MATCH: DOWNSTREAM
///
///       #or
///       SECONDARY_MATCH: UPSTREAM_CST 2
///
///       # which instead of building its own hits and hashing them examines the hits generated in
///       # previous rounds for compatibility with the geometry that you're seeking.
///       # When performing secondary matching, it is not required that you specify all 6 degrees
///       # of freedom.  If you specify only a subset the geometry of only the subset you've specified
///       # will be examined (e.g. if you don't care about torsion_AB, don't specify it
///       # in the CST::BEGIN/CST::END block.
///       # You can perform secondary matching to the ligand (the downstream target) or to an upstream
///       # residue whose geometry was constructed by an earlier geometric constraint.  In the example
///       # above, the geometric constraint pointed to is #2.
///       # The first geometric constraint cannot use secondary matching, it must always use the discrete
///       # classic match algorithm.  (A later geometric constraint may of course perform secondary
///       # matching to the hits produced by the first geometric constraint.  At that point the first
///       # constraint is playing the same role in upstream matching as the ligand plays in downstream matching)
///       # END Secondary Matching comments.
///
///       # Below: chi sample strategies -- "1" is for chi-1
///       # These are the traditional ex?::level options ( ? == 1, 2, 3 or 4 )
///       CHI_STRATEGY:: CHI 1 EX_ONE_STDDEV
///       CHI_STRATEGY:: CHI 1 EX_ONE_HALF_STEP_STDDEV
///       CHI_STRATEGY:: CHI 1 EX_TWO_FULL_STEP_STDDEVS
///       CHI_STRATEGY:: CHI 1 EX_TWO_HALF_STEP_STDDEVS
///       CHI_STRATEGY:: CHI 1 EX_FOUR_HALF_STEP_STDDEVS
///       CHI_STRATEGY:: CHI 1 EX_THREE_THIRD_STEP_STDDEVS
///       CHI_STRATEGY:: CHI 1 EX_SIX_QUARTER_STEP_STDDEVS
///
///       # Below: The "AA" field, followed by the 3-letter AA code gives the sub-specification for
///       # a block that contains multiple amino acids.  If this block allowed both an ASN and a GLN, then
///       # the first line would apply to GLN rotamers only and the second line to ASN rotamers only.
///       # "AA" fields may be included with
///       # any of the CHI_STRATEGY options listed here.
///       CHI_STRATEGY:: AA GLN CHI 2 EX_SIX_QUARTER_STEP_STDDEVS
///       CHI_STRATEGY:: AA ASN CHI 2 EX_ONE_STDDEV
///
///       # Everything below: additional chi sample strategies unique to the matcher.
///
///       CHI_STRATEGY:: CHI 1 STEP_WITHIN_SD_RANGE STEP 3.0
///       # the above line says "step by 3 degrees within a single standard deviation
///       # to the left and right of the mean for chi 1"
///
///       CHI_STRATEGY:: CHI 1 STEP_WITHIN_SD_RANGE STEP 3.0 SD_RANGE 2.5
///       # the above line says "step by 3 degrees within 2.5 standard deviations to
///       # the left and the right of the mean for chi 1"
///
///       CHI_STRATEGY:: CHI 2 AA HIS NON_ROTAMERIC_CHI_EXPANSION N_SAMPLES 3 REQUISIT_PROBABILITY 0.10
///       # the line above says that "for histidine chi 2 (which is non-rotameric in the 2010 Dunbrak library)
///       # take 3 extra samples to the left and right of the median chi value inside the non-rotameric chi well
///       # for rotamer wells that have a probability at least 0.1".  The 2010 Dunbrack library must be active.
///       # (the -dun10 flag should be on the command line).
///
///       CHI_STRATEGY:: CHI 2 AA HIS NON_ROTAMERIC_CHI_EXPANSION N_SAMPLES 3
///       # the line above is similar for the previous command, but leaves off the REQUISIT_PROBABILITY
///       # which sets a default REQUISIT_PROBABILITY of 2 / #chi-bins.  For HIS, this is
///       # 2 / 36th, since there are 12 bins for chi 2 and 3 bins for chi1 (36 chi bins total).
///
///       #or
///       DESYMMETERIZE
///       # the line above instructs the DunbrackSCSampler to expand the sampling for the symmetric amino acids
///       # around the final (symmetric) chi angle for ASP/GLU/PHE/TYR.  For example, it would treat ASP's sampling
///       # as if there's a difference between OD1 and OD2.
///
///   ALGORITHM_INFO::END
///
void Matcher::initialize_from_file(
	toolbox::match_enzdes_util::EnzConstraintIO const & enz_data,
	MatcherTask const & mtask
)
{
	//std::cout << "APL DEBUG Matcher.cc::initialize_from_file begin" << std::endl;
	using namespace toolbox::match_enzdes_util;

	runtime_assert( upstream_pose_ != 0 );
	runtime_assert( downstream_pose_ != 0 );
	TR << "Matcher::initialize_from_file: n geometric constraints: " << enz_data.mcfi_lists_size() << std::endl;

	set_n_geometric_constraints( enz_data.mcfi_lists_size() );

	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		TR << "Begin constraint" << ii << std::endl;

		MatchConstraintFileInfoListCOP constraint_list = enz_data.mcfi_list( ii );

		utility::vector1< core::chemical::ResidueTypeCOP > const & upres( constraint_list->upstream_restypes() );

		for ( Size jj = 1; jj <= upres.size(); ++jj ) {
			add_upstream_restype_for_constraint( ii, upres[ jj ] );

			utility::vector1< MatchConstraintFileInfoCOP > const & jj_mcfis(
				constraint_list->mcfis_for_upstream_restype( upres[ jj ] ));
			TR << " Upstream residue type " << upres[ jj ]->name() << " for geometric constraint #" << ii << std::endl;

			for ( Size kk = 1; kk <= jj_mcfis.size(); ++kk ) {

				//Author: Kui Chan Date:101309
				std::string sec_match_str("");

				bool secondary_matching( false );
				bool secondary_match_upstream_residue( false );
				Size secondary_match_upstream_geomcst_id( 0 );

				utility::vector1< bool > chi_sample_data_in_file( upres[ jj ]->nchi(), false );
				//utility::vector1< SampleStrategyData > file_chi_sample_data( upres[ jj ]->nchi() );

				/// Process ALGORITHM_INFO data
				std::map< std::string, utility::vector1< std::string > > const &
					alg_info( jj_mcfis[ kk ]->algorithm_inputs() );
				if ( alg_info.find( "match" ) != alg_info.end() ) {
					utility::vector1< std::string > const & info( alg_info.find( "match" )->second );

					/// Line by line.  Currently, there are no checks for consistency across multiple lines.
					for ( Size ll = 1; ll <= info.size(); ++ll ) {
						std::string llstr = info[ ll ];
						std::istringstream llstream( llstr );
						std::string first;
						llstream >> first;
						if( first == "MAX_DUNBRACK_ENERGY" ){
							core::Real cutoff;
							llstream >> cutoff;
							TR << "Setting dunbrack energy cutoff for restype " << upres[jj]->name() << " in constraint " << ii << " to " << cutoff << "." << std::endl;
							set_fa_dun_cutoff_for_constraint( ii, upres[jj], cutoff );
						} else if ( first == "IGNORE_UPSTREAM_PROTON_CHI" ) {
							/// iterate across the proton chi, set their sample strategy to no_samples.
							upstream::SampleStrategyData nosamps; nosamps.set_strategy( upstream::no_samples );
							for ( Size mm = 1; mm <= upres[ jj ]->n_proton_chi(); ++mm ) {
								Size mmchi = upres[ jj ]->proton_chi_2_chi( mm );
								if ( chi_sample_data_in_file[ mmchi ] ) {
									TR << "  WARNING:: Already encountered chi sampling strategy for proton chi " << mmchi << " and will NOT ignore this proton chi" << std::endl;
								} else {
									set_sample_startegy_for_constraint(	ii, upres[ jj ], mmchi, nosamps );
									TR << "  ALGORITHM_INFO:: match -- Ignoring proton chi for " << upres[ jj ]->name() << " chi # " << mmchi << std::endl;
									chi_sample_data_in_file[ mmchi ] = true;
								}
							}
						} else if ( first == "CHI_STRATEGY::" ) {
							std::string chi_string;
							llstream >> chi_string;
							if ( llstream.bad() ) {
								utility_exit_with_message( "Expected 'CHI' or 'AA' following CHI_STRATEGY:: on line '" + llstr + "'" );
							}
							if ( chi_string == "AA" ) {
								std::string aa3;
								llstream >> aa3;
								//std::map< std::string, core::chemical::AA >::const_iterator iter = core::chemical::name2aa().find( aa3 );
								//if ( iter == core::chemical::name2aa().end() ) {
								//	utility_exit_with_message( "Expected amino acid 3-letter code following 'CHI_STRATEGY:: AA ' but read " + aa3 );
								//}
								//core::chemical::AA aa = *iter;
								core::chemical::AA aa = core::chemical::aa_from_name( aa3 );
								bool aa_matches_any = false;
								for ( Size mm = 1; mm <= upres.size(); ++mm ) {
									if ( upres[ mm ]->aa() == aa ) {
										aa_matches_any = true;
										break;
									}
								}
								if ( ! aa_matches_any ) {
									std::cerr << "ERROR: amino acid " << aa3 << " on line\n" << llstr << "\nis not accessible for this geometric constraint." << std::endl;
									std::cerr << "Available amino acids:";
									for ( Size mm = 1; mm <= upres.size(); ++mm ) {
										std::cerr << " " << upres[ mm ]->name();
									}
									std::cerr << std::endl;
									utility_exit_with_message( "Amino acid restriction in CHI_STRATEGY:: block is invalid" );
								}
								if ( aa != upres[ jj ]->aa() ) {
									//TR << "  Ignoring line '" << llstr << "' in processing amino acid " << upres[ jj ]->name() << std::endl;
									continue;
								}

								llstream >> chi_string;
							}

							if ( chi_string != "CHI" ) {
								utility_exit_with_message( "Expected 'CHI' following CHI_STRATEGY:: on line '" + llstr + "'" );
							}
							Size which_chi;
							llstream >> which_chi;
							if ( which_chi > upres[ jj ]->nchi() ) {
								TR <<  "WARNING: Ignoring rotamer sampling strategy data for chi # "
									<< which_chi
									 << " for residue type " << upres[ jj ]->name() << " because it only has " << upres[ jj ]->nchi() << " chi angles." <<  std::endl;
								continue;
							}
							if ( chi_sample_data_in_file[ which_chi ] ) {
								TR << "  WARNING:: Repeat chi info for chi " << which_chi << " is being ignored!" << std::endl;
								TR << "  WARNING:: Ignoring: '" << llstr << "'" << std::endl;
								continue;
							}
							chi_sample_data_in_file[ which_chi ] = true;
							if ( llstream.bad() ) {
								utility_exit_with_message( "Error parsing CHI_STRATEGY.  Expected an "
									"integer following 'CHI_STRATEGY::' on line '" + llstr );
							}
							std::string strategy;
							llstream >> strategy;
							if ( core::pack::task::is_rot_sample_name( strategy ) ) {
								core::pack::task::ExtraRotSample sample_level =  core::pack::task::rot_sample_from_name( strategy );
								upstream::SampleStrategyData stratdat; stratdat.set_strategy( upstream::rotameric_chi_mimic_EX_flags );
								stratdat.set_sample_level( sample_level );
								set_sample_startegy_for_constraint(	ii, upres[ jj ], which_chi, stratdat );

								TR << "  ALGORITHM_INFO:: match -- chi sampling strategy " << strategy << " for chi # " << which_chi << std::endl;

							} else if ( strategy == "STEP_WITHIN_SD_RANGE" ) {
								/// PARSE SD_RANGE
								Real step_size( 1.0 ), sd_range( 1.0 );
								std::string step, range;
								llstream >> step;
								if ( step != "STEP" ) {
									std::cerr << "ERROR:: Bad line in CHI_STRATEGY: \"" << llstr << "\"" << std::endl;
									utility_exit_with_message( "While parsing CHI_STRATEGY:: STEP_WITHIN_SD_RANGE, expected \"STEP\" after \"STEP_WITHIN_SD_RANGE\"" );
								}
								if ( llstream.bad() ) {
									std::cerr << "ERROR:: Bad line in CHI_STRATEGY: \"" << llstr << "\"" << std::endl;
									utility_exit_with_message( "While parsing CHI_STRATEGY:: STEP_WITHIN_SD_RANGE, unexpected EOF" );
								}
								llstream >> step_size;
								if ( llstream.bad() ) {
									std::cerr << "ERROR:: Bad line in CHI_STRATEGY: \"" << llstr << "\"" << std::endl;
									utility_exit_with_message( "While parsing CHI_STRATEGY:: STEP_WITHIN_SD_RANGE.  Could not read step size" );
								}
								if ( step_size == 0.0 ) {
									std::cerr << "ERROR:: Bad line in CHI_STRATEGY: \"" << llstr << "\"" << std::endl;
									utility_exit_with_message( "While parsing CHI_STRATEGY:: STEP_WITHIN_SD_RANGE.  Invalid step size of 0.0" );
								}
								llstream >> range;
								if ( range != "" ) {
									if ( range != "SD_RANGE" ) {
										std::cerr << "ERROR:: Bad line in CHI_STRATEGY: \"" << llstr << "\"" << std::endl;
										utility_exit_with_message( "While parsing CHI_STRATEGY:: STEP_WITHIN_SD_RANGE.  Expected to read SD_RANGE to specify the standard deviation range" );
									}
									llstream >> sd_range;
									if ( llstream.bad() ) {
										std::cerr << "ERROR:: Bad line in CHI_STRATEGY: \"" << llstr << "\"" << std::endl;
										utility_exit_with_message( "While parsing CHI_STRATEGY:: STEP_WITHIN_SD_RANGE.  Could not read standard-deviation range" );
									}
									if ( sd_range <= 0.0 ) {
										std::cerr << "ERROR:: Bad line in CHI_STRATEGY: \"" << llstr << "\"" << std::endl;
										utility_exit_with_message( "While parsing CHI_STRATEGY:: STEP_WITHIN_SD_RANGE.  Invalid standard-deviation range (must be positive). Read: " + utility::to_string( sd_range ) );
									}
								} else {
									/// Implicit sd_range of 1.0
									sd_range = 1.0;
								}
								upstream::SampleStrategyData stratdat;
								stratdat.set_strategy( upstream::rotameric_chi_step_wi_sd_range );
								stratdat.set_step_size( step_size );
								stratdat.set_sd_range( sd_range );
								set_sample_startegy_for_constraint(	ii, upres[ jj ], which_chi, stratdat );

								TR << "  ALGORITHM_INFO:: match -- chi sampling strategy STEP_WITHIN_SD_RANGE with step_size= " << step_size << " degrees across " << sd_range << " standard deviations for chi # " << which_chi << std::endl;
							} else if ( strategy == "NON_ROTAMERIC_CHI_EXPANSION" ) {
								Size npossiblerots( 0 );
								/// Expand the number of rotamers for a non-rotameric chi.
								/// 1st check to make sure that this chi is nonrotameric!
								{/// SCOPE

								using namespace core::scoring;
								using namespace core::pack::dunbrack;
								RotamerLibrary const & rotlib( * core::pack::dunbrack::RotamerLibrary::get_instance() );
								SingleResidueRotamerLibraryCOP res_rotlib( rotlib.get_rsd_library( *upres[ jj ] ) );

								if ( res_rotlib != 0 ) {

									SingleResidueDunbrackLibraryCOP dun_rotlib(
										utility::pointer::dynamic_pointer_cast< core::pack::dunbrack::SingleResidueDunbrackLibrary const > ( res_rotlib ));

									if ( dun_rotlib == 0 ) {
										utility_exit_with_message( "Failed to retrieve a Dunbrack rotamer library for AA: " +
											utility::to_string( upres[ jj ]->aa() ) +  " named " +  upres[ jj ]->name() );
									}
									npossiblerots = dun_rotlib->n_rotamer_bins();
									if ( npossiblerots == 0 ) {
										std::cerr << "Error while reading line: " << llstr << std::endl;
										utility_exit_with_message( "Rotamer library for " + upres[ jj ]->name() + " says it contains no rotamers" );
									}
									if ( which_chi != dun_rotlib->nchi() ) {
										utility_exit_with_message( "Cannot treat chi " + utility::to_string( which_chi ) +
											" on residue " + upres[ jj ]->name() +
											" as non-rotameric since there are " + utility::to_string( dun_rotlib->nchi() ) +
											" chi in the library, and the last chi is the only one that could be non-rotameric" );
									}
									bool failed_cast = false;
									switch ( dun_rotlib->nchi() ) {
										case 2: {
											switch ( dun_rotlib->nbb() ) {
												case 1: {
													SemiRotamericSingleResidueDunbrackLibrary< ONE, ONE > const * sr2 =
													dynamic_cast< SemiRotamericSingleResidueDunbrackLibrary< ONE, ONE > const * >
													( dun_rotlib.get() );
													failed_cast = sr2 == 0;
												}
												case 2: {
													SemiRotamericSingleResidueDunbrackLibrary< ONE, TWO > const * sr2 =
													dynamic_cast< SemiRotamericSingleResidueDunbrackLibrary< ONE, TWO > const * >
													( dun_rotlib.get() );
													failed_cast = sr2 == 0;
												}
												case 3: {
													SemiRotamericSingleResidueDunbrackLibrary< ONE, THREE > const * sr2 =
													dynamic_cast< SemiRotamericSingleResidueDunbrackLibrary< ONE, THREE > const * >
													( dun_rotlib.get() );
													failed_cast = sr2 == 0;
												}
												case 4: {
													SemiRotamericSingleResidueDunbrackLibrary< ONE, FOUR > const * sr2 =
													dynamic_cast< SemiRotamericSingleResidueDunbrackLibrary< ONE, FOUR > const * >
													( dun_rotlib.get() );
													failed_cast = sr2 == 0;
												}
												default: {
													utility_exit_with_message( "While parsing CHI_STRATEGY::NON_ROTAMERIC_CHI_EXPANSION\n"
																			  "All semi-rotameric libraries have 1 - 4 bb, but the library for "+
																			  upres[ jj ]->name() + " has " + utility::to_string( dun_rotlib->nbb() ) + " bb." );
												}
											}
										} break;
										case 3: {
											switch ( dun_rotlib->nbb() ) {
												case 1: {
													SemiRotamericSingleResidueDunbrackLibrary< TWO, ONE > const * sr2 =
													dynamic_cast< SemiRotamericSingleResidueDunbrackLibrary< TWO, ONE > const * >
													( dun_rotlib.get() );
													failed_cast = sr2 == 0;
												}
												case 2: {
													SemiRotamericSingleResidueDunbrackLibrary< TWO, TWO > const * sr2 =
													dynamic_cast< SemiRotamericSingleResidueDunbrackLibrary< TWO, TWO > const * >
													( dun_rotlib.get() );
													failed_cast = sr2 == 0;
												}
												case 3: {
													SemiRotamericSingleResidueDunbrackLibrary< TWO, THREE > const * sr2 =
													dynamic_cast< SemiRotamericSingleResidueDunbrackLibrary< TWO, THREE > const * >
													( dun_rotlib.get() );
													failed_cast = sr2 == 0;
												}
												case 4: {
													SemiRotamericSingleResidueDunbrackLibrary< TWO, FOUR > const * sr2 =
													dynamic_cast< SemiRotamericSingleResidueDunbrackLibrary< TWO, FOUR > const * >
													( dun_rotlib.get() );
													failed_cast = sr2 == 0;
												}
												default: {
													utility_exit_with_message( "While parsing CHI_STRATEGY::NON_ROTAMERIC_CHI_EXPANSION\n"
																			  "All semi-rotameric libraries have 1 - 4 bb, but the library for "+
																			  upres[ jj ]->name() + " has " + utility::to_string( dun_rotlib->nbb() ) + " bb." );
												}
											}
										} break;
										default: {
											utility_exit_with_message( "While parsing CHI_STRATEGY::NON_ROTAMERIC_CHI_EXPANSION\n"
												"All semi-rotameric libraries have 2 or 3 chi, but the library for "+
												upres[ jj ]->name() + " has " + utility::to_string( dun_rotlib->nchi() ) + " chi." );
										}
									}
									if ( failed_cast ) {
										utility_exit_with_message( "While parsing CHI_STRATEGY::NON_ROTAMERIC_CHI_EXPANSION\n"
											"Failed to find a semi-rotameric rotamer library for " + upres[ jj ]->name() +
											" (Did you forget the -dun10 flag?)\n"
											"The following amino acids define semi-rotameric rotamer libraries: DEFHNQ" );
									}
								} else {
									utility_exit_with_message( "While parsing CHI_STRATEGY::NON_ROTAMERIC_CHI_EXPANSION\n"
										"Failed to find a rotamer library for " + upres[ jj ]->name() );
								}
								} // scope to check we're looking at a semi-rotameric rotamer library.
								std::string nsamps;
								Size nsamples;
								if ( ! llstream.good() ) {
									std::cerr << "Error reading line: " << llstr << std::endl;
									utility_exit_with_message( "While parsing CHI_STRATEGY::NON_ROTAMERIC_CHI_EXPANSION\n"
										"Expected to read N_SAMPLES after NON_ROTAMERIC_CHI_EXPANSIONbut reached an end of line");
								}
								llstream >> nsamps;
								if ( nsamps != "N_SAMPLES" ) {
									std::cerr << "Error reading line: " << llstr << std::endl;
									utility_exit_with_message( "While parsing CHI_STRATEGY::NON_ROTAMERIC_CHI_EXPANSION\n"
										"Expected to read N_SAMPLES after NON_ROTAMERIC_CHI_EXPANSION but found '" + nsamps + "'");
								}
								llstream >> nsamples;
								if ( !llstream ) {
									std::cerr << "Error reading line: " << llstr << std::endl;
									utility_exit_with_message( "While parsing CHI_STRATEGY::NON_ROTAMERIC_CHI_EXPANSION\n"
										"Expected to read an integer following N_SAMPLES but could not" );
								}
								std::string minprob;
								Real minprobability( 2.0 / npossiblerots );
								llstream >> minprob;
								if ( minprob != "" ) {
									if ( minprob != "REQUISIT_PROBABILITY" ) {
										std::cerr << "Error reading line: " << llstr << std::endl;
										utility_exit_with_message( "While parsing CHI_STRATEGY::NON_ROTAMERIC_CHI_EXPANSION\n"
											"Expected end-of-line or 'REQUISIT_PROBABILITY' after reading the number of samples, but found '" +  minprob + "'" );
									}
									if ( ! llstream.good() ) {
										std::cerr << "Error reading line: " << llstr << std::endl;
										utility_exit_with_message( "While parsing CHI_STRATEGY::NON_ROTAMERIC_CHI_EXPANSION\n"
											"Expected to read a probability following the 'REQUISIT_PROBABILITY' string but found an end-of-line" );
									}
									llstream >> minprobability;
									if ( llstream.bad() ) {
										std::cerr << "Error reading line: " << llstr << std::endl;
										utility_exit_with_message( "While parsing CHI_STRATEGY::NON_ROTAMERIC_CHI_EXPANSION\n"
											"Expected to read a probability following the 'REQUISIT_PROBABILITY' string but could not" );
									}
								}
								// OK -- Now set the chi-sample strategy.
								upstream::SampleStrategyData stratdat;
								stratdat.set_strategy( upstream::nonrotameric_chi_sample_wi_nrchi_bin );
								stratdat.set_n_samples_per_side_of_nrchi_bin( nsamples );
								stratdat.set_nrchi_prob_minimum_for_extra_samples( minprobability );
								set_sample_startegy_for_constraint(	ii, upres[ jj ], which_chi, stratdat );
								TR << "  ALGORITHM_INFO:: match -- chi sampling strategy NON_ROTAMERIC_CHI_EXPANSION with nsteps= " << nsamples << " for rotamers with a probability better than " << minprobability << std::endl;
							} else {
								utility_exit_with_message( "While parsing CHI_STRATEGY:: unsupported sample strategy: " + strategy  + " for chi " + utility::to_string( which_chi ) );
							}
						} else if ( first == "SECONDARY_MATCH:" ) {
							if ( ii == 1 ) {
								std::cerr << "ERROR Reading line " << llstr << " " << " for geometric constraint " << ii << std::endl;
								utility_exit_with_message( "Seconary matching cannot be chosen for the first geometric constraint!" );
							}
							if ( ! llstream ) {
								utility_exit_with_message( "While parsing SECONDARY_MATCH: line, exptected to read 'UPSTREAM_CST <int>' but reached an end-of-line" );
							}
							std::string second;
							llstream >> second;
							if ( second != "UPSTREAM_CST" && second != "DOWNSTREAM" ) {
								std::cerr << "Error reading line: " << llstr << std::endl;
								utility_exit_with_message( "While parsing SECONDARY_MATCH: line, expected 'UPSTREAM_CST' or 'DOWNSTREAM' but encountered '" + second +"'." );
							}
							if ( second == "UPSTREAM_CST" ) {
								Size cst_id;
								llstream >> cst_id;
								if ( ! llstream ) {
									std::cerr << "Error reading line: " << llstr << std::endl;
									utility_exit_with_message( "While parsing SECONDARY_MATCH: line, read 'UPSTREAM_CST' and expected to read an integer following, but did not find one." );
								}
								if ( cst_id < 1 || cst_id >= ii ) {
									std::cerr << "Error reading line: '" << llstr << "' for geometric constraint " << ii << std::endl;
									utility_exit_with_message( "Secondary match algorithm requested to an upstream residue "
										"produced by an invalid geometric constraint " + utility::to_string( cst_id ) +
										"\nThe geometric constraint ID must be less than the current geometric constraint id"
										"and greater than 0" );
								}
								secondary_match_upstream_residue = true;
								secondary_match_upstream_geomcst_id = cst_id;
							}
							secondary_matching = true;

						//Author: Kui Chan Date:101309
						//Description: secondary scoring by score term(s).
						} else if ( first == "SCORING_SECMATCH::" ){
							if( !secondary_matching ){
								utility_exit_with_message( "SCORING_SECMATCH line detected without previous SECONDARY_MATCH specifier.");
							}
							sec_match_str += llstr + "\n";
						} else if ( first == "DESYMMETERIZE" ) {
							desymmeterize_upstream_restype_for_constraint( ii );
						} else {
							utility_exit_with_message( "While parsing ALGORITHM:: match data. Command '" + first +"' not supported on line '" + llstr + "'");
						}
					}
				} else {
					TR << "  Did not locate ALGORITHM_INFO:: match for constraint " << ii << std::endl;
				}


				Size const upstream_id = jj_mcfis[ kk ]->upstream_res();
				Size const downstream_id = jj_mcfis[ kk ]->downstream_res();

				toolbox::match_enzdes_util::ExternalGeomSamplerCOP exgs;
				if ( !secondary_matching ) {
					exgs = jj_mcfis[ kk ]->create_exgs();
					if ( exgs.get() == 0 ) {
						utility_exit_with_message( "ERROR: could not define external geometry between upstream and downstream residues.  All 6 parameters must be defined." );
					}
				}

				if ( jj_mcfis[ kk ]->allowed_restypes( downstream_id ).size() == 0 ) {
					utility_exit_with_message( "Input file lists no residue types for template residue " +
						utility::to_string( downstream_id ) + " for geometric constraint " +
						utility::to_string( ii ) + ".  There must be at least one." );
				}

				for ( Size ll = 1; ll <= jj_mcfis[ kk ]->allowed_restypes( downstream_id ).size(); ++ll ) {

					core::chemical::ResidueTypeCOP ll_downres( jj_mcfis[ kk ]->allowed_restypes( downstream_id )[ ll ] );

					utility::vector1< std::string > upstream_launch_atoms( 3 );
					//utility::vector1< std::string > downstream_launch_atoms( 3 ); /// TEMP HACK.  Assume single ligand downstream
					utility::vector1< core::id::AtomID > downstream_3atoms( 3 );
					for ( Size mm = 1; mm <= 3; ++mm ) downstream_3atoms[ mm ].rsd() = 1;


					utility::vector1< utility::vector1< Size > > up_ats( 3 ), down_ats( 3 );
					for ( Size mm = 1; mm <= 3; ++mm ) up_ats[ mm ] = jj_mcfis[ kk ]->template_atom_inds( upstream_id, mm, *upres[ jj ] );
					for ( Size mm = 1; mm <= 3; ++mm ) down_ats[ mm ] = jj_mcfis[ kk ]->template_atom_inds( downstream_id, mm, *ll_downres );

					runtime_assert( up_ats[ 1 ].size() == up_ats[ 2 ].size() );
					runtime_assert( up_ats[ 1 ].size() == up_ats[ 3 ].size() );

					utility::vector1< Size > up_down_ncombs( 2 );
					up_down_ncombs[ 1 ] = up_ats[ 1 ].size();
					up_down_ncombs[ 2 ] = down_ats[ 2 ].size();

					utility::LexicographicalIterator lex( up_down_ncombs );

					while ( ! lex.at_end() ) {
						for ( Size mm = 1; mm <= 3; ++mm ) upstream_launch_atoms[ mm ] = upres[ jj ]->atom_name( up_ats[ mm ][ lex[ 1 ] ] );
						for ( Size mm = 1; mm <= 3; ++mm ) downstream_3atoms[ mm ].atomno() = down_ats[ mm ][ lex[ 2 ] ];

						TR << "   " << upres[ jj ]->name() << " " << ll_downres->name() << std::endl;
						TR << "   ";
						TR << " U3: " << upstream_launch_atoms[ 3 ];
						TR << " U2: " << upstream_launch_atoms[ 2 ];
						TR << " U1: " << upstream_launch_atoms[ 1 ];
						TR << " D1: " << ll_downres->atom_name( down_ats[ 1 ][ lex[ 2 ] ] );
						TR << " D2: " << ll_downres->atom_name( down_ats[ 2 ][ lex[ 2 ] ] );
						TR << " D3: " << ll_downres->atom_name( down_ats[ 3 ][ lex[ 2 ] ] );
						TR << std::endl;

						if ( secondary_matching ) {

							if ( ii == 2 && mtask.build_round1_hits_twice() ) {
								utility_exit_with_message( "When using the build-round1-hits-twice algorithm, round 2 must also use the classic matching algorithm" );
							}

							utility::vector1< Size > candidate_atids( 3 );
							utility::vector1< Size > target_atids( 3 );
							for ( Size nn = 1; nn <= 3; ++nn ) {
								candidate_atids[ nn ] = up_ats[ nn ][ lex[ 1 ]];
								target_atids[ nn ]    = down_ats[ nn ][ lex[ 2 ]];
							}
							if ( secondary_match_upstream_residue ) {
								add_secondary_upstream_match_geometry_for_constraint(
									ii, secondary_match_upstream_geomcst_id, upres[ jj ],
									ll_downres, candidate_atids, target_atids,
									jj_mcfis[ kk ], sec_match_str, *upstream_pose_ );
							} else {
								add_secondary_downstream_match_geometry_for_constraint(
									ii, upres[ jj ], ll_downres,
									candidate_atids, target_atids, jj_mcfis[ kk ], sec_match_str, *upstream_pose_,
									jj_mcfis[kk]->is_covalent() );
							}
						} else {
							add_external_geometry_samples_for_constraint(
								ii, upres[ jj ], upstream_launch_atoms, downstream_3atoms, *exgs,
								jj_mcfis[ kk ]->index(),
								mtask.enumerate_ligand_rotamers(),
								jj_mcfis[kk]->is_covalent(),
								mtask.build_round1_hits_twice() );
						}
						++lex;
					} //lex loop
				} // ll loop over downstream residue types
			} //kk loop over mcfis
		} //jj loop over restypes for constraint
	} //ii loop over geometric constraints

	//std::cout << "APL DEBUG Matcher.cc::initialize_from_file end" << std::endl;
} //initialize_from_file function


/// @brief Main worker function
bool Matcher::find_hits()
{
	if( !initialize_scaffold_build_points() ) return false;
	//std::cout << "APL DEBUG Matcher.cc::find_hits 1" << std::endl;
	initialize_bump_grids();
	//std::cout << "APL DEBUG Matcher.cc::find_hits 2" << std::endl;
	initialize_active_site_grid();
	//std::cout << "APL DEBUG Matcher.cc::find_hits 3" << std::endl;
	initialize_occupied_space_hash();
	//std::cout << "APL DEBUG Matcher.cc::find_hits 4" << std::endl;
	initialize_downstream_algorithms();
	//std::cout << "APL DEBUG Matcher.cc::find_hits 5" << std::endl;

	return generate_hits();
}


void
Matcher::process_matches( output::MatchProcessor & processor ) const
{
	processor.begin_processing();
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		representative_downstream_algorithm_[ ii ]->prepare_for_match_enumeration( *this );
	}
	for( std::list< downstream::DownstreamBuilderOP >::const_iterator ds_it( all_downstream_builders_.begin() ),
				 ds_end( all_downstream_builders_.end()); ds_it != ds_end; ++ds_it ){
		if( (*ds_it)->hits_potentially_incompatible() ){
			check_potential_dsbuilder_incompatibility_ = true;
			break;
		}
	}

	if ( (! output_matches_as_singular_downstream_positioning_) || check_potential_dsbuilder_incompatibility_ ) {
		if( output_matches_as_singular_downstream_positioning_ ){
			TR << "Potential DownstreamBuilder hit incompatibilities have been detected. All possible hit combinations will be enumerated. This overrides the MatcherTask instruction to output matches with only a single downstream position." << std::endl;
		}
		process_matches_main_loop_enumerating_all_hit_combos( processor );
	} else {
		process_matches_where_one_geomcst_defines_downstream_location( processor );
	}
	processor.end_processing();
}

/*utility::vector1< std::list< Hit > > const &
Matcher::hits() const
{
	return hits_;
}
*/

upstream::ScaffoldBuildPointCOP
Matcher::build_point( Size index ) const
{
	return all_build_points_[ index ];
}


core::pose::PoseCOP
Matcher::upstream_pose() const
{
	return upstream_pose_;
}

core::pose::PoseCOP
Matcher::downstream_pose() const
{
	return downstream_pose_;
}

upstream::UpstreamBuilderCOP
Matcher::upstream_builder( Size cst_id ) const
{
	return upstream_builders_[ cst_id ];
}

downstream::DownstreamBuilderCOP
Matcher::downstream_builder( Size geom_cst ) const
{
	if ( downstream_builders_[ geom_cst ].empty() ) {
		return 0;
	} else {
		return *(downstream_builders_[ geom_cst ].begin());
	}
}

std::list< downstream::DownstreamAlgorithmCOP >
Matcher::downstream_algorithms( Size cst_id ) const
{
	std::list< downstream::DownstreamAlgorithmCOP > dsalgs;
	for ( std::list< downstream::DownstreamAlgorithmOP >::const_iterator
			iter = downstream_algorithms_[ cst_id ].begin(),
			iter_end = downstream_algorithms_[ cst_id ].end();
			iter != iter_end; ++iter ) {
		dsalgs.push_back( *iter );
	}
	return dsalgs;
}

downstream::DownstreamAlgorithmCOP
Matcher::representative_downstream_algorithm( Size cst_id ) const
{
	return representative_downstream_algorithm_[ cst_id ];
}


Matcher::HitList const &
Matcher::hits( Size cst_id ) const
{
	return hits_[ cst_id ];
}


OccupiedSpaceHashCOP
Matcher::occ_space_hash() const {
	return occ_space_hash_;
}

utility::vector1< upstream::ScaffoldBuildPointCOP > const &
Matcher::per_constraint_build_points( Size cst_id ) const
{
	return per_constraint_build_points_[ cst_id ];
}


upstream::ScaffoldBuildPointOP
Matcher::build_point( Size index )
{
	return all_build_points_[ index ];
}

upstream::UpstreamBuilderOP
Matcher::upstream_builder( Size cst_id )
{
	return upstream_builders_[ cst_id ];
}

bool
Matcher::has_upstream_only_geomcsts() const
{
	for( core::Size i =1; i<= n_geometric_constraints_; ++i){
		if( geomcst_is_upstream_only_[i] ) return true;
	}
	return false;
}

downstream::DownstreamBuilderOP
Matcher::downstream_builder( Size cst_id )
{
	runtime_assert( ! downstream_builders_[ cst_id ].empty() );
	return *(downstream_builders_[ cst_id ].begin());
}

std::list< downstream::DownstreamBuilderOP > const &
Matcher::downstream_builders( Size cst_id ) const
{
	return downstream_builders_[ cst_id ];
}

std::list< downstream::DownstreamAlgorithmOP > const &
Matcher::nonconst_downstream_algorithms( Size cst_id )
{
	return downstream_algorithms_[ cst_id ];
}


OccupiedSpaceHashOP
Matcher::occ_space_hash() {
	return occ_space_hash_;
}

Matcher::HitListIterator
Matcher::hit_list_begin( Size geom_cst_id )
{
	return hits_[ geom_cst_id ].begin();
}

Matcher::HitListIterator
Matcher::hit_list_end( Size geom_cst_id )
{
	return hits_[ geom_cst_id ].end();
}

void
Matcher::erase_hit(
	downstream::DownstreamAlgorithm const & dsalg,
	Size geom_cst_id_for_hit,
	HitListIterator const & iter
)
{
	hits_[ geom_cst_id_for_hit ].erase( iter );
	if ( geom_cst_id_for_hit != dsalg.geom_cst_id() ) {
		/// A downstream algorithm is making a primary modification to
		/// the hit list for some other geometric constraint.

		// ASSUMPTION: hits for geom_cst i may not depend on hits for geom_cst j, for j > i.
		// Under this assumption, the downstream algorithm for geom_cst i is unable to
		// direct the deletion of hits from geom_cst j.  Maybe the downstream algorithm
		// for geom_cst j deletes its own hits during a call to respond_to_peripheral_hitlist_change,
		// as a result of round i hits disappearing; however, round i may not delete round j's hits
		// directly.
 		runtime_assert( dsalg.geom_cst_id() > geom_cst_id_for_hit );

		note_primary_change_to_geom_csts_hitlist( geom_cst_id_for_hit );
	}
}


bool Matcher::generate_hits() {
	//std::cout << "APL DEBUG Matcher.cc::generate_hits() begin" << std::endl;
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		//std::cout << "APL DEBUG Matcher.cc::generate_hits() ii=" << ii << std::endl;
		prepare_for_hit_generation_for_constraint( ii );
		generate_hits_for_constraint( ii );
		if( !finish_hit_generation_for_constraint( ii ) ) {
			//std::cout << "APL DEBUG Matcher.cc::generate_hits() early exit ii=" << ii << std::endl;
			return false;
		}
	}
	//std::cout << "APL DEBUG Matcher.cc::generate_hits() end" << std::endl;
	return true;
}

/// @details Update the per_constraint_build_points_ array for the given constraint before
/// beginning hit generation.  This follows one of two paths depending on whether
/// same_build_resids_for_all_csts_ is active or not.  I cannot see any reason that this
/// function could not be called earlier than immediately before hit generation for a particular
/// constraint.
void Matcher::prepare_for_hit_generation_for_constraint( Size cst_id )
{
	assert( pose_build_resids_.size() == all_build_points_.size() ); // this function assumes the entries in these two vectors correspond to each other

	if ( same_build_resids_for_all_csts_ ) {
		per_constraint_build_points_[ cst_id ].reserve( all_build_points_.size() );
		for ( Size ii = 1; ii <= all_build_points_.size(); ++ii ) {
			/// if ( logic ) goes here to decide if all build points
			/// are apropriate for a pariticular geometric constraint.
			per_constraint_build_points_[ cst_id ].push_back( all_build_points_[ ii ] );
		}
	} else {
		per_constraint_build_points_[ cst_id ].reserve( per_cst_build_resids_[ cst_id ].size() );

		/// "merge sort" inspired algorithm; advance through two sorted lists.
		Size counter = 1;
		for ( Size ii = 1; ii <= pose_build_resids_.size(); ++ii ) {
			if ( per_cst_build_resids_[ cst_id ][ counter ] == pose_build_resids_[ ii ] ) {
				per_constraint_build_points_[ cst_id ].push_back( all_build_points_[ ii ] );
				++counter;
				if ( counter > per_cst_build_resids_[ cst_id ].size() ) break;
			}
		}
	}

}

//Author: Kui Chan
//access function to pose_build_resids_
//Reason: Use to update the SecondaryMatcherToUpstreamResidue hit.second()
utility::vector1< core::Size > const &
Matcher::get_pose_build_resids() const{
	return pose_build_resids_;
}

void Matcher::generate_hits_for_constraint( Size cst_id )
{
	/// At the conclusion of hit generation, there will be new hits for this geometric constraint;
	/// put this constraint ID at the front of the list of geom-csts to trigger primary-change responses.
	note_primary_change_to_geom_csts_hitlist( cst_id );

	std::list< Hit > hits = representative_downstream_algorithm_[ cst_id ]->build_hits_at_all_positions( *this );
	hits_[ cst_id ].splice( hits_[ cst_id ].end(), hits );

}

/// @details just like generate_hits_for_constraint, but without calling not_primary_change_to_geom_csts_hitlist.
/// Used with the "build round-1 hits twice" scheme.
void Matcher::regenerate_round1_hits()
{
	core::Size const cst_id = 1;
	std::list< Hit > hits = representative_downstream_algorithm_[ cst_id ]->build_hits_at_all_positions( *this );
	hits_[ cst_id ].splice( hits_[ cst_id ].end(), hits );
}

bool
Matcher::finish_hit_generation_for_constraint( Size cst_id )
{

	/// During the primary and peripheral hit-list prunings, new elements may be pushed-back
	/// into the hit-lists-modified-by-primary-deletion list.
	for ( std::list< Size >::iterator iter = geom_csts_with_primary_hitlist_modificiations_.begin();
			iter != geom_csts_with_primary_hitlist_modificiations_.end(); /* no increment */ ) {

		representative_downstream_algorithm_[ *iter ]->respond_to_primary_hitlist_change( *this, cst_id );

		/// There is an implicit dependency DAG between geometric constraints:
		/// round i cannot be dependent on round j, if j > i.
		/// Therefore a very simple "topological sort" may be performed offline to come up with
		/// the appropriate order in which to update peripheral hits: go from 1 to n.
		for ( Size jj = 1; jj <= cst_id; ++jj ) {
			if ( jj == *iter ) continue;
			representative_downstream_algorithm_[ jj ]->respond_to_peripheral_hitlist_change( *this );
		}

		std::list< Size >::iterator iter_next( iter );
		++iter_next;
		geom_cst_has_primary_modification_[ *iter ] = false;
		geom_csts_with_primary_hitlist_modificiations_.erase( iter );
		iter = iter_next;
	}

	//if there are no hits in the list for this constraint, we can abort preemptively
	bool hits_found( hits_[cst_id].size() != 0 );

	if ( cst_id == 2 && build_round1_hits_twice_ ) {
		regenerate_round1_hits();
		TR << "Regenerated " << hits_[ 1 ].size() << " hits for round 1" << std::endl;
	}

	if ( ( cst_id == n_geometric_constraints_ ) || (!hits_found) ) {
		/// We're completely done with hit generation; clean up.
		/// Delete the occ_space_hash_
		for ( std::list< downstream::DownstreamBuilderOP >::const_iterator
				iter = all_downstream_builders_.begin(),
				iter_end = all_downstream_builders_.end();
				iter != iter_end; ++iter ) {
			(*iter)->set_occupied_space_hash( 0 );
		}
		occ_space_hash_.reset();
	}
	return hits_found;
}

/// @details this function returns false in case
/// there are no build points for a certain cst.
/// this can happen if MPMs in the MatcherTask
/// ruled out all user defined build points
bool
Matcher::initialize_scaffold_build_points()
{
	runtime_assert( upstream_pose_.get() );
	if( same_build_resids_for_all_csts_){
		if( pose_build_resids_.size() == 0 ) {
			TR << "WARNING: No build points were set in the matcher, could not initialize scaffold build points." << std::endl;
			return false;
		}
	}
	else{
		for( core::Size i =1; i <= per_cst_build_resids_.size(); ++i ){
			if( per_cst_build_resids_[i].size() == 0 ){
				TR << "WARNING: No build points for geomcst " << i << " were set in the matcher, could not initialize scaffold build points." << std::endl;
				return false;
			}
		}
	}

	all_build_points_.resize( pose_build_resids_.size() );
	for ( Size ii = 1; ii <= pose_build_resids_.size(); ++ii ) {
		runtime_assert_msg( pose_build_resids_[ ii ] <= upstream_pose_->n_residue(),
		                    "pos file contains position outside of valid range.");
		all_build_points_[ ii ] = protocols::match::upstream::ScaffoldBuildPointOP( new upstream::OriginalBackboneBuildPoint(
			upstream_pose_->residue( pose_build_resids_[ ii ] ), ii ) );
	}
	return true;
}


void
Matcher::initialize_bump_grids()
{
	runtime_assert( upstream_pose_.get() );

	clock_t starttime = clock();
	TR << "Initializing BumpGrids... " << std::endl;

	if ( ! bb_grid_ ) {
		bb_grid_ = bump_grid_to_enclose_pose( *upstream_pose_ );
	}

	/// This code is fixed-backbone only... it needs to be expanded.
	original_scaffold_residue_bump_grids_.resize( upstream_pose_->total_residue() );
	for ( Size ii = 1; ii <= upstream_pose_->total_residue(); ++ii ) {
		BumpGridOP resbgop = bump_grid_to_enclose_residue_backbone( upstream_pose_->residue( ii ), *bb_grid_ );
		fill_grid_with_backbone_heavyatom_spheres( upstream_pose_->residue( ii ), *resbgop );
		bb_grid_->or_with( *resbgop );
		original_scaffold_residue_bump_grids_[ ii ] = resbgop;
	}


	// Inform everyone of the grid
	for ( std::list< downstream::DownstreamBuilderOP >::const_iterator
			iter = all_downstream_builders_.begin(),
			iter_end = all_downstream_builders_.end();
			iter != iter_end; ++iter ) {
		(*iter)->set_bb_grid( bb_grid_ );
	}
	for ( std::list< downstream::DownstreamAlgorithmOP >::const_iterator
			iter = all_downstream_algorithms_.begin(),
			iter_end = all_downstream_algorithms_.end();
			iter != iter_end; ++iter ) {
		(*iter)->set_bb_grid( bb_grid_ );
	}
	for ( utility::vector1< upstream::UpstreamBuilderOP >::const_iterator
			iter = upstream_builders_.begin(),
			iter_end = upstream_builders_.end();
			iter != iter_end; ++iter ) {
		(*iter)->set_bb_grid( bb_grid_ );
	}

	TR << "...done" << std::endl;
	clock_t stoptime = clock();
	TR << " TIMING: Bump grids took " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " seconds to compute" << std::endl;
}

void
Matcher::initialize_active_site_grid()
{
	/* TEMP
	if ( downstream_atoms_required_inside_active_site_.empty() ) return;

	if ( upstream_resids_and_radii_defining_active_site_.empty() ) {
		utility_exit_with_message( "ERROR: Active site undefined, yet downstream atoms are required to be in active site" );
	}*/

	if ( read_gridlig_file_ ) {
		active_site_grid_ = downstream::ActiveSiteGridOP( new downstream::ActiveSiteGrid );
		active_site_grid_->initialize_from_gridlig_file( gridlig_fname_ );

		/*Bool3DGridKinemageWriter writer;
		writer.set_write_facets( true );
		std::ofstream ostr( "active_site_gridlig.kin" );
		writer.set_line_color( "green" );
		writer.write_grid_to_kinemage( ostr, "act_site", active_site_grid_->grid() );*/
	} else {
		active_site_grid_ = downstream::ActiveSiteGridOP( new downstream::ActiveSiteGrid );
		active_site_grid_->set_bin_width( 0.25 ); /// Same resolution as the bump grid!

		for ( std::list< std::pair< Size, Real > >::const_iterator
				iter = upstream_resids_and_radii_defining_active_site_.begin(),
				iter_end = upstream_resids_and_radii_defining_active_site_.end();
				iter != iter_end; ++iter ) {
			runtime_assert( iter->first <= upstream_pose_->total_residue() );
			active_site_grid_->enlargen_to_capture_volume_within_radius_of_residue(
				upstream_pose_->residue( iter->first ), iter->second );
		}

		for ( std::list< std::pair< Size, Real > >::const_iterator
				iter = upstream_resids_and_radii_defining_active_site_.begin(),
				iter_end = upstream_resids_and_radii_defining_active_site_.end();
				iter != iter_end; ++iter ) {
			active_site_grid_->or_within_radius_of_residue(
				upstream_pose_->residue( iter->first ), iter->second );
		}

		// Inform everyone of the grid
		for ( std::list< downstream::DownstreamBuilderOP >::const_iterator
				iter = all_downstream_builders_.begin(),
				iter_end = all_downstream_builders_.end();
				iter != iter_end; ++iter ) {
			(*iter)->set_active_site_grid( active_site_grid_ );
		}
		for ( std::list< downstream::DownstreamAlgorithmOP >::const_iterator
				iter = all_downstream_algorithms_.begin(),
				iter_end = all_downstream_algorithms_.end();
				iter != iter_end; ++iter ) {
			(*iter)->set_active_site_grid( active_site_grid_ );
		}

		/// Inform dowsntream builders of their atoms that must be contained inside the
		/// active site grid.
		for ( std::list< downstream::DownstreamBuilderOP >::const_iterator
				iter = all_downstream_builders_.begin(),
				iter_end = all_downstream_builders_.end();
				iter != iter_end; ++iter ) {
			for ( std::list< core::id::AtomID >::const_iterator
					atid_iter = downstream_atoms_required_inside_active_site_.begin(),
					atid_iter_end = downstream_atoms_required_inside_active_site_.end();
					atid_iter != atid_iter_end; ++atid_iter ) {
				(*iter)->require_atom_to_reside_in_active_site( *atid_iter );
			}
		}

		/*std::cout << "Writing active-site kinemage" << std::endl;

		Bool3DGrid copy_active = active_site_grid_->grid();
		copy_active.subtract( bb_grid_->grid( C_ALA ) );

		Bool3DGridKinemageWriter writer;
		writer.set_write_facets( true );
		std::ofstream ostr( "active_site_grid.kin" );
		writer.set_line_color( "green" );
		writer.write_grid_to_kinemage( ostr, "act_site", copy_active );
		*/
	}
}

void
Matcher::initialize_occupied_space_hash()
{
	occ_space_hash_ = OccupiedSpaceHashOP( new OccupiedSpaceHash );
	occ_space_hash_->set_bounding_box( occ_space_bounding_box_ );
	occ_space_hash_->set_xyz_bin_widths( euclidean_bin_widths_ );
	occ_space_hash_->set_euler_bin_widths( euler_bin_widths_ );

	occ_space_hash_->initialize();

}

void
Matcher::initialize_downstream_algorithms()
{
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		geomcst_is_upstream_only_[ ii ] = representative_downstream_algorithm_[ ii ]->upstream_only();
	}
}

downstream::DownstreamBuilderOP
Matcher::create_ds_builder(
	Size const cst_id,
	core::chemical::ResidueTypeCOP restype,
	utility::vector1< std::string >  const & upstream_launch_atoms,
	utility::vector1< core::id::AtomID > const & downstream_3atoms,
	bool enumerate_ligand_rotamers,
	bool catalytic_bond
)
{
	runtime_assert( upstream_launch_atoms.size() == 3 );
	runtime_assert( downstream_3atoms.size() == 3 );

	runtime_assert( build_set_id_for_restype_[ cst_id ].find( restype->name() )
		!= build_set_id_for_restype_[ cst_id ].end() );

	//Size build_set_id = build_set_id_for_restype_[ cst_id ][ restype->name() ];
	runtime_assert( dynamic_cast< upstream::ProteinUpstreamBuilder * > ( upstream_builders_[ cst_id ].get() ) );
	upstream::ProteinUpstreamBuilderOP prot_sc_builder( utility::pointer::static_pointer_cast< upstream::ProteinUpstreamBuilder > ( upstream_builders_[ cst_id ] ));

	upstream::BuildSet & build_set = prot_sc_builder->build_set( restype );

	runtime_assert( build_set.has_restype() );

	/// Only supports rigid-ligand builders for now... This code will expand in the future.
	runtime_assert( downstream_pose_ != 0 );
	runtime_assert( downstream_pose_->total_residue() == 1 );

	for ( Size ii = 1; ii <= 3; ++ii ) {
		runtime_assert( downstream_3atoms[ ii ].rsd() == 1 );
		runtime_assert( downstream_3atoms[ ii ].atomno() <= downstream_pose_->residue( 1 ).natoms() );
	}

	downstream::DownstreamBuilderOP builder;
	if ( ! enumerate_ligand_rotamers ) {

		downstream::RigidLigandBuilderOP rigid_builder( new downstream::RigidLigandBuilder );
		rigid_builder->ignore_h_collisions( true );
		rigid_builder->initialize_from_residue(
			downstream_3atoms[ 1 ].atomno(),
			downstream_3atoms[ 2 ].atomno(),
			downstream_3atoms[ 3 ].atomno(),
			downstream_orientation_atoms_[ 1 ].atomno(),
			downstream_orientation_atoms_[ 2 ].atomno(),
			downstream_orientation_atoms_[ 3 ].atomno(),
			downstream_pose_->residue(1) );

		if ( catalytic_bond ) {
			using namespace core::scoring::etable::count_pair;
			utility::vector1< std::pair< Size, Size > > bond_list;

			Size upstream_atom_id = build_set.restype().atom_index( upstream_launch_atoms[ 1 ] );
			Size downstream_atom_id = downstream_3atoms[ 1 ].atomno();

			bond_list.push_back( std::make_pair( upstream_atom_id, downstream_atom_id ) );

			CountPairGenericOP cpgen( new CountPairGeneric(
				build_set.restype(),
				downstream_pose_->residue_type(1),
				bond_list ) );
			/// Unclear what the xover value should be... 3 ignores collisions for
			/// atoms that are 3 bonds apart
			cpgen->set_crossover( 3 );

			rigid_builder->initialize_upstream_residue( build_set.restype().get_self_ptr(), cpgen );
		} else {
			rigid_builder->initialize_upstream_residue( build_set.restype().get_self_ptr() );
		}
		builder = rigid_builder;
	} else {
		downstream::LigandConformerBuilderOP ligand_rotamer_builder( new downstream::LigandConformerBuilder );
		ligand_rotamer_builder->ignore_h_collisions( true );
		ligand_rotamer_builder->initialize_from_residue(
			downstream_3atoms[ 1 ].atomno(),
			downstream_3atoms[ 2 ].atomno(),
			downstream_3atoms[ 3 ].atomno(),
			downstream_orientation_atoms_[ 1 ].atomno(),
			downstream_orientation_atoms_[ 2 ].atomno(),
			downstream_orientation_atoms_[ 3 ].atomno(),
			downstream_pose_->residue(1) );

		//note: if the relevant downstream atoms have been set,
		//this means that the downsteam conformer builder should split
		//up it's rotamer library accordingly
		if( relevant_downstream_atoms_.size() != 0 ){
			//std::cerr << "determining redundant conformers " << std::endl;
			utility::vector1< core::Size > relevant_atom_indices;
			for( core::Size i = 1; i <= relevant_downstream_atoms_.size(); ++i){
				relevant_atom_indices.push_back( relevant_downstream_atoms_[i].atomno() );
			}
			ligand_rotamer_builder->determine_redundant_conformer_groups( relevant_atom_indices );
		}

		/// Refactor this!
		if ( catalytic_bond ) {
			using namespace core::scoring::etable::count_pair;
			utility::vector1< std::pair< Size, Size > > bond_list;

			Size upstream_atom_id = build_set.restype().atom_index( upstream_launch_atoms[ 1 ] );
			Size downstream_atom_id = downstream_3atoms[ 1 ].atomno();

			bond_list.push_back( std::make_pair( upstream_atom_id, downstream_atom_id ) );

			CountPairGenericOP cpgen( new CountPairGeneric(
				build_set.restype(),
				downstream_pose_->residue_type(1),
				bond_list ) );
			/// Unclear what the xover value should be... 3 ignores collisions for
			/// atoms that are 3 bonds apart
			cpgen->set_crossover( 3 );

			ligand_rotamer_builder->initialize_upstream_residue( build_set.restype().get_self_ptr(), cpgen );
		} else {
			ligand_rotamer_builder->initialize_upstream_residue( build_set.restype().get_self_ptr() );
		}

		builder = ligand_rotamer_builder;
	}

	downstream_builders_[ cst_id ].push_back( builder );
	all_downstream_builders_.push_back( builder );
	return builder;
}


/// @details Subsample all the available hits for a single voxel in 6D to select
/// 10 or fewer hits for each geometric constraint.  If there are already 10 or
/// fewer hits for geom-cst ii, then it takes them all.  Otherwise, it selects
/// 10 representatives for each upstream conformation.  For example,
/// if all the hits in this voxel come from a single rotamer of the upsteram hit,
/// then exactly 10 upstream hits will be chosen.  The purpose of this code is
/// to speed up match enumeration when it might otherwise get very bogged down
/// by the combinatorics
void
Matcher::select_hit_representatives(
	utility::vector1< utility::vector1< Hit const * > > const & hit_vectors,
	utility::vector1< Size > & n_hits_per_geomcst,
	utility::vector1< utility::vector1< Size > > & reps
) const
{

	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii  ) {
		if ( ! dynamic_grid_refinement_ && n_hits_per_geomcst[ ii ] > 10 ) { // 10 is arbitrary
			Size max_hits_per_us_hit(10); // 10 is arbitrary -- consider a smaller number here

			std::map< upstream_hit, std::set< Size > > us_hit_map;

			for( Size jj = 1; jj <= n_hits_per_geomcst[ ii ]; ++jj ) {
				upstream_hit this_us_hit( *hit_vectors[ii][jj] );

				std::map< upstream_hit, std::set< Size > >::iterator hitmap_it(us_hit_map.find( this_us_hit ));
				if( hitmap_it == us_hit_map.end() ){
					std::set< Size > hits_this_us_hit;
					hits_this_us_hit.insert( jj );
					us_hit_map.insert( std::pair< upstream_hit, std::set< Size > >( this_us_hit, hits_this_us_hit ) );
				}
				else{
					hitmap_it->second.insert( jj );
				}
			} //loop over all hits

			//Size unique_us_hits = us_hit_map.size();
			n_hits_per_geomcst[ ii ] = 0;
			reps[ii].clear();
			for( std::map< upstream_hit, std::set< Size> >::const_iterator map_it( us_hit_map.begin() ), map_end( us_hit_map.end() );
						 map_it != map_end; ++map_it ){
				Size counter(0);
				for( std::set< Size >::const_iterator set_it( map_it->second.begin() ), set_end( map_it->second.end() );
						 set_it != set_end; ++set_it ){
					++counter;
					++n_hits_per_geomcst[ ii ];
					reps[ii].push_back( *set_it );
					if( counter >= max_hits_per_us_hit ) break;
				}
			}
		} else { // n_hits_per_geomcst[ ii ] < 10
			reps[ ii ].resize( n_hits_per_geomcst[ ii ] );
			for ( Size jj = 1; jj <= n_hits_per_geomcst[ ii ]; ++jj ) {
				reps[ ii ][ jj ] = jj;
			}
		}
	}
}

/// @brief Returns false if all non-upstream-only hits are compatible.
/// Returns true if any non-upstream-only hits are incompatible, and increments
/// the lexicographical iterator at the most-significant dimension possible
bool
Matcher::check_non_upstream_only_hit_incompatibility(
	match_dspos1 const & m1,
	utility::LexicographicalIterator & lex,
	output::MatchProcessor const & processor
) const
{
	/// Before descending into the secondary matches, check that none of the non-upstream-only
	/// hits are incompatible with each other.  Once we're iterating over the upstream-only
	/// hits, we'll assume that all the non-upstream-only hits are compatible with each other.
	/// Structure this loop so that the "earliest" incompatibility is found.
	/// ii will iterate from 2 to n_geometrict_constraints and look at whether
	/// it is incomaptible with any of the hits from 1 to ii-1; if it is,
	/// then we can increment the lexicographical iterator at position ii.
	/// This guarantees the most-significant dimension producing an incompatibility
	/// will be advanced, skipping the largest possible number of (incompatible)
	/// hit combinations.
	for ( Size ii = 2; ii <= n_geometric_constraints_; ++ii ) {
		if ( geomcst_is_upstream_only_[ ii ] ) continue; // ignore upstream-only hits

		Size ii_bp( m1.upstream_hits[ ii ].scaffold_build_id() );
		for ( Size jj = 1; jj < ii; ++jj ) {
			if ( geomcst_is_upstream_only_[ jj ] ) continue; // ignore upstream-only hits

			Size jj_bp( m1.upstream_hits[ jj ].scaffold_build_id() );
			if ( ! all_build_points_[ ii_bp ]->compatible( *all_build_points_[ jj_bp ]) ) {
				lex.continue_at_dimension( ii );
				//std::cout << "Incompatible at 1153 " << ii << " " << jj << " " << ii_bp << " " << jj_bp << std::endl;
			}
			if ( ! upstream_builders_[ ii ]->compatible( fake_hit( m1.upstream_hits[ ii ] ),
					*all_build_points_[ ii_bp ], *upstream_builders_[ jj ],
					fake_hit( m1.upstream_hits[ jj ] ), *all_build_points_[ jj_bp ] ))  {
				lex.continue_at_dimension( ii );
				//std::cout << "Incompatible at 1159 " << ii << " " << jj << " " << ii_bp << " " << jj_bp << std::endl;
				return true;
			}
			if( !processor.up_coll_filt()->passes_hardsphere_filter( ii, jj, fake_hit( m1.upstream_hits[ ii ] ), fake_hit( m1.upstream_hits[ jj ] ) ) ){
				lex.continue_at_dimension( ii );
				return true;
			}
		}
	}
	//std::cout << "Compatible at 1377" << std::endl;
	return false; // no incompatibility
}


/// @brief very similar to above function, the difference being
/// that in checks whether all the downstream builders agree that
/// their hits are compatible with each other
/// flo sep'13 there was a problem with this in that for secondary
/// downstream matching the mather doesn't have downstream builders
/// i.e., compatibility can't be checked.
/// bugfix is to also as the downstream algorithm for its dsbuilder
/// another potential fix could be to set the dsbuilders in the matcher
/// itself when the algorithms are created, but i'm not sure if this would
/// interfere with other things?
bool
Matcher::check_downstream_hit_incompatibility(
	match const & m,
	utility::LexicographicalIterator & lex
) const
{
	for ( Size ii = 2; ii <= n_geometric_constraints_; ++ii ) {
		if ( geomcst_is_upstream_only_[ ii ] ) continue; // ignore upstream-only hits

		downstream::DownstreamBuilderCOP ii_dsbuilder( downstream_builders_[ii].size() == 0 ? representative_downstream_algorithm_[ ii ]->get_dsbuilder() : *(downstream_builders_[ii].begin()) );
		if( !ii_dsbuilder) continue; // can't check compatibility if there are no downstream builders

		for ( Size jj = 1; jj < ii; ++jj ) {
			if ( geomcst_is_upstream_only_[ jj ] ) continue; // ignore upstream-only hits
			downstream::DownstreamBuilderCOP jj_dsbuilder( downstream_builders_[jj].size() == 0 ? representative_downstream_algorithm_[ jj ]->get_dsbuilder() : *(downstream_builders_[jj].begin()) );
			if( !jj_dsbuilder) continue; // can't check compatibility if there are no downstream builders

			if( ! ii_dsbuilder->compatible( m[ii], *jj_dsbuilder ,m[jj] ) ){
				lex.continue_at_dimension( ii );
				return true;
			}
		} // jj loop
	} //ii loop
	return false;
}

/// @details returns false if all upstream-only hits in a particular combination
/// are compatible with each other and with the non-upstream-only hits.  Returns
/// true if any are incompatible.  Ignores compatibility between non-upstream-only hits.
/// If there is an incompatibility, the upstream_only_hit_iterators are advanced.
/// In the event that this increment beyond the incompatible upstream-only-hit combinations
/// advance the most-significant upstream-only geometric-constraint's iterator
/// to its list end, then this function sets the value of
/// last_upstream_only_geomcst_advanced to zero.
/// update flo nov 10:
/// this function now also checks whether upstream only csts clash with any of the downstream
/// objects. it's faster to do this here than in the filters
bool
Matcher::test_upstream_only_hit_incompatibility(
	match const & m,
	utility::vector1< HitPtrListCOP > const & upstream_only_hits,
	utility::vector1< std::list< Hit const * >::const_iterator > & upstream_only_hit_iterators,
	Size & last_upstream_only_geomcst_advanced,
	output::MatchProcessor const & processor
) const
{
	/// Determine if any of the secondary hits are incompatible and therefore
	/// could not produce a match.
	/// If we do find an incompatibility, then "continue" at the most-significant
	/// dimension.
	bool incompatible( false ), outstanding_list_increment( false );
	for ( Size ii = 2; ii <= n_geometric_constraints_; ++ii ) {
		Size ii_bp( m[ ii ].scaffold_build_id() );
		for ( Size jj = 1; jj < ii; ++jj ) {
			Size jj_bp( m[ jj ].scaffold_build_id() );

			/// We have already checked for compatibility between non-upstream-only matches.
			if ( ! geomcst_is_upstream_only_[ ii ] && ! geomcst_is_upstream_only_[ jj ] ) continue;

			if ( ! all_build_points_[ ii_bp ]->compatible( *all_build_points_[ jj_bp ] )) {
				incompatible = true;
			}
			if ( ! incompatible && ! upstream_builders_[ ii ]->compatible(
					m[ ii ], *all_build_points_[ ii_bp ],
					*upstream_builders_[ jj ],  m[ jj ], *all_build_points_[ jj_bp ] ))  {
				incompatible = true;
			}

			//clash check
			//better here than later to alleviate combinatorics
			if( !incompatible && geomcst_is_upstream_only_[ii] && representative_downstream_algorithm_[jj]->generates_primary_hits() ){
				if( !processor.up_down_filt()->passes_hardsphere_filter( ii, jj, m[ii], m[jj] ) ) incompatible = true;
			}

			if( !incompatible && ( geomcst_is_upstream_only_[ii] || geomcst_is_upstream_only_[jj] ) ){
				if( !processor.up_coll_filt()->passes_hardsphere_filter( ii, jj, m[ii], m[jj] ) ) incompatible = true;
			}
			//clash check over

			if ( incompatible ) {
				/// Increment the most-significant geom-cst to alleviate the incompatibility
				/// either ii or jj represents an upstream-only geom-cst.
				/// if it's not ii, then we start our incrementing at jj.
				for ( Size kk = ( ! geomcst_is_upstream_only_[ ii ] ? jj : ii ); kk >= 1; --kk ) {
					if ( geomcst_is_upstream_only_[ kk ] ) {
						++upstream_only_hit_iterators[ kk ];
						last_upstream_only_geomcst_advanced = kk;
						if ( upstream_only_hit_iterators[ kk ] == upstream_only_hits[ kk ]->val().end() ) {
							outstanding_list_increment = true;
						} else {
							outstanding_list_increment = false;
							break;
						}
					}
				}
				//std::cout << "Incompatible at 1222 " << ii << " " << jj << std::endl;
				if ( outstanding_list_increment ) {
					// indicate to the calling function that we've visited all of the upstream-only hit combinations.
					last_upstream_only_geomcst_advanced = 0;
				}
				return true;
			} // if incompatible
		}
	}
	//std::cout << "Compatible at 1231" << std::endl;
	return false;
}

/// @details same as above
bool
Matcher::test_upstream_only_hit_incompatibility(
	match_dspos1 const & m1,
	utility::vector1< HitPtrListCOP > const & upstream_only_hits,
	utility::vector1< std::list< Hit const * >::const_iterator > & upstream_only_hit_iterators,
	Size & last_upstream_only_geomcst_advanced,
  output::MatchProcessor const & processor
) const
{
	/// Determine if any of the secondary hits are incompatible and therefore
	/// could not produce a match.
	/// If we do find an incompatibility, then "continue" at the most-significant
	/// dimension.
	bool incompatible( false ), outstanding_list_increment( false );
	for ( Size ii = 2; ii <= n_geometric_constraints_; ++ii ) {
		Size ii_bp( m1.upstream_hits[ ii ].scaffold_build_id() );
		for ( Size jj = 1; jj < ii; ++jj ) {
			Size jj_bp( m1.upstream_hits[ jj ].scaffold_build_id() );

			/// We have already checked for compatibility between non-upstream-only matches.
			if ( ! geomcst_is_upstream_only_[ ii ] && ! geomcst_is_upstream_only_[ jj ] ) continue;

			if ( ! all_build_points_[ ii_bp ]->compatible( *all_build_points_[ jj_bp ] )) {
				incompatible = true;
			}
			if ( ! incompatible && ! upstream_builders_[ ii ]->compatible(
					fake_hit( m1.upstream_hits[ ii ] ), *all_build_points_[ ii_bp ],
					*upstream_builders_[ jj ], fake_hit( m1.upstream_hits[ jj ] ), *all_build_points_[ jj_bp ] ))  {
				incompatible = true;
			}

						//upstream downstream clash check
			//better here than later to alleviate combinatorics
			if( !incompatible && geomcst_is_upstream_only_[ii] && (jj == m1.originating_geom_cst_for_dspos) ){
				if( !processor.up_down_filt()->passes_hardsphere_filter( ii, jj, fake_hit(m1.upstream_hits[ii]), full_hit(m1) ) ) incompatible = true;
			}

			if( !incompatible && ( geomcst_is_upstream_only_[ii] || geomcst_is_upstream_only_[jj] ) ){
				if( !processor.up_coll_filt()->passes_hardsphere_filter( ii, jj, fake_hit( m1.upstream_hits[ii]), fake_hit(m1.upstream_hits[jj]) ) ) incompatible = true;
			}
			//clash check over

			if ( incompatible ) {
				/// Increment the most-significant geom-cst to alleviate the incompatibility
				/// either ii or jj represents an upstream-only geom-cst.
				/// if it's not ii, then we start our incrementing at jj.
				for ( Size kk = ( ! geomcst_is_upstream_only_[ ii ] ? jj : ii ); kk >= 1; --kk ) {
					if ( geomcst_is_upstream_only_[ kk ] ) {
						++upstream_only_hit_iterators[ kk ];
						last_upstream_only_geomcst_advanced = kk;
						if ( upstream_only_hit_iterators[ kk ] == upstream_only_hits[ kk ]->val().end() ) {
							outstanding_list_increment = true;
						} else {
							outstanding_list_increment = false;
							break;
						}
					}
				}
				//std::cout << "Incompatible at 1222 " << ii << " " << jj << std::endl;
				if ( outstanding_list_increment ) {
					// indicate to the calling function that we've visited all of the upstream-only hit combinations.
					last_upstream_only_geomcst_advanced = 0;
				}
				return true;
			} // if incompatible
		}
	}
	//std::cout << "Compatible at 1231" << std::endl;
	return false;
}

/// @details returns true if more upstream-only hit combinations remain,
/// and false if there are no upstream-only hit cominbations remaining
bool
Matcher::increment_upstream_only_hit_combination(
	utility::vector1< HitPtrListCOP > const & upstream_only_hits,
	Size starting_point,
	utility::vector1< std::list< Hit const * >::const_iterator > & upstream_only_hit_iterators,
	Size & last_upstream_only_geomcst_advanced
) const
{
	for ( Size ii = starting_point; ii >= 1; --ii ) {
		if ( geomcst_is_upstream_only_[ ii ] ) {
			++upstream_only_hit_iterators[ ii ];
			last_upstream_only_geomcst_advanced = ii;
			if ( upstream_only_hit_iterators[ ii ] != upstream_only_hits[ ii ]->val().end() ) {
				return true; // more upstream-only-hit-combinations remain
			}
		}
	}
	return false; // we did not find an upstream-only geom-cst that was not at the end of its hit list.
}


/// @details This needs a major refactoring.
void
Matcher::process_matches_main_loop_enumerating_all_hit_combos( output::MatchProcessor & processor ) const
{

	/// TEMP!  These variables need to be incremented by the process_matches_all_hit_combos_given_subsets routine
	//core::Size num_potential_matches(0), num_sent_to_proc(0),num_non_up_only_incompatible(0),num_up_only_incompatible(0), num_considered_muliple_origins(0), all_lex_states(0),num_ds_hit_incompatible(0), num_empty_uplist(0);
	MatcherOutputStats output_stats;

	utility::vector1< HitNeighborFinder > finders( hits_.size() );
	for ( Size ii = 1; ii <= hits_.size(); ++ii ) {
		if ( geomcst_is_upstream_only_[ ii ] ) continue; // don't create HitNeighborFinders for upstream-only geometric constraints.
		finders[ii].set_bounding_box( occ_space_bounding_box_ );
		finders[ii].set_xyz_bin_widths( euclidean_bin_widths_ );
		finders[ii].set_euler_bin_widths( euler_bin_widths_ );
		finders[ii].initialize();
		finders[ii].add_hits( hits_[ ii ] );
	}
	clock_t starttime = clock();
	utility::vector1< std::list< Hit const * > > hit_ccs = finders[ 1 ].connected_components();
	clock_t stoptime = clock();
	TR << "Found " << hit_ccs.size() << " connected component" << ( hit_ccs.size() != 1 ? "s" : "" ) << " in the hit neighbor graph in " << ((double) stoptime - starttime ) / CLOCKS_PER_SEC << " seconds." << std::endl;
	//TR << "CONNECTED COMPONENTS: " << hit_ccs.size() << std::endl;
	// -- this doesn't work -- don't turn it on -- #pragma omp parallel for
	for ( Size ii = 1; ii <= hit_ccs.size(); ++ii ) {
		//Size n_combos = hit_ccs[ ii ].size();
		//TR << "CC " << ii << " num neighbors: 1: " << hit_ccs[ ii ].size();
		utility::vector1< std::list< Hit const * > > ii_neighbor_hits( n_geometric_constraints_ );
		ii_neighbor_hits[ 1 ] = hit_ccs[ ii ]; // convenience: copy the list of hits in this CC.
		for ( Size jj = 2; jj <= hits_.size(); ++jj ) {
			if ( geomcst_is_upstream_only_[ jj ] ) continue; // no hits for upstream-only geometric constraints
			ii_neighbor_hits[ jj ] = finders[ jj ].neighbor_hits( hit_ccs[ ii ] );
		}

		Vector ii_euclidean_bin_widths( euclidean_bin_widths_ ), ii_euler_bin_widths( euler_bin_widths_ );
		if ( dynamic_grid_refinement_ ) {
			ii_neighbor_hits = refine_grid_and_subsample_for_hit_subsets( ii_euclidean_bin_widths, ii_euler_bin_widths, ii_neighbor_hits );
		}

		/// Create the hit hasher and then proceed to enumerate all match combos.
		HitHasher hit_hasher;
		hit_hasher.set_bounding_box( occ_space_bounding_box_ );
		hit_hasher.set_xyz_bin_widths( ii_euclidean_bin_widths );
		hit_hasher.set_euler_bin_widths( ii_euler_bin_widths );
		hit_hasher.set_nhits_per_match( n_geometric_constraints_ );
		hit_hasher.initialize();

		MatcherOutputStats ii_outstats = process_matches_all_hit_combos_for_hit_subsets( processor, hit_hasher, ii_neighbor_hits );
		output_stats += ii_outstats;

		/// Iterate across all 64 definitions of the origin to find all matches for this connected component.
	}

	TR << "Match enumeration statistics: ";
	TR << " num_potential_matches: " << output_stats.num_potential_matches <<
		";    all_lex_states " << output_stats.all_lex_states <<
		";   num_considered_muliple_origins " << output_stats.num_considered_muliple_origins <<
		";    num_non_up_only_incompatible " << output_stats.num_non_up_only_incompatible <<
		";    num_ds_hit_incompatible " << output_stats.num_ds_hit_incompatible <<
		";    num_up_only_incompatible " << output_stats.num_up_only_incompatible <<
		";    num_sent_to_proc " << output_stats.num_sent_to_proc << std::endl;
	TR << output_stats.num_empty_uplist << " empty uplists were observed." << std::endl;

}

utility::vector1< std::list< Hit const * > >
Matcher::refine_grid_and_subsample_for_hit_subsets(
	Vector & good_euclidean_bin_widths,
	Vector & good_euler_bin_widths,
	utility::vector1< std::list< Hit const * > > const & neighbor_hits
) const
{
	Size const acceptable( 500000 ); // "acceptible" must be larger than "take_it"
	Size const take_it( 300000 );
	Size const bare_minimum( 20 );

	Size n_combos = predict_n_matches_for_hit_subsets(	euclidean_bin_widths_, euler_bin_widths_, neighbor_hits, take_it );
	Size const initial_n_combos = n_combos;

	if ( n_combos > take_it ) {
		utility::vector1< std::list< Hit const * > > good_subsamples = subsample_hits( euclidean_bin_widths_, euler_bin_widths_, neighbor_hits );
		utility::vector1< std::list< Hit const * > > subsamples;

		Size n_combos = predict_n_matches_for_hit_subsets( euclidean_bin_widths_, euler_bin_widths_, good_subsamples, take_it );
		Size last_n_combos = n_combos;
		TR << "subsampling would produce " << n_combos << " matches down from " << initial_n_combos << " matches." << std::endl;

		Size count_refinement( 0 );
		bool found_acceptible( false );
		while ( n_combos > take_it ) {

			Vector test_euclidean_bin_widths( good_euclidean_bin_widths ), test_euler_bin_widths( good_euler_bin_widths );
			test_euclidean_bin_widths *= 0.75;
			test_euler_bin_widths *= 0.75;
			subsamples = subsample_hits( test_euclidean_bin_widths, test_euler_bin_widths, neighbor_hits ); // subsample the ORIGINAL set of hits
			n_combos = predict_n_matches_for_hit_subsets( test_euclidean_bin_widths, test_euler_bin_widths, good_subsamples, take_it );

			TR << "Grid refinement #" << count_refinement + 1 << " predicts " << n_combos << " matches." << std::endl;

			if ( n_combos > take_it ) {
				good_euclidean_bin_widths = test_euclidean_bin_widths;
				good_euler_bin_widths = test_euler_bin_widths;
				good_subsamples = subsamples;
				last_n_combos = n_combos;
				++count_refinement;
				if ( n_combos < acceptable	) found_acceptible = true;
			} else if ( ! found_acceptible && n_combos > bare_minimum ) {
				// Imagine a case where we had 2 Billion matches from the last round, and 20 matches in this round.
				// go ahead and take this refinement
				good_euclidean_bin_widths = test_euclidean_bin_widths;
				good_euler_bin_widths = test_euler_bin_widths;
				good_subsamples = subsamples; // expensive copy, avoid if we just made the grid smaller
				last_n_combos = n_combos;
				++count_refinement;
			} else if ( ! found_acceptible ) {
				// num combos < bare_minimum; try one last time with a slightly larger grid size:
				// NOTE we don't want to iterate through the outer while loop another time after we make this
				// refinement or we could get stuck in an infinite loop.
				Vector test_euclidean_bin_widths2( good_euclidean_bin_widths ), test_euler_bin_widths2( good_euler_bin_widths );
				test_euclidean_bin_widths2 *= 0.9;
				test_euler_bin_widths2 *= 0.9;
				subsamples = subsample_hits( test_euclidean_bin_widths2, test_euler_bin_widths2, neighbor_hits ); // subsample the ORIGINAL set of hits
				n_combos = predict_n_matches_for_hit_subsets( test_euclidean_bin_widths2, test_euler_bin_widths2, subsamples, take_it );

				TR << "(Partial) Grid refinement #" << count_refinement + 1 << " predicts " << n_combos << " matches." << std::endl;
				if ( n_combos > bare_minimum ) {
					good_euclidean_bin_widths = test_euclidean_bin_widths2;
					good_euler_bin_widths = test_euler_bin_widths2;
					good_subsamples = subsamples;
					last_n_combos = n_combos;
					++count_refinement;
				} else {
					TR << "Failed to refine grid to acceptible number of matches" << std::endl;
				}
				break; // we must exit this loop now.
			} else {
				TR << "Too fine a grid: " << n_combos << " matches predicted." << std::endl;
			}
		}
		TR << "Dynamic grid refinement predicts " << last_n_combos << " matches (down from an initial number of " << initial_n_combos << " matches) after " << count_refinement << " grid refinements" << std::endl;

		return good_subsamples;
	} else {
		return neighbor_hits;
	}

	return neighbor_hits;
}

utility::vector1< std::list< Hit const * > >
Matcher::subsample_hits(
	Vector const & euclidean_bin_widths,
	Vector const & euler_bin_widths,
	utility::vector1< std::list< Hit const * > > const & neighbor_hits
) const
{
	/// Take 1 hit per geomcst per halfbin
	utility::vector1< std::list< Hit const * > > subsamples( n_geometric_constraints_ );

	Vector half_width_euclid = 0.5 * euclidean_bin_widths;
	Vector half_width_euler  = 0.5 * euler_bin_widths;


	HitHasher hit_hasher;
	hit_hasher.set_bounding_box( occ_space_bounding_box_ );
	hit_hasher.set_xyz_bin_widths( half_width_euclid );
	hit_hasher.set_euler_bin_widths( half_width_euler );
	hit_hasher.set_nhits_per_match( n_geometric_constraints_ );
	hit_hasher.initialize();

	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		if ( geomcst_is_upstream_only_[ ii ] ) continue;
		for ( std::list< Hit const * >::const_iterator
				hit_iter = neighbor_hits[ ii ].begin(),
				hit_iter_end = neighbor_hits[ ii ].end();
				hit_iter != hit_iter_end; ++hit_iter ) {
			hit_hasher.insert_hit( 1, ii, *hit_iter );
		}
	}

	for ( HitHasher::HitHash::const_iterator halfbin_iter = hit_hasher.hit_hash_begin( 1 );
			halfbin_iter != hit_hasher.hit_hash_end( 1 ); ++halfbin_iter ) {
		HitHasher::MatchSet const & matches = halfbin_iter->second;
		for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
			if ( geomcst_is_upstream_only_[ ii ] ) continue;
			if ( matches[ ii ].begin() == matches[ ii ].end() ) continue;

			// pick one hit at random for each hit coming from a different upstream build point
			// first -- check if there are multiple upstream build points.
			Size build_point = 0;
			bool all_same( true );
			for ( std::list< Hit const * >::const_iterator
					hit_iter = matches[ ii ].begin(), hit_iter_end = matches[ ii ].end();
					hit_iter != hit_iter_end; ++hit_iter ) {
				if ( build_point == 0 ) {
					build_point = (*hit_iter)->first()[ 1 ];
				} else if ( build_point != (*hit_iter)->first()[ 1 ] ) {
					all_same = false;
					break;
				}
			}
			if ( all_same ) {
				Size len = matches[ ii ].size();
				if ( len == 1 ) {
					subsamples[ ii ].push_back( * matches[ ii ].begin() );
				} else {
					// pick a random number between 0 and len-1 and increment through the list of hits
					// that many times to pick a single random hit.
					std::list< Hit const * >::const_iterator hit_iter = matches[ ii ].begin();
					Size nsteps = static_cast< Size > ( numeric::random::rg().uniform() * len );
					for ( Size jj = 1; jj <= nsteps; ++jj ) ++hit_iter;
					subsamples[ ii ].push_back( * hit_iter );
				}
			} else {
				/// insert the hits into an STL map based on the build position; then grab one hit at random from each
				std::map< Size, std::list< Hit const * > > buildpos_hitmap;
				for ( std::list< Hit const * >::const_iterator
						hit_iter = matches[ ii ].begin(), hit_iter_end = matches[ ii ].end();
						hit_iter != hit_iter_end; ++hit_iter ) {
					buildpos_hitmap[ (*hit_iter)->first()[ 1 ] ].push_back( *hit_iter );
				}
				for ( std::map< Size, std::list< Hit const * > >::const_iterator
						builditer = buildpos_hitmap.begin(), builditer_end = buildpos_hitmap.end();
						builditer != builditer_end; ++builditer ) {
					Size len = builditer->second.size();
					if ( len == 1 ) {
						subsamples[ ii ].push_back( * builditer->second.begin() );
					} else {
						// pick a random number between 0 and len-1 and increment through the list of hits
						// that many times to pick a single random hit.
						std::list< Hit const * >::const_iterator hit_iter = builditer->second.begin();
						Size nsteps = static_cast< Size > ( numeric::random::rg().uniform() * len );
						for ( Size jj = 1; jj <= nsteps; ++jj ) ++hit_iter;
						subsamples[ ii ].push_back( * hit_iter );
					}

				}
			}
		}
	}
	return subsamples;
}

Matcher::Size
Matcher::predict_n_matches_for_hit_subsets(
	Vector const & euclidean_bin_widths,
	Vector const & euler_bin_widths,
	utility::vector1< std::list< Hit const * > > const & neighbor_hits,
	Size accuracy_threshold
) const
{
	assert( ! neighbor_hits[ 1 ].empty() );

	/// 1. First approximation: just add the log of the number of hits in all bins.
	/// If this is less than the log of the accuracy threshold, then don't bother computing a more accurate estimate.
	Real log_n_combos = std::log( (double) neighbor_hits[ 1 ].size() );
	TR << "predict_n_matches_for_hit_subsets with #hits for each geomcst: " << neighbor_hits[ 1 ].size();
	for ( Size ii = 2; ii <= n_geometric_constraints_; ++ii ) {
		if ( ! geomcst_is_upstream_only_[ ii ] ) {
			Size ii_nhits = neighbor_hits[ ii ].size();
			if ( ii_nhits == 0 ) { TR << std::endl; return 0; } // Quit now, no hits for this non-upstream-only geometric constraint
			log_n_combos += std::log( (double) ii_nhits ) ;
			TR << " " << ii_nhits;
		}
		//std::cout << "blah 1 " << std::endl;
	}
	TR << std::endl;
	if ( log_n_combos < std::log( (double) accuracy_threshold ) ) return static_cast< Size > ( exp( log_n_combos )) + 1; // good enough approximation


	MatchCounter match_counter;
	match_counter.set_bounding_box( occ_space_bounding_box_ );
	match_counter.set_xyz_bin_widths( euclidean_bin_widths );
	match_counter.set_euler_bin_widths( euler_bin_widths );

	Size count_non_upstream_only_hits( 0 );
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) if ( ! geomcst_is_upstream_only_[ ii ] ) ++count_non_upstream_only_hits;
	match_counter.set_n_geometric_constraints( count_non_upstream_only_hits );

	match_counter.initialize();

	count_non_upstream_only_hits = 0;
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		if ( ! geomcst_is_upstream_only_[ ii ] ) {
			++count_non_upstream_only_hits;
			match_counter.add_hits( count_non_upstream_only_hits, neighbor_hits[ ii ] );
		}
	}
	return match_counter.count_n_matches(); // the heavy lifting happens inside here
}

/// This loop iterates across all 64 origin definitions (loop "ii"), inserts the "neighbor hits" into
/// each of the hashes.  It then iterates across all the bins in the hash map (loop "iter"), and
/// then enumerates all combinations of hits (loop "lex") for the non-upstream-only geometric constraints.
///  For each combination of hits (a partial match if there are any upstream-only geometric constraints),
/// it then retrieves all upstream-only hits, and iterates across all cominbations of upstream-only
/// hits for this match (loop "true").  Then, for every fully-constructed match,
MatcherOutputStats
Matcher::process_matches_all_hit_combos_for_hit_subsets(
	output::MatchProcessor & processor,
	HitHasher & hit_hasher,
	utility::vector1< std::list< Hit const * > > const & neighbor_hits
) const
{
	/// TEMP: declare these variables here -- eventually, return them
	///core::Size num_potential_matches(0), num_sent_to_proc(0),
	///num_non_up_only_incompatible(0),num_up_only_incompatible(0),
	///num_considered_muliple_origins(0), all_lex_states(0),
	///num_ds_hit_incompatible(0), num_empty_uplist(0);
	MatcherOutputStats ostats;

	//MatchOutputTracker tracker; // APL new space saving algorithm does not require the match output tracker!
	match        m(     n_geometric_constraints_ );
	match_lite   mlite( n_geometric_constraints_ );
	match_dspos1 m1(    n_geometric_constraints_ ); // for finding upstream-only hits

	for ( Size ii = 1; ii <= 64; ++ii ) {

		/// Initialize the hit-hasher's hash map with the ii'th definition of the origin.
		/// Significant memory savings by deleting the ii'th hash map at the conclusion
		/// of the ii'th iteration, instead of having all 64 hash maps in memory at the same time.
		for ( Size jj = 1; jj <= n_geometric_constraints_; ++jj ) {
			if ( ! geomcst_is_upstream_only_[ jj ] ) {
				for ( std::list< Hit const * >::const_iterator iter = neighbor_hits[ jj ].begin(),
						iter_end = neighbor_hits[ jj ].end(); iter != iter_end; ++iter ) {
					hit_hasher.insert_hit( ii, jj, (*iter) );
					//hit_scaff_build_pts[ jj ][ iter->first()[ 1 ] ] = true;
				}
			}
			/// else { noop }, do not hash upstream-only hits; such hits are enumerated separately
		}

		for ( HitHasher::HitHash::const_iterator
				iter = hit_hasher.hit_hash_begin( ii ),
				iter_end = hit_hasher.hit_hash_end( ii );
				iter != iter_end; ++iter ) {

			/// MATCHES!
			HitHasher::MatchSet const & match_set( iter->second );

			utility::vector1< Size > n_hits_per_geomcst( n_geometric_constraints_ );
			utility::vector1< utility::vector1< Hit const * > > hit_vectors( n_geometric_constraints_ );
			// representative hits, if we select a subset of hits for each geometric constraint
			// enumerating all combinations of hits can be very very expensive.
			utility::vector1< utility::vector1< Size > > reps( n_geometric_constraints_ );
			/// hits from upstream-only downstream algorithms
			utility::vector1< HitPtrListCOP > upstream_only_hits( n_geometric_constraints_ );
			utility::vector1< std::list< Hit const * >::const_iterator > upstream_only_hit_iterators( n_geometric_constraints_ );

			/// First check that there is at least one hit per geometric constraint;
			/// go ahead and initialize the hit_vectors array at the same time.
			bool any_size_zero( false );
			for ( Size jj = 1; jj <= n_geometric_constraints_; ++jj ) {
				if ( ! geomcst_is_upstream_only_[ jj ] ) {
					n_hits_per_geomcst[ jj ] = match_set[ jj ].size();
					any_size_zero |= n_hits_per_geomcst[ jj ] == 0;
					hit_vectors[ jj ].resize( n_hits_per_geomcst[ jj ] );
					std::copy( match_set[ jj ].begin(), match_set[ jj ].end(), hit_vectors[ jj ].begin());
				} else {
					n_hits_per_geomcst[ jj ] = 1; // for the lex, indicate there's only a single value for geomcst jj
				}
			}
			if ( any_size_zero ) continue; // no matches possible in this voxel.

			select_hit_representatives( hit_vectors, n_hits_per_geomcst, reps );

			/// Prepare to iterate across all combinations of hits.
			utility::LexicographicalIterator lex( n_hits_per_geomcst );
			ostats.all_lex_states += lex.num_states_total();

			while ( ! lex.at_end() ) {
				++ostats.num_potential_matches;

				/// Assemble the (partial) match
				for ( Size jj = 1; jj <= n_geometric_constraints_; ++jj ) {
					if ( ! geomcst_is_upstream_only_[ jj ] ) {
						mlite[ jj ] =   hit_vectors[ jj ][ reps[jj][lex[ jj ]] ];
						m[ jj ]     = *(hit_vectors[ jj ][ reps[jj][lex[ jj ]] ] );
						m1.upstream_hits[ jj ].copy_hit( m[ jj ] );
					}
					else mlite[jj] = NULL; //null pointer to ensure output tracking works
				}

				/// if any of the non-upstream-only hits are incompatible, increment the lex
				/// and proceed to the next combination of hits
				if ( check_non_upstream_only_hit_incompatibility( m1, lex, processor ) ){
					ostats.num_non_up_only_incompatible++;
					continue;
				}
				if( check_potential_dsbuilder_incompatibility_ &&
						check_downstream_hit_incompatibility( m, lex ) ){
					ostats.num_ds_hit_incompatible++;
					continue;
				}

				numeric::geometry::hashing::Bin6D lower_halfbin_spanned( 1 );
				for ( Size jj = 1; jj <= n_geometric_constraints_; ++jj ) {
					if ( ! geomcst_is_upstream_only_[ jj ] ) {
						numeric::geometry::hashing::Bin6D hit_halfbin = hit_hasher.binner( ii ).halfbin6( m[ jj ].second() );
						for ( Size kk = 1; kk <= 6; ++kk ) {
							lower_halfbin_spanned[ kk ] *= hit_halfbin[ kk ]; // once it goes to zero, it stays at zero.
						}
					}
				}
				bool go_to_next_match( false );
				for ( Size jj = 1; jj <= 6; ++jj ) {
					if ( lower_halfbin_spanned[ jj ] != 0 ) {
						/// We've already seen this match or we'll see it again in a later context; don't process now.
						go_to_next_match = true;
						break;
					}
				}
				if ( go_to_next_match ) {
					++lex;
					continue;
				}

				/*if( tracker.match_has_been_output( mlite ) ){
					ostats.num_considered_muliple_origins++;
					++lex;
					continue;
				}
				tracker.note_output_match( mlite );*/
				/// Now descend into the upstream-only hits (and if there are none, output the singular
				/// match combination from the set of non-upstream-only hits)
				// "inner lex" traversed here -- we're enumerating all the combinations of upstream-only hits
				Size last_upstream_only_geomcst_advanced( 0 );
				while ( true ) {
					if ( last_upstream_only_geomcst_advanced != 0 ) {
						Size const lgca = last_upstream_only_geomcst_advanced; // brevity
						mlite[ lgca ] =  *upstream_only_hit_iterators[ lgca ];
						m[ lgca ]     = **upstream_only_hit_iterators[ lgca ];
						m1.upstream_hits[ lgca ].copy_hit( m[ lgca ] );
					}
					Size empty_hitlist_id( 0 );
					for ( Size jj = last_upstream_only_geomcst_advanced + 1; jj <= n_geometric_constraints_; ++jj ) {
						if ( geomcst_is_upstream_only_[ jj ] ) {
							upstream_only_hits[ jj ] = representative_downstream_algorithm_[ jj ]->
								hits_to_include_with_partial_match( m1 );
							if ( upstream_only_hits[ jj ] == 0 ) {
								empty_hitlist_id = jj;
								break;
							}
							upstream_only_hit_iterators[ jj ] = upstream_only_hits[ jj ]->val().begin();
							mlite[ jj ] =  *upstream_only_hit_iterators[ jj ];
							m[ jj ]     = **upstream_only_hit_iterators[ jj ];
							m1.upstream_hits[ jj ].copy_hit( m[ jj ] );
							last_upstream_only_geomcst_advanced = jj;
						}
					}
					if ( empty_hitlist_id != 0 ) {
						ostats.num_empty_uplist++;
						// advance the upstream_only iterators;
						if( ! increment_upstream_only_hit_combination( upstream_only_hits, empty_hitlist_id - 1,
								upstream_only_hit_iterators, last_upstream_only_geomcst_advanced ) ) break;
						continue;
					}

					if ( test_upstream_only_hit_incompatibility(
								m, upstream_only_hits, upstream_only_hit_iterators,
								last_upstream_only_geomcst_advanced, processor ) ) {
						// we have an incompatibility; break if we're at the end of the upstream-only hit combos
						ostats.num_up_only_incompatible++;
						if ( last_upstream_only_geomcst_advanced == 0 ) break;
						continue;
					}

					ostats.num_sent_to_proc++;
					processor.process_match( m );

					if ( ! increment_upstream_only_hit_combination(
							upstream_only_hits,
							n_geometric_constraints_,
							upstream_only_hit_iterators,
							last_upstream_only_geomcst_advanced )) {
						break;
					}
				} //while (true ), iteration over upstream only hits
				++lex;
			} // while ( ! lex.at_end() )
		} //iteration over hit hashes

		hit_hasher.clear_hash_map( ii );
	} //loop over 64 definitions of 6D hash origin

	TR << "Processed " << ostats.num_sent_to_proc << " matches" << std::endl;
	return ostats;
}

void
Matcher::process_matches_where_one_geomcst_defines_downstream_location(
	output::MatchProcessor & processor
) const
{
	typedef utility::fixedsizearray1< Size, 2 > Size2;
	typedef utility::OrderedTuple< Size2 > Size2Tuple;
	typedef std::map< Size2Tuple, Hit const * > UpstreamRotamerRepresentativeMap;
	typedef UpstreamRotamerRepresentativeMap::const_iterator UpRotRepMapConstIterator;

	utility::vector1< UpstreamRotamerRepresentativeMap > upstream_rotamer_representatives( n_geometric_constraints_ );

	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		if ( ! geomcst_is_upstream_only_[ ii ] ) {
			for ( std::list< Hit >::const_iterator hititer = hits_[ ii ].begin(),
					hititer_end = hits_[ ii ].end(); hititer != hititer_end; ++hititer ) {
				Size2 rotid;
				rotid[ 1 ] = hititer->scaffold_build_id();
				rotid[ 2 ] = hititer->upstream_conf_id();
				UpRotRepMapConstIterator representative = upstream_rotamer_representatives[ ii ].find( rotid );
				if ( representative == upstream_rotamer_representatives[ ii ].end() ) {
					upstream_rotamer_representatives[ ii ][ rotid ] = & (*hititer);
				}
			}
		}
	}

	utility::vector1< std::map< Size, bool > > hit_scaff_build_pts( n_geometric_constraints_ );


	TR << "Begining match enumeration:" << std::endl;
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		TR << "Geometric constraint " << ii << " produced " << hits_[ ii ].size() << " hits" << std::endl;
	}

	TR << "Begining examination of each of the 64 definitions of the 6D-hash-grid origin:" << std::endl;
	//TR << "process_matches_where_one_geomcst_defines_downstream_location" << std::endl;

	utility::vector1< HitNeighborFinder > finders( hits_.size() );
	for ( Size ii = 1; ii <= hits_.size(); ++ii ) {
		finders[ii].set_bounding_box( occ_space_bounding_box_ );
		finders[ii].set_xyz_bin_widths( euclidean_bin_widths_ );
		finders[ii].set_euler_bin_widths( euler_bin_widths_ );
		finders[ii].initialize();
		finders[ii].add_hits( hits_[ ii ] );
	}
	utility::vector1< std::list< Hit const * > > hit_ccs = finders[ 1 ].connected_components();
	TR << "CONNECTED COMPONENTS: " << hit_ccs.size() << std::endl;
	for ( Size ii = 1; ii <= hit_ccs.size(); ++ii ) {
		Size n_combos = hit_ccs[ ii ].size();
		std::cout << "CC " << ii << " num neighbors: 1: " << hit_ccs[ ii ].size();
		utility::vector1< std::list< Hit const * > > ii_neighbor_hits( n_geometric_constraints_ );
		ii_neighbor_hits[ 1 ] = hit_ccs[ ii ]; // convenience: copy the list of hits in this CC.
		for ( Size jj = 2; jj <= hits_.size(); ++jj ) {
			ii_neighbor_hits[ jj ] = finders[ jj ].neighbor_hits( hit_ccs[ ii ] );
			Size jj_nhits = ii_neighbor_hits[ jj ].size();
			if ( ! geomcst_is_upstream_only_[ jj ] ) n_combos *= jj_nhits;
			std::cout << " " << jj << "= " << jj_nhits ;
		}
		std::cout << " ncombos: " << n_combos << std::endl;

		HitHasher hit_hasher;
		hit_hasher.set_bounding_box( occ_space_bounding_box_ );
		hit_hasher.set_xyz_bin_widths( euclidean_bin_widths_ );
		hit_hasher.set_euler_bin_widths( euler_bin_widths_ );

		hit_hasher.set_nhits_per_match( n_geometric_constraints_ );
		hit_hasher.initialize();

		MatchOutputTracker tracker;
		//match          m( n_geometric_constraints_ );
		match_dspos1  m1( n_geometric_constraints_ );
		match_lite mlite( n_geometric_constraints_ );

		/// The match lite where each hit but one has been replaced by it's
		/// upstream-rotamer-representative pointer.  This match-lite prevents repitition
		/// of the match_dspos1 outputs.
		match_lite mliteprime( n_geometric_constraints_ );


		for ( Size jj = 1; jj <= 64; ++jj ) {
			for ( Size kk = 1; kk <= n_geometric_constraints_; ++kk ) {
				if ( ! geomcst_is_upstream_only_[ kk ] ) {
					for ( std::list< Hit >::const_iterator iter = hits_[ kk ].begin(), iter_end = hits_[ kk ].end();
							iter != iter_end; ++iter ) {
						hit_hasher.insert_hit( jj, kk, & (*iter) );
						hit_scaff_build_pts[ kk ][ iter->first()[ 1 ] ] = true;
					}
				} /// else noop, do not hash upstream-only hits; such hits are enumerated separately
			}

			for ( HitHasher::HitHash::const_iterator
					iter = hit_hasher.hit_hash_begin( jj ),
					iter_end = hit_hasher.hit_hash_end( jj );
					iter != iter_end; ++iter ) {

				/// MATCHES!
				HitHasher::MatchSet const & match_set( iter->second );

				utility::vector1< Size > n_hits_per_geomcst( n_geometric_constraints_ );
				utility::vector1< utility::vector1< Hit const * > > hit_vectors( n_geometric_constraints_ );
				//utility::vector1< utility::vector1< Size > > reps( n_geometric_constraints_ ); //representatives

				utility::vector1< Size > n_unique_upstream_rots( n_geometric_constraints_ );
				utility::vector1< utility::vector1< Hit const * > > unique_hit_representatives( n_geometric_constraints_ );
				/// hits from upstream-only downstream algorithms
				utility::vector1< HitPtrListCOP > upstream_only_hits( n_geometric_constraints_ );
				utility::vector1< std::list< Hit const * >::const_iterator > upstream_only_hit_iterators( n_geometric_constraints_ );

				bool any_size_zero( false );
				for ( Size kk = 1; kk <= n_geometric_constraints_; ++kk ) {
					if ( ! geomcst_is_upstream_only_[ kk ] ) {
						n_hits_per_geomcst[ kk ] = match_set[ kk ].size();
						if ( n_hits_per_geomcst[ kk ] == 0 ) {
							any_size_zero = true;
							break;
						}
						hit_vectors[ kk ].resize( n_hits_per_geomcst[ kk ] );
						std::copy( match_set[ kk ].begin(), match_set[ kk ].end(), hit_vectors[ kk ].begin());

						std::list< Hit const * > hit_representatives;
						for ( Size ll = 1; ll <= match_set[ kk ].size(); ++ll ) {
							Size2 upstream_rotamer;
							upstream_rotamer[ 1 ] = hit_vectors[ kk ][ ll ]->scaffold_build_id();
							upstream_rotamer[ 2 ] = hit_vectors[ kk ][ ll ]->upstream_conf_id();
							hit_representatives.push_back( upstream_rotamer_representatives[ kk ][ upstream_rotamer ] );
						}
						hit_representatives.sort();
						hit_representatives.unique();
						n_unique_upstream_rots[ kk ] = hit_representatives.size();
						unique_hit_representatives[ kk ].resize( n_unique_upstream_rots[ kk ] );
						std::copy( hit_representatives.begin(), hit_representatives.end(), unique_hit_representatives[ kk ].begin() );
					} else {
						n_unique_upstream_rots[ kk ] = 1; // to molify the lex
					}
				}
				if ( any_size_zero ) continue;

				utility::vector1< Size > n_hits_to_enumerate( n_geometric_constraints_ );

				for ( Size kk = 1; kk <= n_geometric_constraints_; ++kk ) {
					if ( ! output_match_dspos1_for_geomcst_[ kk ] ) continue;
					runtime_assert( ! geomcst_is_upstream_only_[ kk ] ); // THIS DOESN'T BELONG HERE. FIND IT A HOME!
					m1.originating_geom_cst_for_dspos = kk;

					/// For the geometric constraint being examined, construct the counts
					/// of unique matches where the location of the downstream partner is
					/// irrelevant for all other geometric constraints -- only the rotamer id of the
					/// upstream residue is relevant.
					for ( Size ll = 1; ll <= n_geometric_constraints_; ++ll ) {
						if ( kk == ll ) {
							n_hits_to_enumerate[ ll ] = n_hits_per_geomcst[ ll ];
						} else {
							n_hits_to_enumerate[ ll ] = n_unique_upstream_rots[ ll ];
						}
					}

					utility::LexicographicalIterator lex( n_hits_to_enumerate );

					while ( ! lex.at_end() ) {
						/// Assemble the match
						for ( Size ll = 1; ll <= n_geometric_constraints_; ++ll ) {
							if ( ! geomcst_is_upstream_only_[ ll ] ) {
								if ( ll == kk ) {
									mlite[ ll ] =   hit_vectors[ ll ][ lex[ll] ];
									//m[ ll ]     = *(hit_vectors[ ll ][ lex[ll] ] );
									m1.upstream_hits[ ll ] = upstream_hit( *(hit_vectors[ ll ][ lex[ll] ] ));
									m1.downstream_conf_id  = hit_vectors[ ll ][ lex[ll] ]->downstream_conf_id();
									m1.dspos               = hit_vectors[ ll ][ lex[ll] ]->second();
								} else {
									mlite[ ll ] =   unique_hit_representatives[ ll ][ lex[ll] ];
									//m[ ll ]     = *(unique_hit_representatives[ ll ][ lex[ll] ]);
									m1.upstream_hits[ ll ].copy_hit( *(unique_hit_representatives[ ll ][ lex[ll] ]) );
								}
							}
						}

						/// if any of the non-upstream-only hits are incompatible, increment the lex
						/// and proceed to the next combination of hits
						if ( check_non_upstream_only_hit_incompatibility( m1, lex, processor ) ) continue;

						/// Now descend into the upstream-only hits (and if there are none, output the singular
						/// match combination from the set of non-upstream-only hits)
						// "inner lex" traversed here -- we're enumerating all the combinations of upstream-only hits
						Size last_upstream_only_geomcst_advanced( 0 );
						while ( true ) {
							if ( last_upstream_only_geomcst_advanced != 0 ) {
								Size const lgca = last_upstream_only_geomcst_advanced; // brevity
								mlite[ lgca ] =  *upstream_only_hit_iterators[ lgca ];
								//m[ lgca ]     = **upstream_only_hit_iterators[ lgca ];
								m1.upstream_hits[ lgca ].copy_hit( *mlite[ lgca ] );
							}
							Size empty_hitlist_id( 0 );
							for ( Size kk = last_upstream_only_geomcst_advanced + 1; kk <= n_geometric_constraints_; ++kk ) {
								if ( geomcst_is_upstream_only_[ kk ] ) {
									upstream_only_hits[ kk ] = representative_downstream_algorithm_[ kk ]->
										hits_to_include_with_partial_match( m1 );
									if ( upstream_only_hits[ kk ] == 0 ) {
										empty_hitlist_id = kk;
										break;
									}
									upstream_only_hit_iterators[ kk ] = upstream_only_hits[ kk ]->val().begin();
									mlite[ kk ] =  *upstream_only_hit_iterators[ kk ];
									//m[ kk ]     = **upstream_only_hit_iterators[ kk ];
									m1.upstream_hits[ kk ].copy_hit( *mlite[ kk ] );
									last_upstream_only_geomcst_advanced = kk;
								}
							}
							if ( empty_hitlist_id != 0 ) {
								// advance the upstream_only iterators;
								increment_upstream_only_hit_combination( upstream_only_hits, empty_hitlist_id - 1,
									upstream_only_hit_iterators, last_upstream_only_geomcst_advanced );
							}

							if ( test_upstream_only_hit_incompatibility(
									m1, upstream_only_hits, upstream_only_hit_iterators,
										last_upstream_only_geomcst_advanced, processor ) ) {
								// we have an incompatibility; break if we're at the end of the upstream-only hit combos
								if ( last_upstream_only_geomcst_advanced == 0 ) break;
								continue;
							}

							/// If we've already seen this match for a previous value of jj (some alternate
							/// definition of the hasher origin) then avoid outputting it a second time.
							if ( ! tracker.match_has_been_output( mlite ) ) {
								/// Record that we have seen this combination of hits.
								tracker.note_output_match( mlite );
								processor.process_match( m1 );
							}

							if ( ! increment_upstream_only_hit_combination(
									upstream_only_hits,
									n_geometric_constraints_,
									upstream_only_hit_iterators,
									last_upstream_only_geomcst_advanced )) {
								break;
							}
						}
						++lex;
					}
				}
			}

			hit_hasher.clear_hash_map( jj );
			std::cout << "." << std::flush;
			if ( jj % 20 == 0 ) std::cout << std::endl;
		}
	}
	std::cout << std::endl;
}


void
Matcher::note_primary_change_to_geom_csts_hitlist( Size geomcst_id )
{
	if ( ! geom_cst_has_primary_modification_[ geomcst_id ] ) {
		geom_csts_with_primary_hitlist_modificiations_.push_back( geomcst_id );
		geom_cst_has_primary_modification_[ geomcst_id ] = true;
	}
}

/*
	core::pose::PoseOP upstream_pose_;
	core::pose::PoseOP downstream_pose_;

	utility::vector1< Size > pose_build_resids_;
	utility::vector1< ScaffoldBuildPointOP > all_build_points_;
	utility::vector1< utility::vector1< ScaffoldBuildPointOP > > per_constraint_build_points_;

	utility::vector1< std::list< Hit > > hits_;

	Size n_geometric_constraints_;
	utility::vector1< UpstreamBuilderOP   >                       upstream_builders_;
	utility::vector1< std::map< std::string, Size > >             build_set_id_for_restype_;
	utility::vector1< utility::vector1< DownstreamAlgorithmOP > > downstream_algorithms_;

	utility::vector1< DownstreamBuilderOP >   all_downstream_builders_;
	utility::vector1< DownstreamAlgorithmOP > all_downstream_algorithms_;

	BumpGridOP bb_grid_;
	utility::vector1< BumpGridOP > original_scaffold_residue_bump_grids_;

	OccupiedSpaceHashOP occ_space_hash_;

*/

}
}
