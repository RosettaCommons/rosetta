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
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov
/// @author Modified by Ryan Pavlovicz

// Unit headers
#include <protocols/simple_moves/PackPwatRotamersMover.hh>
#include <protocols/simple_moves/PackPwatRotamersMoverCreator.hh>

// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/prepack_pwat_rotamers.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/graph/Graph.hh>

// AUTO-REMOVED #include <utility/string_util.hh> // string_split

// option key includes
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <fstream>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <iostream>
#include <string>

namespace protocols {
namespace simple_moves {

using namespace core;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace scoring;

using basic::Warning;
using basic::t_warning;
static basic::Tracer TR("protocols.simple_moves.PackPwatRotamersMover");

bool sortDwell(const PointDwell &i, const PointDwell &j) { return i.dwell > j.dwell; }

// PackPwatRotamersMover

std::string
PackPwatRotamersMover::mover_name() {
	return PackPwatRotamersMoverCreator::mover_name();
}

std::string
PackPwatRotamersMoverCreator::keyname() const
{
	return PackPwatRotamersMoverCreator::mover_name();
}

protocols::moves::MoverOP
PackPwatRotamersMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new PackPwatRotamersMover );
}

std::string
PackPwatRotamersMoverCreator::mover_name()
{
	return "PackPwatRotamersMover";
}

void PackPwatRotamersMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PackPwatRotamersMover::provide_xml_schema( xsd );
}

PackPwatRotamersMover::PackPwatRotamersMover() :
	protocols::moves::Mover("PackPwatRotamersMover"),
	scorefxn_(0),
	task_(0),
	nloop_( option[ OptionKeys::packing::ndruns ].value() ),
	task_factory_(0),
	rotamer_sets_( RotamerSetsOP( new rotamer_set::RotamerSets ) ),
	ig_(0),
	pack_temp_(1.0),
	cluster_results_(true),
	limit_waters_(false),
	lkb_pwat_(false),
	exclude_exposed_(true),
	use_average_(true)
{}

PackPwatRotamersMover::PackPwatRotamersMover( std::string const & type_name ) :
	protocols::moves::Mover( type_name ),
	scorefxn_(0),
	task_(0),
	nloop_( option[ OptionKeys::packing::ndruns ].value() ),
	task_factory_(0),
	rotamer_sets_( RotamerSetsOP( new rotamer_set::RotamerSets ) ),
	ig_(0),
	pack_temp_(1.0),
	cluster_results_(true),
	limit_waters_(false),
	lkb_pwat_(false),
	exclude_exposed_(true),
	use_average_(true)
{}

// constructors with arguments
PackPwatRotamersMover::PackPwatRotamersMover(
	ScoreFunctionCOP scorefxn,
	PackerTaskCOP task,
	Size nloop
) :
	protocols::moves::Mover("PackPwatRotamersMover"),
	scorefxn_( scorefxn ),
	task_( task ),
	nloop_( nloop ),
	task_factory_(0),
	rotamer_sets_( RotamerSetsOP( new rotamer_set::RotamerSets ) ),
	ig_(0),
	pack_temp_(1.0),
	cluster_results_(true),
	limit_waters_(false),
	lkb_pwat_(false),
	exclude_exposed_(true),
	use_average_(true)
{}

PackPwatRotamersMover::~PackPwatRotamersMover(){}

PackPwatRotamersMover::PackPwatRotamersMover( PackPwatRotamersMover const & other ) :
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( other )
{
	scorefxn_ = other.score_function();
	task_ = other.task();
	nloop_ = other.nloop();
	task_factory_ = other.task_factory();
	rotamer_sets_ = RotamerSetsOP( new rotamer_set::RotamerSets );
	ig_ = 0;
	pack_temp_ = other.pack_temp();
	cluster_results_ = other.cluster_results();
	limit_waters_ = other.limit_waters();
	lkb_pwat_ = other.lkb_pwat();
	exclude_exposed_ = other.exclude_exposed();
	use_average_ = other.use_average();
}

// find dwell-weighted centroid
// return as PointDwell: xyz = centroid coords; dwell = sum of dwell times
PointDwell
PackPwatRotamersMover::centroid_MW_PointDwell( utility::vector1< PointDwell > point_group )
{
	Real weight_sum(0);
	Vector cent_MW(0);
	for ( Size i = 1; i <= point_group.size(); ++i ) {
		weight_sum += point_group[i].dwell;
		cent_MW += point_group[i].xyz * point_group[i].dwell;
	}
	cent_MW /= weight_sum;
	PointDwell centroid;
	centroid.xyz = cent_MW;
	centroid.dwell = weight_sum;
	return centroid;
}

PointDwell
PackPwatRotamersMover::most_visited( utility::vector1< PointDwell > point_group )
{
	for ( Size i=1; i<= point_group.size(); ++i ) {
		TR << " (before sort) point " << i << " of " << point_group.size() << " has dwell of " << point_group[i].dwell << std::endl;
	}

	std::sort(point_group.begin(),point_group.end(),sortDwell);

	for ( Size i=1; i<= point_group.size(); ++i ) {
		TR << " (after sort) point " << i << " of " << point_group.size() << " has dwell of " << point_group[i].dwell << std::endl;
	}

	return point_group[1];
}

// same as centroid_MW_PointDwell, but only returns a Vector (no dwell time information)
Vector
PackPwatRotamersMover::centroid_MW( utility::vector1< PointDwell > point_group )
{
	Real weight_sum(0);
	Vector cent_MW(0);
	for ( Size i = 1; i <= point_group.size(); ++i ) {
		weight_sum += point_group[i].dwell;
		cent_MW += point_group[i].xyz * point_group[i].dwell;
	}
	cent_MW /= weight_sum;
	return cent_MW;
}


utility::vector1< utility::vector1< PointDwell > >
PackPwatRotamersMover::cluster_rotset( utility::vector1< PointDwell > cutoff_set )
{
	Real clust_radius = option[ OptionKeys::corrections::water::cluster_radius ].value();
	Real clust_cutoff = option[ OptionKeys::corrections::water::cluster_cutoff ].value();

	utility::vector1< utility::vector1< PointDwell > > clusters;
	utility::vector1< PointDwell > newcluster; // assign first rotamer to first cluster
	newcluster.push_back(cutoff_set[1]);
	clusters.push_back(newcluster);
	TR << "Clustering group of " << cutoff_set.size() << " PWAT rotamers " << std::endl;
	bool added(false);
	for ( Size i = 2; i <= cutoff_set.size(); ++i ) {
		added = false;
		//    TR << i+1 << " --  distance from existing clusters " << std::endl;
		for ( Size j = 1; j <= clusters.size(); ++j ) {
			Vector clust_cent = centroid_MW(clusters[j]);
			//      TR << "centroid " << j+1 << " = " << clust_cent << std::endl;
			//      TR << "  distance of " << cutoff_set[i].xyz << " to centroid = " << clust_cent.distance(cutoff_set[i].xyz) << std::endl;
			if ( clust_cent.distance(cutoff_set[i].xyz) <= clust_radius ) {
				clusters[j].push_back(cutoff_set[i]);
				//        TR << "  added cutoff_set member " << i+1 << " to cluster " << j+1 << std::endl;
				added = true;
				break; // use break to prevent same rotamer from being added to multiple clusters
			}
		}
		if ( added ) continue;
		//    TR << "  rotamer not added to any existing clusters" << std::endl;
		utility::vector1< PointDwell > newcluster;
		newcluster.push_back(cutoff_set[i]);
		clusters.push_back(newcluster);
		//    TR << "  new cluster created for rotamer " << i+1 << std::endl;
	}

	TR << "Total number of initial clusters = " << clusters.size() << " created from " << cutoff_set.size() << " PWAT rotamers " << std::endl;

	// check cumulative dwell-time for clusters, removing those that don't meet clust_cutoff
	utility::vector1< utility::vector1< PointDwell > > final_clusters;
	Real cut_clusters(0);
	for ( Size i = 1; i <= clusters.size(); ++i ) {
		Real dwell_sum(0);
		for ( Size j = 1; j <= clusters[i].size(); ++j ) {
			dwell_sum += clusters[i][j].dwell;
		}
		if ( dwell_sum >= clust_cutoff ) {
			final_clusters.push_back(clusters[i]);
		} else {
			cut_clusters += 1;
		}
	}
	TR << "Total number of clusters removed = " << cut_clusters << std::endl;

	for ( Size k = 1; k <= final_clusters.size(); ++k ) {
		TR << "  Cluster " << k << " has " << final_clusters[k].size() << " members: " << std::endl;
		for ( Size kk = 1; kk <= final_clusters[k].size(); ++kk ) {
			TR << "    " << final_clusters[k][kk].xyz.to_string() << " " << final_clusters[k][kk].dwell << std::endl;
		}
	}
	return final_clusters;
}

void
PackPwatRotamersMover::apply( Pose & pose )
{
	using namespace core::conformation;

	Real dwell_cutoff = option[ OptionKeys::corrections::water::dwell_cutoff ].value();

	core::chemical::ResidueTypeSetCOP rsd_set( pose.residue_type_set_for_pose( core::chemical::FULL_ATOM_t ) );

	this->setup( pose );

	TR << "pwat packing temperature set to " << pack_temp_ << std::endl;

	utility::vector1< PointDwell > all_rot;

	Pose working_pose = pose;
	utility::vector0< int > rot_to_pack;

	this->run( working_pose, rot_to_pack, all_rot );

	// truncate rotamer list to those with dwell times within the defined cutoff
	utility::vector1< PointDwell > rot_cutoff;
	for ( Size x = 1; x <= all_rot.size(); ++x ) {
		if ( all_rot[x].dwell >= dwell_cutoff ) {
			rot_cutoff.push_back(all_rot[x]);
		}
	}

	std::sort(rot_cutoff.begin(),rot_cutoff.end(),sortDwell);
	TR << "length of rotamer set within cutoff of " << dwell_cutoff << " = " << rot_cutoff.size() << std::endl;
	for ( Size x = 1; x <= rot_cutoff.size(); ++x ) {
		TR << "   " << x << " " << rot_cutoff[x].xyz.to_string() << " " << rot_cutoff[x].dwell << std::endl;
	}

	// remove all point waters from the pose before adding cluster centroids
	// these are PWATs that were added with the WaterBox mover
	for ( Size ires=pose.total_residue(); ires >= 1; --ires ) {
		std::string resname = pose.residue(ires).name();
		if ( resname == "PWAT" || resname == "PWAT_V" ) {
			pose.conformation().delete_residue_slow( ires );
		}
	}

	Size waternum( 0 );
	Real watlim_scale = option[ OptionKeys::corrections::water::watlim_scale ].value();
	if ( limit_waters_ ) {
		Size nres = working_pose.total_residue();
		for ( int i=1; i<=(int)nres; ++i ) {
			//core::conformation::Residue const &res = pose.residue(i);
			if ( task_->pack_residue(i) || task_->design_residue(i) ) {
				waternum++;
			}
		}
	}

	if ( cluster_results_ && rot_cutoff.size() > 1 ) {
		// make cluster from PWAT rotamer cloud that falls within dwell time cutoff
		utility::vector1< utility::vector1< PointDwell > > clusters = cluster_rotset(rot_cutoff);

		// create dwell-weighted centroids from PWAT rotamer clusters
		// and add centroids to pose as rotatable waters
		utility::vector1< PointDwell > centroids;
		for ( Size x = 1; x <= clusters.size(); ++x ) {
			//      centroids.push_back(centroid_MW_PointDwell(clusters[x]));
			// try to place waters at most visited pwat position
			// instead of dwell-weighted centroid position
			//centroids.push_back( most_visited(clusters[x]) );

			// the members of the clusters should already be sorted by descending
			// dwell time, so just keep the first memeber of each cluster
			centroids.push_back(clusters[x][1]);
		}
		// sort final centroids by cumulative dwell time
		// in case a hard cutoff is used for number of waters
		// added to a pose -- this ensures the most relevant
		// waters are placed first
		std::sort(centroids.begin(),centroids.end(),sortDwell);

		// add centroids to pose as rotatable waters
		Size attach_to, ncent;
		if ( limit_waters_ ) {
			ncent = std::ceil( waternum *  watlim_scale );
			if ( ncent > centroids.size() ) {
				TR  << "watlim_scale allows for " << ncent << " waters, while only " << centroids.size() << " possible after cutoffs and/or clustering" << std::endl;
				TR  << "limiting number of waters to " << centroids.size() << std::endl;
				ncent = centroids.size();
			} else {
				TR << "limiting number of waters to " << ncent << std::endl;
			}
		} else {
			ncent = centroids.size();
		}
		for ( Size x = 1; x <= ncent; ++x ) {
			TR << "  centroid #" << x << ": " << centroids[x].xyz.to_string() << " with cumulative dwell time of " << centroids[x].dwell << std::endl;
			ResidueOP vrt_wat = ResidueFactory::create_residue( rsd_set->name_map("HOH") );
			Vector const OH1( vrt_wat->xyz("H1") - vrt_wat->xyz("O") );
			Vector const OH2( vrt_wat->xyz("H2") - vrt_wat->xyz("O") );
			ResidueOP new_res = ResidueOP( new Residue( *vrt_wat ) );
			new_res->set_xyz("O", centroids[x].xyz);
			new_res->set_xyz("H1", centroids[x].xyz+OH1);
			new_res->set_xyz("H2", centroids[x].xyz+OH2);
			attach_to = find_closest( pose, centroids[x].xyz );
			pose.append_residue_by_jump( *new_res, attach_to );
			pose.pdb_info()->set_resinfo( pose.total_residue(), pose.pdb_info()->chain(attach_to) , pose.total_residue(), ' ');
		}
	} else {
		// no clustering -- rotamers within cutoff converted directly to HOH
		Size attach_to, nrot;
		if ( limit_waters_ ) {
			nrot = std::ceil( waternum * watlim_scale );
			if ( nrot > rot_cutoff.size() ) {
				TR  << "watlim_scale allows for " << nrot << " waters, while only " << rot_cutoff.size() << " possible after cutoffs and/or clustering" << std::endl;
				TR  << "limiting number of waters to " << rot_cutoff.size() << std::endl;
				nrot = rot_cutoff.size();
			} else {
				TR << "limiting number of waters to " << nrot << std::endl;
			}
		} else {
			nrot = rot_cutoff.size();
		}
		for ( Size x = 1; x <= nrot; ++x ) {
			ResidueOP vrt_wat = ResidueFactory::create_residue( rsd_set->name_map("HOH") );
			Vector const OH1( vrt_wat->xyz("H1") - vrt_wat->xyz("O") );
			Vector const OH2( vrt_wat->xyz("H2") - vrt_wat->xyz("O") );
			ResidueOP new_res = ResidueOP( new Residue( *vrt_wat ) );
			new_res->set_xyz("O", rot_cutoff[x].xyz);
			new_res->set_xyz("H1", rot_cutoff[x].xyz+OH1);
			new_res->set_xyz("H2", rot_cutoff[x].xyz+OH2);
			attach_to = find_closest( pose, rot_cutoff[x].xyz );
			pose.append_residue_by_jump( *new_res, attach_to );
			pose.pdb_info()->set_resinfo( pose.total_residue(), pose.pdb_info()->chain(attach_to) , pose.total_residue(), ' ');
		}
	}
	//guaruntees proper scoring if this mover is used as a protocol (as in fixbb)
	(*scorefxn_)(pose);
}

Size
PackPwatRotamersMover::find_closest( Pose const & pose, Vector Ocoord )
{
	Real mindist = 1e6;
	Size attach_to = 0;
	for ( Size i(1); i <= pose.total_residue(); ++i ) {
		conformation::Residue const &res = pose.residue(i);
		if ( res.is_water() ) continue;
		Real dist = ( Ocoord-res.atom( res.nbr_atom() ).xyz() ).length();
		if ( dist < mindist ) {
			mindist = dist;
			attach_to = i;
		}
	}
	TR << "attaching new water to residue " << attach_to << " with distance = " << mindist << std::endl;
	return attach_to;
}

std::string
PackPwatRotamersMover::get_name() const {
	return PackPwatRotamersMoverCreator::mover_name();
}

void
PackPwatRotamersMover::show(std::ostream & output) const
{
	Mover::show(output);
	if ( score_function() != 0 ) {
		output << "Score function: " << score_function()->get_name() << std::endl;
	} else { output << "Score function: none" << std::endl; }
}

///@brief when the PackerTask was not generated locally, verify compatibility with pose
///@details the pose residue types must be equivalent to the ones used to generate the ResidueLevelTasks, because of the way that prevent_repacking and its associated flags work
bool
PackPwatRotamersMover::task_is_valid( Pose const & pose ) const
{
	if ( task_->total_residue() != pose.total_residue() ) return false;
	for ( Size i(1); i <= pose.total_residue(); ++i ) {
		chemical::ResidueTypeCOP r = pose.residue_type(i).get_self_ptr();
		if ( ! task_->residue_task(i).is_original_type( r ) ) return false;
	}
	return true;
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
PackPwatRotamersMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose
)
{
	if ( tag->hasOption("nloop") ) {
		nloop_ = tag->getOption<Size>("nloop",1);
		runtime_assert( nloop_ > 0 );
	}
	parse_score_function( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );

	pack_temp_ = tag->getOption<Real>("pack_temp",1.0);
	cluster_results_ = tag->getOption<bool>("cluster_results", true);
	limit_waters_ = tag->getOption<bool>("limit_waters", false);
	lkb_pwat_ = tag->getOption<bool>("lkb_pwat", false);
	exclude_exposed_ = tag->getOption<bool>("exclude_exposed", true);
	use_average_ = tag->getOption<bool>("use_average", true);
}

///@brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
void
PackPwatRotamersMover::parse_score_function(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	ScoreFunctionOP new_score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	if ( new_score_function == 0 ) return;
	score_function( new_score_function );
}

///@brief parse "task_operations" XML option (can be employed virtually by derived Packing movers)
void
PackPwatRotamersMover::parse_task_operations(
	TagCOP const tag,
	basic::datacache::DataMap const & datamap,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const &
)
{
	TaskFactoryOP new_task_factory( protocols::rosetta_scripts::parse_task_operations( tag, datamap ) );
	if ( new_task_factory == 0 ) return;
	task_factory( new_task_factory );
}

///@brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PackPwatRotamersMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new PackPwatRotamersMover );
}

///@brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PackPwatRotamersMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::simple_moves::PackPwatRotamersMover( *this ) );
}

///@brief update task to include all waters
core::pack::task::PackerTaskOP
PackPwatRotamersMover::update_task(
	Pose const & pose,
	core::pack::task::PackerTaskCOP & packer_task,
	bool force_include_current,
	bool only_water
) {
	utility::vector1< bool > residues_allowed_to_be_packed;
	Size nres = pose.total_residue();
	for ( Size i=1; i <= nres; ++i ) {
		std::string resname = pose.residue(i).name();
		if ( resname == "PWAT" || resname == "PWAT_V" || resname == "HOH" || resname == "HOH_V" ) {
			residues_allowed_to_be_packed.push_back(true);
		} else if ( only_water ) {
			residues_allowed_to_be_packed.push_back(false);
		} else if ( packer_task->pack_residue(i) ) {
			residues_allowed_to_be_packed.push_back(true);
		} else {
			residues_allowed_to_be_packed.push_back(false);
		}
	}

	core::pack::task::PackerTaskOP updated_task( new core::pack::task::PackerTask_( pose ) );
	updated_task->restrict_to_residues(residues_allowed_to_be_packed);
	updated_task->restrict_to_repacking();
	updated_task->initialize_from_command_line();
	// copy extra_chi information from original task to new task
	for ( Size i=1; i <= packer_task->total_residue(); ++i ) {
		if ( i > nres ) continue; // protection in case of smaller pose than what entered the WaterDdg mover
		updated_task->nonconst_residue_task( i ).or_ex1( packer_task->residue_task( i ).ex1() );
		updated_task->nonconst_residue_task( i ).or_ex2( packer_task->residue_task( i ).ex2() );
		updated_task->nonconst_residue_task( i ).or_ex3( packer_task->residue_task( i ).ex3() );
		updated_task->nonconst_residue_task( i ).or_ex4( packer_task->residue_task( i ).ex4() );
		updated_task->nonconst_residue_task( i ).or_ex1aro( packer_task->residue_task( i ).ex1aro() );
		updated_task->nonconst_residue_task( i ).or_ex2aro( packer_task->residue_task( i ).ex2aro() );
		updated_task->nonconst_residue_task( i ).or_ex1aro_exposed( packer_task->residue_task( i ).ex1aro_exposed() );
		updated_task->nonconst_residue_task( i ).or_ex2aro_exposed( packer_task->residue_task( i ).ex2aro_exposed() );
		updated_task->nonconst_residue_task( i ).and_extrachi_cutoff( packer_task->residue_task( i ).extrachi_cutoff() );
		if ( force_include_current ) {
			updated_task->nonconst_residue_task( i ).or_include_current( true );
		} else {
			updated_task->nonconst_residue_task( i ).or_include_current( packer_task->residue_task( i ).include_current() );
		}
	}

	return updated_task;
}

///@brief get rotamers, energies. Also performs lazy initialization of ScoreFunction, PackerTask.
void PackPwatRotamersMover::setup( Pose & pose )
{
	using namespace core::pack::interaction_graph;

	// jec update_residue_neighbors() required to update EnergyGraph (ensures graph_state == GOOD) when calling Interface.cc
	pose.update_residue_neighbors();
	// guarantee of valid ScoreFunction and PackerTask postponed until now
	if ( scorefxn_ == 0 ) {
		Warning() << "undefined ScoreFunction -- creating a default one" << std::endl;
		scorefxn_ = get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
	}

	// if present, task_factory_ always overrides/regenerates task_
	if ( task_factory_ != 0 ) {
		task_ = task_factory_->create_task_and_apply_taskoperations( pose );
	} else if ( task_ == 0 ) {
		Warning() << "undefined PackerTask -- creating a default one" << std::endl;
		task_ = TaskFactory::create_packer_task( pose );
	} else runtime_assert( task_is_valid( pose ) );
	// in case PackerTask was not generated locally, verify compatibility with pose

	note_packertask_settings( pose );

	// if using lkbridge / backbone pwat solvation, first get rotamer sets
	if ( lkb_pwat_ ) {
		SetofSets new_pwat_rotsets;
		lkb_pwat_rotamers_setup( pose, *scorefxn_, task_, rotamer_sets_, new_pwat_rotsets, exclude_exposed_, use_average_ );

		// remove existing bb_pwat residue, and add new pwat residues, one for each new rotset
		rotamer_sets_->build_pwat_rotsets( pose, new_pwat_rotsets );

		// update pose and task to reflect different number of waters
		pose.update_residue_neighbors();
		PackerTaskOP updated_task( update_task( pose, task_, false, false ) );

		rotamer_sets_->set_task( updated_task );
		rotamer_sets_->prepare_sets_for_packing( pose, *scorefxn_ );

		TR << "number of final rotamers = " << rotamer_sets_->nrotamers() << std::endl;
		TR << "number of final moltenres = " << rotamer_sets_->nmoltenres() << std::endl;

		scorefxn_->setup_for_packing( pose, updated_task->repacking_residues(), updated_task->designing_residues() );
		utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, *scorefxn_, updated_task );

		time_t const anneal_start = time(NULL);
		ig_ = InteractionGraphFactory::create_and_initialize_annealing_graph(*updated_task, *rotamer_sets_, pose, *scorefxn_, packer_neighbor_graph );
		TR.Trace << "time spent building interaction graph with lkb pwat = " << time(NULL) - anneal_start << " seconds -- runtime" << std::endl;
	} else {
		// normal pwat packing using solvation sites built from database
		pack_rotamers_setup( pose, *scorefxn_, task_, rotamer_sets_, ig_ );
	}
}

core::PackerEnergy PackPwatRotamersMover::run( Pose & pose, utility::vector0< int > rot_to_pack, utility::vector1< PointDwell > & all_rot ) const
{
	return pack_pwat_rotamers_run( pose, rotamer_sets_, ig_, rot_to_pack, all_rot );
}

///@brief note PackerTask's packable and designable residues as string info
void PackPwatRotamersMover::note_packertask_settings( Pose const & pose )
{
	std::ostringstream packable, designable;
	packable << "REMARK PackingRes";
	designable << "REMARK DesignRes";

	for ( Size i(1), end( task_->total_residue() ); i <= end; ++i ) {
		ResidueLevelTask const & rtask( task_->residue_task(i) );
		if ( rtask.being_designed() ) {
			designable << ", ";
			if ( pose.pdb_info() ) {
				designable << pose.pdb_info()->number(i) << " " << pose.pdb_info()->chain(i);
			} else {
				designable << i;
			}
		} else if ( rtask.being_packed() ) {
			packable << ", ";
			if ( pose.pdb_info() ) {
				packable << pose.pdb_info()->number(i) << " " << pose.pdb_info()->chain(i);
			} else {
				packable << i;
			}
		}
	}
	info().clear();
	info().push_back( packable.str() );
	info().push_back( designable.str() );
}

// setters
void PackPwatRotamersMover::score_function( ScoreFunctionCOP sf )
{
	runtime_assert( sf );
	scorefxn_ = sf;
}

void PackPwatRotamersMover::task( task::PackerTaskCOP t ) { task_ = t; }

void PackPwatRotamersMover::task_factory( TaskFactoryCOP tf )
{
	runtime_assert( tf );
	task_factory_ = tf;
}

void PackPwatRotamersMover::nloop( Size nloop_in ) { nloop_ = nloop_in; }

utility::tag::XMLSchemaComplexTypeGeneratorOP
PackPwatRotamersMover::complex_type_generator_for_pack_pwat_rotamers_mover( utility::tag::XMLSchemaDefinition & )
{

	using namespace utility::tag;
	AttributeList attributes;

	attributes + XMLSchemaAttribute::attribute_w_default(  "nloop", xsct_non_negative_integer, "Equivalent to \"-ndruns\"."
		"Number of complete packing runs before an output (best score) is produced.",  "1"  );

	attributes + XMLSchemaAttribute::attribute_w_default(
		"cluster_results", xsct_rosetta_bool,
		"If true, cluster after point-water packing to reduce number of three-point waters (HOH) added to pose",
		"true");

	attributes + XMLSchemaAttribute::attribute_w_default(
		"limit_waters", xsct_rosetta_bool,
		"If true, limit number of final three-point waters (HOH) to number of packable protein residues",
		"false");

	attributes + XMLSchemaAttribute::attribute_w_default(
		"lkb_pwat", xsct_rosetta_bool,
		"If true, pwat rotamer cloud will be based on intesection of all possible lkball sites",
		"false");

	attributes + XMLSchemaAttribute::attribute_w_default(
		"exclude_exposed", xsct_rosetta_bool,
		"If true, possible hydration sites with fewer than 16 CB (or CA for GLY) within 10A will be excluded",
		"true");

	attributes + XMLSchemaAttribute::attribute_w_default(
		"use_average", xsct_rosetta_bool,
		"If using lkball intersection method, use average position of intersecting sites prior to subsequent clustering",
		"true");

	rosetta_scripts::attributes_for_parse_score_function( attributes );
	rosetta_scripts::attributes_for_parse_task_operations( attributes );

	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );
	ct_gen->complex_type_naming_func( & moves::complex_type_name_for_mover )
		.add_attributes( attributes )
		.add_optional_name_attribute();
	return ct_gen;
}

void PackPwatRotamersMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = complex_type_generator_for_pack_pwat_rotamers_mover( xsd );
	ct_gen->element_name( mover_name() )
		.description( "Repacks sidechains with user-supplied options, including TaskOperations." )
		.write_complex_type_to_schema( xsd );
}

// accessors
ScoreFunctionCOP PackPwatRotamersMover::score_function() const { return scorefxn_; }
PackerTaskCOP PackPwatRotamersMover::task() const { return task_; }
TaskFactoryCOP PackPwatRotamersMover::task_factory() const { return task_factory_; }
rotamer_set::RotamerSetsCOP PackPwatRotamersMover::rotamer_sets() const { return rotamer_sets_; }
interaction_graph::AnnealableGraphBaseCOP PackPwatRotamersMover::ig() const { return ig_; }
Real PackPwatRotamersMover::pack_temp() const { return pack_temp_; }
bool PackPwatRotamersMover::cluster_results() const { return cluster_results_; }
bool PackPwatRotamersMover::limit_waters() const { return limit_waters_; }
bool PackPwatRotamersMover::lkb_pwat() const { return lkb_pwat_; }
bool PackPwatRotamersMover::exclude_exposed() const { return exclude_exposed_; }
bool PackPwatRotamersMover::use_average() const { return use_average_; }

std::ostream &operator<< (std::ostream &os, PackPwatRotamersMover const &mover)
{
	mover.show(os);
	return os;
}

} // moves
} // protocols


