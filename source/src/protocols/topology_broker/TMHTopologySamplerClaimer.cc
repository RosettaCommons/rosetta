// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/topology_broker/TMHTopologySamplerClaimer.hh
/// @brief source file for TMHTopologySamplerClaimer protocol
/// @detailed implementation of sampling protocol that treats transmembrane helices as rigid bodies and moves them around to improve
/// 	sampling of membrane protein topologies
///
/// @author Stephanie H. DeLuca (stephanie.h.deluca@vanderbilt.edu)

// Unit Headers
#include <protocols/topology_broker/TMHTopologySamplerClaimer.hh>
#include <protocols/topology_broker/RigidBodyRandomTMHMover.hh>
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/FragmentClaimer.hh>
#include <core/pose/util.hh>
// Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>
#include <core/scoring/MembraneTopology.hh>
#include <core/scoring/MembranePotential.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/BBClaim.hh>
#include <protocols/topology_broker/claims/CutClaim.hh>
#include <protocols/topology_broker/claims/JumpClaim.hh>
#include <protocols/topology_broker/claims/LegacyRootClaim.hh>

// Utility headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>

#include <core/types.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RollMover.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>

static basic::Tracer tr("protocols.topo_broker.TMHTopologySampler");

namespace protocols{
namespace topology_broker{

TMHTopologySamplerClaimer::TMHTopologySamplerClaimer()
	: rotation_mag_(0),
	  translation_mag_(0),
	  rb_mover_stage1_weight_(0),
	  nres_(0),
	  njumps_(0),
	  jump_array_(2,0),
	  topology_root_res_(0),
	  tmhelix_(0)
{
	set_defaults();
}

TMHTopologySamplerClaimer::TMHTopologySamplerClaimer(TopologyBrokerOP /*broker*/)
	: rotation_mag_(0),
	  translation_mag_(0),
	  rb_mover_stage1_weight_(0),
	  nres_(0),
	  njumps_(0),
	  jump_array_(2,0),
	  topology_root_res_(0),
	  tmhelix_(0)
{
	set_defaults();
}

TMHTopologySamplerClaimer::~TMHTopologySamplerClaimer() {}

using namespace basic::options;
using namespace basic::options::OptionKeys;

//@brief register cmd-line options in option system ( call before core::init )
void
TMHTopologySamplerClaimer::register_options()
{
	option.add_relevant(basic::options::OptionKeys::rigid::rotation);
	option.add_relevant(basic::options::OptionKeys::rigid::translation);
	option.add_relevant(basic::options::OptionKeys::broker::rb_mover_stage1_weight);
}

///@brief read tag from topology broker file (setup.tpb)
bool
TMHTopologySamplerClaimer::read_tag( std::string tag, std::istream& is)
{
	Parent::read_tag(tag,is);
	return true;
}

core::scoring::MembraneTopologyOP
TMHTopologySamplerClaimer::get_membrane_topology(core::pose::Pose& pose)
{
	//get the membrane_topology
	core::scoring::MembraneTopologyOP membrane_topology;
	if (pose.data().has( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY ) )
	{
		membrane_topology = static_cast< core::scoring::MembraneTopology * >
		( pose.data().get_ptr( core::pose::datacache::CacheableDataType::MEMBRANE_TOPOLOGY) () );
	}else{
		utility_exit_with_message("Must have MembraneTopology!");
	}
	return membrane_topology;
}

core::scoring::MembraneEmbed
TMHTopologySamplerClaimer::get_membrane_embed(core::pose::Pose& pose)
{
	//membrane embedding not computed until ScoringManager::get_instance() called
	core::scoring::MembranePotential membrane_potential(core::scoring::ScoringManager::get_instance()->get_MembranePotential());
	membrane_potential.compute_membrane_embedding(pose);
	core::scoring::MembraneEmbed membrane_embed = core::scoring::nonconst_MembraneEmbed_from_pose(pose);

	return membrane_embed;
}

//getter function to retrieve this claimer's current_pose (to get after moving spans)
core::pose::PoseOP
TMHTopologySamplerClaimer::get_pose_from_claimer()
{
	return new core::pose::Pose(current_pose_);
}

///@brief the broker checks if the claimer builds its own fold tree to figure out if needs to build one itself
bool
TMHTopologySamplerClaimer::claimer_builds_own_fold_tree()
{
	return true;
}

//@brief called by constructor ---  calls all set_default_XXX methods
void
TMHTopologySamplerClaimer::set_defaults()
{
	set_label("DEFAULT");
	if(option[basic::options::OptionKeys::rigid::rotation].value())
	{
		rotation_mag_ = option[basic::options::OptionKeys::rigid::rotation].value();
	}
	else{
		rotation_mag_ = 3.0;
	}

	if(option[basic::options::OptionKeys::rigid::translation].value())
	{
		translation_mag_ = option[basic::options::OptionKeys::rigid::translation].value();
	}
	else{
		translation_mag_ = 8.0;
	}

	if(option[basic::options::OptionKeys::broker::rb_mover_stage1_weight].value())
	{
		rb_mover_stage1_weight_ = option[basic::options::OptionKeys::broker::rb_mover_stage1_weight].value();
	}else{
		rb_mover_stage1_weight_ = 5.0;
	}
}

///@brief get the pose from the boker and set it as this object's pose
void
TMHTopologySamplerClaimer::set_pose_from_broker(core::pose::Pose& pose)
{
	current_pose_ = pose;
}

///@brief read in the pose's spans via the MembraneTopology stored in the pose, determine jumps, TMHs, loops, etc.
void
TMHTopologySamplerClaimer::pre_process(core::pose::Pose& pose)
{
	core::scoring::MembraneTopologyOP membrane_topology = get_membrane_topology(pose);

	assert(membrane_topology);
	tmhelix_ = membrane_topology->tmhelix();
	core::Size nspan = tmhelix_;
	core::Size num_cut_loops = nspan;
	njumps_ = nspan;
	nres_ = pose.total_residue();

	if(tr.Trace.visible())
	{
		tr.Trace << "fold tree at beginning of TMHTopologySampler::build_fold_tree()\n" << pose.fold_tree();
	}
	tr.Info << "will set up fold_tree(nres) with " << nres_ << " residues, " << nspan << " spans, and " << njumps_ << " jump(s) in a total of " << num_cut_loops << " loops to cut\n";

	// exit if spanfile not formatted correctly
	if(nspan==0)
	{
		utility_exit_with_message("bad format for spanfile total_tmhelix==0");
	}

	//exit if nspan is not more than njumps and if # of cuttable loops != # jumps
	assert(nspan>=njumps_);
	assert(njumps_==num_cut_loops);

	//put all residues in spans in a one-dimensional array that is nres long, and 0=not in span, 1-3=is in span
	ObjexxFCL::FArray1D_int span(nres_,0);
	for(int i=1;i<=nres_;++i)
	{
		for(core::Size j=1;j<=nspan;++j)
		{
			if(static_cast<core::Size>(i) >= membrane_topology->span_begin(j) && static_cast<core::Size>(i) <= membrane_topology->span_end(j))
			{
				span(i)=j;
			}
		}
	}

	//Get info to make sure membrane span information is being read in correctly
	if(tr.Trace.visible())
	{
		tr.Trace << "size of span array:  " << span.size() << "\tnumber of residues:  " << nres_ << std::endl;
		for(core::Size i=1;i<=span.size();++i)
		{
			tr.Trace << "residue:  " << i << " span:  " << span(i) << std::endl;
		}
	}

	core::Size topology_root_jump_res = 0;

	topology_root_res_ = pose.fold_tree().root();
	topology_root_jump_res = pose.fold_tree().begin()->stop();
	tr.Info << "topology_root:  " << topology_root_res_ << std::endl;
	tr.Info << "topology_root_jump:  " << topology_root_jump_res << std::endl;

	//Make an array of start and end points of loops where we can cut. Needed by FoldTree::random_tree_from_jump_points();
	// initialize some variables for setting up loops array
	if(tr.Trace.visible())
	{
		tr.Trace << "fold tree before setting jumps\n" << pose.fold_tree();
		tr.Trace << "fold_tree_root:  " << pose.fold_tree().root() << std::endl;
	}

	if(tr.Trace.visible())
	{
		for(core::Size atom = 1; atom<=pose.residue(topology_root_res_).natoms(); ++atom)
		{
			tr.Trace << "atom in root residue:  " << pose.residue(topology_root_res_).atom_name(atom) << std::endl;
		}
	}

	core::Size span_index = 1;
	ObjexxFCL::FArray1D_int previous_span_begin((nspan-1),0);
	ObjexxFCL::FArray1D_int previous_span_end((nspan-1),0);
	ObjexxFCL::FArray1D_int span_begin((nspan-1),0);
	ObjexxFCL::FArray1D_int span_end((nspan-1),0);
	ObjexxFCL::FArray1D_int jump_begin(njumps_,0);
	ObjexxFCL::FArray1D_int jump_end(njumps_,0);
	ObjexxFCL::FArray1D_int loop_begin(num_cut_loops,0);
	ObjexxFCL::FArray1D_int loop_end(num_cut_loops,0);

	jump_array_ = ObjexxFCL::FArray2D_int(2,njumps_);

	//set up cuts array
	core::Size ncuts(0);
	ncuts = njumps_;
	cuts_ = ObjexxFCL::FArray1D_int(ncuts,0);

	// set up cut loops array
	for(span_index = 1; span_index <= nspan-1;++span_index)
	{
		//need to know the beginning and end of the TMs to determine where the loops are
		previous_span_begin(span_index) = membrane_topology->span_begin(span_index);
		previous_span_end(span_index) = membrane_topology->span_end(span_index);
		span_begin(span_index) = membrane_topology->span_begin(span_index+1);
		span_end(span_index) = membrane_topology->span_end(span_index+1);

		tr.Info << "span_index:  " << span_index << " previous_span_begin:  " << previous_span_begin(span_index) << " previous_span_end: "
				<< previous_span_end(span_index) << std::endl;
		tr.Info << "span_index:  " << span_index << " span_begin:  " << span_begin(span_index) << " span_end: "
				<< span_end(span_index) << std::endl;

		loop_begin(span_index) = previous_span_end(span_index) + 1;
		loop_end(span_index) = span_begin(span_index) - 1;

		//if the predicted loop (that is, not span) is shorter than 3 res
		if(loop_end(span_index) - loop_begin(span_index) < 2)
		{
			loop_begin(span_index) -= 1;
			loop_end(span_index) += 1;
		}

		assert(loop_begin(span_index) != 0);
		assert(loop_end(span_index) != 0);

		//put a cut at a random number in this loop
		core::Size random_number = numeric::random::RG.random_range(loop_begin(span_index),loop_end(span_index));
		cuts_(span_index) = random_number;

		if(tr.Debug.visible())
		{
			tr.Debug << "LOOP:  " << loop_begin(span_index) << " " << loop_end(span_index) << std::endl;
			tr.Debug << "random number:  " << random_number << std::endl;
		}

		tr.Info << "span_index:  " << span_index << " cuts(span_index):  " << cuts_(span_index) << std::endl;

		//make a 2D array of jumps, where it looks like:
		//res1 res2
		//res1 res2
		//...for njumps.  The res1 and res2 are the two residues that are connected by the jump edge
		jump_end(span_index) = topology_root_res_;
		jump_begin(span_index) = core::pose::residue_center_of_mass(pose,previous_span_begin(span_index),previous_span_end(span_index));
		//special case, last span
		core::Size last_span_jump_end = 0;
		last_span_jump_end = core::pose::residue_center_of_mass(pose,membrane_topology->span_begin(nspan),membrane_topology->span_end(nspan));

		for(core::Size jump_array_index = 1;jump_array_index <= njumps_; ++jump_array_index)
		{
			if(jump_array_index==njumps_)
			{
				jump_array_(1,jump_array_index) = last_span_jump_end;
				jump_array_(2,jump_array_index) = topology_root_res_;
			}else{
				jump_array_(1,jump_array_index) = jump_begin(jump_array_index);
				jump_array_(2,jump_array_index) = jump_end(jump_array_index);
			}
		}
	}//for span_index

	assert(core::Size(nres_)==pose.total_residue());
	assert(pose.fold_tree().root()==nres_);

	//make the last cut the last residue before the root residue
	cuts_(ncuts) = nres_-1;
	for(core::Size cut_index = 1; cut_index <= ncuts; ++cut_index)
	{
		tr.Info << "cuts array:  " << cut_index << " " << cuts_(cut_index) << std::endl;
	}

	//Make extended chain an idealized helix
	set_pose_torsions(pose);
}

//@brief this claimer builds its own radial fold tree based on read-in spanfile
void
TMHTopologySamplerClaimer::build_fold_tree(core::pose::Pose& pose, core::kinematics::FoldTree& fold_tree_in)
{
	// Generate the fold tree
	tr.Debug << fold_tree_in << std::endl;
	core::kinematics::FoldTree new_fold_tree = fold_tree_in;
	new_fold_tree.tree_from_jumps_and_cuts(nres_,njumps_,jump_array_,cuts_,topology_root_res_);
	pose.fold_tree(new_fold_tree);

	if(tr.Debug.visible())
	{
		tr.Debug << "TMHTopologySampler finished making fold tree" << std::endl;
		tr.Debug << new_fold_tree << std::endl;
		tr.Debug << "pose has " << pose.total_residue() << " residues.  Fold_tree_begin:  " << pose.fold_tree().begin()->start()
				<< " " << pose.fold_tree().begin()->stop() << "\tfold_tree_end:  " << pose.fold_tree().end()->start() <<
				" " << pose.fold_tree().end()->stop() << std::endl;
		tr.Debug << pose.annotated_sequence() << std::endl;
	}

	for(Size jump_num = 1; jump_num <= new_fold_tree.num_jump(); ++jump_num){
		tr.Info << "final_jump_array:  " << jump_array_(1,jump_num) << " "
			<< jump_array_(2,jump_num) << std::endl;
		new_fold_tree.set_jump_atoms(jump_num, "X", "N" );
	}

	// output fold tree
	if(new_fold_tree.empty())
	{
		tr << "Warning!  fold_tree is empty!!" << std::endl;
	}
	else{
		tr << new_fold_tree;
	}
}
//getter function to retrieve this claimer's fold tree for the broker
core::kinematics::FoldTreeOP
TMHTopologySamplerClaimer::get_fold_tree(core::pose::Pose& pose)
{
	return new core::kinematics::FoldTree(pose.fold_tree());
}

///@brief generate DoF Claims to be read in by broker
void
TMHTopologySamplerClaimer::generate_claims(claims::DofClaims &dof_claims)
{
	for(Size jump_num = 1; jump_num <= njumps_; ++jump_num)
	{
		dof_claims.push_back(new claims::JumpClaim(this,jump_array_(1,jump_num),jump_array_(2,jump_num),claims::DofClaim::CAN_INIT));
	}
	for(int i=1;i<=nres_;++i)
	{
		for(core::Size cut_index = 1; cut_index <= cuts_.size(); ++cut_index)
		{
			if(cuts_(cut_index)==i)
			{
					dof_claims.push_back(new claims::CutClaim(this,std::make_pair(TopologyClaimer::label(),i),claims::DofClaim::CAN_INIT));
			}
			if(i == topology_root_res_)
			{
				dof_claims.push_back(new claims::BBClaim(this,i,claims::DofClaim::CAN_INIT));
				dof_claims.push_back(new claims::LegacyRootClaim(this,i,claims::DofClaim::CAN_INIT));
			}
		}
	}
	current_pose_ = broker().current_pose();
}

///@brief make move_map and add the DoF claims from generate_claims() to the movemap.  Now we can move certain parts with certain DOFs.
void
TMHTopologySamplerClaimer::initialize_dofs( core::pose::Pose& pose, claims::DofClaims const& init_dofs, claims::DofClaims& /*failed_to_init */)
{
	core::kinematics::MoveMapOP init_map = new core::kinematics::MoveMap;
	init_map->set_jump( false );

	for ( claims::DofClaims::const_iterator it = init_dofs.begin(), eit = init_dofs.end(); it != eit; ++it )
	{
		if ( (*it)->owner()==this && (*it)->str_type() == "JUMP") {
			if(tr.Trace.visible())
			{
				(*it)->show(tr.Trace);
				tr.Trace << std::endl;
			}
			(*it)->toggle( *init_map, true );
		}else{
			if(tr.Trace.visible())
			{
				(*it)->show(tr.Trace);
				tr.Trace << std::endl;
				tr.Trace << "No need to init dof in MoveMap" << std::endl;
			}
		}
	}
	//Move TMHs to origin
	move_spans(pose);
}

//@brief Make extended chain an idealized helix
void
TMHTopologySamplerClaimer::set_pose_torsions(core::pose::Pose& pose)
{
	for(core::Size i=1;i<=pose.total_residue();i++)
	{
		if(pose.residue(i).is_virtual_residue() || pose.residue(i).name3() == "XXX")
		{
			continue;
		}
		else{
			//set phi and psi to make "ideal" helix
			pose.set_phi(i,-58.0);
			pose.set_psi(i,-47.0);
		}
	}
}

//@brief move helices closer together
void
TMHTopologySamplerClaimer::move_spans(core::pose::Pose& pose)
{
	//calculate and output membrane center and normal
	core::Vector membrane_vector = output_membrane_vector(pose);

	//get the membrane_topology
	core::scoring::MembraneTopologyOP membrane_topology = get_membrane_topology(pose);

	//get the grid_points vector
	utility::vector1<core::Vector> grid_points = pre_compute_grid_points(pose);

	//jump_index = span#  loop through all spans and move each one at a time
	for(core::Size jump_index = 1; jump_index <= pose.num_jump(); jump_index++)
	{
		//Since I've added residues to pose by jump, have more jumps than helices
		if(!(pose.residue(pose.fold_tree().jump_edge(jump_index).stop()).is_virtual_residue()) &&
				!(pose.residue(pose.fold_tree().jump_edge(jump_index).stop()).name3()=="XXX"))
		{
			if(tr.Debug.visible())
			{
				tr.Debug << "jump " << jump_index << "\t" << pose.fold_tree().jump_edge(jump_index).start() << " " <<
						pose.fold_tree().jump_edge(jump_index).stop() << std::endl;
			}

			//Compute a helix vector between a CA and n+4 CA
			core::Size current_span_CoM(pose.fold_tree().jump_edge(jump_index).stop());
			core::Vector rb_centroid(pose.residue(current_span_CoM).xyz("CA"));
			core::Vector rb_centroid_second(0,0,0);
			if(pose.residue(current_span_CoM+4).seqpos() >= pose.residue(membrane_topology->span_end(jump_index)).seqpos())
			{
				rb_centroid_second = pose.residue(current_span_CoM-4).xyz("CA");
			}else{
				rb_centroid_second = pose.residue(current_span_CoM+4).xyz("CA");
			}
			core::Vector current_helix_vector(rb_centroid - rb_centroid_second);

			//Compute the angle between membrane vector and helix vector
			// Why the hell doesn't xyzVector::angle_of() compile?  copy/pasting code here:
			core::Real angle(0);
			core::Real const mag = membrane_vector.length() *current_helix_vector.length();
			angle =  (mag > 0.0) ? std::acos( numeric::sin_cos_range( membrane_vector.dot(current_helix_vector  ) / mag ) ) :  0  ;
			angle = numeric::conversions::degrees(angle);
			if(tr.Trace.visible()){tr.Trace << "angle between membrane_vector and helix_vector:  " << angle << std::endl;};
			if(membrane_topology->helix_id(jump_index) % 2 == 0)
			{
				angle+=180;
				if(tr.Trace.visible()){tr.Trace << "angle between membrane_vector and helix_vector:  " << angle << std::endl;}
			}

			//compute a spin axis by taking cross product of helix vector and membrane vector.
			core::Vector current_spin_axis = current_helix_vector.cross_product(membrane_vector);

			//use spinmover to rotate helix vector by calculated angle along the spin axis
			protocols::rigid::RigidBodyDeterministicSpinMover spin_mover(jump_index,current_spin_axis,rb_centroid,angle);
			spin_mover.apply(pose);

			//Get the vector of grid points and make sure it's the same size as the number of jumps	core::Size grid_point(0);
			core::Size grid_point(0);
			core::Size random_gp = numeric::random::RG.random_range(1,grid_points.size());
			core::Vector desired_centroid(0.0,0.0,0.0);

			grid_point = random_gp;

			desired_centroid = grid_points[grid_point];
			if(tr.Trace.visible()){
				tr.Debug << "jump:  " << jump_index << "; grid_point:  " << grid_point << " " << grid_points[grid_point].x() << " " << grid_points[grid_point].y() << " " << grid_points[grid_point].z() << std::endl;
				tr.Debug << "grid_points.size() before erase:  " << grid_points.size() << std::endl;
			}
			grid_points.erase(grid_points.begin()+(grid_point-1));
			if(tr.Trace.visible()){tr.Trace << "grid_points.size() after erase:  " << grid_points.size() << std::endl;}

			//Move helices to origin (for now)
			core::Vector trans_vec = desired_centroid - rb_centroid;
			core::Real trans_len = trans_vec.length();
			if (trans_len > 1e-3) { // otherwise we get NaNs
				protocols::rigid::RigidBodyTransMover mover(pose, jump_index);
				mover.step_size(trans_len);
				mover.trans_axis(trans_vec);
				mover.apply(pose);
			}
		}
	}
	current_pose_ = pose;
	if (tr.Debug.visible()){
		pose.dump_pdb("test_move_spans.pdb","test_move_spans");
	}
}

core::Vector
TMHTopologySamplerClaimer::output_membrane_vector(core::pose::Pose& pose)
{
	core::scoring::MembraneEmbed membrane_embed = get_membrane_embed(pose);
	//get the center and normal from membrane embed.
	//compute center-normal to get the membrane vector
	//this can be done outside the loop here
	core::Vector membrane_center(membrane_embed.center());
	core::Vector membrane_normal(membrane_embed.normal());
	core::Vector membrane_vector = (membrane_center - membrane_normal);

	//virtual residues stuff
	bool fullatom = pose.is_fullatom();
	core::chemical::ResidueTypeSetCAP const &residue_set( core::chemical::ChemicalManager::get_instance()->residue_type_set
		( fullatom ? core::chemical::FA_STANDARD : core::chemical::CENTROID ));
	core::chemical::ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map("VRT") );
	core::conformation::ResidueOP membrane_center_res = core::conformation::ResidueFactory::create_residue( *(rsd_type_list[1]) ) ;
	core::conformation::ResidueOP membrane_normal_res = core::conformation::ResidueFactory::create_residue( *(rsd_type_list[1]) ) ;

	if(option[basic::options::OptionKeys::out::file::output_virtual].user())
	{
		for ( Size j=1; j<= membrane_center_res->natoms(); ++j ) {
			membrane_center_res->atom(j).xyz( membrane_center_res->atom(j).xyz() + membrane_center);
		}
		for ( Size j=1; j<= membrane_normal_res->natoms(); ++j ) {
			membrane_normal_res->atom(j).xyz( membrane_normal_res->atom(j).xyz() + membrane_normal);
		}

		//add virtual residues to pose so can output when dump pose to pdb
		core::Size root = pose.residue(pose.fold_tree().root()).seqpos();
		pose.append_residue_by_jump( *membrane_center_res,root);
		pose.append_residue_by_jump( *membrane_normal_res,root);
	}
	return membrane_vector;
}

///@brief pre-compute grid points where you want to move the helices
utility::vector1<core::Vector>
TMHTopologySamplerClaimer::pre_compute_grid_points(core::pose::Pose& pose)
{
	utility::vector1<core::Vector> grid_points;
	utility::vector1<core::Vector> used_grid_points;
	core::Vector origin(get_membrane_embed(pose).center());
	core::Vector new_vector(0,0,0);
	core::Vector add_vector(0,0,0);
	core::Size sequence_counter(1);
	while(grid_points.size()<tmhelix_)
	{
		if(sequence_counter>6)
		{
			origin = *grid_points.rbegin();
			sequence_counter = 1;
		}
		switch(sequence_counter)
		{
			tr.Debug << "sequence_counter:  " << sequence_counter << std::endl;
			case 1:
				add_vector = core::Vector(15,0,0);
				break;
			case 2:
				add_vector = core::Vector(7.5,sqrt(168.75),0); // this will overlap in 2nd round
				break;
			case 3:
				add_vector = core::Vector(-7.5,sqrt(168.75),0);
				break;
			case 4:
				add_vector = core::Vector(-15,0,0); // this will overlab in second round
				break;
			case 5:
				add_vector = core::Vector(-7.5,-(sqrt(168.75)),0);
				break;
			case 6:
				add_vector = core::Vector(7.5,-(sqrt(168.75)),0);
				break;
		}
		bool gp_exists(false);
		for(utility::vector1<core::Vector>::iterator used_gp_it=used_grid_points.begin(); used_gp_it!=used_grid_points.end(); used_gp_it++)
		{
			if(*used_gp_it == origin+add_vector)
			{
				if(tr.Trace.visible())
				{
					tr.Trace << "*used_gp_it:  " << used_gp_it->x() << " " << used_gp_it->y() << " " << used_gp_it->z() << std::endl;
					tr.Trace << "origin+add_vector:  " << (origin+add_vector).x() << " " << (origin+add_vector).y() << " " << (origin+add_vector).z() << std::endl;
					tr.Trace << "grid point already exists, try another!" << std::endl;
				}
				gp_exists = true;
			}
		}
		if(gp_exists==true)
		{
			tr.Debug << "gp_exists is TRUE" << std::endl;
			sequence_counter++;
			continue;
		}else{
			new_vector = origin + add_vector;
			grid_points.push_back(new_vector);
			used_grid_points.push_back(new_vector);
			if(tr.Debug.visible())
			{
				tr.Debug << "gp_exists is FALSE" << std::endl;
				tr.Debug << "add_vector xyz():  x=" << add_vector.x() << "; y=" <<
						add_vector.y() << "; z=" << add_vector.z() << std::endl;
				tr.Debug << "grid point " << sequence_counter << " new_vector xyz():  x=" << new_vector.x() << "; y=" <<
						new_vector.y() << "; z=" << new_vector.z() << std::endl;
				tr.Debug << "grid_points.size() after:  " << grid_points.size() << std::endl;
				tr.Debug << "used_grid_points.size() after:  " << used_grid_points.size() << std::endl;
			}
			sequence_counter++;
		}
	}
	return grid_points;
}

///@brief claimers can add movers to the RandomMover (Container).
/// add your moves, make it dependent on stage if you want to. So far this is called only by abinitio...
/// if you don't want to do anything special --- don't overload this method!
/// default: adds mover given by virtual call get_mover()  with stage-dependent weight given by abinitio_mover_weight_
void
TMHTopologySamplerClaimer::add_mover(
	moves::RandomMover& random_mover,
	core::pose::Pose const& /*pose*/,
	abinitio::StageID stageID,
	core::scoring::ScoreFunction const& /*scorefxn*/,
	core::Real /*progress*/
)
{
	core::Real weight(0.0);
	if(stageID == 1)
	{
		weight = rb_mover_stage1_weight_;
	}
//	}else if(stageID == 2)
//	{
//		weight = 1.0;
//	}
//	if(stageID == 1 || stageID == 2)
	if(stageID == 1)
	{
		if(tr.Trace.visible()) { tr.Trace << "stageID is " << stageID << " so adding RBMover!" << std::endl; }
		tr.Debug << "rotation_mag:  " << rotation_mag_ << " translation_mag:  " << translation_mag_ << " num_jump_moved:  " << tmhelix_ << std::endl;
		//Use the RigidBodyRandomTMHMover to move TMHs.  The mover chooses a random jump from the pose and does a RB perturb on it given the rotation and translation
		//If the helix moves too far away from starting position, moves it back to initial position so it doesnt escape out into space
		core::Real max_trans = (6*tmhelix_)+10;
		RigidBodyRandomTMHMoverOP RBMover = new RigidBodyRandomTMHMover(max_trans,rotation_mag_,translation_mag_,tmhelix_,this);
		//Add RBMover to the random_mover, which is then passed on to TopologyBroker::mover(pose, stageID, scorefxn, progress)
		random_mover.add_mover(RBMover,weight);
	}else{
		if(tr.Trace.visible())
		{
			tr.Trace << "stageID is " << stageID << std::endl;
		}
	}
}

} // topology_broker
} // protocols




