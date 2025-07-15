// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com), adapted from the ResfileReader code
/// by Steven Lewis (smlewi@gmail.com) and Andrew Leaver-Fay

// Unit Headers
#include <protocols/ligand_docking/MinimizeBackbone.hh>
#include <protocols/ligand_docking/MinimizeBackboneCreator.hh>
#include <protocols/ligand_docking/InterfaceBuilder.hh>
#include <protocols/ligand_docking/ligand_options/Interface.hh>
#include <protocols/loops/Loop.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

//Project Headers
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/pose/Pose.hh>


#include <core/chemical/VariantType.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

// Utility Headers
#include <numeric/random/random.hh>
#include <utility>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>


//STL headers
#include <algorithm> // for std::min, find_if

#include <core/pose/variant_util.hh>
#include <core/pose/chains_util.hh>
#include <protocols/ligand_docking/LigandArea.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/tree/Atom.hh>

// Boost Headers
#include <utility/excn/Exceptions.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace ligand_docking {

static basic::Tracer TR( "protocols.ligand_docking.ligand_options.MinimizeBackbone" );




MinimizeBackbone::MinimizeBackbone():
	protocols::moves::Mover(),
	interface_builder_()
{}

MinimizeBackbone::MinimizeBackbone(InterfaceBuilderOP interface_builder):
	protocols::moves::Mover(),
	interface_builder_(std::move(interface_builder))
{}

MinimizeBackbone::MinimizeBackbone(MinimizeBackbone const & /*that*/) = default;

MinimizeBackbone::~MinimizeBackbone() = default;

protocols::moves::MoverOP MinimizeBackbone::clone() const {
	return utility::pointer::make_shared< MinimizeBackbone >( *this );
}

protocols::moves::MoverOP MinimizeBackbone::fresh_instance() const {
	return utility::pointer::make_shared< MinimizeBackbone >();
}


//@brief parse XML (specifically in the context of the parser/scripting scheme)
void
MinimizeBackbone::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	if ( tag->getName() != "MinimizeBackbone" ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This should be impossible");

	if ( ! tag->hasOption("interface") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'HighResDocker' requires 'interface' tag (comma separated chains to dock)");

	std::string interface_name= tag->getOption<std::string>("interface");
	interface_builder_= datamap.get_ptr<protocols::ligand_docking::InterfaceBuilder>( "interface_builders", interface_name);
}

/// setting up backbone_minimization foldtree
void
MinimizeBackbone::apply( core::pose::Pose & pose ){
	debug_assert(pose.fold_tree().check_edges_for_atom_info());
	debug_assert(interface_builder_);// make sure the pointer points
	constraints_.clear(); // Reset constraints for new pose
	cutpoints_.clear(); // Reset cutpoints for new pose.
	ligand_options::Interface interface= interface_builder_->build(pose);
	restrict_to_protein_residues(interface, pose);
	TR.Debug << "interface for fold_tree: "<< interface << std::endl;
	debug_assert(pose.fold_tree().check_edges_for_atom_info());
	reorder_foldtree_around_mobile_regions(interface, pose);
	restrain_protein_Calphas(interface, pose);
	debug_assert(pose.fold_tree().check_edges_for_atom_info());
}

void
MinimizeBackbone::remove_constraints( core::pose::Pose & pose ) {
	for ( core::Size ii(1); ii <= constraints_.size(); ++ii ) {
		pose.remove_constraint(constraints_[ii],true);
	}
}

void
MinimizeBackbone::remove_cutpoints( core::pose::Pose & pose ) {
	for ( core::Size ii(1); ii <= cutpoints_.size(); ++ii ) {
		core::Size cutpt( cutpoints_[ii] );
		if ( cutpt+1 > pose.total_residue() ) {
			TR.Warning << "Warning: pose no longer has one of the residues listed as a surrounding a cutpoint." << std::endl;
			continue;
		}
		core::pose::remove_variant_type_from_pose_residue(pose, core::chemical::CUTPOINT_LOWER, cutpt);
		core::pose::remove_variant_type_from_pose_residue(pose, core::chemical::CUTPOINT_UPPER, cutpt + 1);
	}
}

void MinimizeBackbone::reorder_foldtree_around_mobile_regions(
	ligand_options::Interface const & interface,
	core::pose::Pose & pose
) {
	debug_assert(pose.fold_tree().check_edges_for_atom_info());
	core::kinematics::FoldTree const & fold_tree_copy = pose.fold_tree();
	debug_assert(fold_tree_copy.check_edges_for_atom_info());
	TR.Debug << "Initial foldtree " << fold_tree_copy << std::endl;
	core::kinematics::FoldTreeOP better_ligand_jumps = create_fold_tree_with_ligand_jumps_from_attach_pts(fold_tree_copy, interface, pose);
	debug_assert(better_ligand_jumps->check_edges_for_atom_info());
	TR.Debug << "foldtree with ligand jumps from attach pts" << *better_ligand_jumps << std::endl;
	core::kinematics::FoldTreeOP with_cutpoints = create_fold_tree_with_cutpoints(better_ligand_jumps, interface, pose);
	TR.Debug  << "foldtree with cutpoints" << *with_cutpoints << std::endl;
	debug_assert(with_cutpoints->check_edges_for_atom_info());
	with_cutpoints->delete_extra_vertices();
	debug_assert(with_cutpoints->check_edges_for_atom_info());
	TR.Debug << "foldtree with less vertices" << *with_cutpoints << std::endl;
	debug_assert(with_cutpoints->check_edges_for_atom_info());
	reorder_with_first_non_mobile_as_root(with_cutpoints, interface, pose);
	TR.Debug << "Final loops foldtree " << *with_cutpoints << std::endl;
	if ( !with_cutpoints->check_fold_tree() ) {
		utility_exit_with_message("Invalid fold tree after trying to set up for minimization!");
	}
	debug_assert(with_cutpoints->check_edges_for_atom_info());
	pose.fold_tree(*with_cutpoints);
	debug_assert(pose.fold_tree().check_edges_for_atom_info());
}

void reorder_with_first_non_mobile_as_root(
	core::kinematics::FoldTreeOP f,
	ligand_options::Interface const & interface,
	core::pose::Pose & pose
) {
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue(i).is_polymer() && interface[i].type == ligand_options::InterfaceInfo::non_interface ) {
			f->reorder(i);
			return;
		}
	}
	utility_exit_with_message("Every residue is mobile!  We need a non-mobile residue");
}

core::kinematics::FoldTreeOP MinimizeBackbone::create_fold_tree_with_cutpoints(
	core::kinematics::FoldTreeCOP f,
	ligand_options::Interface const & interface,
	core::pose::Pose & pose
) {
	core::kinematics::FoldTreeOP f_new( new core::kinematics::FoldTree() );// Deleting edges is bad so add them to a new foldtree

	if ( pose.conformation().is_membrane() ) {
		f_new = utility::pointer::make_shared< core::kinematics::FoldTree >( pose.fold_tree() );
		return f_new;
	}

	int new_jump = f->num_jump();

	for ( auto const & edge_itr : *f ) {

		if ( !edge_itr.is_polymer() ) {
			f_new->add_edge(edge_itr);
			continue; // only subdivide "peptide" edges of the fold tree
		}

		utility::vector1< protocols::loops::Loop> loops = add_cut_points( edge_itr, interface, pose);

		const int e_start = edge_itr.start();
		const int e_stop = edge_itr.stop();

		int last_rigid = e_start;
		for ( core::Size i = 1; i <= loops.size(); ++i ) {
			protocols::loops::Loop const & l = loops[i];
			int first = std::max(int(e_start), int(l.start() - 1));
			int last = std::min(int(e_stop), int(l.stop() + 1));
			if ( first != last_rigid ) {
				f_new->add_edge(last_rigid, first, core::kinematics::Edge::PEPTIDE);
			}
			f_new->add_edge(first, l.cut(), core::kinematics::Edge::PEPTIDE);
			f_new->add_edge(core::kinematics::Edge(first, last, ++new_jump, "CA", "CA", false ));
			f_new->add_edge(last, l.cut() + 1, core::kinematics::Edge::PEPTIDE);

			last_rigid = last;
		}

		if ( last_rigid != int(e_stop) ) {
			f_new->add_edge(last_rigid, e_stop, core::kinematics::Edge::PEPTIDE);
		}
	}
	return f_new;
}

utility::vector1< protocols::loops::Loop> MinimizeBackbone::add_cut_points(
	core::kinematics::Edge const & edge,
	ligand_options::Interface const & interface,
	core::pose::Pose & pose
) {
	if ( edge.start() > edge.stop() ) {
		utility_exit_with_message("Not prepared to deal with backwards fold tree edges!");
	}
	utility::vector1< protocols::loops::Loop> loops;

	core::Size start= interface.find_first_interface_residue(edge.start(), edge.stop());
	core::Size stop= 0;

	while ( start ) {
		stop= interface.find_stop_of_this_interface_region(start, edge.stop());

		TR.Debug << "start,stop:"<< start << ','<< stop << std::endl;
		TR.Debug << "edge.start,edge.stop:"<< edge.start() << ','<< edge.stop() << std::endl;
		runtime_assert( stop >= start );
		if ( (stop - start + 1) < 4 ) {
			TR.Warning
				<< "for backbone minimization to work properly, a stretch of at least 4 residues needs to be allowed to move. Stretch between "
				<< start << " and " << stop << " is too short.  This should never happen after extending the interface"
				<< std::endl;
			continue;
		}

		// Cutpoint should fall between start and stop, but not if the rsd on either side is a terminus.
		// Because this happens within one peptide edge, no "internal" residue should ever be a terminus.
		core::Size const cut_start = (pose.residue(start).is_terminus() ? start + 1 : start);
		core::Size const cut_end = ( pose.residue(stop).is_terminus() ? stop - 2 : stop - 1);

		runtime_assert( cut_start <= cut_end );

		auto cutpt = core::Size(numeric::random::rg().random_range(cut_start, cut_end)); // cut is made between cutpt and cutpt+1

		loops.push_back(protocols::loops::Loop(start, stop, cutpt));
		// also need to set up residue variants so chainbreak score works correctly!
		cutpoints_.push_back( cutpt ); // So we can remove the variants later
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::CUTPOINT_LOWER, cutpt);
		core::pose::add_variant_type_to_pose_residue(pose, core::chemical::CUTPOINT_UPPER, cutpt + 1);

		if ( stop == edge.stop() ) break; // cast the stop as an int
		start= interface.find_start_of_next_interface_region(stop+1,  edge.stop());
	}
	return loops;
}

core::kinematics::FoldTreeOP
MinimizeBackbone::create_fold_tree_with_ligand_jumps_from_attach_pts(
	core::kinematics::FoldTree const & f_const,
	ligand_options::Interface const & interface,
	core::pose::Pose & pose

) const {
	std::map<core::Size, core::Size> jump_to_attach = find_attach_pts(interface, pose);

	core::kinematics::FoldTreeOP new_fold_tree( new core::kinematics::FoldTree() );

	if ( pose.conformation().is_membrane() ) {
		new_fold_tree = utility::pointer::make_shared< core::kinematics::FoldTree >( pose.fold_tree() );
	} else {

		for ( auto const & e : f_const ) {
			std::map<core::Size, core::Size>::const_iterator const jump_attach = jump_to_attach.find(e.label());
			if ( jump_attach != jump_to_attach.end() ) {
				core::Size const jump_id = jump_attach->first;
				core::Size const attach_pt = jump_attach->second;
				core::Size const ligand_residue_id = pose.fold_tree().downstream_jump_residue(jump_id);
				new_fold_tree->add_edge(attach_pt, ligand_residue_id, jump_id);
			} else {
				core::Size const attach_pt= find_peptide_attach_pt(e.start(), e.stop(), jump_to_attach);

				if ( e.label() == core::kinematics::Edge::PEPTIDE && attach_pt != 0 ) {
					new_fold_tree->add_edge(e.start(), attach_pt, core::kinematics::Edge::PEPTIDE);
					new_fold_tree->add_edge(attach_pt, e.stop(), core::kinematics::Edge::PEPTIDE);
				} else if ( e.label() == core::kinematics::Edge::CHEMICAL ) {
					new_fold_tree->add_edge(e.start(), e.stop(), e.start_atom(), e.stop_atom());
				} else { // add a jump edge
					new_fold_tree->add_edge(e.start(), e.stop(), e.label());
				}
			}
		}
	}

	if ( !new_fold_tree->check_fold_tree() ) {
		utility_exit_with_message("Fold tree did not pass check!");
	}
	if ( !new_fold_tree->check_edges_for_atom_info() ) {
		utility_exit_with_message("Fold tree has chemical edge problems!");
	}
	if ( f_const.nres() != new_fold_tree->nres() ) {
		utility_exit_with_message("Number of residues changed?!");
	}
	if ( f_const.num_jump() != new_fold_tree->num_jump() ) {
		utility_exit_with_message("Number of jumps changed?!");
	}
	if ( !new_fold_tree->connected() ) {
		utility_exit_with_message("Fold tree not connected?!");
	}
	return new_fold_tree;
}

/// @detail: residues that are near the interface residues will be allowed to move with restraint
void MinimizeBackbone::restrain_protein_Calphas(
	ligand_options::Interface const & interface,
	core::pose::Pose & pose
) {
	core::id::AtomID fixed_pt(pose.atom_tree().root()->atom_id());
	TR.Debug << interface << std::endl;
	for ( core::Size residue_id = 1; residue_id <= interface.size(); ++residue_id ) { // didn't use an iterator because pose needs a #
		if ( interface[residue_id].type != ligand_options::InterfaceInfo::non_interface ) {
			restrain_protein_Calpha(interface, pose, residue_id, fixed_pt);
		}
	}
}

void MinimizeBackbone::restrain_protein_Calpha(
	ligand_options::Interface const & interface,
	core::pose::Pose & pose,
	core::Size residue_id,
	core::id::AtomID const & fixed_pt
) {
	core::conformation::Residue const & residue = pose.residue(residue_id);
	if ( !residue.is_protein() ) {
		//Skip non-protein residues even if they are part of an interface
		return;
	}
	std::string const & ligand_chain= interface[residue_id].chain;
	std::map<std::string, LigandAreaOP> const & ligand_areas= interface_builder_->get_ligand_areas();
	auto found= ligand_areas.find(ligand_chain);
	debug_assert( found != ligand_areas.end() );// this shouldn't be possible
	LigandAreaOP const ligand_area= found->second;
	core::id::AtomID const atom_ID(residue.atom_index("CA"), residue_id);
	core::scoring::func::FuncOP const harmonic_function( new core::scoring::func::HarmonicFunc(0, ligand_area->Calpha_restraints_) );
	core::scoring::constraints::ConstraintOP const constraint( new core::scoring::constraints::CoordinateConstraint(
		atom_ID,
		fixed_pt,
		residue.xyz("CA"),
		harmonic_function
		) )
		;
	TR.Debug << "Restraining C-alpha of residue " << residue_id << " with "<< ligand_area->Calpha_restraints_<< std::endl;

	pose.add_constraint(constraint);
	constraints_.push_back(constraint);
}

std::string MinimizeBackbone::get_name() const {
	return mover_name();
}

std::string MinimizeBackbone::mover_name() {
	return "MinimizeBackbone";
}

void MinimizeBackbone::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("interface", xs_string, "Comma separated chains to dock.");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Perform backbone mimimization of a docked structure.", attlist );
}

std::string MinimizeBackboneCreator::keyname() const {
	return MinimizeBackbone::mover_name();
}

protocols::moves::MoverOP
MinimizeBackboneCreator::create_mover() const {
	return utility::pointer::make_shared< MinimizeBackbone >();
}

void MinimizeBackboneCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MinimizeBackbone::provide_xml_schema( xsd );
}


std::map<core::Size, core::Size> MinimizeBackbone::find_attach_pts(
	ligand_options::Interface const & interface,
	core::pose::Pose const & pose
) const {
	std::map<core::Size, core::Size> jumpToAttach;
	std::map<std::string, LigandAreaOP> const ligand_areas= interface_builder_->get_ligand_areas();
	//std::map<char, LigandAreaOP>::const_iterator index = ligand_areas.begin();
	for ( LigandAreas::value_type const & ligand_area_pair : ligand_areas ) {
		//for (; index != ligand_areas.end(); ++index) { // these are just the ligand chains to dock
		utility::vector1<core::Size> jump_ids = core::pose::get_jump_ids_from_chain(ligand_area_pair.first, pose);
		for ( core::Size const jump_id : jump_ids ) {
			jumpToAttach[jump_id] = find_attach_pt(jump_id, interface, pose);
		}
	}
	return jumpToAttach;
}

//std::map<char, LigandArea>
//MinimizeBackbone::get_ligand_areas(){
// return interface_builder_->get_ligand_areas();
//}

core::Size find_peptide_attach_pt (
	int const & start,
	int const & stop,
	std::map<core::Size, core::Size> const & jump_to_attach
){
	auto index =
		jump_to_attach.begin();
	std::map<core::Size, core::Size>::const_iterator const end =
		jump_to_attach.end();
	for ( ; index != end; ++index ) {
		auto const attach_pt = (int) index->second;

		if ( start < attach_pt && attach_pt < stop ) {
			return index->second;
		}
		if ( stop < attach_pt && attach_pt < start ) {
			return index->second;
		}
	}
	return 0;
}


core::Size find_attach_pt(
	core::Size const jump_id,
	ligand_options::Interface const & interface,
	core::pose::Pose const & pose
){
	core::Size attach_pt = 1; // for now, attach ligand to residue 1
	core::Real shortest_dist2 = 1e99;
	//core::Vector center= protocols::geometry::downstream_centroid_by_jump(pose, jump_id);
	const core::Size residue_id = pose.fold_tree().downstream_jump_residue(jump_id);
	const core::conformation::Residue & residue = pose.residue(residue_id);
	const core::Vector lig_nbr_vector = residue.nbr_atom_xyz();
	//for(core::Size i = 1; i <= pose.size(); ++i) { // only look at upstream residues, else multi-residue ligand will attach to itself
	for ( core::Size i = 1; i < residue_id; ++i ) {
		if ( !pose.residue(i).is_protein() && interface[i].type != ligand_options::InterfaceInfo::non_interface ) {
			core::Vector res_coords;
			if ( pose.residue(i).has("CA") ) {
				res_coords = pose.residue(i).xyz("CA");
			} else {
				res_coords = pose.residue(i).nbr_atom_xyz();
			}
			core::Real const new_dist2 = lig_nbr_vector.distance_squared(
				res_coords);
			if ( new_dist2 < shortest_dist2 ) {
				shortest_dist2 = new_dist2;
				attach_pt = i;
			}
		}
	}

	return attach_pt;
}

void restrict_to_protein_residues(
	ligand_options::Interface & interface,
	core::pose::Pose const & pose
){
	for ( core::Size i=1; i <= interface.size(); ++i ) {
		if ( interface[i].type == ligand_options::InterfaceInfo::is_interface
				&& pose.residue(i).is_ligand()
				) {
			interface[i].type = ligand_options::InterfaceInfo::non_interface;
		}
	}
}

} //namespace ligand_docking
} //namespace protocols
