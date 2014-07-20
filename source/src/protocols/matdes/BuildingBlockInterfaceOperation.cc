// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/matdes/BuildingBlockInterfaceOperation2.0.hh
/// @brief  Restrict design to only residues at inter-building block interfaces
/// @author Will Sheffler (willsheffler@gmail.com) Jacob Bale (balej@uw.edu)

// Unit Headers
#include <protocols/matdes/BuildingBlockInterfaceOperation.hh>
#include <protocols/matdes/BuildingBlockInterfaceOperationCreator.hh>

// Project Headers
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>


static basic::Tracer TR("protocols.matdes.BuildingBlockInterfaceOperation" );

namespace protocols {
namespace matdes {

core::pack::task::operation::TaskOperationOP
BuildingBlockInterfaceOperationCreator::create_task_operation() const
{
	return new BuildingBlockInterfaceOperation;
}


BuildingBlockInterfaceOperation::BuildingBlockInterfaceOperation( core::Size nsub_bblock, std::string sym_dof_names, core::Real contact_dist /* = 10*/, core::Real bblock_dist /*= 5 */, core::Real fa_rep_cut /* = 3.0 */, bool filter_intrabb, bool intrabb_only, bool multicomponent ):
	nsub_bblock_(nsub_bblock),
	sym_dof_names_(sym_dof_names),
	contact_dist_(contact_dist),
	bblock_dist_(bblock_dist),
	fa_rep_cut_(fa_rep_cut),
	filter_intrabb_(filter_intrabb),
	intrabb_only_(intrabb_only),
	multicomponent_(multicomponent)
{}

BuildingBlockInterfaceOperation::~BuildingBlockInterfaceOperation() {}

core::pack::task::operation::TaskOperationOP BuildingBlockInterfaceOperation::clone() const
{
	return new BuildingBlockInterfaceOperation( *this );
}

void
BuildingBlockInterfaceOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace core;
	using namespace basic;
	using namespace pose;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;
	using namespace scoring;
	using namespace utility;
	typedef vector1<Size> Sizes;

	utility::vector1<std::string> sym_dof_name_list;
	if( sym_dof_names_ == "" ) {
		sym_dof_name_list = sym_dof_names( pose );
	} else {
		sym_dof_name_list = utility::string_split( sym_dof_names_ , ',' );
	}

	Sizes intra_subs1, intra_subs2;
	if( multicomponent_ ) {
		runtime_assert (sym_dof_name_list.size() == 2);  //fpd  multicomponent code assumes this holds
		intra_subs1 = get_jump_name_to_subunits(pose,sym_dof_name_list[1]);
		intra_subs2 = get_jump_name_to_subunits(pose,sym_dof_name_list[2]);
	}

	// Find out which positions are near the inter-subunit interfaces
	// These will be further screened below, then passed to design()
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	vector1<bool> indy_resis = sym_info->independent_residues();
	Real const contact_dist_sq = contact_dist_ * contact_dist_;
	Sizes design_pos; 
	std::set<Size> filtered_design_pos;
	Sizes comp_chains;
	std::string select_interface_pos("select interface_pos, resi ");
	std::string select_comp1_chains("select comp1_chains, chain ");
	std::string select_comp2_chains("select comp2_chains, chain ");
	for(Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
		if(sym_info->subunit_index(ir) != 1) continue;
		std::string atom_i = (pose.residue(ir).name3() == "GLY") ? "CA" : "CB";
		for(Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
			std::string atom_j = (pose.residue(jr).name3() == "GLY") ? "CA" : "CB";

			//If one component, then check for clashes between all residues in primary subunit and subunits with indices > nsub_bb
			if( !multicomponent_ && sym_info->subunit_index(jr) <= nsub_bblock_ ) continue;

			//If two component, then check for clashes between all residues in primary subunitA and other building blocks, and all resis in primary subB and other building blocks. 
			if( multicomponent_ ) {
				Sizes const & isubs( get_component_of_residue(pose,ir)=='A'?intra_subs1:intra_subs2);
				if (find(comp_chains.begin(),comp_chains.end(),pose.chain(jr))==comp_chains.end()) {
					if (get_component_of_residue(pose,jr)=='A') {
						select_comp1_chains.append("\"" + std::string(1, pose.pdb_info()->chain( jr )) +
								"\"+"); // TODO: Jacob
						comp_chains.push_back(pose.chain(jr));
					} else if (get_component_of_residue(pose,jr)!='A') {
						select_comp2_chains.append("\"" + std::string(1, pose.pdb_info()->chain( jr )) +
								"\"+"); // TODO: Jacob
						comp_chains.push_back(pose.chain(jr));
					}
				}	
				if(get_component_of_residue(pose,ir)==get_component_of_residue(pose,jr)&&find(isubs.begin(),isubs.end(),sym_info->subunit_index(jr))!=isubs.end()) continue;
			}

			if(pose.residue(ir).xyz(atom_i).distance_squared(pose.residue(jr).xyz(atom_j)) <= contact_dist_sq) {
				design_pos.push_back(ir);
				TR.Debug << ir << std::endl;
				core::Size output_resi = ir;
				if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
					output_resi = pose.pdb_info()->number( ir );
				}
				select_interface_pos.append(ObjexxFCL::string_of(output_resi) + "+");   
				break;
			}
		}
	}
	TR << select_interface_pos << std::endl;
	TR.Debug << select_comp1_chains << std::endl;
	TR.Debug << select_comp2_chains << std::endl;
	
	// Here we filter the residues that we are selecting for design
	// to get rid of those that make intra-building block interactions
	Pose scored_pose( pose );
	get_score_function()->score(scored_pose);
	Real bblock_dist_sq = bblock_dist_ * bblock_dist_;
	std::string select_filtered_interface_pos("select filtered_interface_pos, resi ");
	bool contact;

	for(Size iip=1; iip<=design_pos.size(); iip++) {
		Size ir = design_pos[iip];
		if ( filter_intrabb_ ) {
			TR.Debug << "Filtering: Checking resi: " << ir << std::endl;
			contact = true;
			for(Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
				if( !multicomponent_ ) {
					if(sym_info->subunit_index(ir) > nsub_bblock_ || sym_info->subunit_index(jr) > nsub_bblock_) continue;
				}	else {
					Sizes const & intra_subs(get_component_of_residue(pose,ir)=='A'?intra_subs1:intra_subs2);
					if(get_component_of_residue(pose,ir)!=get_component_of_residue(pose,jr)) continue;
					if(find(intra_subs.begin(), intra_subs.end(), sym_info->subunit_index(jr)) == intra_subs.end()) continue;
				}

				if(sym_info->subunit_index(jr) == 1) continue;

				for(Size ia = 1; ia<=pose.residue(ir).nheavyatoms(); ia++) {
					for(Size ja = 1; ja<=pose.residue(jr).nheavyatoms(); ja++) {
						if(pose.residue(ir).xyz(ia).distance_squared(pose.residue(jr).xyz(ja)) <= bblock_dist_sq)	{
							if ( intrabb_only_ ) { contact = false; break; } 
							else {
								// However, if the residue in question is clashing badly (usually with a
								// residue from another building block), it needs to be designed.
								core::scoring::EnergyMap em1 = scored_pose.energies().residue_total_energies(ir);
								Real resi_fa_rep = em1[core::scoring::fa_rep];
								TR.Debug << "resi_fa_rep: " << resi_fa_rep << " fa_rep_cut_: " << fa_rep_cut_ << std::endl;
								if (resi_fa_rep < fa_rep_cut_) { contact = false; TR.Debug << "Filtered out resi: " << ir << std::endl; break; }
							}
						}
					}
					if( contact == false ) break;
				}
				if( contact == false ) break;
			}
			if( (contact && !intrabb_only_) || ((contact == false) && (intrabb_only_ == true ))) {
				filtered_design_pos.insert(ir);
				core::Size output_resi = ir;
				if ( !basic::options::option[ basic::options::OptionKeys::out::file::renumber_pdb ]() ) {
					output_resi = pose.pdb_info()->number( ir );
				}
				select_filtered_interface_pos.append(ObjexxFCL::string_of(output_resi) + "+");
			}
		} else {
			filtered_design_pos.insert(ir);
		}
	}
	TR << select_filtered_interface_pos << std::endl;
	// Now prevent_repacking at all positions that are not defined filtered design positions:
	std::string output = "design_pos ";
	for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
		if (filtered_design_pos.find(ir) != filtered_design_pos.end()) {
			output += ObjexxFCL::string_of(ir)+"+";
		} else {
			TR.Debug << "resi " << ir << " will not be designed" << std::endl;
			task.nonconst_residue_task(ir).prevent_repacking();
		}
	}
	TR.Debug << output << std::endl;
	//core::pack::make_symmetric_PackerTask_by_truncation(pose, task); // Does this need to be fixed or omitted?
}

void
BuildingBlockInterfaceOperation::parse_tag( TagCOP tag , DataMap & )
{
  nsub_bblock_ = tag->getOption<core::Size>("nsub_bblock", 1);
	sym_dof_names_ = tag->getOption< std::string >( "sym_dof_names", "" );
	contact_dist_ = tag->getOption<core::Real>("contact_dist", 10.0);
	bblock_dist_ = tag->getOption<core::Real>("bblock_dist", 5.0);
	fa_rep_cut_ = tag->getOption<core::Real>("fa_rep_cut", 3.0);
	filter_intrabb_ = tag->getOption< bool >("filter_intrabb", 1);
	intrabb_only_ = tag->getOption< bool >("intrabb_only", 0);
	multicomponent_ = tag->getOption< bool >("multicomp", 0);
}

void
BuildingBlockInterfaceOperation::parse_def( utility::lua::LuaObject const & def)
{
	nsub_bblock_ = def["nsub_bblock"] ? def["nsub_bblock"].to<core::Size>() : 1;
	contact_dist_ = def["contact_dist"] ? def["contact_dist"].to<core::Real>() : 10.0;
	bblock_dist_ = def["bblock_dist"] ? def["bblock_dist"].to<core::Real>() : 5.0;
	fa_rep_cut_ = def["fa_rep_cut"] ? def["fa_rep_cut"].to<core::Real>() : 3.0;
}

} //namespace matdes
} //namespace protocols
