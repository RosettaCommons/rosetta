// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/matdes/BuildingBlockInterfaceOperation2.0.hh
/// @brief  Restrict design to only residues at inter-building block interfaces
/// @author Will Sheffler (willsheffler@gmail.com) Jacob Bale (balej@uw.edu)

// Unit Headers
#include <protocols/matdes/BuildingBlockInterfaceOperation.hh>
#include <protocols/matdes/BuildingBlockInterfaceOperationCreator.hh>

// Project Headers
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <numeric/HomogeneousTransform.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/UniformRotationSampler.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <utility>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <utility/vector1.hh>


static basic::Tracer TR( "protocols.matdes.BuildingBlockInterfaceOperation" );

namespace protocols {
namespace matdes {

using namespace core::pack::task::operation;
using namespace utility::tag;

//fd  GIVEN: a symmetric pose and a symmetric SUBCOMPLEX
//fd  RETURN: the index of all subcomplex subunits in the pose
//fd  NOTE that both symmetries might be unbounded (e.g., a 1D fiber in a 3D lattice)
//fd    so we need to consider that the subcomplex is not "extended" as far as the pose
utility::vector1< core::Size >
get_matching_subunits( core::pose::Pose const &pose, core::pose::Pose const &subpose ) {
	runtime_assert( subpose.pdb_info() != nullptr );
	runtime_assert( core::pose::symmetry::is_symmetric( pose ) );

	// split reference pose into chains
	std::map< char, utility::vector1< numeric::xyzVector<core::Real> > > cas_by_chain;
	core::Size nres = subpose.total_residue();
	for ( core::Size i=1; i<=nres; ++i ) {
		if ( subpose.residue(i).is_protein() ) {
			char chn = subpose.pdb_info()->chain(i);
			numeric::xyzVector< core::Real > x_i = subpose.residue(i).atom(" CA ").xyz();
			cas_by_chain[chn].push_back( x_i );
		}
	}

	// "main" chain (in _ref_ pose) is 'A' (if it exists) or the first chain
	char refchain = 'A';
	if ( cas_by_chain.count('A')==0 ) {
		refchain = subpose.pdb_info()->chain(1);
	}

	core::Size nres_asu = cas_by_chain[refchain].size();
	numeric::xyzVector< core::Real > ref_com(0.,0.,0.);
	ObjexxFCL::FArray2D< core::Real > ref_coords( 3, nres_asu );
	for ( core::Size i=1; i<=nres_asu; ++i ) {
		ref_com += cas_by_chain[refchain][i];
		for ( int j=0; j<3; ++j ) ref_coords(j+1,i) = cas_by_chain[refchain][i][j];
	}
	ref_com /= nres_asu;
	for ( core::Size i=1; i<=nres_asu; ++i ) {
		for ( int j=0; j<3; ++j ) ref_coords(j+1,i) -= ref_com[j];
	}


	// get all transforms
	utility::vector1< numeric::HomogeneousTransform< core::Real > > transforms;
	for ( auto iter = cas_by_chain.begin(); iter != cas_by_chain.end(); ++iter ) {
		//char tgtchain = iter->first;
		//if (tgtchain == refchain) continue;

		// get COM
		runtime_assert( nres_asu == iter->second.size() );
		numeric::xyzVector< core::Real > tgt_com(0.,0.,0.);
		ObjexxFCL::FArray2D< core::Real > tgt_coords( 3, nres_asu );
		for ( core::Size i=1; i<=nres_asu; ++i ) {
			tgt_com += iter->second[i];
			for ( int j=0; j<3; ++j ) tgt_coords(j+1,i) = iter->second[i][j];
		}
		tgt_com /= nres_asu;
		for ( core::Size i=1; i<=nres_asu; ++i ) {
			for ( int j=0; j<3; ++j ) tgt_coords(j+1,i) -= tgt_com[j];
		}

		// get optimal superposition
		// rotate >init< to >final<
		ObjexxFCL::FArray1D< numeric::Real > ww( nres_asu, 1.0 );
		ObjexxFCL::FArray2D< numeric::Real > uu( 3, 3, 0.0 );
		numeric::Real ctx;
		float rms;

		numeric::model_quality::findUU( ref_coords, tgt_coords, ww, nres_asu, uu, ctx );
		numeric::model_quality::calc_rms_fast( rms, ref_coords, tgt_coords, ww, nres_asu, ctx );
		runtime_assert(rms < 1.0);

		numeric::xyzMatrix<core::Real> R;
		R.xx( uu(1,1) ); R.xy( uu(2,1) ); R.xz( uu(3,1) );
		R.yx( uu(1,2) ); R.yy( uu(2,2) ); R.yz( uu(3,2) );
		R.zx( uu(1,3) ); R.zy( uu(2,3) ); R.zz( uu(3,3) );
		numeric::xyzVector<core::Real> T = R*(-ref_com)+tgt_com;
		transforms.push_back( numeric::HomogeneousTransform< core::Real >(R, T) );
	}

	// compare transforms from ref to target
	utility::vector1< core::Size > ref_subunits;

	core::pose::Pose posecopy = pose;
	core::conformation::symmetry::SymmetricConformation & symmconf =
		dynamic_cast< core::conformation::symmetry::SymmetricConformation & >( posecopy.conformation() );
	core::Size nsubs = symmconf.Symmetry_Info()->subunits();
	for ( core::Size i=1; i<= nsubs; ++i ) {
		numeric::HomogeneousTransform< core::Real > ht_i =
			symmconf.get_transformation(i,true) *
			symmconf.get_transformation(1,true).inverse() ;
		for ( core::Size j=1; j<=transforms.size(); ++j ) {
			numeric::HomogeneousTransform< core::Real > &ht_j = transforms[j];
			core::Real angle_ij = numeric::urs_R2ang(
				ht_i.rotation_matrix().transpose() * ht_j.rotation_matrix() );
			core::Real distance_ij = ( ht_i.point() - ht_j.point() ).length();
			if ( angle_ij <= 20 && distance_ij <= 4 ) { // ? not sure on cutoffs
				ref_subunits.push_back(i);
				break;
			}
		}
	}

	TR.Debug << "Found " << ref_subunits.size() << " reference chains" << std::endl;
	return ref_subunits;
}

TaskOperationOP
BuildingBlockInterfaceOperationCreator::create_task_operation() const
{
	return utility::pointer::make_shared< BuildingBlockInterfaceOperation >();
}

void BuildingBlockInterfaceOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BuildingBlockInterfaceOperation::provide_xml_schema( xsd );
}

std::string BuildingBlockInterfaceOperationCreator::keyname() const
{
	return BuildingBlockInterfaceOperation::keyname();
}

BuildingBlockInterfaceOperation::BuildingBlockInterfaceOperation( core::Size nsub_bblock, std::string const & sym_dof_names, core::Real contact_dist /* = 10*/, core::Real bblock_dist /*= 5 */, core::Real fa_rep_cut /* = 3.0 */, bool filter_intrabb, bool intrabb_only, bool multicomponent ):
	nsub_bblock_(nsub_bblock),
	sym_dof_names_(sym_dof_names),
	contact_dist_(contact_dist),
	bblock_dist_(bblock_dist),
	fa_rep_cut_(fa_rep_cut),
	filter_intrabb_(filter_intrabb),
	intrabb_only_(intrabb_only),
	multicomponent_(multicomponent)
{}

BuildingBlockInterfaceOperation::~BuildingBlockInterfaceOperation() = default;

core::pack::task::operation::TaskOperationOP BuildingBlockInterfaceOperation::clone() const
{
	return utility::pointer::make_shared< BuildingBlockInterfaceOperation >( *this );
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
	using Sizes = vector1<core::Size>;

	utility::vector1<std::string> sym_dof_name_list;
	if ( sym_dof_names_ == "" ) {
		sym_dof_name_list = sym_dof_names( pose );
	} else {
		sym_dof_name_list = utility::string_split( sym_dof_names_ , ',' );
	}

	Sizes intra_subs1, intra_subs2;
	if ( multicomponent_ ) {
		runtime_assert (sym_dof_name_list.size() == 2);  //fpd  multicomponent code assumes this holds
		intra_subs1 = get_jump_name_to_subunits(pose,sym_dof_name_list[1]);
		intra_subs2 = get_jump_name_to_subunits(pose,sym_dof_name_list[2]);
	}

	utility::vector1<core::Size> include_subunits;

	if ( bblock_reference_pose_ ) {
		runtime_assert (!multicomponent_); // ref pdb & multicomponent incompatible
		if ( nsub_bblock_ != 1 ) {
			TR.Error << "WARNING! bblock_reference_pdb provided; nsub_bblock is ignored!";
		}
		include_subunits = get_matching_subunits(pose, *bblock_reference_pose_);

	} else {
		for ( core::Size i=1; i<=nsub_bblock_; ++i ) {
			include_subunits.push_back(i);
		}
	}

	// Find out which positions are near the inter-subunit interfaces
	// These will be further screened below, then passed to design()
	SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	vector1<bool> indy_resis = sym_info->independent_residues();
	Real const contact_dist_sq = contact_dist_ * contact_dist_;
	Sizes design_pos;
	std::set<core::Size> filtered_design_pos;
	Sizes comp_chains;
	std::string select_interface_pos("select interface_pos, resi ");
	std::string select_comp1_chains("select comp1_chains, chain ");
	std::string select_comp2_chains("select comp2_chains, chain ");
	for ( core::Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++ ) {
		if ( sym_info->subunit_index(ir) != 1 ) continue;
		std::string atom_i = (pose.residue(ir).name3() == "GLY") ? "CA" : "CB";
		for ( core::Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++ ) {
			std::string atom_j = (pose.residue(jr).name3() == "GLY") ? "CA" : "CB";

			//If one component, then check for clashes between all residues in primary subunit and subunits with indices > nsub_bb
			//if ( !multicomponent_ && sym_info->subunit_index(jr) <= nsub_bblock_ ) continue;
			if ( !multicomponent_ &&
					std::find( include_subunits.begin(), include_subunits.end(), sym_info->subunit_index(jr)
					) != include_subunits.end() ) continue;

			//If two component, then check for clashes between all residues in primary subunitA and other building blocks, and all resis in primary subB and other building blocks.
			if ( multicomponent_ ) {
				Sizes const & isubs( get_component_of_residue(pose,ir)=='A'?intra_subs1:intra_subs2);
				if ( find(comp_chains.begin(),comp_chains.end(),pose.chain(jr))==comp_chains.end() ) {
					if ( get_component_of_residue(pose,jr)=='A' ) {
						select_comp1_chains.append("\"" + std::string(1, pose.pdb_info()->chain( jr )) +
							"\"+"); // TODO: Jacob
						comp_chains.push_back(pose.chain(jr));
					} else if ( get_component_of_residue(pose,jr)!='A' ) {
						select_comp2_chains.append("\"" + std::string(1, pose.pdb_info()->chain( jr )) +
							"\"+"); // TODO: Jacob
						comp_chains.push_back(pose.chain(jr));
					}
				}
				if ( get_component_of_residue(pose,ir)==get_component_of_residue(pose,jr)&&find(isubs.begin(),isubs.end(),sym_info->subunit_index(jr))!=isubs.end() ) continue;
			}

			if ( pose.residue(ir).xyz(atom_i).distance_squared(pose.residue(jr).xyz(atom_j)) <= contact_dist_sq ) {
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

	for ( core::Size iip=1; iip<=design_pos.size(); iip++ ) {
		core::Size ir = design_pos[iip];
		if ( filter_intrabb_ ) {
			TR.Debug << "Filtering: Checking resi: " << ir << std::endl;
			contact = true;
			for ( core::Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++ ) {
				if ( !multicomponent_ ) {
					if ( sym_info->subunit_index(ir) > nsub_bblock_ || sym_info->subunit_index(jr) > nsub_bblock_ ) continue;
				} else {
					Sizes const & intra_subs(get_component_of_residue(pose,ir)=='A'?intra_subs1:intra_subs2);
					if ( get_component_of_residue(pose,ir)!=get_component_of_residue(pose,jr) ) continue;
					if ( find(intra_subs.begin(), intra_subs.end(), sym_info->subunit_index(jr)) == intra_subs.end() ) continue;
				}

				if ( sym_info->subunit_index(jr) == 1 ) continue;

				for ( core::Size ia = 1; ia<=pose.residue(ir).nheavyatoms(); ia++ ) {
					for ( core::Size ja = 1; ja<=pose.residue(jr).nheavyatoms(); ja++ ) {
						if ( pose.residue(ir).xyz(ia).distance_squared(pose.residue(jr).xyz(ja)) <= bblock_dist_sq ) {
							if ( intrabb_only_ ) { contact = false; break; }
							else {
								// However, if the residue in question is clashing badly (usually with a
								// residue from another building block), it needs to be designed.
								core::scoring::EnergyMap em1 = scored_pose.energies().residue_total_energies(ir);
								Real resi_fa_rep = em1[core::scoring::fa_rep];
								TR.Debug << "resi_fa_rep: " << resi_fa_rep << " fa_rep_cut_: " << fa_rep_cut_ << std::endl;
								if ( resi_fa_rep < fa_rep_cut_ ) { contact = false; TR.Debug << "Filtered out resi: " << ir << std::endl; break; }
							}
						}
					}
					if ( contact == false ) break;
				}
				if ( contact == false ) break;
			}
			if ( (contact && !intrabb_only_) || ((contact == false) && (intrabb_only_ == true )) ) {
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
	for ( core::Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++ ) {
		if ( filtered_design_pos.find(ir) != filtered_design_pos.end() ) {
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
	filter_intrabb_ = tag->getOption< bool >("filter_intrabb", true);
	intrabb_only_ = tag->getOption< bool >("intrabb_only", false);
	multicomponent_ = tag->getOption< bool >("multicomp", false);

	if ( tag->hasOption("bblock_reference_pdb") ) {
		std::string bblock_ref_pdb = tag->getOption< std::string >("bblock_reference_pdb");
		bblock_reference_pose_ = utility::pointer::make_shared< core::pose::Pose >();
		core::import_pose::pose_from_file( *bblock_reference_pose_, bblock_ref_pdb , core::import_pose::PDB_file);
	}
}

void BuildingBlockInterfaceOperation::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "nsub_bblock", xsct_non_negative_integer, "The number of subunits in the symmetric building block (e.g., 3 for a trimer). This option is not needed for multicomponent systems.",  "1"  )
		+ XMLSchemaAttribute::attribute_w_default(  "bblock_reference_pdb", xs_string, "A reference PDB for defining the bblock interface. NOT COMPATIBLE with multicomponent symmetry.",  "1"  )
		+ XMLSchemaAttribute::attribute_w_default(  "sym_dof_names", xs_string, "Names of the sym_dofs corresponding to the symmetric building blocks. (Eventually replace the need for this option by having is_singlecomponent or is_multicomponent utility functions). If no sym_dof_names are specified, then they will be extracted from the pose.",  "XRW TO DO"  )
		+ XMLSchemaAttribute::attribute_w_default(  "contact_dist", xsct_real, "Residues with beta carbons not within this distance of any beta carbon from another building block are prevented from repacking.",  "10.0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "bblock_dist", xsct_real, "The all-heavy atom cutoff distance used to specify residues that are making inter-subunit contacts within the building block. Because these residues are making presumably important intra-building block interactions, they are prevented from repacking unless they are clashing.",  "5.0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "fa_rep_cut", xsct_real, "The cutoff used to determine whether residues making inter-subunit contacts within the building block are clashing.",  "3.0"  )
		+ XMLSchemaAttribute::attribute_w_default(  "filter_intrabb", xsct_rosetta_bool, "Filter intra-building block interactions?",  "true"  )
		+ XMLSchemaAttribute::attribute_w_default(  "intrabb_only", xsct_rosetta_bool, "Only intra-building block residues?",  "false"  )
		+ XMLSchemaAttribute::attribute_w_default(  "multicomp", xsct_rosetta_bool, "If true, system has more than one component. If false, it is a single component system.",  "false"  );

	task_op_schema_w_attributes( xsd, keyname(), attributes, "For use when designing with symmetric building blocks. Prevents repacking at residues that are: 1) distant from the inter-building block interface, or 2) near the inter-building block interface, but also make intra-building block interface contacts that are not clashing." );
}


} //namespace matdes
} //namespace protocols
