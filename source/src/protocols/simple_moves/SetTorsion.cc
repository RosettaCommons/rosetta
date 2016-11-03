// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SetTorsion.c
/// @brief Sets the value of a desired torsion.
/// @author Modified 4 June 2015 by Vikram K. Mulligan (vmullig@uw.edu), Baker Laboratory,
/// to add perturb torsion option (I didn't write this file, though).

// Unit headers
#include <protocols/simple_moves/SetTorsion.hh>
#include <protocols/simple_moves/SetTorsionCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
//parsing
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/filters/Filter.fwd.hh> //Filters_map
#include <protocols/rosetta_scripts/util.hh>
#include <basic/Tracer.hh>

#include <core/pose/selection.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>

#include <boost/lexical_cast.hpp>


// Utility Headers

// Unit Headers

// C++ headers


namespace protocols {
namespace simple_moves {

using namespace core;
using namespace core::chemical;
using namespace std;
using namespace numeric::conversions;

using core::pose::Pose;
using core::conformation::Residue;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.SetTorsion" );

std::string
SetTorsionCreator::keyname() const
{
	return SetTorsionCreator::mover_name();
}

protocols::moves::MoverOP
SetTorsionCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetTorsion );
}

std::string
SetTorsionCreator::mover_name()
{
	return "SetTorsion";
}

SetTorsion::~SetTorsion() = default;

/// @brief default ctor
SetTorsion::SetTorsion() :
	parent(),
	random_set_(false),
	angle_(),
	residues_(),
	torsion_name_(),
	custom_rama_map_(),
	extending_(),
	torsion_atoms_(),
	perturbation_type_(),
	perturbation_magnitude_(),
	fold_tree_root_(0)
{}

/// @brief Actually get the value that the torsion will be set to.
/// @details Depending on settings, this will look up a value, generate a random value, or perturb an input value.
core::Real
SetTorsion::angle(
	core::Size const iset,
	core::Real const &old_angle
) const {
	if ( angle_[iset] == "random" ) {
		return 360.0*numeric::random::rg().uniform()-180.0;
	} else if ( angle_[iset] == "perturb" ) {
		core::Real returnval( old_angle );
		if ( perturbation_type(iset)==perturbtorsion_uniform ) {
			returnval += (numeric::random::rg().uniform()-0.5)*perturbation_magnitude(iset);
		} else if ( perturbation_type(iset)==perturbtorsion_gaussian ) {
			returnval += numeric::random::rg().gaussian()*perturbation_magnitude(iset);
		} else {
			utility_exit_with_message("Error in protocols::simple_moves::SetTorsion::angle(): Perturbation type not recognized!");
		}
		return returnval;
	} else {
		return boost::lexical_cast<core::Real>(angle_[iset]);
	}
}

utility::vector1<core::Size> SetTorsion::residue_list(core::Size iset, core::pose::Pose const & pose) {
	utility::vector1<core::Size> residue_numbers;
	if ( residues_[iset] == "ALL" ) {
		for ( Size ires=1; ires<=pose.size(); ++ires ) {
			residue_numbers.push_back(ires);
		}
	} else if ( residues_[iset] == "pick_atoms" ) {
		// shouldn't get here, do nothing
	} else if ( residues_[iset] == "random" ) {
		Size res_num = numeric::random::rg().random_range(1, pose.size());
		residue_numbers.push_back(res_num);
		if ( extending_[iset] != 0 ) {
			for ( int ishift = 1; ishift <= (int)extending_[iset]; ++ishift ) {
				if ( res_num - ishift >= 1 ) {
					residue_numbers.push_back(res_num - ishift);
				}
				if ( res_num + ishift <= pose.size() ) {
					residue_numbers.push_back(res_num + ishift);
				}
			}
		}
	} else {
		utility::vector1<std::string> buff = utility::string_split( residues_[iset], ',' );
		for ( std::string const & field : buff ) {
			Size const value = std::atoi( field.c_str() );
			residue_numbers.push_back(value);
		}
	}

	return residue_numbers;
}
void SetTorsion::apply( Pose & pose ) {

	//Store the old FoldTree in case we're rerooting:
	core::kinematics::FoldTree saved_ft( pose.fold_tree() );

	//Reroot the FoldTree:
	if ( fold_tree_root_ > 0 && fold_tree_root_ <=pose.size() ) {
		TR << "Old FoldTree:" << std::endl;
		pose.fold_tree().show(TR);
		core::kinematics::FoldTree ft_copy;
		ft_copy.clear();
		core::Size jumpcount=1;
		core::Size const chaincount = pose.conformation().num_chains();
		bool firstchainrerooted = false;
		for ( core::Size chain=1; chain<=chaincount; ++chain ) { //Loop through all chains
			if ( fold_tree_root_ <= pose.conformation().chain_end(chain) && fold_tree_root_ >= pose.conformation().chain_begin(chain) ) { //If the new root is in the current chain
				if ( chain>1 ) {
					ft_copy.add_edge( pose.conformation().chain_begin(1), fold_tree_root_, jumpcount++);
					firstchainrerooted=false;
				} else {
					firstchainrerooted=true;
				}
				if ( fold_tree_root_ < pose.conformation().chain_end(chain) ) ft_copy.add_edge( fold_tree_root_, pose.conformation().chain_end(chain), -1);
				if ( fold_tree_root_ > pose.conformation().chain_begin(chain) ) ft_copy.add_edge( fold_tree_root_, pose.conformation().chain_begin(chain), -1);
			} else { //Otherwise, just add the chain.
				if ( chain>1 ) ft_copy.add_edge( (firstchainrerooted ? fold_tree_root_ : pose.conformation().chain_begin(1) ) , pose.conformation().chain_begin(chain), jumpcount++);
				ft_copy.add_edge( pose.conformation().chain_begin(chain), pose.conformation().chain_end(chain), -1);
			}
		}
		//ft_copy.delete_self_edges();
		if ( firstchainrerooted ) ft_copy.reorder(fold_tree_root_);
		else ft_copy.reorder(1);
		TR << "New FoldTree:" << std::endl;
		ft_copy.show(TR);
		pose.fold_tree( ft_copy );
		TR.flush();
	}

	core::scoring::Ramachandran & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran_nonconst(); //Must be nonconst to allow lazy loading of custom rama tables.
	Size picked_set(1);
	if ( random_set_ ) {
		picked_set = numeric::random::rg().random_range(1, n_torsion_sets());
	}
	for ( core::Size iset=1; iset<=n_torsion_sets(); ++iset ) {
		if ( random_set_ ) {
			if ( iset != picked_set ) continue;
		}
		if ( residues_[iset] == "pick_atoms" ) {
			core::Real angle_in( //The original torsion angle value
				pose.conformation().torsion_angle(
				id::AtomID(pose.residue(torsion_atoms_[iset][1].rsd()).atom_index(torsion_atoms_[iset][1].atom()), torsion_atoms_[iset][1].rsd()),
				id::AtomID(pose.residue(torsion_atoms_[iset][2].rsd()).atom_index(torsion_atoms_[iset][2].atom()), torsion_atoms_[iset][2].rsd()),
				id::AtomID(pose.residue(torsion_atoms_[iset][3].rsd()).atom_index(torsion_atoms_[iset][3].atom()), torsion_atoms_[iset][3].rsd()),
				id::AtomID(pose.residue(torsion_atoms_[iset][4].rsd()).atom_index(torsion_atoms_[iset][4].atom()), torsion_atoms_[iset][4].rsd())
				)
			);
			pose.conformation().set_torsion_angle(
				id::AtomID(pose.residue(torsion_atoms_[iset][1].rsd()).atom_index(torsion_atoms_[iset][1].atom()), torsion_atoms_[iset][1].rsd()),
				id::AtomID(pose.residue(torsion_atoms_[iset][2].rsd()).atom_index(torsion_atoms_[iset][2].atom()), torsion_atoms_[iset][2].rsd()),
				id::AtomID(pose.residue(torsion_atoms_[iset][3].rsd()).atom_index(torsion_atoms_[iset][3].atom()), torsion_atoms_[iset][3].rsd()),
				id::AtomID(pose.residue(torsion_atoms_[iset][4].rsd()).atom_index(torsion_atoms_[iset][4].atom()), torsion_atoms_[iset][4].rsd()),
				radians( angle(iset, angle_in) )
			);
		} else {
			for ( core::Size ires=1; ires<=residue_list(iset, pose).size(); ++ires ) {
				Size resnum = residue_list(iset, pose)[ires];

				if ( torsion_name(iset) == "phi" ) {
					if ( pose.residue(resnum).type().is_alpha_aa()
							|| pose.residue(resnum).type().is_beta_aa()
							|| pose.residue(resnum).type().is_peptoid() ) {
						pose.set_phi( resnum, angle(iset, pose.phi(resnum)) );
					}
				} else if ( torsion_name(iset) == "theta" ) {
					if ( pose.residue(resnum).type().is_beta_aa() ) {
						pose.set_theta( resnum, angle(iset, pose.theta(resnum)) );
					}
				} else if ( torsion_name(iset) == "psi" ) {
					if ( pose.residue(resnum).type().is_alpha_aa()
							|| pose.residue(resnum).type().is_beta_aa()
							|| pose.residue(resnum).type().is_peptoid() ) {
						pose.set_psi( resnum, angle(iset, pose.psi(resnum)) );
					}
				} else if ( torsion_name(iset) == "omega" ) {
					if ( pose.residue(resnum).type().is_alpha_aa()
							|| pose.residue(resnum).type().is_beta_aa()
							|| pose.residue(resnum).type().is_peptoid() ) {
						pose.set_omega( resnum, angle(iset, pose.omega(resnum)) );
					} //TODO -- beta-amino acids.
				} else if ( torsion_name(iset) == "rama" ) {
					if ( pose.residue(resnum).type().is_alpha_aa() ) {
						if ( angle_[iset] == "rama_biased" ) {
							Real phi(0), psi(0);
							if ( custom_rama_map_[iset] == "" ) {
								assert(pose.residue(resnum).aa()!=core::chemical::aa_unk);
								rama.random_phipsi_from_rama( pose.residue_type(resnum).aa(), phi, psi); //TODO -- use backbone_aa
							} else { // custom_ramna_map_[iset] != ""
								rama.draw_random_phi_psi_from_extra_cdf( rama.get_ramatable_type_by_name(custom_rama_map_[iset]) , phi, psi);
							}
							pose.set_phi( resnum, phi );
							pose.set_psi( resnum, psi );
						} else {
							pose.set_phi( resnum, angle(iset, pose.phi(resnum)) );
							pose.set_psi( resnum, angle(iset, pose.psi(resnum)) );
							if ( pose.residue(resnum).type().is_beta_aa() ) pose.set_theta( resnum, angle(iset, pose.theta(resnum)) );
						}
					}
				}
				TR.Debug <<"Set " << resnum <<"'s "<<torsion_name(iset)<<" to " << angle_[iset]<<std::endl;
			}
		}
	}

	if ( fold_tree_root_ > 0 ) pose.fold_tree( saved_ft ); //Reset the foldtree

	pose.update_residue_neighbors();
}

std::string
SetTorsion::get_name() const {
	return SetTorsionCreator::mover_name();
}

void SetTorsion::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & //pose
)
{
	utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
	utility::vector1< utility::tag::TagCOP >::const_iterator tag_it;

	random_set_ = tag->getOption< bool >( "random", false );

	fold_tree_root_=tag->getOption< core::Size >( "foldtree_root", 0 ); //Get the residue index for the fold tree root for this operation.

	for ( tag_it = branch_tags.begin(); tag_it != branch_tags.end(); ++tag_it ) {
		if ( (*tag_it)->getName() == "Torsion" ) {
			utility::vector1< id::NamedAtomID > atoms;
			atoms.resize(4, id::BOGUS_NAMED_ATOM_ID);

			angle_.push_back((*tag_it)->getOption< std::string >( "angle" ));
			add_perturbation_type( (*tag_it)->getOption< std::string >("perturbation_type", "gaussian") );
			add_perturbation_magnitude( (*tag_it)->getOption< core::Real >("perturbation_magnitude", 1.0) );
			residues_.push_back((*tag_it)->getOption< std::string >( "residue" ));
			torsion_name_.push_back((*tag_it)->getOption< std::string >( "torsion_name", ""));
			custom_rama_map_.push_back( (*tag_it)->getOption< std::string >("custom_rama_table", "") );
			extending_.push_back((*tag_it)->getOption< Size >( "extending", 0 )); // expanding picked residue

			utility::vector1< utility::tag::TagCOP > const sub_branch_tags( (*tag_it)->getTags() );
			utility::vector1< utility::tag::TagCOP >::const_iterator sub_tag_it;

			for ( sub_tag_it = sub_branch_tags.begin(); sub_tag_it != sub_branch_tags.end(); ++sub_tag_it ) {
				if ( (*sub_tag_it)->getName() == "Atom1" ) {
					atoms[1] = id::NamedAtomID( (*sub_tag_it)->getOption< string >( "atom" ), (*sub_tag_it)->getOption< Size >( "residue" ) );
				}
				if ( (*sub_tag_it)->getName() == "Atom2" ) {
					atoms[2] = id::NamedAtomID( (*sub_tag_it)->getOption< string >( "atom" ), (*sub_tag_it)->getOption< Size >( "residue" ) );
				}
				if ( (*sub_tag_it)->getName() == "Atom3" ) {
					atoms[3] = id::NamedAtomID( (*sub_tag_it)->getOption< string >( "atom" ), (*sub_tag_it)->getOption< Size >( "residue" ) );
				}
				if ( (*sub_tag_it)->getName() == "Atom4" ) {
					atoms[4] = id::NamedAtomID( (*sub_tag_it)->getOption< string >( "atom" ), (*sub_tag_it)->getOption< Size >( "residue" ) );
				}
			}
			torsion_atoms_.push_back(atoms);
		}
	}
}


} // moves
} // protocols
