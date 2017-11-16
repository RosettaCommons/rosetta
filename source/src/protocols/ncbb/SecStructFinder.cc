// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeSet.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/chemical/rings/RingConformerSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/ncbb/SecStructMinimizeMultiFunc.hh>
#include <protocols/ncbb/SecStructMinimizeMover.hh>
#include <protocols/ncbb/SecStructFinder.hh>
#include <protocols/ncbb/SecStructFinderCreator.hh>
#include <protocols/ncbb/util.hh>

#include <numeric/random/random.hh>

//Basic headers
#include <basic/basic.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <string>
#include <sstream>
#include <cmath>

// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace func;
using namespace constraints;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace core::pack::task;
using namespace core::id;


// tracer - used to replace cout
static basic::Tracer TR( "SecStructFinder" );

namespace protocols {
namespace ncbb {

SecStructFinder::SecStructFinder():
	residue_( "ALA" ),
	min_length_( 5 ),
	max_length_( 5 ),
	bin_size_( 10 ),
	dump_threshold_( -10000 ),
	dihedral_min_( -180 ),
	dihedral_max_(  180 ),
	dihedral_pattern_( expand_pattern_to_fit( "A", min_length_ ) ),
	alpha_beta_pattern_( expand_pattern_to_fit( "A", min_length_ ) ),
	min_everything_( false ),
	cart_( false ),
	constrain_( false ),
	dissimilarity_( 10 )
{
	Mover::type("SecStructFinder");

	score_fxn_ = get_score_function();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
}

SecStructFinder::SecStructFinder(
	std::string residue,
	Size min_length,
	Size max_length,
	Real bin_size,
	Real dissimilarity,
	Real dihedral_min,
	Real dihedral_max,
	Real dump_threshold,
	std::string dihedral_pattern,
	std::string alpha_beta_pattern,
	bool min_everything,
	bool cart,
	bool constrain
):
	residue_(std::move( residue )),
	min_length_( min_length ),
	//max_length_( max_length ),
	bin_size_( bin_size ),
	dump_threshold_( dump_threshold ),
	dihedral_min_( dihedral_min ),
	//dihedral_max_( ( dihedral_min > dihedral_max ) ? dihedral_min : dihedral_max ),
	//dihedral_pattern_( dihedral_pattern ),
	//alpha_beta_pattern_( expand_pattern_to_fit( alpha_beta_pattern, min_length_ ) ),
	min_everything_( min_everything ),
	cart_( cart ),
	constrain_( constrain ),
	dissimilarity_( ( bin_size > dissimilarity ) ? bin_size : dissimilarity )
{
	max_length_ = ( min_length_ > max_length ) ? min_length_ : max_length;
	dihedral_pattern_ = expand_pattern_to_fit( dihedral_pattern, min_length_ );
	alpha_beta_pattern_ = expand_pattern_to_fit( alpha_beta_pattern, min_length_ );
	dihedral_max_ = ( dihedral_min_ > dihedral_max ) ? dihedral_min_ : dihedral_max;

	Mover::type("SecStructFinder");

	score_fxn_ = get_score_function();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn_);
}

std::string
SecStructFinder::alpha_to_beta( std::string alpha ) {
	if ( alpha == "ALA" ) {
		return "B3A";
	}
	if ( alpha == "CYS" ) {
		return "B3C";
	}
	if ( alpha == "ASP" ) {
		return "B3D";
	}
	if ( alpha == "GLU" ) {
		return "B3E";
	}
	if ( alpha == "PHE" ) {
		return "B3F";
	}
	if ( alpha == "GLY" ) {
		return "B3G";
	}
	if ( alpha == "HIS" ) {
		return "B3H";
	}
	if ( alpha == "ILE" ) {
		return "B3I";
	}
	if ( alpha == "LYS" ) {
		return "B3K";
	}
	if ( alpha == "LEU" ) {
		return "B3L";
	}
	if ( alpha == "MET" ) {
		return "B3M";
	}
	if ( alpha == "ASN" ) {
		return "B3N";
	}
	if ( alpha == "GLN" ) {
		return "B3Q";
	}
	if ( alpha == "ARG" ) {
		return "B3R";
	}
	if ( alpha == "SER" ) {
		return "B3S";
	}
	if ( alpha == "THR" ) {
		return "B3T";
	}
	if ( alpha == "VAL" ) {
		return "B3V";
	}
	if ( alpha == "TRP" ) {
		return "B3W";
	}
	if ( alpha == "TYR" ) {
		return "B3Y";
	}
	if ( alpha == "PRO" ) {
		return "B3P";
	}

	return "B3A";
}

std::string
SecStructFinder::expand_pattern_to_fit( std::string pattern, Size length ) {
	if ( length < pattern.length() ) {
		return pattern.substr( length );
	} else if ( length == pattern.length() ) {
		return pattern;
	} else {
		std::string expanded_pattern;
		while ( length > pattern.length() ) {
			expanded_pattern += pattern;
			length -= pattern.length();
		}
		expanded_pattern += pattern.substr( 0, length );
		return expanded_pattern;
	}
}

bool
SecStructFinder::uniq_refers_to_beta (
	char uniq
) {
	for ( Size i = 1; i <= dihedral_pattern_.length(); ++i ) {
		if ( dihedral_pattern_[i] == uniq ) {
			if ( alpha_beta_pattern_[i] == 'A' || alpha_beta_pattern_[i] == 'P' ) return false;
			else return true;
		}
	}
	return false;
}

void
SecStructFinder::initialize_rtype_vector( utility::vector1< ResidueType > & restypes ) {

	core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

	for ( Size resi = 0; resi < min_length_; ++resi ) {
		std::string name;
		std::string npatch;
		std::string cpatch;

		if ( alpha_beta_pattern_[ resi ] == 'B' ) {
			name = alpha_to_beta( residue_ );
			npatch = ":AcetylatedNtermProteinFull";
			cpatch = ":MethylatedCtermProteinFull";
		} else if ( alpha_beta_pattern_[ resi ] == 'A' ) {
			name = residue_;
			npatch = ":AcetylatedNtermProteinFull";
			cpatch = ":MethylatedCtermProteinFull";
		} else {
			name = "201";
			npatch = ":AcetylatedNtermPeptoidFull";
			cpatch = ":CtermPeptoidFull"; // can't constrain dihedral with peptoid final res
		}

		if ( resi == 0 ) {
			name += npatch;
		}

		// todo
		if ( resi == min_length_-1 ) {
			name += cpatch;
		}

		restypes.push_back( residue_set_cap->name_map( name ) );
	}
}

std::string
SecStructFinder::make_filename (
	Size number_dihedrals,
	utility::vector1< Real > dihedrals
) {
	std::stringstream filename;
	filename << residue_ << "_" << dihedral_pattern_ << "_" << min_length_;
	for ( Size i = 1; i <= number_dihedrals; ++i ) {
		filename << "_" << dihedrals[i];
	}
	filename << ".pdb";
	return filename.str();
}


bool
SecStructFinder::too_similar( Size i, Size j, utility::vector1< Real > dihedrals ) {
	bool similar = false;
	if ( alpha_beta_pattern_[i] == alpha_beta_pattern_[j] ) { //  don't compare dissimilars
		if ( alpha_beta_pattern_[i] == 'A' || alpha_beta_pattern_[i] == 'P' ) {
			if ( std::abs(basic::periodic_range( dihedrals[i] - dihedrals[j], 360 )) < dissimilarity_ && std::abs(basic::periodic_range( dihedrals[i+1] - dihedrals[j+1], 360 ) ) < dissimilarity_ ) {
				similar = true;
			}
		} else if ( alpha_beta_pattern_[i] == 'B' ) {
			if ( std::abs(basic::periodic_range( dihedrals[i] - dihedrals[j], 360 )) < dissimilarity_ && std::abs(basic::periodic_range( dihedrals[i+1] - dihedrals[j+1], 360 ))  < dissimilarity_ && std::abs(basic::periodic_range( dihedrals[i+2] - dihedrals[j+2], 360 ) )  < dissimilarity_ ) {
				similar = true;
			}
		}
	}
	return similar;
}

void
SecStructFinder::show_current_dihedrals( Size number_dihedral_sets, utility::vector1< char > uniqs, utility::vector1< Real > dihedrals ) {
	TR << "Presently ";
	for ( Size i = 1; i <= number_dihedral_sets; ++i ) {
		TR << uniqs[i] << " ( ";
		Size take_from_here = give_dihedral_index( i, uniqs, dihedral_pattern_, alpha_beta_pattern_ );
		if ( uniq_refers_to_beta( uniqs[i] ) ) {
			TR << dihedrals[ take_from_here ] << ", " << dihedrals[ take_from_here+1 ] << ", " << dihedrals[ take_from_here+2 ] << " ), ";
		} else {
			TR << dihedrals[ take_from_here ] << ", " << dihedrals[ take_from_here+1 ] << " ), ";
		}
	}
}

void
SecStructFinder::apply( Pose & pose )
{
	Size number_dihedral_sets;
	utility::vector1< char > uniqs;
	count_uniq_char( dihedral_pattern_, number_dihedral_sets, uniqs );

	Size number_dihedrals = get_number_dihedrals( uniqs, dihedral_pattern_, alpha_beta_pattern_ );//2 * number_dihedral_sets; // plus one per beta AA!
	TR << "Investigating dihedral pattern " << dihedral_pattern_ << " with " << number_dihedral_sets << " uniques " << std::endl;

	utility::vector1< ResidueType > restypes;
	initialize_rtype_vector( restypes );
	TR << "Restype vector: ";
	for ( Size ii = 1; ii <= restypes.size(); ++ii ) {
		TR << restypes[ ii ].name() << ", ";
	}
	TR << std::endl;


	core::conformation::ResidueOP new_rsd( nullptr );
	new_rsd = conformation::ResidueFactory::create_residue( restypes[1] );
	pose.append_residue_by_jump( *new_rsd, 1 );
	for ( Size i = 2; i <= min_length_; ++i ) {
		new_rsd = conformation::ResidueFactory::create_residue( restypes[ i ] );
		pose.append_residue_by_bond( *new_rsd, true );
	}

	kinematics::MoveMapOP min_mm( new kinematics::MoveMap );
	min_mm->set_bb( true );
	min_mm->set_chi( true );

	protocols::simple_moves::MinMover minmover( min_mm, score_fxn_, "lbfgs_armijo_nonmonotone", 0.0001, true );
	minmover.cartesian( cart_ );
	utility::vector1< Pose > poses_to_min;
	utility::vector1< utility::vector1< Real > > dihedrals_to_min;
	utility::vector1< std::string > minned_filenames;

	// now we have our pose
	utility::vector1< Real > dihedrals( number_dihedrals+1, dihedral_min_ );

	// while loop idiom for looping through every layer of a vector of indices
	Size p = 1;
	while ( dihedrals[ number_dihedrals+1 ] == dihedral_min_ ) {

		// skip if the dihedrals for A are the same as B, etc.
		// break after first example found
		bool skip = false;
		for ( Size i = 1; i <= number_dihedral_sets; i += 2 ) {
			for ( Size j = i+1; j <= number_dihedral_sets; j++ ) {
				skip = too_similar( i, j, dihedrals );
				if ( skip ) break;
			}
			if ( skip ) break;
		}

		if ( !skip ) {

			TR << "Trying dihedral DOF vector: ";
			for ( Size ii = 1; ii <= number_dihedrals; ++ii ) {
				TR << dihedrals[ii] << ", ";
			} TR << std::endl;

			for ( Size resi = 1; resi <= pose.size(); ++resi ) {
				Size take_from_here = 1;
				for ( Size i = 1; i <= number_dihedral_sets; ++i ) {
					if ( dihedral_pattern_[ resi-1 ] == uniqs[ i ] ) { // -1 to sync indexing
						take_from_here = give_dihedral_index( i, uniqs, dihedral_pattern_, alpha_beta_pattern_  );
					}
				}

				//TR << "For res " << resi << " ( " << pose.residue_type( resi ).name() << " )... " << std::endl;
				if ( pose.residue( resi ).type().is_beta_aa() ) {
					//TR << "Assigning DOF vector items " << take_from_here << " to " << take_from_here+2 << " to this BAA" << std::endl;
					pose.set_torsion( TorsionID( resi, BB, 1 ), dihedrals[ take_from_here ] );
					pose.set_torsion( TorsionID( resi, BB, 2 ), dihedrals[ take_from_here+1 ] );
					pose.set_torsion( TorsionID( resi, BB, 3 ), dihedrals[ take_from_here+2 ] );
					pose.set_torsion( TorsionID( resi, BB, 4 ), 180 );
				} else {
					//TR << "Assigning DOF vector items " << take_from_here << " to " << take_from_here+1 << " to this AAA" << std::endl;
					pose.set_phi( resi, dihedrals[ take_from_here ] ); //phi_vec[ vec_index ] );
					pose.set_psi( resi, dihedrals[ take_from_here+1 ] ); //psi_vec[ vec_index ] );
					pose.set_omega( resi, 180 );
				}
			}

			show_current_dihedrals( number_dihedral_sets, uniqs, dihedrals );

			Real score = ( *score_fxn_ ) ( pose );
			TR << "score is " << score << std::endl;

			if ( min_everything_ ) {

				if ( score_fxn_->has_zero_weight( core::scoring::dihedral_constraint ) ) {
					score_fxn_->set_weight( core::scoring::dihedral_constraint, 1.0 );
				}

				Pose minpose = pose;

				SecStructMinimizeMoverOP ssmm( new SecStructMinimizeMover() );
				ssmm->set_constrain( constrain_ );
				ssmm->set_scorefunction( score_fxn_ );
				ssmm->set_alpha_beta_pattern( alpha_beta_pattern_ );
				ssmm->set_dihedral_pattern( dihedral_pattern_ );
				ssmm->set_dihedrals( dihedrals );
				ssmm->apply( minpose );

				score_fxn_->set_weight( core::scoring::dihedral_constraint, 0.0 );
				Real score = ( ( *score_fxn_ ) ( minpose ) );

				TR << " and minned is " << score << std::endl;
				if ( score <= dump_threshold_ ) {
					TR << "Thus, dumping ";
					std::string filename = make_filename( number_dihedrals, dihedrals );

					minpose.dump_scored_pdb( filename.c_str(), ( *score_fxn_ ) );
					if ( minpose.residue_type( 1 ).is_beta_aa() ) {
						TR << " ( " << minpose.torsion( TorsionID( 1, BB, 1 ) ) << ", "<< minpose.torsion( TorsionID( 1, BB, 2 ) ) << ", " << minpose.torsion( TorsionID( 1, BB, 3 ) ) << " ), ";
					} else {
						TR << " ( " << minpose.torsion( TorsionID( 1, BB, 1 ) ) << ", " << minpose.psi( 1 ) << " ), ";
					}

					for ( Size resi = 2; resi < pose.size(); ++resi ) {
						if ( minpose.residue_type( resi ).is_beta_aa() ) {
							TR << " ( " << minpose.torsion( TorsionID( resi, BB, 1 ) ) << ", "<< minpose.torsion( TorsionID( resi, BB, 2 ) ) << ", " << minpose.torsion( TorsionID( resi, BB, 3 ) ) << " ), ";
						} else {
							TR << " ( " << minpose.phi( resi ) << ", " << minpose.psi( resi ) << " ), ";
						}
					}

					if ( minpose.residue_type( pose.size() ).is_beta_aa() ) {
						TR << " ( " << minpose.torsion( TorsionID(pose.size(), BB, 1 ) ) << ", "<< minpose.torsion( TorsionID( pose.size(), BB, 2 ) ) << ", "
							<< minpose.torsion( TorsionID( pose.size(), BB, 3 ) ) << " ), ";
					} else {
						TR << " ( " << minpose.phi( pose.size() ) << ", " << minpose.torsion( TorsionID( pose.size(), BB, 2 ) ) << " ), " << std::endl;
					}
				}
			} else {
				if ( score <= dump_threshold_ ) {
					TR << " thus, dumping ";
					std::string filename = make_filename( number_dihedrals, dihedrals );

					pose.dump_scored_pdb( filename.c_str(), ( *score_fxn_ ) );
					TR << filename.c_str() << " " << score << std::endl;
					poses_to_min.push_back( pose );
					utility::vector1< Real > dihedrals_for_minimization;
					for ( Size index = 1; index <= number_dihedrals; ++index ) {
						dihedrals_for_minimization.push_back( dihedrals[ index ] );
					}
					dihedrals_to_min.push_back( dihedrals_for_minimization );
					minned_filenames.push_back( filename );
				}
			}
		}

		dihedrals[ 1 ] += bin_size_;
		while ( dihedrals[ p ] == dihedral_max_ ) {
			dihedrals[ p ] = dihedral_min_;
			dihedrals[ ++p ] += bin_size_;
			if ( dihedrals[ p ] != dihedral_max_ ) p = 1;
		}
	}

	// now min all our poses to min
	if ( poses_to_min.size() == 0 ) return;

	for ( Size ii = 1, sz = poses_to_min.size(); ii < sz; ++ii ) {

		Pose minpose = poses_to_min[ii];

		SecStructMinimizeMoverOP ssmm( new SecStructMinimizeMover() );
		ssmm->set_constrain( constrain_ );
		ssmm->set_scorefunction( score_fxn_ );
		ssmm->set_alpha_beta_pattern( alpha_beta_pattern_ );
		ssmm->set_dihedral_pattern( dihedral_pattern_ );
		ssmm->set_dihedrals( dihedrals );
		ssmm->apply( minpose );

		score_fxn_->set_weight( core::scoring::dihedral_constraint, 0.0 );
		Real score = ( ( *score_fxn_ ) ( minpose ) );
		TR << "for stored pose with dihedrals ";
		for ( Size resi = 1, n = poses_to_min[ii].size(); resi <= n; ++resi ) {
			if ( alpha_beta_pattern_[ resi-1 ] == 'A' || alpha_beta_pattern_[ resi-1 ] == 'P' ) {
				TR << "( " << poses_to_min[ii].phi( resi ) << ", " << poses_to_min[ii].psi( resi ) << " ), ";
			} else {
				Real phi = poses_to_min[ii].torsion( TorsionID( resi, BB, 1 )  );
				Real tht = poses_to_min[ii].torsion( TorsionID( resi, BB, 2 )  );
				Real psi = poses_to_min[ii].torsion( TorsionID( resi, BB, 3 )  );
				TR << "( " << phi << ", " << tht << ", " << psi << " ), ";
			}
		}
		TR << std::endl << " now ";
		for ( Size resi = 1, n = minpose.size(); resi <= n; ++resi ) {
			if ( alpha_beta_pattern_[ resi-1 ] == 'A' || alpha_beta_pattern_[ resi-1 ] == 'P' ) {
				TR << "( " << minpose.phi( resi ) << ", " << minpose.psi( resi ) << " ), ";
			} else {
				Real phi = minpose.torsion( TorsionID( resi, BB, 1 )  );
				Real tht = minpose.torsion( TorsionID( resi, BB, 2 )  );
				Real psi = minpose.torsion( TorsionID( resi, BB, 3 )  );
				TR << "( " << phi << ", " << tht << ", " << psi << " ), ";
			}
		}

		TR << "minned score is " << score << std::endl;
		minned_filenames[ii] = "minned_" + minned_filenames[ii];
		minpose.dump_scored_pdb( minned_filenames[ii].c_str(), ( * score_fxn_ ) );
	}
}//apply

protocols::moves::MoverOP SecStructFinder::clone() const {
	return SecStructFinderOP( new SecStructFinder(
		residue_,
		min_length_,
		max_length_,
		bin_size_,
		dissimilarity_,
		dihedral_min_,
		dihedral_max_,
		dump_threshold_,
		dihedral_pattern_,
		alpha_beta_pattern_,
		min_everything_,
		cart_ ) );
}

void SecStructFinder::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & ) {

	if ( tag->hasOption( "residue" ) ) {
		this->residue_ = tag->getOption<std::string>("residue", residue_);
	} else {
		residue_ = "ALA";
	}

	if ( tag->hasOption( "min_length" ) ) {
		this->min_length_ = tag->getOption<Size>("min_length", min_length_);
	} else {
		min_length_ = 5;
	}

	if ( tag->hasOption( "max_length" ) ) {
		this->max_length_ = tag->getOption<Size>("max_length", max_length_);
	} else {
		max_length_ = 5;
	}

	if ( max_length_ < min_length_ ) max_length_ = min_length_;

	if ( tag->hasOption( "bin_size" ) ) {
		this->bin_size_ = tag->getOption<Real>("bin_size", bin_size_);
	} else {
		bin_size_ = 10;
	}

	if ( tag->hasOption( "dissimilarity" ) ) {
		this->dissimilarity_ = tag->getOption<Real>("dissimilarity", dissimilarity_);
	} else {
		dissimilarity_ = 10;
	}

	if ( dissimilarity_ < bin_size_ ) dissimilarity_ = bin_size_;

	if ( tag->hasOption( "dump_threshold" ) ) {
		this->dump_threshold_ = tag->getOption<Real>("dump_threshold", dump_threshold_);
	} else {
		dump_threshold_ = 10;
	}

	if ( tag->hasOption( "dihedral_pattern" ) ) {
		this->dihedral_pattern_ = tag->getOption<std::string>("dihedral_pattern", dihedral_pattern_);
	} else {
		dihedral_pattern_ = "A";
	}

	dihedral_pattern_ = expand_pattern_to_fit( dihedral_pattern_, min_length_ );

	if ( tag->hasOption( "alpha_beta_pattern" ) ) {
		this->alpha_beta_pattern_ = tag->getOption<std::string>("alpha_beta_pattern", alpha_beta_pattern_);
	} else {
		alpha_beta_pattern_ = "A";
	}

	alpha_beta_pattern_ = expand_pattern_to_fit( alpha_beta_pattern_, min_length_ );

	if ( tag->hasOption( "min_everything" ) ) {
		this->min_everything_ = tag->getOption<bool>("min_everything", min_everything_);
	} else {
		min_everything_ = false;
	}

	if ( tag->hasOption( "cart" ) ) {
		this->cart_ = tag->getOption<bool>("cart", cart_);
	} else {
		cart_ = false;
	}

}

// MoverCreator
moves::MoverOP SecStructFinderCreator::create_mover() const {
	return SecStructFinderOP( new SecStructFinder() );
}

std::string SecStructFinderCreator::keyname() const {
	return SecStructFinderCreator::mover_name();
}

std::string SecStructFinderCreator::mover_name(){
	return "SecStructFinder";
}

}
}
