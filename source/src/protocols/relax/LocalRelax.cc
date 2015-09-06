// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/relax/LocalRelax.cc
/// @brief A relax protocol that iteratively cart relaxes clustered subsets of residues
/// @author Frank DiMaio


#include <protocols/relax/LocalRelax.hh>
#include <protocols/relax/LocalRelaxCreator.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/util.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/CrystInfo.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/selection.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/constraints/util.hh>
#include <core/types.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/constants.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <basic/datacache/DataMap.hh>

#include <boost/foreach.hpp>
#define foreach_ BOOST_FOREACH

// option includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <fstream>
#include <iostream>
#include <math.h>

#include <sstream>
#include <string>
#include <queue>


namespace protocols {
namespace relax {


using namespace core;

static basic::Tracer TR("LocalRelax");


std::string
LocalRelaxCreator::keyname() const {
	return LocalRelaxCreator::mover_name();
}

protocols::moves::MoverOP
LocalRelaxCreator::create_mover() const {
	return protocols::moves::MoverOP( new LocalRelax() );
}

std::string
LocalRelaxCreator::mover_name() {
	return "LocalRelax";
}


LocalRelax::LocalRelax() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	NCYC_ = option[ basic::options::OptionKeys::relax::default_repeats ](); // n relax cycles
	NEXP_ = 2; // n expansions
	K_ = 16; // CB dist cut
	max_iter_ = 200;
	verbose_ = false;
	ramp_cart_ = false;

	ramp_schedule_.push_back(0.02);
	ramp_schedule_.push_back(0.25);
	ramp_schedule_.push_back(0.55);
	ramp_schedule_.push_back(1.0);

	pack_sfxn_ = core::scoring::get_score_function();
	min_sfxn_ = core::scoring::get_score_function();

	if ( option[ edensity::mapfile ].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *pack_sfxn_ );
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *min_sfxn_ );
	}

	if ( option[ OptionKeys::constraints::cst_fa_weight].user() ) {
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *pack_sfxn_ );
		core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *min_sfxn_ );
	} else {
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *pack_sfxn_ );
		core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *min_sfxn_ );
	}
}


void
LocalRelax::optimization_loop(
	Pose & pose,
	core::pack::task::PackerTaskOP ptask,
	core::kinematics::MoveMapOP mm,
	core::Real fa_rep_scale,
	core::Real min_tol)
{
	using namespace core::scoring;

	// minpack+cartmin
	core::optimization::CartesianMinimizer minimizer;

	core::scoring::ScoreFunctionOP local_pack_sf = pack_sfxn_->clone();
	core::scoring::ScoreFunctionOP local_min_sf = min_sfxn_->clone();

	local_pack_sf->set_weight( fa_rep, fa_rep_scale*pack_sfxn_->get_weight( fa_rep ) );
	local_min_sf->set_weight( fa_rep, fa_rep_scale*min_sfxn_->get_weight( fa_rep ) );

	if ( ramp_cart_ ) {
		core::Real cart_scale = std::max( fa_rep_scale, 0.1 );
		local_pack_sf->set_weight( cart_bonded, cart_scale*pack_sfxn_->get_weight( cart_bonded ) );
		local_min_sf->set_weight( cart_bonded, cart_scale*min_sfxn_->get_weight( cart_bonded ) );
		local_pack_sf->set_weight( cart_bonded_angle, cart_scale*pack_sfxn_->get_weight( cart_bonded_angle ) );
		local_min_sf->set_weight( cart_bonded_angle, cart_scale*min_sfxn_->get_weight( cart_bonded_angle ) );
		local_pack_sf->set_weight( cart_bonded_length, cart_scale*pack_sfxn_->get_weight( cart_bonded_length ) );
		local_min_sf->set_weight( cart_bonded_length, cart_scale*min_sfxn_->get_weight( cart_bonded_length ) );
		local_pack_sf->set_weight( cart_bonded_torsion, cart_scale*pack_sfxn_->get_weight( cart_bonded_torsion ) );
		local_min_sf->set_weight( cart_bonded_torsion, cart_scale*min_sfxn_->get_weight( cart_bonded_torsion ) );

	}

	core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", min_tol, true, false, false );
	options.max_iter(max_iter_);
	core::pack::pack_rotamers( pose, *local_pack_sf, ptask );
	minimizer.run( pose, *mm, *local_min_sf, options );

	// AMW: cppcheck flags this but unnecessarily so
	static int dump_idx=1;
	if ( verbose_ ) {
		std::string name = "opt_"+utility::to_string( dump_idx++ )+".pdb";
		TR << "Write " << name << std::endl;
		pose.dump_pdb( name );
	}
}

void
LocalRelax::get_neighbor_graph(
	Pose & pose,
	utility::vector1< utility::vector1<bool> > &neighbor) {
	using namespace core;
	using namespace core::scoring;


	// grab symminfo (if defined) from the pose
	core::conformation::symmetry::SymmetryInfoCOP symminfo=NULL;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation &>(pose.conformation()).Symmetry_Info();
	}

	core::Real K=K_;
	core::Real b=0.28;

	core::Size nres = pose.total_residue();

	// make pose polyA
	Pose pose_working = pose;
	utility::vector1< Size > protein_residues;
	for ( Size i=1, i_end = nres; i<= i_end; ++i ) {
		if ( pose.residue( i ).is_protein() ) {
			protein_residues.push_back( i );
		}
	}
	protocols::toolbox::pose_manipulation::construct_poly_XXX_pose( "ALA", pose_working, protein_residues, false, false, true );

	neighbor.clear();
	neighbor.resize( nres );

	for ( Size i=1, i_end = nres; i<= i_end; ++i ) {
		if ( symminfo && !symminfo->bb_is_independent( i ) ) continue;
		neighbor[i].resize(nres, false);
		neighbor[i][i] = true;

		conformation::Residue const & rsd1( pose_working.residue( i ) );
		for ( Size j=1, j_end = nres; j<= j_end; ++j ) {
			conformation::Residue const & rsd2( pose_working.residue( j ) );

			if ( i==j ) continue;
			if ( !rsd1.is_protein() || !rsd2.is_protein() ) continue;

			core::Real dist = (rsd1.atom("CB").xyz() - rsd2.atom("CB").xyz()).length();
			core::Real angle1 = numeric::angle_degrees(rsd1.atom("CA").xyz(), rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz() ) ;
			core::Real angle2 = numeric::angle_degrees(rsd1.atom("CB").xyz(), rsd2.atom("CB").xyz(), rsd2.atom("CA").xyz() ) ;

			core::Real angle_tgt = K*exp(b*dist);

			if ( angle_tgt < 180 && angle1 > angle_tgt && angle2 > angle_tgt ) {
				core::Size j_asu=j;
				if ( symminfo && !symminfo->bb_is_independent( j ) ) {
					j_asu = symminfo->bb_follows( j );
				}
				neighbor[i][j_asu] = true;
			}
		}
	}
}


void
LocalRelax::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/
) {
	using namespace basic::options;
	using namespace core::scoring;

	// scorefxns
	if ( tag->hasOption( "scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "scorefxn" ) );
		pack_sfxn_ = (data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_name ));
		min_sfxn_ = (data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_name ));
	}
	if ( tag->hasOption( "pack_scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "pack_scorefxn" ) );
		pack_sfxn_ = (data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_name ));
	}
	if ( tag->hasOption( "min_scorefxn" ) ) {
		std::string const scorefxn_name( tag->getOption<std::string>( "min_scorefxn" ) );
		min_sfxn_ = (data.get_ptr< ScoreFunction >( "scorefxns", scorefxn_name ));
	}

	if ( tag->hasOption( "ncyc" ) ) {
		NCYC_ = tag->getOption< int >( "ncyc" );
	}
	if ( tag->hasOption( "nexp" ) ) {
		NEXP_ = tag->getOption< int >( "nexp" );
	}
	if ( tag->hasOption( "K" ) ) {
		K_ = tag->getOption< int >( "K" );
	}

	if ( tag->hasOption( "max_iter" ) ) {
		max_iter_ = tag->getOption< int >( "max_iter" );
	}


	if ( tag->hasOption( "ramp_schedule" ) ) {
		std::string ramp_schedule_str = tag->getOption< std::string >( "ramp_schedule" );
		utility::vector1< std::string > ramp_schedule_strs = utility::string_split ( ramp_schedule_str, ',' );
		ramp_schedule_.clear();
		for ( core::Size i=1; i<= ramp_schedule_strs.size(); ++i ) {
			ramp_schedule_.push_back( atoi(ramp_schedule_strs[i].c_str()) );
		}
		runtime_assert( ramp_schedule_.size() >= 1);
	}

	verbose_ = tag->getOption< bool >( "verbose" , false );
	ramp_cart_ = tag->getOption< bool >( "ramp_cart" , false );
}


void
LocalRelax::apply( core::pose::Pose & pose) {
	core::Size nres = pose.total_residue();
	core::Size nres_asu = nres;

	// set up symm
	core::conformation::symmetry::SymmetryInfoCOP symminfo=NULL;
	if ( core::pose::symmetry::is_symmetric(pose) )  {
		symminfo = dynamic_cast<const core::conformation::symmetry::SymmetricConformation &>(pose.conformation()).Symmetry_Info();
		nres_asu = symminfo->num_independent_residues();
	}

	// set up packer task [to do: make this RS selectible]
	core::pack::task::TaskFactoryOP task ( new core::pack::task::TaskFactory );
	task->push_back( core::pack::task::operation::TaskOperationOP( new core::pack::task::operation::InitializeFromCommandline ) );
	task->push_back( core::pack::task::operation::TaskOperationOP( new core::pack::task::operation::RestrictToRepacking ) );
	task->push_back( core::pack::task::operation::TaskOperationOP( new core::pack::task::operation::IncludeCurrent ) );
	core::pack::task::PackerTaskOP ptask_resfile = task->create_task_and_apply_taskoperations( pose );

	// for each residue
	utility::vector1< utility::vector1<bool> > neighbor;
	for ( Size cyc = 1; cyc <= NCYC_; ++cyc ) {
		for ( Size innercyc = 1; innercyc <= ramp_schedule_.size(); ++innercyc ) {
			get_neighbor_graph( pose, neighbor );

			// "priority list" on residues
			//   - sort by connectedness
			utility::vector1< core::Size > neighborcounts(nres, 0);
			for ( Size i=1, i_end = nres; i<= i_end; ++i ) {
				if ( !symminfo || symminfo->bb_is_independent(i) ) {
					for ( Size j=1; j<=nres; ++j ) if ( neighbor[i][j] ) neighborcounts[j]++;
				}
			}

			utility::vector1<bool> shell0, shell1, visited(nres_asu, false);

			// mark non-packable as visited
			for ( Size i=1; i<=nres; ++i ) {
				if ( !ptask_resfile->pack_residue( i ) ) visited[i] = true;
			}

			// main loop
			while ( true ) {

				// find most connected residue
				Size maxneighb=0;
				Size currres=0;
				for ( Size i=1; i<=nres; ++i ) {
					if ( !visited[i] && neighborcounts[i] > maxneighb ) {
						maxneighb = neighborcounts[i];
						currres = i;
					}
				}

				if ( maxneighb==0 ) {
					// all done
					break;
				} else if ( maxneighb < 2 ) {
					TR << "PACK SURFACE" << std::endl;

					utility::vector1<bool> neigh_merge(nres, false);
					for ( Size i=1; i<=nres; ++i ) {
						if ( visited[i] ) continue;
						if ( symminfo && !symminfo->bb_is_independent(i) ) continue;
						for ( Size j=1; j<=nres; ++j ) {
							neigh_merge[j] = neigh_merge[j] || neighbor[i][j];
						}
					}

					// "surface pack" << generally lots of surface residues in small clusters.  Pack them all at once
					shell0 = neigh_merge;
					shell1 = shell0;
					for ( Size j=1; j<=nres; ++j ) {
						if ( shell0[j] ) {
							for ( Size k=1; k<=nres; ++k ) {
								if ( !shell0[k] && neigh_merge[k] ) shell1[k] = true;
							}
						}
					}
				} else {
					shell1 = neighbor[currres];
					for ( Size i=1; i<=NEXP_; ++i ) {
						shell0 = shell1;
						for ( Size j=1; j<=nres; ++j ) {
							if ( shell0[j] ) {
								for ( Size k=1; k<=nres; ++k ) {
									if ( !shell0[k] && neighbor[j][k] ) shell1[k] = true;
								}
							}
						}
					}
				}

				// build 1 residue packer task
				core::pack::task::PackerTaskOP ptask_working (core::pack::task::TaskFactory::create_packer_task( pose ));
				ptask_working->restrict_to_residues(shell1);
				ptask_working->or_include_current(true);

				for ( Size j=1, j_end = nres; j<= j_end; ++j ) {
					if ( shell0[j] ) {
						visited[j] = true;
						dynamic_cast<core::pack::task::ResidueLevelTask_&>
							(ptask_working->nonconst_residue_task(j)).update_commutative( ptask_resfile->nonconst_residue_task(j) );
					} else if ( shell1[j] ) {
						ptask_working->nonconst_residue_task(j).restrict_to_repacking();
						ptask_working->nonconst_residue_task(j).or_include_current(true);
					}
				}

				// set up movemap
				core::kinematics::MoveMapOP mm(new core::kinematics::MoveMap);
				mm->set_jump(true);
				for ( Size j=1, j_end = nres; j<= j_end; ++j ) {
					if ( shell1[j] ) {
						mm->set_bb(j, false); mm->set_chi(j, true);
					} else {
						mm->set_bb(j, false); mm->set_chi(j, false);
					}
				}

				for ( Size j=1, j_end = nres; j<= j_end; ++j ) {
					if ( shell0[j] ) {
						// allow a window of bb movement around each central residue
						mm->set_bb(j, true);
						mm->set_chi(j, true);
						if ( j<nres ) { mm->set_bb(j+1, true); mm->set_chi(j+1, true); }
						if ( j>1 )    { mm->set_bb(j-1, true); mm->set_chi(j-1, true); }
					}
				}

				if ( core::pose::symmetry::is_symmetric(pose) )  {
					core::pose::symmetry::make_symmetric_movemap( pose, *mm );
				}

				Size nvis=0,nsh0=0,nsh1=0;
				for ( Size j=1, j_end = nres; j<= j_end; ++j ) {
					if ( visited[j] ) nvis++;
					if ( shell0[j] ) nsh0++;
					if ( shell1[j] ) nsh1++;
				}

				// optimize
				optimization_loop( pose, ptask_working, mm,  ramp_schedule_[innercyc], 1e-4 );
				TR << "[" << cyc << "." << innercyc << "] res " << currres << " [" << nsh0 << "/" << nsh1 << "] ("
					<< nvis << "/" << nres_asu << ")  E=" << (*min_sfxn_)(pose)
					<< "  ramp=" << ramp_schedule_[innercyc] << std::endl;
			}
		}
	}
}

}
}

