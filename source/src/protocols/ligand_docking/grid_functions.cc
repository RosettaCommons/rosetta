// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/grid_functions.cc
///
/// @brief
/// @author Ian W. Davis


#include <protocols/ligand_docking/grid_functions.hh>

#include <core/graph/Graph.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/SingleResidueRotamerLibrary.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random_permutation.hh>
#include <utility/vector1.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <utility/vector0.hh>


namespace protocols {
namespace ligand_docking {

typedef utility::pointer::shared_ptr<core::grid::CartGrid<int> > CartGridIntOP;

static thread_local basic::Tracer TR( "protocols.ligand_docking.grid_functions", basic::t_debug );


/// @details If score exceeds max_score, stop counting and return (faster).
/// Atoms that fall outside the grid are scored as zero.
//template<typename T>
int grid_score(
	core::grid::CartGrid<int> const & grid,
	core::conformation::Residue const & rsd,
	int max_score /*= 9999*/
)
{
	int score = 0;
	for(core::Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end && score < max_score; ++i) {
		core::Vector const & atom = rsd.xyz(i);
		if( grid.is_in_grid( atom.x(), atom.y(), atom.z() ) ) {
			score += grid.getValue( atom.x(), atom.y(), atom.z() );
		}
	}
	TR << "Grid score: " << score << std::endl;
	return score;
}


/// @details If rep exceeds max_rep, atr is likely to be incorrect.
void grid_score_atr_rep(
	core::grid::CartGrid<int> const & grid,
	core::conformation::Residue const & rsd,
	int & atr_out,
	int & rep_out,
	int max_rep /*= 9999*/
)
{
	for(core::Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end && rep_out <= max_rep; ++i) {
		core::Vector const & atom = rsd.xyz(i);
		if( grid.is_in_grid( atom.x(), atom.y(), atom.z() ) ) {
			int val = grid.getValue( atom.x(), atom.y(), atom.z() );
			if( val > 0 ) rep_out += val;
			else atr_out += val;
		}
	}
	TR << "Grid score  rep=" << rep_out << "  atr=" << atr_out << std::endl;
}

/// @brief Sum the grid values for all heavy atoms in the residue
void rb_grid_score_atr_rep(
	core::grid::CartGrid<int> const & grid,
	core::pose::Pose const & pose,
	core::Size begin,
	core::Size const end,
	int & atr_out, //< sum of negative grid values
	int & rep_out, //< sum of positive grid values
	int max_rep
){
	for(; begin <= end; ++begin){
		core::conformation::Residue const & r= pose.residue(begin);
		grid_score_atr_rep(grid, r, atr_out, rep_out, max_rep);
	}
}

std::pair<int, int> get_rb_atr_and_rep_scores(
		core::grid::CartGrid<int> const & grid,
		core::pose::Pose const & pose,
		core::Size begin,
		core::Size end
){
	int atr=0;
	int rep=0;
	rb_grid_score_atr_rep(grid, pose, begin, end, atr, rep);
	return std::pair<int, int>(atr,rep);
}

/// @details min_score exists so you can search for the rotamer that clashes with the grid
/// a minimal ammount instead of none, to disfavor free-floating, non-interacting poses.
void grid_rotamer_trials(
	core::grid::CartGrid<int> const & grid,
	core::pose::Pose & pose,
	core::Size rsd_no,
	int const min_score /*= 0*/
)
{
	using core::conformation::ResidueOP;

	// Retrieve list of conformers
	utility::vector1< ResidueOP > conformers;
	rotamers_for_trials(pose, rsd_no, conformers);
	if( conformers.empty() ) return;
	// Conformers are already aligned onto original residue coords at this point.

	// Select best conformer, stopping early if score == 0
	int best_score = 9999;
	core::Size best_i = 0;
	for(core::Size i = 1, i_end = conformers.size(); i <= i_end; ++i) {
		ResidueOP rsd = conformers[i];
		int const new_score = grid_score(grid, *rsd, best_score);
		if( new_score < best_score && new_score >= min_score ) {
			best_score = new_score;
			best_i = i;
			if(best_score == min_score) break;
		}
	}

	// Replace current residue with best residue
	if( best_i > 0 ) {
		TR << "Best fit is conformer " << best_i << " with score " << best_score << std::endl;
		// Residue library has already superimpose residues appropriately, so don't orient again
		pose.replace_residue(rsd_no, *(conformers[best_i]), false /*orient backbone*/);
	}
}


void grid_rotamer_trials_atr_rep(
	core::grid::CartGrid<int> const & grid,
	core::pose::Pose & pose,
	core::Size rsd_no
)
{
	using core::conformation::ResidueOP;

	// Retrieve list of conformers
	utility::vector1< ResidueOP > conformers;
	rotamers_for_trials(pose, rsd_no, conformers);
	if( conformers.empty() ) return;
	TR<< "conformers.size: "<< conformers.size()<< std::endl;
	// Conformers are already aligned onto original residue coords at this point.

	// Select best conformer, stopping early if score == 0
	int const perfect_rep = 0, perfect_atr = -(pose.residue(rsd_no).nheavyatoms());
	int best_rep = 9999, best_atr = 0;
	core::Size best_i = 0;
	for(core::Size i = 1, i_end = conformers.size(); i <= i_end; ++i) {
		ResidueOP rsd = conformers[i];
		int new_atr(0), new_rep(0);
		grid_score_atr_rep(grid, *rsd, new_atr, new_rep, best_rep);
		if( new_rep < best_rep || (new_rep == best_rep && new_atr < best_atr) ) {
			best_rep = new_rep;
			best_atr = new_atr;
			best_i = i;
			if(best_rep == perfect_rep && best_atr == perfect_atr) break;
		}
	}

	// Replace current residue with best residue
	if( best_i > 0 ) {
		TR << "Best fit is conformer " << best_i << " with score rep=" << best_rep << " atr=" << best_atr << std::endl;
		// Residue library has already superimpose residues appropriately, so don't orient again
		pose.replace_residue(rsd_no, *(conformers[best_i]), false /*orient backbone*/);
	}
}

void rb_grid_rotamer_trials_atr_rep(
		core::grid::CartGrid<int> const & grid,
		core::pose::Pose & pose,
		core::Size begin,
		core::Size end
){
	// Residues don't interact with each other -- no reason to do this more than
	// once with the same pose. Likewise, ligand residues don't interact.
	for(; begin <= end; ++begin){
			TR<< "now performing rotamer trials on "<< begin << std::endl;
			grid_rotamer_trials_atr_rep(grid, pose, begin);
		}
}

void rotamers_for_trials(
	core::pose::Pose & pose,
	core::Size rsd_no,
	utility::vector1< core::conformation::ResidueOP > & conformers_out
)
{
	using core::conformation::ResidueOP;
	using namespace core::pack::task;

	// Dummy parameters that the ligand rotamer library doesn't use:
	core::scoring::ScoreFunction dummy_sfxn;
	PackerTaskOP dummy_task = TaskFactory::create_packer_task(pose); // actually used, in a trivial way

	core::graph::GraphCOP empty_graph( new core::graph::Graph() );
	// Retrieve conformers
	core::pack::dunbrack::SingleResidueRotamerLibraryCOP reslib = core::pack::dunbrack::RotamerLibrary::get_instance().get_rsd_library( pose.residue_type(rsd_no) ).lock();
	if( reslib.get() == NULL ) return;

	core::chemical::ResidueType const & res_type =  pose.residue_type(rsd_no);
	utility::vector1< utility::vector1< core::Real > > empty_extra_chi_steps( res_type.nchi(), utility::vector1< core::Real >() );

	// Retrieve list of conformers
	reslib->fill_rotamer_vector(
		pose,
		dummy_sfxn,
		*dummy_task,
		empty_graph,
		pose.residue_type(rsd_no).get_self_ptr(), //ResidueTypeCOP
		pose.residue(rsd_no),
		empty_extra_chi_steps,
		true /* sure, let's pretend it's buried */,
		conformers_out // output appended here
	);
	// Conformers are already aligned onto original residue coords at this point.

	// Because scoring is so "coarse", with many conformers possibly getting the same score,
	// it's important to try them in random order to avoid biases.
	numeric::random::random_permutation(conformers_out, numeric::random::rg());
}


/// @details Make a bounding box around the sphere, and visit all grid points
/// that the box intersects.  If the grid point center is within the sphere,
/// fill that grid space with the specified value.
void set_sphere(
	core::grid::CartGrid<int> & grid,
	core::Vector const & center,
	core::Real radius,
	int value
)
{
	typedef core::grid::CartGrid<int>::GridPt GridPt;
	using core::Vector;
	using namespace std; // min, max

	core::Real radius2 = radius*radius;
	int nx(0), ny(0), nz(0); // grid points in each dimension
	grid.getNumberOfPoints(nx, ny, nz);
	Vector vec_rad( radius );
	GridPt grid_min = grid.gridpt( center - vec_rad );
	GridPt grid_max = grid.gridpt( center + vec_rad );
	for(int i = max(0, grid_min.x()), i_end = min(nx-1, grid_max.x()); i <= i_end; ++i) {
		for(int j = max(0, grid_min.y()), j_end = min(ny-1, grid_max.y()); j <= j_end; ++j) {
			for(int k = max(0, grid_min.z()), k_end = min(nz-1, grid_max.z()); k <= k_end; ++k) {
				GridPt grid_pt(i, j, k);
				//std::cout << "Checking grid pt " << i << " " << j << " " << k << std::endl;
				Vector box_ctr = grid.coords( grid_pt );
				if( box_ctr.distance_squared(center) <= radius2 ) grid.setValue(grid_pt, value);
			}
		}
	}
}

void set_repulsive_bb_cores( utility::pointer::shared_ptr<core::grid::CartGrid<int> >grid, core::pose::Pose const & pose, core::Real const rep_rad){
	// Set repulsive core around each backbone heavy atom (including CB)
	for(core::Size r = 1, r_end = pose.total_residue(); r <= r_end; ++r) {
		core::conformation::Residue const & rsd = pose.residue(r);
		if( !rsd.is_protein() ) continue;
		if( rsd.has("CB") ) set_sphere(*grid, rsd.xyz("CB"), rep_rad, 1);
		if( rsd.has("N") ) set_sphere(*grid, rsd.xyz("N"), rep_rad, 1);
		if( rsd.has("CA") ) set_sphere(*grid, rsd.xyz("CA"), rep_rad, 1);
		if( rsd.has("C") ) set_sphere(*grid, rsd.xyz("C"), rep_rad, 1);
		if( rsd.has("O") ) set_sphere(*grid, rsd.xyz("O"), rep_rad, 1);
	}
}

/// @detail this function assumes there is only one ligand so it only considers protein residues
utility::pointer::shared_ptr<core::grid::CartGrid<int> > make_atr_rep_grid(
	core::pose::Pose const & pose,
	core::Vector const & center
)
{
	TR << "make_atr_rep_grid"<<std::endl;
	using namespace core;
	int const num_pts = 160;
	Real const resol = 0.25;
	Real const grid_halfwidth = (num_pts * resol) / 2.0;
	Real const atr_rad = 4.75; // optimal contact distance <= 4A for most atom pairs (but we're ignoring H)
	Real const rep_rad = 2.25; // don't want to exclude H-bonds (~2.8A heavy-heavy) or make clefts too narrow

	// Designer of this class did not believe in RAII:
	CartGridIntOP grid( new core::grid::CartGrid<int>() );
	grid->setBase(
		center.x() - grid_halfwidth,
		center.y() - grid_halfwidth,
		center.z() - grid_halfwidth
	);
	grid->setDimensions(num_pts, num_pts, num_pts, resol, resol, resol);
	grid->setupZones();
	grid->zero();

	// Set attractive zones around all heavy atoms -- assume most ligand atoms
	// will be near *something*, and most sidechains will stay put.
	for(Size r = 1, r_end = pose.total_residue(); r <= r_end; ++r) {
		conformation::Residue const & rsd = pose.residue(r);
		if( !rsd.is_protein() ) continue;
		for(Size a = 1, a_end = rsd.nheavyatoms(); a <= a_end; ++a)
		{
			set_sphere(*grid, rsd.xyz(a), atr_rad, -1);
		}
	}

	// Set neutral core around each sidechain heavy atom, as MOST of these stay put.
	for(Size r = 1, r_end = pose.total_residue(); r <= r_end; ++r) {
		conformation::Residue const & rsd = pose.residue(r);
		if( !rsd.is_protein() ) continue;
		for(Size a = rsd.first_sidechain_atom(), a_end = rsd.nheavyatoms(); a <= a_end; ++a)
		{
			set_sphere(*grid, rsd.xyz(a), rep_rad, 0);
		}
	}

	set_repulsive_bb_cores(grid, pose, rep_rad);

	return grid;
}

/// @detail this function assumes excludes one ligand from the grid
utility::pointer::shared_ptr<core::grid::CartGrid<int> > make_atr_rep_grid_without_ligand(
	core::pose::Pose const & pose,
	core::Vector const & center,
	core::Size const & ligand_chain_id_to_exclude
)
{
	utility::vector1<core::Size> ligand_chain_ids_to_exclude;
	ligand_chain_ids_to_exclude.push_back(ligand_chain_id_to_exclude);
	return make_atr_rep_grid_without_ligands(pose, center, ligand_chain_ids_to_exclude);
}

/// @brief Make a grid around the specified point with attractive (negative)
/// and repulsive (positive) values for all heavy atoms not in ligand_chain_id_to_exclude
utility::pointer::shared_ptr<core::grid::CartGrid<int> > make_atr_rep_grid_without_ligands(
	core::pose::Pose const & pose,
	core::Vector const & center,
	utility::vector1<core::Size> ligand_chain_ids_to_exclude
){
	TR << "make_atr_rep_grid_without_ligands"<<std::endl;
	using namespace core;
	int const num_pts = 160;
	Real const resol = 0.25;
	Real const grid_halfwidth = (num_pts * resol) / 2.0;
	Real const atr_rad = 4.75; // optimal contact distance <= 4A for most atom pairs (but we're ignoring H)
	Real const rep_rad = 2.25; // don't want to exclude H-bonds (~2.8A heavy-heavy) or make clefts too narrow

	// Designer of this class did not believe in RAII:
	CartGridIntOP grid( new grid::CartGrid<int>() );
	grid->setBase(
		center.x() - grid_halfwidth,
		center.y() - grid_halfwidth,
		center.z() - grid_halfwidth
	);
	grid->setDimensions(num_pts, num_pts, num_pts, resol, resol, resol);
	grid->setupZones();
	grid->zero();

	// Set attractive zones around all heavy atoms -- assume most ligand atoms
	// will be near *something*, and most sidechains will stay put.
	for(Size r = 1, r_end = pose.total_residue(); r <= r_end; ++r) {
		conformation::Residue const & rsd = pose.residue(r);
		if( find(
				ligand_chain_ids_to_exclude.begin(),
				ligand_chain_ids_to_exclude.end(),
				rsd.chain()
			) !=  ligand_chain_ids_to_exclude.end()
		) {
			continue;
		}
		for(Size a = 1, a_end = rsd.nheavyatoms(); a <= a_end; ++a)
		{
			set_sphere(*grid, rsd.xyz(a), atr_rad, -1);
		}
	}

	// Set neutral core around each sidechain heavy atom, as MOST of these stay put.
	for(Size r = 1, r_end = pose.total_residue(); r <= r_end; ++r) {
		conformation::Residue const & rsd = pose.residue(r);
		if( find(
				ligand_chain_ids_to_exclude.begin(),
				ligand_chain_ids_to_exclude.end(),
				rsd.chain()
			) !=  ligand_chain_ids_to_exclude.end()
		) {
			continue;
		}
		if( rsd.is_protein() ){
			for(Size a = rsd.first_sidechain_atom(), a_end = rsd.nheavyatoms(); a <= a_end; ++a)
			{
				set_sphere(*grid, rsd.xyz(a), rep_rad, 0);//ligand can run into side-chains, they'll be repacked
			}
		}else{
			for(Size a = 1, a_end = rsd.nheavyatoms(); a <= a_end; ++a)
			{
				set_sphere(*grid, rsd.xyz(a), rep_rad, 1); //ligand shouldn't run into other ligands
			}
		}
		// else don't add this ligand to the grid
	}

	set_repulsive_bb_cores(grid, pose,rep_rad);

	return grid;
}





} // namespace ligand_docking
} // namespace protocols
