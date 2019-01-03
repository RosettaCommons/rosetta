// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRDummySpinlabelEnsemble.hh
/// @brief   Implementation of classes NMRDummySpinlabelConformer and NMRDummySpinlabelEnsemble
/// @details last Modified: 28/11/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.hh>

// Package headers
#include <core/scoring/nmr/NMRDummySpinlabelVoxelGrid.hh>
#include <core/io/nmr/util.hh>
#include <core/scoring/nmr/util.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/scoring/rms_util.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/io/util.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/constants.hh>
#include <numeric/random/random.hh>
#include <numeric/HomogeneousTransform.hh>

// C++ headers
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <fstream>
#include <iterator>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>

// Boost headers
#include <boost/algorithm/string.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <boost/unordered/unordered_set.hpp>
#include <boost/functional/hash.hpp>

namespace core {
namespace scoring {
namespace nmr {

static basic::Tracer TR( "core.scoring.nmr.NMRDummySpinlabelEnsemble" );

//////////////////////////////////////////////////
/// Class NMRDummySpinlabelConformer
//////////////////////////////////////////////////

/// @brief Construct from ID, number of observations, frequency and residue
NMRDummySpinlabelConformer::NMRDummySpinlabelConformer(
	Size const id,
	Size const nobs,
	Real const freq,
	conformation::Residue const & residue
) :
	utility::pointer::ReferenceCount(),
	id_(id),
	nobs_(nobs),
	frequency_(freq),
	clash_(false),
	clash_score_(0.0),
	residue_( new core::conformation::Residue( residue ) )
{
	for ( Size n(1); n <= residue.natoms(); ++n ) {
		atom_table_.insert(std::make_pair(residue.atom_name(n), NMRDummySpinlabelAtom(residue.atom_name(n), id, residue.atom(n).xyz(), *this)));
	}
}

/// @brief Construct from ID and number of observations (frequency)
///        as well as atom names and xyz coordinates
NMRDummySpinlabelConformer::NMRDummySpinlabelConformer(
	Size const id,
	Size const nobs,
	Real const freq,
	chemical::ResidueType const & restype,
	utility::vector1< std::string > const & names,
	utility::vector1< Vector > const & coords
) :
	utility::pointer::ReferenceCount(),
	id_(id),
	nobs_(nobs),
	frequency_(freq),
	clash_(false),
	clash_score_(0.0)
{
	runtime_assert_msg(names.size() == coords.size(),
		"ERROR during instantiation of NMRDummySpinlabelConformer. Vector of atom names and xyz coordinates must have the same size.");

	// Build residue from atom names and coordinates
	// Fill table of NMRSpinlabelAtoms
	residue_ = conformation::ResidueFactory::create_residue( restype );
	std::set< std::string > skipped_atom_names;
	for ( Size i(1), i_end(names.size()); i <= i_end; ++i ) {
		atom_table_.insert(std::make_pair(names[i], NMRDummySpinlabelAtom(names[i], id, coords[i], *this)));
		if ( residue_->has( names[i] ) ) {
			residue_->set_xyz( names[i], coords[i] );
		} else if ( skipped_atom_names.count(names[i]) == 0 ) {
			TR.Warning << "Skipping unrecognized atom '" << names[i] << "' in building NMRDummySpinlabelConformer of type " << restype.name() << std::endl;
			skipped_atom_names.insert(names[i]);
		}
	}
	// TODO function for filling missing atoms

	// Torsion angles are not automatically calculated from the coordinates.
	// We should manually assign chi() and mainchain_torsions() for each conformer.
	conformation::set_chi_according_to_coordinates( *residue_ );
}

NMRDummySpinlabelConformer::NMRDummySpinlabelConformer(
	Size const id,
	Size const nobs,
	Real const freq,
	chemical::ResidueType const & restype,
	AtomNamePosMap const & names_coords
) :
	utility::pointer::ReferenceCount(),
	id_(id),
	nobs_(nobs),
	frequency_(freq),
	clash_(false),
	clash_score_(0.0)
{
	// Build residue from atom names and coordinates
	// Fill table of NMRSpinlabelAtoms
	residue_ = conformation::ResidueFactory::create_residue( restype );
	std::set< std::string > skipped_atom_names;
	for ( auto const & atm : names_coords ) { // Pair of atom names and xyz vector
		atom_table_.insert(std::make_pair(atm.first, NMRDummySpinlabelAtom(atm.first, id, atm.second, *this)));
		if ( residue_->has( atm.first ) ) {
			residue_->set_xyz( atm.first, atm.second );
		} else if ( skipped_atom_names.count(atm.first) == 0 ) {
			TR.Warning << "Skipping unrecognized atom '" << atm.first << "' in building NMRDummySpinlabelConformer of type " << restype.name() << std::endl;
			skipped_atom_names.insert(atm.first);
		}
	}
	// TODO function for filling missing atoms

	// Torsion angles are not automatically calculated from the coordinates.
	// We should manually assign chi() and mainchain_torsions() for each conformer.
	conformation::set_chi_according_to_coordinates( *residue_ );
}

/// @brief copy constructor
NMRDummySpinlabelConformer::NMRDummySpinlabelConformer(NMRDummySpinlabelConformer const & other) :
	utility::pointer::ReferenceCount( other ),
	id_(other.id_),
	nobs_(other.nobs_),
	frequency_(other.frequency_),
	clash_(other.clash_),
	clash_score_(other.clash_score_),
	residue_( new conformation::Residue( *(other.residue_) ) )
{
	atom_table_.clear();
	atom_table_.reserve(other.atom_table_.size());
	for ( auto const & atm : other.atom_table_ ) {
		atom_table_.insert(std::make_pair(atm.first, NMRDummySpinlabelAtom(atm.second.atom_name(), id_, atm.second.get_coordinates(), *this)));
	}
}

/// @brief assignment operator
NMRDummySpinlabelConformer &
NMRDummySpinlabelConformer::operator=(NMRDummySpinlabelConformer const & rhs) {
	if ( this != &rhs ) {
		id_ = rhs.id_;
		nobs_= rhs.nobs_;
		frequency_ = rhs.frequency_;
		clash_ = rhs.clash_;
		clash_score_ = rhs.clash_score_;
		residue_ = ResidueOP( new conformation::Residue( *(rhs.residue_) ) );
		atom_table_.clear();
		atom_table_.reserve(rhs.atom_table_.size());
		for ( auto const & atm : rhs.atom_table_ ) {
			atom_table_.insert(std::make_pair(atm.first, NMRDummySpinlabelAtom(atm.second.atom_name(), id_, atm.second.get_coordinates(), *this)));
		}
	}
	return *this;
}

/// @brief destructor
NMRDummySpinlabelConformer::~NMRDummySpinlabelConformer() { }

//////////////////////////////////////////////////
/// Class NMRDummySpinlabelEnsemble
//////////////////////////////////////////////////

/// @brief Construct from database file
NMRDummySpinlabelEnsemble::NMRDummySpinlabelEnsemble(std::string const & database_file, chemical::ResidueType const & restype) :  /* throw(utility::excn::FileNotFound, utility::excn::Exception) */
	utility::pointer::ReferenceCount(),
	grid_(nullptr),
	elaborate_clash_check_(false)
{
	init_from_database_file(database_file, restype);

	// Define the local origin of the coordinate frame of the ensemble around which we
	// will build a voxel grid used to quickly identify neighbors for the clash score calculation.
	define_ensemble_frame();
	register_options();
	init_from_cml();
}

/// @brief destructor
NMRDummySpinlabelEnsemble::~NMRDummySpinlabelEnsemble() { }

/// @brief utility function used for initialization from database files
void
NMRDummySpinlabelEnsemble::init_from_database_file(std::string const & database_file, chemical::ResidueType const & restype)
{
	// Open database file
	utility::vector1< std::string > lines_database_file(utility::io::get_lines_from_file_data(database_file));
	TR.Debug << "Creating NMRDummySpinlabelEnsemble from database file " << database_file << std::endl;

	// Read lines from spin-label conformer database file (in PDB format)
	// and store atom names and xyz coordinates temporarily in AtomNamePosMap
	AtomNamePosMap name_pos_map;
	Size sl_id(1);
	std::map< Size, AtomNamePosMap > id_to_conformer_atom_map;
	Size n_lines_database_file(lines_database_file.size());
	for ( Size i = 1; i <= n_lines_database_file; ++i ) {
		if ( lines_database_file[i].find("HETATM") != std::string::npos ||
				lines_database_file[i].find("ATOM")   != std::string::npos ) {

			std::string atom_name = utility::strip(lines_database_file[i].substr(12,4)); // Strip whitespaces from the atom name
			Vector atom_pos;
			atom_pos.x() = std::atof( lines_database_file[i].substr(30,8).c_str() );
			atom_pos.y() = std::atof( lines_database_file[i].substr(38,8).c_str() );
			atom_pos.z() = std::atof( lines_database_file[i].substr(46,8).c_str() );

			if ( name_pos_map.count( atom_name ) != 0 ) {
				TR.Warning << "Found multiple lines for atom " << atom_name << " for spinlabel conformer with ID " << sl_id << ". Use the latest atom entry." << std::endl;
			}
			name_pos_map[ atom_name ] = atom_pos;
		} else if ( lines_database_file[i].find("TER") != std::string::npos ) {
			if ( name_pos_map.size() > 0 ) {
				if ( id_to_conformer_atom_map.count( sl_id ) != 0 ) {
					TR.Warning << "Found multiple entries for spinlabel with ID " << sl_id << " in spinlabel conformers file. Use the latest entry." << std::endl;
				}
				id_to_conformer_atom_map[ sl_id ] = name_pos_map;
				name_pos_map.clear();
			} else {
				TR.Warning << "No HETATM entries found for spinlabel conformer with ID " << sl_id << ". Skip spinlabel " << sl_id << "." << std::endl;
			}
			++sl_id;
		}
	}
	ensemble_size_ = id_to_conformer_atom_map.size();
	if ( ensemble_size_ == 0 ) {
		utility_exit_with_message("ERROR during creation of NMRDummySpinlabelEnsemble. Number of conformers is zero.");
	}
	conformer_table_.reserve(ensemble_size_);

	// Create NMRDummySpinlabelConformers from information in AtomNamePosMap
	// and add them to vector in NMRDummySpinlabelEnsemble
	std::map< Size, AtomNamePosMap >::iterator conformer_map_iter, conformer_map_end;
	for ( conformer_map_iter = id_to_conformer_atom_map.begin(), conformer_map_end = id_to_conformer_atom_map.end();
			conformer_map_iter != conformer_map_end; ++conformer_map_iter ) {

		// Align every conformer with its backbone coordinates to origin
		HT conformer_frame = define_local_frame((*conformer_map_iter).second.at("CA"),
			(*conformer_map_iter).second.at("CB"), (*conformer_map_iter).second.at("C"));

		for ( AtomNamePosMap::iterator atom_map_iter( (*conformer_map_iter).second.begin() ),
				atom_map_end( (*conformer_map_iter).second.end() ); atom_map_iter != atom_map_end; ++atom_map_iter ) {
			(*atom_map_iter).second =  conformer_frame.inverse() * (*atom_map_iter).second;
		}
		TR.Debug << "Adding conformer to NMRDummySpinlabelEnsemble ..." << std::endl;
		conformer_table_.push_back( NMRDummySpinlabelConformerOP(
			new NMRDummySpinlabelConformer((*conformer_map_iter).first, 1, 1.0/ensemble_size_, restype, (*conformer_map_iter).second ) ) );
	}
	TR.Debug << "Created NMRDummySpinlabelEnsemble with " << conformer_table_.size() << " conformers." << std::endl;

	// Calculate RMSD matrix (this is after the conformers have been aligned to the backbone)
	TR.Debug << "Calculating RMSDs between conformers of NMRDummySpinlabelEnsemble" << std::endl;
	rmsd_mat_.resize(conformer_table_.size());
	for ( Size i=1; i <= rmsd_mat_.size(); ++i ) {
		rmsd_mat_[i].resize(rmsd_mat_.size());
		std::fill(rmsd_mat_[i].begin(),rmsd_mat_[i].end(),0.0);
	}
	for ( Size i=1; i < conformer_table_.size(); ++i ) {
		for ( Size j=i+1; j <= conformer_table_.size(); ++j ) {
			ResidueOP rsd1 = conformer_table_[i]->get_residue();
			ResidueOP rsd2 = conformer_table_[j]->get_residue();
			Real conf_rmsd = automorphic_rmsd(*rsd1,*rsd2,false);
			rmsd_mat_[i][j] = conf_rmsd;
			rmsd_mat_[j][i] = conf_rmsd;
			TR.Debug << "RMSD [ Conformer " << i << " ][ Conformer " << j << " ] = " << conf_rmsd << std::endl;
		}
	}
}

void
NMRDummySpinlabelEnsemble::register_options() {
	using namespace basic::options;
	option.add_relevant(OptionKeys::nmr::spinlabel::elaborate_rotamer_clash_check);
}
void
NMRDummySpinlabelEnsemble::init_from_cml() {
	using namespace basic::options;
	elaborate_clash_check_ = option[ basic::options::OptionKeys::nmr::spinlabel::elaborate_rotamer_clash_check ]();
}

/// @brief calculate the homogeneous transform that relates the target residue
///        (i.e. in our case the residue that represents the spinlabel site)
///        with the coordinate frame of the dummy spinlabel ensemble. In this
///        way we don't need to move the spinlabel ensemble which contains many
///        more coordinates but can leave it static.
HT
NMRDummySpinlabelEnsemble::coordinate_transform_from_target_site(
	pose::Pose const & pose,
	Size const resid
) const
{
	std::map< std::string, Vector > bb_coo;
	bb_coo.insert(std::make_pair("CA", Vector(0.0,0.0,0.0)));
	bb_coo.insert(std::make_pair("CB", Vector(0.0,0.0,0.0)));
	bb_coo.insert(std::make_pair("C",  Vector(0.0,0.0,0.0)));
	bb_coo.insert(std::make_pair("N",  Vector(0.0,0.0,0.0)));

	for ( auto const & atom : bb_coo ) {
		std::string atom_name = atom.first;
		if ( pose.residue(resid).has(atom_name) ) {
			bb_coo[atom_name] = pose.residue(resid).xyz(atom_name);
		} else if ( atom_name == "CB" && pose.residue(resid).has("2HA") ) { // special case for Glycine
			bb_coo[atom_name] = pose.residue(resid).xyz("2HA");
		} else if ( atom_name == "CB" && !pose.residue( resid ).has("CB") &&
				!pose.residue( resid ).has("2HA") ) {
			utility_exit_with_message("ERROR in coordinate transformation of NMRDummySpinlabelEnsemble. Target protein residue has no CB or 2HA atom.");
		} else {
			utility_exit_with_message("ERROR in coordinate transformation of NMRDummySpinlabelEnsemble. Target protein residue has no " + atom_name + " atom.");
		}
	}
	// local frame the spinlabel would have in the protein
	HT protein_frame = define_local_frame(bb_coo["CA"], bb_coo["CB"], bb_coo["C"]);
	HT transformation = ensemble_origin_ * protein_frame.inverse();
	return transformation;
}

/// @brief returns the transformation matrix which transforms a xyz coordinate from
///        the frame of the spinlabel ensemble into the frame of the target site residue.
HT
NMRDummySpinlabelEnsemble::coordinate_transform_onto_target_site(
	pose::Pose const & pose,
	Size const resid
) const
{
	// Since transformation proceeds into the opposite direction (SL frame -> target residue frame)
	// we simply return the inverse of the homogeneous transform this time.
	return coordinate_transform_from_target_site(pose, resid).inverse();
}

/// @brief define the origin of the xyz coordinate frame of the ensemble.
///        this needs to be done only once when we initialize the object.
void
NMRDummySpinlabelEnsemble::define_ensemble_frame() {
	if ( !(ensemble_size_ > 0) ) {
		utility_exit_with_message("ERROR in defining the local coordinate frame of the NMRDummySpinlabelEnsemble. Ensemble size must be > 0.");
	}

	NMRDummySpinlabelAtomTable const & conf_1_atom_table_ref = (*(conformer_table_.begin()))->get_atom_table();

	std::map< std::string, Vector > focus_bb_coo;
	focus_bb_coo.insert(std::make_pair("CA", Vector(0.0,0.0,0.0)));
	focus_bb_coo.insert(std::make_pair("CB", Vector(0.0,0.0,0.0)));
	focus_bb_coo.insert(std::make_pair("C",  Vector(0.0,0.0,0.0)));
	focus_bb_coo.insert(std::make_pair("N",  Vector(0.0,0.0,0.0)));

	for ( auto const & atom : focus_bb_coo ) {
		std::string atom_name = atom.first;
		if ( conf_1_atom_table_ref.find(atom_name) != conf_1_atom_table_ref.end() ) {
			focus_bb_coo[atom_name] = conf_1_atom_table_ref.at(atom_name).get_coordinates();
		} else {
			utility_exit_with_message("ERROR in defining the local coordinate frame of the NMRDummySpinlabelEnsemble. Spinlabel residue has no " + atom_name + " atom.");
		}
	}
	// Calculate the homogenous transform that defines the static local frame of the ensemble
	ensemble_origin_ = define_local_frame(focus_bb_coo["CA"], focus_bb_coo["CB"], focus_bb_coo["C"]);
}

/// @brief Perform a clash score calculation for every spinlabel conformer in the ensemble by
///        calculating pairwise distances to neighborhood residues which are within a given radius.
///        Specifically, the distance between every side-chain heavy atom in the
///        spinlabel conformer and the neighbor atom of a pose residue is calculated.
///        If the distance is smaller than the pose residue's neighbor radius, is is assumed
///        that the spinlabel conformer will fall within the side-chain radius of the
///        neighborhood residue and thus it will be labeled as clashing otherwise it will not.
void
NMRDummySpinlabelEnsemble::clash_check(
	pose::Pose const & pose,
	Size const target_resid,
	Real const radius
)
{
	using namespace select::residue_selector;

	// Set up neighborhood residue selector to select the spinlabel environment.
	// The third argument in the constructor means that the focus is not included in the selection.
	utility::vector1< bool > focus(pose.total_residue());
	focus[target_resid] = true;
	NeighborhoodResidueSelectorOP neighborhood_selector = NeighborhoodResidueSelectorOP( new NeighborhoodResidueSelector(focus, radius, false) );
	utility::vector1< bool > neighbor_residue_selection = neighborhood_selector->apply( pose );

	// Reset clash score boolean for all spinlabel conformers
	for ( Size i(1); i <= conformer_table_.size(); ++i ) {
		conformer_table_[i]->clash_off();
		conformer_table_[i]->clash_score() = 0.0;
	}

	TR.Debug << "Performing clash filter of NMRDummySpinlabelEnsemble by neighbor count method. Searching for neighbors within radius "
		<< radius << " around spinlabel site " << target_resid << "." << std::endl;
	// return if there are no neighbors
	auto nbr_resi_sum_calc = [](Size const previous, bool current) { return previous+Size(current); };
	if ( std::accumulate(neighbor_residue_selection.begin(), neighbor_residue_selection.end(), 0, nbr_resi_sum_calc) == 0 ) {
		TR.Debug << "No neighbors detected within " << radius << " angstrom. No NMRDummySpinlabel conformers were labeled as clashing." << std::endl;
		return;
	}

	// create voxel grid
	// resolution of the grid is the maximal NBR_RADIUS of the neighborhood residues
	Real resolution(0.0);
	for ( Size i(1), i_end(neighbor_residue_selection.size()); i <= i_end; ++i ) {
		if ( neighbor_residue_selection[i] && !pose.residue(i).is_virtual_residue() ) {
			if ( pose.residue(i).nbr_radius() > resolution ) {
				resolution = pose.residue(i).nbr_radius();
			}
		}
	}

	// If grid does not exist or if resolution is higher than previous value, create a new grid.
	// Note that we don't create a new grid if the resolution is lower to avoid reallocation
	// at every scoring step.
	if ( !grid_ || grid_->GetResolution() < resolution ) {
		TR.Debug << "Reallocating NMRSpinlabelVoxelGrid." << std::endl;
		utility::vector1< VoxelGridPoint const * > points;
		for ( NMRDummySpinlabelConformerTableIter conformer_table_iter( conformer_table_.begin() ), conformer_table_end( conformer_table_.end() );
				conformer_table_iter != conformer_table_end; ++conformer_table_iter ) {
			for ( NMRDummySpinlabelAtomTableCOIter sl_atom_table_iter( (*conformer_table_iter)->get_atom_table().begin() ),
					sl_atom_table_end( (*conformer_table_iter)->get_atom_table().end() ); sl_atom_table_iter != sl_atom_table_end; ++sl_atom_table_iter ) {
				// Use only SL side chain atoms for clash filter
				// Use only heavy atoms
				Size atmid = (*conformer_table_iter)->get_residue()->atom_index(sl_atom_table_iter->first);
				bool ishydrogen = (*conformer_table_iter)->get_residue()->atom_is_hydrogen(atmid);
				if ( sl_atom_table_iter->first != "N" && sl_atom_table_iter->first != "C" &&
						sl_atom_table_iter->first != "O" && sl_atom_table_iter->first != "CA" && !ishydrogen ) {
					points.push_back( &(sl_atom_table_iter->second) );
				}
			}
		}
		grid_ = NMRDummySpinlabelVoxelGridOP( new NMRDummySpinlabelVoxelGrid(resolution, points) );
	}

	// Coordinate transformation for amino acid neighbor atoms into the frame of the grid
	HT ht_to_ensemble_origin = coordinate_transform_from_target_site(pose, target_resid);

	std::string sl_atom, aa_atom;
	Size num_clash_conformers(0);

	// Iterate over neighbor residues
	for ( Size i(1), i_end(pose.total_residue()); i <= i_end; ++i ) {

		// early break if all conformers have been marked as clashing and we don't keep track
		// of how many clashes each individual conformer has
		if ( num_clash_conformers == ensemble_size_ && !elaborate_clash_check_ ) { break; }

		if ( neighbor_residue_selection[i] && !pose.residue(i).is_virtual_residue() ) { // exclude virtual residues

			TR.Debug << "Searching for NMRDummySpinlabel conformers that are within the neighbor radius of residue " << pose.residue(i).name3() << i << std::endl;
			Real neighborhood(pose.residue(i).nbr_radius());
			// For direct sequence neighbors of the spin-label, the max observed O-CB distance in ubiquitin (3.4473 A)
			// is used as cutoff distance. This is to avoid that this simple distance check is too over-sensitive for
			// AAs with large NBR_RADIUS (especially ARG) which are direct sequence neighbors of the spin-label.
			neighborhood = ( ( std::max(target_resid,i) - std::min(target_resid,i) ) <= 3 ) ? 3.4473 : neighborhood;
			Real tolerance(1.0e-9);
			if ( neighborhood > tolerance ) {
				// The residue's nbr atom gets the test atom to detect any clashes,
				// i.e. pairwise distances between the nbr atom and every spinlabel
				// heavy atom are calculated
				aa_atom = pose.residue(i).atom_name(pose.residue(i).nbr_atom());
				Vector aa_coords_in_grid_frame = ht_to_ensemble_origin * pose.residue(i).nbr_atom_xyz();

				// Make a grid point object and search for its neighbors within radius <= neighborhood
				VoxelGridPoint_AA aa_point(aa_atom, i, aa_coords_in_grid_frame, pose.residue(i));
				bool check_only_when_relevant(!elaborate_clash_check_);  // skip spinlabel atoms in conformers that have been labeled as clashing before
				utility::vector1< std::pair< VoxelGridPoint const *, Real > > neighbor_atoms = grid_->GetNeighbors(aa_point, neighborhood, check_only_when_relevant);

				// all spinlabel atoms in the neighbor list are within NBR_RADIUS of the amino acid residue
				// mark these spinlabels as clashing and add the number of atoms within NBR_RADIUS to the clash score
				for ( auto const & neighbor : neighbor_atoms ) {
					Size conformer_key(neighbor.first->id());
					TR.Debug << "Clash detected for NMRDummySpinlabel conformer " << conformer_key << std::endl;
					NMRDummySpinlabelConformer & conformer_ref = *(conformer_table_[ conformer_key ]);
					if ( !conformer_ref.has_clash() ) {
						TR.Debug << "Mark NMRDummySpinlabel conformer " << conformer_key << " as clashing." << std::endl;
						conformer_ref.clash_on();
						num_clash_conformers += 1;
					}
					// Keep track of the number of clashes in case we want to determine the 'best' conformer
					if ( elaborate_clash_check_ ) conformer_ref.clash_score() += 1.0;
				}
			}
		}
	} // neighborhood residues

	// If all spinlabel conformers clash with the protein, turn the conformer with the
	// lowest score (or a random conformer) on, so that we are still able to calculate the PRE/PCS
	if ( num_clash_conformers == ensemble_size_ ) {
		TR.Debug << "No conformer in NMRDummySpinlabelEnsemble passed neighbor count clash filter. Return emergency conformer instead." << std::endl;
		if ( elaborate_clash_check_ ) {
			auto compare_func = []( NMRDummySpinlabelConformerOP a, NMRDummySpinlabelConformerOP b) { return a->get_clash_score() < b->get_clash_score(); };
			(*std::min_element(conformer_table_.begin(), conformer_table_.end(), compare_func))->clash_off();
		} else {
			Size dummy = numeric::random::random_range(1, ensemble_size_);
			conformer_table_[ dummy ]->clash_off();
		}
	}

	if ( TR.Debug.visible() ) {
		auto sl_sum_calc = []( Size const previous, NMRDummySpinlabelConformerOP current)
			{ return previous + Size(current->has_clash()); };
		Size sum = std::accumulate(conformer_table_.begin(), conformer_table_.end(), 0, sl_sum_calc);
		TR.Debug << "Size of NMRDummySpinlabelEnsemble before neighbor count clash filter = " << ensemble_size_
			<< ". Size of NMRDummySpinlabelEnsemble after neighbor count clash filter = " << (ensemble_size_ - sum) << "." << std::endl;
		utility::vector1< NMRDummySpinlabelConformerOP >::reverse_iterator lastNonclashing =
			std::find_if(conformer_table_.rbegin(), conformer_table_.rend(), [](NMRDummySpinlabelConformerOP current) {return !current->has_clash();});
		TR.Debug << "Non-clashing NMRDummySpinlabel conformers: [ ";
		if ( lastNonclashing != conformer_table_.rend() ) {
			for ( Size i(1), i_end(std::distance(lastNonclashing,conformer_table_.rend())); i<i_end; ++i ) {
				if ( !conformer_table_[i]->has_clash() ) TR.Debug << conformer_table_[i]->get_id() << ", ";
			}
			TR.Debug << (*lastNonclashing)->get_id();
		}
		TR.Debug << " ]" << std::endl;
	}
}

/// @brief copy constructor
NMRDummySpinlabelEnsemble::NMRDummySpinlabelEnsemble(NMRDummySpinlabelEnsemble const & other) :
	utility::pointer::ReferenceCount( other ),
	ensemble_size_(other.ensemble_size_),
	ensemble_origin_(other.ensemble_origin_),
	rmsd_mat_(other.rmsd_mat_),
	elaborate_clash_check_(other.elaborate_clash_check_)
{
	conformer_table_.resize(other.conformer_table_.size());
	for ( Size i(1); i <= other.conformer_table_.size(); ++i ) {
		conformer_table_[i] =  NMRDummySpinlabelConformerOP( new NMRDummySpinlabelConformer( *(other.conformer_table_[i]) ) );
	}
	if ( grid_ ) {
		grid_->Clear();
		grid_.reset();
	}
}

/// @brief assignment operator
NMRDummySpinlabelEnsemble &
NMRDummySpinlabelEnsemble::operator=(NMRDummySpinlabelEnsemble const & rhs) {
	if ( this != &rhs ) {
		ensemble_size_ = rhs.ensemble_size_;
		conformer_table_.resize(rhs.conformer_table_.size());
		for ( Size i(1); i <= rhs.conformer_table_.size(); ++i ) {
			conformer_table_[ i ] = NMRDummySpinlabelConformerOP ( new NMRDummySpinlabelConformer( *(rhs.conformer_table_[i]) ) );
		}
		ensemble_origin_ = rhs.ensemble_origin_;
		rmsd_mat_ = rhs.rmsd_mat_;
		elaborate_clash_check_ = rhs.elaborate_clash_check_;
		if ( grid_ ) {
			grid_->Clear();
			grid_.reset();
		}
	}
	return *this;
}

} // nmr
} // scoring
} // core
