#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/annotated_sequence.hh>  //make_pose_from_sequence
#include <core/pose/util.hh>

#include <core/import_pose/import_pose.hh>

#include <core/conformation/Residue.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/util/SwitchResidueTypeSet.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

// minimize pose into density
#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>


#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/fragment/util.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>

//#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>

#include <iostream>
#include <string>
#include <list>
#include <algorithm>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <apps/pilot/rayyrw/min_pack_min.hh>
#include <apps/pilot/rayyrw/util.hh>


static basic::Tracer TR("place_fragment_into_density");
static basic::Tracer func("density_grid_search_func");

// options declaration
class ThisApplication  {
public:
	ThisApplication();
	static void register_options();
private:
};

ThisApplication::ThisApplication()
{}

OPT_KEY( File, fragfile )
OPT_KEY( Integer, num_frags )
OPT_KEY( IntegerVector, pos )
OPT_KEY( IntegerVector, designated_rank )
OPT_KEY( Integer, keep )
OPT_KEY( Integer, movestep )
OPT_KEY( Integer, skip_edge )
OPT_KEY( Boolean, min_pack_min )
OPT_KEY( Integer, cycles )
OPT_KEY( Boolean, min_bb )
OPT_KEY( Boolean, median_radius )
OPT_KEY( Real, no_density_score )
OPT_KEY( Boolean, debug )
OPT_KEY( File, debug_outfn )

using namespace ObjexxFCL::format;

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( fragfile, "fragment file name", "" );
	NEW_OPT( num_frags, "number of top fragments to use for each residue position", 25 );
	NEW_OPT( pos, "only calculate the designate position", 0 );
	NEW_OPT( designated_rank, "only calculate the designate rank", 0 );
	NEW_OPT( keep, "how many you want keep for each search", 10 );
	NEW_OPT( movestep,  "how big a move for grid search", 2 );
	NEW_OPT( skip_edge, "skip the edge region for search", 0 );
	NEW_OPT( min_pack_min, "MinPackMin", false );
	NEW_OPT( cycles, "how many cycles of min_pack_min", 1 );
	NEW_OPT( min_bb, "backbone minimization - add phipsi minimization", false );
	NEW_OPT( median_radius, "get median radius from fragments set given a position", false );
	NEW_OPT( no_density_score, "no_density_score to prevent fitting into no density region s-den_min/den_stdv", 3.0 );
	//NEW_OPT( no_density_score, "no_density_score to prevent fitting into no density region", 0.0 );
	NEW_OPT( debug, "print out debug trans/rotat matrix", false );
	NEW_OPT( debug_outfn, "debug output file", "debug_outfile.txt" );
}

class Results_Keeper{
public:
	core::Real score_;
	numeric::xyzVector<core::Real> pre_trans_;
	numeric::xyzVector<core::Real> translation_;
	numeric::xyzMatrix<core::Real> rotation_;
	numeric::xyzVector<core::Real> mycenterofmass_;
	double map_s_, map_s2_, pose_s_, pose_s2_;

	Results_Keeper();
	Results_Keeper( core::Real, numeric::xyzVector<core::Real>, numeric::xyzVector<core::Real>, numeric::xyzMatrix<core::Real>, numeric::xyzVector<core::Real>, double, double, double, double );

private:
};

Results_Keeper::Results_Keeper(){} // default constructor

Results_Keeper::Results_Keeper(
	core::Real score,
	numeric::xyzVector<core::Real> pre_trans,
	numeric::xyzVector<core::Real> translation,
	numeric::xyzMatrix<core::Real> rotation,
	numeric::xyzVector<core::Real> mycenterofmass,
	double map_s,
	double map_s2,
	double pose_s,
	double pose_s2
){
	score_          = score;
	pre_trans_      = pre_trans;
	translation_    = translation;
	rotation_       = rotation;
	mycenterofmass_ = mycenterofmass;
	map_s_          = map_s;
	map_s2_         = map_s2;
	pose_s_         = pose_s;
	pose_s2_        = pose_s2;
	//std::cout << "non-default constructor gets called: " << std::endl;
}

/////////////////////////////    end of results_keeper class definition   ////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
// get_radius: calculate the radius given a protein
//    note: set default value ATOM_MASK to be 3.2, according to the basic/options/options_rosetta.py
//    nRsteps: # of shells
core::Size
get_radius(
	core::pose::Pose & pose
){
	//core::Real ATOM_MASK = 3.2;
	core::Real ATOM_MASK = basic::options::option[ basic::options::OptionKeys::edensity::atom_mask ](); // the default value is 3.2
	core::Size nRsteps;
	core::Real delRsteps=2.0;

	utility::vector1< numeric::xyzVector< core::Real > > atmList;
	numeric::xyzVector< core::Real > massSum(0.0,0.0,0.0);
	for ( core::Size i=1; i<= pose.size(); ++i ) {
		core::conformation::Residue const & rsd( pose.residue(i) );
		if ( rsd.aa() == core::chemical::aa_vrt ) continue;
		// use every atom ...
		for ( core::Size j=1; j<= rsd.nheavyatoms(); ++j ) {
			core::conformation::Atom const & atom( rsd.atom(j) );
			atmList.push_back(atom.xyz());
			massSum += atom.xyz();
		}
	}
	int nAtms=atmList.size();
	massSum /= nAtms;
	core::Real maxSize = 0;
	for ( core::Size i=1; i<=atmList.size(); ++i ) {
		maxSize = std::max( maxSize, massSum.distance_squared(atmList[i]) );
	}
	maxSize = std::sqrt( maxSize )+std::min(ATOM_MASK, 5.0); // this the radius
	nRsteps = (int)std::floor(maxSize/delRsteps + 0.5);

	return nRsteps;
}


// calculate median given a vector of integers
core::Size
get_median(
	utility::vector1< core::Size > & v
){
	std::size_t n = v.size() / 2;
	std::nth_element(v.begin(), v.begin()+n, v.end());
	return v[n];
}


// a function that loop over each fragments, finding the median radius for spherical harmonic to calculate map
// the number (how many fragments are going to be used in a position) of input fragments has beeb controlled upstream
core::Size
get_radius_in_frag_set(
	core::fragment::FragSetCOP fragments,
	std::string & sequence,
	core::Size given_seq_pos
){
	utility::vector1< core::Size > container;
	//core::Size the_smallest_radius = 500;
	for ( core::fragment::ConstFrameIterator frame=fragments->begin(), eframe=fragments->end();
			frame != eframe;
			++frame ) {

		core::Size seq_pos = frame->start();
		if ( seq_pos == given_seq_pos ) {

			core::Size Nmers_size = frame->length();
			core::pose::Pose frag_pose;
			std::string frag_seq = sequence.substr( seq_pos-1, Nmers_size );
			core::pose::make_pose_from_sequence( frag_pose, frag_seq, *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )) );

			for ( core::Size j=1; j<=frame->nr_frags(); ++j ) {
				frame->fragment(j).apply( frag_pose, 1, Nmers_size );
				core::Size radius = get_radius( frag_pose ); // call the get radius by giving the frag_pose
				TR << "radius" << RJ(6, Nmers_size) << RJ(6,seq_pos) << RJ(6,j) << RJ(6,radius) << std::endl;
				container.push_back( radius );
			} // loop over "the frame" given a seq_pos
		}
	}
	return get_median( container );
}



// 130131 adding a new variable "keep": how many results you are going to keep for one round of grid search
void
density_grid_search(
	core::Real const & no_density_score,
	core::scoring::electron_density::ElectronDensity & density,
	core::pose::Pose & pose,
	core::Size nRsteps,
	double & map_s,
	double & map_s2,
	double & pose_s,
	double & pose_s2,
	core::Size list_size,
	core::Real movestep,
	core::Size edge,
	std::list< Results_Keeper > & keeper_list
){
	func << "density_grid_search function gets called: " << std::endl;
	// hacky way to get private data from ElectronDensity Class
	//numeric::xyzVector<core::Real> centerOfMass   = density.call_getCoM();
	//numeric::xyzVector<core::Real> origin         = density.call_getOrigin();
	//func << "center of Mass: " << centerOfMass << std::endl;
	//func << "origin        : " << origin << std::endl;

	numeric::xyzVector<int> grid           = density.get_grid();
	func << "grid: " << grid << std::endl;

	// define the MyCenterOfMass as the origin
	numeric::xyzVector<core::Real> MyCenterOfMass;
	//func << "MyCenterOfMass: " << MyCenterOfMass << std::endl;

	// get parameters that are going to be used in search
	numeric::xyzVector< core::Real > pre_trans;
	ObjexxFCL::FArray3D< double > sigR, poseCoefR, poseCoefI;

	// inputs: nRsteps, poseCoefR, poseCoefI
	// outputs: pre_trans( center of mass ) sigR
	func << "get parameters from map and pose - call pose_spherical_harmonic " << std::endl;
	pre_trans = density.poseSHT( pose, nRsteps, sigR, poseCoefR, poseCoefI ); // get center of mass

	// begin the grid search
	func << "-------------- grid search begins --------------" << std::endl;
	for ( core::Real i=1+edge; i<=grid[0]-edge; i+=movestep ) {
		for ( core::Real j=1+edge; j<=grid[1]-edge; j+=movestep ) {
			for ( core::Real k=1+edge; k<=grid[2]-edge; k+=movestep ) {

				MyCenterOfMass[0] = i;
				MyCenterOfMass[1] = j;
				MyCenterOfMass[2] = k;

				assert( i > 0 ); assert( j > 0 ); assert( k > 0 );

				func << "Searching at MyCenterOfMass (" << i << "," << j << "," << k << ")" << std::endl;
				numeric::xyzMatrix< core::Real > rotation;
				numeric::xyzVector< core::Real > post_trans;

				// get the correlation score: between pose and map at that grid point
				// outputs are passed by reference, and this function will change the variables for you
				//   input: nRstep, MyCenterOfMass, poseCoefR, poseCoefI,
				//   outputs: map_s, map_s2, pose_s, pose_s2, sigR, rotation, pre_trans, post_trans,
				//core::Real correlation_score = 0.5;
				core::Real correlation_score = density.mapSHT( no_density_score,
					nRsteps,
					map_s, map_s2, pose_s, pose_s2,
					poseCoefR, poseCoefI,
					rotation, pre_trans, post_trans,
					MyCenterOfMass ); // pre_trans - const

				// if the correlation_score < 0 -> the map area didn't pass the threshold just go to next grid point
				if ( correlation_score == 0.0 ) {
					continue;
				}

				// the keeper_list is designed to be the smallest number at "front" and the biggest number at "back":
				// sort based on correlation_score
				// front | 1 | 4 | 100 | 123 | back
				// 130327 Yuan helped add a check to prevent potential problem caused by uninitialized keeper_list
				//
				std::list< Results_Keeper >::iterator back_element_it  = keeper_list.end();
				std::list< Results_Keeper >::iterator front_element_it = keeper_list.begin();

				Results_Keeper results = Results_Keeper( correlation_score, pre_trans, post_trans, rotation, MyCenterOfMass, map_s, map_s2, pose_s, pose_s2 );

				if ( back_element_it==front_element_it || correlation_score >= (*(--back_element_it)).score_ ) { // if the score is higher than the highest one, just adds to the back directly
					keeper_list.push_back( results );
					//if ( keeper_list.size() > list_size ){ keeper_list.pop_front(); }

				} else if ( correlation_score >  (*front_element_it).score_ ) { // if the score is higher than the lowest, but lower than the highest one - try to find a place to insert the results
					for ( std::list< Results_Keeper >::iterator it=keeper_list.begin();
							it != keeper_list.end(); ++it ) {
						if ( correlation_score < (*it).score_ ) {
							keeper_list.insert( it, results );
							//if ( keeper_list.size() > list_size ){ keeper_list.pop_front(); } //
							break; // break the iteration once insert the value
						} // find the best place to insert "results"
					} // iterator

				} else if ( correlation_score < (*front_element_it).score_ && keeper_list.size() < list_size ) {
					keeper_list.push_front( results );
				} // keep adding items if not reach the max size
				if ( keeper_list.size() > list_size ) { keeper_list.pop_front(); } //
			}
		}
	}
	//return keeper_list;
}



// move the pose given rotation/translation matrix
void
move_it(
	core::pose::Pose & pose,
	numeric::xyzVector<core::Real> & pre_trans,
	numeric::xyzVector<core::Real> & translation,
	numeric::xyzMatrix<core::Real> & rotation
){
	// then each atom x can be transformed to this optimal configuration by:
	utility::vector1< core::id::AtomID > ids;
	utility::vector1< numeric::xyzVector< core::Real > > positions;
	numeric::xyzVector< core::Real > xyz_rot;

	for ( core::Size irsd = 1; irsd <= pose.size(); ++irsd ) {
		for ( core::Size iatom = 1; iatom <= pose.residue_type( irsd ).natoms(); ++iatom ) {
			numeric::xyzVector< core::Real > atom_xyz = pose.xyz( core::id::AtomID( iatom, irsd ) );
			ids.push_back( core::id::AtomID( iatom, irsd ) );

			xyz_rot = rotation*( atom_xyz + pre_trans ) + translation;

			positions.push_back( xyz_rot );
		}
	}
	pose.batch_set_xyz( ids, positions );
}




// main
int main(int argc, char* argv[]) {
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::pose;

		ThisApplication::register_options();
		devel::init(argc, argv);

		// read in map file and get some data from map
		//std::string mapfile = option[ edensity::mapfile ]();

		// create ElectronDensity object:
		// this is going to call "init()" in the class
		core::scoring::electron_density::ElectronDensity &density = core::scoring::electron_density::getDensityMap();
		//density.readMRCandResize( mapfile );
		core::Real map_stdev = density.getStdev();
		//core::Real map_mean  = density.getMean();
		core::Real map_min   = density.getMin();
		numeric::xyzVector<int> grid = density.get_grid();
		TR << "main grid: " << grid << std::endl;
		//TR << "mapfile: " << mapfile << " loaded." << std::endl;



		//// Do we have a map?
		//if (!
		//    utility_exit_with_message("Density scoring function called but no map loaded.");
		//}

		// get sequence from fasta or native pdb if provided
		std::string sequence;
		core::pose::Pose native_pose;
		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_file( native_pose, option[ in::file::native ]().name() , core::import_pose::PDB_file);
			sequence = native_pose.sequence();
		} else {
			sequence = core::sequence::read_fasta_file( option[ in::file::fasta ]()[1])[1]->sequence();
		}

		// read in fragments file
		// use top25 fragments for each position
		// annotated = true; to return pdbid
		core::fragment::FragSetCOP fragments;
		fragments = core::fragment::FragmentIO( basic::options::option[ num_frags ], 1, true ).read_data( basic::options::option[ fragfile ]() );

		utility::vector1<core::Size> const designated_positions = basic::options::option[ pos ]();
		if ( designated_positions[1] == 0 ) {
			utility_exit_with_message("ERROR: you have to provide residue positions eg. -pos: 1 2 3 4 5");
		} else {
			TR << "positions choose: " << designated_positions << std::endl;
		}

		utility::vector1<core::Size> const designated_picker_rank = basic::options::option[ designated_rank ]();

		core::scoring::ScoreFunctionOP scorefxn_dens( new core::scoring::ScoreFunction() );
		scorefxn_dens->set_weight( core::scoring::elec_dens_whole_structure_allatom, 1.0);

		// dump fragment pdb as silent file
		std::string silent_fn = option[ out::file::silent ]();
		core::io::silent::SilentFileData silent_file_data( silent_fn, false, false, "binary" ); //true to store argv in silent file

		// debug print
		//if ( option[ debug ].user() ){
		utility::io::ozstream out( option[ debug_outfn ] );
		//}

		//
		MinPackMin mpm;

		for ( core::fragment::ConstFrameIterator frame=fragments->begin(), eframe=fragments->end();
				frame != eframe;
				++frame ) {

			core::Size seq_pos = frame->start();
			if ( std::find( designated_positions.begin(), designated_positions.end(), seq_pos ) !=designated_positions.end() ) {

				core::Size Nmers_size = frame->length();
				std::string frag_seq = sequence.substr( seq_pos-1, Nmers_size );

				TR << "sequence: " << sequence << std::endl;
				TR << "frag_seq: " << RJ( seq_pos+frag_seq.size()-1, frag_seq ) << std::endl;

				// create native_frag_pose from native_pose if native are provided
				core::pose::Pose native_frag_pose;
				if ( option[ in::file::native ].user() ) {
					utility::vector1< core::Size > positions;
					core::kinematics::FoldTree fold_tree( Nmers_size );
					for ( core::Size irsd=seq_pos; irsd<seq_pos+Nmers_size; ++irsd ) { positions.push_back( irsd ); } // get the residues numbers to make native fragment
					core::pose::create_subpose( native_pose, positions, fold_tree, native_frag_pose ); // make native_frag_pose
				}

				// 130320 bug caught:
				//   WAS: universal approximated radius coming from calculating radius for the whole fragments set
				//   NOW: should be local universal radius from a fragment set given a position
				//        the thought is to make the correlation scores from a fragment set ( given a position) comparable
				// 130320_comment:
				//  this is actually doesn't make that much sense since different shape should give different sampling space to calculate spherical harmonic
				core::Size nRsteps = 0;
				if ( option[ median_radius ].user() ) {
					nRsteps = get_radius_in_frag_set( fragments, sequence, seq_pos );
					TR << "median nRsteps is set from fragments given position: " << nRsteps << std::endl;
				}

				// each frag candidate at a designated position
				core::pose::Pose frag_pose;
				core::pose::make_pose_from_sequence( frag_pose, frag_seq, *( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID )) );
				for ( core::Size j=1; j<=frame->nr_frags(); ++j ) {


					core::Size picker_rank = j;
					frame->fragment( picker_rank ).apply( frag_pose, 1, Nmers_size );
					std::string const pdbid ( frame->fragment( picker_rank ).pdbid() );

					std::string frag_fn = "frags."+int2str(seq_pos)+"."+int2str(picker_rank)+".pdb";

					// syd has an evict issue
					if ( utility::file::file_exists( frag_fn ) ) {
						TR << frag_fn << " exists!" << std::endl;
						continue;
					}

					if ( option[ designated_rank ].user() ) {
						if ( ! ( std::find( designated_picker_rank.begin(), designated_picker_rank.end(), picker_rank ) != designated_picker_rank.end() ) ) {
							TR << "skip " << picker_rank << std::endl;
							continue;
						}
					}

					TR << RJ(6, Nmers_size) << RJ(6,seq_pos) << RJ(6,picker_rank) << std::endl;

					// set input variables
					//core::Real score;
					double map_s=0, map_s2=0, pose_s=0, pose_s2=0;
					core::Size edge = basic::options::option[ skip_edge ]();
					TR << edge << " is going to be skipped at the edge." << std::endl;

					if ( ! option[ median_radius ].user() ) {
						nRsteps = get_radius( frag_pose );
						TR << "get nRsteps from individule frag_pose: " << nRsteps << std::endl;
					}
					assert( nRsteps != 0 );

					std::list< Results_Keeper > results_list;

					// density grid search, pass the results_list into as a reference
					density_grid_search( option[ no_density_score ](), density, frag_pose, nRsteps, map_s, map_s2, pose_s, pose_s2, basic::options::option[ keep ](), basic::options::option[ movestep ](), edge, results_list  );

					// iterate the results_list
					core::Size cor_rank = results_list.size(); // correlation score ranking

					// exit if no results passed the filter no_density_score
					if ( cor_rank < 1 ) { utility_exit_with_message("you have no results passed the filter you set."); }

					for ( std::list< Results_Keeper >::iterator it=results_list.begin();
							it != results_list.end(); ++it ) {

						// retrieve data from results_list
						core::Real cor_score                          = (*it).score_ ;
						numeric::xyzVector<core::Real> pre_trans      = (*it).pre_trans_;
						numeric::xyzVector<core::Real> translation    = (*it).translation_;
						numeric::xyzMatrix<core::Real> rotation       = (*it).rotation_;
						numeric::xyzVector<core::Real> mycenterofmass = (*it).mycenterofmass_;
						double map_s                                  = (*it).map_s_;
						double map_s2                                 = (*it).map_s2_;
						double pose_s                                 = (*it).pose_s_;
						double pose_s2                                = (*it).pose_s2_;

						// frag_pose are all the same!
						core::pose::Pose temp_pose = frag_pose;

						// rotate and translate to the position
						move_it( temp_pose, pre_trans, translation, rotation );

						// watch out for this - without native pose it would just print out 0.000
						core::Real before_super_rmsd        = 0.0;
						core::Real before_no_super_rmsd     = 0.0;
						core::Real after_super_rmsd         = 0.0;
						core::Real after_no_super_rmsd      = 0.0;
						core::Real after_density_score_fast = 0.0;

						// before refinement rmsd
						if ( option[ in::file::native ].user() ) {
							before_super_rmsd    = core::scoring::CA_rmsd( native_frag_pose, temp_pose );
							before_no_super_rmsd = core::scoring::rmsd_no_super( temp_pose, native_frag_pose, core::scoring::is_protein_CA );
							TR << "calculate fragfile rmsd to native: " << before_no_super_rmsd << std::endl;
						}

						// density score before minimization -> the scorefxn instance has existsed
						core::pose::addVirtualResAsRoot( temp_pose );
						core::Real before_density_score = (*scorefxn_dens)( temp_pose );

						// refine it!
						if ( option[ min_pack_min ].user() ) {
							for ( core::Size i=1; i<= core::Size( option[ cycles ]() ); i++ ) {

								// Step 1. rigid body minimization
								mpm.rigid_body_minimization( temp_pose );

								// Step 2. pack sidechains
								mpm.pack_sidechains( temp_pose );

								// Step 3. rigid body minimization
								mpm.rigid_body_minimization( temp_pose );

								// Step 4. backbone torsion minimization
								if ( option[ min_bb ].user() ) {
									mpm.backbone_minimization( temp_pose );
								}
							}
						}

						// density score before minimization
						core::pose::addVirtualResAsRoot( temp_pose );
						core::Real after_density_score = (*scorefxn_dens)( temp_pose );

						// after refinement rmsd
						if ( option[ in::file::native ].user() ) {
							remove_all_virtual_residues( temp_pose );
							after_super_rmsd    = core::scoring::CA_rmsd( native_frag_pose, temp_pose );
							after_no_super_rmsd = core::scoring::rmsd_no_super( temp_pose, native_frag_pose, core::scoring::is_protein_CA );
							TR << "calculate no_super_rmsd: " << after_no_super_rmsd <<std::endl;
						}

						std::string tag = "after_rotation_frags."+int2str( Nmers_size )+"."+int2str( seq_pos )+"."+int2str( picker_rank )+"."+int2str( cor_rank )+"."+pdbid;

						core::io::silent::BinarySilentStruct silent_stream( temp_pose, tag );
						silent_stream.add_energy( "Nmers_size",               Nmers_size);
						silent_stream.add_energy( "seq_pos",                  seq_pos);
						silent_stream.add_energy( "picker_rank",              picker_rank );
						silent_stream.add_energy( "cor_rank",                 cor_rank );
						silent_stream.add_energy( "cor_score",                cor_score );
						silent_stream.add_energy( "before_density_score",     before_density_score );
						silent_stream.add_energy( "before_super_rmsd",        before_super_rmsd );
						silent_stream.add_energy( "before_no_super_rmsd",     before_no_super_rmsd );
						silent_stream.add_energy( "after_density_score",      after_density_score );
						silent_stream.add_energy( "after_super_rmsd",         after_super_rmsd );
						silent_stream.add_energy( "after_no_super_rmsd",      after_no_super_rmsd );
						silent_stream.add_energy( "after_density_score_fast", after_density_score_fast );
						silent_stream.add_energy( "map_s",                    map_s );
						silent_stream.add_energy( "map_s_divide_by_stdev",    map_s/map_stdev );
						silent_stream.add_energy( "map_z_score",              (map_s-map_min)/map_stdev );
						silent_stream.add_energy( "map_s2",                   map_s2 );
						silent_stream.add_energy( "pose_s",                   pose_s );
						silent_stream.add_energy( "pose_s2",                  pose_s2 );

						silent_file_data.write_silent_struct( silent_stream, silent_fn );

						--cor_rank;

						if ( option[ out::file::scorefile ].user() ) {
							silent_file_data.write_silent_struct( silent_stream, option[ out::file::scorefile ](), true );
						}

						if ( option[ debug ].user() ) {
							out << tag
								<< " pre_trans: ( "
								<< F(10, 4, pre_trans[0]) << ", "
								<< F(10, 4, pre_trans[1]) << ", "
								<< F(10, 4, pre_trans[2]) << " ) "

								<< " translation: ( "
								<< F(10, 4, translation[0]) << ", "
								<< F(10, 4, translation[1]) << ", "
								<< F(10, 4, translation[2]) << " ) "

								<< " mycenterofmass: ( "
								<< RJ(6, mycenterofmass[0]) << ", "
								<< RJ(6, mycenterofmass[1]) << ", "
								<< RJ(6, mycenterofmass[2]) << " )"

								<< " after_no_super_rmsd: "
								<< after_no_super_rmsd

								<< std::endl;
						}


					} // results_list iterator
					frag_pose.dump_pdb( frag_fn );

				} // positional candidate check
			} // position check -> only work on some positions
		} // residual fragments check
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
