// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   demo/ian_test/ligand_dock.cc
///
/// @brief
/// @author Ian Davis (ian.w.davis@gmail.com)

// must be here to avoid VC++ ambiguous symbol w/ ObjexxFCL::byte
// for boinc builds - dek
#include <protocols/jobdist/JobDistributors.hh>


#include <devel/init.hh>
#include <core/types.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>

#include <protocols/enzdes/EnzConstraintIO.hh> //for addding constraints if demanded by user
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/ligand_docking/LigandDockProtocol.hh>
#include <protocols/moves/Mover.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/tools/make_vector1.hh>

#include <ctime>
#include <fstream>
#include <sstream>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>




//////////////////////////////////////////////////////////////////////////////


void print_node(
	std::ostream & out,
	int residue_num,
	int atom_num,
	core::conformation::Conformation const & conf,
	std::string extras = "" //< P for points, color, width/radius, etc.
)
{
	// atom_num is often 0 in fold tree, means no specific atom.
	// might as well use the first one:
	if (atom_num == 0) atom_num = 1;
	core::conformation::Residue const & res = conf.residue(residue_num);
	core::chemical::ResidueType const & res_type = conf.residue_type(residue_num);
	core::conformation::Atom const & atom = res.atom(atom_num);
	core::chemical::AtomType const & atom_type = res.atom_type(atom_num);
	// This info appears when you click on the point
	out << "{" << res_type.name3() << " " << res.seqpos()
		<< " " << res.atom_name(atom_num) << " (" << atom_type.name() << ")"
		<< "}";
	// Color, width, etc. followed by coordinates
	//out << "col" << residue_num << " ";
	out << extras;
	out << " " << atom.xyz().x() << " " << atom.xyz().y() << " " << atom.xyz().z() << "\n";
}


void print_node(
	std::ostream & out,
	int residue_num,
	std::string atom_name,
	core::conformation::Conformation const & conf,
	std::string extras = "" //< P for points, color, width/radius, etc.
)
{
	// atom_num is often 0 in fold tree, means no specific atom.
	// might as well use the first one:
	core::conformation::Residue const & res = conf.residue(residue_num);

	int atom_num;
	if (atom_name == "") {
		atom_num = 1;
	}	else {
		atom_num = res.atom_index( atom_name );
	}
	print_node( out, residue_num, atom_num, conf, extras );
}


void print_interres_bond(
	std::ostream & out,
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::conformation::Conformation const & conf
)
{
	print_node(out, rsd1.seqpos(), rsd1.connect_atom(rsd2), conf, "P");
	print_node(out, rsd2.seqpos(), rsd2.connect_atom(rsd1), conf);
}


void dump_residue_kinemage(
	std::ostream & out,
	core::conformation::Residue const & rsd,
	core::conformation::Conformation const & conf,
	std::string tag = ""
)
{
	// intra-residue connections
	// do residues in different (~random) colors to help distinguish them
	int const num_colors = 6;
	std::string colors[num_colors] = {"pinktint", "peachtint", "yellowtint", "greentint", "bluetint", "lilactint"};
	std::string color = colors[ rsd.seqpos() % num_colors ];
	out << "@vectorlist {} color= " << color << " master= {intra-res}\n";
	for(core::Size atom_i = 1; atom_i <= rsd.natoms(); ++atom_i) {
		core::conformation::Residue::AtomIndices const & nbrs = rsd.nbrs(atom_i);
		for(core::conformation::Residue::AtomIndices::const_iterator j = nbrs.begin(), end_j = nbrs.end(); j != end_j; ++j) {
			core::Size atom_j = *j;
			if(atom_j <= atom_i) continue; // so we draw each bond just once, not twice
			bool const is_H = rsd.atom_is_hydrogen(atom_j) || rsd.atom_is_hydrogen(atom_i);
			std::string const ptmaster = ( is_H ? " 'h'" : "" );
			print_node(out, rsd.seqpos(), atom_i, conf, tag+"P"+ptmaster);
			print_node(out, rsd.seqpos(), atom_j, conf, tag+ptmaster);
		}
	}
	// inter-residue connections
	// there *has* to be a better way of getting next/prev residue...
	out << "@vectorlist {} color= gray master= {inter-res}\n";
	core::chemical::ResidueType const & res_type = rsd.type();
	if (rsd.seqpos() > 1 && rsd.is_bonded( conf.residue(rsd.seqpos()-1) )) {
		print_interres_bond(out, rsd, conf.residue(rsd.seqpos()-1), conf);
	}
	if ((core::Size)rsd.seqpos() < conf.size() && rsd.is_bonded( conf.residue(rsd.seqpos()+1) )) {
		print_interres_bond(out, rsd, conf.residue(rsd.seqpos()+1), conf);
	}
	for(core::Size i = 1; i <= res_type.n_residue_connections(); ++i) {
		print_interres_bond(out, rsd, conf.residue( rsd.residue_connection_partner(i) ), conf);
	}
}


void dump_structure_kinemage(
	std::ostream & out,
	core::conformation::Conformation const & conf
)
{
	out << "@subgroup {by residue} dominant\n";
	for(core::Size i = 1; i <= conf.size(); ++i) {
		dump_residue_kinemage(out, conf.residue(i), conf);
	}
}


void dump_energy_kinemage(
	std::ostream & out,
	core::pose::Pose const & pose,
	utility::vector1< core::scoring::ScoreType > scoretypes,
	std::string tag,
	bool animate = true
)
{
	using namespace core::scoring;
	core::Size const ligid = pose.total_residue();
	utility::vector1< core::Real > scores;
	core::Real sum_scores = 0;
	for( core::Size i = 1, end_i = pose.total_residue(); i <= end_i; ++i ) {
		EnergyEdge const * e = pose.energies().energy_graph().find_energy_edge(ligid, i);
		scores.push_back(0);
		if( !e ) continue;
		typedef utility::vector1< core::scoring::ScoreType > ST;
		for( ST::const_iterator j = scoretypes.begin(), end_j = scoretypes.end(); j != end_j; ++j ) {
			scores[i] += e->energy_map()[ *j ];
		}
		sum_scores += scores[i];
	}
	core::Real const smin = utility::min(scores), smax = utility::max(scores), absmax = std::max(std::abs(smin), std::abs(smax));
	for( core::Size i = 1, end_i = pose.total_residue(); i <= end_i; ++i ) {
		if( i == ligid ) {
			out << "@colorset {" << tag << i << "} green\n";
			continue;
		}
		out << "@hsvcolor {" << tag << i << "} " << (scores[i] > 0 ? "20" : "275") << " " << 100*std::abs(scores[i])/absmax << " 90\n";
	}

	out << "@group {" << tag << "} dominant";
	if(animate) out << " animate";
	out << "\n";
	for(core::Size i = 1; i <= pose.total_residue(); ++i) {
		std::ostringstream s;
		s << tag << i << ' ';
		dump_residue_kinemage(out, pose.residue(i), pose.conformation(), s.str());
	}
	out << "@labellist {} color= magenta master= {labels}\n";
	for(core::Size i = 1; i <= pose.total_residue(); ++i) {
		if( std::abs(scores[i]) < 0.01 ) continue;
		out << "{  ";
		std::streamsize p = out.precision();  out.precision(2);  out << scores[i]; out.precision(p);
		out << "}" << tag << i
			<< " " << pose.residue(i).nbr_atom_xyz().x()
			<< " " << pose.residue(i).nbr_atom_xyz().y()
			<< " " << pose.residue(i).nbr_atom_xyz().z()
			<< "\n";
	}
	out << "{  ";
	std::streamsize p = out.precision();  out.precision(3);  out << sum_scores; out.precision(p);
	out << "}" << tag << ligid
		<< " " << pose.residue(ligid).nbr_atom_xyz().x()
		<< " " << pose.residue(ligid).nbr_atom_xyz().y()
		<< " " << pose.residue(ligid).nbr_atom_xyz().z()
		<< "\n";

	out << "@pointmaster 'h' {Hs} off\n";
}


//////////////////////////////////////////////////////////////////////////////
/// @details Assumes option system has already been initialized!
int
main( int argc, char * argv [] )
{
	try {

	// Parses command line options and inits RNG.
	// Doesn't seem to hurt to do it again if already done once (?)
	// Except in unit testing mode, where it wipes out e.g. -database
	devel::init(argc, argv);

	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using namespace protocols::ligand_docking;
	basic::Tracer TR("ligand_dock.main");

	// Build overall docking protocol Mover
	LigandDockProtocolOP dockingProtocol = new LigandDockProtocol();

	time_t overall_start_time = time(NULL);
	bool const use_silent_in = option[ in::file::silent ].active();
	utility::vector1< BasicJobOP > input_jobs;
	core::io::silent::AtomTreeDiffOP atdiff;
	if( use_silent_in ) {
		using core::io::silent::AtomTreeDiff;
		int const nstruct = std::max( 1, option[ out::nstruct ]() );
		atdiff = new AtomTreeDiff( *(option[ in::file::silent ]().begin()) );
		if( (option [ in::file::silent ]()).size() > 1 ) utility_exit_with_message("ligand_dock.main can only handle one input silent file at a time!");
		if( option[ in::file::tags ].active() ) {
			// Do only the ones specified in tags.  If it's a number, try it as an index too.
			utility::vector1< std::string > tags = option[ in::file::tags ]();
			for(core::Size i = 1; i <= tags.size(); ++i) {
				core::Size tag_as_int = 0;
				if( is_int(tags[i]) ) tag_as_int = int_of(tags[i]);
				if( atdiff->has_tag( tags[i] ) ) {
					input_jobs.push_back(new BasicJob( tags[i], tags[i], nstruct ));
				} else if( 0 < tag_as_int && tag_as_int <= atdiff()->scores().size() ) {
					input_jobs.push_back(new BasicJob( atdiff()->scores()[tag_as_int].first, atdiff()->scores()[tag_as_int].first, nstruct ));
				} else {
					utility_exit_with_message("Can't find tag "+tags[i]);
				}
			}
		} else {
			// Just do them all!
			for(core::Size i = 1; i <= atdiff()->scores().size(); ++i) {
				input_jobs.push_back(new BasicJob( atdiff()->scores()[i].first, atdiff()->scores()[i].first, nstruct ));
			}
		}
	} else {
		input_jobs = load_s_and_l();
		// Reduce read contention between processes by randomizing the order in which structures are processed
		// Don't want to do this with silent file input -- slows access and screws up reference structure tracking.
		numeric::random::random_permutation( input_jobs, numeric::random::RG );
	}
	std::string outfile = "silent.out";
	if( option[ out::file::silent].user() ) outfile = option[ out::file::silent];
	AtomTreeDiffJobDistributor< BasicJobOP > jobdist( input_jobs, outfile );
	jobdist.set_precision(6,3,1); // makes silent file much smaller, ~3x vs. default 6,4,2

	protocols::enzdes::EnzConstraintIOOP constraint_io = NULL;
	if( option[basic::options::OptionKeys::enzdes::cstfile].user() ){
		//we need the residue type set, assuming FA standard is used
		core::chemical::ResidueTypeSetCAP restype_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		option[basic::options::OptionKeys::run::preserve_header ].value(true);
		constraint_io = new protocols::enzdes::EnzConstraintIO( restype_set );
		constraint_io->read_enzyme_cstfile(basic::options::option[basic::options::OptionKeys::enzdes::cstfile]);
		core::scoring::ScoreFunctionOP scorefunc = dockingProtocol->get_scorefxn();
		if (scorefunc->has_zero_weight( core::scoring::coordinate_constraint ) ) scorefunc->set_weight(core::scoring::coordinate_constraint, 1.0 ) ;
		if (scorefunc->has_zero_weight( core::scoring::atom_pair_constraint ) ) scorefunc->set_weight(core::scoring::atom_pair_constraint, 1.0 ) ;
		if (scorefunc->has_zero_weight( core::scoring::angle_constraint ) ) scorefunc->set_weight(core::scoring::angle_constraint, 1.0 ) ;
		if (scorefunc->has_zero_weight( core::scoring::dihedral_constraint ) ) scorefunc->set_weight(core::scoring::dihedral_constraint, 1.0 ) ;
	}

	BasicJobOP curr_job, prev_job;
	int curr_nstruct = 0, num_structures_processed = 0;
	core::pose::PoseOP input_pose; // starts NULL, coords *never* modified!
	jobdist.startup();
	while( jobdist.next_job(curr_job, curr_nstruct) ) {
		time_t pdb_start_time = time(NULL);
		TR << "Starting " << curr_job->output_tag(curr_nstruct) << " ..." << std::endl;

		// we read each PDB just once to save on disk I/O
		if( curr_job.get() != prev_job.get() || input_pose.get() == NULL ) {
			input_pose = new core::pose::Pose();
			if( use_silent_in ) atdiff->read_pose( curr_job->input_tag(), *input_pose );
			else core::import_pose::pose_from_pdb( *input_pose, curr_job->input_tag() );

			//if constraints are requested
			if( option[basic::options::OptionKeys::enzdes::cstfile].user() ){
				constraint_io->add_constraints_to_pose( *input_pose, dockingProtocol->get_scorefxn(), false );
			}
		}

		// Make a modifiable copy of the pose read from disk
		core::pose::PoseOP the_pose = new core::pose::Pose( *input_pose );
		the_pose->data().set(basic::JOBDIST_OUTPUT_TAG, new basic::CacheableString(curr_job->output_tag(curr_nstruct)));

		//*** dockingProtocol->apply( *the_pose ); ***
		//std::ostream & out = std::cout;
		std::ofstream out("tmp.kin");
		using namespace core::scoring;
		core::scoring::ScoreFunctionOP sfxn = dockingProtocol->get_scorefxn();
		(*sfxn)( *the_pose );
		/// Now handled automatically.  sfxn->accumulate_residue_total_energies( *the_pose );
		out << "@kinemage {" << curr_job->input_tag() << "}\n";
		out << "@title {" << curr_job->input_tag() << "}\n";
		core::Size const ligid = the_pose->total_residue();
		core::Vector const & ctr = the_pose->residue(ligid).nbr_atom_xyz();
		out << "@1center " << ctr.x() << " " << ctr.y() << " " << ctr.z() << "\n";
		out << "@1span 25\n";
		dump_energy_kinemage(out, *the_pose, utility::tools::make_vector1(fa_atr, fa_rep, fa_sol, hbond_sc, hbond_bb_sc, hbond_lr_bb, hbond_sr_bb, fa_elec), "a+r+s+h+e");
		dump_energy_kinemage(out, *the_pose, utility::tools::make_vector1(fa_atr, fa_rep), "a+r");
		dump_energy_kinemage(out, *the_pose, utility::tools::make_vector1(fa_sol, hbond_sc, hbond_bb_sc, hbond_lr_bb, hbond_sr_bb, fa_elec), "s+h+e");
		dump_energy_kinemage(out, *the_pose, utility::tools::make_vector1(fa_atr), "atr");
		dump_energy_kinemage(out, *the_pose, utility::tools::make_vector1(fa_rep), "rep");
		dump_energy_kinemage(out, *the_pose, utility::tools::make_vector1(fa_sol), "sol");
		dump_energy_kinemage(out, *the_pose, utility::tools::make_vector1(hbond_sc, hbond_bb_sc, hbond_lr_bb, hbond_sr_bb), "hbond");
		dump_energy_kinemage(out, *the_pose, utility::tools::make_vector1(fa_elec), "elec");
		//dump_energy_kinemage(out, *the_pose, utility::tools::make_vector1(coordinate_constraint, atom_pair_constraint, angle_constraint, dihedral_constraint), "constr"); // these values not cached?
		out.close();

		// Score new structure and add to silent file
		//std::map< std::string, core::Real > scores;
		//core::io::silent::map_of_weighted_scores(*the_pose, *(dockingProtocol->get_scorefxn()), scores);
		//if( rms_native_pose.get() == NULL ) { // input serves as native
		//	dockingProtocol->append_ligand_docking_scores(*input_pose, *the_pose,
		//		dockingProtocol->get_scorefxn(), scores, constraint_io);
		//} else {
		//	dockingProtocol->append_ligand_docking_scores(*rms_native_pose, *the_pose,
		//		dockingProtocol->get_scorefxn(), scores, constraint_io);
		//}
		//// Want to recycle the reference poses from the input silent file, or else every entry becomes a new reference pose!
		//if( use_silent_in ) {
		//	//std::cout << "Ref addr " << (atdiff->ref_pose_for( curr_job->input_tag() ))() << std::endl;
		//	jobdist.dump_pose( curr_job->output_tag(curr_nstruct), scores, *(atdiff->ref_pose_for( curr_job->input_tag() )), *the_pose );
		//}
		//else jobdist.dump_pose( curr_job->output_tag(curr_nstruct), scores, *input_pose, *the_pose );

		prev_job = curr_job; // pointer assignment, not a copy op
		num_structures_processed += 1;
		time_t pdb_end_time = time(NULL);
		TR << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;
	} // loop over jobs and nstructs

	time_t overall_end_time = time(NULL);
	TR << "Finished all " << num_structures_processed << " structures in " << (long)(overall_end_time - overall_start_time) << " seconds." << std::endl;
	if ( num_structures_processed == 0 )
		basic::Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them?" << std::endl;
	jobdist.shutdown(); // under BOINC, this will cause program exit!

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}

