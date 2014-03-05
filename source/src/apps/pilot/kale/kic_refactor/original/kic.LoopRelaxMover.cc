
if ( refine() == "refine_kic" ) {
	LoopMover_Refine_KIC refine_kic( loops, fa_scorefxn_ );
	refine_kic.set_native_pose( new core::pose::Pose ( native_pose ) );
	refine_kic.apply( pose );
}
else if ( refine() == "refine_kic_refactor" ) {
	LoopProtocolOP protocol = new LoopProtocol;
	LoopMoverOP mover = new LoopMover;
	LoggerOP acceptance_logger = new AcceptanceRates;
	loops::Loop loop = (*loops)[1];

	Size repack_period = 20;
	if (option[OptionKeys::loops::repack_period].user()) {
		repack_period = option[OptionKeys::loops::repack_period]();
	}

	mover->add_task(new KicSampler(acceptance_logger));
	mover->add_task(new PeriodicTask(new RepackingRefiner, repack_period));
	mover->add_task(new RotamerTrialsRefiner);
	mover->add_task(new LocalMinimizationRefiner);

	Size outer_cycles = 3;
	Size inner_cycles = 10 * loop.size();

	if (option[OptionKeys::loops::outer_cycles].user()) {
		outer_cycles = option[ OptionKeys::loops::outer_cycles ]();
	}
	if (option[OptionKeys::loops::max_inner_cycles].user()) {
		Size max_cycles = option[OptionKeys::loops::max_inner_cycles]();
		inner_cycles = std::max(inner_cycles, max_cycles);
	}
	if (option[OptionKeys::loops::fast]) {
		outer_cycles = 3;
		inner_cycles = 12;
	}

	protocol->set_mover(mover);
	protocol->set_loop(loop);
	protocol->set_score_function(fa_scorefxn_);
	protocol->set_iterations(outer_cycles, inner_cycles, 2);
	protocol->add_logger(new ProgressBar);
	protocol->add_logger(new ScoreVsRmsd(native_pose, loop));
	protocol->add_logger(acceptance_logger);

	protocol->apply(pose);
}
