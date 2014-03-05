/// @namespace protocols::loop_modeling::openers
///
/// @brief Unused algorithms for created open loop conformations.
/// 
/// @details This namespace was originally intended to provide a general set of 
/// algorithms for sampling new backbone torsions.  The idea is that these 
/// algorithms not concern themselves with whether or not the backbone was 
/// closed so that they could sample freely.  Then an algorithm like KIC or CCD 
/// would be used to ensure closure.
///
/// The problem is that these methods ended up being too slow for kinematic 
/// closure, presumably because they were causing the AtomTree to recalculate 
/// its coordinates too often.  (I'm not totally sure about that, though.)  So 
/// instead, KIC uses custom perturber classes that operate directly on its 
/// internal matrix representation of the closure problem.
///
/// So these algorithms would only be used by CCD, but CCD has not yet been 
/// ported to the new loop modeling framework.  For now, then, these algorithms 
/// are not being used.

