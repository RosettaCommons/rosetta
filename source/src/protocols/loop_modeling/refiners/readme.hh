/// @namespace protocols::loop_modeling::refiners
///
/// @brief Algorithms for lowering loop region scores during sampling.
/// 
/// @details the classes in the namespace are broadly supposed to be 
/// responsible for lowering the score in loop regions sampled using the 
/// algorithms in the samplers namespace.  The primary refinement algorithms 
/// are RepackingRefiner, RotamerTrialsRefiner, and LocalMinimizationRefiner.  
/// These algorithms together account for the lion's share of the total running 
/// time in most loop modeling simulations.
