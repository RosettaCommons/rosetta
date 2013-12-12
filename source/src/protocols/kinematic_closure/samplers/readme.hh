// #error Don't include; for documentation only.

/// @namespace protocols::kinematic_closure::samplers
///
/// @brief Simple interfaces to the kinematic closure algorithm.
///
/// @details The sampler classes provided by this namespace are meant to be 
/// used by other protocols that want to use kinematic closure functionality.  
/// In general, the classes in this namespace provide apply() functions which 
/// sample new backbone conformations and a series of helper methods which 
/// customize how apply() behaves.
///
/// Feel free to create new sampler classes for specific applications.  For 
/// example, if you decide to write custom pivot picking and perturbing 
/// algorithms that are meant to be used together, make a new sampler class 
/// that wraps around KicSampler and uses your custom algorithms by default.
///
/// @note In the future, I might decide to make everything in this namespace 
/// inherit from protocols::moves::Mover.  The only difference right now is 
/// that the apply() methods all take a loop object, but it might be better to 
/// delegate the responsibility for choosing a loop to the pivot picker 
/// classes.
