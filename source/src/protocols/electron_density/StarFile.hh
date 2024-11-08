// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/electron_density/StarFile.hh
/// @brief A holder for all of the crap that a starfile can be
/// @author Danny Farrell
/// @note I'm pretty sure this is the most unsafe file that's ever been written.

#ifndef INCLUDED_protocols_electron_density_StarFile_hh
#define INCLUDED_protocols_electron_density_StarFile_hh

#include <utility/vector1.hh>
#include <core/types.hh>

// std c++ headers
#include <iostream>
#include <iosfwd>
#include <fstream>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <map>
#include <memory>

namespace protocols {
namespace electron_density {

class StarOptions {
public:
	StarOptions() {}
	virtual ~StarOptions() = default;

	friend std::ostream& operator<<( std::ostream& o, const StarOptions& c_in ) {
		return c_in.print(o);
	}


	template<class T, class U> void set( U u );

	template<class T> bool is_set();

	template<class T> T get();

private:
	virtual std::ostream& print(std::ostream & ) const = 0;
};


template< class T >
class StarOption : public StarOptions {
public:
	StarOption( T val ) { val_ = val; }
	StarOption() {
		set_ = false;
	}

	void
	set( T val ) {
		val_ = val;
		set_ = true;
	}

	bool
	is_set() {
		return set_;
	}

	T
	get() {
		return val_;
	}

private:
	std::ostream& print(std::ostream& o) const {
		//return o << utility::to_string(val_);
		return o << val_;
	}
	T val_;
	bool set_;
};

// TODO check for bad casts? to exit gracefully...
// TODO check if is an option
// TODO check if option is set
//
template<class T, class U> void StarOptions::set( U u )
{ return dynamic_cast<StarOption<T>&>(*this).set(u); }


template<class T> bool StarOptions::is_set()
{ return dynamic_cast<StarOption<T>&>(*this).is_set(); }

template<class T> T StarOptions::get()
{ return dynamic_cast<StarOption<T>&>(*this).get(); }


class StarFile {
private:
	utility::vector1< std::map< std::string, std::unique_ptr<StarOptions> > > all_image_info_;
	utility::vector1< std::string > set_values_;
public:
	
	StarFile( StarFile const & ) = delete;
        StarFile & operator=( StarFile const & ) = delete;

	template< typename F > void
	set(core::Size i, std::string option_name, F val_in ) {
		if ( i > all_image_info_.size() ) throw std::runtime_error("set StarFile out of range");
		all_image_info_[i][option_name]->set<F>(val_in);

		set_values_.push_back(option_name);
		std::sort( set_values_.begin(), set_values_.end() );
		set_values_.erase( std::unique( set_values_.begin(), set_values_.end() ), set_values_.end() );
	}

	template< typename F > bool
	is_set(core::Size i, std::string option_name) {
		return all_image_info_[i][option_name]->is_set<F>();
	}

	template< typename F > void
	get(core::Size i, std::string option_name) {
		return all_image_info_[i][option_name]->get<F>();
	}

	//std::unique_ptr<StarOptions>
	//get_StarOption(core::Size i, std::string option_name ) {
	// return all_image_info_[i][option_name];
	//}

	//std::ostream& print(core::Size i, std::string option_name,  std::ostream& o) const {
	// return os << all_image_info_[i][option_name]->print();
	//}


	//std::string
	//write_to_string() {
	// std::string ret = "data_\nloop_\n";
	// for ( core::Size i = 1; i <= set_values_.size(); ++i ) {
	//  ret += "_" + set_values_[i] + "\n";
	// }
	// for ( core::Size i = 1; i <= all_image_info_.size(); ++i ) {
	//  for ( core::Size j = 1; i <= set_values_.size(); ++j ) {
	//   ret += utility::to_string( StarFile::get(i, set_values_[j]) ) + " ";
	//  }
	//  ret += "\n";
	// }
	// return ret;
	//}

	/* Write StarFile information to filename
	*
	* @param filename
	* @param starfile_information A vector of StarFiles
	*/
	void
	writeStarFile( std::string filename ) {
		std::ofstream out_star( filename.c_str() );
		out_star << "data_\nloop_\n";
		for ( core::Size i = 1; i <= set_values_.size(); ++i ) {
			out_star << "_" << set_values_[i] << "\n";
		}
		for ( core::Size i = 1; i <= all_image_info_.size(); ++i ) {
			for ( core::Size j = 1; j <= set_values_.size(); ++j ) {
				out_star << *all_image_info_[i][set_values_[j]] << " ";
			}
			out_star << std::endl;
		}
	}


	core::Size
	size() {
		return all_image_info_.size();
	}


	StarFile() {}

	/* I'm so sorry
	*
	*/
	void
	append_one() {
		std::map< std::string, std::unique_ptr<StarOptions> > star_map;
		//rlnAccuracyRotations (RFLOAT) : Estimated accuracy (in degrees) with which rotations can be assigned
		star_map["rlnAccuracyRotations"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAdaptiveOversampleFraction (RFLOAT) : Fraction of the weights that will be oversampled in a second pass of the adaptive oversampling strategy
		star_map["rlnAdaptiveOversampleFraction"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAdaptiveOversampleOrder (int)    : Order of the adaptive oversampling (0=no oversampling, 1= 2x oversampling; 2= 4x oversampling, etc)
		star_map["rlnAdaptiveOversampleOrder"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		//rlnAmplitudeContrast (RFLOAT) : Amplitude contrast (as a fraction, i.e. 10% = 0.1)
		star_map["rlnAmplitudeContrast"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAmplitudeCorrelationMaskedMaps (RFLOAT) : Correlation coefficient between amplitudes in Fourier shells of masked maps
		star_map["rlnAmplitudeCorrelationMaskedMaps"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAmplitudeCorrelationUnmaskedMaps (RFLOAT) : Correlation coefficient between amplitudes in Fourier shells of unmasked maps
		star_map["rlnAmplitudeCorrelationUnmaskedMaps"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAnglePsi (RFLOAT) : Third Euler, or in-plane angle (psi, in degrees)
		star_map["rlnAnglePsi"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAnglePsiFlipRatio (RFLOAT) : Flip ratio of bimodal psi prior (0~0.5, 0 means an ordinary prior, 0.5 means a perfect bimodal prior)
		star_map["rlnAnglePsiFlipRatio"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAnglePsiPrior (RFLOAT) : Center of the prior (in degrees) on the third Euler angle (psi)
		star_map["rlnAnglePsiPrior"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAngleRot (RFLOAT) : First Euler angle (rot, in degrees)
		star_map["rlnAngleRot"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAngleRotPrior (RFLOAT) : Center of the prior (in degrees) on the first Euler angle (rot)
		star_map["rlnAngleRotPrior"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAngleTilt (RFLOAT) : Second Euler angle (tilt, in degrees)
		star_map["rlnAngleTilt"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAngleTiltPrior (RFLOAT) : Center of the prior (in degrees) on the second Euler angle (tilt)
		star_map["rlnAngleTiltPrior"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAngstromResolution (RFLOAT) : Resolution (in Angstroms)
		star_map["rlnAngstromResolution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAreaId (long)   : ID (i.e. a unique number) of an area (i.e. field-of-view)
		star_map["rlnAreaId"] = std::unique_ptr< StarOptions >(new StarOption<long>(0.0));
		//rlnAreaName (string) : Name of an area (i.e. field-of-view)
		star_map["rlnAreaName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		//rlnAutoLocalSearchesHealpixOrder (int)    : Healpix order (before oversampling) from which autosampling procedure will use local angular searches
		star_map["rlnAutoLocalSearchesHealpixOrder"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		//rlnAutopickFigureOfMerit (RFLOAT) : Autopicking FOM for a particle
		star_map["rlnAutopickFigureOfMerit"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAvailableMemory (RFLOAT) : Available memory per computing node (i.e. per MPI-process)
		star_map["rlnAvailableMemory"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAverageNrOfFrames (int)    : Number of movie frames that one averages over upon extraction of movie-particles
		star_map["rlnAverageNrOfFrames"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		//rlnAveragePmax (RFLOAT) : Average value (over all images) of the maxima of the probability distributions
		star_map["rlnAveragePmax"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnAverageValue (RFLOAT) : Average value for the pixels in an image
		star_map["rlnAverageValue"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnBeamTiltGroupName (string) : Name of a group (of images) with assumedly identical beam-tilts
		star_map["rlnBeamTiltGroupName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		//rlnBeamTiltX (RFLOAT) : Beam tilt in the X-direction (in mrad)
		star_map["rlnBeamTiltX"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnBeamTiltY (RFLOAT) : Beam tilt in the Y-direction (in mrad)
		star_map["rlnBeamTiltY"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnBestResolutionThusFar (RFLOAT) : The highest resolution that has been obtained in this optimization thus far
		star_map["rlnBestResolutionThusFar"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnBfactorUsedForSharpening (RFLOAT) : Applied B-factor in the sharpening of the map
		star_map["rlnBfactorUsedForSharpening"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnChangesOptimalClasses (RFLOAT) : The number of particles that changed their optimal clsas assignment in the last iteration
		star_map["rlnChangesOptimalClasses"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnChangesOptimalOffsets (RFLOAT) : The average change in optimal translation in the last iteration (in pixels)
		star_map["rlnChangesOptimalOffsets"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnChangesOptimalOrientations (RFLOAT) : The average change in optimal orientation in the last iteration (in degrees)
		star_map["rlnChangesOptimalOrientations"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnChromaticAberration (RFLOAT) : Chromatic aberration (in millimeters)
		star_map["rlnChromaticAberration"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnClassDistribution (RFLOAT) : Probability Density Function of the different classes (i.e. fraction of images assigned to each class)
		star_map["rlnClassDistribution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnClassNumber (int)    : Class number for which a particle has its highest probability
		star_map["rlnClassNumber"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		//rlnClassPriorOffsetX (RFLOAT) : Prior in the X-offset for a class (in pixels)
		star_map["rlnClassPriorOffsetX"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnClassPriorOffsetY (RFLOAT) : Prior in the Y-offset for a class (in pixels)
		star_map["rlnClassPriorOffsetY"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCoarseImageSize (int)    : Current size of the images to be used in the first pass of the adaptive oversampling strategy (may be smaller than the original image size)
		star_map["rlnCoarseImageSize"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		//rlnClassPriorOffsetX (RFLOAT) : Prior in the X-offset for a class (in pixels)
		star_map["rlnClassPriorOffsetX"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnClassPriorOffsetY (RFLOAT) : Prior in the Y-offset for a class (in pixels)
		star_map["rlnClassPriorOffsetY"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCoarseImageSize (int)    : Current size of the images to be used in the first pass of the adaptive oversampling strategy (may be smaller than the original image size)
		star_map["rlnCoarseImageSize"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		//rlnComment (string) : A metadata comment (This is treated in a special way)
		star_map["rlnComment"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		//rlnConvergenceCone (RFLOAT) : Convergence cone (in mrad)
		star_map["rlnConvergenceCone"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCoordinateX (RFLOAT) : X-Position of an image in a micrograph (in pixels)
		star_map["rlnCoordinateX"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCoordinateY (RFLOAT) : Y-Position of an image in a micrograph (in pixels)
		star_map["rlnCoordinateY"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCoordinateZ (RFLOAT) : Z-Position of an image in a 3D micrograph, i.e. tomogram (in pixels)
		star_map["rlnCoordinateZ"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps (RFLOAT) : FSC value after masking of the randomized-phases maps
		star_map["rlnCorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCorrelationFitGuinierPlot (RFLOAT) : The correlation coefficient of the fitted line through the Guinier-plot
		star_map["rlnCorrelationFitGuinierPlot"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCtfBfactor (RFLOAT) : B-factor (in A^2) that describes CTF power spectrum fall-off
		star_map["rlnCtfBfactor"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCtfDataAreCtfPremultiplied (bool)   : Flag to indicate that the input images have been premultiplied with their CTF
		star_map["rlnCtfDataAreCtfPremultiplied"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		//rlnCtfDataArePhaseFlipped (bool)   : Flag to indicate that the input images have been phase-flipped
		star_map["rlnCtfDataArePhaseFlipped"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		//rlnCtfFigureOfMerit (RFLOAT) : Figure of merit for the fit of the CTF (not used inside relion_refine)
		star_map["rlnCtfFigureOfMerit"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCtfImage (string) : Name of an image with all CTF values
		star_map["rlnCtfImage"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		//rlnCtfMaxResolution (RFLOAT) : Estimated maximum resolution (in A) of significant CTF Thon rings
		star_map["rlnCtfMaxResolution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCtfScalefactor (RFLOAT) : Linear scale-factor on the CTF (values between 0 and 1)
		star_map["rlnCtfScalefactor"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCtfValidationScore (RFLOAT) : Gctf-based validation score for the quality of the CTF fit
		star_map["rlnCtfValidationScore"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		//rlnCtfValue (RFLOAT) : Value of the Contrast Transfer Function
		star_map["rlnCtfValue"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnCurrentImageSize (int)    : Current size of the images used in the refinement
		star_map["rlnCurrentImageSize"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnCurrentIteration (int)    : The number of the current iteration
		star_map["rlnCurrentIteration"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnCurrentResolution (RFLOAT) : Current resolution where SSNR^MAP drops below 1 (in 1/Angstroms)
		star_map["rlnCurrentResolution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnDataDimensionality (int)    : Dimensionality of the data (2D/3D)
		star_map["rlnDataDimensionality"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnDataType (int)    : Type of data stored in an image (e.g. int, RFLOAT etc)
		star_map["rlnDataType"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnDefocusAngle (RFLOAT) : Angle between X and defocus U direction (in degrees)
		star_map["rlnDefocusAngle"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnDefocusU (RFLOAT) : Defocus in U-direction (in Angstroms, positive values for underfocus)
		star_map["rlnDefocusU"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnDefocusV (RFLOAT) : Defocus in V-direction (in Angstroms, positive values for underfocus)
		star_map["rlnDefocusV"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnDetectorPixelSize (RFLOAT) : Pixel size of the detector (in micrometers)
		star_map["rlnDetectorPixelSize"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnDiff2RandomHalves (RFLOAT) : Power of the differences between two independent reconstructions from random halves of the data
		star_map["rlnDiff2RandomHalves"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnDifferentialPhaseResidualMaskedMaps (RFLOAT) : Differential Phase Residual in Fourier shells of masked maps
		star_map["rlnDifferentialPhaseResidualMaskedMaps"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnDifferentialPhaseResidualUnmaskedMaps (RFLOAT) : Differential Phase Residual in Fourier shells of unmasked maps
		star_map["rlnDifferentialPhaseResidualUnmaskedMaps"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnDoAutoRefine (bool)   : Flag to indicate that 3D auto-refine procedure is being used
		star_map["rlnDoAutoRefine"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoCorrectCtf (bool)   : Flag to indicate that CTF-correction should be performed
		star_map["rlnDoCorrectCtf"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoCorrectMagnification (bool)   : Flag to indicate that (per-group) magnification correction should be performed
		star_map["rlnDoCorrectMagnification"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoCorrectNorm (bool)   : Flag to indicate that (per-image) normalisation-error correction should be performed
		star_map["rlnDoCorrectNorm"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoCorrectScale (bool)   : Flag to indicate that internal (per-group) intensity-scale correction should be performed
		star_map["rlnDoCorrectScale"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoHelicalRefine (bool)   : Flag to indicate that helical refinement should be performed
		star_map["rlnDoHelicalRefine"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoIgnoreCtfUntilFirstPeak (bool)   : Flag to indicate that the CTFs should be ignored until their first peak
		star_map["rlnDoIgnoreCtfUntilFirstPeak"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoMapEstimation (bool)   : Flag to indicate that MAP estimation should be performed (otherwise ML estimation)
		star_map["rlnDoMapEstimation"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoOnlyFlipCtfPhases (bool)   : Flag to indicate that CTF-correction should only comprise phase-flipping
		star_map["rlnDoOnlyFlipCtfPhases"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoRealignMovies (bool)   : Flag to indicate that individual frames of movies are being re-aligned
		star_map["rlnDoRealignMovies"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoSkipAlign (bool)   : Flag to indicate that orientational (i.e. rotational and translational) searches will be omitted from the refinement, only marginalisation over classes will take place
		star_map["rlnDoSkipAlign"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoSkipRotate (bool)   : Flag to indicate that rotational searches will be omitted from the refinement, only marginalisation over classes and translations will take place
		star_map["rlnDoSkipRotate"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoSolventFlattening (bool)   : Flag to indicate that the references should be masked to set their solvent areas to a constant density
		star_map["rlnDoSolventFlattening"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoSplitRandomHalves (bool)   : Flag to indicate that the data should be split into two completely separate, random halves
		star_map["rlnDoSplitRandomHalves"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoStochasticGradientDescent (bool)   : Flag to indicate that SGD-optimisation should be performed (otherwise expectation maximisation)
		star_map["rlnDoStochasticGradientDescent"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnDoZeroMask (bool)   : Flag to indicate that the surrounding solvent area in the experimental particles will be masked to zeros (by default random noise will be used
		star_map["rlnDoZeroMask"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnEnabled (bool)   : Not used in RELION, only included for backward compatibility with XMIPP selfiles
		star_map["rlnEnabled"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnEnergyLoss (RFLOAT) : Energy loss (in eV)
		star_map["rlnEnergyLoss"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnEstimatedResolution (RFLOAT) : Estimated resolution (in A) for a reference
		star_map["rlnEstimatedResolution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnExperimentalDataStarFile (string) : STAR file with metadata for the experimental images
		star_map["rlnExperimentalDataStarFile"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnFinalResolution (RFLOAT) : Final estimated resolution after postprocessing (in Angstroms)
		star_map["rlnFinalResolution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnFittedInterceptGuinierPlot (RFLOAT) : The fitted intercept of the Guinier-plot
		star_map["rlnFittedInterceptGuinierPlot"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnFittedSlopeGuinierPlot (RFLOAT) : The fitted slope of the Guinier-plot
		star_map["rlnFittedSlopeGuinierPlot"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnFixSigmaNoiseEstimates (bool)   : Flag to indicate that the estimates for the power spectra of the noise should be kept constant
		star_map["rlnFixSigmaNoiseEstimates"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnFixSigmaOffsetEstimates (bool)   : Flag to indicate that the estimates for the stddev in the origin offsets should be kept constant
		star_map["rlnFixSigmaOffsetEstimates"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnFixTauEstimates (bool)   : Flag to indicate that the estimates for the power spectra of the signal (i.e. the references) should be kept constant
		star_map["rlnFixTauEstimates"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnFourierCompleteness (RFLOAT) : Fraction of Fourier components (per resolution shell) with SNR>1
		star_map["rlnFourierCompleteness"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnFourierShellCorrelation (RFLOAT) : FSC value (of unspecified type, e.g. masked or unmasked)
		star_map["rlnFourierShellCorrelation"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnFourierShellCorrelationCorrected (RFLOAT) : Final FSC value: i.e. after correction based on masking of randomized-phases maps
		star_map["rlnFourierShellCorrelationCorrected"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnFourierShellCorrelationMaskedMaps (RFLOAT) : FSC value after masking of the original maps
		star_map["rlnFourierShellCorrelationMaskedMaps"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnFourierShellCorrelationUnmaskedMaps (RFLOAT) : FSC value before masking of the original maps
		star_map["rlnFourierShellCorrelationUnmaskedMaps"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnFourierSpaceInterpolator (int)    : The kernel used for Fourier-space interpolation (NN=0, linear=1)
		star_map["rlnFourierSpaceInterpolator"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnGoldStandardFsc (RFLOAT) : Fourier shell correlation between two independent reconstructions from random halves of the data
		star_map["rlnGoldStandardFsc"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnGroupName (string) : The name of a group of images (e.g. all images from a micrograph)
		star_map["rlnGroupName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnGroupNrParticles (long)   : Number particles in a group of images
		star_map["rlnGroupNrParticles"] = std::unique_ptr< StarOptions >(new StarOption<long>(0.0));
		// rlnGroupNumber (long)   : The number of a group of images
		star_map["rlnGroupNumber"] = std::unique_ptr< StarOptions >(new StarOption<long>(0.0));
		// rlnGroupScaleCorrection (RFLOAT) : Intensity-scale correction for a group of images
		star_map["rlnGroupScaleCorrection"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHasConverged (bool)   : Flag to indicate that the optimization has converged
		star_map["rlnHasConverged"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnHasHighFscAtResolLimit (bool)   : Flag to indicate that the FSC at the resolution limit is significant
		star_map["rlnHasHighFscAtResolLimit"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnHasLargeSizeIncreaseIterationsAgo (int)    : How many iterations have passed since the last large increase in image size
		star_map["rlnHasLargeSizeIncreaseIterationsAgo"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnHealpixOrder (int)    : Healpix order for the sampling of the first two Euler angles (rot, tilt) on the 3D sphere
		star_map["rlnHealpixOrder"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnHelicalCentralProportion (RFLOAT) : Only expand this central fraction of the Z axis when imposing real-space helical symmetry
		star_map["rlnHelicalCentralProportion"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalKeepTiltPriorFixed (bool)   : Flag to indicate that helical tilt priors are kept fixed (at 90 degrees) in global angular searches
		star_map["rlnHelicalKeepTiltPriorFixed"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnHelicalMaskTubeInnerDiameter (RFLOAT) : Inner diameter of helical tubes in Angstroms (for masks of helical references and particles)
		star_map["rlnHelicalMaskTubeInnerDiameter"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalMaskTubeOuterDiameter (RFLOAT) : Outer diameter of helical tubes in Angstroms (for masks of helical references and particles)
		star_map["rlnHelicalMaskTubeOuterDiameter"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalOffsetStep (RFLOAT) : Step size for the searches of offsets along helical axis (in Angstroms)
		star_map["rlnHelicalOffsetStep"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalRise (RFLOAT) : The helical rise (translation per subunit) in Angstroms
		star_map["rlnHelicalRise"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalRiseInitial (RFLOAT) : The initial helical rise (translation per subunit) in Angstroms before refinement
		star_map["rlnHelicalRiseInitial"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalRiseInitialStep (RFLOAT) : Initial step of helical rise search (in Angstroms)
		star_map["rlnHelicalRiseInitialStep"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalRiseMax (RFLOAT) : Maximum helical rise (in Angstroms)
		star_map["rlnHelicalRiseMax"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalRiseMin (RFLOAT) : Minimum helical rise (in Angstroms)
		star_map["rlnHelicalRiseMin"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalSigmaDistance (RFLOAT) : Sigma of distance along the helical tracks
		star_map["rlnHelicalSigmaDistance"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalSymmetryLocalRefinement (bool)   : Flag to indicate that local refinement of helical parameters should be performed
		star_map["rlnHelicalSymmetryLocalRefinement"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnHelicalTrackLength (RFLOAT) : Distance from the position of this helical segment to the starting point of the tube
		star_map["rlnHelicalTrackLength"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalTubeID (int)    : Helical tube ID for a helical segment
		star_map["rlnHelicalTubeID"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnHelicalTubePitch (RFLOAT) : Corss-over distance for a helical segment (A)
		star_map["rlnHelicalTubePitch"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalTwist (RFLOAT) : The helical twist (rotation per subunit) in degrees
		star_map["rlnHelicalTwist"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalTwistInitial (RFLOAT) : The intial helical twist (rotation per subunit) in degrees before refinement
		star_map["rlnHelicalTwistInitial"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalTwistInitialStep (RFLOAT) : Initial step of helical twist search (in degrees)
		star_map["rlnHelicalTwistInitialStep"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalTwistMax (RFLOAT) : Maximum helical twist (in degrees, + for right-handedness)
		star_map["rlnHelicalTwistMax"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHelicalTwistMin (RFLOAT) : Minimum helical twist (in degrees, + for right-handedness)
		star_map["rlnHelicalTwistMin"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHighresLimitExpectation (RFLOAT) : High-resolution-limit (in Angstrom) for the expectation step
		star_map["rlnHighresLimitExpectation"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnHighresLimitSGD (RFLOAT) : High-resolution-limit (in Angstrom) for Stochastic Gradient Descent
		star_map["rlnHighresLimitSGD"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnIgnoreHelicalSymmetry (bool)   : Flag to indicate that helical symmetry is ignored in 3D reconstruction
		star_map["rlnIgnoreHelicalSymmetry"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnImageDimensionality (int)    : Dimensionality of data stored in an image (i.e. 2 or 3)
		star_map["rlnImageDimensionality"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnImageId (long)   : ID (i.e. a unique number) of an image
		star_map["rlnImageId"] = std::unique_ptr< StarOptions >(new StarOption<long>(0.0));
		// rlnImageName (string) : Name of an image
		star_map["rlnImageName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnImageOriginalName (string) : Original name of an image
		star_map["rlnImageOriginalName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnImageSize (int)    : Size of an image (in pixels)
		star_map["rlnImageSize"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnImageSizeX (int)    : Size of an image in the X-direction (in pixels)
		star_map["rlnImageSizeX"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnImageSizeY (int)    : Size of an image in the Y-direction (in pixels)
		star_map["rlnImageSizeY"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnImageSizeZ (int)    : Size of an image in the Z-direction (in pixels)
		star_map["rlnImageSizeZ"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnImageWeight (RFLOAT) : Relative weight of an image
		star_map["rlnImageWeight"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnIncrementImageSize (int)    : Number of Fourier shells to be included beyond the resolution where SSNR^MAP drops below 1
		star_map["rlnIncrementImageSize"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnIs3DSampling (bool)   : Flag to indicate this concerns a 3D sampling
		star_map["rlnIs3DSampling"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnIs3DTranslationalSampling (bool)   : Flag to indicate this concerns a x,y,z-translational sampling
		star_map["rlnIs3DTranslationalSampling"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnIsFlip (bool)   : Flag to indicate that an image should be mirrored
		star_map["rlnIsFlip"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnIsHelix (bool)   : Flag to indicate that helical refinement should be performed
		star_map["rlnIsHelix"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnJoinHalvesUntilThisResolution (RFLOAT) : Resolution (in Angstrom) to join the two random half-reconstructions to prevent their diverging orientations (for C-symmetries)
		star_map["rlnJoinHalvesUntilThisResolution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnKullbackLeibnerDivergence (RFLOAT) : Kullback-Leibner divergence for a particle
		star_map["rlnKullbackLeibnerDivergence"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnKurtosisExcessValue (RFLOAT) : Kurtosis excess (4th moment - 3) for the pixel values in an image
		star_map["rlnKurtosisExcessValue"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnLensStability (RFLOAT) : Lens stability (in ppm)
		star_map["rlnLensStability"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnLocalSymmetryFile (string) : Local symmetry description file containing list of masks and their operators
		star_map["rlnLocalSymmetryFile"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnLogAmplitudesIntercept (RFLOAT) : Y-value for Guinier plot: the fitted plateau of the logarithm of the radially averaged amplitudes
		star_map["rlnLogAmplitudesIntercept"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnLogAmplitudesMTFCorrected (RFLOAT) : Y-value for Guinier plot: the logarithm of the radially averaged amplitudes after MTF correction
		star_map["rlnLogAmplitudesMTFCorrected"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnLogAmplitudesOriginal (RFLOAT) : Y-value for Guinier plot: the logarithm of the radially averaged amplitudes of the input map
		star_map["rlnLogAmplitudesOriginal"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnLogAmplitudesSharpened (RFLOAT) : Y-value for Guinier plot: the logarithm of the radially averaged amplitudes after sharpening
		star_map["rlnLogAmplitudesSharpened"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnLogAmplitudesWeighted (RFLOAT) : Y-value for Guinier plot: the logarithm of the radially averaged amplitudes after FSC-weighting
		star_map["rlnLogAmplitudesWeighted"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnLogLikeliContribution (RFLOAT) : Contribution of a particle to the log-likelihood target function
		star_map["rlnLogLikeliContribution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnLogLikelihood (RFLOAT) : Value of the log-likelihood target function
		star_map["rlnLogLikelihood"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnLongitudinalDisplacement (RFLOAT) : Longitudinal displacement (in Angstroms)
		star_map["rlnLongitudinalDisplacement"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMagnification (RFLOAT) : Magnification at the detector (in times)
		star_map["rlnMagnification"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMagnificationCorrection (RFLOAT) : Magnification correction value for an image
		star_map["rlnMagnificationCorrection"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMagnificationSearchRange (RFLOAT) : Search range for magnification correction
		star_map["rlnMagnificationSearchRange"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMagnificationSearchStep (RFLOAT) : Step size  for magnification correction
		star_map["rlnMagnificationSearchStep"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMaskName (string) : Name of an image that contains a [0,1] mask
		star_map["rlnMaskName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnMatrix_1_1 (RFLOAT) : Matrix element (1,1) of a 3x3 matrix
		star_map["rlnMatrix_1_1"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMatrix_1_2 (RFLOAT) : Matrix element (1,2) of a 3x3 matrix
		star_map["rlnMatrix_1_2"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMatrix_1_3 (RFLOAT) : Matrix element (1,3) of a 3x3 matrix
		star_map["rlnMatrix_1_3"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMatrix_2_1 (RFLOAT) : Matrix element (2,1) of a 3x3 matrix
		star_map["rlnMatrix_2_1"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMatrix_2_2 (RFLOAT) : Matrix element (2,1) of a 3x3 matrix
		star_map["rlnMatrix_2_2"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMatrix_2_3 (RFLOAT) : Matrix element (2,1) of a 3x3 matrix
		star_map["rlnMatrix_2_3"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMatrix_3_1 (RFLOAT) : Matrix element (3,1) of a 3x3 matrix
		star_map["rlnMatrix_3_1"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMatrix_3_2 (RFLOAT) : Matrix element (3,1) of a 3x3 matrix
		star_map["rlnMatrix_3_2"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMatrix_3_3 (RFLOAT) : Matrix element (3,1) of a 3x3 matrix
		star_map["rlnMatrix_3_3"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMaxNumberOfPooledParticles (int)    : Maximum number particles that are processed together to speed up calculations
		star_map["rlnMaxNumberOfPooledParticles"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnMaxValueProbDistribution (RFLOAT) : Maximum value of the (normalised) probability function for a particle
		star_map["rlnMaxValueProbDistribution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMaximumCoarseImageSize (int)    : Maximum size of the images to be used in the first pass of the adaptive oversampling strategy (may be smaller than the original image size)
		star_map["rlnMaximumCoarseImageSize"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnMaximumValue (RFLOAT) : Maximum value for the pixels in an image
		star_map["rlnMaximumValue"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMicrographId (long)   : ID (i.e. a unique number) of a micrograph
		star_map["rlnMicrographId"] = std::unique_ptr< StarOptions >(new StarOption<long>(0.0));
		// rlnMicrographMovieName (string) : Name of a micrograph movie stack
		star_map["rlnMicrographMovieName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnMicrographName (string) : Name of a micrograph
		star_map["rlnMicrographName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnMicrographNameNoDW (string) : Name of a micrograph without dose weighting
		star_map["rlnMicrographNameNoDW"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnMicrographTiltAngle (RFLOAT) : Tilt angle (in degrees) used to collect a micrograph
		star_map["rlnMicrographTiltAngle"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMicrographTiltAxisDirection (RFLOAT) : Direction of the tilt-axis (in degrees) used to collect a micrograph
		star_map["rlnMicrographTiltAxisDirection"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMicrographTiltAxisOutOfPlane (RFLOAT) : Out-of-plane angle (in degrees) of the tilt-axis used to collect a micrograph (90=in-plane)
		star_map["rlnMicrographTiltAxisOutOfPlane"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnMinRadiusNnInterpolation (int)    : Minimum radius for NN-interpolation (in Fourier pixels), for smaller radii linear int. is used
		star_map["rlnMinRadiusNnInterpolation"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnMinimumValue (RFLOAT) : Minimum value for the pixels in an image
		star_map["rlnMinimumValue"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnModelStarFile (string) : STAR file with metadata for the model that is being refined
		star_map["rlnModelStarFile"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnModelStarFile2 (string) : STAR file with metadata for the second model that is being refined (from random halves of the data)
		star_map["rlnModelStarFile2"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnMovieFrameNumber (int)    : Number of a movie frame
		star_map["rlnMovieFrameNumber"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnMovieFramesRunningAverage (int)    : Number of movie frames inside the running average that will be used for movie-refinement
		star_map["rlnMovieFramesRunningAverage"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnMtfValue (RFLOAT) : Value of the detectors modulation transfer function (between 0 and 1)
		star_map["rlnMtfValue"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnNormCorrection (RFLOAT) : Normalisation correction value for an image
		star_map["rlnNormCorrection"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnNormCorrectionAverage (RFLOAT) : Average value (over all images) of the normalisation correction values
		star_map["rlnNormCorrectionAverage"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnNrBodies (int)    : The number of independent rigid bodies to be refined in multi-body refinement
		star_map["rlnNrBodies"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnNrClasses (int)    : The number of references (i.e. classes) to be used in refinement
		star_map["rlnNrClasses"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnNrGroups (int)    : The number of different groups of images (each group has its own noise spectrum, and intensity-scale correction)
		star_map["rlnNrGroups"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnNrHelicalAsymUnits (int)    : How many new helical asymmetric units are there in each box
		star_map["rlnNrHelicalAsymUnits"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnNrOfFrames (int)    : Number of movie frames that were collected for this particle
		star_map["rlnNrOfFrames"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnNrOfSignificantSamples (int)    : Number of orientational/class assignments (for a particle) with sign.probabilities in the 1st pass of adaptive oversampling
		star_map["rlnNrOfSignificantSamples"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnNumberOfIterWithoutChangingAssignments (int)    : Number of iterations that have passed without large changes in orientation and class assignments
		star_map["rlnNumberOfIterWithoutChangingAssignments"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnNumberOfIterWithoutResolutionGain (int)    : Number of iterations that have passed without a gain in resolution
		star_map["rlnNumberOfIterWithoutResolutionGain"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnNumberOfIterations (int)    : Maximum number of iterations to be performed
		star_map["rlnNumberOfIterations"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnOffsetRange (RFLOAT) : Search range for the origin offsets (in Angstroms)
		star_map["rlnOffsetRange"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnOffsetStep (RFLOAT) : Step size for the searches in the origin offsets (in Angstroms)
		star_map["rlnOffsetStep"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnOrientSamplingStarFile (string) : STAR file with metadata for the orientational sampling
		star_map["rlnOrientSamplingStarFile"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnOrientationDistribution (RFLOAT) : Probability Density Function of the orientations  (i.e. fraction of images assigned to each orient)
		star_map["rlnOrientationDistribution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnOrientationalPriorMode (int)    : Mode for prior distributions on the orientations (0=no prior; 1=(rot,tilt,psi); 2=(rot,tilt); 3=rot; 4=tilt; 5=psi)
		star_map["rlnOrientationalPriorMode"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnOrientationsID (long)   : ID (i.e. a unique number) for an orientation
		star_map["rlnOrientationsID"] = std::unique_ptr< StarOptions >(new StarOption<long>(0.0));
		// rlnOriginX (RFLOAT) : X-coordinate (in pixels) for the origin of rotation
		star_map["rlnOriginX"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnOriginXPrior (RFLOAT) : Center of the prior on the X-coordinate (in pixels) for the origin of rotation
		star_map["rlnOriginXPrior"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnOriginY (RFLOAT) : Y-coordinate (in pixels) for the origin of rotation
		star_map["rlnOriginY"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnOriginYPrior (RFLOAT) : Center of the prior on the X-coordinate (in pixels) for the origin of rotation
		star_map["rlnOriginYPrior"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnOriginZ (RFLOAT) : Z-coordinate (in pixels) for the origin of rotation
		star_map["rlnOriginZ"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnOriginZPrior (RFLOAT) : Center of the prior on the X-coordinate (in pixels) for the origin of rotation
		star_map["rlnOriginZPrior"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnOriginalImageSize (int)    : Original size of the images (in pixels)
		star_map["rlnOriginalImageSize"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnOriginalParticleName (string) : Original name for a particles
		star_map["rlnOriginalParticleName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnOutputRootName (string) : Rootname for all output files (this may include a directory structure, which should then exist)
		star_map["rlnOutputRootName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnOverallAccuracyRotations (RFLOAT) : Overall accuracy of the rotational assignments (in degrees)
		star_map["rlnOverallAccuracyRotations"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnOverallAccuracyTranslations (RFLOAT) : Overall accuracy of the translational assignments (in pixels)
		star_map["rlnOverallAccuracyTranslations"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnOverallFourierCompleteness (RFLOAT) : Fraction of all Fourier components up to the current resolution with SNR>1
		star_map["rlnOverallFourierCompleteness"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnPaddingFactor (RFLOAT) : Oversampling factor for Fourier transforms of the references
		star_map["rlnPaddingFactor"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnParticleDiameter (RFLOAT) : Diameter of the circular mask to be applied to all experimental images (in Angstroms)
		star_map["rlnParticleDiameter"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnParticleFigureOfMerit (RFLOAT) : Developmental FOM for a particle
		star_map["rlnParticleFigureOfMerit"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnParticleId (long)   : ID (i.e. a unique number) for a particle
		star_map["rlnParticleId"] = std::unique_ptr< StarOptions >(new StarOption<long>(0.0));
		// rlnParticleName (string) : Name for a particles
		star_map["rlnParticleName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnParticleSelectZScore (RFLOAT) : Sum of Z-scores from particle_select. High Z-scores are likely to be outliers.
		star_map["rlnParticleSelectZScore"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnPerFrameCumulativeWeight (RFLOAT) : Sum of the resolution-dependent relative weights from the first frame until the given frame
		star_map["rlnPerFrameCumulativeWeight"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnPerFrameRelativeWeight (RFLOAT) : The resolution-dependent relative weights for a given frame
		star_map["rlnPerFrameRelativeWeight"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnPhaseShift (RFLOAT) : Phase-shift from a phase-plate (in degrees)
		star_map["rlnPhaseShift"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnPipeLineEdgeFromNode (string) : Name of the origin of an edge
		star_map["rlnPipeLineEdgeFromNode"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnPipeLineEdgeProcess (string) : Name of the destination of an edge
		star_map["rlnPipeLineEdgeProcess"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnPipeLineEdgeToNode (string) : Name of the to-Node in an edge
		star_map["rlnPipeLineEdgeToNode"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnPipeLineJobCounter (int)    : Number of the last job in the pipeline
		star_map["rlnPipeLineJobCounter"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnPipeLineNodeName (string) : Name of a Node in the pipeline
		star_map["rlnPipeLineNodeName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnPipeLineNodeType (int)    : Type of a Node in the pipeline
		star_map["rlnPipeLineNodeType"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnPipeLineProcessAlias (string) : Alias of a Process in the pipeline
		star_map["rlnPipeLineProcessAlias"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnPipeLineProcessName (string) : Name of a Process in the pipeline
		star_map["rlnPipeLineProcessName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnPipeLineProcessStatus (int)    : Status of a Process in the pipeline (running, scheduled, finished or cancelled)
		star_map["rlnPipeLineProcessStatus"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnPipeLineProcessType (int)    : Type of a Process in the pipeline
		star_map["rlnPipeLineProcessType"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnPixelSize (RFLOAT) : Size of the pixels in the references and images (in Angstroms)
		star_map["rlnPixelSize"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnPsiStep (RFLOAT) : Step size (in degrees) for the sampling of the in-plane rotation angle (psi)
		star_map["rlnPsiStep"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnRadiusMaskExpImages (int)    : Radius of the circular mask to be applied to all experimental images (in Angstroms)
		star_map["rlnRadiusMaskExpImages"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnRadiusMaskMap (int)    : Radius of the spherical mask to be applied to all references (in Angstroms)
		star_map["rlnRadiusMaskMap"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnRandomSeed (int)    : Seed (i.e. a number) for the random number generator
		star_map["rlnRandomSeed"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnRandomSubset (int)    : Random subset to which this particle belongs
		star_map["rlnRandomSubset"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnReconstructImageName (string) : Name of an image to be used for reconstruction only
		star_map["rlnReconstructImageName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnReferenceDimensionality (int)    : Dimensionality of the references (2D/3D)
		star_map["rlnReferenceDimensionality"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnReferenceImage (string) : Name of a reference image
		star_map["rlnReferenceImage"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnReferenceSigma2 (RFLOAT) : Spherical average of the estimated power in the noise of a reference
		star_map["rlnReferenceSigma2"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnReferenceSpectralPower (RFLOAT) : Spherical average of the power of the reference
		star_map["rlnReferenceSpectralPower"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnReferenceTau2 (RFLOAT) : Spherical average of the estimated power in the signal of a reference
		star_map["rlnReferenceTau2"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnRefsAreCtfCorrected (bool)   : Flag to indicate that the input references have been CTF-amplitude corrected
		star_map["rlnRefsAreCtfCorrected"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnResolution (RFLOAT) : Resolution (in 1/Angstroms)
		star_map["rlnResolution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnResolutionInversePixel (RFLOAT) : Resolution (in 1/pixel, Nyquist = 0.5)
		star_map["rlnResolutionInversePixel"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnResolutionSquared (RFLOAT) : X-value for Guinier plot: squared resolution in 1/Angstrom^2
		star_map["rlnResolutionSquared"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSGDGradientImage (string) : Name of image containing the SGD gradient
		star_map["rlnSGDGradientImage"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnSamplingPerturbFactor (RFLOAT) : Factor for random perturbation on the orientational sampling (between 0 no perturbation and 1 very strong perturbation)
		star_map["rlnSamplingPerturbFactor"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSamplingPerturbInstance (RFLOAT) : Random instance of the random perturbation on the orientational sampling
		star_map["rlnSamplingPerturbInstance"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSamplingRate (RFLOAT) : Sampling rate of an image (in Angstrom/pixel)
		star_map["rlnSamplingRate"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSamplingRateX (RFLOAT) : Sampling rate in X-direction of an image (in Angstrom/pixel)
		star_map["rlnSamplingRateX"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSamplingRateY (RFLOAT) : Sampling rate in Y-direction of an image (in Angstrom/pixel)
		star_map["rlnSamplingRateY"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSamplingRateZ (RFLOAT) : Sampling rate in Z-direction of an image (in Angstrom/pixel)
		star_map["rlnSamplingRateZ"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSelected (bool)   : Flag whether an entry in a metadatatable is selected in the viewer or not
		star_map["rlnSelected"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnSgdMaxSubsets (long)   : Stop SGD after doing this many subsets (possibly spanning more than 1 iteration)
		star_map["rlnSgdMaxSubsets"] = std::unique_ptr< StarOptions >(new StarOption<long>(0.0));
		// rlnSgdMuFactor (RFLOAT) : The mu-parameter that controls the momentum of the SGD gradients
		star_map["rlnSgdMuFactor"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSgdNextSubset (int)    : Number of the next subset to restart this run with
		star_map["rlnSgdNextSubset"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnSgdSigma2FudgeHalflife (long)   : After processing this many particles the multiplicative factor for the noise variance will have halved
		star_map["rlnSgdSigma2FudgeHalflife"] = std::unique_ptr< StarOptions >(new StarOption<long>(0.0));
		// rlnSgdSigma2FudgeInitial (RFLOAT) : The variance of the noise will initially be multiplied with this value (larger than 1)
		star_map["rlnSgdSigma2FudgeInitial"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSgdStepsize (RFLOAT) : Stepsize in SGD updates)
		star_map["rlnSgdStepsize"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSgdSubsetSize (long)   : The number of particles in the random subsets for SGD
		star_map["rlnSgdSubsetSize"] = std::unique_ptr< StarOptions >(new StarOption<long>(0.0));
		// rlnSgdWriteEverySubset (int)    : Every this many subsets the model is written to disk
		star_map["rlnSgdWriteEverySubset"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnSigma2Noise (RFLOAT) : Spherical average of the standard deviation in the noise (sigma)
		star_map["rlnSigma2Noise"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSigmaOffsets (RFLOAT) : Standard deviation in the origin offsets (in Angstroms)
		star_map["rlnSigmaOffsets"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSigmaPriorPsiAngle (RFLOAT) : Standard deviation of the prior on the psi (i.e. third Euler) angle
		star_map["rlnSigmaPriorPsiAngle"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSigmaPriorRotAngle (RFLOAT) : Standard deviation of the prior on the rot (i.e. first Euler) angle
		star_map["rlnSigmaPriorRotAngle"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSigmaPriorTiltAngle (RFLOAT) : Standard deviation of the prior on the tilt (i.e. second Euler) angle
		star_map["rlnSigmaPriorTiltAngle"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSignalToNoiseRatio (RFLOAT) : Spectral signal-to-noise ratio for a reference
		star_map["rlnSignalToNoiseRatio"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSkewnessValue (RFLOAT) : Skewness (3rd moment) for the pixel values in an image
		star_map["rlnSkewnessValue"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSmallestChangesClasses (int)    : Smallest changes thus far in the optimal class assignments (in numer of particles).
		star_map["rlnSmallestChangesClasses"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnSmallestChangesOffsets (RFLOAT) : Smallest changes thus far in the optimal offset assignments (in pixels).
		star_map["rlnSmallestChangesOffsets"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSmallestChangesOrientations (RFLOAT) : Smallest changes thus far in the optimal orientation assignments (in degrees).
		star_map["rlnSmallestChangesOrientations"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSolventMask2Name (string) : Name of a secondary solvent mask (e.g. to flatten density inside an icosahedral virus)
		star_map["rlnSolventMask2Name"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnSolventMaskName (string) : Name of an image that contains a (possibly soft) mask for the solvent area (values=0 for solvent, values =1 for protein)
		star_map["rlnSolventMaskName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnSortedIndex (long)   : Index of a metadata entry after sorting (first sorted index is 0).
		star_map["rlnSortedIndex"] = std::unique_ptr< StarOptions >(new StarOption<long>(0.0));
		// rlnSpectralIndex (int)    : Spectral index (i.e. distance in pixels to the origin in Fourier space)
		star_map["rlnSpectralIndex"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		// rlnSpectralOrientabilityContribution (RFLOAT) : Spectral SNR contribution to the orientability of individual particles
		star_map["rlnSpectralOrientabilityContribution"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSphericalAberration (RFLOAT) : Spherical aberration (in millimeters)
		star_map["rlnSphericalAberration"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnSsnrMap (RFLOAT) : Spectral signal-to-noise ratio as defined for MAP estimation (SSNR^MAP)
		star_map["rlnSsnrMap"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnStandardDeviationValue (RFLOAT) : Standard deviation for the pixel values in an image
		star_map["rlnStandardDeviationValue"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnStarFileMovieParticles (string) : Filename of a STAR file with movie-particles in it
		star_map["rlnStarFileMovieParticles"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnSymmetryGroup (string) : Symmetry group (e.g., C1, D7, I2, I5, etc.)
		star_map["rlnSymmetryGroup"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnTau2FudgeFactor (RFLOAT) : Regularisation parameter with which estimates for the power in the references will be multiplied (T in original paper)
		star_map["rlnTau2FudgeFactor"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnTauSpectrumName (string) : Name of a STAR file that holds a tau2-spectrum
		star_map["rlnTauSpectrumName"] = std::unique_ptr< StarOptions >(new StarOption<std::string>(""));
		// rlnTiltAngleLimit (RFLOAT) : Values to which to limit the tilt angles (positive for keeping side views, negative for keeping top views)
		star_map["rlnTiltAngleLimit"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnTransversalDisplacement (RFLOAT) : Transversal displacement (in Angstroms)
		star_map["rlnTransversalDisplacement"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnUseTooCoarseSampling (bool)   : Flag to indicate that the angular sampling on the sphere will be one step coarser than needed to speed up calculations
		star_map["rlnUseTooCoarseSampling"] = std::unique_ptr< StarOptions >(new StarOption<bool>(false));
		// rlnVoltage (RFLOAT) : Voltage of the microscope (in kV)
		star_map["rlnVoltage"] = std::unique_ptr< StarOptions >(new StarOption<float>(0.0));
		// rlnWidthMaskEdge (int)    : Width (in pixels) of the soft edge for spherical/circular masks to be used for solvent flattening
		star_map["rlnWidthMaskEdge"] = std::unique_ptr< StarOptions >(new StarOption<int>(0));
		all_image_info_.push_back(std::move(star_map));
	}
};

}  // electron_density
}  // protocols

#endif
