//
//  LinearJoeyMoogSC.hpp
//
//  Linear Moog Voltage Control Filter Model.
//
//  Created by Joey Hook
//  
//

#ifndef LinearJoeyMoogSC_hpp
#define LinearJoeyMoogSC_hpp

#include <iostream>
#include <cmath>
/**
 MoogVCF class: Initialised with sample rate (getSampleRate() in JUCE)
 @version version number
 @author Joey Hook
 */
class JoeyMoogSC
{
public:
    //==============================================================================
    /** Constructor */
    JoeyMoogSC(){};
    JoeyMoogSC(double extSampRate)
    {
        init(extSampRate);
    };
    /** Destructor */
    ~JoeyMoogSC(){};
    //==============================================================================
    /**
     apply Moog Voltage Controlled Filter to incoming audio samples
     
     @param sample input audio sample as double
     
     @returns sample processed through MoogVCF algorithm
     */
    double filter (const double sample, const double resonanceSideChain/*sideChain*/, const double cutoffSideChain);
    //==============================================================================
    
    /**
     prints the current values of all matrices and vectors to standard character out
     */
    void printMatsAndVects();
    
    /** initialiase filter settings and coefficients
     @param extSampRate sample rate of environment
     */
    void init(double extSampRate)
    {
        setSampleRate(extSampRate);
        initCoefMatrices();
    };
    
private:
    //==============================================================================
    /**
     Sets the internal sample rate and time step values
     
     @param extSampRate sample rate from external environment (e.g. JUCE)
     */
    void setSampleRate(double extSampRate);
    //==============================================================================
    /**
     set Coefficient matrices to initial values
     */
    void initCoefMatrices();
    
    /**
     Trim variable to equal max or min if outwith this range
     
     @param var variable to be trimmed
     @param max maximum value
     @param min minimum value
     */
    void trimRange(double &var, double max, double min);
    //==============================================================================
    
    /**
     sets the inverse matrix of Im [I-(kA/2)] to the correct values: This is hard coded as the
     solution as been precomputed to save on computation
     
     @param coefD coefficient to multiply matrix with.
     */
    void setInvMatrix(double coefD);
    
    /**
     sets the Im matrix (2D array):

     @param normFreq normalised cutoff frequency
     */
    void setIm(double normFreq);
    
    /**
     propagates the IpkAx array
     @param normFreq normalised cutoff frequency
     */
    void setIpkAx(double normFreq);

private:
    //==============================================================================
    /** resonance: is in range [0,1] 
        @attention unstable when > .98729
     */
    double resonance = 0.5;
    /** cutoff: is in range [0,sampleRate/2]
     @attention Maxing out when > .745
     */
    double cutoff = .7;
    //==============================================================================
    /** internal pi constant */
    const double pi = 3.1415926536;
    /** resonance parameter */
    double r = 0.9; // resonance variable
    /** cutoff frequency: w0 is in the range 2*pi*[20,20k] */
    double w0 = 2*pi*20*pow(2,10*(0.9f));
    /** internal sampleRate */
    double sampleRate;
    /** internal time step (1/sampleRate) */
    double timeStep;
    
    /** VCF ladder output tap */
    double x[4]     = {0,0,0,0};
    /** The left term in the full output calculation: inv(I-kA/2)*(I+kA/2)*x */
    double left[4]  = {0,0,0,0};
    /** the right term in the full output calculation: inv(I-kA/2)k*(I+kA/2)*B*input */
    double right[4] = {0,0,0,0};
    /** (I plus half(k times A)) times x*/
    double IpkAx[4] = {0,0,0,0};
    /** k*(I+kA/2)*B*in */
    double kIpkABin[4] = {0,0,0,0};
    /** I-(kA/2) */
    double Im[4][4] = {0};
    /** The inverse of Im: I-(kA/2) */
    double INV[4][4] = {0};

};

#endif /* LinearJoeyMoogSC_hpp */
