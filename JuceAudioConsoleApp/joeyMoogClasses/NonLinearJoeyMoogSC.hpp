#ifndef NonLinearJoeyMoogSC_hpp
#define NonLinearJoeyMoogSC_hpp

#include <iostream>
#include <cmath>
/**
 MoogVCF class: Initialised with sample rate (getSampleRate() in JUCE)
 @version version number
 @author Joey Hook
 */
class JoeyNonLinearMoogSC
{
public:
    //==========================================================================
    /** Constructor */
    JoeyNonLinearMoogSC(){};
    JoeyNonLinearMoogSC(double extSampRate)
    {
        init(extSampRate);
    };
    /** Destructor */
    ~JoeyNonLinearMoogSC(){};
    //==========================================================================
    /**
     apply Moog VCF filter to incoming audio samples
     
     @param sample input audio sample as double
     
     @returns sample processed through MoogVCF algorithm
     */
    double filter (const double sample, const double sideChain/*sideChain*/, const double secSideChain);
    //==========================================================================
    /** initialiase filter settings and coefficients
     @param extSampRate sample rate of environment
     */
    void init(double extSampRate)
    {
        setSampleRate(extSampRate);
        initCoefMatrices();
    }
private:
    /**
     Sets the internal sample rate and time step values
     
     @param extSampRate sample rate from external environment (e.g. JUCE)
     */
    void setSampleRate(double extSampRate);
    //==========================================================================
    /**
     <#Description#>
     */
    void initCoefMatrices();
    
    /**
     <#Description#>
     @param coefD <#coefD description#>
     */
    void setInvMatrix(double coefD);
    //==========================================================================
    /**
     Trim variable to equal max or min if out
     of range
     
     @param var variable to be trimmed
     @param max maximum value
     @param min minimum value
     */
    void trimRange(double &var, double max, double min);
    
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
    //==========================================================================
public:
    
    /** <#Description#> */
    double resonance = 0;
    /** cutoff: is in range [1,0] */
    double cutoff;
    //==========================================================================
private:
    /** internal pi constant */
    const double pi = 3.1415926536;
    /** resonance parameter */
    double r = 0.9; // resonance variable
    /** <#Description#> */
    double rho;
    /** cutoff frequency: cutoff is in range [1,0]
     
     @important w0 is 2*pi*[20,20k]
     */
    double w0 = 2*pi*20*pow(2,10*(0.9f));// cutoff frequency (*cutoff is in range [1,0], w0 is 2*pi*[20,20k])
    /** internal sampleRate */
    double sampleRate;
    /** internal time step (1/sampleRate) */
    double timeStep;
    
    /** VCF ladder output tap */
    double x[4]     = {0,0,0,0};
    /** <#Description#> */
    double tanhx[4]    = {0,0,0,0};
    /** <#Description#> inv(I-kA/2)*(I+kA/2)*x */
    double left[4]  = {0,0,0,0};
    /** <#Description#> inv(I-kA/2)k*(I+kA/2)*B*input */
    double right[4] = {0,0,0,0};
    /** <#Description#> */
    double IpkAx[4] = {0,0,0,0};
    /** <#Description#> */
    double kIpkABin[4] = {0,0,0,0};
    
    /** <#Description#> */
    double Im[4][4] = {0};
    /** <#Description#> */
    double INV[4][4] = {0};
    //==========================================================================
};

#endif /* NonLinearJoeyNonLinearMoogSC_hpp */
