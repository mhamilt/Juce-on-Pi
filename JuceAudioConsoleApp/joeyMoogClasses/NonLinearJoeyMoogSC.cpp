#include "NonLinearJoeyMoogSC.hpp"

void JoeyNonLinearMoogSC::setSampleRate(double extSampRate)
{
    sampleRate = extSampRate;
    timeStep = 1/sampleRate;
}
//==============================================================================
double JoeyNonLinearMoogSC::filter(const double sample, const double resonanceSideChain, const double cutoffSideChain)
{
    //==========================================================================
    const double in = sample;
    resonance = resonanceSideChain;
    cutoff = cutoffSideChain;
    //==========================================================================
    r = resonance;
    trimRange(r, 0., .9873);
    w0 = 2*pi*20*pow(2,10*(cutoff));
    trimRange(w0, 0., .5*sampleRate);
    //==========================================================================
    double mu = 1.0;
    if (in!=0.0)
        mu = tanh(in)/in;
    rho = r;
    if (x[3]!=0.0)
        rho = tanh(4.0*r*x[3])/(4.0*tanh(x[3]));
    //==========================================================================
    
    const double k = timeStep;
    const double wk = w0*k;     // normalised cutoff frequency
    
    for (int i = 0; i < 4; i++)
        tanhx[i] = tanh(x[i]);
    
    //==========================================================================
    
    setIm(wk);       //  I-kA/2
    setIpkAx(wk);
    // determinant = pow(ophkw0,4) + 0.25f*r*pow(wk,4);
    const double D = 1.0/pow(1.0f+0.5f*wk,4) + 0.25f*rho*pow(wk,4); // reciprocal of the determinate of Im:
    setInvMatrix(D); // inverse of I-kA/2
    
    //==========================================================================
    const double tanhrx3 = tanh(r*x[3]);
    const double B = mu*(1.0-tanhrx3*tanhrx3)/(1.0-tanh(in)*tanhrx3);
    // k*(I+kA/2)*B*input
    const double termA = B*(wk-0.5*wk*wk);
    const double termB = 0.5*wk*wk*B;
    //==========================================================================
    // START HERE IF SIDE CHAIN IS STATIC
    kIpkABin[0] = in*termA;
    kIpkABin[1] = in*termB;
    
    for (int i = 0; i < 4; i++)
    {
        left[i] = 0;
        right[i] = 0;
        for (int j = 0; j < 4; j++)
        {
            left[i]  += IpkAx[j]*INV[i][j];      // inv(I-kA/2)*(I+kA/2)*x
            right[i] += kIpkABin[j]*INV[i][j];  // k*(I+kA/2)*B*input
        }
        x[i] = left[i] + right[i];              // x = inv(I-kA/2)*(I+kA/2)*x + inv(I-kA/2)*k*(I+kA/2)*B*input
    }
    
    return x[3]; // output
}
//==============================================================================
void JoeyNonLinearMoogSC::initCoefMatrices()
{
    for (int i = 0; i < 4; i++)
        std::fill(Im[i], Im[i]+4, 0);
}

void JoeyNonLinearMoogSC::setIm(double normFreq)
{
    //    double Im[4][4] =
    //
    //      { { ophkw0,      0,      0, 2*r*wk },
    //        {  mhkw0, ophkw0,      0,      0 },
    //        {      0,  mhkw0, ophkw0,      0 },
    //        {      0,      0,  mhkw0, ophkw0 } };
    
    //  I is 4x4 identity matrix,
    //      A is:               B is:
    // w0*[-1  0  0 -4r           [ w0
    //      1 -1  0  0               0
    //      0  1 -1  0               0
    //      0  0  1 -1 ]             0 ]
    const double diagCoef = 1.0+0.5*normFreq;
    const double offDiagCoef = -0.5*normFreq;
    
    for(int i = 0; i < 4; i++)
    {
        Im[i][i] = diagCoef;
    }
    
    for(int i = 0; i < 3; i++)
    {
        Im[i+1][i] = offDiagCoef;
    }
    
    Im[0][3] = 2*rho*normFreq;
}

//==============================================================================

void JoeyNonLinearMoogSC::setIpkAx(double normFreq /* w0*k */)
{
    const double a = 1.0-0.5*normFreq;
    const double b = 0.5*normFreq;
    IpkAx[0] = tanhx[0]*a  - tanhx[3]*2*rho*normFreq;
    IpkAx[1] = tanhx[0]*b  + tanhx[1]*a;
    IpkAx[2] = tanhx[1]*b  + tanhx[2]*a;
    IpkAx[3] = tanhx[2]*b  + tanhx[3]*a;
}

//==============================================================================

void JoeyNonLinearMoogSC::setInvMatrix(double D)
{
    INV[0][0] =  D*Im[1][1]*Im[2][2]*Im[3][3];
    INV[0][1] = -D*Im[0][3]*Im[2][1]*Im[3][2];
    INV[0][2] =  D*Im[0][3]*Im[1][1]*Im[3][2];
    INV[0][3] = -D*Im[0][3]*Im[1][1]*Im[2][2];
    
    INV[1][0] = -D*Im[1][0]*Im[2][2]*Im[3][3];
    INV[1][1] =  D*Im[0][0]*Im[2][2]*Im[3][3];
    INV[1][2] = -D*Im[0][3]*Im[1][0]*Im[3][2];
    INV[1][3] =  D*Im[0][3]*Im[1][0]*Im[2][2];
    
    INV[2][0] =  D*Im[1][0]*Im[2][1]*Im[3][3];
    INV[2][1] = -D*Im[0][0]*Im[2][1]*Im[3][3];
    INV[2][2] =  D*Im[0][0]*Im[1][1]*Im[3][3];
    INV[2][3] = -D*Im[0][3]*Im[1][0]*Im[2][1];
    
    INV[3][0] = -D*Im[1][0]*Im[2][1]*Im[3][2];
    INV[3][1] =  D*Im[0][0]*Im[2][1]*Im[3][2];
    INV[3][2] = -D*Im[0][0]*Im[1][1]*Im[3][2];
    INV[3][3] =  D*Im[0][0]*Im[1][1]*Im[2][2];
}

//==============================================================================

void JoeyNonLinearMoogSC::trimRange(double &var, const double min, const double max)
{
    if (var < min){var = min;}
    if (var > max){var = max;}
}
