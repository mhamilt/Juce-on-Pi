#include "LinearJoeyMoogSC.hpp"

void JoeyMoogSC::setSampleRate(double extSampRate)
{
    sampleRate = extSampRate;
    timeStep = 1/sampleRate;
}
//==============================================================================
double JoeyMoogSC::filter(const double sample, const double resonanceSideChain/*sideChain*/, const double cutoffSideChain)
{
    //==========================================================================
    // signals coming in
    const double in = sample;
    resonance = resonanceSideChain;
    cutoff = cutoffSideChain;
    //==========================================================================
    r = resonance;
    trimRange(r, 0., .9873);
    w0 = 2*pi*20*pow(2,10*(cutoff));
    trimRange(w0, 0., .5*sampleRate);
    
    //==========================================================================
    const double k = timeStep;
    const double wk = w0*k;    // normalised frequency cutoff
    
    setIm(wk);     //  I-(kA/2)
    setIpkAx(wk);  // (I+(kA/2))x
    
    // determinant = pow(ophkw0,4) + 0.25f*r*pow(wk,4);
    const double D = 1.0/pow((1+ (wk*.5)),4) + 0.25*r*pow(wk,4); // reciprocal of the determinate of Im:
    setInvMatrix(D); // inverse of Im
    
    //==========================================================================
    // START HERE IF SIDE CHAIN IS STATIC
    //    kIpkABin[4] = { in*(wk-0.5*wk*wk), in*(0.5*wk*wk), 0, 0 }; // input here
    kIpkABin[0] = in*(wk-0.5*wk*wk);
    kIpkABin[1] = in*(0.5*wk*wk);
    
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
void JoeyMoogSC::initCoefMatrices()
{
    for (int i = 0; i < 4; i++)
    {
        std::fill(Im[i], Im[i]+4, 0);
    }
}

void JoeyMoogSC::setIm(double normFreq /*w0*k*/)
{
    //    double Im[4][4] =
    //
    //      { { diagCoef,    0,           0,          2*r*w0*k },
    //        { offDiagCoef, diagCoef,    0,          0        },
    //        { 0,           offDiagCoef, diagCoef,   0        },
    //        { 0,           0,           offDiagCoef,diagCoef } };
    
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
    
    Im[0][3] = 2*r*normFreq;
}

//==============================================================================

void JoeyMoogSC::setIpkAx(double normFreq /*w0*k*/)
{
    const double a = 1.0-0.5*normFreq;
    const double b = 0.5*normFreq;
    IpkAx[0] = x[0]*a  - x[3]*2*r*normFreq;
    IpkAx[1] = x[0]*b  + x[1]*a;
    IpkAx[2] = x[1]*b  + x[2]*a;
    IpkAx[3] = x[2]*b  + x[3]*a;
}

//==============================================================================

void JoeyMoogSC::setInvMatrix(double D)
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

void JoeyMoogSC::trimRange(double &var, const double min, const double max)
{
    if (var < min){var = min;}
    if (var > max){var = max;}
}


//==============================================================================

void JoeyMoogSC::printMatsAndVects()
{
    printf("x:      left:   right:  IpkAx:  kIpkABin:\n");
    for (int i = 0; i < 4; i++) {
        printf("%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", x[i], left[i], right[i], IpkAx[i], kIpkABin[i]);
    }
    printf("\n");
    
    printf("Im:\n");
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            printf("%.4f\t",Im[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    
    printf("INV:\n");
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            printf("%.4f\t",INV[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
