//
//  AudioProcessing.cpp
//  TestConsoleApp - ConsoleApp
//
//  Created by mhamilt7 on 29/10/2018.
//

#include "AudioProcessing.hpp"

//==============================================================================
AudioProcessing::AudioProcessing()
{
    deviceManager.initialise(0,2, nullptr, true);
    deviceManager.addAudioCallback(this);
    
}
AudioProcessing::~AudioProcessing()
{
    deviceManager.removeAudioCallback(this);
}

//==============================================================================
/// This is where your sample by sample processing goes

void AudioProcessing::audioDeviceIOCallback(const float** inputChannelData,
                                            int numInputChannels,
                                            float** outputChannelData,
                                            int numOutputChannels,
                                            int numSamples)
{
    
    for (int i = 0; i < numSamples; i++)
    {
        const float modSig = (sin(currRad) + 1) * .5f;
        outputChannelData[0][i] = moogNonLin.filter(random.nextFloat()* 0.5 - 0.25,
                                                    modSig,
                                                    .4);
        currRad += radsPerSec;
    }
    
    if (numOutputChannels > 1)
    {
        for (int channel = 1; channel < numOutputChannels; ++channel)
            outputChannelData[channel] = outputChannelData[0];
    }
}
//==============================================================================
void AudioProcessing::audioDeviceAboutToStart(AudioIODevice* device)
{
    double sampleRate = device->getCurrentSampleRate();
    moogNonLin.init(sampleRate);
    radsPerSec = 2 * float_Pi * modFreq / sampleRate;
}
//==============================================================================
void AudioProcessing::audioDeviceStopped()
{
}

