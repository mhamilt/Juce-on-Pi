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
    deviceManager.initialise(1, 1, nullptr, true);
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
    for (int channel = 0; channel < numOutputChannels; channel++)
    {
        for (int i = 0; i < numSamples; i++)
        {
            outputChannelData[0][i] = random.nextFloat()* 0.5 - 0.25;
        }
    }
}
//==============================================================================
void AudioProcessing::audioDeviceAboutToStart(AudioIODevice* device)
{
}
//==============================================================================
void AudioProcessing::audioDeviceStopped()
{
}
 
