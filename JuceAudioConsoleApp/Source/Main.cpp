/*
 ==============================================================================
 
 This file was auto-generated!
 
 It contains the basic startup code for a JUCE application.
 
 ==============================================================================
 */

#include "../JuceLibraryCode/JuceHeader.h"
#include "AudioProcessing.hpp"
#include <wiringPi.h>
//==============================================================================
int main (int argc, char* argv[])
{
    AudioProcessing processor;
    wiringPiSetup () ;
    pinMode (0, OUTPUT) ;
    
    while (true)
    {
        digitalWrite (0, HIGH) ; delay (500) ;
        printf("on\n");
        digitalWrite (0,  LOW) ; delay (500) ;
        printf("off\n");
    }
    return 0;
}
//==============================================================================

