#ifndef __ZERO_CROSSING_H__
#define __ZERO_CROSSING_H__

#include "AudioStream.h"

class AudioAnalyzeZeroCrossing: public AudioStream
{
    public:
        AudioAnalyzeZeroCrossing(void) : AudioStream(1, inputQueueArray){
            crossings = 0;
        }

        virtual void update(void);
        uint32_t getNumCrossings(void){return crossings;};
        void resetCounter(void) {crossings = 0;};

    private:
        audio_block_t *inputQueueArray[1];
        uint32_t crossings;
};

void AudioAnalyzeZeroCrossing::update() {
    audio_block_t *block  = receiveReadOnly(0); // get the read-only memory block
    if (!block || !active) return;
    // start at second sample so we can always back-check
    for (int i = 1; i < AUDIO_BLOCK_SAMPLES; i++) {
        if (block->data[i] >= 0 && block->data[i-1] < 0) {
            crossings++;
        } else if (block->data[i] <= 0 && block->data[i-1] > 0) {
            crossings++;
        }
    }
    release(block);// releases the block from memory to prevent leaks
}
#endif // __ZERO_CROSSING_H__
