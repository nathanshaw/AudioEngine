#ifndef __CLIP_COUNTER_H__
#define __CLIP_COUNTER_H__

#include "AudioStream.h"

class AudioAnalyzeClipCounter : public AudioStream
{
    public:
        AudioAnalyzeClipCounter(void) : AudioStream(1, inputQueueArray){
            clips = 0;
            threshold = 32767;    // the max value possible (same for min)
            clip_debounce = 1000; // equal to one ms 
            active = true;
        }

        virtual void update(void);
        uint32_t getNumClips(void){return clips;};

        void resetCounter(void) {clips = 0;};
        void setThreshold(uint16_t t) {threshold = t;};
        void setDebounceTime(uint32_t t) {clip_debounce = t;};

        void deactivate(){active = false;};
        void activate(){active = true;};

    private:
        audio_block_t *inputQueueArray[1];
        elapsedMicros last_clip;
        uint32_t clips;
        uint16_t threshold;
        uint32_t clip_debounce; // 1000 equals one ms

        bool active = false; // used to turn off the analyser for when the bot is actuating
};

void AudioAnalyzeClipCounter::update() {
    audio_block_t *block  = receiveReadOnly(0); // get the read-only memory block
    if (!block || !active) return;

    // go through each sample to see if a clip is registered
    for (int i = 0; i < AUDIO_BLOCK_SAMPLES; i++) {
        if (last_clip > clip_debounce) {
            if (block->data[i] >= threshold) {
                clips++;
                last_clip = 0;
            } else if (block->data[i] <= threshold * -1) {
                clips++;
                last_clip = 0;
            }
        }
    }
    release(block);// releases the block from memory to prevent leaks
}

#endif // __CLIP_COUNTER_H__
