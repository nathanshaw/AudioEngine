#ifndef __FEATURE_COLLECTOR_H__
#define __FEATURE_COLLECTOR_H__

#include <Audio.h>
#include "../NeopixelManager/NeopixelManager.h"
#include "../ValueTracker/ValueTrackerDouble.h"


#ifndef PEAK_LOG_RESET_MIN
#define PEAK_LOG_RESET_MIN 2000
#endif

class FeatureCollector {
    public:
        FeatureCollector(uint8_t channels, String _id);
        //////////////////////////////////////////////////
        /////////////////// Configuration ////////////////
        //////////////////////////////////////////////////
        /////////////// General //////////////////////////
        String name = "";
        bool microphone_active[2] = {true, true};
        elapsedMillis last_update_timer;

        /////////////////// Microphone Testing /////////////////
        bool testMicrophone();

        /////////////////// Printing ///////////////////////////
        void printFeatures();

        //////////////////////////////////////////////////
        /////////////////// Runtime //////////////////////
        //////////////////////////////////////////////////
        //
        ////////////////// Getting Functions ///////////////////
        String getName() {return name;};
        bool isActive();
        uint8_t getNumChannels(){Serial.print(F("num_channels: " ));Serial.println(num_channels);return num_channels;};

        //////////////// Gain Tracking /////////////////////////
        void   setGain(double g, int channel);
        double getGain(int channel);
        void   linkAmplifier(AudioAmplifier * amp, double low, double high, double _max_gain_adj);

        //////////////////////////////////////////////
        //////////////// RMS /////////////////////////
        void   linkRMS(AudioAnalyzeRMS *r, bool print);
        bool   getRMSActive(){return rms_active;};

        double getRMS(int channel);
        double getDominateRMS(){return rms_val[dominate_channel];}; // get the RMS reading from the dominate microphone

        double getRMSPosDelta(int channel){return rms_tracker[channel].getPosDelta();};
        double getDominatePosDelta(){return getRMSPosDelta(dominate_channel);};

        double getRMSAvg(int channel){return rms_tracker[channel].getAvg();};
        double getDominateRMSAvg(){return getRMSAvg(dominate_channel);};
        bool   resetRMSAvg(int c);
        void   resetDominateRMSAvg(){resetRMSAvg(dominate_channel);};

        void   printRMSVals();

        ///////////////////////////////////////////////
        //////////////// Peak /////////////////////////
        void   linkPeak(AudioAnalyzePeak *r, bool print);
        bool   getPeakActive(){ return peak_active;};

        double getPeak(int channel);
        double getDominatePeak(){return peak_val[dominate_channel];};

        ////////////////// Pos Delta
        double getPeakPosDelta(int channel) {return peak_tracker[channel].getPosDelta();};
        double getDominatePeakPosDelta(){return peak_tracker[dominate_channel].getPosDelta();};

        ////////////////// Avg
        double getPeakAvg(int channel, bool print);

        double getDominatePeakAvg(){return getPeakAvg(dominate_channel, false);};
        void   resetPeakAvg(int channel);
        void   resetDominatePeakAvg(){return resetPeakAvg(dominate_channel);};

        ////////////////// Clipping
        uint32_t getClips(int channel){return clips[channel];};
        void     resetClips(int channel){clips[channel] = 0;};

        ////////////////// Min / Max
        double getPeakMin(int channel){return peak_tracker[channel].getMin();};
        double getPeakMax(int channel){return peak_tracker[channel].getMax();};
        double getDominatePeakMin(){return peak_tracker[dominate_channel].getMin();};
        double getDominatePeakMax(){return peak_tracker[dominate_channel].getMax();};

        void   printPeakVals();

        //////////////// General ///////////////////////
        bool update(FFTManager1024 _fft[]);
        void setDominateChannel(int c){dominate_channel = c;};
        uint8_t calculateDominateChannel(FFTManager1024 _fft[]);
        uint8_t getDominateChannel() {return dominate_channel;};

        int getNumRMSAnas() {return num_rms_anas;};
        int getNumPeakAnas() {return num_peak_anas;};

        // change the update rate for the FC
        void setUpdateRate(uint16_t t) {fc_update_rate = t;};

        //////////////// Printing ///////////////////////
        void setPrintRMS(bool p) {print_rms = p;};
        void setPrintRMSDebug(bool p) {print_rms_debug = p;};
        void setPrintPeak(bool p) {print_peak = p;};
        void setPrintDominateChannel(bool p) {print_dom = p;};
        void setPrintAutogain(bool p) {print_autogain = p;};
        // void setPrintPeakDebug(bool p) {print_peak_debug = p;};
        //
        ////////////////// Autogain //////////////////////
        bool calculateAutogain(FFTManager1024 _fft[]);
        bool calculateAutogain(int channel, FFTManager1024 _fft[]);
        void activateAutogain(uint32_t update_rate, uint32_t initial_delay);
        double calculateAutogainCost(double val, int idx);

        void autogainTrackAvgRMS(double target, double tolerance, double _min, double _max);
        void autogainTrackOnRatio(double target, double tolerance, double _min, double _max, NeoGroup *_n);
        void autogainTrackAvgPeak(double target, double tolerance, double _min, double _max);
        void autogainTrackClipping(double target, double tolerance, double _min, double _max);

        bool updateAutogain(FFTManager1024 _fft[]);

    private:
        /////////// Autogain Active //////////////////////
        bool autogain_active = false; 
        double ag_target_vals[4];
        double ag_low_thresh[4];
        double ag_high_thresh[4];
        double ag_min_thresh[4];
        double ag_max_thresh[4];

        elapsedMillis last_autogain;

        uint32_t autogain_initial_delay = 60000 * 5;
        uint32_t autogain_update_rate   = 60000 * 3;

        bool gain_good[2] =                    {false, false};

        int ag_rms_idx;
        bool ag_tracking_rms =                 false;

        int ag_on_ratio_idx;
        bool ag_tracking_on_ratio =            false;

        int ag_peak_idx;
        bool ag_tracking_peak  =               false;

        int ag_clipping_idx;
        bool ag_tracking_clipping =            false;

        uint8_t num_tracked_autogain_features =      0;

        ////////////////// General ///////////////////////
        // this is the channel which will be used when the "dominate"
        // feature is desired
        int dominate_channel    = 0;
        uint8_t num_channels    = 0;
        // how often will the feature collecter update?, if 
        // the feature collector is updated more frequently than 
        // this number is ms, the update will be skipped
        uint16_t fc_update_rate = 10;

        ////////////////// Printing //////////////////////
        bool print_rms          = false;
        bool print_rms_debug    = false;
        bool print_peak         = false;
        bool print_peak_debug   = false;
        bool print_z_crossings  = false;
        bool print_rms_delta    = false;
        bool print_peak_delta   = false;
        bool print_dom          = false;
        bool print_autogain     = false;

        /////////////// Testing Functions ////////////////
        bool testMicrophoneRMS() {return testMicrophoneRMS(5000);};
        bool testMicrophoneRMS(int _dur);
        bool testMicrophonePeak() {return testMicrophonePeak(5000);};
        bool testMicrophonePeak(int _dur);

        //////////////// Gain Tracking ///////////////
        AudioAmplifier *amp_ana[4];
        double gain[4];
        // ValueTrackerDouble gain_tracker[2] = {ValueTrackerDouble("front_gain",&gain[0], 0.5), 
        //        ValueTrackerDouble("rear_gain", &gain[1], 0.5)};
        double min_gain[4];
        double max_gain[4];
        // 25% of the current gain
        double max_gain_adj[4];
        uint8_t audio_amp_add_idx = 0; // how many current gain objects are linked?

        //////////////// RMS /////////////////////////
        AudioAnalyzeRMS *rms_ana[2];
        uint8_t num_rms_anas = 0;
        bool rms_active = false;

        void calculateRMS(int channel);

        double rms_val[2];
        ValueTrackerDouble rms_tracker[2] = {ValueTrackerDouble("front_rms", &rms_val[0], 0.5), ValueTrackerDouble("rear_rms", &rms_val[1], 0.5)};

        // double rms_totals[2];
        // unsigned long rms_readings[2];

        elapsedMillis last_rms_reset[2];
        uint16_t rms_log_reset_min = 2000;

        //////////////// Peak /////////////////////////
        AudioAnalyzePeak *peak_ana[2];
        uint8_t num_peak_anas = 0;
        bool peak_active = false;

        double peak_val[2];
        ValueTrackerDouble peak_tracker[2] = {ValueTrackerDouble("front_peak", &peak_val[0], 0.5), ValueTrackerDouble("rear_peak", &peak_val[1], 0.5)};

        void calculatePeak(int channel);

        //////////////// Clipping /////////////////////////
        // how many times the peak value is >= 1.0
        uint32_t clips[2] = {0, 0};
        // ValueTrackerDouble clip_tracker[2] = {ValueTrackerDouble(&clips[0], 0.5), ValueTrackerDouble(&clips[1], 0.5)};

        elapsedMillis last_peak_reset[2];
        uint16_t peak_log_reset_min = 2000;

        //////////////// Z-crossings //////////////////////
        // TODO
        // uint32_t z_crossings = {0, 0};
        // ValueTrackerDouble z_crossing_tracker[2] = {ValueTrackerDouble(&z_crossings[0], 0.5), ValueTrackerDouble(&z_crossings[1], 0.5)};
        NeoGroup *neopixel_manager;
};

FeatureCollector::FeatureCollector(uint8_t channels, String _id) {
    name = _id;
    num_channels = channels;
    Serial.print(F("Initalising feature collector with "));
    Serial.print(num_channels);
    Serial.print(F(" channels of audio and an id of :"));
    Serial.println(name);
}

void FeatureCollector::linkPeak(AudioAnalyzePeak *r, bool print) {
    print_peak = print;
    peak_ana[num_peak_anas] = r;
    num_peak_anas++;
    peak_active = true;
};

void FeatureCollector::linkRMS(AudioAnalyzeRMS *r, bool print) {
    print_rms = print;
    rms_ana[num_rms_anas] = r;
    num_rms_anas++;
    rms_active = true;
};

double FeatureCollector::getPeakAvg(int channel, bool print) {
    double avg;
    avg = peak_tracker[channel].getAvg();
    if (print){
       Serial.print(F("getPeakAvg() is returning a peak average of :"));
       Serial.println(avg);
    }
    return avg;
}

bool FeatureCollector::updateAutogain(FFTManager1024 _fft[]) {
    // first thing is to see if it is time to actually update the gain or to
    // exit the function. This logic is different depending of if the
    // gain is expected to be good or bad.

    // the gain will be false if the auto-gain has not been run before
    // or if the auto_gain has not been effecting and the tracked features
    // are outside of the established values

    // no matter what if less time than the initial_update_rate has
    // passed since the last update no update will occur, so the
    // function exits
    if (last_autogain < autogain_initial_delay) {
        return false;
    }

    // if we made it this far then we are updating the gain
    // this is done here as it is more efficient due to the various
    // exiting conditions for the function
    dprintMinorDivide(print_autogain);
    dprint(print_autogain, F("AUTO GAIN INITIALISED"));

    for (int channel = 0; channel < 2; channel++) {
        // if the gain is within the accepted bounds, and less time than
        // the update_rate has passed, the function will exit
        if (gain_good[channel] == true && last_autogain < autogain_update_rate) {
            return false;
        } if (gain_good[channel] == false && last_autogain < autogain_update_rate * 0.2) {
            return false;// if the gain is not good, we update the autogain five times faster than usual
        }
        dprint(print_autogain, F("starting auto_gain for channel: "));
        dprintln(print_autogain, channel);
        calculateAutogain(channel, _fft);
    }
    last_autogain = 0;
    return true;
}

void FeatureCollector::activateAutogain(uint32_t update_rate, uint32_t initial_delay){
    autogain_active = true;
    autogain_update_rate = update_rate;
    autogain_initial_delay = initial_delay;
    last_autogain = 0;
}

bool FeatureCollector::calculateAutogain(FFTManager1024 _fft[]) {
    bool good = false;
    for (int i = 0; i < num_channels; i++) {
        if (calculateAutogain(i, _fft)){
            good = true;
        }
    }
    return good;
}

double FeatureCollector::calculateAutogainCost(double val, int idx) {
    /* Thee value is the new value for the feature and idx is which index
     * that feature has (according to when the features are added.
     *
     * This function will compare the value to the low and high thresholds
     * if the value does not exceed these values it will return a cost of 
     * 0.0
     *
     * If it is found that the value is outside of the low and high
     * thresholds the function will calculate a cost based on the appropiate
     * low or high threshold and the min and max values for that feature.
     * If it is found that the feature exceeds the min or max value a cost
     * of -1.0 or 1.0 will be returned accordingly. If the value is between
     * the appropiate low/high and min/max the cost will be calculated linearly
     * according to where the values lies within this range and a value of between
     * -1.0 and 1.0 will be returned accordingly.
    */
    double cost = 0.0;
    // if the value is within the low and high values return a 0.0 as no gain
    // adjustment is needed
    if (val > ag_low_thresh[idx] && val < ag_high_thresh[idx]) {
        return 0.0;
    }

    // if the value is below the min threshold. this results in the max cost of 1.0
    if (val < ag_min_thresh[idx]) {
        dprint(print_autogain, F("Increasing autogain by the max amount "));
        return 1.0;
    }

    // if the value is above the max threshold. this results in the min cost of -1.0
    if (val > ag_max_thresh[idx]) {
        dprint(print_autogain, F("Decreasing autogain by the max amount ")); 
        return -1.0;
    }

    // if the value is between the low and min or high and max then return a scaled
    // cost according to how to lies within that range
    if (val < ag_low_thresh[idx]) {
        // if thee gain needs to be raised
        cost = (val - ag_min_thresh[idx]) / (ag_low_thresh[idx] - ag_min_thresh[idx]);
        dprint(print_autogain, F("returning a cost of : "));
        dprint(print_autogain, cost);
        return cost;
    }

    if (val > ag_high_thresh[idx]) {
        cost = (1.0 - ((val - ag_high_thresh[idx]) / (ag_max_thresh[idx]- ag_high_thresh[idx]))) * -1;
        dprint(print_autogain, F("returning a cost of : "));
        dprint(print_autogain, cost);
        return cost;
    }
    // the program should never reach this point...
    Serial.println(F("WARNING - there is something wrong with the calculate feature function, it should be returning something before the end of the function."));
    return 0.0;
}

bool FeatureCollector::calculateAutogain(int channel, FFTManager1024 _fft[]) {
    // now that we have established that enough time has passed for the
    // auto_gain to act, the auto_gain must act. The first order of business
    // is to update the tracked features.
    double cost = 0.0;
    double new_gain = gain[channel];
    dprintMinorDivide(print_autogain);
    ////////////////////////// RMS //////////////////////////////////
    if (ag_tracking_rms == true) {
        // IMPORTANT - for the RMS feature, if the value is half
        // or less than the expected value then the gain needs to be
        // halfed or doubled accordingly.
        // dprint(print_autogain, F("rms adjusted the cost to: "));
        double _avg = getRMSAvg(channel);
        cost += calculateAutogainCost(_avg, ag_rms_idx);
        // need to reset the RMS average calculations for next time
        // dprintln(print_autogain, cost);
        // we aree halfing or doubling the signal
        // so this would equate to a max adjustment of half the signal when 
        // decreasing or double the signal when increasing
        dprintln(print_autogain, F("Now resetting the RMS avg for the channel"));
        resetRMSAvg(channel);
        dprintln(print_autogain, F("RMS Average has been reset"));
    }
    ////////////////////// On Ratio ////////////////////////////////
    if (ag_tracking_on_ratio == true) {
        dprint(print_autogain, F("On Ratio is : "));
        cost += calculateAutogainCost(neopixel_manager->getOnRatio(), ag_on_ratio_idx);
        dprintln(print_autogain, neopixel_manager->getOnRatio());
        dprint(print_autogain, F("on_ratio adjusted the cost to: "));
        dprintln(print_autogain, cost);
        neopixel_manager->resetOnRatio();
    }
    ////////////////////// On Ratio ////////////////////////////////
    if (ag_tracking_peak == true) {
        dprint(print_autogain, F("Peak Average is                   : "));
        cost += calculateAutogainCost(getPeakAvg(channel, false), ag_peak_idx);
        dprintln(print_autogain, getPeakAvg(channel, false));
        dprint(print_autogain, F("Peak average adjusted the cost to : "));
        dprintln(print_autogain, cost);
        resetPeakAvg(channel);
    }
    dprint(print_autogain, F("current gain is : "));
    dprintln(print_autogain, gain[channel]);

    dprint(print_autogain, F("Applying the cost and calculating the new gain, the adjustment is "));
    if (cost > 0.0) {
        cost = (gain[channel] * cost * max_gain_adj[channel]);
    } else if (cost < 0.0) {
        cost = (gain[channel] * cost * max_gain_adj[channel]);
    }
    dprintln(print_autogain, cost);

    dprint(print_autogain, F("ensuring gain is within min/max "));
    dprint(print_autogain, new_gain);
    dprint(print_autogain, F("\t"));
    dprint(print_autogain, F(" min/max "));
    dprint(print_autogain, min_gain[channel]);
    dprint(print_autogain, F("/"));
    dprintln(print_autogain, max_gain[channel]);
    dprint(print_autogain, F("current gain is : "));
    dprintln(print_autogain, gain[channel]);

    new_gain = gain[channel] + cost;
    dprint(print_autogain, F("new gain before clipping to min/max is : "));
    dprintln(print_autogain, new_gain);

    new_gain = maxf(minf(max_gain[channel], new_gain), min_gain[channel]);
    dprintln(print_autogain, new_gain);

    // only update the gain if it is different from the current gain
    if (new_gain != gain[channel] && new_gain > 0.0) {
        // this cost value tells us if gain adjustment is needed, if the value
        // is below 0.0 then the gain needs to be lowered, if the value is above
        // 1.0 then the gain needs to be raised.
        dprint(print_autogain, F("updating gain from "));dprint(print_autogain, gain[channel]);dprint(print_autogain, " to ");
        dprint(print_autogain, new_gain);
        gain[channel] = new_gain;
        dprint(print_autogain, F(" with a cost of "));dprintln(print_autogain, cost);// dprint(print_autogain, " auto_gain timer");dprint(print_autogain, last_autogain);
        // dprintln(print_autogain, " / ");dprintln(print_autogain, update_rate);
        // if there is a cost, then the gain needs to be adjusted
        setGain(new_gain, channel);
        // if the gain needed to be adjusted then it is not good
        gain_good[channel] = false;
        return true;
    } else {
        // if the new gain is the same as the old gain then our
        // gain is good and we dont have to check as often
        gain_good[channel] = true;
    }
    // if we got to this point and have not updated the gain we will
    return false;
}

uint8_t FeatureCollector::calculateDominateChannel(FFTManager1024 _fft[]){
    // first ensure that the gain is good for both channels before determining a dominate
    // second compare the RMS averages from each channel
    // then compare the Peak averages from each channel
    // subtract the RMS from the Peak to get the dynamic range metric
    //
    // for the second metric, noise calculate the average centroid for each channel
    // assume that the channel with the lower centroid has less noise
    //
    // for the third metric, look at the channel which is clipping less
    //
    // TODO -- feature collector needs to be tracking clips
    dprintMinorDivide(print_dom);
    dprintln(print_dom, F("Entering calculateDominateChannel()"));

    /////////////////////////// Clipping //////////////////////////////////////
    dprintMinorDivide(print_dom);
    dprint(print_dom, F("Starting to calculate the clipping feature for "));
    dprint(print_dom, num_channels);
    dprintln(print_dom, F(" channels"));
    uint32_t clips = 0;
    int clip_chan = -1;
    for (int i = 0; i < num_channels; i++) {
        dprint(print_dom, F("channel "));
        dprint(print_dom, i);
        dprint(print_dom, F(" has clipped "));
        dprint(print_dom, getClips(i));
        dprintln(print_dom, F(" times since the last reset"));
        if (i == 0) {
            clips = getClips(i);
            clip_chan = 0;
        } else if (clips > getClips(i)){
            clip_chan = i;
        } else if (clips == getClips(i)){
            clip_chan = -1;
        }
        resetClips(i);
    }

    if (clip_chan == -1) {
        dprintln(print_dom, F("Both microphones recorded the same number of clipped waveforms since the last reset"));
    } else {
        dprint(print_dom, F("channel "));
        dprint(print_dom, clip_chan);
        dprintln(print_dom, F(" clipped less than other channels since last feature reset"));
    }

    ///////////////////////////// Dynamic Range ///////////////////////////////////////
    dprintMinorDivide(print_dom);
    // dprintln(print_dom, F("Now calculating the Dynamic range of each channel"));
    uint8_t dr_chan = 0;
    double max_diff = 0.0;
    for (int i = 0; i < num_channels; i++) {
        dprint(print_dom, F("channel "));
        dprint(print_dom, i);
        dprint(print_dom, F(" has a Peak avg of "));
        double _peak = getPeakAvg(i, false);
        dprint(print_dom, _peak, 6);
        resetPeakAvg(i);
        dprint(print_dom, F(" and a RMS avg of "));
        double _rms = getRMSAvg(i);
        dprint(print_dom, _rms, 6);
        resetRMSAvg(i);

        dprint(print_dom, F(" The channel has a diff of : "));
        dprintln(print_dom, _peak - _rms, 6);
        if ((_peak - _rms) > max_diff) {
            max_diff = _peak - _rms;
            dr_chan = i;
        }
    }
    if (dr_chan == 255) {
        dprintln(print_dom, F("No channel holds an advantage for the dynamic range"));

    } else {
        dprint(print_dom, F("channel "));
        dprint(print_dom, dr_chan);
        dprintln(print_dom, F(" has the highest difference between the rms avg and the peak avg"));
    }

    ///////////////////////////// Spectral Centroid ///////////////////////////////////////
    dprintMinorDivide(print_dom);
    // dprintln(print_dom, F("Now calculating the Avg Centroid of each channel"));
    uint8_t sc_chan = -1;
    double _cent = 100000;
    for (int i = 0; i < 2; i++) {
        // dprintln(print_dom, "WARNING - THE FFTManager does not currently switch channels");
        dprint(print_dom, F("channel "));
        dprint(print_dom, i);
        dprint(print_dom, F(" has a centroid of "));
        dprintln(print_dom, _fft[i].getAvgCentroid());
        if (_fft[i].getAvgCentroid() < _cent) {
            _cent = _fft[i].getAvgCentroid();
            sc_chan = i;
        }
        _fft[i].resetAvgCentroid();
    }
    dprint(print_dom, F("channel "));
    dprint(print_dom, sc_chan);
    dprintln(print_dom, F("has the lowest centroid and is assumed to be less noisy"));

    ///////////////////////////////////////////////////////////////////////////
    uint8_t dom_channel = 0;
    uint8_t  total = 0;

    if (clip_chan >= 0) {
        total += clip_chan;
    }
    if (dr_chan < 255) {
        total += dr_chan;
    }
    if (sc_chan >= 0) {
        total  += sc_chan;
    }
    if (total  > 1) {
        dom_channel = 1;
    } else {
        // need different logic as 255 equates to no advantage detected
        dom_channel = 0;
    }
    dprint(print_dom, "channel ");
    dprint(print_dom, dom_channel);
    dprint(print_dom, F(" is the dominate channel"));
    printMinorDivide();
    dominate_channel = dom_channel;
    return dom_channel;
}

bool FeatureCollector::isActive() {
 return microphone_active[0] || microphone_active[1];
}

void FeatureCollector::linkAmplifier(AudioAmplifier * amp, double low, double high, double _max_gain_adj) { 
    max_gain[audio_amp_add_idx] = high;
    min_gain[audio_amp_add_idx] = low;
    max_gain_adj[audio_amp_add_idx] = _max_gain_adj;
    if (audio_amp_add_idx < 4) {
        Serial.print(F("Linked an audio amplifier "));Serial.print(audio_amp_add_idx);printTab();
        amp_ana[audio_amp_add_idx] = amp;
        audio_amp_add_idx = audio_amp_add_idx + 1;
        Serial.println(audio_amp_add_idx);
    }
    else {
        Serial.println(F("ERROR, can't link audio amplifier, there are not enough available slots"));
    }
}


void FeatureCollector::setGain(double g, int channel) {
    Serial.print(F("setGain called for channel "));
    Serial.print(channel);
    Serial.print(F(" and a gain of "));
    Serial.println(g, 6);
    gain[channel] = g;
    Serial.print(F("currently the min and max gain for this channel are: "));
    Serial.print(min_gain[channel]);
    Serial.print("\t");
    Serial.println(max_gain[channel]);

    if (gain[channel] > max_gain[channel]) {
        Serial.println(F("gain is greater than max gain, changing the gain to max gain"));
        gain[channel] = max_gain[channel];
    }
    if (gain[channel] < min_gain[channel]) {
        Serial.println(F("gain is less than min gain, changing the gain to min gain"));
        gain[channel] = min_gain[channel];
    }
    // the add_idx is the number of objects
    if (audio_amp_add_idx > channel) {
        Serial.println(F("gain and channel is within range, now sending gain to the gain objects"));
        amp_ana[channel]->gain(g);
    } else {
        Serial.println(F("WARNING THE CHANNEL IS OUTOF RANGE, EXITNG"));
    }
    Serial.print(F("the current gain for channel "));
    Serial.print(channel);
    Serial.print(F(" is "));
    Serial.println(gain[channel]);
}

double FeatureCollector::getGain(int channel) {
    return gain[channel];
}

bool FeatureCollector::testMicrophoneRMS(int _dur) {
    // go through and gather 10 features from each channel and make sure it is picking up audio
    for (int i = 0; i < num_rms_anas; i++) {
        uint8_t readings = 0;
        double values = 0.0;
        unsigned long a_time = millis();
        Serial.print(F("Testing "));Serial.print(i);
        Serial.println(F(" Microphone using RMS"));
        while (readings < 10 
                && millis() < a_time + _dur) {
            if (rms_ana[i]->available()) {
                values += rms_ana[i]->read();
                readings++;
                Serial.print(". ");
                delay(20);
            }
        }
        if (values > 0) {
            Serial.println();
            Serial.print(name);
            Serial.println(F(" Microphone is good"));
            microphone_active[i] = true;
            return true;
        } 
        Serial.println(F("\nERROR, "));
        Serial.print(name);Serial.println(F(" Microphone does not work"));
        printDivideLn();
        microphone_active[i] = false;
    }
    return false;
}

bool FeatureCollector::testMicrophonePeak(int _dur) {
    // go through and gather 10 features from each channel and make sure it is picking up audio
    for (int i = 0; i < num_peak_anas; i++) {
        uint8_t readings = 0;
        double values = 0.0;
        unsigned long a_time = millis();
        printDivide();
        Serial.print(F("Testing "));Serial.print(name);Serial.println(F(" Microphone using Peak"));
        while (readings < 10 && millis() < a_time + _dur) {
            if (peak_ana[i]->available()) {
                values += peak_ana[i]->read();
                readings++;
                Serial.print(". ");
                delay(20);
            }
        }
        if (values > 0) {
            Serial.println();
            Serial.print(name);
            Serial.println(F(" Microphone is good"));
            microphone_active[i] = true;
            return true;
        } 
        Serial.println(F("\nERROR, "));
        Serial.print(name);Serial.println(F(" Microphone does not work"));
        printMinorDivideLn();
        microphone_active[i] = false;
}
    return false;
}

bool FeatureCollector::testMicrophone () {
    bool done = false;
    if (rms_active) {
        testMicrophoneRMS();
        done = true;
    }
    if (peak_active) {
        testMicrophonePeak();
        done = true;
    } 
    if (done == true) {
        return true;
    }
    else {
        Serial.println(F("Sorry unable to test microphone as neither the RMS or Peak feature is active"));
        return false;
    }
}


//////////////// Update Functions ///////////////////////////////
/////////////////////////////////////////////////////////////////
void FeatureCollector::calculatePeak(int channel) {
    bool avail = peak_ana[channel]->available();
    if (peak_active && avail) {
        peak_val[channel] = peak_ana[channel]->read();
        peak_tracker[channel].update();
        if (print_peak_debug){
            peak_tracker[channel].print();
        }
        dprint(print_peak_debug, peak_tracker[channel].getName());
        dprint(print_peak_debug, F(" Peaks (normal, pos_delta):\t"));
        dprint(print_peak_debug, peak_val[channel]);
        dprint(print_peak_debug, "\t");
        dprintln(print_peak_debug, peak_tracker[channel].getPosDelta());
    }
    if (peak_val[channel] >= 1.0) {
        dprintln(print_peak_debug, F("Clip identified"));
        clips[channel]++;
    }
}

void FeatureCollector::resetPeakAvg(int channel) {
    if (last_peak_reset[channel] > PEAK_LOG_RESET_MIN) {
        peak_tracker[channel].getAvg(true);// this resets the peak value
        last_peak_reset[channel] = 0;
    }
}

void FeatureCollector::calculateRMS(int channel) {
    if (rms_active  && (rms_ana[channel]->available())) {
            double _rms = rms_ana[channel]->read();
        if (_rms > 0.0) {
            rms_val[channel] = _rms;
            rms_tracker[channel].update();
        } else {
            dprint(print_rms_debug, F("WARNING RMS is equal to 0 for channel: "));
            dprint(print_rms_debug, num_rms_anas);
        }
    }
}

bool FeatureCollector::resetRMSAvg(int c) {
    if (c < 2) {
        if (last_rms_reset[c] > rms_log_reset_min) {
            rms_tracker[c].getAvg(true);
            last_rms_reset[c] = 0;
            return true;
        }
    }
    return false;
}

///////////////////// Getter functions ///////////////////////////////
double FeatureCollector::getRMS(int channel) {
    if (rms_active) {
        return rms_val[channel];
    }
    Serial.println(F("ERROR  - RMS IS NOT AN ACTIVE AUDIO FEATURE : ")); Serial.println(name);
    return -1.0;
}

double FeatureCollector::getPeak(int channel) {
    if (peak_active) {
        return peak_val[channel];
    }
    Serial.println(F("ERROR  - Peak IS NOT AN ACTIVE AUDIO FEATURE : ")); Serial.println(name);
    return -1.0;
}

//////////////////////////////// Print Functions /////////////////////////////////////////
void FeatureCollector::printFeatures() {
    if (rms_active && print_rms) {
        printRMSVals();
    };
    if (peak_active && print_peak) {
        printPeakVals();
    };
}

void FeatureCollector::printRMSVals() {
    if (rms_active > 0) {
        for (int channel = 0; channel < num_rms_anas; channel++){
            Serial.print(name); Serial.print(F(" RMS vals\t"));
            Serial.print(rms_val[channel]);printTab();
            Serial.print(F("delta\t"));Serial.print(rms_tracker[channel].getDelta());
            Serial.print(F(" average\t"));Serial.println(rms_tracker[channel].getAvg());
        }
    }
}

void FeatureCollector::printPeakVals() {
    if (peak_active > 0) {
        for (int channel = 0; channel < num_peak_anas; channel++){
            Serial.print(name); Serial.print(" Peak vals\t");
            Serial.print(peak_val[channel]);printTab();
            Serial.print(F("delta\t"));Serial.print(peak_tracker[channel].getDelta());
            Serial.print(F(" average\t"));Serial.println(peak_tracker[channel].getAvg());
        }
    }
}

/////////////////////////////////// UPDATE / INIT //////////////////////////////////////
bool FeatureCollector::update(FFTManager1024 _fft[]) {
    if (autogain_active){updateAutogain(_fft);};
    if (microphone_active[0] + microphone_active[1] > 0) {
        if (last_update_timer > fc_update_rate) {
            for (int i = 0; i < num_rms_anas; i++) {
                calculateRMS(i);
            }
            for (int i = 0; i < num_peak_anas; i++) {
                calculatePeak(i);
            }
            printFeatures();
            last_update_timer = 0;
            return true;
        }
    }
    else {
        if (last_update_timer > 3000){
            Serial.print(name);
            Serial.println(F(" sorry the microphone does not work, not updating the feature collector"));
            last_update_timer = 0;
        }
    }
    return false;
}

////////////////////////////////////////////////////////////////////////
///////////////////////////// Tracking Features ////////////////////////
////////////////////////////////////////////////////////////////////////
void FeatureCollector::autogainTrackAvgRMS(double target, double tolerance, double _min, double _max) {
    ag_tracking_rms = true;
    ag_rms_idx = num_tracked_autogain_features;
    ag_target_vals[num_tracked_autogain_features] = target;
    ag_low_thresh[num_tracked_autogain_features] = target - (target * tolerance);
    ag_high_thresh[num_tracked_autogain_features] = target + (target * tolerance);
    ag_min_thresh[num_tracked_autogain_features] = _min;
    ag_max_thresh[num_tracked_autogain_features] = _max;
    num_tracked_autogain_features++;
}

void FeatureCollector::autogainTrackAvgPeak(double target, double tolerance, double _min, double _max) {
    ag_tracking_peak  = true;
    ag_peak_idx = num_tracked_autogain_features;
    ag_target_vals[num_tracked_autogain_features] = target;
    ag_low_thresh[num_tracked_autogain_features] = target - (target * tolerance);
    ag_high_thresh[num_tracked_autogain_features] = target + (target * tolerance);
    ag_min_thresh[num_tracked_autogain_features] = _min;
    ag_max_thresh[num_tracked_autogain_features] = _max;
    num_tracked_autogain_features++;
}

void FeatureCollector::autogainTrackOnRatio(double target, double tolerance, double _min, double _max, NeoGroup *_n) {
    ag_tracking_on_ratio = true;
    ag_on_ratio_idx = num_tracked_autogain_features;
    ag_target_vals[num_tracked_autogain_features] = target;
    ag_low_thresh[num_tracked_autogain_features] = target - (target * tolerance);
    ag_high_thresh[num_tracked_autogain_features] = target + (target * tolerance);
    ag_min_thresh[num_tracked_autogain_features] = _min;
    ag_max_thresh[num_tracked_autogain_features] = _max;
    neopixel_manager = _n;
    num_tracked_autogain_features++;
}

void FeatureCollector::autogainTrackClipping(double target, double tolerance, double _min, double _max) {
    ag_tracking_clipping = true;
    ag_clipping_idx = num_tracked_autogain_features;
    ag_target_vals[num_tracked_autogain_features] = target;
    ag_low_thresh[num_tracked_autogain_features] = target - (target * tolerance);
    ag_high_thresh[num_tracked_autogain_features] = target + (target * tolerance);
    ag_min_thresh[num_tracked_autogain_features] = _min;
    ag_max_thresh[num_tracked_autogain_features] = _max;
    num_tracked_autogain_features++;
}
#endif
