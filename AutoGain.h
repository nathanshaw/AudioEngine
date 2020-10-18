#ifndef __AUDIOENGINE_H__
#define __AUDIOENGINE_H__

#include "FeatureCollector.h"
#include "FFTManager1024.h"
#include <ValueTracker.h>

#ifndef P_AUTO_GAIN
#define P_AUTO_GAIN 0
#endif

class AutoGain {
  public:
    // Any feature collector or FFTManager which is passed into the init
    // functions will be linked to the feature collector
    // normal init
    AutoGain(FeatureCollector *, double min, double max, double max_gain_a);
    // one feature collector and one FFTManager
    AutoGain(FeatureCollector *, FFTManager1024 *, double min, double max, double max_gain_a);

    // for linked additional external feature collectors and FFTManagers
    void linkFeatureCollector(FeatureCollector *_fc );
    void linkFFTManager(FFTManager1024 * f);

    ///////////////////////////////////////////////////////////////////////////
    // for changing how often the auto-gain updates once the gain value is 
    // determined to be good
    void setUpdateRate(uint32_t r) {update_rate = r;};
    void setInitialUpdateRate(uint32_t d) {initial_update_rate = d;};
    // for setting up which features the AutoGain will be evaluating when an 
    // update is called
    //
    void trackAvgRMS(double target, double tolerance, double min, double max);
    void trackOnRatio(double target, double tolerance, double min, double max);
    // TODO - add track clips

    ////////////////////////////////////////////////////////////////////////////
    // the general update function that will look to see if the AutoGain is
    // ready to perform a gain adjustment. If a gain adjustment is determined to
    // be needed the update function will carry out all operations. In theory, 
    // after all the initial settings are passed into the AutoGain this is the
    // only function that needs to be called
    bool update();

    int calculateDominateChannel();
    void forceCalibartion();

  private:
    /////////////////// Linkages and Tracking /////////////////////
    // this will be set to true when the tracked values conform 
    // with the needed levels
    bool calculate(int channel);
    bool gain_good =                    false;
    uint8_t num_feature_collectors =    0;
    uint8_t num_fft_managers       =    0;

    bool fc_linked  =                   false;
    bool fft_linked =                   false;

    int rms_idx;
    bool tracking_rms =                 false;

    int on_ratio_idx;
    bool tracking_on_ratio =            false;

    uint8_t num_tracked_features =      0;

    //////////////////// Thresholds ///////////////////////////////
    // for the adjustment features
    double target_vals[2];
    double low_thresh[2];
    double high_thresh[2];
    double min_thresh[2];
    double max_thresh[2];

    // for the gain adjustments
    double min_gain;
    double max_gain;
    double max_gain_adj;

    //////////////////// Cost /////////////////////////////////////
    double calculateCost(double val, int idx);

    /////////////////// When to engage the routine ///////////////
    elapsedMillis last_autogain;

    // default update rate of 10 minutes
    uint32_t update_rate = 1000 * 60 * 10;

    // default initial update rate of 2 minutes
    uint32_t initial_update_rate = 1000 * 60 * 2;

    FeatureCollector        *fc;
    FFTManager1024          *fft;
};

float minf(float l, float b) {
    if (l < b){
        return l;
    }
    return b;
}

float maxf(float l, float b) {
    if (l > b){
        return l;
    }
    return b;
}

////////////////////////////////////////////////////////////////////////
///////////////////////////// Initialisation ///////////////////////////
////////////////////////////////////////////////////////////////////////

AutoGain::AutoGain(FeatureCollector *_fc, double min, double max, double max_gain_a) {
    if (_fc->isActive() == false) {
        Serial.println("WARNING - \
                the passed feature collector does not have amplifiers active");
    }
    min_gain     = min;
    max_gain     = max;
    max_gain_adj = max_gain_a;
    // this will add the feature collector to and
    // increase the num_feature_collectors by one
    linkFeatureCollector(_fc);
}

AutoGain::AutoGain(FeatureCollector *_fc, FFTManager1024 *_fft, double min, double max, double max_gain_a) {
    fft = _fft;
    AutoGain(_fc, min, max, max_gain_a);
}

void AutoGain::linkFeatureCollector(FeatureCollector *_fc ) {
    fc = _fc;
    fc_linked = true;
}

void AutoGain::linkFFTManager(FFTManager1024 * f) {
    fft = f;
    fft_linked = true;
}

void AutoGain::forceCalibartion() {
    // TODO
}

////////////////////////////////////////////////////////////////////////
///////////////////////////// Updates //////////////////////////////////
////////////////////////////////////////////////////////////////////////

bool AutoGain::calculate(int channel) {
    // now that we have established that enough time has passed for the
    // auto_gain to act, the auto_gain must act. The first order of business
    // is to update the tracked features.
    double cost = 0.0;
    double new_gain = 0.0;
    ////////////////////////// RMS //////////////////////////////////
    if (tracking_rms == true) {
        // IMPORTANT - for the RMS feature, if the value is half
        // or less than the expected value then the gain needs to be
        // halfed or doubled accordingly.
        cost += calculateCost(fc->getRMSAvg(channel), rms_idx);
        // need to reset the RMS average calculations for next time
        fc->resetRMSAvg(channel);
        dprint(P_AUTO_GAIN, "rms adjusted the cost to: ");
        dprintln(P_AUTO_GAIN, cost);
        // we aree halfing or doubling the signal
        // so this would equate to a max adjustment of half the signal when 
        // decreasing or double the signal when increasing
        if (cost > 0.0) {
            new_gain = fc->getGain(channel) + (fc->getGain(channel) * cost);
        } else if (cost < 0.0) {
            new_gain = fc->getGain(channel) - (fc->getGain(channel) * cost * 0.5) * -1.0;
        }
    }
    ////////////////////// On Ratio ////////////////////////////////
    if (tracking_on_ratio == true) {
        Serial.println("WARNING - tracking_on_ratio is not an implemented feature at this time, please use RMS instead");
        // IMPORTANT - for the on-ratio what we need is unknown
        // cost += calculateCost();
        // dprint(P_AUTO_GAIN, "on_ratio adjusted the cost to: ");
        // dprintln(P_AUTO_GAIN, cost);
    }
    new_gain = maxf(minf(max_gain, new_gain), min_gain);
    // only update the gain if it is different from the current gain
    if (new_gain != fc->getGain(channel)) {
        // this cost value tells us if gain adjustment is needed, if the value
        // is below 0.0 then the gain needs to be lowered, if the value is above
        // 1.0 then the gain needs to be raised.
        Serial.print("updating gain from ");Serial.print(fc->getGain(channel));Serial.print(" to ");
        Serial.print(new_gain);
        Serial.print(" with a cost of ");Serial.println(cost);// Serial.print(" auto_gain timer");Serial.print(last_autogain);
        // Serial.println(" / ");Serial.println(update_rate);
        // if there is a cost, then the gain needs to be adjusted
        fc->setGain(new_gain, channel);
        return true;
    }
    // if we got to this point and have not updated the gain we will
    return false;
}

bool AutoGain::update() {
    // first thing is to see if it is time to actually update the gain or to
    // exit the function. This logic is different depending of if the
    // gain is expected to be good or bad.

    // the gain will be false if the auto-gain has not been run before
    // or if the auto_gain has not been effecting and the tracked features
    // are outside of the established values

    // no matter what if less time than the initial_update_rate has
    // passed since the last update no update will occur, so the
    // function exits
    if (last_autogain < initial_update_rate) {
        return false;
    }

    // if the gain is within the accepted bounds, and less time than
    // the update_rate has passed, the function will exit
    if (gain_good == true && last_autogain < update_rate) {
        return false;
    }
    // if we made it this far then we are updating the gain
    // this is done here as it is more efficient due to the various
    // exiting conditions for the function
    dprintMinorDivide(P_AUTO_GAIN);
    dprintln(P_AUTO_GAIN, "\t\tAUTO GAIN INITIALISED");
    last_autogain = 0;

    for (int channel = 0; channel < num_tracked_features; channel++) {
        dprint(P_AUTO_GAIN, "starting auto_gain for channel: ");
        dprintln(P_AUTO_GAIN, channel);
        calculate(channel);
    }
    return true;
}

int AutoGain::calculateDominateChannel(){
    // first ensure that the gain is good for both channels before determining a dominate
    // second compare the RMS averages from each channel
    // then compare the Peak averages from each channel
    // subtract the RMS from the Peak to get the dynamic range metric
    //
    // for the second metric, noise calculate the average centroid for each channel
    // assume that the channel with the lower centroid has less noise
    //
    // for the third metric, look at the channel which is clipping less
    // TODO -- feature collector needs to be tracking clips
    int dom = 0;
    double _max = 0.0;

    for (int i = 0; i < fc->getNumRMSAnas(); i++) {
        double _temp = fc->getRMSAvg(i);
        fc->resetRMSAvg(i);
        if (_temp > _max) {
            dom = i;
            _max = _temp;
        }
    }
    Serial.println("WARNING!!!!!!!!!!!!!!!!!!!!!!");
    Serial.println("calculateDominateChannel is not fully implemented");

    return dom;
}
////////////////////////////////////////////////////////////////////////
///////////////////////////// Cost ...//////////////////////////////////
////////////////////////////////////////////////////////////////////////

double AutoGain::calculateCost(double val, int idx) {
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
    if (val > low_thresh[idx] && val < high_thresh[idx]) {
        return 0.0;
    }

    // if the value is below the min threshold. this results in the max cost of 1.0
    if (val < min_thresh[idx]) {
        dprintln(P_AUTO_GAIN, "Increasing autogain by the max amount");
        return 1.0;
    }

    // if the value is above the max threshold. this results in the min cost of -1.0
    if (val > max_thresh[idx]) {
        dprint(P_AUTO_GAIN, "Decreasing autogain by the max amount"); 
        return -1.0;
    }

    // if the value is between the low and min or high and max then return a scaled
    // cost according to how to lies within that range
    if (val < low_thresh[idx]) {
        // if thee gain needs to be raised
        cost = (val - min_thresh[idx]) / (low_thresh[idx] - min_thresh[idx]);
        dprint(P_AUTO_GAIN,"returning a cost of : ");
        dprintln(P_AUTO_GAIN, cost);
        return cost;
    }

    if (val > high_thresh[idx]) {
        cost = (1.0 - ((val - high_thresh[idx]) / (max_thresh[idx]- high_thresh[idx]))) * -1;
        dprint(P_AUTO_GAIN,"returning a cost of : ");
        dprintln(P_AUTO_GAIN, cost);
        return cost;
    }
    // the program should never reach this point...
    Serial.println("WARNING - there is something wrong with the calculate feature function, it should be returning something before the end of the function.");
    return 0.0;
}

/*
bool checkSongAutoGain() {
  adjustSongGainLedOnRatio();
  bool success = false;
  for (int i = 0; i < num_tracked_features; i++) {
    ///////////////////////////////////////////////////////////////
    // second check is to see if the song gain needs to be adjusted
    ///////////////////////////////////////////////////////////////
    // calculate the average peak values since the last auto-gain adjust
    double avg_song_peak = total_song_peaks[i] / num_song_peaks[i];
    double cost; // our cost variable
    dprint(P_AUTO_GAIN, "\n--------- song "); dprint(P_AUTO_GAIN, i); dprintln(P_AUTO_GAIN, " -------------");
    dprint(P_AUTO_GAIN, "total_song_peaks ");
    dprintln(P_AUTO_GAIN, total_song_peaks[i]);
    dprint(P_AUTO_GAIN, "num_song_peaks ");
    dprintln(P_AUTO_GAIN, (long) num_song_peaks[i]);
    // if the avg value is more than the max...
    if (avg_song_peak > MAX_SONG_PEAK_AVG) {
      // calculate cost between 0 and 1 with higher cost resulting in higher gain amplification
      cost = 1.0 - (MAX_SONG_PEAK_AVG / avg_song_peak);
      // calculate what the new song_gain will be
      double change = fc[i].gain * MAX_GAIN_ADJ* cost;
      fc[i].gain -= change;
      // ensure that what we have is not less than the min
      fc[i].gain = max(fc[i].gain, MIN_SONG_GAIN);
      dprint(P_AUTO_GAIN, "song gain decreased by ");
      dprint(P_AUTO_GAIN, change);
      dprint(P_AUTO_GAIN, " ");
      success = true;
    }
    // if the average value is less than the min....
    else if (avg_song_peak < MIN_SONG_PEAK_AVG) {
      dprintln(P_AUTO_GAIN);
      dprint(P_AUTO_GAIN, "avg_song_peak lower than MIN_SONG_PEAK_AVG ");
      dprintln(P_AUTO_GAIN, avg_song_peak);
      // calculate cost between 0 and 1 with higher cost resulting in higher gain attenuation
      cost = 1.0 - (MIN_SONG_PEAK_AVG / avg_song_peak);
      dprint(P_AUTO_GAIN, "cost : ");
      dprintln(P_AUTO_GAIN, cost);
      // calculate the new song gain
      double change = fc[i].gain * MAX_GAIN_ADJ* cost;
      fc[i].gain += change;
      // ensure what we have is not less than the max...
      fc[i].gain = min(fc[i].gain, MAX_SONG_GAIN);
      dprint(P_AUTO_GAIN, "song gain increased by ");
      dprint(P_AUTO_GAIN, change);
      dprint(P_AUTO_GAIN, " ");
      dprintln(P_AUTO_GAIN);
      success = true;
    }
    // now check the avg on time? todo

    ///////////////////////////////////////////////////////////////
    // last thing to do is reset the last_auto_gain_adjustment timer
    ///////////////////////////////////////////////////////////////
    total_song_peaks[i] = 0;
    num_song_peaks[i] = 0;
  }
  if (success) {
    return 1;
  } else {
    return 0;
  }
}
*/

////////////////////////////////////////////////////////////////////////
///////////////////////////// Tracking Features ////////////////////////
////////////////////////////////////////////////////////////////////////
void AutoGain::trackAvgRMS(double target, double tolerance, double _min, double _max) {
    tracking_rms = true;
    rms_idx = num_tracked_features;
    target_vals[num_tracked_features] = target;
    low_thresh[num_tracked_features] = target - (target * tolerance);
    high_thresh[num_tracked_features] = target + (target * tolerance);
    min_thresh[num_tracked_features] = _min;
    max_thresh[num_tracked_features] = _max;
    num_tracked_features++;
}

void AutoGain::trackOnRatio(double target, double tolerance, double _min, double _max) {
    tracking_on_ratio = true;
    on_ratio_idx = num_tracked_features;
    target_vals[num_tracked_features] = target;
    low_thresh[num_tracked_features] = target - (target * tolerance);
    high_thresh[num_tracked_features] = target + (target * tolerance);
    min_thresh[num_tracked_features] = _min;
    max_thresh[num_tracked_features] = _max;
    num_tracked_features++;
}

#endif // __AUDIO_ENGINE_H__
