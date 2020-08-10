#ifndef __FEATURE_COLLECTOR_H__
#define __FEATURE_COLLECTOR_H__

#include "audio_utils.h"
#include <Audio.h>

#ifndef RMS_LOG_RESET_MIN
#define RMS_LOG_RESET_MIN 2000
#endif

#ifndef PEAK_LOG_RESET_MIN
#define PEAK_LOG_RESET_MIN 2000
#endif

#ifndef MICROPHONE_TEST_DURATION
#define MICROPHONE_TEST_DURATION 2000
#endif

#ifndef FC_UPDATE_RATE
#define FC_UPDATE_RATE 33
#endif

#ifndef PRINT_PEAK_DEBUG
#define PRINT_PEAK_DEBUG 0
#endif

#ifndef PRINT_RMS_DEBUG
#define PRINT_RMS_DEBUG 0
#endif

class FeatureCollector {
    public:
        FeatureCollector(String _id);
        /////////////// General //////////////////////////
        String name = "";
        bool microphone_active[2] = {true, true};
        elapsedMillis last_update_timer;

        /////////////////// Microphone Testing /////////////////
        bool testMicrophone();

        /////////////////// Printing ///////////////////////////
        void printFeatures();

        ////////////////// Getting Functions ///////////////////
        String getName() {return name;};
        bool isActive();

        //////////////// Gain Tracking /////////////////////////
        void   setGain(double g, int channel);
        double getGain(int channel);
        bool   ampActive() {return amp_active;};
        void   linkAmplifier(AudioAmplifier * amp, double low, double high);

        //////////////// RMS /////////////////////////
        void   linkRMS(AudioAnalyzeRMS *r, bool print);

        double getRMS(int channel);
        double getDominateRMS(){return rms_val[dominate_channel];};

        bool   isRMSActive(){return rms_active;};

        double getRMSPosDelta(int channel){return rms_pos_delta[channel];};

        double getRMSAvg(int channel);
        double getDominateRMSAvg(){return getRMSAvg(dominate_channel);};

        void   resetRMSAvg(int channel);
        void   resetDominateRMSAvg(){return resetRMSAvg(dominate_channel);};

        void   printRMSVals();

        //////////////// Peak /////////////////////////
        void   linkPeak(AudioAnalyzePeak *r, bool print);
        bool   isPeakActive(){ return peak_active;};
        double getDominatePeak(){return peak_val[dominate_channel];};

        double getPeak(int channel);
        double getPeakPosDelta(int channel) {calculatePeak(channel);return peak_pos_delta[channel];};
        double getPeakAvg(int channel);
        double getDominatePeakAvg(){return getPeakAvg(dominate_channel);};

        double getClips(int channel){return clips[channel];};
        void   resetClips(int channel){clips[channel] = 0;};

        double getPeakMin(int channel){return peak_min[channel];};
        double getPeakMax(int channel){return peak_max[channel];};

        void   resetPeakAvg(int channel);
        void   resetDominatePeakAvg(){return resetPeakAvg(dominate_channel);};

        void   printPeakVals();

        //////////////// General ///////////////////////
        void update();
        void setDominateChannel(int c){dominate_channel = c;};

        int getNumRMSAnas() {return num_rms_anas;};
        int getNumPeakAnas() {return num_peak_anas;};

        //////////////// Printing ///////////////////////
        void autoPrintRMS(bool p) {print_rms = p;};
        void autoPrintPeak(bool p) {print_peak= p;};

    private:
        ////////////////// General ///////////////////////
        // this is the channel which will be used when the "dominate"
        // feature is desired
        int dominate_channel = 0;

        ////////////////// Printing //////////////////////
        bool print_rms      = false;
        bool print_peak     = false;

        /////////////// Testing Functions ////////////////
        bool testMicrophoneRMS();
        bool testMicrophonePeak();

        //////////////// Gain Tracking ///////////////
        AudioAmplifier *amp_ana[4];

        double gain[2] = {1.0, 1.0};
        double min_gain[2];
        double max_gain[2];

        bool amp_active = false;

        uint8_t audio_amp_add_idx = 0;

        //////////////// RMS /////////////////////////
        AudioAnalyzeRMS *rms_ana[2];
        uint8_t num_rms_anas = 0;
        bool rms_active = false;

        void calculateRMS(int channel);

        double rms_val[2];
        double rms_totals[2];
        unsigned long rms_readings[2];

        double rms_pos_delta[2];
        elapsedMillis last_rms_reset[2];

        //////////////// Peak /////////////////////////
        AudioAnalyzePeak *peak_ana[2];
        uint8_t num_peak_anas = 0;
        bool peak_active = false;

        void calculatePeak(int channel);

        double peak_val[2];
        double peak_min[2] = {9999.9, 9999.9};
        double peak_max[2] = {0.0, 0.0};

        // how many times the peak value is == 1.0
        uint32_t clips[2] = {0, 0};

        double peak_totals[2];
        unsigned long peak_readings[2];
        double peak_pos_delta[2];

        elapsedMillis last_peak_reset[2];
};

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

double FeatureCollector::getGain(int channel) {
    return gain[channel];
}

bool FeatureCollector::isActive() {
 return microphone_active[0] || microphone_active[1];
}

void FeatureCollector::linkAmplifier(AudioAmplifier * amp, double low, double high) { 
    max_gain[audio_amp_add_idx] = high;
    min_gain[audio_amp_add_idx] = low;
    if (audio_amp_add_idx < 4) {
        Serial.print("Linked an audio amplifier ");Serial.print(audio_amp_add_idx);printTab();
        amp_ana[audio_amp_add_idx] = amp;
        audio_amp_add_idx = audio_amp_add_idx + 1;
        Serial.println(audio_amp_add_idx);
    }
    else {
        Serial.println("ERROR, can't link audio amplifier, there are not enough available slots");
    }
}

FeatureCollector::FeatureCollector(String _id) {
    name = _id;
}


void FeatureCollector::setGain(double g, int channel) {
    gain[channel] = g;
    if (gain[channel] > max_gain[channel]) {
        max_gain[channel] = gain[channel];
    }
    if (gain[channel] < min_gain[channel]) {
        min_gain[channel] = gain[channel];
    }
    if (gain[channel] > max_gain[channel]) {
        gain[channel] = max_gain[channel];
    }
    else if (gain[channel] < min_gain[channel]) {
        gain[channel] = min_gain[channel];
    }
    amp_ana[channel]->gain(g);
}

bool FeatureCollector::testMicrophoneRMS() {
    // go through and gather 10 features from each channel and make sure it is picking up audio
    for (int i = 0; i < num_rms_anas; i++) {
        uint8_t readings = 0;
        double values = 0.0;
        unsigned long a_time = millis();
        Serial.print("Testing ");Serial.print(i);
        Serial.println(" Microphone using RMS");
        while (readings < 10 
                && millis() < a_time + MICROPHONE_TEST_DURATION) {
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
            Serial.println(" Microphone is good");
            microphone_active[i] = true;
            return true;
        } 
        Serial.println("\nERROR, ");
        Serial.print(name);Serial.println(" Microphone does not work");
        printDivideLn();
        microphone_active[i] = false;
    }
    return false;
}

bool FeatureCollector::testMicrophonePeak() {
    // go through and gather 10 features from each channel and make sure it is picking up audio
    for (int i = 0; i < num_peak_anas; i++) {
        uint8_t readings = 0;
        double values = 0.0;
        unsigned long a_time = millis();
        printDivide();
        Serial.print("Testing ");Serial.print(name);Serial.println(" Microphone using Peak");
        while (readings < 10 && millis() < a_time + MICROPHONE_TEST_DURATION) {
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
            Serial.println(" Microphone is good");
            microphone_active[i] = true;
            return true;
        } 
        Serial.println("\nERROR, ");
        Serial.print(name);Serial.println(" Microphone does not work");
        printMinorDivideLn();
        microphone_active[i] = false;
}
    return false;
}

bool FeatureCollector::testMicrophone () {
    if (rms_active) {
        return testMicrophoneRMS();
    } else if (peak_active) {
        return testMicrophonePeak();
    } else {
        Serial.println("Sorry unable to test microphone as neither the RMS or Peak feature is active");
        return false;
    }
}


//////////////// Update Functions ///////////////////////////////
void FeatureCollector::calculatePeak(int channel) {
    double last = peak_val[channel];
    double _peak_val = 0.0;
    for (int i = 0; i < num_peak_anas; i++) {
        bool avail = peak_ana[i]->available();
        if (peak_active && avail) {
            double _p = peak_ana[i]->read();
            if (_peak_val < _p){
                _peak_val = _p;
            }
        }
    }
    peak_val[channel] = _peak_val;
    if (peak_val[channel] >= 1.0) {
        clips[channel]++;
    }
    dprint(PRINT_PEAK_DEBUG, name);
    dprint(PRINT_PEAK_DEBUG, " Peaks (normal, pos_delta):\t");
    dprint(PRINT_PEAK_DEBUG, peak_val[channel]);
    peak_pos_delta[channel] = getPosDelta(last, peak_val[channel]);
    dprint(PRINT_PEAK_DEBUG, "\t");
    dprintln(PRINT_PEAK_DEBUG, peak_pos_delta[channel]);
    peak_totals[channel] += peak_val[channel];
    peak_readings[channel]++;
    if (peak_val[channel] > peak_max[channel]) {
        peak_max[channel] = peak_val[channel];
    } else if (peak_val[channel] < peak_min[channel]) {
        peak_min[channel] = peak_val[channel];
    }
}

void FeatureCollector::resetPeakAvg(int channel) {
    if (last_peak_reset[channel] > PEAK_LOG_RESET_MIN) {
        peak_totals[channel] = 0.0;
        peak_readings[channel] = 0;
        last_peak_reset[channel] = 0;
    }
}

void FeatureCollector::calculateRMS(int channel) {
    if (rms_active  && (rms_ana[channel]->available())) {
            double _rms = rms_ana[channel]->read();
        if (_rms > 0.0) {
            double temp = rms_val[channel];
            rms_val[channel] = _rms;
            rms_pos_delta[channel] = getPosDelta(temp, rms_val[channel]);
            rms_totals[channel] += rms_val[channel];
            rms_readings[channel]++;
        } else {
            dprint(PRINT_RMS_DEBUG, "WARNING RMS is equal to 0 for channel: ");
            dprint(PRINT_RMS_DEBUG, num_rms_anas);
        }
    }
}

void FeatureCollector::resetRMSAvg(int channel) {
    if (last_rms_reset[channel] > RMS_LOG_RESET_MIN) {
        rms_totals[channel] = 0.0;
        rms_readings[channel] = 0;
        last_rms_reset[channel] = 0;
    }
}

///////////////////// Getter functions ///////////////////////////////

double FeatureCollector::getRMS(int channel) {
    if (rms_active) {
        return rms_val[channel];
    }
    Serial.println("ERROR  - RMS IS NOT AN ACTIVE AUDIO FEATURE : "); Serial.println(name);
    return -1.0;
}

double FeatureCollector::getPeak(int channel) {
    if (peak_active) {
        return peak_val[channel];
    }
    Serial.println("ERROR  - Peak IS NOT AN ACTIVE AUDIO FEATURE : "); Serial.println(name);
    return -1.0;
}

double FeatureCollector::getPeakAvg(int channel) {
    if (peak_readings[channel] > 0 && peak_totals[channel] > 0) {
        return ((double)peak_totals[channel] / (double)peak_readings[channel]);
    }
    return peak_val[channel];
}

double FeatureCollector::getRMSAvg(int channel) {
    if (rms_readings[channel] > 0 && rms_totals[channel] > 0) {
        return ((double)rms_totals[channel] / (double)rms_readings[channel]);
    }
    return rms_val[channel];
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
            Serial.print(name); Serial.print(" RMS vals\t");
            Serial.print(rms_val[channel]);printTab();
            Serial.print("delta\t");Serial.print(rms_pos_delta[channel]);
            Serial.print(" average\t");Serial.println(getRMSAvg(channel));
        }
    }
}

void FeatureCollector::printPeakVals() {
    if (peak_active > 0) {
        for (int channel = 0; channel < num_peak_anas; channel++){
            Serial.print(name); Serial.print(" Peak vals\t");
            Serial.print(peak_val[channel]);printTab();
            Serial.print("delta\t");Serial.print(peak_pos_delta[channel]);
            Serial.print(" average\t");Serial.println(getPeakAvg(channel));
        }
    }
}

/////////////////////////////////// UPDATE / INIT //////////////////////////////////////
void FeatureCollector::update() {
    if (microphone_active[0] + microphone_active[1] > 0) {
        if (last_update_timer > FC_UPDATE_RATE) {
            last_update_timer = 0;
            for (int i = 0; i < num_rms_anas; i++) {
                calculateRMS(i);
            }
            for (int i = 0; i < num_peak_anas; i++) {
                calculatePeak(i);
            }
            printFeatures();
        }
    }
    else {
        if (last_update_timer > 3000){
            Serial.print(name);
            Serial.println(" sorry the microphone does not work, not updating the feature collector");
            last_update_timer = 0;
        }
    }
}
#endif
