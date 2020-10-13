#ifndef __FEATURE_COLLECTOR_H__
#define __FEATURE_COLLECTOR_H__

#include <Audio.h>
#include <ValueTrackerDouble.h>

#ifndef RMS_LOG_RESET_MIN
#define RMS_LOG_RESET_MIN 2000
#endif

#ifndef PEAK_LOG_RESET_MIN
#define PEAK_LOG_RESET_MIN 2000
#endif

class FeatureCollector {
    public:
        FeatureCollector(String _id);
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

        //////////////// Gain Tracking /////////////////////////
        void   setGain(double g, int channel);
        double getGain(int channel);
        void   linkAmplifier(AudioAmplifier * amp, double low, double high);

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
        void   resetRMSAvg(int channel);
        void   resetDominateRMSAvg(){return resetRMSAvg(dominate_channel);};

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
        double getPeakAvg(int channel){return peak_tracker[channel].getAvg();};
        double getDominatePeakAvg(){return getPeakAvg(dominate_channel);};
        void   resetPeakAvg(int channel);
        void   resetDominatePeakAvg(){return resetPeakAvg(dominate_channel);};

        ////////////////// Clipping
        double getClips(int channel){return clips[channel];};
        void   resetClips(int channel){clips[channel] = 0;};

        ////////////////// Min / Max
        double getPeakMin(int channel){return peak_tracker[channel].getMin();};
        double getPeakMax(int channel){return peak_tracker[channel].getMax();};
        double getDominatePeakMin(){return peak_tracker[dominate_channel].getMin();};
        double getDominatePeakMax(){return peak_tracker[dominate_channel].getMax();};

        void   printPeakVals();

        //////////////// General ///////////////////////
        bool update();
        void setDominateChannel(int c){dominate_channel = c;};

        int getNumRMSAnas() {return num_rms_anas;};
        int getNumPeakAnas() {return num_peak_anas;};

        // change the update rate for the FC
        void setUpdateRate(uint16_t t) {fc_update_rate = t;};

        //////////////// Printing ///////////////////////
        void setPrintRMS(bool p) {print_rms = p;};
        void setPrintRMSDebug(bool p) {print_rms_debug = p;};
        void setPrintPeak(bool p) {print_peak= p;};
        // void setPrintPeakDebug(bool p) {print_peak_debug = p;};

    private:
        ////////////////// General ///////////////////////
        // this is the channel which will be used when the "dominate"
        // feature is desired
        int dominate_channel    = 0;
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

        /////////////// Testing Functions ////////////////
        bool testMicrophoneRMS() {return testMicrophoneRMS(5000);};
        bool testMicrophoneRMS(int _dur);
        bool testMicrophonePeak() {return testMicrophonePeak(5000);};
        bool testMicrophonePeak(int _dur);

        //////////////// Gain Tracking ///////////////
        AudioAmplifier *amp_ana[4];

        double gain[2] = {1.0, 1.0};
        ValueTrackerDouble gain_tracker[2] = {ValueTrackerDouble("front_gain",&gain[0], 0.5), 
                ValueTrackerDouble("rear_gain", &gain[1], 0.5)};
        float min_gain[2];
        float max_gain[2];


        uint8_t audio_amp_add_idx = 0;

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

bool FeatureCollector::testMicrophoneRMS(int _dur) {
    // go through and gather 10 features from each channel and make sure it is picking up audio
    for (int i = 0; i < num_rms_anas; i++) {
        uint8_t readings = 0;
        double values = 0.0;
        unsigned long a_time = millis();
        Serial.print("Testing ");Serial.print(i);
        Serial.println(" Microphone using RMS");
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

bool FeatureCollector::testMicrophonePeak(int _dur) {
    // go through and gather 10 features from each channel and make sure it is picking up audio
    for (int i = 0; i < num_peak_anas; i++) {
        uint8_t readings = 0;
        double values = 0.0;
        unsigned long a_time = millis();
        printDivide();
        Serial.print("Testing ");Serial.print(name);Serial.println(" Microphone using Peak");
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
/////////////////////////////////////////////////////////////////
void FeatureCollector::calculatePeak(int channel) {
    for (int i = 0; i < num_peak_anas; i++) {
        bool avail = peak_ana[i]->available();
        if (peak_active && avail) {
            peak_val[channel] = peak_ana[i]->read();;
        }
    }
    if (peak_val[channel] >= 1.0) {
        clips[channel]++;
    }
    dprint(print_peak_debug, name);
    dprint(print_peak_debug, " Peaks (normal, pos_delta):\t");
    dprint(print_peak_debug, peak_val[channel]);
    dprint(print_peak_debug, "\t");
    dprintln(print_peak_debug, peak_tracker[channel].getPosDelta());
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
            dprint(print_rms_debug, "WARNING RMS is equal to 0 for channel: ");
            dprint(print_rms_debug, num_rms_anas);
        }
    }
}

void FeatureCollector::resetRMSAvg(int channel) {
    if (last_rms_reset[channel] > RMS_LOG_RESET_MIN) {
        rms_tracker[channel].getAvg(true);
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
            Serial.print("delta\t");Serial.print(rms_tracker[channel].getDelta());
            Serial.print(" average\t");Serial.println(rms_tracker[channel].getAvg());
        }
    }
}

void FeatureCollector::printPeakVals() {
    if (peak_active > 0) {
        for (int channel = 0; channel < num_peak_anas; channel++){
            Serial.print(name); Serial.print(" Peak vals\t");
            Serial.print(peak_val[channel]);printTab();
            Serial.print("delta\t");Serial.print(peak_tracker[channel].getDelta());
            Serial.print(" average\t");Serial.println(peak_tracker[channel].getAvg());
        }
    }
}

/////////////////////////////////// UPDATE / INIT //////////////////////////////////////
bool FeatureCollector::update() {
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
            Serial.println(" sorry the microphone does not work, not updating the feature collector");
            last_update_timer = 0;
        }
    }
    return false;
}
#endif
