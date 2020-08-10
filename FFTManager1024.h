#ifndef __FFTManager1024_H__
#define __FFTManager1024_H__

#include "audio_utils.h"

#define NUM_FFT_BINS 512

#ifndef P_FLUX_VALS 
#define P_FLUX_VALS 0
#endif

#ifndef P_FFT_DEBUG
#define P_FFT_DEBUG 0
#endif

#ifndef P_CENTROID_VALS
#define P_CENTROID_VALS 0
#endif

#ifndef SMOOTH_CENTROID 
#define SMOOTH_CENTROID 1
#endif

#ifndef P_FFT_VALS
#define P_FFT_VALS 0 
#endif

class FFTManager1024 {
    public:
        //////////// init ///////////////
        FFTManager1024(String _name);
        void linkFFT(AudioAnalyzeFFT1024*r, bool w);

        // printers
        void   printFFTVals();

        // getters
        float getFFTRangeByIdx(uint16_t s, uint16_t e);
        float getFFTRangeByFreq(uint32_t s, uint32_t e);
        uint16_t    getHighestEnergyIdx(int start, int end);
        uint16_t    getHighestEnergyIdx();
        uint16_t    getHighestEnergyIdx(float array[], int start, int end);
        float getRelativeEnergy(uint16_t);
        float getFFTTotalEnergy();
        float getRelativeBinPos() {return relative_bin_pos;};

        float getCentroid();
        float getLastCentroid(){return last_centroid;};;
        float getCentroid(uint16_t min, uint16_t max);
        float getCentroidDelta();
        float getCentroidPosDelta();
        float getCentroidNegDelta();

        float getFlux();

        // setters
        void setupCentroid(bool v, float min, float max);
        void setCalculateFlux(bool v) {calculate_flux= v;};
        void setFluxActive(bool s) { calculate_flux = s;}; 

        void calculateSpectralFlux();
        void updateFFT(){calculateFFT();};

    private:
        String name = "";
        bool fft_active = true;
        AudioAnalyzeFFT1024*fft_ana;

        float raw_fft_vals[NUM_FFT_BINS];
        float fft_max_vals[NUM_FFT_BINS];
        float fft_vals[NUM_FFT_BINS];
        float last_fft_vals[NUM_FFT_BINS];

        void calculateFFT();

        float fft_tot_energy = 0.0;
        float relative_bin_pos = 0.0;

        uint16_t max_bin = NUM_FFT_BINS;// what is the highest index bin that we care about?
        uint16_t min_bin = 0;// what is the lowest index bin that we care about?

        bool calculate_centroid = false;
        float calculateCentroid();
        float centroid = 0.0;
        float last_centroid = 0.0;
        int centroid_min_bin = 0;
        int centroid_max_bin = NUM_FFT_BINS;

        bool calculate_flux = false;
        float calculateFlux();
        float flux = 0.0;
        elapsedMillis last_fft_reading;
        bool whitening_active = false;
        float whitening_floor = 0.0;
};

FFTManager1024::FFTManager1024(String _id) {
    name = _id;
}

void FFTManager1024::setupCentroid(bool v, float min, float max) {
    calculate_centroid = v;
    centroid_min_bin = uint16_t(min / 43);
    centroid_max_bin = uint16_t(max / 43);
    Serial.print("Now calculating the centroid for energy in bins ");
    Serial.print(centroid_min_bin);
    Serial.print(" through ");
    Serial.println(centroid_max_bin);
}

uint32_t getBinsMidFreq1024(int bin) {
    return (uint32_t)(bin * 43 + 22);
}

void printFreqRangeOfBin256(int idx) {
    Serial.print(idx * 172);
    Serial.print(" - ");
    Serial.println((idx + 1) * 172);
}

void printFreqRangeOfBin1024(int idx) {
    Serial.print(idx * 43);
    Serial.print(" - ");
    Serial.println((idx + 1) * 43);
}

void FFTManager1024::printFFTVals() {
    if (fft_active && (P_FFT_VALS || P_FLUX_VALS || P_CENTROID_VALS)) {
        Serial.print(name); Serial.println(" FFT vals\t");
        if (P_FFT_VALS) {
            // if (USE_SCALED_FFT) {Serial.print("Scaled ");}
            uint8_t w = 12;
            for (int l  = 0; l < w; l++) {
                Serial.println();
                Serial.print(l+min_bin); Serial.print("\t");
                for (int i = l; i < NUM_FFT_BINS; i = i + w) {
                    if (i != l) {
                        Serial.print(", ");
                        Serial.print(i);
                        Serial.print(":");
                    };
                    Serial.print(fft_vals[i], 5);
                }
            }
        }
    }
    if (calculate_flux == true && P_FLUX_VALS) {
        Serial.print("flux:\t");Serial.print(flux);Serial.println();
    }
    if (calculate_centroid == true && P_CENTROID_VALS) {
        Serial.print("centroid:\t");Serial.println(centroid);
    }
}

void FFTManager1024::linkFFT(AudioAnalyzeFFT1024*r, bool w) {
    fft_ana = r;
    fft_ana->averageTogether(4); // average together the readings from 10 FFT's before available() returns true, normally produces over 300 fft per second, this will be closer to 30
    fft_active = true;
    whitening_active = true;
};

float FFTManager1024::getRelativeEnergy(uint16_t idx) {
    // calculateFFT();
    if (fft_tot_energy > 0) {
        float val = 0.0;
        val = fft_vals[idx] / fft_tot_energy;
        // Serial.println(val);
        return val;
    }
    return 0.0;
}

float FFTManager1024::getFFTTotalEnergy() {
    // calculateFFT();
    if (fft_active) {
        return fft_tot_energy;
    }
    Serial.println("ERROR  - FFT IS NOT AN ACTIVE AUDIO FEATURE : "); Serial.println(name);
    return 0.0;
}


float FFTManager1024::getFFTRangeByIdx(uint16_t s, uint16_t e) {
    // calculateFFT();
    if (fft_active) {
        float sum = 0.0;
        for (int i = s; i <= e; i++){
            sum += fft_vals[i];
        }
        return sum / (e - s + 1);
    }
    return 0.0;
}

float FFTManager1024::getFFTRangeByFreq(uint32_t s, uint32_t e) {
    // calculateFFT();
    if (fft_active) {
        uint16_t start_idx = (uint16_t)(s / 43);
        uint16_t end_idx = (uint16_t)(e / 43);
        Serial.print(start_idx);
        Serial.print(" - - ");
        Serial.println(end_idx);
        return fft_ana->read(start_idx, end_idx);
    }
    return -1.0;
}

float FFTManager1024::calculateFlux() {
    float f = 0.0;
    if (last_fft_vals[0] != 0) {
        for (int i = 2; i < NUM_FFT_BINS; i++) {
            f += pow((fft_vals[i] - last_fft_vals[i]), 2);
        }
        if (f != 0.0) {
            flux = f;
        }
    }
    else {
        dprintln(P_FFT_DEBUG, "last_fft_vals[0] is equal to zero");
    }
    return flux;
}

float FFTManager1024::getFlux() {
    // calculateFFT();
    dprint(P_FLUX_VALS, "flux: ");
    dprintln(P_FLUX_VALS, flux);
    return flux;
}

/////////////// Calculate Features //////////////////////////////
float FFTManager1024::calculateCentroid() {
    float mags = 0.0;
    for (int i = centroid_min_bin; i < centroid_max_bin; i++) {
        // take the magnitude of all the bins
        // and multiply if by the mid frequency of the bin
        // then all it to the total cent value
        // we also have to remove the effect the FTT_SCALER has on the bins stored energy
        mags += raw_fft_vals[i] * getBinsMidFreq1024(i);
    }
    last_centroid = centroid;
    centroid = mags;
    dprint(P_CENTROID_VALS, "centroid : ");
    dprintln(P_CENTROID_VALS, centroid);
    return centroid;
}

float FFTManager1024::getCentroidDelta() {
    return centroid - last_centroid;
}

float FFTManager1024::getCentroidPosDelta() {
    if (last_centroid < centroid) {
        return centroid - last_centroid;
    }
    return 0.0;
}

float FFTManager1024::getCentroidNegDelta() {
    if (last_centroid > centroid) {
        return last_centroid - centroid;
    }
    return 0.0;
}

float FFTManager1024::getCentroid() {
    // calculateFFT();
    return centroid;
}

float FFTManager1024::getCentroid(uint16_t min, uint16_t max) {
    // calculateFFT();
    float mags = 0.0;
    for (int i = min; i < max; i++) {
        // take the magnitude of all the bins
        // and multiply if by the mid frequency of the bin
        // then all it to the total cent value
        mags += fft_vals[i] * getBinsMidFreq1024(i);
    }
    dprint(P_CENTROID_VALS, "centroid : ");
    dprintln(P_CENTROID_VALS, mags);
    return mags;
}

void FFTManager1024::calculateFFT() {
    if (fft_active && fft_ana->available() == true) {
        dprint(P_FFT_DEBUG, name);
        dprint(P_FFT_DEBUG, " FFT Available, ");
        dprint(P_FFT_DEBUG, " last FFT reading was ");
        dprint(P_FFT_DEBUG, (uint32_t)last_fft_reading);
        dprintln(P_FFT_DEBUG, " ms ago\n");
        fft_tot_energy = 0.0;
        if (calculate_flux == true) { // only do this if we need to in order to save time
            for (int i = 0; i < NUM_FFT_BINS; i++) {
                last_fft_vals[i] = fft_vals[i];
            }
        }
        for (int i = 0; i < NUM_FFT_BINS; i++) {
            fft_vals[i] = fft_ana->read(i);
            raw_fft_vals[i] = fft_vals[i];
            if (whitening_active) {
                if (fft_vals[i] > fft_max_vals[i]) {
                    fft_max_vals[i] = fft_vals[i];
                } else {
                    // ensure that the max values decay over time to prevent issues
                    fft_max_vals[i] = fft_max_vals[i] * 0.99;
                    // then make sure that the adjusted value is not less than the floor
                    if (fft_max_vals[i] < whitening_floor) {
                        fft_max_vals[i] = whitening_floor;
                    }
                }
                fft_vals[i] = fft_vals[i] / fft_max_vals[i];
            }
        }
        for (int i = 0; i < NUM_FFT_BINS; i++) {
            // fft_vals[i] *= fft_scaler;
            fft_tot_energy += fft_vals[i];
        }
        if (calculate_centroid == true) {calculateCentroid();};
        if (calculate_centroid == true && SMOOTH_CENTROID == true) {centroid = (centroid += last_centroid) * 0.5;};
        if (calculate_flux == true) {calculateFlux();};
        if (P_FFT_DEBUG) {
          printFFTVals();
        }
        last_fft_reading = 0;
    }
}


uint16_t FFTManager1024::getHighestEnergyIdx(float array[], int start, int end) {
    int highest = -1;
    float h_val = 0.0;
    for (int i = start; i <= end ; i++) {
        if (array[i] > h_val){
            highest = i;
            h_val = array[i];
        }
    }
    return highest;
}

uint16_t FFTManager1024::getHighestEnergyIdx(int start, int end){
    return getHighestEnergyIdx(fft_vals, start, end);
}

uint16_t FFTManager1024::getHighestEnergyIdx(){
    return getHighestEnergyIdx(fft_vals, min_bin, max_bin);
}

#endif // feature_Collector_h
