#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <complex.h>
#include <sys/stat.h>

const int samp_rate = 240000;
const int output_rate = 48000;

const float decimate_transition_bw = 800;
const float ssb_bw = 3000;
const float amfm_bw = 12000;

#define hamming(x) (0.54-0.46*(2*M_PI*(0.5+(x/2.))));

int main(int argc, char *argv[]) { //syntax: rx <center_freq> <rx_freq> <l[sb]|u[sb]|a[m]|f[m]>
    if(argc != 4) {
        fprintf(stderr, "Usage: %s <center_freq> <rx_freq> <l[sb]|u[sb]|a[m]|f[m]>\n", argv[0]);
        return 1;
    }
    float center_freq = atof(argv[1]);
    float rx_freq = atof(argv[2]);
    char modulation = argv[3][0]; // l[sb], u[sb], a[m], f[m]

    const int decimate_factor = samp_rate / output_rate;
    const int decimate_taps_length = (int)(4.0f/(decimate_transition_bw/samp_rate)) | 1; //(should be odd)
    const int decimate_taps_middle = decimate_taps_length/2;
    const float dshift = 2*M_PI*(rx_freq - center_freq)/samp_rate;
    fprintf(stderr, "decimate_taps_len: %d\n", decimate_taps_length);
    //calculate filter taps
    complex float* decimate_taps = malloc(sizeof(complex float)*decimate_taps_length);
    const float decimate_cutoff_rate = (modulation == 'u'|| modulation == 'l') ? (ssb_bw/2.)/samp_rate : (amfm_bw/2.)/samp_rate;
    decimate_taps[decimate_taps_middle] = 2 * M_PI * decimate_cutoff_rate * hamming(0);
    for(int i=1; i<=decimate_taps_middle; i++) {
        decimate_taps[decimate_taps_middle-i] = decimate_taps[decimate_taps_middle+i] = (sin(2*M_PI*decimate_cutoff_rate*i)/i) * hamming((float)i/decimate_taps_middle);
    }
    float shift = 0;
    // make the filter asymmetric in case of SSB
    const complex float decimate_dshift = (modulation == 'u' ? 1:-1) * ((ssb_bw/2.)/samp_rate)*2*M_PI; // USB or LSB
    if(1 && (modulation == 'u' || modulation == 'l')) for(int i=0; i<decimate_taps_length; i++) {
        decimate_taps[i] *= sinf(shift) + I*cosf(shift);
        shift += decimate_dshift;
        if (shift > 2*M_PI) shift -= 2*M_PI;
    }
    // normalize filter
    float decimate_taps_sum = 0;
    for(int i=0; i<decimate_taps_length; i++) decimate_taps_sum += cabsf(decimate_taps[i]); 
    for(int i=0; i<decimate_taps_length; i++) decimate_taps[i] /= decimate_taps_sum;

#ifdef PRINTFREQZ
    FILE* filterfile = fopen("test.m", "w"); // <optional> <-- print the filter taps to file
    fprintf(filterfile, "#!/usr/bin/octave\n\nfreqz([\n");
    for(int i=0; i<decimate_taps_length; i++) fprintf(filterfile, " %f+(%f*i)\n", creal(decimate_taps[i]), cimag(decimate_taps[i]));
    fprintf(filterfile, "]);\ninput(\"\");\n");
    fchmod(fileno(filterfile), 0755);
    fclose(filterfile);                      // </optional>
#endif

    shift = 0;
    complex float* samplebuf = calloc(sizeof(complex float), decimate_taps_length);
    int decimate_counter = 0;
    float last_phi = 0;

    while(1) {
        // read 8bit[2] input from stdin   and  convert to complex float
        int in_a = getchar();
        int in_b = getchar();
        if (in_a == EOF || in_b == EOF) break;
        complex float sample = (float)in_a/(UCHAR_MAX/2.0)-1.0   +  I*( (float)in_b/(UCHAR_MAX/2.0)-1.0 );

        // oscillator & mixer
        shift += dshift; if (shift > 2*M_PI) shift-=2*M_PI;
        samplebuf[decimate_taps_length - decimate_factor+decimate_counter] = (sinf(shift) + I*cosf(shift)) * sample; // mix & buf

        // decimator: every N. sample   and   demodulate
        if(++decimate_counter >= decimate_factor) {
            // decimate
            complex float decim = 0;
            for(int i=0; i<=decimate_taps_length; i++) decim += samplebuf[i] * decimate_taps[i];
            memmove(samplebuf, samplebuf+decimate_factor, (decimate_taps_length-decimate_factor)*sizeof(complex float));
            decimate_counter = 0;

            // demodulate
            short soundout;
            switch (modulation) {
                case 'f': { // <-- fmdemod
                    float phi = cargf(decim);
                    float dphi = phi-last_phi;
                    last_phi = phi;
                    while (dphi < -M_PI) dphi += 2*M_PI;
                    while (dphi >  M_PI) dphi -= 2*M_PI;
                    soundout = (SHRT_MAX-1)*(dphi/M_PI);
                } break;
                case 'a': // <-- amdemod
                    soundout = cabsf(decim) * SHRT_MAX;
                break;
                default:  // <-- ssbdemod
                    soundout = crealf(decim) * SHRT_MAX;
                break;
            }
            // write demodulated sound to stdout
            fwrite(&soundout, sizeof(short), 1, stdout);
        }
    } // while(1)
}
