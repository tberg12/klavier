package io;

import tberg.murphy.arrays.a;
import learning.Preprocessing.SpectralAtomsAndEnvelopes;

/**
 * @author jda
 */
public class SynthIO {

    public static SpectralAtomsAndEnvelopes buildSpectralAtomsAndEnvelopes() {
        final int N_NOTES = PitchEventIO.N_MIDI_PITCH_IDS;
        final int N_BUCKETS = 1024; // TODO(jda) get from somewhere
        final int MIN_FRAMES = 10;
        final int MAX_FRAMES = 100;
        
        final double EQUAL_TEMP_MULTIPLIER = 1.059463;
        final double A0_FREQ = 27.5 * .093;
    	
        //return new float[0][];
        //return a.onesFloat(N_NOTES, N_BUCKETS);
        float[][] spectralAtoms = a.zerosFloat(N_NOTES, N_BUCKETS);
        for (int i = 0; i < N_NOTES; i++) {
            double freq = A0_FREQ * Math.pow(EQUAL_TEMP_MULTIPLIER, i);
            System.out.println(freq);
            int bin = (int)(freq);
            if (bin < spectralAtoms[i].length) {
                spectralAtoms[i][bin] = 1;
                spectralAtoms[i][bin+1] = .5f;
                spectralAtoms[i][bin-1] = .5f;
            }
        }
        
        //return a.onesFloat(N_NOTES, MAX_FRAMES);
        //return base.getEnvelopes(instrName);
        float[][] envelopes = a.zerosFloat(N_NOTES, MAX_FRAMES);
        for (int i = 0; i < N_NOTES; i++) {
            for (int j = 0; j < MIN_FRAMES; j++) {
                envelopes[i][j] = 1;
            }
        }
    	
        return new SpectralAtomsAndEnvelopes(spectralAtoms, envelopes, MIN_FRAMES, MAX_FRAMES);
    }
}
