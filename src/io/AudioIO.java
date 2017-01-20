package io;

public class AudioIO {
	
	public static class Wave implements java.io.Serializable {
		private static final long serialVersionUID = 1L;
		public final float[] amplitudes;
		public final int sampleRateHz;
		public final float lengthSec;
		public Wave(float[] amplitudes, int sampleRateHz, float lengthSec) {
			this.amplitudes = amplitudes;
			this.sampleRateHz = sampleRateHz;
			this.lengthSec = lengthSec;
		}
	}
	
	public static Wave readWaveFile(String path) {
		WaveFile wave = new WaveFile(path);
		short[] amplShort = wave.getSampleAmplitudes();
		int numChannels = wave.getWaveHeader().getChannels();
		float[] ampl = new float[amplShort.length/numChannels];
		for (int i=0; i<ampl.length; ++i) {
			float avg = 0.0f;
			for (int j=0; j<numChannels; ++j) {
				avg += ((float) amplShort[i*numChannels+j]) / numChannels;
			}
			ampl[i] = avg;
		}
		return new Wave(ampl, wave.getWaveHeader().getSampleRate(), (float) wave.length());
	}
	
	public static Wave takePrefix(Wave wave, float prefixLengthSec) {
		if (prefixLengthSec < wave.lengthSec) {
			int newNumSamples = (int) Math.floor(prefixLengthSec * wave.sampleRateHz);
			float[] ampl = new float[newNumSamples];
			for (int i=0; i<newNumSamples; ++i) {
				ampl[i] = wave.amplitudes[i];
			}
			return new Wave(ampl, wave.sampleRateHz, newNumSamples * (1.0f / wave.sampleRateHz));
		} else {
			return wave;
		}
	}

}
