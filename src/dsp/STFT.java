package dsp;

import java.util.Arrays;

import main.Main;
import tberg.murphy.threading.BetterThreader;

public class STFT {
	
	public static enum WindowType {RECTANGULAR, BARTLETT, HANNING, HANNINGSQRT, HAMMING, BLACKMAN};
	
	public static float[][] rstftMagnitudes(final float[] realSignal, final int windowSize, WindowType windowType, final int hopSize, final int frequencyPrefixSize, final int frequencyOffset, int numThreads) {
		int numFrames = (int) Math.ceil(((float) realSignal.length) / ((float) hopSize));
		final float[][] magSpect = new float[numFrames][frequencyPrefixSize];
		final float[] windowWeights = makeWindow(windowSize, windowType);
		final int halfWindowSize = (int) Math.floor(windowSize / 2.0);
		BetterThreader.Function<Integer,Object> func = new BetterThreader.Function<Integer,Object>(){public void call(Integer t, Object ignore){
			float[] amplitudesWindow = new float[windowSize];
			for (int i=0; i<windowSize; ++i) {
				if (t*hopSize + (i - halfWindowSize) >= 0 && t*hopSize + (i - halfWindowSize) < realSignal.length) {
					amplitudesWindow[i] = realSignal[t*hopSize + (i - halfWindowSize)] * windowWeights[i];
				}
			}
			float[] mags = FFT.rdftMagnitudes(amplitudesWindow);
			for (int i=0; i<frequencyPrefixSize; ++i) {
				magSpect[t][i] = mags[i+frequencyOffset];
			}
		}};
		BetterThreader<Integer,Object> threader = new BetterThreader<Integer,Object>(func, Main.numThreads);
		for (int t=0; t<numFrames; ++t) threader.addFunctionArgument(t);
		threader.run();
		
		return magSpect;
	}
	
	public static float[] makeWindow(int nSamples, WindowType windowType) {
		// generate nSamples window function values
		// for index values 0 .. nSamples - 1
		int m = nSamples / 2;
		float r;
		float pi = (float) Math.PI;
		float[] w = new float[nSamples];
		switch (windowType) {
		case BARTLETT: // Bartlett (triangular) window
			for (int n = 0; n < nSamples; n++)
				w[n] = 1.0f - Math.abs(n - m) / m;
			break;
		case HANNING: // Hanning window
			r = pi / (m + 1);
			for (int n = -m; n < m; n++)
				w[m + n] = 0.5f + 0.5f * (float) Math.cos(n * r);
			break;
		case HANNINGSQRT: // Hanning window sqrt
			r = pi / (m + 1);
			for (int n = -m; n < m; n++)
				w[m + n] = (float) Math.sqrt(0.5f + 0.5f * Math.cos(n * r));
			break;
		case HAMMING: // Hamming window
			r = pi / m;
			for (int n = -m; n < m; n++)
				w[m + n] = 0.54f + 0.46f * (float) Math.cos(n * r);
			break;
		case BLACKMAN: // Blackman window
			r = pi / m;
			for (int n = -m; n < m; n++)
				w[m + n] = 0.42f + 0.5f * (float) Math.cos(n * r) + 0.08f
						* (float) Math.cos(2 * n * r);
			break;
		default: // Rectangular window function
			for (int n = 0; n < nSamples; n++)
				w[n] = 1.0f;
		}
		return w;
	}
	
	public static float[][] downsample(float[][] matrix, int hopSize) {
		int numDownsampledFrames = (int) Math.ceil(((float) matrix.length) / ((float) hopSize));
		int halfWindowSize = (int) Math.floor(hopSize / 2.0);
		float[][] result = new float[numDownsampledFrames][matrix[0].length];
		for (int t=0; t<numDownsampledFrames; ++t) {
			Arrays.fill(result[t], Float.NEGATIVE_INFINITY);
			for (int i=0; i<hopSize; ++i) {
				if (t*hopSize + (i - halfWindowSize) >= 0 && t*hopSize + (i - halfWindowSize) < matrix.length) {
					for (int j=0; j<result[t].length; ++j) {
						result[t][j] = Math.max(result[t][j], matrix[t*hopSize + (i - halfWindowSize)][j]);
					}
				}
			}
		}
		return result;
	}
	
	public static boolean[][] downsample(boolean[][] matrix, int hopSize) {
		int numDownsampledFrames = (int) Math.ceil(((float) matrix.length) / ((float) hopSize));
		int halfWindowSize = (int) Math.floor(hopSize / 2.0);
		boolean[][] result = new boolean[numDownsampledFrames][matrix[0].length];
		for (int t=0; t<numDownsampledFrames; ++t) {
			Arrays.fill(result[t], false);
			for (int i=0; i<hopSize; ++i) {
				if (t*hopSize + (i - halfWindowSize) >= 0 && t*hopSize + (i - halfWindowSize) < matrix.length) {
					for (int j=0; j<result[t].length; ++j) {
						result[t][j] = result[t][j] || matrix[t*hopSize + (i - halfWindowSize)][j];
					}
				}
			}
		}
		return result;
	}
	
}
