package learning;

import learning.DenseSemiMarkovDP.Model;
import main.Main;
import tberg.murphy.arrays.a;

public class EnvelopeActivationModel implements Model {
	
	public static final float MIN_VOLUME = 0.0f;
	public static final float MAX_VOLUME = 1.0f;

    public static enum State {ON, OFF};

	float[] activations;
	float[] envelope;
	float[] envelopePrefixSqrNorms;
	int minOnWidth;
	int maxOnWidth;
	int[] allowedOffWidths;
	int[] allowedOnWidths;
	float[] logAllowedOnWidthsScores;
	float[][] logTransScores;
	float invEmissionVar;

	public EnvelopeActivationModel(float[] activations, float[] envelope, int minOnWidth, int maxOnWidth, float[] logAllowedOnWidthsScores, float[][] logTransScores) {
		this.activations = activations;
		this.envelope = envelope;
		this.envelopePrefixSqrNorms = new float[envelope.length];
		float currentSqrNorm = 0.0f;
		for (int i=0; i<envelope.length; ++i) {
			currentSqrNorm += envelope[i]*envelope[i];
			envelopePrefixSqrNorms[i] = currentSqrNorm;
		}
		this.allowedOffWidths = new int[] {1};
		this.minOnWidth = minOnWidth;
		this.maxOnWidth = maxOnWidth;
		this.allowedOnWidths = a.enumerate(minOnWidth, maxOnWidth+1);
		this.logAllowedOnWidthsScores = logAllowedOnWidthsScores;
		this.logTransScores = logTransScores;
		this.invEmissionVar = 1.0f / (float) Main.modelActivationVar;
	}

	public int length() {
		return activations.length;
	}

	public int numStates() {
		return State.values().length;
	}

	public int[] allowedWidths(int s) {
		if (s == State.OFF.ordinal()) {
			return allowedOffWidths;
		} else {
			return allowedOnWidths;
		}
	}

	public float logEdgePotential(int t, int prevS, int s) {
        return logTransScores[prevS][s];
	}
	
	public static float bestVolume(int t, int w, float[] activations, float[] envelope) {
		float envelopePrefixSqrNorm = 0.0f;
		for (int i=0; i<w; ++i) {
			envelopePrefixSqrNorm += envelope[i] * envelope[i];
		}
		float volume = 0.0f;
		for (int i=0; i<w; ++i) {
			volume += activations[t+i] * envelope[i];
		}
		volume /= envelopePrefixSqrNorm;
		volume = Math.max(MIN_VOLUME, volume);
		volume = Math.min(MAX_VOLUME, volume);
		return volume;
	}
	
	private float bestVolume(int t, int w) {
		float volume = 0.0f;
		for (int i=0; i<w; ++i) {
			volume += activations[t+i] * envelope[i];
		}
		volume /= envelopePrefixSqrNorms[w-1];
		volume = Math.max(MIN_VOLUME, volume);
		volume = Math.min(MAX_VOLUME, volume);
		return volume;
	}
	
	public float logNodePotential(int t, int w, int s) {
		float score = 0.0f;
		if (s == State.OFF.ordinal()) {
			float sumSqrDiff = 0.0f;
			for (int i=0; i<w; ++i) {
				float a = activations[t+i];
				sumSqrDiff += a * a;
			}
			score += -0.5f * sumSqrDiff * invEmissionVar;
		} else {
			score += logAllowedOnWidthsScores[w-minOnWidth];
			float volume = bestVolume(t, w);
			float sumSqrDiff = 0.0f;
			for (int i=0; i<w; ++i) {
				float diff = activations[t+i] - volume * envelope[i];
				sumSqrDiff += diff*diff;
			}
			score += -0.5f * sumSqrDiff * invEmissionVar;
		}
		return score;
	}
	
}
