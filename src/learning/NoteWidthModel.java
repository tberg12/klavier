package learning;

import learning.Preprocessing.SpectralAtomsAndEnvelopes;
import tberg.murphy.arrays.a;

public class NoteWidthModel {
	
	public static final float WIDTH_MODEL_EXPONENT = 1.0f;

	public static interface WidthScoresFactory {
		public int getMinWidth(int pitchId);
		public int getMaxWidth(int pitchId);
		public double[] getLogAllowedWidthsScores(int pitchId);
	}
	
	public static class SimpleWidthScoresFactory implements WidthScoresFactory {
		int minWidth;
		int maxWidth;
		double[] logAllowedWidthsScores;
		public SimpleWidthScoresFactory(SpectralAtomsAndEnvelopes spectralAtomsAndEnvelopes) {
			this.minWidth = spectralAtomsAndEnvelopes.getMinEnvelopeWidth();
			this.maxWidth = spectralAtomsAndEnvelopes.getMaxEnvelopeWidth();
			this.logAllowedWidthsScores = new double[(maxWidth - minWidth)+1];
			for (int w=minWidth; w<=maxWidth; ++w) {
				logAllowedWidthsScores[w-minWidth] = 1.0;
			}
			a.normalizei(logAllowedWidthsScores);
			a.logi(logAllowedWidthsScores);
			a.scalei(logAllowedWidthsScores, WIDTH_MODEL_EXPONENT);
		}
		public int getMinWidth(int pitchId) {
			return minWidth;
		}
		public int getMaxWidth(int pitchId) {
			return maxWidth;
		}
		public double[] getLogAllowedWidthsScores(int pitchId) {
			return logAllowedWidthsScores;
		}
	}
	
}
