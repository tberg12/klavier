package learning;

import io.AudioIO.Wave;
import io.DatasetIO;
import io.DatasetIO.Datum;
import io.PitchEventIO;

import java.util.ArrayList;
import java.util.List;

import tberg.murphy.threading.BetterThreader;
import tberg.murphy.tuple.Pair;
import main.Main;
import nmf.NMFUtilOpenCL;
import main.Main.NMFType;
import tberg.murphy.arrays.a;
import dsp.STFT;
import eval.PitchEventUtil;

public class Preprocessing {
	public static final boolean SPECT_POW = false;
	public static final boolean SPECT_DIVIDE_BY_MAX = true;
	public static final boolean SPECT_STRETCH_TO_UNIT = false;
	public static final boolean SPECT_LOG = false;
	public static final boolean SPECT_SCALE = false;
	public static final boolean SPECT_BUCKET = false;
	
	public static final boolean ACTIVATIONS_POW = false;
	public static final boolean ACTIVATIONS_DIVIDE_BY_MAX = false;
	public static final boolean ACTIVATIONS_STRETCH_TO_UNIT = false;
	public static final boolean ACTIVATIONS_LOG = false;
	public static final boolean ACTIVATIONS_SCALE = false;
	public static final boolean ACTIVATIONS_BUCKET = false;
	
	public static final float POW = 2.0f;
	public static final float SCALE = 1e2f;
	public static final float SCALE_BEFORE_LOG = 1e1f;
	public static final float ADD_BEFORE_LOG = 1e0f;
	public static final float CLIP_MIN = 1e-15f;
	
	public static class ProcDatum implements java.io.Serializable {
		private static final long serialVersionUID = 1L;
		public final Datum rawDatum;
		public final float[][] spect;
		public final float[][] activations;
		public float secondsPerFrame;
		public ProcDatum(Datum rawDatum, float[][] spect, float[][] activations, float secondsPerFrame) {
			System.out.println("Length sec: "+rawDatum.wave.lengthSec);
			System.out.println("Num frames: "+activations.length);
//			System.out.println("Frames extrap length sec: "+activations.length * secondsPerFrame);
			this.rawDatum = rawDatum;
			this.spect = spect;
			this.activations = activations;
			this.secondsPerFrame = secondsPerFrame;
		}
	}
	
	public static class SpectralAtomsAndEnvelopes implements java.io.Serializable  {
			private static final long serialVersionUID = 1L;
			float[][] spectralAtoms;
			float[][] envelopes;
			int minWidth;
			int maxWidth;
			public SpectralAtomsAndEnvelopes(float[][] spectralAtoms, float[][] envelopes, int minWidth, int maxWidth) {
				this.spectralAtoms = spectralAtoms;
				this.envelopes = JaggedMatrixAverager.padToNumColumnsLastValue(envelopes, maxWidth);
				this.minWidth = minWidth;
				this.maxWidth = maxWidth;
			}
			public float[][] getSpectralAtoms() {
				return spectralAtoms;
			}
			public float[][] getEnvelopes() {
				return envelopes;
			}
			public int getMinEnvelopeWidth() {
				return minWidth;
			}
			public int getMaxEnvelopeWidth() {
				return maxWidth;
			}
	}
	
	public static class JaggedMatrixAverager {
		float[][] average;
		float[][] count;
		public JaggedMatrixAverager(int numRows) {
			this.average = new float[numRows][0];
			this.count = new float[numRows][0];
		}
		public void observe(float[][] mat) {
			observe(mat, 1.0f);
		}
		public void observe(float[][] mat, float n) {
			for (int r=0; r<mat.length; ++r) {
				observe(r, mat[r], n);
			}
		}
		public void observe(int r, float[] row) {
			observe(r, row, 1.0f);
		}
		public void observe(int r, float[] row, float n) {
			if (row.length > average[r].length) {
				{
					float[] longer = new float[row.length];
					System.arraycopy(average[r], 0, longer, 0, average[r].length);
					average[r] = longer;
				}
				{
					float[] longer = new float[row.length];
					System.arraycopy(count[r], 0, longer, 0, count[r].length);
					count[r] = longer;
				}
			}
			for (int c=0; c<row.length; ++c) {
				average[r][c] = ((count[r][c] * average[r][c])/(count[r][c] + n)) + ((n * row[c]) / (count[r][c] + n));
				count[r][c] += n;
			}
		}
		public float[][] getAverage() {
			return average;
		}
	    public static float[][] pad(float[][] mat, float val) {
	    	int maxLength = 0;
	    	for (int i=0; i<mat.length; ++i) {
	    		if (mat[i].length > maxLength) {
	    			maxLength = mat[i].length; 
	    		}
	    	}
	    	return padToNumColumns(mat, val, maxLength);
	    }
	    public static float[][] padLastValue(float[][] mat) {
	    	int maxLength = 0;
	    	for (int i=0; i<mat.length; ++i) {
	    		if (mat[i].length > maxLength) {
	    			maxLength = mat[i].length; 
	    		}
	    	}
	    	return padToNumColumnsLastValue(mat, maxLength);
	    }
		public static float[][] padToNumColumnsLastValue(float[][] mat, int numCols) {
			float[][] result = new float[mat.length][numCols];
			for (int i=0; i<mat.length; ++i) {
				for (int j=0; j<numCols; ++j) {
					if (j < mat[i].length) {
						result[i][j] = mat[i][j];
					} else {
						result[i][j] = mat[i][ mat[i].length-1];
					}
				}
			}
			return result;
		}
		public static float[][] padToNumColumns(float[][] mat, float val, int numCols) {
			float[][] result = new float[mat.length][numCols];
			for (int i=0; i<mat.length; ++i) {
				for (int j=0; j<numCols; ++j) {
					if (j < mat[i].length) {
						result[i][j] = mat[i][j];
					} else {
						result[i][j] = val;
					}
				}
			}
			return result;
		}
	}
	
	public static SpectralAtomsAndEnvelopes buildInitialSpectralAtomsAndEnvelopes(final String notesBasePath, final String notesPaths, List<String> instrNames) {
		final JaggedMatrixAverager spectralAtomsAvg = new JaggedMatrixAverager(PitchEventIO.N_MIDI_PITCH_IDS);
		final JaggedMatrixAverager envelopesAvg = new JaggedMatrixAverager(PitchEventIO.N_MIDI_PITCH_IDS);
		final float[] secondsPerFrame = new float[] {0.0f};
		BetterThreader.Function<String,Object> func = new BetterThreader.Function<String,Object>(){public void call(String instrName, Object ignore){
			Pair<Pair<float[][],float[][]>,Float> spectralAtomsAndEnvelopeAndSecondsPerFrame = buildSpectralAtomsAndEnvelopes(notesBasePath, notesPaths, instrName);
			Pair<float[][],float[][]> spectralAtomsAndEnvelope = spectralAtomsAndEnvelopeAndSecondsPerFrame.getFirst();
			synchronized (spectralAtomsAvg) {
				secondsPerFrame[0] = spectralAtomsAndEnvelopeAndSecondsPerFrame.getSecond();
				spectralAtomsAvg.observe(spectralAtomsAndEnvelope.getFirst());
				envelopesAvg.observe(spectralAtomsAndEnvelope.getSecond());
			}
		}};
		BetterThreader<String,Object> threader = new BetterThreader<String,Object>(func, Main.numThreads);
		for (String instrName : instrNames) threader.addFunctionArgument(instrName);
		threader.run();
		int minWidth = (int) Math.ceil(Main.modelMinNoteLengthMs * 1e-3f / secondsPerFrame[0]);
		int maxWidth = (int) Math.ceil(Main.modelMaxNoteLengthMs * 1e-3f / secondsPerFrame[0]);
		return new SpectralAtomsAndEnvelopes(spectralAtomsAvg.getAverage(), envelopesAvg.getAverage(), minWidth, maxWidth);
	}
	
	private static Pair<Pair<float[][],float[][]>,Float> buildSpectralAtomsAndEnvelopes(String notesBasePath, String notesPaths, String instrName) {
		JaggedMatrixAverager spectralAtomsAvg = new JaggedMatrixAverager(PitchEventIO.N_MIDI_PITCH_IDS);
		JaggedMatrixAverager envelopesAvg = new JaggedMatrixAverager(PitchEventIO.N_MIDI_PITCH_IDS);
		float secondsPerFrame = 0.0f;
		for (String notesPath : notesPaths.trim().split(":")) {
			List<Datum> data = DatasetIO.readLabeledData(notesBasePath+"/"+instrName+"/"+notesPath, Float.POSITIVE_INFINITY);
			System.out.println(notesBasePath+"/"+instrName+"/"+notesPath+" num: "+data.size());
			for (Datum datum : data) {
				Pair<float[][],Float> spectAndSecondsPerFrame = computeSpect(datum.wave);
				float[][] spect = spectAndSecondsPerFrame.getFirst();
				secondsPerFrame = spectAndSecondsPerFrame.getSecond();
				for (int i=0; i<datum.events.size(); ++i) {
					int noteIndex = datum.events.get(i).noteIndex;
					if (noteIndex >= 0 && noteIndex < PitchEventIO.N_MIDI_PITCH_IDS) {

						List<float[]> subSpectList = new ArrayList<float[]>();
						for (int t : PitchEventUtil.framesByOverlap(datum.events.get(i).onsetSec + ((float) Main.initSpectEnvNotesOffsetMs * 1e-3f), datum.events.get(i).offsetSec + ((float) Main.initSpectEnvNotesOffsetMs * 1e-3f), secondsPerFrame, 0, spect.length-1)) {
//						for (int t : PitchEventUtil.framesByCenterWithinInterval(datum.events.get(i).onsetSec + ((float) Main.notesOffsetMs * 1e-3f), datum.events.get(i).offsetSec + ((float) Main.notesOffsetMs * 1e-3f), secondsPerFrame, 0, spect.length-1)) {
							subSpectList.add(spect[t]);
						}
						if (subSpectList.isEmpty()) continue;
						float[][] subSpect = subSpectList.toArray(new float[0][]);

						Pair<float[][],float[][]> atomAndEnvelope = null;
						if (Main.nmfType == NMFType.KL) {
							atomAndEnvelope = NMFUtilOpenCL.nmfKL(subSpect, 1, 20, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps);
						} else if (Main.nmfType == NMFType.Beta) {
							atomAndEnvelope = NMFUtilOpenCL.nmfBeta(subSpect, 1, 20, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfBeta);
						} else if (Main.nmfType == NMFType.LogNormal) {
							int iters = 2000;
							float startStepSize = 1e-8f;
							float endStepSize = 1e-8f;
							float init = 1e-5f;
							atomAndEnvelope = NMFUtilOpenCL.nmfLogNormalExpGrad(a.scale(a.onesFloat(1, subSpect[0].length), init), true, a.scale(a.onesFloat(subSpect.length, 1), init), true, subSpect, startStepSize, endStepSize, iters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfC);
							while (a.hasinf(atomAndEnvelope.getFirst()) || a.hasnan(atomAndEnvelope.getFirst()) || a.hasinf(atomAndEnvelope.getSecond()) || a.hasnan(atomAndEnvelope.getSecond())) {
								startStepSize *= 0.8;
								endStepSize *= 0.8;
								atomAndEnvelope = NMFUtilOpenCL.nmfLogNormalExpGrad(a.scale(a.onesFloat(1, subSpect[0].length), init), true, a.scale(a.onesFloat(subSpect.length, 1), init), true, subSpect, startStepSize, endStepSize, iters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfC);
							}
						} else {
							atomAndEnvelope = NMFUtilOpenCL.nmfL2(subSpect, 1, 20, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps);
						}
						float[] atom = atomAndEnvelope.getFirst()[0];
						float[] envelope = a.transpose(atomAndEnvelope.getSecond())[0];

//						a.scalei(atom, 1.0f / a.max(atom));
//						a.scalei(atom, 1.0f / a.sum(atom));
						a.scalei(atom, 1.0f / (float) Math.sqrt(a.sum(a.sqr(atom))));
						spectralAtomsAvg.observe(noteIndex, atom);

						a.scalei(envelope, 1.0f / a.max(envelope));
						envelopesAvg.observe(noteIndex, envelope);
					}
				}
			}
		}
		return Pair.makePair(Pair.makePair(spectralAtomsAvg.getAverage(), envelopesAvg.getAverage()), secondsPerFrame);
	}
	
	public static List<ProcDatum> preprocessData(final List<Datum> data, final SpectralAtomsAndEnvelopes spectralAtomsAndEnvelopes) {
		final ProcDatum[] resultArray = new ProcDatum[data.size()];
		BetterThreader.Function<Integer,Object> func = new BetterThreader.Function<Integer,Object>(){public void call(Integer i, Object ignore){
			resultArray[i] = preprocessDatum(data.get(i), spectralAtomsAndEnvelopes.getSpectralAtoms());
		}};
		BetterThreader<Integer,Object> threader = new BetterThreader<Integer,Object>(func, 1);
		for (int i=0; i<data.size(); ++i) threader.addFunctionArgument(i);
		threader.run();
		List<ProcDatum> result = new ArrayList<ProcDatum>();
		for (ProcDatum procDatum : resultArray) result.add(procDatum);
		return result;
	}
	
	private static ProcDatum preprocessDatum(Datum datum, float[][] spectralAtoms) {
		Pair<float[][],Float> spectAndSecondsPerFrame = computeSpect(datum.wave);
		float[][] spect = spectAndSecondsPerFrame.getFirst();
		float secondsPerFrame = spectAndSecondsPerFrame.getSecond();
		float[][] activations = null;
		if (Main.nmfType == NMFType.KL) {
			activations = (Main.useGpu ? NMFUtilOpenCL.nmfKLGPU(spectralAtoms, false, a.onesFloat(spect.length, spectralAtoms.length), true, spect, Main.nmfNumIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps).getSecond() : NMFUtilOpenCL.nmfKL(spectralAtoms, false, a.onesFloat(spect.length, spectralAtoms.length), true, spect, Main.nmfNumIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps).getSecond());
		} else if (Main.nmfType == NMFType.Beta) {
			activations = (Main.useGpu ? NMFUtilOpenCL.nmfBetaGPU(spectralAtoms, false, a.onesFloat(spect.length, spectralAtoms.length), true, spect, Main.nmfNumIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfBeta).getSecond() : NMFUtilOpenCL.nmfBeta(spectralAtoms, false, a.onesFloat(spect.length, spectralAtoms.length), true, spect, Main.nmfNumIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfBeta).getSecond());
		} else if (Main.nmfType == NMFType.LogNormal) {
			float startStepSize = 1e-7f;
			float endStepSize = 1e-7f;
			float init = 1e-5f;
			activations = (Main.useGpu ? NMFUtilOpenCL.nmfLogNormalExpGradGPU(spectralAtoms, false, a.scale(a.onesFloat(spect.length, spectralAtoms.length), init), true, spect, startStepSize, endStepSize, Main.nmfNumIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfC).getSecond() : NMFUtilOpenCL.nmfLogNormalExpGrad(spectralAtoms, false, a.scale(a.onesFloat(spect.length, spectralAtoms.length), init), true, spect, startStepSize, endStepSize, Main.nmfNumIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfC).getSecond());
		}  else {
			activations = (Main.useGpu ? NMFUtilOpenCL.nmfL2GPU(spectralAtoms, false, a.onesFloat(spect.length, spectralAtoms.length), true, spect, Main.nmfNumIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps).getSecond() : NMFUtilOpenCL.nmfL2(spectralAtoms, false, a.onesFloat(spect.length, spectralAtoms.length), true, spect, Main.nmfNumIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps).getSecond());
		}
		
		if (ACTIVATIONS_POW) {
			a.powi(activations, POW);
		}
		if (ACTIVATIONS_DIVIDE_BY_MAX) {
			divideByMax(activations);
		}
		if (ACTIVATIONS_STRETCH_TO_UNIT) {
			stretchToUnit(activations);
		}
		if (ACTIVATIONS_LOG) {
			logSpace(activations);
			stretchToUnit(activations);
		}
		if (ACTIVATIONS_SCALE) {
			a.scalei(activations, SCALE);
		}
		if (ACTIVATIONS_BUCKET) {
			bucket(activations);
		}
		return new ProcDatum(datum, spect, activations, secondsPerFrame);
	}
	
	public static Pair<float[][],Float> computeSpect(Wave wave) {
		long start = System.nanoTime();
		float[][] spect = STFT.downsample(STFT.rstftMagnitudes(wave.amplitudes, Main.stftWindowSize, Main.stftWindowType, Main.stftHopSize, Main.stftFrequencyPrefixSize, Main.stftFrequencyOffset, Main.numThreads), Main.stftDownsampleSpectrogramHopSize);
		if (SPECT_POW) {
			a.powi(spect, POW);
		}
		if (SPECT_DIVIDE_BY_MAX) {
			divideByMax(spect);
		}
		if (SPECT_STRETCH_TO_UNIT) {
			stretchToUnit(spect);
		}
		if (SPECT_LOG) {
			logSpace(spect);
			stretchToUnit(spect);
		}
		if (SPECT_SCALE) {
			a.scalei(spect, SCALE);
		}
		if (SPECT_BUCKET) {
			bucket(spect);
		}
		float secondsPerSample = (1.0f / wave.sampleRateHz);
		float samplesPerFrame = Main.stftHopSize * Main.stftDownsampleSpectrogramHopSize;
		float secondsPerFrame = secondsPerSample * samplesPerFrame;
		long end = System.nanoTime();
		System.out.println("Compute spect time: "+(end - start)/(1e9)+"s");
		return Pair.makePair(spect, secondsPerFrame);
	}
	
	public static void divideByMax(float[][] activations) {
		a.scalei(activations, 1.0f / a.max(a.max(activations)));
	}
	
	public static void divideByValAndClip(float[][] activations, float val) {
		for (int t=0; t<activations.length; ++t) {
			for (int p=0; p<activations[t].length; ++p) {
				activations[t][p] = Math.min(activations[t][p], val);
			}
		}
		a.scalei(activations, 1.0f / val);
		for (int t=0; t<activations.length; ++t) {
			for (int p=0; p<activations[t].length; ++p) {
				activations[t][p] = Math.max(activations[t][p], CLIP_MIN);
			}
		}
	}

	public static void logSpace(float[][] activations) {
		a.scalei(activations, SCALE_BEFORE_LOG);
		a.addi(activations, ADD_BEFORE_LOG);
		a.logi(activations);
	}

	public static void stretchToUnit(float[][] activations) {
		float min = a.min(a.min(activations));
		float max = a.max(a.max(activations));
		for (int t=0; t<activations.length; ++t) {
			for (int p=0; p<activations[t].length; ++p) {
				activations[t][p] = (activations[t][p] - min) / (max - min);
				if (Float.isNaN(activations[t][p])) activations[t][p] = 0.0f;
			}
		}
		for (int t=0; t<activations.length; ++t) {
			for (int p=0; p<activations[t].length; ++p) {
				activations[t][p] = Math.max(activations[t][p], CLIP_MIN);
			}
		}
	}
	
	public static void bucket(float[][] activations) {
		for (int t=0; t<activations.length; ++t) {
			for (int p=0; p<activations[t].length; ++p) {
				activations[t][p] = (float) Math.max(SCALE * 1e-4, Math.round(activations[t][p]));
			}
		}
	}

}
