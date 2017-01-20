package learning;

import io.PitchEventIO.Event;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import eval.PitchEventUtil;
import eval.PitchEventUtil.NoteState;
import nmf.NMFUtil;
import tberg.murphy.arrays.a;
import tberg.murphy.threading.BetterThreader;
import tberg.murphy.tuple.Pair;
import learning.NoteTransitionModel.TransScoresFactory;
import learning.NoteTransitionModel.IMSLPTransScoresFactory;
import learning.NoteWidthModel.WidthScoresFactory;
import learning.NoteWidthModel.SimpleWidthScoresFactory;
import learning.Preprocessing.ProcDatum;
import learning.EnvelopeActivationModel;
import learning.EnvelopeActivationModel.State;
import learning.Preprocessing.SpectralAtomsAndEnvelopes;
import main.Main;
import main.Main.NMFType;

public class Model {
	
	SpectralAtomsAndEnvelopes spectralAtomsAndEnvelopes;
	WidthScoresFactory widthScoresFactory;
	TransScoresFactory transScoreFactory;
	
	public Model(SpectralAtomsAndEnvelopes spectralAtomsAndEnvelopes) {
		this.spectralAtomsAndEnvelopes = spectralAtomsAndEnvelopes;
		this.widthScoresFactory = new SimpleWidthScoresFactory(spectralAtomsAndEnvelopes);
		this.transScoreFactory = new IMSLPTransScoresFactory();
	}

	public List<List<Event>> train(List<ProcDatum> data, boolean updateAtoms, boolean updateEnvelopes, int numEMIters, int numApproxIters) {
		float[][] envelopes = a.copy(spectralAtomsAndEnvelopes.getEnvelopes());
		float[][] atoms = a.copy(spectralAtomsAndEnvelopes.getSpectralAtoms());
		float[][] envelopesPrior = a.copy(envelopes);
		float[][] atomsPrior = a.copy(atoms);
		
		List<float[][]> spects = new ArrayList<float[][]>();
		for (ProcDatum datum : data) spects.add(datum.spect);
		List<float[][]> activations = new ArrayList<float[][]>();
		for (ProcDatum datum : data) activations.add(a.copy(datum.activations));
		
		Pair<List<NoteState[][]>,List<float[][]>> noteStatesAndPreActivations = null;
		for (int emIter=0; emIter<numEMIters; ++emIter) {
			System.out.println("Learning iter "+emIter+"..");
			
			if (updateAtoms && emIter > 0) activations = estepActivations(atoms, noteStatesAndPreActivations.getSecond(), activations, spects);
			noteStatesAndPreActivations = estepNotes(activations, envelopes, noteStatesAndPreActivations == null ? null : noteStatesAndPreActivations.getFirst());
			
			for (int approxIter=0; approxIter<numApproxIters; ++approxIter) {
				System.out.println("Update iter "+approxIter+"..");
				activations = estepActivations(atoms, noteStatesAndPreActivations.getSecond(), activations, spects);
				noteStatesAndPreActivations = estepNotes(activations, envelopes, noteStatesAndPreActivations.getFirst());
			}
			
			if (updateEnvelopes) envelopes = mstepEnvelopes(noteStatesAndPreActivations.getFirst(), activations, envelopesPrior);
			if (updateAtoms) atoms = mstepAtoms(atomsPrior, atoms, activations, spects);
		}
		
		this.spectralAtomsAndEnvelopes = new SpectralAtomsAndEnvelopes(atoms, envelopes, spectralAtomsAndEnvelopes.getMinEnvelopeWidth(), spectralAtomsAndEnvelopes.getMaxEnvelopeWidth());
		
		if (updateAtoms && numEMIters > 0) activations = estepActivations(atoms, noteStatesAndPreActivations.getSecond(), activations, spects);
		noteStatesAndPreActivations = estepNotes(activations, envelopes, noteStatesAndPreActivations == null ? null : noteStatesAndPreActivations.getFirst());
		for (int approxIter=0; approxIter<numApproxIters; ++approxIter) {
			System.out.println("Update iter "+approxIter+"..");
			activations = estepActivations(atoms, noteStatesAndPreActivations.getSecond(), activations, spects);
			noteStatesAndPreActivations = estepNotes(activations, envelopes, noteStatesAndPreActivations.getFirst());
		}
    	List<NoteState[][]> noteStates = noteStatesAndPreActivations.getFirst();
    	
    	List<List<Event>> result = new ArrayList<List<Event>>();
    	for (int i=0; i<data.size(); ++i) {
    		result.add(PitchEventUtil.convertToEvents(noteStates.get(i), activations.get(i), envelopes, data.get(0).secondsPerFrame));
    	}
		
		return result;
	}
	
	public SpectralAtomsAndEnvelopes getSpectralAtomsAndEnvelopes() {
		return spectralAtomsAndEnvelopes;
	}

	public List<List<Event>> decode(List<ProcDatum> data, int numApproxIters) {
		float[][] envelopes = a.copy(spectralAtomsAndEnvelopes.getEnvelopes());
		float[][] atoms = a.copy(spectralAtomsAndEnvelopes.getSpectralAtoms());
		
		List<float[][]> spects = new ArrayList<float[][]>();
		for (ProcDatum datum : data) spects.add(datum.spect);
		List<float[][]> activations = new ArrayList<float[][]>();
		for (ProcDatum datum : data) activations.add(a.copy(datum.activations));
		
		Pair<List<NoteState[][]>,List<float[][]>> noteStatesAndPreActivations = estepNotes(activations, envelopes, null);
		for (int approxIter=0; approxIter<numApproxIters; ++approxIter) {
			System.out.println("Update iter "+approxIter+"..");
			activations = estepActivations(atoms, noteStatesAndPreActivations.getSecond(), activations, spects);
			noteStatesAndPreActivations = estepNotes(activations, envelopes, noteStatesAndPreActivations.getFirst());
		}
    	List<NoteState[][]> noteStates = noteStatesAndPreActivations.getFirst();
    	
    	List<List<Event>> result = new ArrayList<List<Event>>();
    	for (int i=0; i<data.size(); ++i) {
    		result.add(PitchEventUtil.convertToEvents(noteStates.get(i), activations.get(i), envelopes, data.get(0).secondsPerFrame));
    	}
		
		return result;
	}
	
	private static List<float[][]> estepActivations(float[][] atoms, List<float[][]> preActivations, List<float[][]> activations, List<float[][]> spects) {
		System.out.println("Update activations..");
		
		Pair<float[][],int[]> totalSpectAndNumRows = appendMatrices(spects);
		float[][] totalSpect = totalSpectAndNumRows.getFirst();
		int[] appendNumRows = totalSpectAndNumRows.getSecond();
		float[][] totalActivation = appendMatrices(activations).getFirst();
		float[][] totalPreActivation = appendMatrices(preActivations).getFirst();
		
		if (Main.nmfType == NMFType.Beta) {
//			int nmfIters = 400;
//			double activationPriorWeight = 5e0;
//			double startStepSize = 5e-4;
//			double endStepSize = 5e-4;
//			int nmfIters = 400;
//			double activationPriorWeight = 5e0;
//			double startStepSize = 2e0;
//			double endStepSize = 2e0;
//			double startStepSize = 2e1;
//			double endStepSize = 2e1;
			
			if (Main.activationsPriorWeight == 0.0) {
				Pair<float[][],float[][]> atomsAndActivations = (Main.useGpu ? NMFUtil.nmfBetaGPU(atoms, false, totalActivation, true, totalSpect, Main.updateActivationsIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfBeta) : NMFUtil.nmfBeta(atoms, false, totalActivation, true, totalSpect, Main.updateActivationsIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfBeta));
				return separateMatrices(atomsAndActivations.getSecond(), appendNumRows);
			} else {
				Pair<float[][],float[][]> atomsAndActivations = (Main.useGpu ? NMFUtil.nmfBetaL2PriorExpGradGPU(atoms, atoms, 0.0f, false, totalPreActivation, totalActivation, (float) Main.activationsPriorWeight, true, totalSpect, (float) Main.updateActivationsStepSize, (float) Main.updateActivationsStepSize, Main.updateActivationsIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfBeta) : NMFUtil.nmfBetaL2PriorExpGrad(atoms, atoms, 0.0f, false, totalPreActivation, totalActivation, (float) Main.activationsPriorWeight, true, totalSpect, (float) Main.updateActivationsStepSize, (float) Main.updateActivationsStepSize, Main.updateActivationsIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfBeta));
				return separateMatrices(atomsAndActivations.getSecond(), appendNumRows);
			}
		} else {
			int nmfIters = 400;
			double activationPriorWeight = 5e0;
			double startStepSize = 1e-1;
			double endStepSize = 1e-1;
			
			if (activationPriorWeight == 0.0) {
				Pair<float[][],float[][]> atomsAndActivations = (Main.useGpu ? NMFUtil.nmfKLGPU(atoms, false, totalActivation, true, totalSpect, nmfIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps) : NMFUtil.nmfKL(atoms, false, totalActivation, true, totalSpect, nmfIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps));
				return separateMatrices(atomsAndActivations.getSecond(), appendNumRows);
			} else {
				Pair<float[][],float[][]> atomsAndActivations = (Main.useGpu ? NMFUtil.nmfKLL2PriorExpGradGPU(atoms, atoms, 0.0f, false, totalPreActivation, totalActivation, (float) activationPriorWeight, true, totalSpect, (float) startStepSize, (float) endStepSize, nmfIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps) : NMFUtil.nmfKLL2PriorExpGrad(atoms, atoms, 0.0f, false, totalPreActivation, totalActivation, (float) activationPriorWeight, true, totalSpect, (float) startStepSize, (float) endStepSize, nmfIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps));
				return separateMatrices(atomsAndActivations.getSecond(), appendNumRows);
			}
		}
	}
	
	private static float[][] mstepAtoms(float[][] atomsPrior, float[][] atoms, List<float[][]> activations, List<float[][]> spects) {
		System.out.println("Update atoms..");
		
		float[][] totalSpect = appendMatrices(spects).getFirst();
		float[][] totalActivation = appendMatrices(activations).getFirst();
		
		if (Main.nmfType == NMFType.Beta) {
//			int nmfIters = 400;
//			double atomsPriorWeight = 1e2;
//			double startStepSize = 1e-5;
//			double endStepSize = 1e-5;
//			int nmfIters = 400;
//			double atomsPriorWeight = 2e-3;
//			double startStepSize = 1e-1;
//			double endStepSize = 1e-1;
//			double startStepSize = 5e-1;
//			double endStepSize = 5e-1;
			
			if (Main.atomsPriorWeight == 0.0) {
				Pair<float[][],float[][]> atomsAndActivations = (Main.useGpu ? NMFUtil.nmfBetaGPU(atoms, true, totalActivation, false, totalSpect, Main.updateAtomsIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfBeta) : NMFUtil.nmfBeta(atoms, true, totalActivation, false, totalSpect, Main.updateAtomsIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfBeta));
				atoms = atomsAndActivations.getFirst();
			} else {
				Pair<float[][],float[][]> atomsAndActivations = (Main.useGpu ? NMFUtil.nmfBetaL2PriorExpGradGPU(atomsPrior, atoms, (float) Main.atomsPriorWeight, true, totalActivation, totalActivation, 0.0f, false, totalSpect, (float) Main.updateAtomsStepSize, (float) Main.updateAtomsStepSize, Main.updateAtomsIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfBeta) : NMFUtil.nmfBetaL2PriorExpGrad(atomsPrior, atoms, (float) Main.atomsPriorWeight, true, totalActivation, totalActivation, 0.0f, false, totalSpect, (float) Main.updateAtomsStepSize, (float) Main.updateAtomsStepSize, Main.updateAtomsIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps, (float) Main.nmfBeta));
				atoms = atomsAndActivations.getFirst();
			}
		} else {
			int nmfIters = 400;
			double atomsPriorWeight = 1e2;
			double startStepSize = 1e-3;
			double endStepSize = 1e-3;
			
			if (atomsPriorWeight == 0.0) {
				Pair<float[][],float[][]> atomsAndActivations = (Main.useGpu ? NMFUtil.nmfKLGPU(atoms, true, totalActivation, false, totalSpect, nmfIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps) : NMFUtil.nmfKL(atoms, true, totalActivation, false, totalSpect, nmfIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps));
				atoms = atomsAndActivations.getFirst();
			} else {
				Pair<float[][],float[][]> atomsAndActivations = (Main.useGpu ? NMFUtil.nmfKLL2PriorExpGradGPU(atomsPrior, atoms, (float) atomsPriorWeight, true, totalActivation, totalActivation, 0.0f, false, totalSpect, (float) startStepSize, (float) endStepSize, nmfIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps) : NMFUtil.nmfKLL2PriorExpGrad(atomsPrior, atoms, (float) atomsPriorWeight, true, totalActivation, totalActivation, 0.0f, false, totalSpect, (float) startStepSize, (float) endStepSize, nmfIters, (float) Main.nmfSilenceEps, (float) Main.nmfMinEps));
				atoms = atomsAndActivations.getFirst();
			}
		}
		return atoms;
	}
	
	private static float[][] mstepEnvelopes(List<NoteState[][]> noteStates, List<float[][]> activations, float[][] envelopesPrior) {
		System.out.println("Update envelopes..");
		
		float envelopesPriorWeight = 0.2f;
		
		float[][] envelopesNum = a.scale(envelopesPrior, envelopesPriorWeight);
		float[][] envelopesDenom = new float[envelopesPrior.length][];
		for (int pitchId=0; pitchId<envelopesPrior.length; ++pitchId) {
			envelopesDenom[pitchId] = new float[envelopesPrior[pitchId].length];
			Arrays.fill(envelopesDenom[pitchId], envelopesPriorWeight);
		}
		int numNotes = activations.get(0)[0].length;
		for (int d=0; d<noteStates.size(); ++d) {
			for (int pitchId=0; pitchId<numNotes; ++pitchId) {
				boolean noteOn = false;
				int start = -1;
				for (int t=0; t<noteStates.get(d).length; ++t) {
					NoteState s = noteStates.get(d)[t][pitchId];
					if (s == NoteState.ONSET) {
						noteOn = true;
						start = t;
					}
					if (noteOn && s == NoteState.OFF) {
						noteOn = false;
						float[] envelopeToken = new float[t-start];
						for (int tt=start; tt<t; ++tt) {
							envelopeToken[tt-start] = activations.get(d)[tt][pitchId];
						}
						float bestVolume = EnvelopeActivationModel.bestVolume(0, envelopeToken.length, envelopeToken, envelopesPrior[pitchId]);
						for (int i=0; i<envelopeToken.length; ++i) {
							envelopesNum[pitchId][i] += bestVolume * envelopeToken[i];
							envelopesDenom[pitchId][i] += bestVolume * bestVolume;
						}
					}
				}
			}
		}

		float[][] envelopes = new float[envelopesPrior.length][];
		for (int pitchId=0; pitchId<envelopesPrior.length; ++pitchId) {
			envelopes[pitchId] = new float[envelopesPrior[pitchId].length];
			for (int i=0; i<envelopesPrior[pitchId].length; ++i) {
				envelopes[pitchId][i] = envelopesNum[pitchId][i] / envelopesDenom[pitchId][i];
			}
		}
		return envelopes;
	}
	
	private Pair<List<NoteState[][]>,List<float[][]>> estepNotes(List<float[][]> activations, float[][] envelopes, List<NoteState[][]> prevNoteStates) {
		System.out.println("Update notes..");
		
		long start = System.nanoTime();
		
		List<NoteState[][]> noteStates = new ArrayList<NoteState[][]>();
		List<float[][]> preActivations = new ArrayList<float[][]>();
		for (int d=0; d<activations.size(); ++d) {
			NoteState[][] prevNoteStatesLocal;
			if (prevNoteStates == null) {
				Pair<NoteState[][],float[][]> prevNoteStateAndPreActivation = viterbi(activations.get(d), envelopes, null);
				prevNoteStatesLocal = prevNoteStateAndPreActivation.getFirst();
			} else {
				prevNoteStatesLocal = prevNoteStates.get(d);
			}
			Pair<NoteState[][],float[][]> noteStateAndPreActivation = viterbi(activations.get(d), envelopes, prevNoteStatesLocal);
			noteStates.add(noteStateAndPreActivation.getFirst());
			preActivations.add(noteStateAndPreActivation.getSecond());
		}
		
		long end = System.nanoTime();
		System.out.println("Compute time: "+(end - start)/(1e9)+"s");
		
		return Pair.makePair(noteStates, preActivations);
	}
	
	private Pair<NoteState[][],float[][]> viterbi(final float[][] activations, final float[][] envelopes, NoteState[][] prevNoteStates) {
		final int numNotes = activations[0].length;
		final NoteState[][] noteStates = new NoteState[activations.length][numNotes];
		final float[][] preActivations = new float[activations.length][numNotes];
		for (int t=0; t<noteStates.length; ++t) {
			Arrays.fill(noteStates[t], NoteState.OFF);
		}
		final float[][] transActivations = a.transpose(activations);
		BetterThreader.Function<Integer,Object> func = new BetterThreader.Function<Integer,Object>(){public void call(Integer pitchId, Object ignore){
			DenseSemiMarkovDP.Model model = new EnvelopeActivationModel(transActivations[pitchId], envelopes[pitchId], widthScoresFactory.getMinWidth(pitchId), widthScoresFactory.getMaxWidth(pitchId), a.toFloat(widthScoresFactory.getLogAllowedWidthsScores(pitchId)), a.toFloat(transScoreFactory.getLogTransScores(pitchId)));
			List<Pair<Integer,Pair<Integer,Integer>>> decode = DenseSemiMarkovDP.viterbiDecode(model, pitchId);
			for (Pair<Integer,Pair<Integer,Integer>> segment : decode) {
				State s = State.values()[segment.getFirst()];
				int startT = segment.getSecond().getFirst();
				int endT = segment.getSecond().getSecond();
				if (s == State.ON) {
					float bestVolume = EnvelopeActivationModel.bestVolume(startT, endT-startT, transActivations[pitchId], envelopes[pitchId]);
					for (int t=startT; t<endT; ++t) {
						noteStates[t][pitchId] = NoteState.SUSTAIN;
						preActivations[t][pitchId] = bestVolume * envelopes[pitchId][t-startT];
					}
					noteStates[startT][pitchId] = NoteState.ONSET;
				}
			}
		}};
		BetterThreader<Integer,Object> threader = new BetterThreader<Integer,Object>(func, Main.numThreads);
		for (int pitchId=0; pitchId<numNotes; ++pitchId) threader.addFunctionArgument(pitchId);
		threader.run();
		return Pair.makePair(noteStates, preActivations);
	}
	
	private static List<float[][]> separateMatrices(float[][] appendedMat, int[] numRows) {
		List<float[][]> mats = new ArrayList<float[][]>();
		int i=0;
		for (int m=0; m<numRows.length; ++m) {
			mats.add(Arrays.copyOfRange(appendedMat, i, i+numRows[m]));
			i += numRows[m];
		}
		return mats;
	}
	
	private static Pair<float[][],int[]> appendMatrices(List<float[][]> mats) {
		int[] numRows = new int[mats.size()];
		for (int i=0; i<mats.size(); ++i) {
			numRows[i] = mats.get(i).length;
		}
		int totalRows = a.sum(numRows);
		float[][] result = new float[totalRows][];
		int i=0;
		for (float[][] mat : mats) {
			for (float[] row : mat) {
				result[i] = row;
				i++;
			}
		}
		return Pair.makePair(result, numRows);
	}
	
}
