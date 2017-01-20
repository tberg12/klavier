package lm;

import tberg.murphy.indexer.HashMapIndexer;
import tberg.murphy.indexer.Indexer;
import io.PitchEventIO;

import java.io.File;
import java.util.*;

import tberg.murphy.counter.Counter;
import eval.PitchEventUtil;
import eval.PitchEventUtil.NoteState;
import tberg.murphy.arrays.a;

public class HeuristicChordAndPitchLM implements ChordAndPitchLM {

	public static final double PITCH_LM_SMOOTHING = 1e-20;
	public static final double CHORD_LM_SMOOTHING = 1e-5;
	public static final double CHORD_PITCH_SMOOTHING = 1e-5;
	
	public static final double SECONDS_PER_CHORD_FRAME = 0.5;
	public static final double CHOMRA_MIN_FRAC = 0.3;
	public static final int PRUNE_CHORDS_BELOW_COUNT = 5;
	
	public static class Chord {
		public final boolean[] chroma;
		public Chord(boolean[] chroma) {
			this.chroma = chroma;
		}
		public boolean equals(Object other) {
		    if (other instanceof Chord) {
		    	Chord that = (Chord) other;
		    	if (!Arrays.equals(this.chroma, that.chroma)) {
		    		return false;
		    	} else {
		    		return true;
		    	}
		    } else {
		    	return false;
		    }
		}
		public int hashCode() {
			return Arrays.hashCode(chroma);
		}
		public String toString() {
			StringBuffer buf = new StringBuffer();
			for (int i=0; i<chroma.length; ++i) {
				buf.append((chroma[i] ? "O" : "."));
			}
			return buf.toString();
		}
	}
	
    public static HeuristicChordAndPitchLM buildFromMIDICorpus(String path, float secondsPerPitchFrame, float extractBeforeOnsetMs, float extractAfterOnsetMs) {
    	File corpusRoot = new File(path);
        assert corpusRoot.isDirectory();
        List<NoteState[][]> songs = new ArrayList<NoteState[][]>();
        //for (File f : corpusRoot.listFiles) {
        {
            File f = new File(corpusRoot.getAbsolutePath() + "/IMSLP00001#1#Smartscore-10.2.1.mid");
            try {
                List<PitchEventIO.Event> events = PitchEventIO.readMIDIFile(f.getAbsolutePath());
                PitchEventIO.Event lastEvent = Collections.max(events, new Comparator<PitchEventIO.Event>() {
                    @Override
                    public int compare(PitchEventIO.Event event, PitchEventIO.Event event2) {
                        return (int)(event.offsetSec - event2.offsetSec);
                    }
                });
                NoteState[][] stateTypes = PitchEventUtil.convertToRoll(events, (int) Math.ceil(lastEvent.offsetSec / secondsPerPitchFrame), secondsPerPitchFrame,
                		extractBeforeOnsetMs, extractAfterOnsetMs);
                songs.add(stateTypes);

                //boolean[][] activations = PitchEventIO.eventsToActivationMatrix(events, 0.01f);
                //List<Object> context = new LinkedList<Object>();
                //for (boolean[] frame : activations) {
                //    Object repr = represent(frame);
                //    if (!scores.containsKey(context)) {
                //        scores.put(context, new HashMap<Object, Double>());
                //    }
                //    if (!scores.get(context).containsKey(repr)) {
                //        scores.get(context).put(repr, 0d);
                //    }
                //    scores.get(context).put(repr, scores.get(context).get(repr) + 1);
                //    context.add(repr);
                //    context.remove(0);
                //}
            } catch(Exception e) {
                e.printStackTrace();
                //assert false;
            }
            //if(count++ < 5) {
            //    break;
            //}
        }
        System.out.println(songs);
        return new HeuristicChordAndPitchLM(songs, secondsPerPitchFrame);
    }

	int pitchFramesPerChord;
	int numPitches;
	Indexer<Chord> chordIndexer;
	double[][][] logPitchLMProbs;
	double[][] logChordLMProbs;
	double[][][] logChordToPitchProbs;

	public HeuristicChordAndPitchLM(List<NoteState[][]> pitches, float secondsPerPitchFrame) {
    	this.pitchFramesPerChord = (int) Math.ceil(SECONDS_PER_CHORD_FRAME / secondsPerPitchFrame);
		this.numPitches = pitches.get(0)[0].length;
		this.chordIndexer = new HashMapIndexer<Chord>();
		indexChords(pitches);
		this.chordIndexer.lock();
		List<int[]> chords = extractChords(pitches);
		this.logPitchLMProbs = buildLogPitchLMProbs(pitches);
		this.logChordLMProbs = buildLogChordLMProbs(numChords(), chords);
		this.logChordToPitchProbs = buildLogChordToPitchProbs(numChords(), numPitches(), pitches, chords);
		
		System.out.println("Chords:");
		for (int c=0; c<numChords(); ++c) {
			System.out.println("chord: "+chordIndexer.getObject(c));
			StringBuffer pitchProbs = new StringBuffer();
			for (int p=0; p<numPitches(); ++p) {
				double onProb = Math.exp(logChordToPitchProbs[p][c][NoteState.ONSET.ordinal()]) + Math.exp(logChordToPitchProbs[p][c][NoteState.SUSTAIN.ordinal()]);
				pitchProbs.append((onProb > 0.33 ? (onProb > 0.66 ? "O" : "o") : ".") + ((p+1) % PitchEventIO.N_CHROMA == 0 ? "|" : ""));
			}
			System.out.println("pitch: "+pitchProbs.toString());
		}
    }

	public int numPitches() {
		return numPitches;
	}
	
	public int numChords() {
		return chordIndexer.size();
	}
	
	public int getPitchFramesPerChord() {
		return pitchFramesPerChord;
	}
	
	public List<int[]> extractChords(List<NoteState[][]> pitches) {
		final List<int[]> chords = new ArrayList<int[]>();
		for (NoteState[][] states : pitches) { 
			int numChordFrames = (int) Math.ceil((double) states.length / pitchFramesPerChord);
			int[] projected = new int[numChordFrames];
			for (int ct=0; ct<numChordFrames; ++ct) {
				int bestIndex = -1;
				double bestHamming = Double.POSITIVE_INFINITY;
				for (int cs=0; cs<chordIndexer.size(); ++cs) {
					double hamming = hammingDist(chordIndexer.getObject(cs).chroma, states, getPitchFrames(ct, states.length, pitchFramesPerChord));
					if (hamming < bestHamming) {
						bestHamming = hamming;
						bestIndex = cs;
					}
				}
				projected[ct] = bestIndex;
			}
			chords.add(projected);
		}
		
		return chords;
	}
	
	public double[][][] getLogPitchLMProbs() {
		return logPitchLMProbs;
	}
	
	public double[][] getLogChordLMProbs() {
		return logChordLMProbs;
	}
	
	public double[][][] getLogChordToPitchProbs() {
		return logChordToPitchProbs;
	}
	
	public static int[] getPitchFrames(int chordFrame, int numPitchFrames, int pitchFramesPerChord) {
		return a.enumerate(chordFrame*pitchFramesPerChord, Math.min(numPitchFrames, (chordFrame+1)*pitchFramesPerChord));
	}
	
	public static int getChordFrame(int pitchFrame, int pitchFramesPerChord) {
		return pitchFrame / pitchFramesPerChord;
	}
	
	private void indexChords(List<NoteState[][]> pitches) {
		Counter<Chord> chordCounter = new Counter<Chord>();
		for (NoteState[][] states : pitches) {
			int numChordFrames = (int) Math.ceil((double) states.length / pitchFramesPerChord);
			for (int ct=0; ct<numChordFrames; ++ct) {
				chordCounter.incrementCount(new Chord(projectToChord(states, getPitchFrames(ct, states.length, pitchFramesPerChord), CHOMRA_MIN_FRAC)), 1.0);
			}
		}
		chordCounter.pruneKeysBelowThreshold(PRUNE_CHORDS_BELOW_COUNT);
		for (Chord state : chordCounter.keySet()) {
			chordIndexer.getIndex(state);
		}
		System.out.println("Indexed "+chordIndexer.size()+" chords");
	}
	
	private static double[][][] buildLogPitchLMProbs(List<NoteState[][]> pitches) {
		int numPitches = pitches.get(0)[0].length;
		double[][][] logProbs = new double[pitches.get(0)[0].length][NoteState.values().length][NoteState.values().length];
		for (int pitchId=0; pitchId<numPitches; ++pitchId) {
			for (NoteState[][] pitchFrames : pitches) {
				for (int t=0; t<pitchFrames.length-1; ++t) {
					NoteState s0 = pitchFrames[t][pitchId];
					NoteState s1 = pitchFrames[t+1][pitchId];
					logProbs[pitchId][s0.ordinal()][s1.ordinal()] += 1.0;
				}
			}
		}
		a.normalizecoli(logProbs);
		a.scalei(logProbs, 1.0 - (PITCH_LM_SMOOTHING * NoteState.values().length));
		a.addi(logProbs, PITCH_LM_SMOOTHING);
		a.logi(logProbs);
		assert !a.hasnan(logProbs);
		assert !a.hasinf(logProbs);
		return logProbs;
	}
	
	private static double[][] buildLogChordLMProbs(int numChords, List<int[]> chords) {
		double[][] logProbs = new double[numChords][numChords];
		for (int[] chordFrames : chords) {
			for (int ct=0; ct<chordFrames.length-1; ++ct) {
				logProbs[chordFrames[ct]][chordFrames[ct+1]] += 1.0;
			}
		}
		a.normalizecoli(logProbs);
		a.scalei(logProbs, 1.0 - (CHORD_LM_SMOOTHING * numChords));
		a.addi(logProbs, CHORD_LM_SMOOTHING);
		a.logi(logProbs);
		assert !a.hasnan(logProbs);
		assert !a.hasinf(logProbs);
		return logProbs;
	}
	
	private double[][][] buildLogChordToPitchProbs(int numChords, int numPitches, List<NoteState[][]> pitches, List<int[]> chords) {
		double[][][] logProbs = new double[numPitches][numChords][NoteState.values().length];
		for (int d=0; d<chords.size(); ++d) {
			int[] chordFrames = chords.get(d);
			NoteState[][] states = pitches.get(d);
			for (int p=0; p<numPitches; ++p) {
				for (int ct=0; ct<chordFrames.length; ++ct) {
					int cs = chordFrames[ct];
					for (int t : getPitchFrames(ct, states.length, pitchFramesPerChord)) {
						logProbs[p][cs][states[t][p].ordinal()] += 1.0;
					}
				}
			}
		}
		a.normalizecoli(logProbs);
		a.scalei(logProbs, 1.0 - (CHORD_PITCH_SMOOTHING * NoteState.values().length));
		a.addi(logProbs, CHORD_PITCH_SMOOTHING);
		a.logi(logProbs);
		assert !a.hasnan(logProbs);
		assert !a.hasinf(logProbs);
		return logProbs;
	}
	
	private static boolean[] projectToChord(NoteState[][] states, int[] positions, double minFrac) {
		double[] chromaProbs = new double[PitchEventIO.N_CHROMA];
		for (int t : positions) {
			for (int p=0; p<states[t].length; ++p) {
				if (states[t][p] != NoteState.OFF) {
					chromaProbs[p % PitchEventIO.N_CHROMA]++;
				}
			}
		}
		a.normalizei(chromaProbs);
		boolean[] chroma = new boolean[PitchEventIO.N_CHROMA];
		for (int cp=0; cp<chroma.length; ++cp) {
			if (chromaProbs[cp] >= minFrac) {
				chroma[cp] = true;
			}
		}
		return chroma;
	}
	
	private static int hammingDist(boolean[] chroma0, boolean[] chroma1) {
		int dist = 0;
		for (int cp=0; cp<chroma0.length; ++cp) {
			if (chroma0[cp] != chroma1[cp]) {
				dist++;
			}
		}
		return dist;
	}
	
	private static double hammingDist(boolean[] chroma, NoteState[][] states, int[] positions) {
		int dist = 0;
		for (int t : positions) {
			boolean[] goldChroma = projectToChord(states, new int[] {t}, CHOMRA_MIN_FRAC);
			dist += hammingDist(chroma, goldChroma);
		}
		return ((double) dist) / positions.length;
	}

}
