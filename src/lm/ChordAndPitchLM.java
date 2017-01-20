package lm;

import java.util.List;

import eval.PitchEventUtil.NoteState;

public interface ChordAndPitchLM {

	public abstract int numPitches();

	public abstract int numChords();

	public abstract int getPitchFramesPerChord();

	public abstract List<int[]> extractChords(List<NoteState[][]> pitches);

	public abstract double[][][] getLogPitchLMProbs();

	public abstract double[][] getLogChordLMProbs();

	public abstract double[][][] getLogChordToPitchProbs();

}