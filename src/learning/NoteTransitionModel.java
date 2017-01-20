package learning;

import io.PitchEventIO;

import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;

import learning.EnvelopeActivationModel.State;
import main.Main;

public class NoteTransitionModel {
	
	public static interface TransScoresFactory {
		public double[][] getLogTransScores(int pitchId);
	}
	
    public static class SimpleTransScoresFactory implements TransScoresFactory {
    	public static final double ON_PROB = 1e-3;
		public double[][] getLogTransScores(int pitchId) {
			double[][] logTransScores = new double[State.values().length][State.values().length];
			logTransScores[State.OFF.ordinal()][State.OFF.ordinal()] = Math.log(1.0 - ON_PROB);
			logTransScores[State.OFF.ordinal()][State.ON.ordinal()] = Math.log(ON_PROB);
			logTransScores[State.ON.ordinal()][State.OFF.ordinal()] = 0.0;
			logTransScores[State.ON.ordinal()][State.ON.ordinal()] = Float.NEGATIVE_INFINITY;
			return logTransScores;
		}
    }
    
    public static class IMSLPTransScoresFactory implements TransScoresFactory {
        double[][][] logTransScores;

        public IMSLPTransScoresFactory() {
            try {
                ObjectInputStream is = new ObjectInputStream(new FileInputStream(new File(Main.imslpTransModelSer)));
                logTransScores = (double[][][])is.readObject();
                is.close();
                
                for (int pitchId=0; pitchId<PitchEventIO.N_MIDI_PITCH_IDS; ++pitchId) {
                	logTransScores[pitchId][State.ON.ordinal()][State.OFF.ordinal()] = 0.0;
                	logTransScores[pitchId][State.ON.ordinal()][State.ON.ordinal()] = Float.NEGATIVE_INFINITY;
                }
            } catch(Exception e) {
                System.err.println("Couldn't read transition scores...");
            }
        }

        public double[][] getLogTransScores(int pitch) {
            return logTransScores[pitch];
        }

    }
 
}
