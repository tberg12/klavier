package lm;

import tberg.murphy.arrays.a;
import io.PitchEventIO;
import learning.EnvelopeActivationModel;
import main.Main;
import tberg.murphy.tuple.Pair;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

import eval.PitchEventUtil;
import eval.PitchEventUtil.NoteState;

/**
 * @author jda
 */
public class IMSLPIndepTransitionReader {

    public static final int N_SONGS = 5000;

    public static final int OFF_STATE = EnvelopeActivationModel.State.OFF.ordinal();
    public static final int ON_STATE = EnvelopeActivationModel.State.ON.ordinal();
    public static final int N_STATES = EnvelopeActivationModel.State.values().length;

    protected List<NoteState[][]> noteStates;
    protected List<List<PitchEventIO.Event>> events;
    double[][][] transitionScores;
    //double[] widths;
    double[] widthMeans;
    double[] widthVars;

    public IMSLPIndepTransitionReader(File corpusRoot) {
        loadSongs(corpusRoot);

        //int[] lengthCounts = new int[1000];

        //for (boolean[][] m : activations) {
        //    int offset = 0;
        //    int counter = 0;
        //    //System.out.println(m.length);
        //    for (int pitch = 0; pitch < m[0].length; pitch++) {
        //        while (offset < m.length && m[offset][pitch]) {
        //            offset++; }
        //        while (offset + counter < m.length && !m[offset + counter][pitch]) {
        //            counter++; }
        //        //while (m[offset++][pitch]) { if (offset >= m.length) break; }
        //        //while (m[offset + counter++][pitch]) { if (offset + counter >= m.length) break; }
        //        counter /= 10;
        //        counter = Math.min(counter, lengthCounts.length-1);
        //        lengthCounts[counter]++;
        //        offset = offset + counter;
        //        counter = 0;
        //    }
        //}

        //System.out.println(a.toString(lengthCounts));
        //System.exit(1);

        transitionScores = computeTransitionScores();
        Pair<double[], double[]> meansAndVars = computeStats();
        widthMeans = meansAndVars.getFirst();
        widthVars = meansAndVars.getSecond();
    }

    protected void loadSongs(File corpusRoot) {
        assert corpusRoot.isDirectory();
        List<NoteState[][]> songs = new ArrayList<NoteState[][]>();
        List<List<PitchEventIO.Event>> allEvents = new ArrayList<List<PitchEventIO.Event>>();
        int counter = 0;
        for (File f : corpusRoot.listFiles()) {
            try {
                List<PitchEventIO.Event> events = PitchEventIO.readMIDIFile(f.getAbsolutePath());
                //boolean[][] activations = PitchEventIO.eventsToActivationMatrix(events, 0.01f);
                float secondsPerSample = (1.0f / 44000);
                float samplesPerFrame = Main.stftHopSize * Main.stftDownsampleSpectrogramHopSize;
                float secondsPerFrame = secondsPerSample * samplesPerFrame;
                NoteState[][] noteStates = PitchEventUtil.convertToRoll(events, secondsPerFrame);
                if (noteStates.length == 0) {
                    continue;
                }
                //System.out.println(activations.length);
                //System.out.println(activations[0].length);
                songs.add(noteStates);
                allEvents.add(events);
                System.out.println("ok song");
                System.out.println(f.getAbsolutePath());
                if (++counter >= N_SONGS) break;
            } catch(UnsupportedOperationException e) {
            } catch(Exception e) {
                e.printStackTrace();
            }
        }
        this.noteStates = songs;
        this.events = allEvents;
    }

    protected double[][][] computeTransitionScores() {
        double[][][] transitionScores = new double[PitchEventIO.N_MIDI_PITCH_IDS][N_STATES][N_STATES];
        double nFrames = 0;
        double[] nOnFrames = new double[transitionScores.length];
        double[] nOffFrames = new double[transitionScores.length];
        double[] onToOn = new double[transitionScores.length];
        double[] onToOff = new double[transitionScores.length];
        double[] offToOn = new double[transitionScores.length];
        double[] offToOff = new double[transitionScores.length];

        for (NoteState[][] noteStates : this.noteStates) {
            nFrames += noteStates.length;
            for (int frame = 1; frame < noteStates.length; frame++) {
                for (int pitchId = 0; pitchId < noteStates[frame].length; pitchId++) {

                    NoteState prevVal = noteStates[frame-1][pitchId];
                    NoteState currVal = noteStates[frame][pitchId];

                    if (prevVal != NoteState.OFF && currVal != NoteState.SUSTAIN) {
                        // this is an ON -> X event
                        nOnFrames[pitchId] += 1;
                        if (currVal == NoteState.ONSET) {
                            onToOn[pitchId] += 1;
                        } else if (currVal == NoteState.OFF) {
                            onToOff[pitchId] += 1;
                        }
                    } else if (prevVal == NoteState.OFF) {
                        // this is an OFF -> X event
                        nOffFrames[pitchId] += 1;
                        if (currVal == NoteState.ONSET) {
                            offToOn[pitchId] += 1;
                        } else if (currVal == NoteState.OFF) {
                            offToOff[pitchId] += 1;
                        } else {
                            assert false;
                        }
                    }

                    //boolean prev = activations[frame-1][pitchId];
                    //boolean curr = activations[frame][pitchId];
                    //if (prev) {
                    //    nOnFrames += 1;
                    //} else {
                    //    nOffFrames += 1;
                    //}
                    //if (!prev && !curr) {
                    //    transitionScores[pitchId][OFF_STATE][OFF_STATE] += 1;
                    //} else if (!prev && curr) {
                    //    transitionScores[pitchId][OFF_STATE][ON_STATE] += 1;
                    //} else if (prev && !curr) {
                    //    transitionScores[pitchId][ON_STATE][OFF_STATE] += 1;
                    //} else {
                    //    transitionScores[pitchId][ON_STATE][ON_STATE] += 1;
                    //}
                }
            }
        }

        // nOffFrames += 2;
        // nOnFrames += 2;

        // System.out.println(nOnFrames);
        // System.out.println(nOffFrames);

        //System.out.println(a.toString(transitionScores));

        //a.scalei(transitionScores, 1.0 / nFrames);
        //a.logi(transitionScores);
        //a.scalei(transitionScores, 1e-2);

        //double scale = 1e-1;
        double scale = 1;

        for (int i = 0; i < transitionScores.length; i++) {

            transitionScores[i][OFF_STATE][OFF_STATE] = Math.log(offToOff[i] / nOffFrames[i]);
            transitionScores[i][OFF_STATE][ON_STATE] = Math.log(offToOn[i] / nOffFrames[i]);
            transitionScores[i][ON_STATE][OFF_STATE] = Math.log(onToOff[i] / nOnFrames[i]);
            transitionScores[i][ON_STATE][ON_STATE] = Math.log(onToOn[i] / nOnFrames[i]);

            // // laplace smoothing
            // //a.addi(transitionScores[i], 1);
            // //a.scalei(transitionScores[i], 1.0 / nFrames);
            // a.scalei(transitionScores[i][OFF_STATE], 1.0 / nOffFrames);
            // a.scalei(transitionScores[i][ON_STATE], 1.0 / nOnFrames);
            // a.logi(transitionScores[i]);
            // //transitionScores[i][OFF_STATE][OFF_STATE] *= 1e-4;
            // transitionScores[i][OFF_STATE][OFF_STATE] *= scale;
            // transitionScores[i][OFF_STATE][ON_STATE] *= scale;
            // transitionScores[i][ON_STATE][OFF_STATE] *= scale;
            // transitionScores[i][ON_STATE][ON_STATE] *= scale;
            // //transitionScores[i][OFF_STATE][OFF_STATE] *= 1e-2;
        }
        //System.out.println(a.toString(transitionScores));

//        for (int pitchId = 0; pitchId < PitchEventIO.N_MIDI_PITCH_IDS; pitchId++) {
//            System.out.println(a.toString(transitionScores[pitchId]));
//        }

        return transitionScores;
    }

    protected Pair<double[], double[]> computeStats() {

        float secondsPerSample = (1.0f / 44000000);
        float samplesPerFrame = Main.stftHopSize * Main.stftDownsampleSpectrogramHopSize;
        float secondsPerFrame = secondsPerSample * samplesPerFrame;

        double[] totals = new double[PitchEventIO.N_MIDI_PITCH_IDS];
        double[] counts = new double[totals.length];
        for (List<PitchEventIO.Event> song : events) {
            for (PitchEventIO.Event event : song) {
                if (event.noteIndex < 0 || event.noteIndex >= totals.length) {
                    continue;
                }
                double width = event.offsetSec - event.onsetSec;
                totals[event.noteIndex] += width;
                counts[event.noteIndex]++;
            }
        }

        double[] means = a.pointwiseDiv(totals, counts);
        //System.out.println(a.toString(r));
        //System.out.println(a.toString(counts));

        double[] vars = new double[means.length];
        for (List<PitchEventIO.Event> song : events) {
            for (PitchEventIO.Event event : song) {
                if (event.noteIndex < 0 || event.noteIndex >= totals.length) {
                    continue;
                }
                double width = event.offsetSec - event.onsetSec;
                vars[event.noteIndex] += Math.pow(width - means[event.noteIndex], 2);
            }
        }

        vars = a.pointwiseDiv(vars, counts);

        return new Pair<double[],double[]>(means, vars);
    }

    public void save(String path) throws IOException {
        ObjectOutputStream tos = new ObjectOutputStream(new FileOutputStream(new File(path + "_transitions.ser")));
        tos.writeObject(transitionScores);
        tos.close();

        ObjectOutputStream wos = new ObjectOutputStream(new FileOutputStream(new File(path + "_widths.ser")));
        wos.writeObject(widthMeans);
        wos.writeObject(widthVars);
        wos.close();
    }

    public static void main(String[] args) throws IOException {
        new IMSLPIndepTransitionReader(new File("/Users/jda/Corpora/imslp")).save("/Users/jda/Code/music_transcription2/data/imslp_indep");
    }
}
