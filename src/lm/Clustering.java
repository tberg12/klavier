package lm;

import tberg.murphy.counter.Counter;
import tberg.murphy.counter.CounterInterface;
//import tberg.murphy.floatsequence.ForwardBackward;
//import tberg.murphy.floatsequence.ForwardBackward.StationaryLattice;
//import tberg.murphy.floatsequence.ForwardBackward.StationaryStateProjector;
import tberg.murphy.indexer.HashMapIndexer;
import tberg.murphy.indexer.Indexer;
import io.PitchEventIO;
import tberg.murphy.opt.*;
import tberg.murphy.sequence.ForwardBackward;
import tberg.murphy.tuple.Pair;
import tberg.murphy.arrays.a;

import java.io.File;
import java.io.IOException;
import java.util.*;

import eval.PitchEventUtil;
import eval.PitchEventUtil.NoteState;

class Chord {

    protected List<Integer> notes;

    public Chord(List<Integer> notes) {
        this.notes = notes;
    }

    public List<Integer> getNotes() {
        return Collections.unmodifiableList(notes);
    }

    public String toString() {
        return notes.toString();
    }

    @Override
    public int hashCode() {
        return notes.hashCode();
    }

    public boolean equals(Object o) {
        if (o.getClass() != Chord.class) {
            return false;
        }
        return notes.equals(((Chord)o).getNotes());
    }
}


// musical initialization

/**
 * @author jda
 */
public class Clustering {

    public static final int CLUSTERS = 24;
    public static final double REGULARIZER = 0.1f;
    //public static final float REGULARIZER = 0f;

    protected Indexer<String> featIndexer;
    protected Indexer<Chord> chordIndexer;

    protected int[][][][] emissionFeatIndices;
    //protected int[][][] emissionFeatIndices;
    protected int[][][] transitionFeatIndices;
    protected double[][][][] emissionFeatValues;
    //protected double]][][] emissionFeatValues;
    protected double[][][] transitionFeatValues;

    protected int[][][] chordFeatIndices;
    protected double[][][] chordFeatValues;

    protected Set<Chord> knownChords;

    public Clustering(File corpusRoot) {
        knownChords = new HashSet<Chord>();
        List<List<Chord>> sequences = loadSongs(corpusRoot);
        featIndexer = new HashMapIndexer<String>();
        chordIndexer = new HashMapIndexer<Chord>();
        System.out.println("loaded activations");

        Pair<int[][][], double[][][]> chordFeatures = computeChordFeatures();
        chordFeatIndices = chordFeatures.getFirst();
        chordFeatValues = chordFeatures.getSecond();
        System.out.println("computed chord features");

        Pair<int[][][][], double[][][][]> emissionFeatures = computeEmissionFeatures(sequences);
        emissionFeatIndices = emissionFeatures.getFirst();
        emissionFeatValues = emissionFeatures.getSecond();
        System.out.println("computed emission features");
        //System.out.println(sequences);
        //System.out.println(emissionFeatIndices.length);
        //System.out.println(emissionFeatIndices[0].length);
        //System.out.println(emissionFeatIndices[0][0].length);
        //System.out.println(emissionFeatIndices[0][0][0].length);

        Pair<int[][][], double[][][]> transitionFeatures = computeTransitionFeatures(sequences);
        transitionFeatIndices = transitionFeatures.getFirst();
        transitionFeatValues = transitionFeatures.getSecond();
        System.out.println("computed transition features");
    }

    protected List<List<Chord>> loadSongs(File corpusRoot) {
       assert corpusRoot.isDirectory();
       List<List<Chord>> songs = new ArrayList<List<Chord>>();
       int counter = 0;
       for (File f : corpusRoot.listFiles()) {
       //{
           //File f = new File("/Users/jda/Corpora/imslp/IMSLP00001#1#Smartscore-10.2.1.mid");
           try {
               List<PitchEventIO.Event> events = PitchEventIO.readMIDIFile(f.getAbsolutePath());
               NoteState[][] noteStates = PitchEventUtil.convertToRoll(events, 0.01f);
               List<Chord> song = new ArrayList<Chord>();
               for (NoteState[] frame : noteStates) {
                   ArrayList<Integer> notes = new ArrayList<Integer>();
                   for (int pitch = 0; pitch < frame.length; pitch++) {
                       if (frame[pitch] != NoteState.OFF) {
                           int chroma = pitch % 12;
                           if (!notes.contains(chroma)) {
                               notes.add(chroma);
                           }
                           //notes.add(pitch);
                       }
                   }
                   Chord chord = new Chord(notes);
                   if (song.size() > 0 && song.get(song.size()-1).equals(chord)) {
                       continue;
                   }
                   song.add(chord);
                   knownChords.add(chord);
               }
               if (song.size() == 0) {
                   continue;
               }
               songs.add(song);
               System.out.println("ok song");
               if (++counter >= 10) break;
           } catch(UnsupportedOperationException e) {
           //} catch(IOException e) {
           } catch(Exception e) {
               e.printStackTrace();
               //assert false;
           }
       }
       return songs;
    }

    protected Pair<int[][][][], double[][][][]> computeEmissionFeatures(List<List<Chord>> sequences) {
        int[][][][] indices = new int[sequences.size()][][][];
        double[][][][] values = new double[sequences.size()][][][];
        for (int d = 0; d < sequences.size(); d++) {
            List<Chord> song = sequences.get(d);
            indices[d] = new int[song.size()][][];
            values[d] = new double[song.size()][][];
            for (int t = 0; t < song.size(); t++) {
                Chord chord = song.get(t);
                int c = chordIndexer.getIndex(chord);
                indices[d][t] = new int[CLUSTERS][];
                values[d][t] = new double[CLUSTERS][];
                for (int state = 0; state < CLUSTERS; state++) {
                    indices[d][t][state] = chordFeatIndices[c][state];
                    values[d][t][state] = chordFeatValues[c][state];
                //    List<String> emissionFeatures = makeEmissionFeatures(chord, state);
                //    indices[d][t][state] = new int[emissionFeatures.size()];
                //    values[d][t][state] = new double]emissionFeatures.size()];
                //    for (int i = 0; i < emissionFeatures.size(); i++) {
                //        String feat = emissionFeatures.get(i);
                //        int featIndex = featIndexer.getIndex(feat);
                //        indices[d][t][state][i] = featIndex;
                //        values[d][t][state][i] = 1f;
                //    }
                }
            }
        }
        return new Pair<int[][][][], double[][][][]>(indices, values);
    }

    protected Pair<int[][][], double[][][]> computeChordFeatures() {
        int[][][] indices = new int[knownChords.size()][][];
        double[][][] values = new double[knownChords.size()][][];
        for (Chord chord : knownChords) {
            int c = chordIndexer.getIndex(chord);
            indices[c] = new int[CLUSTERS][];
            values[c] = new double[CLUSTERS][];
            for (int s = 0; s < CLUSTERS; s++) {
                List<String> emissionFeatures = makeEmissionFeatures(chord, s);
                indices[c][s] = new int[emissionFeatures.size()];
                values[c][s] = new double[emissionFeatures.size()];
                for (int i = 0; i < emissionFeatures.size(); i++) {
                    String feat = emissionFeatures.get(i);
                    int featIndex = featIndexer.getIndex(feat);
                    indices[c][s][i] = featIndex;
                    values[c][s][i] = 1f;
                }
            }
        }
        return new Pair<int[][][], double[][][]>(indices, values);
    }

    protected List<String> makeEmissionFeatures(Chord chord, int state) {
        ArrayList<String> r = new ArrayList<String>(chord.getNotes().size());
        for (int note : chord.getNotes()) {
            r.add("EMIT_" + state + "_" + note);
        }
        return r;
    }

    protected Pair<int[][][], double[][][]> computeTransitionFeatures(List<List<Chord>> sequences) {
        int[][][] indices = new int[CLUSTERS][][];
        double[][][] values = new double[CLUSTERS][][];
        for (int s1 = 0; s1 < CLUSTERS; s1++) {
            indices[s1] = new int[CLUSTERS][];
            values[s1] = new double[CLUSTERS][];
            for (int s2 = 0; s2 < CLUSTERS; s2++) {
                List<String> transitionFeatures = makeTransitionFeatures(s1, s2);
                indices[s1][s2] = new int[transitionFeatures.size()];
                values[s1][s2] = new double[transitionFeatures.size()];
                for (int i = 0; i < transitionFeatures.size(); i++) {
                    String feat = transitionFeatures.get(i);
                    int featIndex = featIndexer.getIndex(feat);
                    indices[s1][s2][i] = featIndex;
                    values[s1][s2][i] = 1f;
                }
            }
        }
        return new Pair<int[][][], double[][][]>(indices, values);
    }

    protected List<String> makeTransitionFeatures(int s1, int s2) {
        return Arrays.asList("TRANS_" + s1 + "_" + s2);
    }

    protected double dot(int[] featIndices, double[] featValues, double[] weights) {
        double r = 0;
        for (int i = 0; i < featIndices.length; i++) {
            r += featValues[i] * weights[featIndices[i]];
        }
        return r;
    }

    protected double norm(double[] weights) {
        double r = 0;
        for (double f : weights) {
            r += f * f;
        }
        return r;
    }

    public void train() {
        double[] weights = a.randDouble(featIndexer.size(), new Random());
        a.addi(weights, -0.5);
        for (int i = 0; i < 12; i++) {
            for (int j = 0; j < 12; j++) {
                String majFeatName = "EMIT_" + i + "_" + j;
                String minFeatName = "EMIT_" + (i+12) + "_" + j;
                int majFeatIndex = featIndexer.getIndex(majFeatName);
                int minFeatIndex = featIndexer.getIndex(minFeatName);
                weights[majFeatIndex] = 0;
                weights[minFeatIndex] = 0;

                if (j == i) {
                    weights[majFeatIndex] = 1;
                    weights[minFeatIndex] = 1;
                }

                if (j == (i + 3) % 12 ) {
                    weights[minFeatIndex] = 1;
                }

                if (j == (i + 4) % 12) {
                    weights[majFeatIndex] = 1;
                }
            }
        }
        //a.scalei(weights, 0.01f);
        System.out.println("started training");

        for (int t = 0; t < 100; t++) {

            // System.out.println("EM iter " + t);

            weights = doEMIter(weights);

            //for (int i = 0; i < featIndexer.size(); i++) {
            //    System.out.println(featIndexer.getObject(i) + ":\t" + weights[i]);
            //}
            //System.out.println();
        }

        System.out.println("EMISSIONS");

        for (int c = 0; c < CLUSTERS; c++) {
            for (int p = 0; p < 12; p++) {
                System.out.printf("%.3f", weights[featIndexer.getIndex("EMIT_" + c + "_" + p)]);
                System.out.print("\t");
            }
            System.out.println();
        }

        System.out.println();
        System.out.println("TRANSITIONS");

        for (int s1 = 0; s1 < CLUSTERS; s1++) {
            for (int s2 = 0; s2 < CLUSTERS; s2++) {
                System.out.printf("%.3f", weights[featIndexer.getIndex("TRANS_" + s1 + "_" + s2)]);
                System.out.print("\t");
            }
            System.out.println();
        }

        sample(weights, 50, 3);

        //DifferentiableFunction objective = makeDirectObjective();

        //EmpiricalGradientTester.test(objective, a.toDouble(weights), 0.01, 0.01, 0.001);
        //System.exit(1);

        //Minimizer minimizer = new LBFGSMinimizer(0.01, 10);
        //double[] finalWeights = minimizer.minimize(objective, a.toDouble(weights), true, null);
    }

    protected void sample(double[] weights, int chords, int notesPerChord) {
        Random r = new Random();
        int cluster = r.nextInt(CLUSTERS);

        double[][] chordScores = computeChordScores(weights);
        double[][] transitionScores = computeTransitionScores(weights).getFirst();

        for (int t = 0; t < chords; t++) {
            System.out.println(cluster);
            //int bestChord = -1;
            //double bestScore = Double.MIN_VALUE;
            //System.out.println(chordScores.length);
            //for (int c = 0; c < chordScores.length; c++) {
            //    double score = chordScores[c][cluster];
            //    System.out.println("> " + score);
            //    if (score > bestScore) {
            //        bestChord = c;
            //        bestScore = score;
            //    }
            //}
            //System.out.print(cluster);
            //System.out.print(" ");
            //System.out.println(chordIndexer.getObject(bestChord));
            double[] noteScores = new double[12];
            for (int i = 0; i < 12; i++) {
                noteScores[i] = Math.exp(weights[featIndexer.getIndex("EMIT_" + cluster + "_" + i)]);
            }
            for (int i = 0; i < notesPerChord; i++) {
                System.out.print(categorialSample(noteScores));
                System.out.print(" ");
            }
            System.out.println();

            cluster = categorialSample(a.exp(transitionScores[cluster]));
        }
    }

    protected int categorialSample(double[] scores) {
        double[] mScores = new double[scores.length];
        System.arraycopy(scores, 0, mScores, 0, scores.length);
        Arrays.sort(mScores);
        for (int i = 0; i < mScores.length; i++) {
            mScores[i] /= a.sum(mScores);
        }
        for (int i = 1; i < mScores.length; i++) {
            mScores[i] += mScores[i-1];
        }
        double a = new Random().nextDouble();
        for (int i = 0; i < mScores.length; i++) {
            if (mScores[i] > a) {
                return i;
            }
        }
        assert false;
        return -1;
    }

    protected double[] doEMIter(double[] weights) {
        long precomp_start = System.currentTimeMillis();
        Pair<double[][][], double[]> emissions = computeEmissionScores(weights);
        double[][][] emissionScores = emissions.getFirst();
        double[] emissionNormalizers = emissions.getSecond();
        double[] emissionLogNormalizers = a.log(emissionNormalizers);
        double[][] transitionScores = computeTransitionScores(weights).getFirst();
        long precomp_total = System.currentTimeMillis() - precomp_start;
        //System.out.println("initial score computation took " + precomp_total);

        //System.out.println(a.toString(emissionScores));
        // System.out.println("emission normalizers");
        // System.out.println(a.toString(emissionNormalizers));
        // System.out.println("transition scores");
        // System.out.println(a.toString(transitionScores));

        ForwardBackward.StationaryLattice lattice = makeLattice(emissionScores, emissionLogNormalizers, transitionScores);
        ForwardBackward.StationaryStateProjector projector = makeStateProjector();


        long fbstart = System.currentTimeMillis();
        Pair<ForwardBackward.NodeMarginals, ForwardBackward.StationaryEdgeMarginals> marginals = ForwardBackward.computeMarginalsLogSpace(lattice, projector, false, 4);
        ForwardBackward.NodeMarginals nodeMarginals = marginals.getFirst();
        ForwardBackward.StationaryEdgeMarginals edgeMarginals = marginals.getSecond();
        long fbtotal = System.currentTimeMillis() - fbstart;
        //System.out.println("F-B took " + fbtotal);

        // System.out.println(a.toString(edgeMarginals.allowedForwardEdgesExpectedCounts(0, 0)));

        System.out.println("L = " + edgeMarginals.logMarginalProb());
        //System.out.println("L = " + nodeMarginals.logMarginalProb());

        final DifferentiableFunction objective = makeObjective(nodeMarginals, edgeMarginals);
        Minimizer minimizer = new LBFGSMinimizer(1e-5, 100);
        //Minimizer minimizer = new LBFGSFortranWrapperMinimizer(0.01, 100);

        //EmpiricalGradientTester.test(objective, weights, 1e-3, 1e-1, 1e-3);
        //System.exit(1);

        //objective.calculate(a.toDouble(weights));
        //double[] newWeights = new double[3];
        //System.exit(1);

        //System.out.println("started minimizer");
        double[] newWeights = minimizer.minimize(objective, weights, false, new Minimizer.Callback() {

            @Override
            public void callback(double[] guess, int iter, double val, double[] grad) {
                //System.out.println("L = " + val);
                //System.out.println(a.sum(a.pow(objective.calculate(guess).getSecond(), 2)));
                //System.out.println(a.toString(objective.calculate(guess).getSecond()));
            }
        });

        //return a.toDouble(newWeights);
        return newWeights;
    }

    protected Pair<double[][][], double[]> computeEmissionScores(double[] weights) {
        double[][][] emissionScores = new double[emissionFeatIndices.length][][];
        for (int d = 0; d < emissionFeatIndices.length; d++) {
            emissionScores[d] = new double[emissionFeatIndices[d].length][];
            for (int t = 0; t < emissionFeatIndices[d].length; t++) {
                emissionScores[d][t] = new double[CLUSTERS];
                for (int s = 0; s < CLUSTERS; s++) {
                    emissionScores[d][t][s] = dot(emissionFeatIndices[d][t][s], emissionFeatValues[d][t][s], weights);
                }
            }
        }

        double[] emissionNormalizers = new double[CLUSTERS];
        for (int c = 0; c < chordFeatIndices.length; c++) {
            for (int s = 0; s < CLUSTERS; s++) {
                emissionNormalizers[s] += Math.exp(dot(chordFeatIndices[c][s], chordFeatValues[c][s], weights));
            }
        }

        return new Pair<double[][][], double[]>(emissionScores, emissionNormalizers);
    }

    protected double[][] computeChordScores(double[] weights) {
        double[][] chordScores = new double[chordFeatIndices.length][CLUSTERS];
        for (int c = 0; c < chordFeatIndices.length; c++) {
            for (int s = 0; s < CLUSTERS; s++) {
                chordScores[c][s] += dot(chordFeatIndices[c][s], chordFeatValues[c][s], weights);
            }
        }
        return chordScores;
    }

    protected Pair<double[][], double[]> computeTransitionScores(double[] weights) {
        double[][] transitionScores = new double[CLUSTERS][CLUSTERS];
        for (int s1 = 0; s1 < CLUSTERS; s1++) {
            for (int s2 = 0; s2 < CLUSTERS; s2++) {
                transitionScores[s1][s2] = dot(transitionFeatIndices[s1][s2], transitionFeatValues[s1][s2], weights);
            }
        }

        double[] transitionNormalizers = new double[CLUSTERS];
        for (int s1 = 0; s1 < CLUSTERS; s1++) {
            for (int s2 = 0; s2 < CLUSTERS; s2++) {
                transitionNormalizers[s1] += Math.exp(dot(transitionFeatIndices[s1][s2], transitionFeatValues[s1][s2], weights));
            }
        }
        return new Pair<double[][], double[]>(transitionScores, transitionNormalizers);
    }

    protected ForwardBackward.StationaryLattice makeLattice(final double[][][] emissionScores, final double[] emissionLogNormalizers, final double[][] transitionScores) {
        return new ForwardBackward.StationaryLattice() {

            @Override
            public int numSequences() {
                return emissionFeatIndices.length;
            }

            @Override
            public int sequenceLength(int d) {
                return emissionFeatIndices[d].length;
            }

            @Override
            public int numStates(int d) {
                return CLUSTERS;
            }

            @Override
            public double nodeLogPotential(int d, int t, int s) {
                double r = emissionScores[d][t][s] - emissionLogNormalizers[s]; //(float)Math.log(emissionNormalizers[s]);
                assert r < 0;
                return r;
            }

            @Override
            public double nodePotential(int d, int t, int s) {
                return (double)Math.exp(nodeLogPotential(d, t, s));
            }

            @Override
            public double[] allowedEdgesLogPotentials(int d, int s, boolean backward) {
                double[] r = new double[CLUSTERS];
                double norm = 0;
                for (int s2 = 0; s2 < CLUSTERS; s2++) {
                    if (!backward) {
                        r[s2] = transitionScores[s][s2];
                        // TODO this can be precomputed
                        //norm += Math.exp(r[s2]);
                    } else {
                        r[s2] = transitionScores[s2][s];
                        //norm += Math.exp(r[s2]);
                    }
                }

                if (!backward) {
                    norm = 0;
                    for (int s2 = 0; s2 < CLUSTERS; s2++) {
                        norm += Math.exp(transitionScores[s][s2]);
                    }
                    a.addi(r, -Math.log(norm));
                } else {
                    for (int fromState = 0; fromState < CLUSTERS; fromState++) {
                        norm = 0;
                        for (int toOtherState = 0; toOtherState < CLUSTERS; toOtherState++) {
                            norm += Math.exp(transitionScores[fromState][toOtherState]);
                        }
                        r[fromState] -= Math.log(norm);
                    }
                }

                //a.addi(r,-Math.log(norm));
                for (int i = 0 ; i < CLUSTERS; i++) {
                    assert r[i] < 0;
                }
                return r;
            }

            @Override
            public double[] allowedEdgesPotentials(int d, int s, boolean backward) {
                return a.exp(allowedEdgesLogPotentials(d, s, backward));
            }

            private int[] ae = a.enumerate(0, CLUSTERS);
            @Override
            public int[] allowedEdges(int d, int s, boolean backward) {
                return ae;
            }
        };
    }

    protected ForwardBackward.StationaryStateProjector makeStateProjector() {
        return new ForwardBackward.StationaryStateProjector() {

            @Override
            public int domainSize(int d, int t) {
                return CLUSTERS;
            }

            @Override
            public int rangeSize(int d) {
                return CLUSTERS;
            }

            @Override
            public int project(int d, int t, int s) {
                return s;
            }
        };
    }

    protected DifferentiableFunction makeObjective(final ForwardBackward.NodeMarginals nodeMarginals,
                                                   final ForwardBackward.StationaryEdgeMarginals edgeMarginals
                                                   //final double[][][] emissionScores,
                                                   //final double[] emissionNormalizers,
                                                   //final double[][] transitionScores
                                                   ) {
        return new DifferentiableFunction() {

            @Override
            public Pair<Double, double[]> calculate(double[] x) {

                Pair<double[][][],double[]> emissions = computeEmissionScores(x);
                double[][][] emissionScores = emissions.getFirst();
                double[] emissionNormalizers = emissions.getSecond();
                double[] emissionLogNormalizers = a.log(emissionNormalizers);

                Pair<double[][],double[]> transitions = computeTransitionScores(x);
                double[][] transitionScores = transitions.getFirst();
                double[] transitionNormalizers = transitions.getSecond();
                double[] transitionLogNormalizers = a.log(transitionNormalizers);

                double[][] chordScores = computeChordScores(x);
                double[][] expectedChordFeatures = computeExpectedChordFeatures(x, chordScores, emissionNormalizers);

                double[][] expectedTransitionFeatures = computeExpectedTransitionFeatures(x, transitionScores, transitionNormalizers);

                double nll = 0;
                double[] gradient = a.zerosDouble(x.length);

                nll -= .5 * REGULARIZER * a.innerProd(x, x);
                a.combi(gradient, 1, x, -REGULARIZER);

                for (int d = 0; d < emissionFeatIndices.length; d++) {
                    for (int t = 0; t < emissionFeatIndices[d].length; t++) {
                        for (int s = 0; s < CLUSTERS; s++) {
                            nll += nodeMarginals.nodeCondProbs(d, t)[s] * (emissionScores[d][t][s] - emissionLogNormalizers[s]);
                            //nll -= nodeMarginals.nodeCondProbs(d, t)[s] * Math.log(nodeMarginals.nodeCondProbs(d,t)[s]);
                            a.combi(gradient, 1, deltaEmission(d, t, s, expectedChordFeatures[s]), nodeMarginals.nodeCondProbs(d, t)[s]);

                            //double[] de = deltaEmission(d, t, s, expectedChordFeatures[s]);
                            //for(int i = 0; i < de.length; i++) {
                            //    gradient[i] += nodeMarginals.nodeCondProbs(d,t)[s] * de[i];
                            //}

                            // nll += nodeMarginals.nodeCondProbs(d, t)[s] * emissionScores[d][t][s];
                            // a.combi(gradient, 1, deltaEmission(d, t, s, expectedChordFeatures[s]), nodeMarginals.nodeCondProbs(d, t)[s]);
                            // nll += nodeMarginals.nodeCondProbs(d, t)[s] * (-Math.log(emissionNormalizers[s]));
                            // a.combi(gradient, 1, deltaEmission(d, t, s, expectedChordFeatures[s]), nodeMarginals.nodeCondProbs(d, t)[s]);
                        }
                    }
                }

                for (int d = 0; d < emissionFeatIndices.length; d++) {
                    for (int s1 = 0; s1 < CLUSTERS; s1++) {
                        for (int s2 = 0; s2 < CLUSTERS; s2++) {
                            nll += edgeMarginals.allowedForwardEdgesExpectedCounts(d, s1)[s2] * (transitionScores[s1][s2] - transitionLogNormalizers[s1]); // Math.log(transitionNormalizers[s1]));
                            //nll -= edgeMarginals.allowedForwardEdgesExpectedCounts(d, s1)[s2] * Math.log(edgeMarginals.allowedForwardEdgesExpectedCounts(d,s1)[s2]);
                            a.combi(gradient, 1, deltaTransition(s1, s2, expectedTransitionFeatures[s1]), edgeMarginals.allowedForwardEdgesExpectedCounts(d, s1)[s2]);
                        }
                    }
                }

                nll *= -1;
                a.scalei(gradient, -1);

                //System.out.println("L = " + nll);
                return new Pair<Double, double[]>(nll, gradient);
            }
        };
    }

    protected double[][] computeExpectedChordFeatures(double[] weights, double[][] chordScores, double[] emissionNormalizers) {
        double[][] r = new double[CLUSTERS][weights.length];
        for (int s = 0; s < CLUSTERS; s++) {
            for (int c = 0; c < chordScores.length; c++) {
                for (int i = 0; i < chordFeatIndices[c][s].length; i++) {
                    r[s][chordFeatIndices[c][s][i]] += Math.exp(chordScores[c][s]) * chordFeatValues[c][s][i];
                    //System.out.println(Math.exp(chordScores[c][s]));
                    //System.out.println(chordFeatIndices.length);
                    //System.out.println(chordScores.length);
                    //System.out.println(c);
                    //System.out.println();
                }
            }
            //System.out.println(emissionNormalizers[s]);
            a.scalei(r[s], 1/emissionNormalizers[s]);
        }
        return r;
    }

    protected double[][] computeExpectedTransitionFeatures(double[] weights, double[][] transitionScores, double[] transitionNormalizers) {
        double[][] r = new double[CLUSTERS][weights.length];
        for (int s1 = 0; s1 < CLUSTERS; s1++) {
            for (int s2 = 0; s2 < CLUSTERS; s2++) {
                for (int i = 0; i < transitionFeatIndices[s1][s2].length; i++) {
                    r[s1][transitionFeatIndices[s1][s2][i]] += Math.exp(transitionScores[s1][s2]) * transitionFeatValues[s1][s2][i];
                }
            }
            a.scalei(r[s1], 1/transitionNormalizers[s1]);
        }
        return r;
    }

    protected double[] deltaEmission(int d, int t, int s, double[] expectedChordFeatures) {
        double[] r = a.zerosDouble(expectedChordFeatures.length);
        for (int i = 0; i < emissionFeatIndices[d][t][s].length; i++) {
            r[emissionFeatIndices[d][t][s][i]] += emissionFeatValues[d][t][s][i];
        }
        a.combi(r, 1, expectedChordFeatures, -1);
        return r;
    }

    protected double[] deltaTransition(int s1, int s2, double[] expectedTransitionFeatures) {
        double[] r = new double[expectedTransitionFeatures.length];
        for (int i = 0; i < transitionFeatIndices[s1][s2].length; i++) {
            r[transitionFeatIndices[s1][s2][i]] += transitionFeatValues[s1][s2][i];
        }
        a.combi(r, 1, expectedTransitionFeatures, -1);
        return r;
    }

    protected DifferentiableFunction makeDirectObjective() {
        return new DifferentiableFunction() {
            @Override
            public Pair<Double,double[]> calculate(double[] weights) {
                //double[] weights = a.toDouble(dweights);
                Pair<double[][][], double[]> emissions = computeEmissionScores(weights);
                double[][][] emissionScores = emissions.getFirst();
                double[] emissionNormalizers = emissions.getSecond();
                double[] emissionLogNormalizers = a.log(emissionNormalizers);
                //float[][] transitionScores = computeTransitionScores(weights).getFirst();
                Pair<double[][],double[]> transitions = computeTransitionScores(weights);
                double[][] transitionScores = transitions.getFirst();
                double[] transitionNormalizers = transitions.getSecond();
                double[] transitionLogNormalizers = a.log(transitionNormalizers);

                ForwardBackward.StationaryLattice lattice = makeLattice(emissionScores, emissionLogNormalizers, transitionScores);
                ForwardBackward.StationaryStateProjector projector = makeStateProjector();

                Pair<ForwardBackward.NodeMarginals, ForwardBackward.StationaryEdgeMarginals> marginals = ForwardBackward.computeMarginalsLogSpace(lattice, projector, false, 4);
                ForwardBackward.NodeMarginals nodeMarginals = marginals.getFirst();
                ForwardBackward.StationaryEdgeMarginals edgeMarginals = marginals.getSecond();

                double[][] chordScores = computeChordScores(weights);
                double[][] expectedChordFeatures = computeExpectedChordFeatures(weights, chordScores, emissionNormalizers);
                double[][] expectedTransitionFeatures = computeExpectedTransitionFeatures(weights, transitionScores, transitionNormalizers);

                double[] gradient = a.zerosDouble(weights.length);

                //a.combi(gradient, 1, dweights, -REGULARIZER);

                for (int d = 0; d < emissionFeatIndices.length; d++) {
                    for (int t = 0; t < emissionFeatIndices[d].length; t++) {
                        for (int s = 0; s < CLUSTERS; s++) {
                            a.combi(gradient, 1, deltaEmission(d, t, s, expectedChordFeatures[s]), nodeMarginals.nodeCondProbs(d, t)[s]);
                        }
                    }
                }

                for (int d = 0; d < emissionFeatIndices.length; d++) {
                    for (int s1 = 0; s1 < CLUSTERS; s1++) {
                        for (int s2 = 0; s2 < CLUSTERS; s2++) {
                            a.combi(gradient, 1, deltaTransition(s1, s2, expectedTransitionFeatures[s1]), edgeMarginals.allowedForwardEdgesExpectedCounts(d, s1)[s2]);
                        }
                    }
                }

                return new Pair<Double,double[]>((double)nodeMarginals.logMarginalProb(), gradient);
            }
        };
    }

    /*
    public void train() {
        final float[] weights = a.zerosFloat(featIndexer.size());
        //float[] init = a.randFloat(featIndexer.size(), new Random());
        float[] init = a.onesFloat(featIndexer.size());
        final float[] emissionNormalizers = new float[CLUSTERS];

        final StationaryLattice lattice = new StationaryLattice() {

            @Override
            public int numSequences() {
                return emissionFeatIndices.length;
            }

            @Override
            public int sequenceLength(int d) {
                return emissionFeatIndices[d].length;
            }

            @Override
            public int numStates(int d) {
                return CLUSTERS;
            }

            @Override
            public float nodeLogPotential(int d, int t, int s) {
                float r = dot(emissionFeatIndices[d][t][s], emissionFeatValues[d][t][s], weights) - (float)Math.log(emissionNormalizers[s]);
                //System.out.println(r);
                return r;
            }

            @Override
            public float[] allowedEdgesLogPotentials(int d, int s, boolean backward) {
                float[] r = new float[CLUSTERS];
                float norm = 0;
                for (int s2 = 0; s2 < CLUSTERS; s2++) {
                    if (!backward) {
                        r[s2] = dot(transitionFeatIndices[s][s2], transitionFeatValues[s][s2], weights);
                    } else {
                        r[s2] = dot(transitionFeatIndices[s2][s], transitionFeatValues[s2][s], weights);
                    }
                    norm += Math.exp(r[s2]);
                }
                for (int s2 = 0; s2 < CLUSTERS; s2++) {
                    r[s2] -= Math.log(norm);
                }
                return r;
            }

            @Override
            public float nodePotential(int d, int t, int s) {
                return (float)Math.exp(nodeLogPotential(d,t,s));
            }

            @Override
            public float[] allowedEdgesPotentials(int d, int s, boolean backward) {
                float[] r = allowedEdgesLogPotentials(d, s, backward);
                for (int i = 0; i < r.length; i++) {
                    r[i] = (float)Math.exp(r[i]);
                }
                return r;
            }

            @Override
            public int[] allowedEdges(int d, int s, boolean backward) {
                return a.enumerate(0, CLUSTERS);
            }
        };

        final StationaryStateProjector projector = new StationaryStateProjector() {

            @Override
            public int domainSize(int d, int t) {
                return CLUSTERS;
            }

            @Override
            public int rangeSize(int d) {
                return CLUSTERS;
            }

            @Override
            public int project(int d, int t, int s) {
                return s;
            }
        };

        //System.out.println(weights.length);
        //System.out.println(init.length);
        System.arraycopy(init, 0, weights, 0, weights.length);
        //System.out.println(Arrays.toString(weights));
        for (int state = 0; state < CLUSTERS; state++) {
            emissionNormalizers[state] = 0;
            for (int c = 0; c < chordFeatIndices.length; c++) {
                emissionNormalizers[state] += dot(chordFeatIndices[c][state], chordFeatValues[c][state], weights);
            }
        }

        for (int iter = 0; iter < 10; iter++) {

            Pair<ForwardBackward.NodeMarginals, ForwardBackward.StationaryEdgeMarginals> marginals = ForwardBackward.computeMarginalsLogSpace(lattice, projector, false, 1);
            final ForwardBackward.NodeMarginals nodeMarginals = marginals.getFirst();
            final ForwardBackward.StationaryEdgeMarginals edgeMarginals = marginals.getSecond();

            DifferentiableFunction objective = new DifferentiableFunction() {

                @Override
                public Pair<Double, double[]> calculate(double[] x) {

                    // System.out.println(x.length);
                    // System.out.println(chordFeatIndices.length);
                    // System.out.println(emissionFeatIndices.length);
                    // System.out.println(emissionFeatIndices[0].length);

                    for (int i = 0; i < x.length; i++) {
                        weights[i] = (float)x[i];
                    }

                    for (int state = 0; state < CLUSTERS; state++) {
                        emissionNormalizers[state] = 0;
                        for (int c = 0; c < chordFeatIndices.length; c++) {
                            emissionNormalizers[state] += Math.exp(dot(chordFeatIndices[c][state], chordFeatValues[c][state], weights));
                        }
                    }

                    //float logLikelihood = nodeMarginals.logMarginalProb() - 0.5f * REGULARIZER * norm(weights);
                    float logLikelihood = -0.5f * REGULARIZER * norm(weights);
                    float[] gradient = new float[weights.length];
                    scaledAddI(gradient, -REGULARIZER, weights);

                    System.out.println();
                    System.out.println(logLikelihood);
                    System.out.println(Arrays.toString(gradient));

                    float[][] negExpectedFeatValues = new float[CLUSTERS][weights.length];
                    for (int state = 0; state < CLUSTERS; state++) {
                        for (int c = 0; c < chordFeatIndices.length; c++) {
                            float chordScore = (float)Math.exp(dot(chordFeatIndices[c][state], chordFeatValues[c][state], weights));
                            for (int i = 0; i < chordFeatIndices[c][state].length; i++) {
                                negExpectedFeatValues[state][chordFeatIndices[c][state][i]] -= Math.exp(chordScore * chordFeatValues[c][state][i]);
                            }
                        }
                        for (int i = 0; i < negExpectedFeatValues.length; i++) {
                            negExpectedFeatValues[state][i] /= emissionNormalizers[state];
                        }
                    }

                    for (int d = 0; d < emissionFeatIndices.length; d++) {
                        for (int t = 0; t < emissionFeatIndices[d].length; t++) {
                            //System.out.println(".." + t);
                            for (int state = 0; state < CLUSTERS; state++) {
                                //for (Chord chord : knownChords) {
                                //    scaledAddI(gradient, nodeMarginals.nodeCondProbs(d, t)[state], deltaEmission(d,t,state,chord,weights,emissionNormalizers[state]));
                                //}
                                logLikelihood += nodeMarginals.nodeCondProbs(d,t)[state] * (dot(emissionFeatIndices[d][t][state], emissionFeatValues[d][t][state], weights) - Math.log(emissionNormalizers[state]));
                                System.out.println(nodeMarginals.nodeCondProbs(d,t)[state]);
                                System.out.println(dot(emissionFeatIndices[d][t][state], emissionFeatValues[d][t][state], weights));

                                //System.out.println(emissionNormalizers[state]);
                                //System.exit(1);
                                scaledAddI(gradient, nodeMarginals.nodeCondProbs(d,t)[state], deltaEmission(d,t,state,weights,negExpectedFeatValues[state]));
                                System.exit(1);
                            }
                        }
                    }
                    //System.exit(1);

                    // for (int d = 0; d < emissionFeatIndices.length; d++) {
                    //     //for (int t = 0; t < emissionFeatIndices[t].length; t++) {
                    //         for (int state1 = 0; state1 < CLUSTERS; state1++) {
                    //             for (int state2 = 0; state2 < CLUSTERS; state2++) {
                    //                 float logLikelihoodNorm = 0;
                    //                 for (int ostate2 = 0; ostate2 < CLUSTERS; ostate2++) {
                    //                     logLikelihoodNorm += Math.exp(dot(transitionFeatIndices[state1][ostate2], transitionFeatValues[state1][ostate2], weights));
                    //                 }
                    //                 float logLikelihoodTerm = edgeMarginals.allowedForwardEdgesExpectedCounts(d, state1)[state2] * (dot(transitionFeatIndices[state1][state2], transitionFeatValues[state1][state2], weights) - (float)Math.log(logLikelihoodNorm));
                    //                 logLikelihood += logLikelihoodTerm;

                    //                 scaledAddI(gradient, edgeMarginals.allowedForwardEdgesExpectedCounts(d, state1)[state2], deltaTransition(state1, state2, weights));
                    //             }
                    //         }
                    //     //}
                    // }

                    a.scalei(gradient,-1);
                    logLikelihood *= -1;
                    System.out.println(logLikelihood);

                    return new Pair<Double,double[]>((double)logLikelihood, a.toDouble(gradient));
                }
            };

            EmpiricalGradientTester.test(objective, a.toDouble(init), 0.01, 0.01, 0.001);

            LBFGSMinimizer minimizer = new LBFGSMinimizer(0.001, 100);
            double[] min = minimizer.minimize(objective, a.toDouble(weights), true, null);
            System.arraycopy(min, 0, weights, 0, weights.length);

            System.exit(1);

            //    @Override
            //    public void callback(double[] guess, int iter, double val, double[] grad) {
            //        System.out.println("L = " + val);
            //    }
            //});
        }
    }

    protected void scaledAddI(float[] destination, float scale, Pair<int[], float[]> source) {
        for (int i = 0; i < source.getFirst().length; i++) {
            destination[source.getFirst()[i]] += scale * source.getSecond()[i];
        }
    }

    protected void scaledAddI(float[] destination, float scale, float[] source) {
        for (int i = 0; i < source.length; i++) {
            destination[i] = scale * source[i];
        }
    }

    protected float[] deltaEmission(int d, int t, int state, float[] weights, float[] negExpectedFeatValue) {
        float[] delta = new float[weights.length];
        System.arraycopy(negExpectedFeatValue, 0, delta, 0, negExpectedFeatValue.length);
        for (int i = 0; i < emissionFeatIndices[d][t][state].length; i++) {
            delta[emissionFeatIndices[d][t][state][i]] = emissionFeatValues[d][t][state][i];

        }
        return delta;
    }

    protected float[] deltaTransition(int state1, int state2, float[] weights) {
        float[] delta = new float[weights.length];
        for (int i = 0; i < transitionFeatIndices[state1][state2].length; i++) {
            delta[transitionFeatIndices[state1][state2][i]] = transitionFeatValues[state1][state2][i];
        }
        float normalizer = 0;
        for (int oState2 = 0; oState2 < CLUSTERS; oState2++) {
            float transScore = dot(transitionFeatIndices[state1][oState2], transitionFeatValues[state1][oState2], weights);
            normalizer += transScore;
            for (int i = 0; i < transitionFeatIndices[state1][state2].length; i++) {
                delta[transitionFeatIndices[state1][oState2][i]] -= transScore * transitionFeatValues[state1][oState2][i];
            }
        }
        for (int i = 0; i < delta.length; i++) {
            delta[i] /= normalizer;
        }
        return delta;
    }
    */

    public static void main(String[] args) {
        Clustering clustering = new Clustering(new File("/Users/jda/Corpora/imslp"));
        clustering.train();
    }
}
