package lm;

import java.io.File;
import java.io.FilenameFilter;
import java.util.*;

import eval.PitchEventUtil;
import eval.PitchEventUtil.NoteState;
import io.PitchEventIO;

/**
 * @author jda
 */
public class NGramLanguageModel implements LanguageModel {

    private HashMap<Object, HashMap<Object, Double>> scores;

    public NGramLanguageModel(File corpusRoot, int n) {
        scores = new HashMap<Object, HashMap<Object, Double>>();
        assert corpusRoot.isDirectory();
        int count = 0;
        //for (File f : corpusRoot.listFiles) {
        {
            File f = new File("/Users/jda/Corpora/imslp/IMSLP00001#1#Smartscore-10.2.1.mid");
            try {
                List<PitchEventIO.Event> events = PitchEventIO.readMIDIFile(f.getAbsolutePath());
                NoteState[][] noteStates = PitchEventUtil.convertToRoll(events, 0.01f);
                List<Object> context = new LinkedList<Object>();
                for (NoteState[] frame : noteStates) {
                    Object repr = represent(frame);
                    if (!scores.containsKey(context)) {
                        scores.put(context, new HashMap<Object, Double>());
                    }
                    if (!scores.get(context).containsKey(repr)) {
                        scores.get(context).put(repr, 0d);
                    }
                    scores.get(context).put(repr, scores.get(context).get(repr) + 1);
                    context.add(repr);
                    context.remove(0);
                }
            } catch(Exception e) {
                e.printStackTrace();
                assert false;
            }
            //if(count++ < 5) {
            //    break;
            //}
        }
        System.out.println(scores);
    }

    public Object represent(NoteState[] frame) {
        HashSet<Integer> o = new HashSet<Integer>();
        for (int i = 0; i < frame.length; i++) {
            if (frame[i] != NoteState.OFF) {
                o.add(i);
            }
        }
        return o;
    }

    @Override
    public Object sample(List context) {
        return null;
    }

    @Override
    public double score(Object entry, List context) {
        return 0;
    }

    public static void main(String[] args) {
        new NGramLanguageModel(new File("/Users/jda/Corpora/imslp"), 2);
    }
}
