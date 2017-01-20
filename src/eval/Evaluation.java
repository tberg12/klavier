package eval;

import io.MatrixVis;
import io.PitchEventIO;
import io.PitchEventIO.Event;

import java.util.*;

import eval.PitchEventUtil.NoteState;
import tberg.murphy.arrays.a;
import main.Main;

public class Evaluation {
	
    public static final float MIREX_EVAL_FRAME_LENGTH_MS = 10.0f;
    //public static final float MIREX_EVAL_FRAME_LENGTH_MS = 23.2f;

	public static class EvalSuffStats {
	    private float truePositive;
	    private float falsePositive;
	    private float falseNegative;

	    public EvalSuffStats(float truePositive, float falsePositive, float falseNegative) {
	        this.truePositive = truePositive;
	        this.falsePositive = falsePositive;
	        this.falseNegative = falseNegative;
	    }

	    public void increment(EvalSuffStats other) {
	    	this.truePositive += other.truePositive;
	    	this.falsePositive += other.falsePositive;
	    	this.falseNegative += other.falseNegative;
	    }
	    
	    public float getPrecision() {
	        return truePositive / (truePositive + falsePositive);
	    }

	    public float getRecall() {
	        return truePositive / (truePositive + falseNegative);
	    }

	    public float getF1() {
            if (getPrecision() + getRecall() == 0) {
                return 0;
            }
	        return 2.0f * getPrecision() * getRecall() / (getPrecision() + getRecall());
	    }
	    
	    public float getAcc() {
	    	return truePositive / (truePositive + falseNegative + falsePositive);
	    }

	    public String toString() {
	    	return String.format("prec: %f, recall: %f, f1: %f, acc: %f", getPrecision(), getRecall(), getF1(), getAcc());
	    }
	}

    public static EvalSuffStats evaluateFrames(List<Event> guessEvents, List<Event> goldEvents, float prefixLengthSec) {
    	guessEvents = PitchEventIO.takePrefix(guessEvents, prefixLengthSec);
		goldEvents = PitchEventIO.takePrefix(goldEvents, prefixLengthSec);
		float maxOffsetSec = Float.NEGATIVE_INFINITY;
		for (Event event : guessEvents) {
			maxOffsetSec = Math.max(event.offsetSec, maxOffsetSec);
		}
		for (Event event : goldEvents) {
			maxOffsetSec = Math.max(event.offsetSec, maxOffsetSec);
		}
    	
        int nFrames = (int)(maxOffsetSec / MIREX_EVAL_FRAME_LENGTH_MS * 1000);

        int[][] goldRoll = new int[nFrames][PitchEventIO.N_MIDI_PITCH_IDS];
        int[][] predRoll = new int[nFrames][PitchEventIO.N_MIDI_PITCH_IDS];

        for (Event event : goldEvents) {
            int startFrame = (int)Math.floor(event.onsetSec / MIREX_EVAL_FRAME_LENGTH_MS * 1000); // round DOWN for start of frame
            int endFrame = (int)Math.ceil(event.offsetSec / MIREX_EVAL_FRAME_LENGTH_MS * 1000); // round UP for end of frame
            for (int f = Math.max(startFrame, 0); f < Math.min(endFrame, nFrames); f++) {
                goldRoll[f][event.noteIndex]++;
            }
        }

        for (Event event : guessEvents) {
            int startFrame = (int)Math.floor(event.onsetSec / MIREX_EVAL_FRAME_LENGTH_MS * 1000); // round DOWN for start of frame
            int endFrame = (int)Math.ceil(event.offsetSec / MIREX_EVAL_FRAME_LENGTH_MS * 1000); // round UP for end of frame
            for (int f = Math.max(startFrame, 0); f < Math.min(endFrame, nFrames); f++) {
                predRoll[f][event.noteIndex]++;
            }
        }

        int tp = 0, fp = 0, fn = 0;

        for (int frame = 0; frame < nFrames; frame++) {
            for (int pitch = 0; pitch < PitchEventIO.N_MIDI_PITCH_IDS; pitch++) {

                int nGold = goldRoll[frame][pitch];
                int nPred = predRoll[frame][pitch];

                // o'hanlon eval (no double-counting frames with both an onset and an offset)
                if (nGold > 0 && nPred > 0) {
                    tp++;
                } else if (nGold > 0) {
                    fn++;
                } else if (nPred > 0) {
                    fp++;
                }

                // (corrected) vincent eval (with double-counting)
//                if (nPred >= nGold) {
//                    tp += nGold;
//                    fp += nPred - nGold;
//                } else {
//                    tp += nPred;
//                    fn += nGold - nPred;
//                }
            }
        }

        return new EvalSuffStats(tp, fp, fn);

    }
	
//    public static EvalSuffStats evaluateFrames(List<Event> guessEvents, List<Event> goldEvents, float prefixLengthSec) {
//    	guessEvents = PitchEventIO.takePrefix(guessEvents, prefixLengthSec);
//		goldEvents = PitchEventIO.takePrefix(goldEvents, prefixLengthSec);
//		float maxOffsetSec = Float.NEGATIVE_INFINITY;
//		for (Event event : guessEvents) {
//			maxOffsetSec = Math.max(event.offsetSec, maxOffsetSec);
//		}
//		for (Event event : goldEvents) {
//			maxOffsetSec = Math.max(event.offsetSec, maxOffsetSec);
//		}
//
//		NoteState[][] guess = PitchEventUtil.convertToRoll(guessEvents, PitchEventUtil.nearestFrameByCenter(maxOffsetSec, MIREX_EVAL_FRAME_LENGTH_MS * 1e-3f)+1, MIREX_EVAL_FRAME_LENGTH_MS * 1e-3f);
//
//        //System.out.println(a.toString(guess));
//        for (int t = 50; t < 85; t++) {
//            System.out.print(t + "\t");
//            for (int i = 0; i < guess[t].length; i++) {
//                if (guess[t][i] == NoteState.OFF) {
//                    System.out.print("  ");
//                } else if (guess[t][i] == NoteState.ONSET) {
//                    System.out.print("! ");
//                } else {
//                    System.out.print(". ");
//                }
//            }
//            System.out.println();
//            //System.out.println(a.toString(guess[t]));
//        }
//        System.exit(1);
//
//		NoteState[][] gold = PitchEventUtil.convertToRoll(goldEvents, PitchEventUtil.nearestFrameByCenter(maxOffsetSec, MIREX_EVAL_FRAME_LENGTH_MS * 1e-3f)+1, MIREX_EVAL_FRAME_LENGTH_MS * 1e-3f);
//
//        //System.out.println(guess.length);
//        //System.out.println(gold.length);
//        //System.exit(1);
//
//		int truePositive = 0;
//		int falsePositive = 0;
//		int falseNegative = 0;
//
//		for (int pitchId=0; pitchId<PitchEventIO.N_MIDI_PITCH_IDS; ++pitchId) {
//			for (int t=0; t<guess.length; ++t) {
//				if (guess[t][pitchId] != NoteState.OFF && gold[t][pitchId] != NoteState.OFF) {
//					truePositive++;
//				} else if (guess[t][pitchId] == NoteState.OFF && gold[t][pitchId] != NoteState.OFF) {
//					falseNegative++;
//				} else if (guess[t][pitchId] != NoteState.OFF && gold[t][pitchId] == NoteState.OFF) {
//					falsePositive++;
//				}
//			}
//		}
//
//		return new EvalSuffStats(truePositive, falsePositive, falseNegative);
//    }

    static interface EvalPredicate {
        public boolean apply(Event guessEvent, Event goldEvent);
    }

    public static EvalSuffStats evaluateOnsets(List<Event> guessEvents, List<Event> goldEvents, float prefixLengthSec, final float toleranceWindowMs) {
        return evaluateNoteTracking(guessEvents, goldEvents, prefixLengthSec, new EvalPredicate() {
            public boolean apply(Event guessEvent, Event goldEvent) {
                return Math.abs(guessEvent.onsetSec - goldEvent.onsetSec) <= (toleranceWindowMs / 2.0f) * 1e-3f;
                //if (Math.abs(guessOnset - goldOnset) <= (toleranceWindowMs / 2.0f) * 1e-3f) {
            }
        });
    }

    public static EvalSuffStats evaluateFullNote(List<Event> guessEvents, List<Event> goldEvents, float prefixLengthSec, final float onsetToleranceMs, final float offsetToleranceFrac) {
        //guessEvents = goldEvents;
        return evaluateNoteTracking(guessEvents, goldEvents, prefixLengthSec, new EvalPredicate() {
            public boolean apply(Event guessEvent, Event goldEvent) {
                float goldWidth = goldEvent.offsetSec - goldEvent.onsetSec;
                return Math.abs(guessEvent.onsetSec - goldEvent.onsetSec) <= (onsetToleranceMs / 2.0f) * 1e-3f &&
                       //Math.abs(guessEvent.offsetSec - goldEvent.offsetSec) <= (goldEvent.onsetSec - goldEvent.offsetSec) * offsetToleranceFrac / 2.0f * 1e-3f
//                       guessEvent.offsetSec >= goldWidth * (1 + offsetToleranceFrac) &&
//                       guessEvent.offsetSec <= goldWidth * (1 - offsetToleranceFrac);
                       Math.abs(guessEvent.offsetSec - goldEvent.offsetSec) <= goldWidth * offsetToleranceFrac;
                //if (Math.abs(guessOnset - goldOnset) <= (toleranceWindowMs / 2.0f) * 1e-3f) {
            }
        });
    }

	public static EvalSuffStats evaluateNoteTracking(List<Event> guessEvents, List<Event> goldEvents, float prefixLengthSec, EvalPredicate evalPredicate) {
		guessEvents = PitchEventIO.takePrefix(guessEvents, prefixLengthSec);	
		goldEvents = PitchEventIO.takePrefix(goldEvents, prefixLengthSec);	

		int truePositive = 0;
		int falsePositive = 0;
		int falseNegative = 0;
		
		Comparator<Event> onsetComp = new Comparator<Event>() {
			public int compare(Event o1, Event o2) {
				if (o1.onsetSec < o2.onsetSec) {
					return -1;
				} else if (o1.onsetSec > o2.onsetSec) {
					return 1;
				} else {
					return 0;
				}
			}
		};
		List<Event> sortedGuessEvents = new ArrayList<Event>(guessEvents);
		Collections.sort(sortedGuessEvents, onsetComp);
		List<Event> sortedGoldEvents = new ArrayList<Event>(goldEvents);
		Collections.sort(sortedGoldEvents, onsetComp);
		
		
		for (int pitchId=0; pitchId<PitchEventIO.N_MIDI_PITCH_IDS; ++pitchId) {
			//List<Float> guessOnsetsList = new ArrayList<Float>();
            List<Event> guessEventsList = new ArrayList<Event>();
			for (Event event : guessEvents) {
				if (event.noteIndex == pitchId) {
					//guessOnsetsList.add(event.onsetSec);
                    guessEventsList.add(event);
				}
			}
			//float[] guessOnsets = a.toFloatArray(guessOnsetsList);
			//boolean[] guessMapped = new boolean[guessOnsets.length];
            boolean[] guessMapped = new boolean[guessEvents.size()];

			//List<Float> goldOnsetsList = new ArrayList<Float>();
            List<Event> goldEventsList = new ArrayList<Event>();
			for (Event event : goldEvents) {
				if (event.noteIndex == pitchId) {
					//goldOnsetsList.add(event.onsetSec);
                    goldEventsList.add(event);
				}
			}
			//float[] goldOnsets = a.toFloatArray(goldOnsetsList);
			//boolean[] goldMapped = new boolean[goldOnsets.length];
            boolean[] goldMapped = new boolean[goldEvents.size()];
			
			//for (int guessi=0; guessi<guessOnsets.length; ++guessi) {
            for (int guessi=0; guessi<guessEventsList.size(); ++guessi) {
				//float guessOnset = guessOnsets[guessi];
				//for (int goldi=0; goldi<goldOnsets.length; ++goldi) {
                for (int goldi=0; goldi<goldEventsList.size(); ++goldi) {
					if (!goldMapped[goldi]) {
						//float goldOnset = goldOnsets[goldi];
                        if (evalPredicate.apply(guessEventsList.get(guessi), goldEventsList.get(goldi))) {
						//if (Math.abs(guessOnset - goldOnset) <= (toleranceWindowMs / 2.0f) * 1e-3f) {
							guessMapped[guessi] = true;
							goldMapped[goldi] = true;
							truePositive++;
						}
					}
				}
				if (!guessMapped[guessi]) falsePositive++;
			}
			//for (int goldi=0; goldi<goldOnsets.length; ++goldi) {
            for (int goldi = 0; goldi < goldEventsList.size(); ++goldi) {
				if (!goldMapped[goldi]) falseNegative++;
			}
		}
		
		return new EvalSuffStats(truePositive, falsePositive, falseNegative);
	}

//  public static EvalSuffStats evaluateFrames(NoteState[][] predicted, List<Event> events, float secondsPerFrame, float prefixLengthSec) {
//	//int numEvalFrames = (int) Math.min(Math.floor(prefixLengthSec / secondsPerFrame), predicted.length);
//    int numEvalMirexFrames = (int) Math.min(Math.floor(prefixLengthSec / FRAME_EVAL_MIREX_FRAME_SIZE_SEC), predicted.length * secondsPerFrame / FRAME_EVAL_MIREX_FRAME_SIZE_SEC);
//
//    //NoteState[][] gold = convertToRoll(events, numEvalMirexFrames, FRAME_EVAL_MIREX_FRAME_SIZE_SEC);
//
//    boolean[][] gold = new boolean[numEvalMirexFrames][predicted[0].length];
//    for (Event e : events) {
//        int t = 0;
//        while (e.onsetSec / FRAME_EVAL_MIREX_FRAME_SIZE_SEC + t <= e.offsetSec / FRAME_EVAL_MIREX_FRAME_SIZE_SEC) {
//            gold[(int)(e.onsetSec / FRAME_EVAL_MIREX_FRAME_SIZE_SEC + t)][e.noteIndex] = true;
//            t++;
////            System.out.println((int)(e.onsetSec / FRAME_EVAL_MIREX_FRAME_SIZE_SEC + t));
//        }
//    }
//
//    //System.out.println(secondsPerFrame);
//    //System.exit(1);
//    // mirex frames are shorter
//    assert FRAME_EVAL_MIREX_FRAME_SIZE_SEC < secondsPerFrame;
//
//    int tp = 0, fp = 0, fn = 0;
//    for (int frame = 0; frame < numEvalMirexFrames; frame++) {
//        int realFrame = (int) (frame * FRAME_EVAL_MIREX_FRAME_SIZE_SEC / secondsPerFrame);
//        //System.out.println(realFrame);
//        for (int pitch = 0; pitch < gold[0].length; pitch++) {
//            boolean predOn = predicted[realFrame][pitch] != NoteState.OFF;
//            boolean goldOn = gold[frame][pitch]; // != NoteState.OFF;
//
//            if (predOn && goldOn) {
//                tp++;
//            } else if (predOn && !goldOn) {
//                fp++;
//            } else if (!predOn && goldOn) {
//                fn++;
//            }
//        }
//    }
//    //for (int t=0; t<numEvalFrames; t++) {
//    //    for (int pitch = 0; pitch < gold[t].length; pitch++) {
//    //        boolean predOn = predicted[t][pitch] != NoteState.OFF;
//    //        boolean goldOn = gold[t][pitch] != NoteState.OFF;
//    //        if (predOn && goldOn) {
//    //            tp++;
//    //        } else if (predOn && !goldOn) {
//    //            fp++;
//    //        } else if (!predOn && goldOn) {
//    //            fn++;
//    //        }
//    //    }
//    //}
//    return new EvalSuffStats(tp, fp, fn);
//}
//
//public static EvalSuffStats evaluateOnsets(NoteState[][] predicted, List<Event> events, float secondsPerFrame, float prefixLengthSec, float toleranceWindowMs) {
//	int numEvalFrames = (int) Math.min(PitchEventUtil.frameNum(prefixLengthSec, secondsPerFrame), predicted.length);
//	float lengthSec = numEvalFrames * secondsPerFrame;
//	
//    int truePositive = 0;
//    int falsePositive = 0;
//    int falseNegative = 0;
//
//    List<Event> sortedEvents = new ArrayList<Event>(events);
//    Collections.sort(sortedEvents, new Comparator<Event>() {
//		public int compare(Event o1, Event o2) {
//			if (o1.onsetSec < o2.onsetSec) {
//				return -1;
//			} else if (o1.onsetSec > o2.onsetSec) {
//				return 1;
//			} else {
//				return 0;
//			}
//		}
//	});
//    
//    float halfToleranceWindowSec = (toleranceWindowMs / 2.0f) * 1e-3f;
//    boolean[][] mappedPredicted = new boolean[numEvalFrames][PitchEventIO.N_MIDI_PITCH_IDS];
//    for (Event event : sortedEvents) {
//    	if (event.onsetSec < lengthSec) {
//    		boolean mapped = false;
//    		for (int t : PitchEventUtil.frameNumsCenterWithinInterval(event.onsetSec - halfToleranceWindowSec, event.onsetSec + halfToleranceWindowSec, secondsPerFrame, 0, numEvalFrames-1)) {
//    			if (predicted[t][event.noteIndex] == NoteState.ONSET && !mappedPredicted[t][event.noteIndex]) {
//    				truePositive += 1.0;
//    				mappedPredicted[t][event.noteIndex] = true;
//    				mapped = true;
//    				break;
//    			}
//    		}
//    		
//    		if (!mapped) {
//    			falseNegative += 1.0;
//    		}
//    	}
//    }
//    
//    for (int t=0; t<numEvalFrames; ++t) {
//    	for (int noteIndex=0; noteIndex<PitchEventIO.N_MIDI_PITCH_IDS; ++noteIndex) {
//    		if (predicted[t][noteIndex] == NoteState.ONSET && !mappedPredicted[t][noteIndex]) {
//    			falsePositive += 1.0;
//    		}
//    	}
//    }
//
//    return new EvalSuffStats(truePositive, falsePositive, falseNegative);
//}
	
}
