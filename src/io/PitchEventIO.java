package io;

import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import tberg.murphy.arrays.a;
import tberg.murphy.fileio.f;
import tberg.murphy.tuple.Pair;

import javax.sound.midi.*;

public class PitchEventIO {

    public static final int MIDI_NOTE_ON = 0x90;
    public static final int MIDI_TEMPO = 0x51;
    public static final int MIDI_NOTE_OFF = 0x80;
    public static final int MIDI_CONTROLLER = 0xb0;
    public static final int MIDI_CONTROLLER_PEDAL = 0x40;

	public static final int MIN_MIDI_PITCH_ID = 21;
	public static final int MAX_MIDI_PITCH_ID = 109;
	
//	public static final int MIN_MIDI_PITCH_ID = 48;
//	public static final int MAX_MIDI_PITCH_ID = 73;
    
	public static final int N_MIDI_PITCH_IDS = MAX_MIDI_PITCH_ID - MIN_MIDI_PITCH_ID;
	public static final int N_CHROMA = 12;
    public static final int[] midiPitchIds = a.enumerate(MIN_MIDI_PITCH_ID, MAX_MIDI_PITCH_ID);

    // P&E
    
//    public static final float MIDI_ONSET_FACTOR = 1;
//    public static final float MIDI_OFFSET_FACTOR = 1;
//    public static final float PITCH_EVENT_BIAS_MS = -20.0f;

    // MAPS
    
    public static final float MIDI_ONSET_FACTOR = 0.996064f; 
    public static final float MIDI_OFFSET_FACTOR = 0.996064f;
    public static final float PITCH_EVENT_BIAS_MS = 0.0f;

	public static class Event implements java.io.Serializable {
		private static final long serialVersionUID = 1L;
		public final int noteIndex;
		public final float onsetSec;
		public final float offsetSec;
		public final int velocity;
		public Event(int noteIndex, float onsetSec, float offsetSec) {
			this.noteIndex = noteIndex;
			this.onsetSec = onsetSec;
			this.offsetSec = offsetSec;
			this.velocity = 80;
		}
		public Event(int noteIndex, float onsetSec, float offsetSec, int velocity) {
			this.noteIndex = noteIndex;
			this.onsetSec = onsetSec;
			this.offsetSec = offsetSec;
			this.velocity = velocity;
		}
        public String toString() {
            return "Event(" + noteIndex + ", " + onsetSec + ":" + offsetSec + ")";
        }
	}
	
	public static List<Event> readPitchEvents(String path) {
		List<Event> events = new ArrayList<Event>();
		List<String> lines = f.readLines(path);
		while (lines.get(0).trim().equals("")) lines.remove(0);
		lines.remove(0);
		for (String line : lines) {
			if (line.trim().equals("")) continue;
			String[] split = line.trim().split("\\s+");
			float onset = Math.max(0, Float.parseFloat(split[0]) + PITCH_EVENT_BIAS_MS*1e-3f);
			float offset = Math.max(0, Float.parseFloat(split[1]) + PITCH_EVENT_BIAS_MS*1e-3f);
			int noteId = Integer.parseInt(split[2]);
			int noteIndex = noteId - MIN_MIDI_PITCH_ID;
			events.add(new Event(noteIndex, onset, offset));
		}
		return events;
	}

    public static void writeMIDIFile(List<Event> events, String path) throws IOException, InvalidMidiDataException {
        Sequence seq = new Sequence(Sequence.PPQ, 1000, 1);
        long currentBPM = 120;
        for (Event event : events) {
            int pitch = event.noteIndex + MIN_MIDI_PITCH_ID;

            ShortMessage shortMessageOnset = new ShortMessage();
            shortMessageOnset.setMessage(ShortMessage.NOTE_ON, 0, pitch, event.velocity);
            MidiEvent midiEventOnset = new MidiEvent(shortMessageOnset, getMidiTick(seq, Math.max(0, event.onsetSec), currentBPM));

            ShortMessage shortMessageOffset = new ShortMessage();
            shortMessageOffset.setMessage(ShortMessage.NOTE_OFF, 0, pitch, event.velocity);
            MidiEvent midiEventOffset = new MidiEvent(shortMessageOffset, getMidiTick(seq, Math.max(0, event.offsetSec), currentBPM));

            // System.out.println(event.noteIndex);
            // System.out.println(event.onsetSec);
            // System.out.println(event.offsetSec);
            // System.out.println();

            seq.getTracks()[0].add(midiEventOnset);
            seq.getTracks()[0].add(midiEventOffset);

            //MidiEvent midiEvent = new MidiEvent()
            //seq.getTracks()[0].
        }
        MidiSystem.write(seq, MidiSystem.getMidiFileTypes(seq)[0], new File(path));
//        System.exit(1);
    }

    public static List<Event> readMIDIFile(String path) {

        Sequence seq = null;
		try {
            //System.out.println(path);
			seq = MidiSystem.getSequence(new File(path));
		} catch (Exception e) {
			//e.printStackTrace();
            throw new UnsupportedOperationException("Couldn't read song");
		}

        //int trackId;
        //if (path.contains("Poliner")) {
        //    trackId = 0;
        //} else {
        //    trackId = 0;
        //    if (seq.getTracks().length != 1) {
        //        throw new UnsupportedOperationException("Multi-track files not supported");
        //    }
        //}

        //System.out.println(seq.getTracks().length);
        //Track track = seq.getTracks()[trackId];

        boolean setTempo = false;

        List<Event> events = new ArrayList<Event>();

        long currentBPM = 120;

        for (Track track : seq.getTracks()) {

            boolean isPedalDown = false;
            List<Pair<Float,Integer>> waitingForKeyUp = new ArrayList<Pair<Float,Integer>>();
            List<Pair<Float,Integer>> waitingForPedalUp = new ArrayList<Pair<Float,Integer>>();

            //System.out.println(path);
            //System.out.println(track.size());
            //System.exit(1);

            for (int i = 0; i < track.size(); ++i) {
                //System.out.println("evt");
                MidiEvent evt = track.get(i);
                MidiMessage msg = evt.getMessage();
                if (msg instanceof MetaMessage) {
                    MetaMessage mm = (MetaMessage)msg;
                    if (mm.getType() != MIDI_TEMPO) {
                        continue;
                    }
                    //System.out.println("METAMESSAGE");
                    //System.out.println("TEMPO MESSAGE");
                    int microsecondsPerQuarter = new BigInteger(mm.getData()).intValue();
                    currentBPM = 60000000 / microsecondsPerQuarter;
                    System.out.println("current BPM is " + currentBPM);
                    //System.out.println(60000000 / microsecondsPerQuarter);
                    //System.out.println(Arrays.toString(mm.getData()));
                    //System.out.println(mm);
                    if (setTempo) {
                        System.err.println("Warning: changing tempo midstream (probably be in the wrong place");
                    }
                    setTempo = true;
                }
                if (!(msg instanceof ShortMessage)) {
                    continue;
                }
                ShortMessage sm = (ShortMessage)msg;
                int midiPitch = sm.getData1();
                int pitch = midiPitch - MIN_MIDI_PITCH_ID;
                float time = getRealTime(seq, evt.getTick(), currentBPM);
                // apparently, NOTE_ON events with a velocity of 0 are treated as NOTE_OFF events
                if (sm.getCommand() == MIDI_NOTE_ON && sm.getData2() > 0) {
                    //onsets[frame][pitch] = true;
                    waitingForKeyUp.add(new Pair<Float, Integer>(time, pitch));
                    //System.out.println("key down " + pitch + " @" + time);
                } else if (sm.getCommand() == MIDI_NOTE_OFF || sm.getCommand() == MIDI_NOTE_ON && sm.getData2() == 0) {
                    //System.out.println("OFFSET");
                    //offsets[frame][pitch] = true;
                    List<Pair<Float,Integer>> toRemove = new ArrayList<Pair<Float,Integer>>();
                    for (Pair<Float,Integer> onset : waitingForKeyUp) {
                        if (onset.getSecond() == pitch) {
                            toRemove.add(onset);
                            //System.out.println("key up " + pitch + " @" + time);
                        }
                    }
                    if (isPedalDown) {
                        waitingForPedalUp.addAll(toRemove);
                    } else {
                        for (Pair<Float,Integer> onset : toRemove) {
                            events.add(new Event(onset.getSecond(), onset.getFirst() * MIDI_ONSET_FACTOR, time * MIDI_OFFSET_FACTOR));
                        }
                    }
                    waitingForKeyUp.removeAll(toRemove);
                } else if (sm.getCommand() == MIDI_CONTROLLER && (sm.getData1() == 0x41 || sm.getData1() == 0x42 || sm.getData1() == 0x43 || sm.getData1() == 0x44)) {
                    //System.out.println("weird pedal " + sm.getData1());
                    //System.out.println("EVT @" + time + " " + sm.toString());
                } else if (sm.getCommand() == MIDI_CONTROLLER && sm.getData1() == MIDI_CONTROLLER_PEDAL) {
                    //System.out.println(sm.getData2());
                    if (sm.getData2() >= 64) {
                        // pedal down
                        isPedalDown = true;
                        //System.out.println("pedal down @" + time);
                    } else {
                        // pedal up
                        //System.out.println("pedal up @" + time);
                        isPedalDown = false;
                        for (Pair<Float,Integer> onset : waitingForPedalUp) {
                            events.add(new Event(onset.getSecond(), onset.getFirst() * MIDI_ONSET_FACTOR, time * MIDI_OFFSET_FACTOR));
                        }
                        waitingForPedalUp.clear();
                    }
                } else {
                    //System.err.println("Unknown command");
                    //System.err.println(sm.getCommand());
                    //System.err.println(sm.getData1());
                    //throw new UnsupportedOperationException("Unknown MIDI command");
                }
            }

            assert waitingForKeyUp.isEmpty();
            assert waitingForPedalUp.isEmpty();

        }

        Collections.sort(events, new Comparator<Event>() {
            public int compare(Event a, Event b) {
                return (int)Math.signum(a.onsetSec - b.onsetSec);
            }
        });

        //System.out.println(events);
        //System.out.println(waitingForKeyUp.size());
        //System.out.println(track.size());
        //System.out.println(events.size());
        
//        return events;
 
        List<Event> eventsBiased = new ArrayList<Event>();
        for (Event event : events) {
        	eventsBiased.add(new Event(event.noteIndex, Math.max(0, event.onsetSec + PITCH_EVENT_BIAS_MS*1e-3f), Math.max(0, event.offsetSec + PITCH_EVENT_BIAS_MS*1e-3f)));
        }
        return eventsBiased;
    }

    public static float getRealTime(Sequence seq, long tick, long bpm) {
        float dt = seq.getDivisionType();
        float ticksPerSecond;
        if (dt == Sequence.PPQ) {
            ticksPerSecond = seq.getResolution() * (bpm / 60.0f);
        } else {
            System.err.println("Warning: weird timing format");
            double framesPerSecond =
                    (dt == Sequence.SMPTE_24 ? 24
                            : (dt == Sequence.SMPTE_25 ? 25
                            : (dt == Sequence.SMPTE_30 ? 30
                            : (dt == Sequence.SMPTE_30DROP ? 29.97 : -1)))); // something is wrong if we got to the end
            ticksPerSecond = seq.getResolution() * ((float) framesPerSecond);
        }
        return 1.0f * tick / ticksPerSecond;
    }

    public static long getMidiTick(Sequence seq, float realTime, long bpm) {
        float dt = seq.getDivisionType();
        float ticksPerSecond;
        if (dt == Sequence.PPQ) {
            ticksPerSecond = seq.getResolution() * (bpm / 60.0f);
        } else {
            System.err.println("Warning: weird timing format");
            double framesPerSecond =
                    (dt == Sequence.SMPTE_24 ? 24
                            : (dt == Sequence.SMPTE_25 ? 25
                            : (dt == Sequence.SMPTE_30 ? 30
                            : (dt == Sequence.SMPTE_30DROP ? 29.97 : -1)))); // something is wrong if we got to the end
            ticksPerSecond = seq.getResolution() * ((float) framesPerSecond);
        }
        return (long) (realTime * ticksPerSecond);
    }
    
    public static List<Event> takePrefix(List<Event> events, float prefixLengthSec) {
    	List<Event> result = new ArrayList<Event>();
    	for (Event event : events) {
    		if (event.onsetSec < prefixLengthSec) {
    			result.add(new Event(event.noteIndex, event.onsetSec, Math.min(event.offsetSec, prefixLengthSec)));
    		}
    	}
    	return result;
    }

}
