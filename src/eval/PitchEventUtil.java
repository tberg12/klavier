package eval;

import io.PitchEventIO;
import io.PitchEventIO.Event;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import learning.EnvelopeActivationModel;
import tberg.murphy.arrays.a;

public class PitchEventUtil {
	
	public static enum NoteState {OFF, SUSTAIN, ONSET};

	public static float frameOnsetSec(int frameNum, float secondsPerFrame) {
		return frameCenterSec(frameNum, secondsPerFrame) - (secondsPerFrame / 2.0f);
	}
	
	public static float frameOffsetSec(int frameNum, float secondsPerFrame) {
		return frameCenterSec(frameNum, secondsPerFrame) + (secondsPerFrame / 2.0f);
	}
	
	public static float frameCenterSec(int frameNum, float secondsPerFrame) {
		return frameNum * secondsPerFrame;
	}

	public static int nextFrameByCenter(float timeSec, float secondsPerFrame) {
		return (int) Math.ceil(timeSec / secondsPerFrame);
	}
	
	public static int prevFrameByCenter(float timeSec, float secondsPerFrame) {
		return (int) Math.floor(timeSec / secondsPerFrame);
	}
	
	public static int nearestFrameByCenter(float timeSec, float secondsPerFrame) {
		return (int) Math.round(timeSec / secondsPerFrame);
	}
	
	public static int nearestFrameByOnset(float timeSec, float secondsPerFrame) {
		return nearestFrameByCenter(timeSec + (secondsPerFrame / 2.0f), secondsPerFrame);
	}
	
	public static int nearestFrameByOffset(float timeSec, float secondsPerFrame) {
		return nearestFrameByCenter(timeSec - (secondsPerFrame / 2.0f), secondsPerFrame);
	}
	
	public static int[] framesByCenterWithinInterval(float startSec, float stopSec, float secondsPerFrame, int minFrameNum, int maxFrameNum) {
		int onsetFrame = Math.max(minFrameNum, Math.min(maxFrameNum, nextFrameByCenter(startSec, secondsPerFrame)));
		int offsetFrame = Math.max(minFrameNum, Math.min(maxFrameNum, prevFrameByCenter(stopSec, secondsPerFrame)));
		if (onsetFrame < offsetFrame+1) {
			return a.enumerate(onsetFrame, offsetFrame+1);
		} else {
			return new int[0];
		}
	}
	
	public static int[] framesByOverlap(float startSec, float stopSec, float secondsPerFrame, int minFrameNum, int maxFrameNum) {
		int onsetFrame = Math.max(minFrameNum, Math.min(maxFrameNum, nearestFrameByOnset(startSec, secondsPerFrame)));
		int offsetFrame = Math.max(minFrameNum, Math.min(maxFrameNum, nearestFrameByOffset(stopSec, secondsPerFrame)));
		if (onsetFrame < offsetFrame+1) {
			return a.enumerate(onsetFrame, offsetFrame+1);
		} else {
			return new int[0];
		}
	}
	
	////////////////////////////////////////////////////////////////////////////////////
	
	public static List<Event> convertToEvents(NoteState[][] noteStates, float secondsPerFrame) {
		List<Event> events = new ArrayList<Event>();
		for (int pitchId=0; pitchId<PitchEventIO.N_MIDI_PITCH_IDS; ++pitchId) {
			boolean on = false;
			float onsetSec = Float.NaN;
			for (int t=0; t<noteStates.length; ++t) {
				if (!on && noteStates[t][pitchId] == NoteState.ONSET) {
					onsetSec = frameOnsetSec(t, secondsPerFrame);
					on = true;
				} else if (on && noteStates[t][pitchId] == NoteState.OFF) {
					float offsetSec = frameOffsetSec(t-1, secondsPerFrame);
					events.add(new Event(pitchId, onsetSec, offsetSec));
					on = false;
				} else if (on && t == noteStates.length-1) {
					float offsetSec = frameOffsetSec(t, secondsPerFrame);
					events.add(new Event(pitchId, onsetSec, offsetSec));
					on = false;
				}
			}
		}
		return events;
	}
	
	public static List<Event> convertToEvents(NoteState[][] noteStates, float[][] activations, float[][] envlopes, float secondsPerFrame) {
		float[][] transActivations = a.transpose(activations);
		List<Event> events = new ArrayList<Event>();
		for (int pitchId=0; pitchId<PitchEventIO.N_MIDI_PITCH_IDS; ++pitchId) {
			boolean on = false;
			int onsetFrame = -1;
			float onsetSec = Float.NaN;
			for (int t=0; t<noteStates.length; ++t) {
				if (!on && noteStates[t][pitchId] == NoteState.ONSET) {
					onsetFrame = t;
					onsetSec = frameOnsetSec(t, secondsPerFrame);
					on = true;
				} else if (on && noteStates[t][pitchId] == NoteState.OFF) {
					float offsetSec = frameOffsetSec(t-1, secondsPerFrame);
					float volume = EnvelopeActivationModel.bestVolume(onsetFrame, t-onsetFrame, transActivations[pitchId], envlopes[pitchId]); 
					events.add(new Event(pitchId, onsetSec, offsetSec, volumeToVelocity(volume, pitchId)));
					on = false;
				} else if (on && t == noteStates.length-1) {
					float offsetSec = frameOffsetSec(t, secondsPerFrame);
					float volume = EnvelopeActivationModel.bestVolume(onsetFrame, t-onsetFrame, transActivations[pitchId], envlopes[pitchId]); 
					events.add(new Event(pitchId, onsetSec, offsetSec, volumeToVelocity(volume, pitchId)));
					on = false;
				}
			}
		}
		return events;
	}
	
	private static int volumeToVelocity(double volume, int pitchId) {
		final double c = 1.0;
		double rescaled = (Math.log(volume+c) - Math.log(EnvelopeActivationModel.MIN_VOLUME+c)) / (Math.log(EnvelopeActivationModel.MAX_VOLUME+c) - Math.log(EnvelopeActivationModel.MIN_VOLUME+c));
		final double d = 0.3;
		double rerescaled = d + (1.0 - d) * rescaled;
		int velocity = (int) Math.floor(rerescaled * 127);
		return velocity;
	}
	
	public static List<Event> filter(List<Event> events, int minNoteIndex, int maxNoteIndex, int minVelocity) {
		List<Event> result = new ArrayList<Event>();
		for (Event event : events) {
			if (event.noteIndex >= minNoteIndex && event.noteIndex <= maxNoteIndex && event.velocity >= minVelocity) {
				result.add(event);
			}
		}
		return result;
	}
	
	////////////////////////////////////////////////////////////////////////////////////
	
	public static NoteState[][] convertToRoll(List<Event> events, int numFrames, float secondsPerFrame, float beforeOnsetToleranceWindowMs, float afterOnsetToleranceWindowMs) {
		return convertToRoll(events, numFrames, secondsPerFrame, true, beforeOnsetToleranceWindowMs, afterOnsetToleranceWindowMs);
	}
	
	public static NoteState[][] convertToRoll(List<Event> events, float secondsPerFrame, float beforeOnsetToleranceWindowMs, float afterOnsetToleranceWindowMs) {
		int maxFrame = -1;
		for (Event event : events) {
			int offsetFrame = nearestFrameByCenter(event.offsetSec, secondsPerFrame);
			maxFrame = Math.max(offsetFrame, maxFrame);
			int onsetFrame = nearestFrameByCenter(event.onsetSec + (afterOnsetToleranceWindowMs*1e-3f), secondsPerFrame);
			maxFrame = Math.max(onsetFrame, maxFrame);
		}
		return convertToRoll(events, maxFrame+1, secondsPerFrame, false, beforeOnsetToleranceWindowMs, afterOnsetToleranceWindowMs);
	}
	
	public static NoteState[][] convertToRoll(List<Event> events, int numFrames, float secondsPerFrame) {
		return convertToRoll(events, numFrames, secondsPerFrame, false, 0, 0);
	}
	
	public static NoteState[][] convertToRoll(List<Event> events, float secondsPerFrame) {
		int maxFrame = -1;
		for (Event event : events) {
			int frame = nearestFrameByCenter(event.offsetSec, secondsPerFrame);
			maxFrame = Math.max(frame, maxFrame);
		}
		return convertToRoll(events, maxFrame+1, secondsPerFrame, false, 0, 0);
	}

    public static float[][] convertToFloatRoll(List<Event> events, float secondsPerFrame) {
        NoteState[][] noteRoll = convertToRoll(events, secondsPerFrame);
        float[][] floatRoll = new float[noteRoll.length][noteRoll[0].length * 3];
        for (int i = 0; i < noteRoll.length; i++) {
            for (int j = 0; j < noteRoll[i].length; j++) {
                if (noteRoll[i][j] == NoteState.OFF) {
                    floatRoll[i][j*3] = 1f;
                    floatRoll[i][j*3+1] = 1f;
                } else if (noteRoll[i][j] == NoteState.SUSTAIN) {
                    floatRoll[i][j*3] = 0.667f;
                    floatRoll[i][j*3+1] = 0.667f;
                } else {
                    floatRoll[i][j*3] = 0f;
                    floatRoll[i][j*3+1] = 0f;
                }
                floatRoll[i][j*3+2] = 1f;
            }
        }
        for (int j = 0; j < noteRoll.length; j++) {
            floatRoll[j][0] = 0;
        }
        return floatRoll;
    }

	private static NoteState[][] convertToRoll(List<Event> events, int numFrames, float secondsPerFrame, boolean onsetWindow, float beforeOnsetToleranceWindowMs, float afterOnsetToleranceWindowMs) {
		NoteState[][] result = new NoteState[numFrames][PitchEventIO.N_MIDI_PITCH_IDS];
		for (int i=0; i<result.length; ++i) Arrays.fill(result[i], NoteState.OFF);
		for (Event event : events) {
			for (int t : framesByOverlap(event.onsetSec, event.offsetSec, secondsPerFrame, 0, numFrames-1)) {
				result[t][event.noteIndex] = NoteState.SUSTAIN;
			}
		}
		for (Event event : events) {
			if (onsetWindow) {
				for (int t : framesByOverlap(event.onsetSec - (beforeOnsetToleranceWindowMs*1e-3f), Math.min(event.offsetSec, event.onsetSec + (afterOnsetToleranceWindowMs*1e-3f)), secondsPerFrame, 0, numFrames-1)) {
					result[t][event.noteIndex] = NoteState.ONSET;
				}
			}
			result[Math.max(0, Math.min(numFrames-1, nearestFrameByOnset(event.onsetSec, secondsPerFrame)))][event.noteIndex] = NoteState.ONSET;
		}
		for (int j=0; j<result[0].length; ++j) {
			for (int i=1; i<result.length; ++i) {
				if (result[i-1][j] == NoteState.OFF && result[i][j] == NoteState.SUSTAIN) {
					assert false;
				}
			}
		}
		return result;
	}

}
