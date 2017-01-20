package io;

import io.AudioIO.Wave;
import io.PitchEventIO.Event;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;

import learning.Preprocessing.SpectralAtomsAndEnvelopes;
import tberg.murphy.fileio.f;

public class DatasetIO {
	
	public static class Datum implements java.io.Serializable {
		private static final long serialVersionUID = 1L;
		public final Wave wave;
		public final List<Event> events;
		public final String baseName;
		public Datum(Wave wave, List<Event> events, String baseName) {
			System.out.println("Length sec: "+wave.lengthSec);
			System.out.println("Sample rate: "+wave.sampleRateHz);
			System.out.println("Num samples: "+wave.amplitudes.length);
//			System.out.println("samples extrap length sec: "+wave.amplitudes.length * (1.0 / wave.sampleRateHz));
			this.wave = wave;
			this.events = events;
			this.baseName = baseName;
		}
	}
	
	public static List<Datum> readLabeledData(String path, float prefixLengthSec) {
		List<Datum> dataset = new ArrayList<Datum>();
		File baseNamesFile = new File(path);
		List<String> baseNames = null;
		String dirPath = null;
		if (baseNamesFile.isDirectory()) {
			baseNames = getLabeledBaseNamesFromDir(baseNamesFile);
			dirPath = path;
		} else {
			baseNames = getBaseNamesFromFile(baseNamesFile);
			dirPath = baseNamesFile.getParent();
		}
		for (String baseName : baseNames) {
			dataset.add(readLabeledDatum(dirPath, baseName, prefixLengthSec));
		}
		return dataset;
	}
	
	public static Datum readLabeledDatum(String dirPath, String baseName, float prefixLengthSec) {
      Wave wave = AudioIO.takePrefix(AudioIO.readWaveFile(dirPath+"/"+baseName+".wav"), prefixLengthSec);
      return new Datum(wave, PitchEventIO.takePrefix(PitchEventIO.readMIDIFile(dirPath+"/"+baseName+".mid"), wave.lengthSec), baseName);
	}
	
	public static Datum readLiveDatum(String filePath) {
		AudioCaptureIO.captureAudio(filePath);
		Wave wave = AudioIO.readWaveFile(filePath);
	    return new Datum(wave, null, filePath);
	}

    public static Datum readLiveDatum(String filePath, long captureDuration) {
        AudioCaptureIO.captureAudio(filePath, captureDuration);
        Wave wave = AudioIO.readWaveFile(filePath);
        return new Datum(wave, null, filePath);
    }
	
	public static Datum readUnlabeledDatum(String filePath) {
		Wave wave = AudioIO.readWaveFile(filePath);
	    return new Datum(wave, null, filePath);
	}
	
	private static List<String> getLabeledBaseNamesFromDir(File baseDir) {
		List<String> baseNames = new ArrayList<String>();
		for (File f : baseDir.listFiles()) {
			String fileName = f.getName();
			if (fileName.endsWith(".mid")) {
				baseNames.add(fileName.substring(0,fileName.length()-4));
			}
		}
		return baseNames;
	}
	
	private static List<String> getBaseNamesFromFile(File baseFile) {
		return f.readLines(baseFile.getAbsolutePath());
	}
	
	public static SpectralAtomsAndEnvelopes readSpectralAtomsAndEnvelopes(String path) {
		SpectralAtomsAndEnvelopes spectralAtomsAndEvelopes = null;
		try {
			FileInputStream fileIn = new FileInputStream(path);
			ObjectInputStream in = new ObjectInputStream(fileIn);
			spectralAtomsAndEvelopes = (SpectralAtomsAndEnvelopes) in.readObject();
			in.close();
			fileIn.close();
		} catch (IOException i) {
			i.printStackTrace();
		} catch (ClassNotFoundException c) {
			System.out.println("SpectralAtomsAndEnvelopes class not found");
			c.printStackTrace();
		}
		return spectralAtomsAndEvelopes;
	}
	
	public static void writeSpectralAtomsAndEnvelopes(SpectralAtomsAndEnvelopes spectralAtomsAndEnvelopes, String path) {
		try {
			FileOutputStream fileOut = new FileOutputStream(path);
			ObjectOutputStream out = new ObjectOutputStream(fileOut);
			out.writeObject(spectralAtomsAndEnvelopes);
			out.close();
			fileOut.close();
		} catch(IOException i) {
			i.printStackTrace();
		}
	}
	
}
