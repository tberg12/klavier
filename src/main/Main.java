package main;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

import io.*;
import io.DatasetIO.Datum;
import learning.Model;
import learning.Preprocessing;
import learning.Preprocessing.ProcDatum;
import learning.Preprocessing.SpectralAtomsAndEnvelopes;
import dsp.STFT.WindowType;
import eval.Evaluation;
import eval.Evaluation.EvalSuffStats;
import eval.PitchEventUtil;
import tberg.murphy.fig.Option;
import tberg.murphy.fig.OptionsParser;
import tberg.murphy.fileio.f;
import tberg.murphy.gpu.JOCLBlasUtil;
import io.PitchEventIO.Event;

public class Main implements Runnable {
	
	public static enum InputType {LIVE, FILE, DATASET}
	
	@Option(gloss = "")
	public static InputType inputType = InputType.FILE;
	@Option(gloss = "")
	public static boolean train = false;
	@Option(gloss = "")
	public static boolean evaluate = false;
	@Option(gloss = "")
	public static boolean buildInitializer = false;
	@Option(gloss = "")
	public static boolean writeInitializer = false;
	
    @Option(gloss = "")
    public static int filterNoteMinIndex = 0;
    @Option(gloss = "")
//    public static int filterNoteMaxIndex = 62;
    public static int filterNoteMaxIndex = 108;
    @Option(gloss = "")
    public static int filterNoteMinVel = 20;
//    public static int filterNoteMinVel = 0;
	
	@Option(gloss = "")
	public static String audioFilePath = "/Users/tberg/Dropbox/demos/klavier/gould2.wav";
//	public static String audioFilePath = "./piano_demo/wav/MAPS_MUS-deb_clai_ENSTDkCl.wav";
	@Option(gloss = "")
	public static String datasetPaths = "/Users/tberg/Dropbox/corpora/music_data/maps/benetos_test.txt";
	@Option(gloss = "")
	public static float datasetPrefixLengthSec = 30;
	@Option(gloss = "")
	public static String outputDir = "/Users/tberg/Desktop/piano_test/";
	@Option(gloss = "")
	public static String initSpectEnvSer = "/Users/tberg/Dropbox/demos/klavier/piano_demo/initSpectEnvFull.ser";
	@Option(gloss = "")
//	public static String spectEnvSer = "./piano_demo/learnedSpecCNMAT2.ser";
//	public static String spectEnvSer = "./piano_demo/learnedSpecCNMAT.ser";
//	public static String spectEnvSer = "./piano_demo/learnedSpecEnvGould.ser";
//	public static String spectEnvSer = "./piano_demo/learnedSpecEnvFullDemo1.ser";
//	public static String spectEnvSer = "./piano_demo/learnedSpecEnvMiddleCPiano12.ser";
//	public static String spectEnvSer = "./piano_demo/learnedSpecEnvFullMaps.ser";
//	public static String spectEnvSer = "/Users/tberg/Dropbox/demos/klavier/piano_demo/learnedSpecEnvGould2Itunes.ser";
	public static String spectEnvSer = "/Users/tberg/Desktop/piano_test/learnedSpecEnvGould2Itunes.ser";
    @Option(gloss = "")
    public static String imslpTransModelSer = "/Users/tberg/Dropbox/demos/klavier/data/imslp_indep_transitions.ser";

    @Option(gloss = "")
    public static int updateActivationsIters = 800;
//    public static int updateActivationsIters = 400;
    @Option(gloss = "")
    public static double updateActivationsStepSize = 5e-1;
//    public static double updateActivationsStepSize = 1e0;
    @Option(gloss = "")
    public static double activationsPriorWeight = 5e0;
//    public static double activationsPriorWeight = 5e0;
    
    @Option(gloss = "")
    public static int updateAtomsIters = 400;
    @Option(gloss = "")
    public static double updateAtomsStepSize = 1e-1;
    @Option(gloss = "")
    public static double atomsPriorWeight = 2e-3;
    
    
	@Option(gloss = "")
	public static String initSpectEnvInstrNamesPath = "/Users/tberg/Dropbox/corpora/music_data/maps/mapsm_instrument_names.txt";
	@Option(gloss = "")
	public static String initSpectEnvNotesBasePath = "/Users/tberg/Dropbox/corpora/music_data/maps/";
	@Option(gloss = "")
	public static String initSpectEnvNotesPaths = "ISOL/NO/";
	@Option(gloss = "")
	public static double initSpectEnvNotesOffsetMs = -40;

	@Option(gloss = "")
	public static double modelActivationVar = 1e-3;
	@Option(gloss = "")
	public static int modelNumEMIters = 5;
	@Option(gloss = "")
	public static int modelNumApproxIters = 0
	;
	@Option(gloss = "")
	public static boolean modelEMUpdateEnvelopes = true;
	@Option(gloss = "")
	public static boolean modelEMUpdateAtoms = true;
	@Option(gloss = "")
	public static double modelMinNoteLengthMs = 200;
	@Option(gloss = "")
	public static double modelMaxNoteLengthMs = 2000;
	
	@Option(gloss = "")
	public static int stftWindowSize = 4096;
	@Option(gloss = "")
	public static int stftFrequencyPrefixSize = 512;
	@Option(gloss = "")
	public static int stftFrequencyOffset = 0;
	@Option(gloss = "")
	public static int stftHopSize = 512;
	@Option(gloss = "")
	public static WindowType stftWindowType = WindowType.HANNINGSQRT;
	@Option(gloss = "")
	public static int stftDownsampleSpectrogramHopSize = 2;
	@Option(gloss = "")
	public static NMFType nmfType = NMFType.Beta;
	@Option(gloss = "")
	public static double nmfBeta = 0.5;
	@Option(gloss = "")
	public static double nmfC = -0.5;
	@Option(gloss = "")
	public static double nmfSilenceEps = 1e-10;
	@Option(gloss = "")
	public static double nmfMinEps = 1e-10;
	@Option(gloss = "")
	public static int nmfNumIters = 800;
	@Option(gloss = "")
	public static boolean useGpu = true;
//	@Option(gloss = "")
//	public static int nmfGpuId = 0;
	@Option(gloss = "")
	public static int numThreads = 8;

	public static enum NMFType {KL, L2, Beta, LogNormal};

	public void run() {
		if (useGpu) JOCLBlasUtil.startup();
		
		EvalSuffStats totalEvalOnsets = new EvalSuffStats(0, 0, 0);
		EvalSuffStats totalEvalFrames = new EvalSuffStats(0, 0, 0);
		
		(new File(outputDir)).mkdirs();
		
		if (train) {
			List<Datum> inputData = new ArrayList<Datum>();
			if (inputType == InputType.LIVE) {
				inputData.add(DatasetIO.readLiveDatum(audioFilePath));
			} else if (inputType == InputType.FILE) {
				inputData.add(DatasetIO.readUnlabeledDatum(audioFilePath));
			} else {
				inputData = DatasetIO.readLabeledData(datasetPaths, datasetPrefixLengthSec);
			}
			
			SpectralAtomsAndEnvelopes spectralAtomsAndEnvelopes = null;
			if (buildInitializer) {
				List<String> spectralAtomsAndEnvelopesInstrNames = f.readLines(initSpectEnvInstrNamesPath);
				spectralAtomsAndEnvelopes = Preprocessing.buildInitialSpectralAtomsAndEnvelopes(initSpectEnvNotesBasePath, initSpectEnvNotesPaths, spectralAtomsAndEnvelopesInstrNames);
				if (writeInitializer) DatasetIO.writeSpectralAtomsAndEnvelopes(spectralAtomsAndEnvelopes, initSpectEnvSer);
			} else {
				spectralAtomsAndEnvelopes = DatasetIO.readSpectralAtomsAndEnvelopes(initSpectEnvSer);
			}
			List<ProcDatum> testSet = Preprocessing.preprocessData(inputData, spectralAtomsAndEnvelopes);

			System.out.println("Milliseconds per frame: " + String.format("%.1f", testSet.get(0).secondsPerFrame * 1e3));
			double numTestFrames = 0;
			for (ProcDatum datum : testSet) numTestFrames += datum.activations.length;
			System.out.println("Num test frames: "+String.format("%.0f", numTestFrames));

			Model model = new Model(spectralAtomsAndEnvelopes);

			List<List<Event>> decode = model.train(testSet, modelEMUpdateAtoms, modelEMUpdateEnvelopes, modelNumEMIters, modelNumApproxIters);
			List<String> midiFilePaths = new ArrayList<String>();
			for (int i=0; i<testSet.size(); ++i) {
				ProcDatum datum = testSet.get(i);
				List<Event> guess = decode.get(i);
				guess = PitchEventUtil.filter(guess, filterNoteMinIndex, filterNoteMaxIndex, filterNoteMinVel);
//				System.out.println(guess.toString());
				try {
					String midiFilePath = outputDir + "/" + datum.rawDatum.baseName.split("/")[datum.rawDatum.baseName.split("/").length-1]+".mid";
					PitchEventIO.writeMIDIFile(guess, midiFilePath);
					midiFilePaths.add(midiFilePath);
				} catch (Exception e) {}
				
				if (evaluate) {
					EvalSuffStats evalOnsets = Evaluation.evaluateOnsets(guess, datum.rawDatum.events, datasetPrefixLengthSec, 100.0f);
					EvalSuffStats evalFrames = Evaluation.evaluateFrames(guess, datum.rawDatum.events, datasetPrefixLengthSec);
					System.out.println(datum.rawDatum.baseName);
					System.out.println("Eval onsets: " + evalOnsets);
					System.out.println("Eval frames: " + evalFrames);
					totalEvalOnsets.increment(evalOnsets);
					totalEvalFrames.increment(evalFrames);
				}
	            
			}
			
			for (String midiFilePath : midiFilePaths) {
//				String command = "/usr/bin/midisheetmusic.mono.exe "+midiFilePath;
//				String command = "/Applications/MidiSheetMusic-2.6.app/Contents/MacOS/MidiSheetMusic "+midiFilePath;
				String command = "open "+midiFilePath;
//				System.out.println("Executing\n"+command);
				try {
					Runtime.getRuntime().exec(command);
//					Process process = Runtime.getRuntime().exec(command);
//					Thread.sleep(5000);
//					process.destroy();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			
			DatasetIO.writeSpectralAtomsAndEnvelopes(model.getSpectralAtomsAndEnvelopes(), spectEnvSer);
		} else {
			List<Datum> inputData = new ArrayList<Datum>();
			if (inputType == InputType.LIVE) {
				inputData.add(DatasetIO.readLiveDatum(audioFilePath));
			} else if (inputType == InputType.FILE) {
				inputData.add(DatasetIO.readUnlabeledDatum(audioFilePath));
			} else {
				inputData = DatasetIO.readLabeledData(datasetPaths, datasetPrefixLengthSec);
			}
			
			SpectralAtomsAndEnvelopes spectralAtomsAndEnvelopes = DatasetIO.readSpectralAtomsAndEnvelopes(spectEnvSer);
			List<ProcDatum> testSet = Preprocessing.preprocessData(inputData, spectralAtomsAndEnvelopes);

			System.out.println("Milliseconds per frame: " + String.format("%.1f", testSet.get(0).secondsPerFrame * 1e3));
			double numTestFrames = 0;
			for (ProcDatum datum : testSet) numTestFrames += datum.activations.length;
			System.out.println("Num test frames: "+String.format("%.0f", numTestFrames));

			Model model = new Model(spectralAtomsAndEnvelopes);

			List<List<Event>> decode = model.decode(testSet, modelNumApproxIters);
			List<String> midiFilePaths = new ArrayList<String>();
			for (int i=0; i<testSet.size(); ++i) {
				ProcDatum datum = testSet.get(i);
				List<Event> guess = decode.get(i);
				guess = PitchEventUtil.filter(guess, filterNoteMinIndex, filterNoteMaxIndex, filterNoteMinVel);
//				System.out.println(guess.toString());
				try {
					String midiFilePath = outputDir + "/" + datum.rawDatum.baseName.split("/")[datum.rawDatum.baseName.split("/").length-1]+".mid";
					PitchEventIO.writeMIDIFile(guess, midiFilePath);
					midiFilePaths.add(midiFilePath);
				} catch (Exception e) {}
				
				if (evaluate) {
					EvalSuffStats evalOnsets = Evaluation.evaluateOnsets(guess, datum.rawDatum.events, datasetPrefixLengthSec, 100.0f);
					EvalSuffStats evalFrames = Evaluation.evaluateFrames(guess, datum.rawDatum.events, datasetPrefixLengthSec);
					System.out.println(datum.rawDatum.baseName);
					System.out.println("Eval onsets: " + evalOnsets);
					System.out.println("Eval frames: " + evalFrames);
					totalEvalOnsets.increment(evalOnsets);
					totalEvalFrames.increment(evalFrames);
				}
	            
			}
			
			for (String midiFilePath : midiFilePaths) {
//				String command = "/usr/bin/midisheetmusic.mono.exe "+midiFilePath;
//				String command = "/Applications/MidiSheetMusic-2.6.app/Contents/MacOS/MidiSheetMusic "+midiFilePath;
				String command = "open "+midiFilePath;
//				System.out.println("Executing\n"+command);
				try {
					Runtime.getRuntime().exec(command);
//					Process process = Runtime.getRuntime().exec(command);
//					Thread.sleep(5000);
//					process.destroy();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			
		}
		
		if (evaluate) {
			System.out.println("Total eval onsets: "+totalEvalOnsets);
			System.out.println("Total eval frames: "+totalEvalFrames);
		}
			
		if (useGpu) JOCLBlasUtil.shutdown();
	}
	
    private static byte[] toByteArray(InputStream inputStream) throws IOException {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        byte buffer[] = new byte[8192];
        while (true)
        {
            int read = inputStream.read(buffer);
            if (read == -1)
            {
                break;
            }
            baos.write(buffer, 0, read);
        }
        return baos.toByteArray();
    }
	
	public static void main(String[] args) {
		Main main = new Main();
		OptionsParser parser = new OptionsParser();
		parser.doRegisterAll(new Object[] {main});
		if (!parser.doParse(args)) System.exit(1);
		main.run();
	}

}
