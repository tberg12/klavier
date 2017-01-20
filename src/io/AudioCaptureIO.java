package io;

import javax.sound.sampled.*;

import java.io.*;
 
/**
 * A sample program is to demonstrate how to record sound in Java
 * author: www.codejava.net
 */
public class AudioCaptureIO {
	
	public static final float SAMPLE_RATE = 44100;
	public static final int SAMPLE_SIZE_BITS = 16;
	public static final int CHANNELS = 2;
	public static final boolean SIGNED = true;
	public static final boolean BIG_ENDIAN = true;
	public static final AudioFileFormat.Type FILE_TYPE = AudioFileFormat.Type.WAVE;
	
	public static class AudioCapture {
		TargetDataLine line;

		void start(String filePath) {
			File wavFile = new File(filePath);
			try {
				AudioFormat format = new AudioFormat(SAMPLE_RATE, SAMPLE_SIZE_BITS, CHANNELS, SIGNED, BIG_ENDIAN);
				DataLine.Info info = new DataLine.Info(TargetDataLine.class, format);

				// checks if system supports the data line
				if (!AudioSystem.isLineSupported(info)) {
					System.out.println("Line not supported");
					System.exit(0);
				}
				line = (TargetDataLine) AudioSystem.getLine(info);
				line.open(format);
				line.start();   // start capturing

				AudioInputStream ais = new AudioInputStream(line);

				System.out.println("Recording...");

				// start recording
				AudioSystem.write(ais, FILE_TYPE, wavFile);

			} catch (LineUnavailableException ex) {
				ex.printStackTrace();
			} catch (IOException ioe) {
				ioe.printStackTrace();
			}
		}

		void finish() {
			line.stop();
			line.close();
			System.out.println("Finished");
		}
	}

	public static void captureAudio(String filePath, final long durationMs) {
        final AudioCapture recorder = new AudioCapture();
 
        // creates a new thread that waits for a specified
        // of time before stopping
        Thread stopper = new Thread(new Runnable() {
            public void run() {
                try {
                    Thread.sleep(durationMs);
                } catch (InterruptedException ex) {
                    ex.printStackTrace();
                }
                recorder.finish();
            }
        });
 
        stopper.start();
 
        // start recording
        recorder.start(filePath);
    }
	
	public static void captureAudio(String filePath) {
		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		try {
			System.out.println("Press return to begin recording.");
			in.readLine();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		final AudioCapture recorder = new AudioCapture();
		
		// creates a new thread that waits for a specified
		// of time before stopping
		Thread stopper = new Thread(new Runnable() {
			public void run() {
				try {
					BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
					in.readLine();
					in.close();
				} catch (Exception ex) {
					ex.printStackTrace();
				}
				recorder.finish();
			}
		});
		
		stopper.start();
		
		// start recording
		recorder.start(filePath);
	}

    public static void main(String[] args) {
//    	captureAudio("/Users/tberg/Desktop/test.wav", 5000);
    	captureAudio("/Users/tberg/Desktop/test.wav");
    }
    
}