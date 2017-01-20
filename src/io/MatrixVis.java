package io;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.WritableRaster;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JScrollPane;

import java.awt.BorderLayout;

import javax.swing.JPanel;
import javax.swing.JButton;

import eval.PitchEventUtil.NoteState;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.File;
import java.io.IOException;
import java.util.concurrent.atomic.AtomicInteger;

import tberg.murphy.arrays.a;

public class MatrixVis {
	
//	public static int maxDim = 1920;
	public static int maxDim = 10000;
	
	public static void visualizeGrayscale(boolean[][] mat) {
		float[][] dmat = new float[mat.length][mat[0].length];
		for (int i=0; i<mat.length; ++i) {
			for (int j=0; j<mat[i].length; ++j) {
				dmat[i][j] = (mat[i][j] ? 1.0f : 0.0f);
			}
		}
		visualize(makeGrayscaleImage(dmat));
	}
	
	public static void visualizeGrayscale(float[] vect) {
		float[][] mat = new float[100][vect.length];
		for (int i=0; i<100; ++i) {
			mat[i] = vect;
		}
		visualize(makeGrayscaleImage(mat));
	}
	
	public static void visualizeGrayscale(NoteState[][] mat) {
		visualize(makeGrayscaleImage(stateTypeToFloatMatrix(mat)));
	}

	public static void visualizeGrayscale(float[][] mat) {
		visualize(makeGrayscaleImage(mat));
	}

    public static void visualizePair(float[][] m1, float[][] m2) {
        visualize(makePairImage(m1, m2));
    }
	
	public static void writeGrayscale(String path, float[][] mat) {
		writeBufferedImageHard(path, makeGrayscaleImage(mat));
	}
	
    public static float[][] stateTypeToFloatMatrix(NoteState[][] mat) {
        float[][] dmat = new float[mat.length][mat[0].length];
        for (int i=0; i<mat.length; ++i) {
            for (int j=0; j<mat[i].length; ++j) {
                dmat[i][j] = (float) mat[i][j].ordinal() / (float) NoteState.values().length-1;
            }
        }
        return dmat;
    }

    public static void visualizeOverlay(float[][] input, NoteState[][] pred, NoteState[][] gold) {
        visualizeOverlay(new float[input.length][input[0].length], input, pred, gold);
    }

    public static void visualizeOverlay(float[][] goldInput, float[][] input, NoteState[][] pred, NoteState[][] gold) {
        visualize(makeOverlayImage(goldInput, input, stateTypeToFloatMatrix(pred), stateTypeToFloatMatrix(gold)));
    }

    public static void writeOverlay(String path, float[][] input, NoteState[][] pred, NoteState[][] gold) {
        writeOverlay(path, new float[input.length][input[0].length], input, pred, gold);
    }

    public static void writeOverlay(String path, float[][] goldInput, float[][] input, NoteState[][] pred, NoteState[][] gold) {
    	writeBufferedImageHard(path, makeOverlayImage(goldInput, input, stateTypeToFloatMatrix(pred), stateTypeToFloatMatrix(gold)));
    }
    
    // public static void visualizeOverlay(float[][] input0, float[][] input1, NoteState[][] gold) {
    // 	visualize(makeOverlayImage(input0, input1, stateTypeToFloatMatrix(gold)));
    // }
    
    //public static void writeOverlay(String path, float[][] input0, float[][] input1, NoteState[][] pred, NoteState[][] gold) {
    //	writeBufferedImageHard(path, makeOverlayImage(input0, input1, stateTypeToFloatMatrix(gold)));
    //}

    //private static BufferedImage makeOverlayImage(float[][] input, float[][] pred, float[][] gold) {
    //    return makeOverlayImage(new float[input.length][input[0].length], input, pred, gold);
    //}

    private static BufferedImage makeOverlayImage(float[][] goldInput, float[][] input, float[][] pred, float[][] gold) {
        final int LIGHTEN = 15;
        final int SCALE = 255 - LIGHTEN;
        input = padAndNormalize(input);
        pred = padAndNormalize(pred);
        gold = padAndNormalize(gold);
        assert input.length == pred.length && input.length == gold.length && input[0].length == pred[0].length && input[0].length == gold[0].length;
        int cols = input.length;
        int rows = input[0].length;
        final BufferedImage img = new BufferedImage(Math.min(cols, maxDim), Math.min(rows * 5, maxDim), BufferedImage.TYPE_INT_RGB);
        WritableRaster raster = img.getRaster();
        for (int i = 0; i < raster.getWidth(); ++i) {
            for (int j = 0; j < raster.getHeight() / 5; ++j) {
                int inv_j = raster.getHeight() - 1 - 5 * j;
                float goldInputInt = goldInput[i][j] * SCALE + LIGHTEN;
                raster.setPixel(i, inv_j, new float[] { goldInputInt, goldInputInt, goldInputInt});
                float inputInt = input[i][j] * SCALE + LIGHTEN;
                raster.setPixel(i, inv_j - 1, new float[] { inputInt, inputInt, inputInt });
                float predInt = pred[i][j] * SCALE + LIGHTEN;
                raster.setPixel(i, inv_j - 2, new float[] { LIGHTEN, LIGHTEN, predInt });
                float goldInt = gold[i][j] * SCALE + LIGHTEN;
                raster.setPixel(i, inv_j - 3, new float[] { LIGHTEN, goldInt, LIGHTEN });
            }
        }
        return img;
    }
    
    public static BufferedImage makeGrayscaleImage(float[][] mat) {
        mat = padAndNormalize(mat);
		final BufferedImage img = new BufferedImage(Math.min(mat.length, maxDim), Math.min(mat[0].length, maxDim), BufferedImage.TYPE_BYTE_GRAY);
		WritableRaster writeableRaster = img.getRaster();
		for (int i=0; i<writeableRaster.getWidth(); ++i) {
			for (int j=0; j<writeableRaster.getHeight(); ++j) {
				writeableRaster.setPixel(i, writeableRaster.getHeight()-1-j, new float[] { mat[i][j] * 255 });
			}
		}
        return img;
    }

    private static BufferedImage makePairImage(float[][] m1, float[][] m2) {
        m1 = padAndNormalize(m1);
        m2 = padAndNormalize(m2);
        assert m1.length == m2.length;
        assert m1[0].length == m2[0].length;
        final BufferedImage img = new BufferedImage(m1.length * 2, m1[0].length, BufferedImage.TYPE_INT_RGB);
        WritableRaster writeableRaster = img.getRaster();
        for (int i = 0; i < writeableRaster.getWidth() / 2; ++i) {
            for (int j = 0; j < writeableRaster.getHeight(); ++j) {
                float v1 = m1[i][j] * 255;
                float v2 = m2[i][j] * 255;
                writeableRaster.setPixel(i*2, writeableRaster.getHeight()-1-j, new float[] { v1, v1/2, v1/2 });
                writeableRaster.setPixel(i*2 + 1, writeableRaster.getHeight()-1-j, new float[] { v2/2, v2, v2/2 });
            }
        }
        return img;
    }
    
    public static float[][] padMatrix(float[][] mat) {
    	int maxLength = 0;
    	for (int i=0; i<mat.length; ++i) {
    		if (mat[i].length > maxLength) {
    			maxLength = mat[i].length; 
    		}
    	}
    	float[][] result = new float[mat.length][maxLength];
    	for (int i=0; i<mat.length; ++i) {
    		System.arraycopy(mat[i], 0, result[i], 0, mat[i].length);
    	}
    	return result;
    }
    
    private static float[][] padAndNormalize(float[][] mat) {
        mat = padMatrix(mat);
        for (int i=0; i<mat.length; ++i) {
            for (int j=0; j<mat[i].length; ++j) {
                if (mat[i][j] == Double.NEGATIVE_INFINITY) mat[i][j] = 0;
            }
        }
        a.addi(mat, -a.min(a.min(mat)));
        a.scalei(mat, 1.0f / a.max(a.max(mat)));
        return mat;
    }
    
	private static void writeBufferedImageHard(String path, BufferedImage image) {
		try {
			File outputfile = new File(path);
			ImageIO.write(image, path.substring(path.lastIndexOf(".") + 1), outputfile);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	private static void visualize(final BufferedImage img) {
		final JFrame frame = new JFrame();
		frame.getContentPane().setLayout(new BorderLayout());
		
        final ImageIcon imageIcon = new ImageIcon(img);
        JLabel imageLabel = new JLabel(imageIcon);
        frame.getContentPane().add(new JScrollPane(imageLabel), BorderLayout.CENTER);
		//frame.getContentPane().add(new JLabel(new ImageIcon(img)));

        JPanel buttonPanel = new JPanel();
        buttonPanel.setLayout(new GridLayout(1,4));
        JButton zoomInHButton = new JButton("+H");
        JButton zoomOutHButton = new JButton("-H");
        JButton zoomInVButton = new JButton("+V");
        JButton zoomOutVButton = new JButton("-V");
        buttonPanel.add(zoomInHButton);
        buttonPanel.add(zoomOutHButton);
        buttonPanel.add(zoomInVButton);
        buttonPanel.add(zoomOutVButton);

        final AtomicInteger zoomX = new AtomicInteger(0);
        final AtomicInteger zoomY = new AtomicInteger(0);

        zoomInHButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                zoomX.set(zoomX.get() + 1);
                refreshViewer(img, zoomX, zoomY, imageIcon, frame);
            }
        });
        zoomOutHButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                zoomX.set(zoomX.get() - 1);
                refreshViewer(img, zoomX, zoomY, imageIcon, frame);
            }
        });
        zoomInVButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                zoomY.set(zoomY.get() + 1);
                refreshViewer(img, zoomX, zoomY, imageIcon, frame);
            }
        });
        zoomOutVButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                zoomY.set(zoomY.get() - 1);
                refreshViewer(img, zoomX, zoomY, imageIcon, frame);
            }
        });

        frame.getContentPane().add(buttonPanel, BorderLayout.SOUTH);

		frame.pack();
		frame.setVisible(true);
        frame.invalidate();
	}

	private static void refreshViewer(Image img, AtomicInteger zoomX, AtomicInteger zoomY, ImageIcon icon, Frame frame) {
        //System.err.println(zoomX);
        //System.err.println(zoomY);
        Image newImage;
        if (zoomX.get() == 1 && zoomY.get() == 1) {
            newImage = img;
        } else {
            newImage = img.getScaledInstance((int)(img.getWidth(frame) * Math.pow(2, zoomX.get())),
                                             (int)(img.getHeight(frame) * Math.pow(2, zoomY.get())),
                                             Image.SCALE_SMOOTH);
        }
        icon.setImage(newImage);
        frame.repaint();
    }
	
	public static void main(String[] args) {
		int n = 1000;
		float[][] image = new float[n][n];
		for (int i=0; i<n; ++i) {
			for (int j=0; j<n; ++j) {
				image[i][j] = (float) Math.exp(-(1.0/(100.0*n))*((i - n/2.0)*(i - n/2.0) + (j - n/2.0)*(j - n/2.0)));
			}
		}
		visualizeGrayscale(image);
	}

}
