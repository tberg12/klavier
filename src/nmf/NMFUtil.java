package nmf;

import tberg.murphy.gpu.CublasUtil;
import tberg.murphy.gpu.CublasUtil.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import tberg.murphy.opt.DifferentiableFunction;
import tberg.murphy.opt.EmpiricalGradientTester;

import org.jblas.FloatMatrix;
import org.jblas.MatrixFunctions;

import tberg.murphy.tuple.Pair;
import tberg.murphy.arrays.a;

public class NMFUtil {

	private static final int PRINT_PERIOD = 50;
	private static final int RAND_SEED = 4;
	private static final float MIN_LOG_EPS = 1e-10f;

	public static void gradientTester(float[][] Wtrans0, float[][] Htrans0, float[][] Xtrans0, float silenceEps) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			Htrans0 = filterSilence(Htrans0, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		final int n = Xtrans0[0].length;
		final int m = Xtrans0.length;
		final int r = Wtrans0.length;

		System.out.println("Gradient Test GPU: ");

		final Matrix X = Matrix.build(a.transpose(Xtrans0));
		float[] flattenedWinit = Matrix.build(a.transpose(Wtrans0)).toArray();
		float[] flattenedHinit = Matrix.build(a.transpose(Htrans0)).toArray();

		DifferentiableFunction obj = new DifferentiableFunction() {
			public Pair<Double, double[]> calculate(double[] xDouble) {
				float[] x = a.toFloat(xDouble);
				Matrix W = Matrix.build(n, r, Arrays.copyOfRange(x, 0, n * r));
				Matrix H = Matrix.build(r, m, Arrays.copyOfRange(x, n * r, n * r + r * m));
				// Matrix Htrans = H.transpose();
				// Matrix Wtrans = W.transpose();
				Matrix WH = W.mmul(H);
				// Matrix logX = X.log();
				// Matrix logWH = WH.log();
				// float c = 1.0f;
				// Matrix logXsubLogWHdivWH = logX.sub(logWH).add(Matrix.ones(n,
				// m).mul(c)).div(WH);

				float beta = 0.5f;
				Matrix WHbetaSub1 = WH.pow(beta - 1.0f);
				Matrix WHbetaSub2 = WH.pow(beta - 2.0f);
				Matrix gH = W.transpose().mmul(WHbetaSub1.sub(X.mul(WHbetaSub2)));
				Matrix gW = (WHbetaSub1.sub(X.mul(WHbetaSub2))).mmul(H.transpose());
				double val = (1.0f / (beta * (beta - 1.0f))) * (X.pow(beta).norm1()
						+ (beta - 1.0f) * WH.pow(beta).norm1() - beta * (X.mul(WH.pow(beta - 1.0f)).norm1()));

				// Matrix gW = logXsubLogWHdivWH.mmul(Htrans).mul(-2.0f);
				// Matrix gH = Wtrans.mmul(logXsubLogWHdivWH).mul(-2.0f);
				// double val = Math.pow(logX.sub(logWH).add(Matrix.ones(n,
				// m).mul(c)).norm2(), 2.0);

				// Matrix gW = (X.mmul(Htrans)).comb(-2.0f, 2.0f,
				// (WH).mmul(Htrans));
				// Matrix gH = (Wtrans.mmul(X)).comb(-2.0f, 2.0f,
				// Wtrans.mmul(WH));
				// double val = Math.pow(X.distance2(WH), 2.0);
				Pair<Double, double[]> result = Pair.makePair(val, a.toDouble(a.append(gW.toArray(), gH.toArray())));
				CublasUtil.freeAllBut(X);
				return result;
			}
		};

		EmpiricalGradientTester.test(obj, a.toDouble(a.append(flattenedWinit, flattenedHinit)), 1e-2, 1.0, 1e-7);
	}

	public static Pair<float[][], float[][]> nmfKLL2PriorExpGradGPU(float[][] priorWtrans, float[][] initWtrans,
			final float priorWeightW, boolean updateW, float[][] priorHtrans, float[][] initHtrans,
			final float priorWeightH, boolean updateH, float[][] Xtrans0, float startStepSize, float endStepSize,
			int iters, float silenceEps, float minEps) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			priorHtrans = filterSilence(priorHtrans, loud);
			initHtrans = filterSilence(initHtrans, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		final int n = Xtrans0[0].length;
		final int m = Xtrans0.length;

		System.out.println("NMF KL Exp Grad GPU: ");

		long start = System.nanoTime();

		Matrix X = Matrix.build(a.transpose(Xtrans0));
		Matrix logX = X.log();
		Matrix Wprior = Matrix.build(a.transpose(priorWtrans));
		Matrix Hprior = Matrix.build(a.transpose(priorHtrans));
		Matrix W = Matrix.build(a.transpose(initWtrans));
		Matrix H = Matrix.build(a.transpose(initHtrans));

		Matrix denomUpdateW = Matrix.ones(n, 1).mmul(Matrix.ones(1, m).mmul(H.transpose()));
		Matrix denomUpdateH = (W.transpose().mmul(Matrix.ones(n, 1))).mmul(Matrix.ones(1, m));
		{
			Matrix WH = W.mmul(H);
			System.out.println(((X.mul(logX.sub(WH.log())).sub(X)).add(WH)).norm1()
					+ (priorWeightW > 0.0 ? priorWeightW * Math.pow(W.distance2(Wprior), 2.0) : 0.0)
					+ (priorWeightH > 0.0 ? priorWeightH * Math.pow(H.distance2(Hprior), 2.0) : 0.0));
		}
		for (int i = 0; i < iters; ++i) {
			float stepSize = startStepSize + ((float) i / iters) * (endStepSize - startStepSize);
			if (updateH) {
				Matrix priorHGrad = H.sub(Hprior).mul(priorWeightH * 2.0f);
				Matrix numerUpdateH = W.transpose().mmul(X.div(W.mmul(H)));
				Matrix gradH = (denomUpdateH.sub(numerUpdateH)).add(priorHGrad);
				H = H.mul((gradH.mul(-stepSize)).exp());
				H.maxi(minEps);
				denomUpdateW = Matrix.ones(n, 1).mmul(Matrix.ones(1, m).mmul(H.transpose()));
			}
			if (updateW) {
				Matrix priorWGrad = W.sub(Wprior).mul(priorWeightW * 2.0f);
				Matrix numerUpdateW = (X.div(W.mmul(H))).mmul(H.transpose());
				Matrix gradW = (denomUpdateW.sub(numerUpdateW)).add(priorWGrad);
				W = W.mul((gradW.mul(-stepSize)).exp());
				W.maxi(minEps);
				denomUpdateH = (W.transpose().mmul(Matrix.ones(n, 1))).mmul(Matrix.ones(1, m));
			}
			if (i % PRINT_PERIOD == 0) {
				Matrix WH = W.mmul(H);
				System.out.println(((X.mul(logX.sub(WH.log())).sub(X)).add(WH)).norm1()
						+ (priorWeightW > 0.0 ? priorWeightW * Math.pow(W.distance2(Wprior), 2.0) : 0.0)
						+ (priorWeightH > 0.0 ? priorWeightH * Math.pow(H.distance2(Hprior), 2.0) : 0.0));
			}
			CublasUtil.freeAllBut(W, Wprior, H, Hprior, X, logX, denomUpdateW, denomUpdateH);
		}

		float[][] resultWtrans = W.transpose().toArray2();
		float[][] resultHtrans = H.transpose().toArray2();
		CublasUtil.freeAll();

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(resultWtrans, insertSilence(resultHtrans, loud));
		} else {
			return Pair.makePair(resultWtrans, resultHtrans);
		}
	}

	public static Pair<float[][], float[][]> nmfKLL2PriorExpGrad(float[][] priorWtrans, float[][] initWtrans,
			final float priorWeightW, boolean updateW, float[][] priorHtrans, float[][] initHtrans,
			final float priorWeightH, boolean updateH, float[][] Xtrans0, float startStepSize, float endStepSize,
			int iters, float silenceEps, float minEps) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			priorHtrans = filterSilence(priorHtrans, loud);
			initHtrans = filterSilence(initHtrans, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		final int n = Xtrans0[0].length;
		final int m = Xtrans0.length;

		System.out.println("NMF KL Exp Grad: ");

		long start = System.nanoTime();

		FloatMatrix X = new FloatMatrix(a.transpose(Xtrans0));
		FloatMatrix logX = MatrixFunctions.log(X);
		FloatMatrix Wprior = new FloatMatrix(a.transpose(priorWtrans));
		FloatMatrix Hprior = new FloatMatrix(a.transpose(priorHtrans));
		FloatMatrix W = new FloatMatrix(a.transpose(initWtrans));
		FloatMatrix H = new FloatMatrix(a.transpose(initHtrans));

		FloatMatrix denomUpdateW = FloatMatrix.ones(n, 1).mmul(FloatMatrix.ones(1, m).mmul(H.transpose()));
		FloatMatrix denomUpdateH = (W.transpose().mmul(FloatMatrix.ones(n, 1))).mmul(FloatMatrix.ones(1, m));
		{
			FloatMatrix WH = W.mmul(H);
			System.out.println(((X.mul(logX.sub(MatrixFunctions.log(WH))).sub(X)).add(WH)).norm1()
					+ (priorWeightW > 0.0 ? priorWeightW * Math.pow(W.distance2(Wprior), 2.0) : 0.0)
					+ (priorWeightH > 0.0 ? priorWeightH * Math.pow(H.distance2(Hprior), 2.0) : 0.0));
		}
		for (int i = 0; i < iters; ++i) {
			float stepSize = startStepSize + ((float) i / iters) * (endStepSize - startStepSize);
			if (updateH) {
				FloatMatrix priorHGrad = H.sub(Hprior).mul(priorWeightH * 2.0f);
				FloatMatrix numerUpdateH = W.transpose().mmul(X.div(W.mmul(H)));
				FloatMatrix gradH = (denomUpdateH.sub(numerUpdateH)).add(priorHGrad);
				H = H.mul(MatrixFunctions.exp(gradH.mul(-stepSize)));
				H.maxi(minEps);
				denomUpdateW = FloatMatrix.ones(n, 1).mmul(FloatMatrix.ones(1, m).mmul(H.transpose()));
			}
			if (updateW) {
				FloatMatrix priorWGrad = W.sub(Wprior).mul(priorWeightW * 2.0f);
				FloatMatrix numerUpdateW = (X.div(W.mmul(H))).mmul(H.transpose());
				FloatMatrix gradW = (denomUpdateW.sub(numerUpdateW)).add(priorWGrad);
				W = W.mul(MatrixFunctions.exp(gradW.mul(-stepSize)));
				W.maxi(minEps);
				denomUpdateH = (W.transpose().mmul(FloatMatrix.ones(n, 1))).mmul(FloatMatrix.ones(1, m));
			}
			if (i % PRINT_PERIOD == 0) {
				FloatMatrix WH = W.mmul(H);
				System.out.println(((X.mul(logX.sub(MatrixFunctions.log(WH))).sub(X)).add(WH)).norm1()
						+ (priorWeightW > 0.0 ? priorWeightW * Math.pow(W.distance2(Wprior), 2.0) : 0.0)
						+ (priorWeightH > 0.0 ? priorWeightH * Math.pow(H.distance2(Hprior), 2.0) : 0.0));
			}
		}

		float[][] resultWtrans = W.transpose().toArray2();
		float[][] resultHtrans = H.transpose().toArray2();

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(resultWtrans, insertSilence(resultHtrans, loud));
		} else {
			return Pair.makePair(resultWtrans, resultHtrans);
		}
	}

	public static Pair<float[][], float[][]> nmfBeta(float[][] Xtrans0, int r, int iters, float silenceEps,
			float minEps, float beta) {
		int n = Xtrans0[0].length;
		int m = Xtrans0.length;
		Random rand = new Random(RAND_SEED);
		return nmfBeta(a.randFloat(r, n, rand), true, a.randFloat(m, r, rand), true, Xtrans0, iters, silenceEps, minEps,
				beta);
	}

	public static Pair<float[][], float[][]> nmfBeta(float[][] Wtrans0, boolean updateW, float[][] Htrans0,
			boolean updateH, float[][] Xtrans0, int iters, float silenceEps, float minEps, float beta) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			Htrans0 = filterSilence(Htrans0, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		System.out.println("NMF Beta: ");

		long start = System.nanoTime();

		FloatMatrix W = new FloatMatrix(a.transpose(Wtrans0));
		FloatMatrix H = new FloatMatrix(a.transpose(Htrans0));
		FloatMatrix X = new FloatMatrix(a.transpose(Xtrans0));

		{
			FloatMatrix WH = W.mmul(H);
			float divergence = 0.0f;
			if (beta == 0.0f) {
				FloatMatrix XdivWH = X.div(WH);
				divergence = XdivWH.sub(MatrixFunctions.log(XdivWH)).norm1() - 1.0f;
			} else if (beta == 1.0f) {
				divergence = ((X.mul((MatrixFunctions.log(X)).sub(MatrixFunctions.log(WH)))).add(WH.sub(X))).norm1();
			} else {
				divergence = (1.0f / (beta * (beta - 1.0f)))
						* (MatrixFunctions.pow(X, beta).norm1() + (beta - 1.0f) * MatrixFunctions.pow(WH, beta).norm1()
								- beta * (X.mul(MatrixFunctions.pow(WH, beta - 1.0f)).norm1()));
			}
			System.out.println(divergence);
		}
		for (int i = 0; i < iters; ++i) {
			if (updateH) {
				FloatMatrix WH = W.mmul(H);
				FloatMatrix WHbetaSub2 = MatrixFunctions.pow(WH, beta - 2.0f);
				FloatMatrix WHbetaSub1 = MatrixFunctions.pow(WH, beta - 1.0f);
				H = H.mul(W.transpose().mmul(WHbetaSub2.mul(X))).div(W.transpose().mmul(WHbetaSub1));
				H.maxi(minEps);
			}
			if (updateW) {
				FloatMatrix WH = W.mmul(H);
				FloatMatrix WHbetaSub2 = MatrixFunctions.pow(WH, beta - 2.0f);
				FloatMatrix WHbetaSub1 = MatrixFunctions.pow(WH, beta - 1.0f);
				W = W.mul((X.mul(WHbetaSub2)).mmul(H.transpose())).div((WHbetaSub1).mmul(H.transpose()));
				W.maxi(minEps);
			}
			if (i % PRINT_PERIOD == 0) {
				FloatMatrix WH = W.mmul(H);
				float divergence = 0.0f;
				if (beta == 0.0f) {
					FloatMatrix XdivWH = X.div(WH);
					divergence = XdivWH.sub(MatrixFunctions.log(XdivWH)).norm1() - 1.0f;
				} else if (beta == 1.0f) {
					divergence = ((X.mul((MatrixFunctions.log(X)).sub(MatrixFunctions.log(WH)))).add(WH.sub(X)))
							.norm1();
				} else {
					divergence = (1.0f / (beta * (beta - 1.0f))) * (MatrixFunctions.pow(X, beta).norm1()
							+ (beta - 1.0f) * MatrixFunctions.pow(WH, beta).norm1()
							- beta * (X.mul(MatrixFunctions.pow(WH, beta - 1.0f)).norm1()));
				}
				System.out.println(divergence);
			}
		}

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(W.transpose().toArray2(), insertSilence(H.transpose().toArray2(), loud));
		} else {
			return Pair.makePair(W.transpose().toArray2(), H.transpose().toArray2());
		}
	}

	public static Pair<float[][], float[][]> nmfBetaGPU(float[][] Xtrans0, int r, int iters, float silenceEps,
			float minEps, float beta) {
		int n = Xtrans0[0].length;
		int m = Xtrans0.length;
		Random rand = new Random(RAND_SEED);
		return nmfBetaGPU(a.randFloat(r, n, rand), true, a.randFloat(m, r, rand), true, Xtrans0, iters, silenceEps,
				minEps, beta);
	}

	public static Pair<float[][], float[][]> nmfBetaGPU(float[][] Wtrans0, boolean updateW, float[][] Htrans0,
			boolean updateH, float[][] Xtrans0, int iters, float silenceEps, float minEps, float beta) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			Htrans0 = filterSilence(Htrans0, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		// System.out.println("NMF Beta GPU: ");
		System.out.print("Fitting model.");

		long start = System.nanoTime();

		Matrix W = Matrix.build(a.transpose(Wtrans0));
		Matrix H = Matrix.build(a.transpose(Htrans0));
		Matrix X = Matrix.build(a.transpose(Xtrans0));

		{
			// Matrix WH = W.mmul(H);
			// float divergence = 0.0f;
			// if (beta == 0.0f) {
			// Matrix XdivWH = X.div(WH);
			// divergence = (XdivWH.sub(XdivWH.log())).norm1() - 1.0f;
			// } else if (beta == 1.0f) {
			// divergence =
			// (X.mul(X.log().sub(WH.log())).add(WH.sub(X))).norm1();
			// } else {
			// divergence = (1.0f / (beta * (beta - 1.0f))) *
			// (X.pow(beta).norm1() + (beta - 1.0f) * WH.pow(beta).norm1() -
			// beta * (X.mul(WH.pow(beta-1.0f)).norm1()));
			// }
			// System.out.println(divergence);
		}
		for (int i = 0; i < iters; ++i) {
			if (updateH) {
				Matrix WH = W.mmul(H);
				Matrix WHbetaSub2 = WH.pow(beta - 2.0f);
				Matrix WHbetaSub1 = WH.pow(beta - 1.0f);
				H = H.mul(W.transpose().mmul(WHbetaSub2.mul(X))).div(W.transpose().mmul(WHbetaSub1));
				H.maxi(minEps);
			}
			if (updateW) {
				Matrix WH = W.mmul(H);
				Matrix WHbetaSub2 = WH.pow(beta - 2.0f);
				Matrix WHbetaSub1 = WH.pow(beta - 1.0f);
				W = W.mul((X.mul(WHbetaSub2)).mmul(H.transpose())).div((WHbetaSub1).mmul(H.transpose()));
				W.maxi(minEps);
			}
			if (i % PRINT_PERIOD == 0) {
				// Matrix WH = W.mmul(H);
				// float divergence = 0.0f;
				// if (beta == 0.0f) {
				// Matrix XdivWH = X.div(WH);
				// divergence = (XdivWH.sub(XdivWH.log())).norm1() - 1.0f;
				// } else if (beta == 1.0f) {
				// divergence =
				// (X.mul(X.log().sub(WH.log())).add(WH.sub(X))).norm1();
				// } else {
				// divergence = (1.0f / (beta * (beta - 1.0f))) *
				// (X.pow(beta).norm1() + (beta - 1.0f) * WH.pow(beta).norm1() -
				// beta * (X.mul(WH.pow(beta-1.0f)).norm1()));
				// }
				// System.out.println(divergence);
				System.out.print(".");
			}
			CublasUtil.freeAllBut(W, H, X);
		}

		System.out.println();

		float[][] resultWtrans = W.transpose().toArray2();
		float[][] resultHtrans = H.transpose().toArray2();
		CublasUtil.freeAll();

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(resultWtrans, insertSilence(resultHtrans, loud));
		} else {
			return Pair.makePair(resultWtrans, resultHtrans);
		}
	}

	public static Pair<float[][], float[][]> nmfBetaL2PriorExpGradGPU(float[][] priorWtrans, float[][] initWtrans,
			final float priorWeightW, boolean updateW, float[][] priorHtrans, float[][] initHtrans,
			final float priorWeightH, boolean updateH, float[][] Xtrans0, float startStepSize, float endStepSize,
			int iters, float silenceEps, float minEps, float beta) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			priorHtrans = filterSilence(priorHtrans, loud);
			initHtrans = filterSilence(initHtrans, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		final int n = Xtrans0[0].length;
		final int m = Xtrans0.length;
		final int r = initWtrans.length;

		// System.out.println("NMF Beta Exp Grad GPU: ");
		System.out.print("Fitting model.");

		long start = System.nanoTime();

		Matrix X = Matrix.build(a.transpose(Xtrans0));
		Matrix Wprior = Matrix.build(a.transpose(priorWtrans));
		Matrix Hprior = Matrix.build(a.transpose(priorHtrans));
		Matrix W = Matrix.build(a.transpose(initWtrans));
		Matrix H = Matrix.build(a.transpose(initHtrans));

		{
			// Matrix WH = W.mmul(H);
			// float divergence = 0.0f;
			// if (beta == 0.0f) {
			// Matrix XdivWH = X.div(WH);
			// divergence = (XdivWH.sub(XdivWH.log())).norm1() - 1.0f;
			// } else if (beta == 1.0f) {
			// divergence =
			// (X.mul(X.log().sub(WH.log())).add(WH.sub(X))).norm1();
			// } else {
			// divergence = (1.0f / (beta * (beta - 1.0f))) *
			// (X.pow(beta).norm1() + (beta - 1.0f) * WH.pow(beta).norm1() -
			// beta * (X.mul(WH.pow(beta-1.0f)).norm1()));
			// }
			//// System.out.println(divergence / (n*m) + (priorWeightW > 0.0 ?
			// priorWeightW * Math.pow(W.distance2(Wprior), 2.0) : 0.0) / (n*r)
			// + (priorWeightH > 0.0 ? priorWeightH *
			// Math.pow(H.distance2(Hprior), 2.0) : 0.0) / (m*r));
			// System.out.println(divergence / (m) + (priorWeightW > 0.0 ?
			// priorWeightW * Math.pow(W.distance2(Wprior), 2.0) : 0.0) +
			// (priorWeightH > 0.0 ? priorWeightH *
			// Math.pow(H.distance2(Hprior), 2.0) : 0.0) / (m));
		}
		for (int i = 0; i < iters; ++i) {
			float stepSize = startStepSize + ((float) i / iters) * (endStepSize - startStepSize);
			if (updateH) {
				Matrix WH = W.mmul(H);
				Matrix WHbetaSub1 = WH.pow(beta - 1.0f);
				Matrix WHbetaSub2 = WH.pow(beta - 2.0f);
				Matrix priorHGrad = H.sub(Hprior).mul(priorWeightH * 2.0f);
				priorHGrad.muli(1.0f / (m));
				Matrix divergeHGrad = W.transpose().mmul(WHbetaSub1.sub(X.mul(WHbetaSub2)));
				divergeHGrad.muli(1.0f / (m));
				Matrix gradH = divergeHGrad.add(priorHGrad);
				H = H.mul((gradH.mul(-stepSize)).exp());
				H.maxi(minEps);
			}
			if (updateW) {
				Matrix WH = W.mmul(H);
				Matrix WHbetaSub1 = WH.pow(beta - 1.0f);
				Matrix WHbetaSub2 = WH.pow(beta - 2.0f);
				Matrix priorWGrad = W.sub(Wprior).mul(priorWeightW * 2.0f);
				Matrix divergeWGrad = (WHbetaSub1.sub(X.mul(WHbetaSub2))).mmul(H.transpose());
				divergeWGrad.muli(1.0f / (m));
				Matrix gradW = divergeWGrad.add(priorWGrad);
				W = W.mul((gradW.mul(-stepSize)).exp());
				W.maxi(minEps);
			}
			if (i % PRINT_PERIOD == 0) {
				// Matrix WH = W.mmul(H);
				// float divergence = 0.0f;
				// if (beta == 0.0f) {
				// Matrix XdivWH = X.div(WH);
				// divergence = (XdivWH.sub(XdivWH.log())).norm1() - 1.0f;
				// } else if (beta == 1.0f) {
				// divergence =
				// (X.mul(X.log().sub(WH.log())).add(WH.sub(X))).norm1();
				// } else {
				// divergence = (1.0f / (beta * (beta - 1.0f))) *
				// (X.pow(beta).norm1() + (beta - 1.0f) * WH.pow(beta).norm1() -
				// beta * (X.mul(WH.pow(beta-1.0f)).norm1()));
				// }
				//// System.out.println(divergence / (n*m) + (priorWeightW > 0.0
				// ? priorWeightW * Math.pow(W.distance2(Wprior), 2.0) : 0.0) /
				// (n*r) + (priorWeightH > 0.0 ? priorWeightH *
				// Math.pow(H.distance2(Hprior), 2.0) : 0.0) / (m*r));
				// System.out.println(divergence / (m) + (priorWeightW > 0.0 ?
				// priorWeightW * Math.pow(W.distance2(Wprior), 2.0) : 0.0) +
				// (priorWeightH > 0.0 ? priorWeightH *
				// Math.pow(H.distance2(Hprior), 2.0) : 0.0) / (m));
				System.out.print(".");
			}
			CublasUtil.freeAllBut(W, Wprior, H, Hprior, X);
		}

		System.out.println();

		float[][] resultWtrans = W.transpose().toArray2();
		float[][] resultHtrans = H.transpose().toArray2();
		CublasUtil.freeAll();

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(resultWtrans, insertSilence(resultHtrans, loud));
		} else {
			return Pair.makePair(resultWtrans, resultHtrans);
		}
	}

	public static Pair<float[][], float[][]> nmfBetaL2PriorExpGrad(float[][] priorWtrans, float[][] initWtrans,
			final float priorWeightW, boolean updateW, float[][] priorHtrans, float[][] initHtrans,
			final float priorWeightH, boolean updateH, float[][] Xtrans0, float startStepSize, float endStepSize,
			int iters, float silenceEps, float minEps, float beta) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			priorHtrans = filterSilence(priorHtrans, loud);
			initHtrans = filterSilence(initHtrans, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		final int n = Xtrans0[0].length;
		final int m = Xtrans0.length;
		final int r = initWtrans.length;

		// System.out.println("NMF Beta Exp Grad: ");
		System.out.print("Fitting model.");

		long start = System.nanoTime();

		FloatMatrix X = new FloatMatrix(a.transpose(Xtrans0));
		FloatMatrix Wprior = new FloatMatrix(a.transpose(priorWtrans));
		FloatMatrix Hprior = new FloatMatrix(a.transpose(priorHtrans));
		FloatMatrix W = new FloatMatrix(a.transpose(initWtrans));
		FloatMatrix H = new FloatMatrix(a.transpose(initHtrans));

		for (int i = 0; i < iters; ++i) {
			float stepSize = startStepSize + ((float) i / iters) * (endStepSize - startStepSize);
			if (updateH) {
				FloatMatrix WH = W.mmul(H);
				FloatMatrix WHbetaSub1 = MatrixFunctions.pow(WH, beta - 1.0f);
				FloatMatrix WHbetaSub2 = MatrixFunctions.pow(WH, beta - 2.0f);
				FloatMatrix priorHGrad = H.sub(Hprior).mul(priorWeightH * 2.0f);
				priorHGrad.muli(1.0f / (m));
				FloatMatrix divergeHGrad = W.transpose().mmul(WHbetaSub1.sub(X.mul(WHbetaSub2)));
				divergeHGrad.muli(1.0f / (m));
				FloatMatrix gradH = divergeHGrad.add(priorHGrad);
				H = H.mul(MatrixFunctions.exp(gradH.mul(-stepSize)));
				H.maxi(minEps);
			}
			if (updateW) {
				FloatMatrix WH = W.mmul(H);
				FloatMatrix WHbetaSub1 = MatrixFunctions.pow(WH, beta - 1.0f);
				FloatMatrix WHbetaSub2 = MatrixFunctions.pow(WH, beta - 2.0f);
				FloatMatrix priorWGrad = W.sub(Wprior).mul(priorWeightW * 2.0f);
				FloatMatrix divergeWGrad = (WHbetaSub1.sub(X.mul(WHbetaSub2))).mmul(H.transpose());
				divergeWGrad.muli(1.0f / (m));
				FloatMatrix gradW = divergeWGrad.add(priorWGrad);
				W = W.mul(MatrixFunctions.exp(gradW.mul(-stepSize)));
				W.maxi(minEps);
			}
			if (i % PRINT_PERIOD == 0) {
				System.out.print(".");
			}
		}

		System.out.println();

		float[][] resultWtrans = W.transpose().toArray2();
		float[][] resultHtrans = H.transpose().toArray2();

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(resultWtrans, insertSilence(resultHtrans, loud));
		} else {
			return Pair.makePair(resultWtrans, resultHtrans);
		}
	}

	public static Pair<float[][], float[][]> nmfLogNormalExpGradGPU(float[][] Xtrans0, int r, float startStepSize,
			float endStepSize, int iters, float silenceEps, float minEps, float c) {
		int n = Xtrans0[0].length;
		int m = Xtrans0.length;
		Random rand = new Random(RAND_SEED);
		return nmfLogNormalExpGradGPU(a.randFloat(r, n, rand), true, a.randFloat(m, r, rand), true, Xtrans0,
				startStepSize, endStepSize, iters, silenceEps, minEps, c);
	}

	public static Pair<float[][], float[][]> nmfLogNormalExpGradGPU(float[][] Wtrans0, boolean updateW,
			float[][] Htrans0, boolean updateH, float[][] Xtrans0, float startStepSize, float endStepSize, int iters,
			float silenceEps, float minEps, float c) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			Htrans0 = filterSilence(Htrans0, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		final int n = Xtrans0[0].length;
		final int m = Xtrans0.length;

		System.out.println("NMF Log Normal Exp Grad GPU: ");

		long start = System.nanoTime();

		Matrix W = Matrix.build(a.transpose(Wtrans0));
		Matrix H = Matrix.build(a.transpose(Htrans0));
		Matrix X = Matrix.build(a.transpose(Xtrans0));
		Matrix logX = X.max(MIN_LOG_EPS).log();
		Matrix C = Matrix.ones(n, m).mul(c);

		{
			System.out.println(Math.pow(logX.sub((W.mmul(H)).max(MIN_LOG_EPS).log()).add(C).norm2(), 2.0));
		}
		for (int i = 0; i < iters; ++i) {
			float stepSize = startStepSize + ((float) i / iters) * (endStepSize - startStepSize);
			if (updateH) {
				Matrix WH = W.mmul(H);
				Matrix logWH = WH.max(MIN_LOG_EPS).log();
				Matrix logXsubLogWHdivWH = logX.sub(logWH).add(C).div(WH);
				Matrix gradH = W.transpose().mmul(logXsubLogWHdivWH).mul(-2.0f);
				H = H.mul((gradH.mul(-stepSize)).exp());
				H.maxi(minEps);
			}
			if (updateW) {
				Matrix WH = W.mmul(H);
				Matrix logWH = WH.max(MIN_LOG_EPS).log();
				Matrix logXsubLogWHdivWH = logX.sub(logWH).add(C).div(WH);
				Matrix gradW = logXsubLogWHdivWH.mmul(H.transpose()).mul(-2.0f);
				W = W.mul((gradW.mul(-stepSize)).exp());
				W.maxi(minEps);
			}
			if (i % PRINT_PERIOD == 0) {
				System.out.println(Math.pow(logX.sub((W.mmul(H)).max(MIN_LOG_EPS).log()).add(C).norm2(), 2.0));
			}
			CublasUtil.freeAllBut(W, H, X, logX, C);
		}

		float[][] resultWtrans = W.transpose().toArray2();
		float[][] resultHtrans = H.transpose().toArray2();
		CublasUtil.freeAll();

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(resultWtrans, insertSilence(resultHtrans, loud));
		} else {
			return Pair.makePair(resultWtrans, resultHtrans);
		}
	}

	public static Pair<float[][], float[][]> nmfLogNormalExpGrad(float[][] Xtrans0, int r, float startStepSize,
			float endStepSize, int iters, float silenceEps, float minEps, float c) {
		int n = Xtrans0[0].length;
		int m = Xtrans0.length;
		Random rand = new Random(RAND_SEED);
		return nmfLogNormalExpGrad(a.randFloat(r, n, rand), true, a.randFloat(m, r, rand), true, Xtrans0, startStepSize,
				endStepSize, iters, silenceEps, minEps, c);
	}

	public static Pair<float[][], float[][]> nmfLogNormalExpGrad(float[][] Wtrans0, boolean updateW, float[][] Htrans0,
			boolean updateH, float[][] Xtrans0, float startStepSize, float endStepSize, int iters, float silenceEps,
			float minEps, float c) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			Htrans0 = filterSilence(Htrans0, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		final int n = Xtrans0[0].length;
		final int m = Xtrans0.length;

		System.out.println("NMF Log Normal Exp Grad: ");

		long start = System.nanoTime();

		FloatMatrix W = new FloatMatrix(a.transpose(Wtrans0));
		FloatMatrix H = new FloatMatrix(a.transpose(Htrans0));
		FloatMatrix X = new FloatMatrix(a.transpose(Xtrans0));
		FloatMatrix logX = MatrixFunctions.log(X.max(MIN_LOG_EPS));
		FloatMatrix C = FloatMatrix.ones(n, m).mul(c);

		{
			System.out.println(Math.pow(logX.sub(MatrixFunctions.log(W.mmul(H).max(MIN_LOG_EPS))).add(C).norm2(), 2.0));
		}
		for (int i = 0; i < iters; ++i) {
			float stepSize = startStepSize + ((float) i / iters) * (endStepSize - startStepSize);
			if (updateH) {
				FloatMatrix WH = W.mmul(H);
				FloatMatrix logWH = MatrixFunctions.log(WH.max(MIN_LOG_EPS));
				FloatMatrix logXsubLogWHdivWH = logX.sub(logWH).add(C).div(WH);
				FloatMatrix gradH = W.transpose().mmul(logXsubLogWHdivWH).mul(-2.0f);
				H = H.mul(MatrixFunctions.exp(gradH.mul(-stepSize)));
				H.maxi(minEps);
			}
			if (updateW) {
				FloatMatrix WH = W.mmul(H);
				FloatMatrix logWH = MatrixFunctions.log(WH.max(MIN_LOG_EPS));
				FloatMatrix logXsubLogWHdivWH = logX.sub(logWH).add(C).div(WH);
				FloatMatrix gradW = logXsubLogWHdivWH.mmul(H.transpose()).mul(-2.0f);
				W = W.mul(MatrixFunctions.exp(gradW.mul(-stepSize)));
				W.maxi(minEps);
			}
			if (i % PRINT_PERIOD == 0) {
				System.out.println(
						Math.pow(logX.sub(MatrixFunctions.log(W.mmul(H).max(MIN_LOG_EPS))).add(C).norm2(), 2.0));
			}
		}

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(W.transpose().toArray2(), insertSilence(H.transpose().toArray2(), loud));
		} else {
			return Pair.makePair(W.transpose().toArray2(), H.transpose().toArray2());
		}
	}

	public static Pair<float[][], float[][]> nmfL2(float[][] Xtrans0, int r, int iters, float silenceEps,
			float minEps) {
		int n = Xtrans0[0].length;
		int m = Xtrans0.length;
		Random rand = new Random(RAND_SEED);
		return nmfL2(a.randFloat(r, n, rand), true, a.randFloat(m, r, rand), true, Xtrans0, iters, silenceEps, minEps);
	}

	public static Pair<float[][], float[][]> nmfL2(float[][] Wtrans0, boolean updateW, float[][] Htrans0,
			boolean updateH, float[][] Xtrans0, int iters, float silenceEps, float minEps) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			Htrans0 = filterSilence(Htrans0, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		System.out.println("NMF L2: ");

		long start = System.nanoTime();

		FloatMatrix W = new FloatMatrix(a.transpose(Wtrans0));
		FloatMatrix H = new FloatMatrix(a.transpose(Htrans0));
		FloatMatrix X = new FloatMatrix(a.transpose(Xtrans0));

		{
			float norm2 = X.distance2(W.mmul(H));
			System.out.println(norm2 * norm2);
		}
		for (int i = 0; i < iters; ++i) {
			if (updateH) {
				H = H.mul(W.transpose().mmul(X)).div(W.transpose().mmul(W.mmul(H)));
				H.maxi(minEps);
			}
			if (updateW) {
				W = W.mul(X.mmul(H.transpose())).div((W.mmul(H)).mmul(H.transpose()));
				W.maxi(minEps);
			}
			if (i % PRINT_PERIOD == 0) {
				float norm2 = X.distance2(W.mmul(H));
				System.out.println(norm2 * norm2);
			}
		}

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(W.transpose().toArray2(), insertSilence(H.transpose().toArray2(), loud));
		} else {
			return Pair.makePair(W.transpose().toArray2(), H.transpose().toArray2());
		}
	}

	public static Pair<float[][], float[][]> nmfL2GPU(float[][] Xtrans0, int r, int iters, float silenceEps,
			float minEps) {
		int n = Xtrans0[0].length;
		int m = Xtrans0.length;
		Random rand = new Random(RAND_SEED);
		return nmfL2GPU(a.randFloat(r, n, rand), true, a.randFloat(m, r, rand), true, Xtrans0, iters, silenceEps,
				minEps);
	}

	public static Pair<float[][], float[][]> nmfL2GPU(float[][] Wtrans0, boolean updateW, float[][] Htrans0,
			boolean updateH, float[][] Xtrans0, int iters, float silenceEps, float minEps) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			Htrans0 = filterSilence(Htrans0, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		System.out.println("NMF L2 GPU: ");

		long start = System.nanoTime();

		Matrix W = Matrix.build(a.transpose(Wtrans0));
		Matrix H = Matrix.build(a.transpose(Htrans0));
		Matrix X = Matrix.build(a.transpose(Xtrans0));

		{
			float norm2 = X.distance2(W.mmul(H));
			System.out.println(norm2 * norm2);
		}
		for (int i = 0; i < iters; ++i) {
			if (updateH) {
				H = H.mul(W.transpose().mmul(X)).div(W.transpose().mmul(W.mmul(H)));
				H.maxi(minEps);
			}
			if (updateW) {
				W = W.mul(X.mmul(H.transpose())).div((W.mmul(H)).mmul(H.transpose()));
				W.maxi(minEps);
			}
			if (i % PRINT_PERIOD == 0) {
				float norm2 = X.distance2(W.mmul(H));
				System.out.println(norm2 * norm2);
			}
			CublasUtil.freeAllBut(W, H, X);
		}

		float[][] resultWtrans = W.transpose().toArray2();
		float[][] resultHtrans = H.transpose().toArray2();
		CublasUtil.freeAll();

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(resultWtrans, insertSilence(resultHtrans, loud));
		} else {
			return Pair.makePair(resultWtrans, resultHtrans);
		}
	}

	public static Pair<float[][], float[][]> nmfKL(float[][] Xtrans0, int r, int iters, float silenceEps,
			float minEps) {
		int n = Xtrans0[0].length;
		int m = Xtrans0.length;
		Random rand = new Random(RAND_SEED);
		return nmfKL(a.randFloat(r, n, rand), true, a.randFloat(m, r, rand), true, Xtrans0, iters, silenceEps, minEps);
	}

	public static Pair<float[][], float[][]> nmfKL(float[][] Wtrans0, boolean updateW, float[][] Htrans0,
			boolean updateH, float[][] Xtrans0, int iters, float silenceEps, float minEps) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			System.out.println("silence eps: " + silenceEps);
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			Htrans0 = filterSilence(Htrans0, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		System.out.println("NMF KL: ");

		long start = System.nanoTime();

		int n = Xtrans0[0].length;
		int m = Xtrans0.length;

		FloatMatrix W = new FloatMatrix(a.transpose(Wtrans0));
		FloatMatrix H = new FloatMatrix(a.transpose(Htrans0));
		FloatMatrix X = new FloatMatrix(a.transpose(Xtrans0));
		FloatMatrix logX = MatrixFunctions.log(X);

		FloatMatrix denomUpdateW = FloatMatrix.ones(n, 1).mmul(FloatMatrix.ones(1, m).mmul(H.transpose()));
		FloatMatrix denomUpdateH = (W.transpose().mmul(FloatMatrix.ones(n, 1))).mmul(FloatMatrix.ones(1, m));
		{
			FloatMatrix WH = W.mmul(H);
			System.out.println(((X.mul(logX.sub(MatrixFunctions.log(WH))).sub(X)).add(WH)).norm1());
		}
		for (int i = 0; i < iters; ++i) {
			if (updateH) {
				H = H.mul(W.transpose().mmul(X.div(W.mmul(H)))).div(denomUpdateH);
				H.maxi(minEps);
				denomUpdateW = FloatMatrix.ones(n, 1).mmul(FloatMatrix.ones(1, m).mmul(H.transpose()));
			}
			if (updateW) {
				W = W.mul((X.div(W.mmul(H))).mmul(H.transpose())).div(denomUpdateW);
				W.maxi(minEps);
				denomUpdateH = (W.transpose().mmul(FloatMatrix.ones(n, 1))).mmul(FloatMatrix.ones(1, m));
			}
			if (i % PRINT_PERIOD == 0) {
				FloatMatrix WH = W.mmul(H);
				System.out.println(((X.mul(logX.sub(MatrixFunctions.log(WH))).sub(X)).add(WH)).norm1());
			}
		}

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(W.transpose().toArray2(), insertSilence(H.transpose().toArray2(), loud));
		} else {
			return Pair.makePair(W.transpose().toArray2(), H.transpose().toArray2());
		}
	}

	public static Pair<float[][], float[][]> nmfKLGPU(float[][] Xtrans0, int r, int iters, float silenceEps,
			float minEps) {
		int n = Xtrans0[0].length;
		int m = Xtrans0.length;
		Random rand = new Random(RAND_SEED);
		return nmfKLGPU(a.randFloat(r, n, rand), true, a.randFloat(m, r, rand), true, Xtrans0, iters, silenceEps,
				minEps);
	}

	public static Pair<float[][], float[][]> nmfKLGPU(float[][] Wtrans0, boolean updateW, float[][] Htrans0,
			boolean updateH, float[][] Xtrans0, int iters, float silenceEps, float minEps) {
		int[] loud = null;
		if (silenceEps > 0.0) {
			Pair<float[][], int[]> filterSilence = filterSilence(Xtrans0, silenceEps);
			Xtrans0 = filterSilence.getFirst();
			loud = filterSilence.getSecond();
			Htrans0 = filterSilence(Htrans0, loud);
		}

		// W : n x r
		// H : r x m
		// X : n x m

		System.out.println("NMF KL GPU: ");

		long start = System.nanoTime();

		int n = Xtrans0[0].length;
		int m = Xtrans0.length;

		Matrix W = Matrix.build(a.transpose(Wtrans0));
		Matrix H = Matrix.build(a.transpose(Htrans0));
		Matrix X = Matrix.build(a.transpose(Xtrans0));
		Matrix logX = X.log();

		Matrix denomUpdateW = Matrix.ones(n, 1).mmul(Matrix.ones(1, m).mmul(H.transpose()));
		Matrix denomUpdateH = (W.transpose().mmul(Matrix.ones(n, 1))).mmul(Matrix.ones(1, m));
		{
			Matrix WH = W.mmul(H);
			System.out.println(((X.mul(logX.sub(WH.log())).sub(X)).add(WH)).norm1());
		}
		for (int i = 0; i < iters; ++i) {
			if (updateH) {
				H = H.mul(W.transpose().mmul(X.div(W.mmul(H)))).div(denomUpdateH);
				H.maxi(minEps);
				denomUpdateW = Matrix.ones(n, 1).mmul(Matrix.ones(1, m).mmul(H.transpose()));
			}
			if (updateW) {
				W = W.mul((X.div(W.mmul(H))).mmul(H.transpose())).div(denomUpdateW);
				W.maxi(minEps);
				denomUpdateH = (W.transpose().mmul(Matrix.ones(n, 1))).mmul(Matrix.ones(1, m));
			}
			if (i % PRINT_PERIOD == 0) {
				Matrix WH = W.mmul(H);
				System.out.println(((X.mul(logX.sub(WH.log())).sub(X)).add(WH)).norm1());
			}
			CublasUtil.freeAllBut(W, H, X, logX, denomUpdateW, denomUpdateH);
		}

		float[][] resultWtrans = W.transpose().toArray2();
		float[][] resultHtrans = H.transpose().toArray2();
		CublasUtil.freeAll();

		long end = System.nanoTime();
		System.out.println("Compute time: " + (end - start) / (1e9) + "s");

		if (silenceEps > 0.0) {
			return Pair.makePair(resultWtrans, insertSilence(resultHtrans, loud));
		} else {
			return Pair.makePair(resultWtrans, resultHtrans);
		}
	}

	private static Pair<float[][], int[]> filterSilence(float[][] X, float thresh) {
		List<float[]> Xs = new ArrayList<float[]>();
		int[] loud = new int[X.length];
		int c = 0;
		for (int i = 0; i < X.length; ++i) {
			float[] f = X[i];
			if (a.sum(f) > thresh) {
				Xs.add(f);
				loud[i] = c;
				c++;
			} else {
				loud[i] = -1;
			}
		}
		return new Pair<float[][], int[]>(Xs.toArray(new float[0][0]), loud);
	}

	private static float[][] filterSilence(float[][] X, int[] loud) {
		List<float[]> Xs = new ArrayList<float[]>();
		for (int i = 0; i < X.length; ++i) {
			if (loud[i] != -1) {
				float[] f = X[i];
				Xs.add(f);
			}
		}
		return Xs.toArray(new float[0][0]);
	}

	private static float[][] insertSilence(float[][] Xs, int[] loud) {
		float[][] X = new float[loud.length][];
		for (int i = 0; i < loud.length; ++i) {
			if (loud[i] >= 0) {
				X[i] = Xs[loud[i]];
			} else {
				X[i] = new float[Xs[0].length];
			}
		}
		return X;
	}

	public static void main(String[] args) {
//		CublasUtil.startup(1);

		int iters = 1000;
		int n = 100;
		int m = 100;
		int r = 50;
//		int n = 10;
//		int m = 10;
//		int r = 5;

		float[][] X = a.randFloat(n, m, new Random());
		// a.addi(X, -0.5f);
		// a.scalei(X, 10.0f);
		// a.expi(X);
		a.scalei(X, 1.0f / a.max(a.max(X)));

		float[][] Wtrans = a.scale(a.randFloat(r, m, new Random()), 3e-2f);
		float[][] Htrans = a.scale(a.randFloat(n, r, new Random()), 3e-2f);
		float[][] WtransPrior = a.scale(a.randFloat(r, m, new Random()), 3e-2f);
		float[][] HtransPrior = a.scale(a.randFloat(n, r, new Random()), 3e-2f);

		// gradientTester(Wtrans, Htrans, X, 1e-6f);

		// nmfLogNormalExpGrad(Wtrans, true, Htrans, true, X, 1e-5f, 1e-5f,
		// iters, 1e-8f, 1e-8f, 0.5f);
		// nmfLogNormalExpGradGPU(Wtrans, true, Htrans, true, X, 1e-5f, 1e-5f,
		// iters, 1e-8f, 1e-8f, 0.5f);

		
		 nmfBeta(Wtrans, true, Htrans, true, X, iters, 1e-8f, 0.0f, 1.0f);
		// result = nmfBetaGPU(Wtrans, true, Htrans, true, X, iters, 1e-8f,
		// 0.0f, 1.0f);
//		nmfBetaL2PriorExpGradGPU(WtransPrior, Wtrans, 1.0f, true, HtransPrior, Htrans, 1.0f, true, X, 1e-3f, 1e-3f,
//				iters, 1e-10f, 1e-20f, 0.5f);

		// nmfKLLBFGSGPU(Wtrans, 1e0f, Htrans, 1e0f, X, iters, 1e-8f);
		// nmfKLLBFGS(Wtrans, 1e0f, Htrans, 1e0f, X, iters, 1e-8f);

		// nmfL2LBFGS(Wtrans, 1e0f, Htrans, 1e0f, X, iters, 1e-8f);
		// nmfL2LBFGSGPU(Wtrans, 1e0f, Htrans, 1e0f, X, iters, 1e-8f);
		// nmfL2(Wtrans, true, Htrans, true, X, iters, 1e-8f, 0.0f);
		// nmfL2GPU(Wtrans, true, Htrans, true, X, iters, 1e-8f, 0.0f);

		// nmfKLLBFGSGPU(Wtrans, 0.0f, Htrans, 0.0f, X, iters, 1e-8f);
		// nmfKLGPU(Wtrans, true, Htrans, true, X, iters, 1e-8f, 0.0f);
		// nmfKLLBFGS(Wtrans, 0.0f, Htrans, 0.0f, X, iters, 1e-8f);
		// nmfKL(Wtrans, true, Htrans, true, X, iters, 1e-8f, 0.0f);

//		CublasUtil.shutdown();
	}

}
