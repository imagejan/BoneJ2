
package org.bonej.ops;

import static org.bonej.ops.CountInterfacesGrid.findGridSize;
import static org.bonej.ops.CountInterfacesGrid.foo;
import static org.bonej.ops.CountInterfacesGrid.plotSamplers;
import static org.bonej.ops.CountInterfacesGrid.plotSamplersF;

import java.util.Random;
import java.util.stream.Stream;

import net.imagej.ImgPlus;
import net.imagej.axis.Axes;
import net.imagej.axis.DefaultLinearAxis;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.FloatArray;
import net.imglib2.type.numeric.real.FloatType;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.Pair;
import org.scijava.vecmath.Vector3d;
import org.scijava.vecmath.Vector3f;

/**
 * @author Richard Domander
 */
public class CountInterfacesGridTest {
	private static Random random = new Random(System.currentTimeMillis());
	/*
	public static void plotSampling(final RandomAccessibleInterval<IntType> stack,
		final long[] bounds, final double gridSize, final long rotations)
	{
		final Stream<ValuePair<Vector3d, Vector3d>> samplers = Stream.generate(
			() -> plotSamplers(gridSize, 100, bounds)).limit(rotations).flatMap(
				s -> s);
		final Stream<Vector3d> points = samplePoints(samplers, 1 / 2.3, gridSize);
		points.map(CountInterfacesGrid::toVoxelCoordinates).filter(
			c -> !outOfBounds(c, bounds)).forEach(c -> sample(stack, c));
	}*/

	/* Plot the average sampling direction */
	public static long plotSamplingOrientation(final double[][][][] orientations,
											   final long[] bounds, final double gridSize, final long rotations)
	{

		final Stream<ValuePair<Vector3d, Vector3d>> samplers = Stream.generate(
				() -> plotSamplers(gridSize, (long) gridSize, bounds)).limit(rotations)
				.flatMap(s -> s);
		final double increment = 1 / 2.3;
		samplers.forEach(sampler -> sampleVector(sampler, gridSize, increment, orientations, bounds));
		return 0; //samplers.count();

	}

	private static long plotSamplingOrientationF(final double[][][][] orientations, final long[] bounds, final double gridSize, final long rotations) {
		final Stream<Pair<Vector3f, Vector3f>> samplers = Stream.generate(() -> plotSamplersF((float) gridSize, (long) gridSize, bounds)).limit(rotations).flatMap(s -> s);
		final float increment = 1.0f / 2.3f;
		samplers.forEach(sampler -> sampleVectorF(sampler, gridSize, increment, orientations, bounds));
		return 0; //samplers.count();

	}

	private static void plotSamplingOrientation3(
		final double[][][][] orientations, final long[] bounds,
		final double gridSize, final long rotations)
	{
		final double increment = 1;
		final Stream<double[][]> samplers = Stream.generate(() -> foo(gridSize,
				(int) gridSize, bounds)).limit(rotations).flatMap(s -> s);
		samplers.forEach(sampler -> sample(sampler, gridSize, increment, bounds,
			orientations));
	}

	private static long plotSamplingOrientation2(final double[][][][] orientations, final long[] bounds, final double gridSize, final long rotations) {
		final double increment = 1;
		long count = 0;
		long foo = (long) gridSize;
		for (int i = 0; i < rotations; i++) {
			final double[][][] samplers = plotSamplers(gridSize, (double) foo, bounds);
			count += samplers.length;
			for (final double[][] sampler : samplers) {
				sample(sampler, gridSize, increment, bounds, orientations);
			}

		}
		return count;
	}

	private static void sample(final double[][] sampler, final double gridSize, final double increment, final long[] bounds, final double[][][][] orientations) {
		final double u = Math.random();
		final double[] jitter = {u * sampler[1][0], u * sampler[1][1], u * sampler[1][2]};
		final long iterations = (long) FastMath.floor(gridSize / increment) + 1;
		for (int i = 0; i < iterations; i++) {
			final double t = i * increment;
			final int x = (int) FastMath.floor(sampler[0][0] + sampler[1][0] * t + jitter[0]);
			if (CountInterfacesGrid.outOfBounds(x, bounds[0])) {
				continue;
			}
			final int y = (int) FastMath.floor(sampler[0][1] + sampler[1][1] * t + jitter[1]);
			if (CountInterfacesGrid.outOfBounds(y, bounds[1])) {
				continue;
			}
			final int z = (int) FastMath.floor(sampler[0][2] + sampler[1][2] * t + jitter[2]);
			if (CountInterfacesGrid.outOfBounds(z, bounds[2])) {
				continue;
			}
			orientations[x][y][z][0]++; //+= sampler[1][0];
			orientations[x][y][z][1]++; //+= sampler[1][1];
			orientations[x][y][z][2]++; //+= sampler[1][2];
		}

	}

	private static void sampleVector(final ValuePair<Vector3d, Vector3d> sampler, final double gridSize, final double increment, final double[][][][] orientations, final long[] bounds) {
		final long iterations = (long) Math.floor(gridSize / increment) + 1;
		final Vector3d jitter = new Vector3d(sampler.b);
		jitter.scale(Math.random());
		for (int i = 0; i < iterations; i++) {
			final double t = i * increment;
			final Vector3d point = new Vector3d(sampler.b);
			point.scale(t);
			point.add(sampler.a);
			point.add(jitter);
			final int[] coordinates = CountInterfacesGrid.toVoxelCoordinates(
					point);
			if (CountInterfacesGrid.outOfBounds(coordinates, bounds)) {
				continue;
			}
				orientations[coordinates[0]][coordinates[1]][coordinates[2]][0]++;
				//+= sampler.b.x;
				orientations[coordinates[0]][coordinates[1]][coordinates[2]][1]++;
				//+= sampler.b.y;
				orientations[coordinates[0]][coordinates[1]][coordinates[2]][2]++;
				//+= sampler.b.z;
		}
	}

	private static void sampleVectorF(final Pair<Vector3f, Vector3f> sampler, final double gridSize, final float increment, final double[][][][] orientations, final long[] bounds) {
		final long iterations = (long) Math.floor(gridSize / increment) + 1;
		final Vector3f jitter = new Vector3f(sampler.getSecond());
		jitter.scale((float)Math.random());
		for (int i = 0; i < iterations; i++) {
			final float t = i * increment;
			final Vector3f point = new Vector3f(sampler.getSecond());
			point.scale(t);
			point.add(sampler.getFirst());
			point.add(jitter);
			final int[] coordinates = CountInterfacesGrid.toVoxelCoordinates(point);
			if (CountInterfacesGrid.outOfBounds(coordinates, bounds)) {
				continue;
			}
			orientations[coordinates[0]][coordinates[1]][coordinates[2]][0]++;
			//+= sampler.b.x;
			orientations[coordinates[0]][coordinates[1]][coordinates[2]][1]++;
			//+= sampler.b.y;
			orientations[coordinates[0]][coordinates[1]][coordinates[2]][2]++;
			//+= sampler.b.z;
		}
	}

	//TODO Bresenham 3D coordinates

	public static ImgPlus<FloatType> orientationsToImage(
		final double[][][][] orientations, final long[] bounds)
	{
		final ArrayImg<FloatType, FloatArray> ints = ArrayImgs.floats(bounds[0],
			bounds[1], 3, bounds[2]);
		final ImgPlus<FloatType> imgPlus = new ImgPlus<>(ints, "",
			new DefaultLinearAxis(Axes.X), new DefaultLinearAxis(Axes.Y), new DefaultLinearAxis(Axes.CHANNEL),
			new DefaultLinearAxis(Axes.Z));
		final RandomAccess<FloatType> access = imgPlus.randomAccess();
		for (int k = 0; k < bounds[2]; k++) {
			access.setPosition(k, 3);
			for (int j = 0; j < bounds[1]; j++) {
				access.setPosition(j, 1);
				for (int i = 0; i < bounds[0]; i++) {
					access.setPosition(i, 0);
					access.setPosition(0, 2);
					access.get().set((float) orientations[i][j][k][0]);
					access.setPosition(1, 2);
					access.get().set((float) orientations[i][j][k][1]);
					access.setPosition(2, 2);
					access.get().set((float) orientations[i][j][k][2]);
				}
			}
		}
		return imgPlus;
	}

	public static void printStatistics(
		final RandomAccessibleInterval<FloatType> stack)
	{
		final IterableInterval<FloatType> iterable = Views.flatIterable(stack);
		final SummaryStatistics statistics = new SummaryStatistics();
		iterable.cursor().forEachRemaining(i -> statistics.addValue(i.get()));
		System.out.println(statistics.toString());
	}



	// TODO make sure that same number of samplers in methods
	public static void main(String... args) {
		// final RandomAccessibleInterval<IntType> stack = ArrayImgs.ints(100, 100,
		// 3, 100);
		final int width = 100;
		final int height = 100;
		final int depth = 100;
		final long[] bounds = { width, height, depth };
		final double gridSize = findGridSize(bounds);
		final long rotations = 100;
		final double[][][][] orientations = new double[width][height][depth][3];

		long startTime = System.nanoTime();
		long count = plotSamplingOrientation(orientations, bounds, gridSize, rotations);
		long endTime = System.nanoTime();
		long durationMs = (endTime - startTime) / 1_000_000;
		System.out.println("Duration A (ms) " + durationMs);
		System.out.println("Samplers " + count);

		startTime = System.nanoTime();
		count = plotSamplingOrientationF(orientations, bounds, gridSize, rotations);
		endTime = System.nanoTime();
		durationMs = (endTime - startTime) / 1_000_000;
		System.out.println("Duration D (ms) " + durationMs);
		System.out.println("Samplers " + count);


		/*startTime = System.nanoTime();
		plotSamplingOrientation2(orientations, bounds, gridSize, rotations);
		endTime = System.nanoTime();
		durationMs = (endTime - startTime) / 1_000_000;
		System.out.println("Duration B (ms) " + durationMs);

		startTime = System.nanoTime();
		plotSamplingOrientation3(orientations, bounds, gridSize, rotations);
		endTime = System.nanoTime();
		durationMs = (endTime - startTime) / 1_000_000;
		System.out.println("Duration C (ms) " + durationMs);*/


		final ImgPlus<FloatType> image = orientationsToImage(orientations, bounds);
		//final ImageJ imageJ = new ImageJ();
		printStatistics(image);
		//imageJ.launch(args);
		//imageJ.ui().show(image);
	}

}
