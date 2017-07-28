
package org.bonej.ops;

import static org.bonej.ops.CountInterfacesGrid.findGridSize;
import static org.bonej.ops.CountInterfacesGrid.findIntersections;
import static org.bonej.ops.CountInterfacesGrid.outOfBounds;
import static org.bonej.ops.CountInterfacesGrid.plotSamplers;

import java.util.Random;
import java.util.stream.Stream;

import net.imagej.ImageJ;
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
import org.scijava.vecmath.Vector3d;

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
			() -> plotSamplers(gridSize, 100, bounds)).limit(rotations)
			.flatMap(s -> s);
		final double increment = 1;
		samplers.forEach(sampler -> sampleVector(sampler, gridSize, increment,
			orientations, bounds));
		return 0; // samplers.count();
	}



	private static void sampleVector(final ValuePair<Vector3d, Vector3d> sampler,
		final double gridSize, final double increment,
		final double[][][][] orientations, final long[] bounds)
	{
		final Vector3d origin = sampler.a;
		final Vector3d direction = sampler.b;
		final Vector3d jitter = new Vector3d(direction);
		jitter.scale(Math.random());
		final ValuePair<Double, Double> tPair = findIntersections(sampler, bounds);
		if (tPair == null) {
			return;
		}
		// TODO find safe tMinMax and remove outOfBounds
		double tMin = FastMath.min(tPair.a, tPair.b);
		double tMax = FastMath.max(tPair.a, tPair.b);
		if (Double.isNaN(tMin) || Double.isNaN(tMax)) {
			return;
		}
		for (double t = tMin; t <= tMax; t += increment) {
			final Vector3d point = new Vector3d(direction);
			point.scale(t);
			point.add(origin);
			point.add(jitter);
			final int x = (int) FastMath.floor(point.x);
			final int y = (int) FastMath.floor(point.y);
			final int z = (int) FastMath.floor(point.z);
			if (outOfBounds(x, y, z, bounds)) {
				continue;
			}
			orientations[x][y][z][0] += direction.x;
			orientations[x][y][z][1] += direction.y;
			orientations[x][y][z][2] += direction.z;
		}
	}

	private static double safeMin(final Vector3d origin, final Vector3d direction, final Vector3d jitter, final double tMin, final double tMax, final double increment, final long[] bounds) {
		for (double t = tMin; t <= tMax; t += increment) {
			final Vector3d point = new Vector3d(direction);
			point.scale(t);
			point.add(origin);
			point.add(jitter);
			final int x = (int) FastMath.floor(point.x);
			final int y = (int) FastMath.floor(point.y);
			final int z = (int) FastMath.floor(point.z);
			if (x >= 0 && x < bounds[0] && y >= 0 && y < bounds[1] && z >= 0 && z < bounds[2]) {
				return t;
			}
		}
		return Double.NaN;
	}

	private static double safeMax(final Vector3d origin, final Vector3d direction, final Vector3d jitter, final double tMin, final double tMax, final double increment, final long[] bounds) {
		for (double t = tMax; t >= tMin; t -= increment) {
			final Vector3d point = new Vector3d(direction);
			point.scale(t);
			point.add(origin);
			point.add(jitter);
			final int x = (int) FastMath.floor(point.x);
			final int y = (int) FastMath.floor(point.y);
			final int z = (int) FastMath.floor(point.z);
			if (x < bounds[0] && x >= 0 && y < bounds[1] && y >= 0 && z < bounds[2] && z >= 0) {
				return t;
			}
		}
		return Double.NaN;
	}

	private static long vectorIterations(final ValuePair<Double, Double> tPair, final double increment) {
		return (long) FastMath.floor(FastMath.abs(tPair.b - tPair.a) / increment) + 1;
	}

	public static ImgPlus<FloatType> orientationsToImage(
		final double[][][][] orientations, final long[] bounds)
	{
		final ArrayImg<FloatType, FloatArray> ints = ArrayImgs.floats(bounds[0],
			bounds[1], 3, bounds[2]);
		final ImgPlus<FloatType> imgPlus = new ImgPlus<>(ints, "",
			new DefaultLinearAxis(Axes.X), new DefaultLinearAxis(Axes.Y),
			new DefaultLinearAxis(Axes.CHANNEL), new DefaultLinearAxis(Axes.Z));
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

	public static void main(String... args) {
		// final RandomAccessibleInterval<IntType> stack = ArrayImgs.ints(100, 100,
		// 3, 100);
		final int width = 100;
		final int height = 100;
		final int depth = 100;
		final long[] bounds = { width, height, depth };
		final double gridSize = findGridSize(bounds);
		final long rotations = 1000;
		final double[][][][] orientations = new double[width][height][depth][3];

		long startTime = System.nanoTime();
		long count = plotSamplingOrientation(orientations, bounds, gridSize,
			rotations);
		long endTime = System.nanoTime();
		long durationMs = (endTime - startTime) / 1_000_000;
		System.out.println("Duration A (ms) " + durationMs);
		System.out.println("Samplers " + count);

		final ImgPlus<FloatType> image = orientationsToImage(orientations, bounds);
		final ImageJ imageJ = new ImageJ();
		printStatistics(image);
		imageJ.launch(args);
		imageJ.ui().show(image);

	}

}
