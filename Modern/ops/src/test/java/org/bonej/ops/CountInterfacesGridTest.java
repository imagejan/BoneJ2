
package org.bonej.ops;

import static org.bonej.ops.CountInterfacesGrid.findGridSize;
import static org.bonej.ops.CountInterfacesGrid.plotSamplers;

import java.util.Random;
import java.util.concurrent.locks.ReadWriteLock;
import java.util.concurrent.locks.ReentrantReadWriteLock;
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
	public static void plotSamplingOrientation(final double[][][][] orientations,
		final long[] bounds, final double gridSize, final long rotations)
	{
		final Stream<ValuePair<Vector3d, Vector3d>> samplers = Stream.generate(
			() -> plotSamplers(gridSize, (long) gridSize, bounds)).limit(rotations)
			.flatMap(s -> s);
		final double increment = 1.0 / 2.3;
		ReadWriteLock lock = new ReentrantReadWriteLock();
		samplers.forEach(sampler -> sampleVector(sampler, gridSize, increment, orientations, bounds, lock));
	}

	private static void sampleVector(final ValuePair<Vector3d, Vector3d> sampler, final double gridSize, final double increment, final double[][][][] orientations, final long[] bounds, final ReadWriteLock lock) {
		final long iterations = (long) Math.floor(gridSize / increment) + 1;
		final Vector3d jitter = new Vector3d(sampler.b);
		jitter.scale(Math.random());
		for (int i = 0; i < iterations; i++) {
			final double t = i * increment;
			final Vector3d point = new Vector3d(sampler.b);
			point.scale(t);
			point.add(sampler.a);
			//point.add(jitter);
			final int[] coordinates = CountInterfacesGrid.toVoxelCoordinates(
					point);
			if (CountInterfacesGrid.outOfBounds(coordinates, bounds)) {
				continue;
			}


			//lock.writeLock().lock();
			try {
				orientations[coordinates[0]][coordinates[1]][coordinates[2]][0] +=
						sampler.b.x;
				orientations[coordinates[0]][coordinates[1]][coordinates[2]][1] +=
						sampler.b.y;
				orientations[coordinates[0]][coordinates[1]][coordinates[2]][2] +=
						sampler.b.z;
			} finally {
				//lock.writeLock().unlock();
			}


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
		final long startTime = System.nanoTime();
		plotSamplingOrientation(orientations, bounds, gridSize, rotations);

		final long endTime = System.nanoTime();
		final long durationMs = (endTime - startTime) / 1_000_000;
		final ImgPlus<FloatType> image = orientationsToImage(orientations, bounds);
		// final int[] permutation = IntStream.range(0, depth).toArray();
		// knuthShuffle(permutation);
		// System.out.println(Arrays.toString(permutation));
		// final RandomAccessibleInterval<IntType> view =
		// Views.permuteCoordinates(image, permutation, 3);
		System.out.println("Duration (ms) " + durationMs);
		final ImageJ imageJ = new ImageJ();
		// final Img<IntType> img =
		// imageJ.op().convert().int32(Views.flatIterable(view));
		// final ImgPlus<IntType> imgPlus = new ImgPlus<>(img, "",
		// new DefaultLinearAxis(Axes.X), new DefaultLinearAxis(Axes.Y),
		// new DefaultLinearAxis(Axes.CHANNEL), new DefaultLinearAxis(Axes.Z));
		printStatistics(image);
		imageJ.launch(args);
		imageJ.ui().show(image);
	}
}
