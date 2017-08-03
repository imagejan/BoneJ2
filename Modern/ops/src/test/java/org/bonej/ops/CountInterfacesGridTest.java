
package org.bonej.ops;

import static org.bonej.ops.CountInterfacesGrid.findGridSize;
import static org.bonej.ops.CountInterfacesGrid.findIntersections;
import static org.bonej.ops.CountInterfacesGrid.plotGrid;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Stream;

import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.ComplexType;
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

	public static <C extends ComplexType<C>> void futureParallel(
			final RandomAccessibleInterval<C> interval, final long[] bounds,
			final double gridSize, final long sections, final double increment,
			final long rotations) throws ExecutionException, InterruptedException
	{
		final List<Future> futures = new ArrayList<>();
		final int processors = 8;
		final ExecutorService executors = Executors.newFixedThreadPool(processors);
		final long limit = rotations / processors;
		for (int i = 0; i < processors; i++) {
			final Runnable runnable = () -> Stream.generate(() -> plotGrid(gridSize,
					sections, bounds)).limit(limit).flatMap(s -> s).forEach(
					sampler -> sampleVector(interval, sampler, increment, bounds));
			final Future future = executors.submit(runnable);
			futures.add(future);
		}

		for (final Future future : futures) {
			future.get();
		}

		executors.shutdown();
	}

	/* Plot the average sampling direction */
	public static <C extends ComplexType<C>> void plotSamplingOrientation(
		RandomAccessibleInterval<C> interval, final long[] bounds,
		final double gridSize, final long rotations)
	{

		final Stream<ValuePair<Vector3d, Vector3d>> samplers = Stream.generate(
			() -> plotGrid(gridSize, 100, bounds)).limit(rotations).flatMap(
				s -> s);
		final double increment = 1;
		samplers.forEach(sampler -> sampleVector(interval, sampler, increment,
			bounds));
	}

	private static <C extends ComplexType<C>> void sampleVector(
		final RandomAccessibleInterval<C> interval,
		final ValuePair<Vector3d, Vector3d> sampler, final double increment,
		final long[] bounds)
	{
		final Vector3d origin = sampler.a;
		final Vector3d direction = sampler.b;
		final Vector3d jitter = new Vector3d(direction);
		final double jitterT = random.nextDouble() * increment;
		jitter.scale(jitterT);
		final ValuePair<Double, Double> tPair = findIntersections(sampler, bounds);
		if (tPair == null) {
			return;
		}
		final double tMin = FastMath.min(tPair.a, tPair.b);
		final double tMax = FastMath.max(tPair.a, tPair.b) - jitterT;
		final Vector3d point = new Vector3d();
		point.scaleAdd(tMin, direction, origin);
		point.add(jitter);
		final Vector3d bit = new Vector3d(direction);
		bit.scale(increment);
		final RandomAccess<C> access = interval.randomAccess();
		for (double t = tMin; t <= tMax; t += increment) {
			addOrientation(access, point, direction);
			point.add(bit);
		}
	}

	private static <C extends ComplexType<C>> void addOrientation(
		final RandomAccess<C> access, final Vector3d position,
		final Vector3d orientation)
	{
		access.setPosition((int) position.x, 0);
		access.setPosition((int) position.y, 1);
		access.setPosition((int) position.z, 3);
		access.setPosition(0, 2);
		access.get().setReal(access.get().getRealDouble() + orientation.x);
		access.setPosition(1, 2);
		access.get().setReal(access.get().getRealDouble() + orientation.y);
		access.setPosition(2, 2);
		access.get().setReal(access.get().getRealDouble() + orientation.z);
	}

	private static void printStatistics(
		final RandomAccessibleInterval<FloatType> stack)
	{
		final IterableInterval<FloatType> iterable = Views.flatIterable(stack);
		final SummaryStatistics statistics = new SummaryStatistics();
		iterable.cursor().forEachRemaining(i -> statistics.addValue(i.get()));
		System.out.println(statistics.toString());
	}

	public static void main(String... args) throws ExecutionException, InterruptedException {
		final int width = 100;
		final int height = 100;
		final int depth = 100;
		final int rotations = 1000;
		final long[] bounds = { width, height, depth };
		final double gridSize = findGridSize(bounds);
		//final Img<FloatType> img = ArrayImgs.floats(width, height, 3, depth);
		// plotSamplingOrientation(img, bounds, gridSize, rotations);


		final RandomAccessibleInterval<BitType> img = ArrayImgs.bits(width, height, depth);
		final int iterations = 25;
		long totalTime = 0;
		for (int i = 0; i < iterations; i++) {
			long startTime = System.nanoTime();
			final long sum = CountInterfacesGrid.futureParallel(img, bounds, gridSize, (long) gridSize, 1.0 / 2.3, rotations);
			long endTime = System.nanoTime();
			final long duration = endTime - startTime;
			totalTime += duration;
			long durationMs = duration / 1_000_000;
			System.out.println("Duration A (ms) " + durationMs);
			System.out.println("Intersections " + sum);
		}
		long totalMs = totalTime / 1_000_000;
		System.out.println("Total time " + totalMs + ", Avg. " + totalMs / iterations);

		/*final ImgPlus<FloatType> image = new ImgPlus<>(img, "", new DefaultLinearAxis(Axes.X), new DefaultLinearAxis(Axes.Y), new DefaultLinearAxis(Axes.CHANNEL), new DefaultLinearAxis(Axes.Z));
		printStatistics(image);
		final ImageJ imageJ = new ImageJ();
		imageJ.launch(args);
		imageJ.ui().show(image);*/
	}
}
