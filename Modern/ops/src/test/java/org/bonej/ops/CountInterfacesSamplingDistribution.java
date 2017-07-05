
package org.bonej.ops;

import static org.bonej.ops.CountInterfaces.Plane;
import static org.bonej.ops.CountInterfaces.createSamplePoints;
import static org.bonej.ops.CountInterfaces.createSegment;
import static org.bonej.ops.CountInterfaces.createStackPlanes;
import static org.bonej.ops.CountInterfaces.findBounds;
import static org.bonej.ops.CountInterfaces.findDirection;
import static org.bonej.ops.CountInterfaces.selectPlanes;
import static org.bonej.ops.CountInterfaces.toVoxelCoordinates;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import net.imagej.ImageJ;
import net.imglib2.Cursor;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.ComplexType;
import net.imglib2.type.numeric.integer.LongType;
import net.imglib2.type.numeric.integer.UnsignedIntType;
import net.imglib2.type.numeric.integer.UnsignedShortType;
import net.imglib2.util.ValuePair;
import net.imglib2.view.Views;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
public class CountInterfacesSamplingDistribution {

	private synchronized static <C extends ComplexType<C>> void sample(
		final RandomAccessibleInterval<C> interval, final long[] bounds, final Vector3d point)
	{
		final long[] coordinates = toVoxelCoordinates(point, bounds);
		final RandomAccess<C> access = interval.randomAccess(interval);
		access.setPosition(coordinates);
		final C voxel = access.get();
		//System.out.println(Arrays.toString(coordinates));
		voxel.setReal(voxel.getRealDouble() + 1.0);
	}

	private static <C extends ComplexType<C>> void plotOrigin(final RandomAccessibleInterval interval) {
		final List<Plane> planes = createStackPlanes(interval);

	}

	public static void main(String... args) throws IOException {
/*
        System.out.println(Arrays.toString(CountInterfaces.getOrientation(new Vector3d(-1, 0, 0))));
        System.out.println(Arrays.toString(CountInterfaces.getOrientation(new Vector3d(1, 0, 0))));
        System.out.println(Arrays.toString(CountInterfaces.getOrientation(new Vector3d(0, -1, 0))));
        System.out.println(Arrays.toString(CountInterfaces.getOrientation(new Vector3d(0, 1, 0))));
        System.out.println(Arrays.toString(CountInterfaces.getOrientation(new Vector3d(0, 0, -1))));
        System.out.println(Arrays.toString(CountInterfaces.getOrientation(new Vector3d(0, 0, 1))));
        */
        final long width = 10;
		final long height = 10;
		final long depth = 10;
		final RandomAccessibleInterval<UnsignedIntType> stack = ArrayImgs.unsignedInts(width, height, depth);
		CountInterfaces.setSeed(0xC0FF33);
		final List<Plane> planes = createStackPlanes(stack);
		final long[] bounds = CountInterfaces.findBounds(stack);

		final long startTime = System.nanoTime();
		final Stream<ValuePair<Plane, Plane>> pairs = Stream.generate(() -> selectPlanes(planes));
		final Stream<Vector3d> samplePoints = CountInterfaces.createSamplePoints(pairs, 1 / 2.3, bounds);

		//TODO random perioidicity?
		samplePoints.limit(1_000_000).sequential().forEach(point -> sample(stack, bounds, point));
		final long endTime = System.nanoTime();
		final long durationMs = (endTime - startTime) / 1_000_000;
		System.out.println("Duration (ms) " + durationMs);

		final IterableInterval<UnsignedIntType> iterable = Views.flatIterable(stack);
		final SummaryStatistics statistics = new SummaryStatistics();
		iterable.cursor().forEachRemaining(i -> statistics.addValue(i.getIntegerLong()));
		System.out.println(statistics.toString());
		final ImageJ imageJ = new ImageJ();
        imageJ.launch(args);
        imageJ.ui().show(stack);
    }
}
