
package org.bonej.ops;

import static org.bonej.ops.CountInterfaces.Plane;
import static org.bonej.ops.CountInterfaces.createSamplePoints;
import static org.bonej.ops.CountInterfaces.createStackPlanes;
import static org.bonej.ops.CountInterfaces.findBounds;
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
		final RandomAccessibleInterval<C> interval, final Vector3d point)
	{
		final long[] coordinates = toVoxelCoordinates(point);
		final RandomAccess<C> access = interval.randomAccess(interval);
		access.setPosition(coordinates);
		final C voxel = access.get();
		voxel.setReal(voxel.getRealDouble() + 1.0);
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
		final long voxels = width * height * depth;
		final RandomAccessibleInterval<UnsignedShortType> stack = ArrayImgs.unsignedShorts(width, height, depth);
		CountInterfaces.setSeed(0xC0FF33);
		final List<Plane> planes = createStackPlanes(stack);
        final Stream<ValuePair<Plane, Plane>> pairs = Stream.generate(
			() -> selectPlanes(planes));
		final long[] bounds = findBounds(stack);
		final Stream<Vector3d> samplePoints = createSamplePoints(pairs, 1.0,
			bounds);

		double[] orientations = new double[3];

		samplePoints.limit(100).sequential().map(CountInterfaces::getOrientation).forEach(o -> {
            System.out.println(Arrays.toString(o));
            for (int i = 0; i < 3; i++) {
                orientations[i] += o[i];
            }
        });

        System.out.println(Arrays.toString(orientations));

		/*samplePoints.limit(1_000_000).sequential().forEach(point -> sample(stack, point));
		final IterableInterval<UnsignedShortType> iterable = Views.flatIterable(stack);
        final long[] voxelSamples = new long[(int) voxels];

		final Cursor<UnsignedShortType> cursor = iterable.cursor();
		int i = 0;
		while (cursor.hasNext()) {
			final long samples = cursor.next().get();
			voxelSamples[i] = samples;
			i++;
		}

        final SummaryStatistics statistics = new SummaryStatistics();
        Arrays.stream(voxelSamples).forEach(statistics::addValue);
        System.out.println(statistics.toString());

        final ImageJ imageJ = new ImageJ();
        imageJ.launch(args);
        imageJ.ui().show(stack);*/
    }
}
