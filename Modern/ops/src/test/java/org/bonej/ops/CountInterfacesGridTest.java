
package org.bonej.ops;

import static org.bonej.ops.CountInterfacesGrid.findBounds;
import static org.bonej.ops.CountInterfacesGrid.findGridSize;
import static org.bonej.ops.CountInterfacesGrid.outOfBounds;
import static org.bonej.ops.CountInterfacesGrid.plotSamplers;
import static org.bonej.ops.CountInterfacesGrid.sample;
import static org.bonej.ops.CountInterfacesGrid.samplePoints;

import java.util.Arrays;
import java.util.stream.Stream;

import net.imagej.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.type.numeric.integer.UnsignedIntType;
import net.imglib2.util.ValuePair;

import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
public class CountInterfacesGridTest {

	public static void main(String... args) {
		final RandomAccessibleInterval<UnsignedIntType> stack = ArrayImgs
			.unsignedInts(2, 2, 2);
		final long[] bounds = findBounds(stack);
		final long gridSize = findGridSize(bounds);

		/*
		final Stream.Builder<ValuePair<Vector3d, Vector3d>> samplingBuilder = Stream
				.builder();
		plotPlane(samplingBuilder, bounds, gridSize, 2, 0, 1);
		plotPlane(samplingBuilder, bounds, gridSize, 2, 0, 2);
		plotPlane(samplingBuilder, bounds, gridSize, 2, 1, 2);
		final double[] angles = {0, 0, 0}; //random.doubles(3, 0, Math.PI).toArray();
		final Vector3d centroid = new Vector3d(bounds[0] - 1, bounds[1] - 1, bounds[2] - 1);
		centroid.scale(0.5);
		samplingBuilder.build().map(p -> rotateSampler(p, angles)).map(p -> { p.a.add(centroid); return p; }).forEach(p -> System.out.println(p.a.toString()));
		    */

		final Stream<ValuePair<Vector3d, Vector3d>> samplers = plotSamplers(
			gridSize, 2L, bounds);
		//samplers.peek(p -> System.out.print(p.a + " " + p.b + " ")).map(p -> CountInterfacesGrid.findIntersections(p, bounds)).forEach(p -> System.out.println(p == null));
		// TODO rotation broken
		// TODO (Normal) vector rounding causes some vectors to not intersect the stack, that is, don't sample from (0,0,0) along (0,0,1)
		final Stream<Vector3d> points = samplePoints(samplers, 1.0, bounds);

		points.peek(System.out::println).map(
			CountInterfacesGrid::toVoxelCoordinates).peek(c -> System.out.println(
				Arrays.toString(c))).filter(c -> !outOfBounds(c, bounds)).forEach(
					c -> sample(stack, c));
		final ImageJ imageJ = new ImageJ();
		imageJ.launch(args);
		imageJ.ui().show(stack);
	}
}
