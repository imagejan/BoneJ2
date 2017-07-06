
package org.bonej.ops;

import static org.bonej.ops.CountInterfacesGrid.findBounds;
import static org.bonej.ops.CountInterfacesGrid.findGridSize;
import static org.bonej.ops.CountInterfacesGrid.outOfBounds;
import static org.bonej.ops.CountInterfacesGrid.plotSamplers;
import static org.bonej.ops.CountInterfacesGrid.sample;
import static org.bonej.ops.CountInterfacesGrid.samplePoints;

import java.util.Arrays;
import java.util.stream.Stream;

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
		final Stream<ValuePair<Vector3d, Vector3d>> samplers = plotSamplers(
			gridSize, 1L, bounds);
		// TODO Rotation broken
		final Stream<Vector3d> points = samplePoints(samplers, 1.0, bounds);
		points.map(CountInterfacesGrid::toVoxelCoordinates).peek(c -> System.out
			.println(Arrays.toString(c))).filter(c -> !outOfBounds(c, bounds))
			.forEach(c -> sample(stack, c));
		/*final ImageJ imageJ = new ImageJ();
		imageJ.launch(args);
		imageJ.ui().show(stack);*/
	}
}
