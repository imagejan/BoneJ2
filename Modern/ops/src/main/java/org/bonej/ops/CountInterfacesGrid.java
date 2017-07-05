package org.bonej.ops;

import net.imagej.ImageJ;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.array.ArrayImg;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.basictypeaccess.array.IntArray;
import net.imglib2.type.BooleanType;
import net.imglib2.type.numeric.ComplexType;
import net.imglib2.type.numeric.integer.UnsignedIntType;
import org.scijava.vecmath.Vector3d;

import java.util.Arrays;
import java.util.Random;
import java.util.stream.Stream;

import static org.bonej.ops.CountInterfaces.toVoxelCoordinates;

/**
 * @author Richard Domander
 */
public class CountInterfacesGrid {
	private static Random random = new Random(0xc0ff33);

	public static <C extends ComplexType<C>> long[] findBounds(final RandomAccessibleInterval<C> interval) {
		final long[] bounds = new long[interval.numDimensions()];
		interval.dimensions(bounds);
		return bounds;
	}

	public static Stream<Vector3d> plotGrid(final long gridSize, final double increment, final long[] bounds) {
		// TODO add random offset
		// TODO plot and order grid starting points and direction vectors to xyz
		final Vector3d centroid = new Vector3d(bounds[0] * 0.5, bounds[1] * 0.5, bounds[2] * 0.5);
		final Vector3d gridStart = new Vector3d(centroid);
		final double[] offsetCoords = random.doubles(3, 0, increment - 1).toArray();
		final Vector3d offset = new Vector3d(offsetCoords);

		gridStart.sub(new Vector3d(gridSize * 0.5, gridSize * 0.5, gridSize * 0.5));
		final Stream.Builder<Vector3d> originBuilder = Stream.builder();
		for (double z = 0.0; z < gridSize; z += increment) {
			for (double y = 0.0; y < gridSize; y += increment) {
				for (double x = 0.0; x < gridSize; x += increment) {
					final Vector3d v = new Vector3d(gridStart);
					v.add(new Vector3d(x, y, z));
					v.add(offset);
					// TODO rotation
					originBuilder.add(v);
				}
			}
		}
		return originBuilder.build().filter(v -> !outOfBounds(v, bounds));
	}

	public static boolean outOfBounds(final Vector3d v, final long[] bounds) {
		return (v.x < 0) || (v.x > (bounds[0])) || (v.y < 0) ||
				(v.y > (bounds[1])) || (v.z < 0) || (v.z > (bounds[2]));
	}

	public static long findGridSize(final long[] bounds) {
		final long sumSquared = Arrays.stream(bounds).map(i -> i * i).sum();
		return (long) Math.ceil(Math.sqrt(sumSquared));
	}

	public static long[] toVoxelCoordinates(final Vector3d v) {
		final double[] coordinates = new double[3];
		v.get(coordinates);
		final long[] voxelCoordinates = new long[3];
		for (int i = 0; i < 3; i++) {
			voxelCoordinates[i] = (long) Math.floor(coordinates[i]);
		}
		return voxelCoordinates;
	}

	public static <C extends ComplexType<C>> void sample(final Vector3d v, final RandomAccessibleInterval<C> interval) {
		final long[] coordinates = toVoxelCoordinates(v);
		final RandomAccess<C> access = interval.randomAccess(interval);
		access.setPosition(coordinates);
		final C voxel = access.get();
		//System.out.println(Arrays.toString(coordinates));
		voxel.setReal(voxel.getRealDouble() + 1.0);
	}

	public static void main(String... args) {
		final RandomAccessibleInterval<UnsignedIntType> stack = ArrayImgs.unsignedInts(10, 10, 10);
		final long[] bounds = findBounds(stack);
		final long gridSize = findGridSize(bounds);
		final Stream<Vector3d> grid = plotGrid(gridSize, 2.0, bounds);
		grid.forEach(v -> sample(v, stack));
		final ImageJ imageJ = new ImageJ();
		imageJ.launch(args);
		imageJ.ui().show(stack);
	}
}
