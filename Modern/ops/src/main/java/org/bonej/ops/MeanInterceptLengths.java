
package org.bonej.ops;

import java.util.Arrays;
import java.util.Collection;
import java.util.Random;
import java.util.stream.Stream;

import net.imagej.ops.AbstractOp;
import net.imagej.ops.Op;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.BooleanType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.util.ValuePair;
import org.apache.commons.math3.random.UnitSphereRandomVectorGenerator;
import org.apache.commons.math3.util.FastMath;
import org.scijava.ItemIO;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.vecmath.Quat4d;
import org.scijava.vecmath.Tuple3d;
import org.scijava.vecmath.Vector3d;

/**
 * Returns a cloud of MIL vectors
 *
 * @author Richard Domander
 * @author Michael Doube
 */
// TODO Implement for 2D
// TODO implements Contingent
@Plugin(type = Op.class)
public class MeanInterceptLengths <B extends BooleanType<B>> extends AbstractOp {

	/**
	 * A random generator for unit quaternions
	 * <p>
	 * To create a uniformly sampled random unit quaternion which represents a
	 * uniformly sampled random rotation, the generator creates an array of four
	 * doubles, where each array element is a sampling of a normal distribution.
	 * It then treats the array as a vector and normalizes its length (a[i] = a[i]
	 * / ||a||).
	 * </p>
	 */
	private static final UnitSphereRandomVectorGenerator generator =
		new UnitSphereRandomVectorGenerator(4);
	private static final Random random = new Random();
	private static final int X = 0;
	private static final int Y = 1;
	private static final int Z = 2;

	@Parameter
	private RandomAccessibleInterval<B> interval;

	/** Number of vectors in the grid in each dimension */
	@Parameter
	private long vectors;

	/** Increment added to the scalar value t used to scale sampling vectors */
	@Parameter
	private double tIncrement = 1.0;

	@Parameter(type = ItemIO.OUTPUT)
	private Collection<Tuple3d> milVectors;

	@Override
	public void run() {
		final long[] bounds = findBounds(interval);
		final double gridSize = findGridSize(bounds);
		plotCuboid(gridSize, vectors, bounds);
	}

	public static <C extends NumericType<C>> long[] findBounds(
			final RandomAccessibleInterval<C> interval)
	{
		final long[] bounds = new long[interval.numDimensions()];
		interval.dimensions(bounds);
		return bounds;
	}

	public static double findGridSize(final long[] bounds) {
		final long sumSquared = Arrays.stream(bounds).map(i -> i * i).sum();
		return FastMath.sqrt(sumSquared);
	}

	public static Stream<ValuePair<Vector3d, Vector3d>> plotCuboid(
			final double gridSize, final long vectors, final long[] bounds)
	{
		final Vector3d centroid = new Vector3d(bounds[X], bounds[Y], bounds[Z]);
		centroid.scale(0.5);
		final double[] quaternion = generator.nextVector();
		final Stream.Builder<ValuePair<Vector3d, Vector3d>> builder = Stream.builder();
		final long origins = vectors * vectors;
		plotPlane(builder, gridSize, origins, X, Y, Z, quaternion, centroid);
		plotPlane(builder, gridSize, origins, X, Z, Y, quaternion, centroid);
		plotPlane(builder, gridSize, origins, Y, Z, X, quaternion, centroid);
		return builder.build();
	}

	public static void plotPlane(
			final Stream.Builder<ValuePair<Vector3d, Vector3d>> builder, final double gridSize,
			final long segments, final int dim0, final int dim1, final int dim2,
			final double[] quaternion, final Tuple3d centroid)
	{
		final int sign = random.nextBoolean() ? 1 : -1;
		final Vector3d normal = createNormal(dim2, sign, quaternion);
		final double planeShift = -0.5 * sign * gridSize;
		final double[] coordinates = new double[3];
		coordinates[dim2] = planeShift;
		final long points = (segments + 1) * (segments + 1);
		for (long i = 0; i < points; i++) {
			coordinates[dim0] = random.nextDouble() * gridSize - 0.5 * gridSize;
			coordinates[dim1] = random.nextDouble() * gridSize - 0.5 * gridSize;
			final Vector3d origin = createOrigin(coordinates, quaternion, centroid);
			builder.add(new ValuePair<>(origin, normal));
		}
	}

	public static Vector3d createOrigin(final double[] coordinates, final double[] quaternion, final Tuple3d centroid)
	{
		final Vector3d v = new Vector3d(coordinates);
		rotate(v, quaternion);
		v.add(centroid);
		return v;
	}

	public static Vector3d createNormal(final int dim2, final int sign,
										 final double[] quaternion)
	{
		final double[] n = new double[3];
		n[dim2] = sign * 1.0;
		final Vector3d normal = new Vector3d(n);
		rotate(normal, quaternion);
		return normal;
	}

	public static void rotate(final Vector3d v, final double[] quaternion) {
		final Quat4d p = new Quat4d();
		final Quat4d qInv = new Quat4d();
		final Quat4d rotated = new Quat4d();
		rotated.set(quaternion);
		p.set(v.x, v.y, v.z, 0.0);
		rotated.mul(p);
		qInv.set(-quaternion[0], -quaternion[1], -quaternion[2], quaternion[3]);
		rotated.mul(qInv);
		v.set(rotated.x, rotated.y, rotated.z);
	}
}
