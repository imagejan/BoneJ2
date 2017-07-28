
package org.bonej.ops;

import java.util.Arrays;
import java.util.Random;
import java.util.stream.Stream;

import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.util.ValuePair;

import org.apache.commons.math3.random.UnitSphereRandomVectorGenerator;
import org.apache.commons.math3.util.FastMath;
import org.scijava.vecmath.Quat4d;
import org.scijava.vecmath.Tuple3d;
import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
public class CountInterfacesGrid {

	private static final UnitSphereRandomVectorGenerator generator =
		new UnitSphereRandomVectorGenerator(4);
	private static final Random random = new Random(System.currentTimeMillis() +
		System.identityHashCode(new Object()));

	public static <C extends NumericType<C>> long[] findBounds(
		final RandomAccessibleInterval<C> interval)
	{
		final long[] bounds = new long[interval.numDimensions()];
		interval.dimensions(bounds);
		return bounds;
	}

	public static Stream<ValuePair<Vector3d, Vector3d>> plotSamplers(
		final double gridSize, final long sections, final long[] bounds)
	{
		final Vector3d centroid = new Vector3d(bounds[0], bounds[1], bounds[2]);
		centroid.scale(0.5);
		// create a four-element vector where each element is a sampling of a normal
		// distribution. Normalize its length and you have a uniformly sampled
		// random unit quaternion which represents a uniformly sampled random
		// rotation.
		final double[] quaternion = generator.nextVector();
		Stream.Builder<ValuePair<Vector3d, Vector3d>> builder = Stream.builder();
		plotPlane(builder, gridSize, sections, 0, 1, 2, quaternion, centroid);
		plotPlane(builder, gridSize, sections, 0, 2, 1, quaternion, centroid);
		plotPlane(builder, gridSize, sections, 1, 2, 0, quaternion, centroid);
		return builder.build();
	}

	public static void rotateXYZ(final Vector3d v, final double[] q) {
		final Quat4d p = new Quat4d();
		final Quat4d qInv = new Quat4d();
		final Quat4d rotated = new Quat4d();
		p.set(v.x, v.y, v.z, 0.0);
		qInv.set(-q[0], -q[1], -q[2], q[3]);
		rotated.set(q);
		rotated.mul(p);
		rotated.mul(qInv);
		v.set(rotated.x, rotated.y, rotated.z);
	}

	public static void plotPlane(
		final Stream.Builder<ValuePair<Vector3d, Vector3d>> builder,
		final double gridSize, final long segments, final int dim0, final int dim1,
		final int dim2, final double[] quaternion, final Tuple3d centroid)
	{
		// TODO fix sign
		final Vector3d normal = createNormal(dim2);
		rotateXYZ(normal, quaternion);
		final double halfGrid = -0.5 * gridSize;
		final double[] coordinates = new double[3];
		coordinates[dim2] = halfGrid;
		final long points = (segments + 1) * (segments + 1);
		for (long i = 0; i < points; i++) {
			final Vector3d v = nextVector(coordinates, dim0, dim1, gridSize, halfGrid,
				quaternion, centroid);
			builder.add(new ValuePair<>(v, normal));
		}
	}

	public static Vector3d nextVector(final double[] coordinates, final int dim0,
		final int dim1, final double gridSize, final double halfGrid,
		final double[] quaternion, final Tuple3d centroid)
	{
		coordinates[dim0] = random.nextDouble() * gridSize + halfGrid;
		coordinates[dim1] = random.nextDouble() * gridSize + halfGrid;
		final Vector3d v = new Vector3d(coordinates);
		rotateXYZ(v, quaternion);
		v.add(centroid);
		return v;
	}

	private static Vector3d createNormal(final int dim2) {
		final double[] n = new double[3];
		n[dim2] = 1.0;
		return new Vector3d(n);
	}

	public static boolean outOfBounds(final int x, final int y, final int z,
									  final long[] bounds)
	{
		return (x < 0) || (x >= (bounds[0])) ||
				(y < 0) || (y >= (bounds[1])) ||
				(z < 0) || (z >= (bounds[2]));
	}

	public static boolean outOfBounds(final int[] coordinates,
		final long[] bounds)
	{
		return (coordinates[0] < 0) || (coordinates[0] >= (bounds[0])) ||
			(coordinates[1] < 0) || (coordinates[1] >= (bounds[1])) ||
			(coordinates[2] < 0) || (coordinates[2] >= (bounds[2]));
	}

	public static double findGridSize(final long[] bounds) {
		final long sumSquared = Arrays.stream(bounds).map(i -> i * i).sum();
		return Math.sqrt(sumSquared);
	}

	public static ValuePair<Double, Double> findIntersections(
		final ValuePair<Vector3d, Vector3d> sampler, final long[] bounds)
	{
		// Make min coordinates slightly negative to make vectors going parallel to
		// stack planes intersect, e.g. origin (0, 0, 0) and direction (0, 0, 1)
		final Vector3d origin = sampler.a;
		final Vector3d direction = sampler.b;
		final double minX = direction.x >= 0.0 ? 0 : bounds[0];
		final double maxX = direction.x >= 0.0 ? bounds[0] : 0;
		final double minY = direction.y >= 0.0 ? 0 : bounds[1];
		final double maxY = direction.y >= 0.0 ? bounds[1] : 0;
		final double minZ = direction.z >= 0.0 ? 0 : bounds[2];
		final double maxZ = direction.z >= 0.0 ? bounds[2]: 0;
		final double tX0 = (minX - origin.x) / direction.x;
		final double tX1 = (maxX - origin.x) / direction.x;
		final double tY0 = (minY - origin.y) / direction.y;
		final double tY1 = (maxY - origin.y) / direction.y;
		final double tZ0 = (minZ - origin.z) / direction.z;
		final double tZ1 = (maxZ - origin.z) / direction.z;
		if (tX0 > tY1 || tY0 > tX1) {
			return null;
		}
		double tMin = maxNan(tX0, tY0);
		double tMax = minNan(tX1, tY1);
		if (tMin > tZ1 || tZ0 > tMax) {
			return null;
		}
		tMin = maxNan(tZ0, tMin);
		tMax = minNan(tZ1, tMax);
		if (Double.isNaN(tMin) || Double.isNaN(tMax)) {
			return null;
		}

		return new ValuePair<>(tMin, tMax);
	}

	public static double maxNan(final double a, final double b) {
		if (Double.isNaN(a) && Double.isNaN(b)) {
			return Double.NaN;
		} else if (Double.isNaN(a)) {
			return b;
		} else if (Double.isNaN(b)) {
			return a;
		} else {
			return FastMath.max(a, b);
		}
	}

	public static double minNan(final double a, final double b) {
		if (Double.isNaN(a) && Double.isNaN(b)) {
			return Double.NaN;
		} else if (Double.isNaN(a)) {
			return b;
		} else if (Double.isNaN(b)) {
			return a;
		} else {
			return FastMath.min(a, b);
		}
	}
}
