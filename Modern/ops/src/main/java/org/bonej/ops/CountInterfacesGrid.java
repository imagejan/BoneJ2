
package org.bonej.ops;

import static org.bonej.ops.CountInterfacesGrid.SamplingPlane.Orientation.XY;
import static org.bonej.ops.CountInterfacesGrid.SamplingPlane.Orientation.XZ;
import static org.bonej.ops.CountInterfacesGrid.SamplingPlane.Orientation.YZ;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.stream.Stream;
import java.util.stream.Stream.Builder;

import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.BooleanType;
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

	/**
	 * Create a four-element vector where each element is a sampling of a normal
	 * distribution. Normalize its length and you have a uniformly sampled random
	 * unit quaternion which represents a uniformly sampled random rotation.
	 */
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

	public static <B extends BooleanType<B>> long sample(
		final RandomAccessibleInterval<B> interval, final long[] bounds,
		final ValuePair<Vector3d, Vector3d> sampler, final double increment)
	{
		final ValuePair<Double, Double> tPair = findIntersections(sampler, bounds);
		if (tPair == null) {
			return 0;
		}
		final double jitterT = Math.random() * increment;
		final double tMin = FastMath.min(tPair.a, tPair.b);
		final double tMax = FastMath.max(tPair.a, tPair.b) - jitterT;
		if (tMin > tMax) {
			// TODO what about just one sample point in the [tMin, tMax] interval?
			return 0;
		}
		final Vector3d origin = sampler.a;
		final Vector3d direction = sampler.b;
		final Vector3d jitter = new Vector3d(direction);
		jitter.scale(jitterT);
		final Vector3d point = new Vector3d();
		point.scaleAdd(tMin, direction, origin);
		point.add(jitter);
		final Vector3d bit = new Vector3d(direction);
		bit.scale(increment);
		final int[] coordinates = new int[3];
		final RandomAccess<B> access = interval.randomAccess();
		long previous = getElement(access, point, coordinates);
		long intersections = 0;
		for (double t = tMin; t <= tMax; t += increment) {
			final long current = getElement(access, point, coordinates);
			intersections += current ^ previous;
			point.add(bit);
			previous = current;
		}
		return intersections;
	}

	public static <B extends BooleanType<B>> long getElement(
		final RandomAccess<B> access, final Vector3d position,
		final int[] coordinates)
	{
		coordinates[0] = (int) position.x;
		coordinates[1] = (int) position.y;
		coordinates[2] = (int) position.z;
		access.setPosition(coordinates);
		return (long) access.get().getRealDouble();
	}

	// TODO Unnecessary, since ops inherit run from Runnable? (Check the forum)
	public static <B extends BooleanType<B>> long futureParallel(
		final RandomAccessibleInterval<B> interval, final long[] bounds,
		final double gridSize, final long sections, final double increment,
		final int rotations) throws ExecutionException, InterruptedException
	{
		final int processors = 5;
		final ExecutorService executors = Executors.newFixedThreadPool(processors);
		final List<Future<Long>> futures = new ArrayList<>();
		final Vector3d centroid = new Vector3d(bounds[0], bounds[1], bounds[2]);
		centroid.scale(0.5);
		final long points = (sections + 1) * (sections + 1);
		for (int i = 0; i < rotations; i++) {
			final double[] quaternion = generator.nextVector();
			final SamplingGrid grid = new SamplingGrid(gridSize, quaternion,
				centroid);
			final Stream<ValuePair<Vector3d, Vector3d>> samplers = grid.getSamplers(
				points);
			final Future<Long> future = executors.submit(() -> samplers.mapToLong(
				sampler -> sample(interval, bounds, sampler, increment)).sum());
			futures.add(future);
		}

		long total = 0;
		for (final Future<Long> future : futures) {
			total += future.get();
		}

		executors.shutdown();
		return total;
	}

	private static final class SamplingGrid {

		public final SamplingPlane xy;
		public final SamplingPlane xz;
		public final SamplingPlane yz;

		public SamplingGrid(final double gridSize, final double[] quaternion,
			final Tuple3d centroid)
		{
			xy = new SamplingPlane(gridSize, XY, quaternion, centroid);
			xz = new SamplingPlane(gridSize, XZ, quaternion, centroid);
			yz = new SamplingPlane(gridSize, YZ, quaternion, centroid);
		}

		private Stream<ValuePair<Vector3d, Vector3d>> getSamplers(final long n) {
			final Stream.Builder<ValuePair<Vector3d, Vector3d>> builder = Stream
				.builder();
			for (long i = 0; i < n; i++) {
				builder.accept(xy.getSampler());
				builder.accept(xz.getSampler());
				builder.accept(yz.getSampler());
			}
			return builder.build();
		}
	}

	public static final class SamplingPlane {

		public enum Orientation {
				XY, XZ, YZ
		}

		private final Vector3d translation;
		private final Vector3d u;
		private final Vector3d v;
		private final Vector3d normal;

		// TODO Fix sign
		public SamplingPlane(final double gridSize, final Orientation orientation,
			final double[] quaternion, final Tuple3d centroid)
		{
			final int direction = random.nextBoolean() ? 1 : -1;
			final double halfGrid = -0.5 * gridSize;
			final double planeShift = direction * halfGrid;
			if (orientation == XY) {
				u = new Vector3d(gridSize, 0, 0);
				v = new Vector3d(0, gridSize, 0);
				normal = new Vector3d(0, 0, direction);
				translation = new Vector3d(halfGrid, halfGrid, planeShift);
			}
			else if (orientation == XZ) {
				u = new Vector3d(gridSize, 0, 0);
				normal = new Vector3d(0, direction, 0);
				v = new Vector3d(0, 0, gridSize);
				translation = new Vector3d(halfGrid, planeShift, halfGrid);
			}
			else if (orientation == YZ) {
				normal = new Vector3d(direction, 0, 0);
				u = new Vector3d(0, gridSize, 0);
				v = new Vector3d(0, 0, gridSize);
				translation = new Vector3d(planeShift, halfGrid, halfGrid);
			} else {
				throw new IllegalArgumentException("Unexpected orientation");
			}

			rotateVectors(quaternion);
			translation.add(centroid);
		}

		public ValuePair<Vector3d, Vector3d> getSampler() {
			final double uScale = random.nextDouble();
			final double vScale = random.nextDouble();
			final double x = uScale * u.x + vScale * v.x;
			final double y = uScale * u.y + vScale * v.y;
			final double z = uScale * u.z + vScale * v.z;
			final Vector3d origin = new Vector3d(x, y, z);
			origin.add(translation);
			return new ValuePair<>(origin, normal);
		}

		private void rotateVectors(final double[] quaternion) {
			rotate(translation, quaternion);
			rotate(u, quaternion);
			rotate(v, quaternion);
			rotate(normal, quaternion);
		}
	}

	public static Stream<ValuePair<Vector3d, Vector3d>> plotGrid(
		final double gridSize, final long sections, final long[] bounds)
	{
		final Vector3d centroid = new Vector3d(bounds[0], bounds[1], bounds[2]);
		centroid.scale(0.5);
		final double[] quaternion = generator.nextVector();
		final Builder<ValuePair<Vector3d, Vector3d>> builder = Stream.builder();
		plotPlane(builder, gridSize, sections, 0, 1, 2, quaternion, centroid);
		plotPlane(builder, gridSize, sections, 0, 2, 1, quaternion, centroid);
		plotPlane(builder, gridSize, sections, 1, 2, 0, quaternion, centroid);
		return builder.build();
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

	public static void plotPlane(
		final Builder<ValuePair<Vector3d, Vector3d>> builder, final double gridSize,
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
			final Vector3d origin = createOrigin(coordinates, dim0, dim1, gridSize,
				quaternion, centroid);
			builder.add(new ValuePair<>(origin, normal));
		}
	}

	public static Vector3d createOrigin(final double[] coordinates,
		final int dim0, final int dim1, final double gridSize,
		final double[] quaternion, final Tuple3d centroid)
	{
		coordinates[dim0] = random.nextDouble() * gridSize - 0.5 * gridSize;
		coordinates[dim1] = random.nextDouble() * gridSize - 0.5 * gridSize;
		final Vector3d v = new Vector3d(coordinates);
		rotate(v, quaternion);
		v.add(centroid);
		return v;
	}

	private static Vector3d createNormal(final int dim2, final int sign,
		final double[] quaternion)
	{
		final double[] n = new double[3];
		n[dim2] = sign * 1.0;
		final Vector3d normal = new Vector3d(n);
		rotate(normal, quaternion);
		return normal;
	}

	public static double findGridSize(final long[] bounds) {
		final long sumSquared = Arrays.stream(bounds).map(i -> i * i).sum();
		return FastMath.sqrt(sumSquared);
	}

	public static ValuePair<Double, Double> findIntersections(
		final ValuePair<Vector3d, Vector3d> sampler, final long[] bounds)
	{
		final Vector3d origin = sampler.a;
		final Vector3d direction = sampler.b;
		// Subtract epsilon from bounds, so that vectors don't intersect the stack
		// at w,h or d which would cause an index out of bounds exception
		final double eps = 1e-12;
		final double minX = direction.x >= 0.0 ? 0 : bounds[0] - eps;
		final double maxX = direction.x >= 0.0 ? bounds[0] - eps : 0;
		final double minY = direction.y >= 0.0 ? 0 : bounds[1] - eps;
		final double maxY = direction.y >= 0.0 ? bounds[1] - eps : 0;
		final double minZ = direction.z >= 0.0 ? 0 : bounds[2] - eps;
		final double maxZ = direction.z >= 0.0 ? bounds[2] - eps : 0;
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
		}
		else if (Double.isNaN(a)) {
			return b;
		}
		else if (Double.isNaN(b)) {
			return a;
		}
		else {
			return FastMath.max(a, b);
		}
	}

	public static double minNan(final double a, final double b) {
		if (Double.isNaN(a) && Double.isNaN(b)) {
			return Double.NaN;
		}
		else if (Double.isNaN(a)) {
			return b;
		}
		else if (Double.isNaN(b)) {
			return a;
		}
		else {
			return FastMath.min(a, b);
		}
	}
}
