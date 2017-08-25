
package org.bonej.ops;

import static org.bonej.ops.MeanInterceptLengths.SamplingPlane.Orientation.XY;
import static org.bonej.ops.MeanInterceptLengths.SamplingPlane.Orientation.XZ;
import static org.bonej.ops.MeanInterceptLengths.SamplingPlane.Orientation.YZ;

import java.util.Arrays;
import java.util.Collection;
import java.util.Objects;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import net.imagej.ops.AbstractOp;
import net.imagej.ops.Op;
import net.imagej.ops.special.function.AbstractUnaryFunctionOp;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.BooleanType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.util.ValuePair;

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
public class MeanInterceptLengths<B extends BooleanType<B>> extends AbstractUnaryFunctionOp<RandomAccessibleInterval<B>, Collection<Tuple3d>> {

	private static final Random random = new Random();

	/** Number of sampling vectors on each plane */
	@Parameter(required = false)
	private Long vectors = null;

	/** Increment added to the scalar value t used to scale sampling vectors */
	@Parameter(required = false)
	private double tIncrement = 1.0 / 2.3;

	@Parameter(type = ItemIO.BOTH, required = false)
	private Quat4d rotation = null;

	@Override
	public Collection<Tuple3d> calculate(final RandomAccessibleInterval<B> interval) {
		final long[] bounds = findBounds(interval);
		final double gridSize = findGridSize(bounds);
		if (vectors == null) {
			vectors = Arrays.stream(bounds).max().orElse(0) + 1;
			vectors *= vectors;
		}
		if (rotation == null) {
			// TODO Test how ItemIO.BOTH works
			rotation = nextQuaternion();
		}
		final Vector3d centroid = findCentroid(bounds);
		final SamplingGrid grid = new SamplingGrid(gridSize, rotation, centroid);
		final Stream<ValuePair<Vector3d, Vector3d>> samplers = grid.getSamplers(
				vectors);
		return samplers.map(sampler -> sample(interval, bounds, sampler,
				tIncrement)).filter(Objects::nonNull).collect(Collectors.toList());
	}

	public static Vector3d findCentroid(final long[] bounds) {
		return new Vector3d(bounds[0] * 0.5, bounds[1] * 0.5, bounds[2] * 0.5);
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
		return Math.sqrt(sumSquared);
	}

	public static <B extends BooleanType<B>> Vector3d sample(
		final RandomAccessibleInterval<B> interval, final long[] bounds,
		final ValuePair<Vector3d, Vector3d> sampler, final double increment)
	{
		// TODO Refactor to findMILVector(sampler) { intersections = intersections ... return milVector;}
		final ValuePair<Double, Double> tPair = findIntersections(sampler, bounds);
		if (tPair == null) {
			return null;
		}
		final double jitterT = random.nextDouble() * increment;
		final double tMin = Math.min(tPair.a, tPair.b);
		final double tMax = Math.max(tPair.a, tPair.b) - jitterT;
		if (tMin > tMax) {
			// interval fits (at most) one sampling point. no intersections
			return null;
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
		final RandomAccess<B> access = interval.randomAccess();
		long previous = getElement(access, point);
		long intersections = 0;
		for (double t = tMin; t <= tMax; t += increment) {
			final long current = getElement(access, point);
			intersections += current ^ previous;
			point.add(bit);
			previous = current;
		}
		return toMILVector(direction, tPair, intersections);
	}

	public static Vector3d toMILVector(
		final Vector3d direction,
		final ValuePair<Double, Double> tPair, final long intersections)
	{
		if (intersections == 0) {
			return direction;
		}
		final double t = Math.abs(tPair.a - tPair.b) / intersections;
		final Vector3d milVector = new Vector3d(direction);
		milVector.scale(t);
		return milVector;
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
			return Math.max(a, b);
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
			return Math.min(a, b);
		}
	}

	public static <B extends BooleanType<B>> long getElement(
		final RandomAccess<B> access, final Vector3d position)
	{
		access.setPosition((int) position.x, 0);
		access.setPosition((int) position.y, 1);
		access.setPosition((int) position.z, 2);
		return (long) access.get().getRealDouble();
	}

	// region -- Helper methods --

	/**
	 * Creates a random unit quaternion which represents a rotation
	 * <p>
	 * The quaternion is first assigned four normally distributed coordinates with
	 * mean {@code 0.0} and standard deviation {@code 1.0}. Then it's normalised.
	 * This method guarantees that the unit quaternions are uniformly sampled, and
	 * they represent uniformly sampled random rotations.
	 * </p>
	 */
	private Quat4d nextQuaternion() {
		final double x = random.nextGaussian();
		final double y = random.nextGaussian();
		final double z = random.nextGaussian();
		final double w = random.nextGaussian();
		return new Quat4d(x, y, z, w);
	}

	// endregion

	// region -- Helper classes --

	public static final class SamplingGrid {

		private final SamplingPlane xy;
		private final SamplingPlane xz;
		private final SamplingPlane yz;

		public SamplingGrid(final double gridSize, final Quat4d quaternion,
			final Tuple3d centroid)
		{
			xy = new SamplingPlane(gridSize, XY, quaternion, centroid);
			xz = new SamplingPlane(gridSize, XZ, quaternion, centroid);
			yz = new SamplingPlane(gridSize, YZ, quaternion, centroid);
		}

		public Stream<ValuePair<Vector3d, Vector3d>> getSamplers(final long n) {
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

		private SamplingPlane(final double gridSize, final Orientation orientation,
			final Quat4d quaternion, final Tuple3d centroid)
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
			}
			else {
				throw new IllegalArgumentException("Unexpected orientation");
			}

			rotateVectors(quaternion);
			translation.add(centroid);
		}

		// TODO Op?
		private static void rotate(final Vector3d v, final Quat4d q) {
			final Quat4d p = new Quat4d();
			final Quat4d qInv = new Quat4d();
			final Quat4d rotated = new Quat4d(q);
			p.set(v.x, v.y, v.z, 0.0);
			rotated.mul(p);
			// No need to divide q^-1 since q is a unit quaternion
			qInv.set(-q.x, -q.y, -q.z, q.w);
			rotated.mul(qInv);
			v.set(rotated.x, rotated.y, rotated.z);
		}

		private ValuePair<Vector3d, Vector3d> getSampler() {
			final double uScale = random.nextDouble();
			final double vScale = random.nextDouble();
			final double x = uScale * u.x + vScale * v.x;
			final double y = uScale * u.y + vScale * v.y;
			final double z = uScale * u.z + vScale * v.z;
			final Vector3d origin = new Vector3d(x, y, z);
			origin.add(translation);
			return new ValuePair<>(origin, normal);
		}

		private void rotateVectors(final Quat4d quaternion) {
			rotate(translation, quaternion);
			rotate(u, quaternion);
			rotate(v, quaternion);
			rotate(normal, quaternion);
		}
	}
	// endregion
}
