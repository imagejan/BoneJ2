
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

import net.imagej.ops.Contingent;
import net.imagej.ops.Op;
import net.imagej.ops.special.function.AbstractUnaryFunctionOp;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.BooleanType;
import net.imglib2.type.numeric.NumericType;
import net.imglib2.util.ValuePair;

import org.apache.commons.math3.complex.Quaternion;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.random.UnitSphereRandomVectorGenerator;
import org.scijava.ItemIO;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

/**
 * Returns a cloud of MIL vectors
 *
 * @author Richard Domander
 * @author Michael Doube
 */
// TODO Implement for 2D
@Plugin(type = Op.class)
public class MeanInterceptLengths<B extends BooleanType<B>> extends
	AbstractUnaryFunctionOp<RandomAccessibleInterval<B>, Collection<Vector3D>>
	implements Contingent
{

	/**
	 * Generates uniformly sampled random unit quaternions, which represent
	 * uniformly sampled random rotations.
	 */
	private static final UnitSphereRandomVectorGenerator generator =
		new UnitSphereRandomVectorGenerator(4);
	private static final Random random = new Random();

	/** Number of sampling vectors on each plane. */
	@Parameter(required = false)
	private Long vectors = null;

	/** Increment added to the scalar value t used to scale sampling vectors. */
	@Parameter(required = false)
	private Double tIncrement = null;

	/**
	 * The unit quaternion used to rotate the sampling grid.
	 * <p>
	 * If left null, then a random rotation is used.
	 * </p>
	 */
	@Parameter(type = ItemIO.BOTH, required = false)
	private double[] rotation = null;

	@Override
	public Collection<Vector3D> calculate(
		final RandomAccessibleInterval<B> interval)
	{
		final long[] bounds = findBounds(interval);
		final double gridSize = findGridSize(bounds);
		fillParameters(bounds);
		final Vector3D centroid = findCentroid(bounds);
		final SamplingGrid grid = new SamplingGrid(gridSize, rotation, centroid);
		final Stream<ValuePair<Vector3D, Vector3D>> samplers = grid.getSamplers(
			vectors);
		return samplers.map(sampler -> sampleMILVector(interval, bounds, sampler,
			tIncrement)).filter(Objects::nonNull).collect(Collectors.toList());
	}

	// TODO Check RAI dimensionality
	@Override
	public boolean conforms() {
		return rotation == null || rotation.length == 4;
	}

	// region -- Helper methods (public static for unit testing) --
	public static <C extends NumericType<C>> long[] findBounds(
		final RandomAccessibleInterval<C> interval)
	{
		final long[] bounds = new long[interval.numDimensions()];
		interval.dimensions(bounds);
		return bounds;
	}

	public static Vector3D findCentroid(final long[] bounds) {
		return new Vector3D(bounds[0] * 0.5, bounds[1] * 0.5, bounds[2] * 0.5);
	}

	public static double findGridSize(final long[] bounds) {
		final long sumSquared = Arrays.stream(bounds).map(i -> i * i).sum();
		return Math.sqrt(sumSquared);
	}

	public static ValuePair<Double, Double> findIntervalIntersections(
		final ValuePair<Vector3D, Vector3D> sampler, final long[] bounds)
	{
		final Vector3D origin = sampler.a;
		final Vector3D direction = sampler.b;
		// Subtract epsilon from bounds, so that vectors don't intersect the stack
		// at w,h or d which would cause an index out of bounds exception
		final double eps = 1e-12;
		final double minX = direction.getX() >= 0.0 ? 0 : bounds[0] - eps;
		final double maxX = direction.getX() >= 0.0 ? bounds[0] - eps : 0;
		final double minY = direction.getY() >= 0.0 ? 0 : bounds[1] - eps;
		final double maxY = direction.getY() >= 0.0 ? bounds[1] - eps : 0;
		final double minZ = direction.getZ() >= 0.0 ? 0 : bounds[2] - eps;
		final double maxZ = direction.getZ() >= 0.0 ? bounds[2] - eps : 0;
		final double tX0 = (minX - origin.getX()) / direction.getX();
		final double tX1 = (maxX - origin.getX()) / direction.getX();
		final double tY0 = (minY - origin.getY()) / direction.getY();
		final double tY1 = (maxY - origin.getY()) / direction.getY();
		final double tZ0 = (minZ - origin.getZ()) / direction.getZ();
		final double tZ1 = (maxZ - origin.getZ()) / direction.getZ();
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

	public static <B extends BooleanType<B>> long getElement(
		final RandomAccess<B> access, final Vector3D position)
	{
		access.setPosition((int) position.getX(), 0);
		access.setPosition((int) position.getY(), 1);
		access.setPosition((int) position.getZ(), 2);
		return (long) access.get().getRealDouble();
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
		return Math.max(a, b);
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
		return Math.min(a, b);
	}

	public static <B extends BooleanType<B>> long sampleIntercepts(
		final RandomAccessibleInterval<B> interval,
		final ValuePair<Vector3D, Vector3D> sampler,
		final ValuePair<Double, Double> intersections, final double increment)
	{
		final double offset = random.nextDouble() * increment;
		final double tMin = Math.min(intersections.a, intersections.b);
		final double tMax = Math.max(intersections.a, intersections.b) - offset;
		if (tMin > tMax) {
			// interval fits (at most) one sampling point. no intercepts
			return -1;
		}
		final Vector3D origin = sampler.a;
		final Vector3D direction = sampler.b;
		final Vector3D samplingGap = direction.scalarMultiply(increment);
		final RandomAccess<B> access = interval.randomAccess();
		final Vector3D point = direction.scalarMultiply(tMin + offset).add(origin);
		final long iterations = (long) Math.max(0, (tMax - tMin) / increment + 1);
		return countIntercepts(access, point, samplingGap, iterations);
	}

	public static <B extends BooleanType<B>> Vector3D sampleMILVector(
		final RandomAccessibleInterval<B> interval, final long[] stack,
		final ValuePair<Vector3D, Vector3D> sampler, final double increment)
	{
		final ValuePair<Double, Double> intersections = findIntervalIntersections(
			sampler, stack);
		if (intersections == null) {
			return null;
		}
		final long intercepts = sampleIntercepts(interval, sampler, intersections,
			increment);
		if (intercepts < 0) {
			return null;
		}
		return toMILVector(sampler.b, intersections, intercepts);
	}

	/**
	 * Creates a mean interception length (MIL) vector.
	 * <p>
	 * The MIL vector is parallel to the direction of the sampling vector. Its
	 * length equals the section of the sampling vector withing the stack
	 * (interval), and its origin is (0, 0, 0).
	 * </p>
	 *
	 * @param direction direction of the sampling vector.
	 * @param intersections scalar values (t) for stack - sampling vector
	 *          intersection points.
	 * @param intercepts number of foreground - background interceptions along
	 *          sampling vector.
	 * @return a new MIL vector.
	 */
	// TODO rename
	public static Vector3D toMILVector(final Vector3D direction,
		final ValuePair<Double, Double> intersections, final long intercepts)
	{
		final double intervalScalar = Math.abs(intersections.a - intersections.b);
		if (intercepts < 2) {
			return direction.scalarMultiply(intervalScalar);
		}
		return direction.scalarMultiply(intervalScalar / intercepts);
	}

	private static <B extends BooleanType<B>> long countIntercepts(
		final RandomAccess<B> access, final Vector3D start,
		final Vector3D samplingGap, final long iterations)
	{
		long intercepts = 0;
		long previous = getElement(access, start);
		Vector3D point = start;
		for (long i = 0; i < iterations; i++) {
			final long current = getElement(access, start);
			intercepts += current ^ previous;
			point = point.add(samplingGap);
			previous = current;
		}
		return intercepts;
	}

	private void fillParameters(final long[] bounds) {
		if (vectors == null) {
			vectors = Arrays.stream(bounds).min().orElse(0);
			vectors *= vectors;
		}
		if (rotation == null) {
			// TODO Test how ItemIO.BOTH works
			rotation = generator.nextVector();
		}
		if (tIncrement == null) {
			tIncrement = 1.0;
		}
	}

	// endregion

	// region -- Helper classes --

	public static final class SamplingGrid {

		private final SamplingPlane xy;
		private final SamplingPlane xz;
		private final SamplingPlane yz;

		public SamplingGrid(final double gridSize, final double[] quaternion,
			final Vector3D centroid)
		{
			xy = new SamplingPlane(gridSize, XY, quaternion, centroid);
			xz = new SamplingPlane(gridSize, XZ, quaternion, centroid);
			yz = new SamplingPlane(gridSize, YZ, quaternion, centroid);
		}

		public Stream<ValuePair<Vector3D, Vector3D>> getSamplers(final long n) {
			final Stream.Builder<ValuePair<Vector3D, Vector3D>> builder = Stream
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

		private Vector3D translation;
		private Vector3D u;
		private Vector3D v;
		private Vector3D normal;

		private SamplingPlane(final double gridSize, final Orientation orientation,
			final double[] quaternion, final Vector3D centroid)
		{
			createVectors(orientation, gridSize);
			rotateVectors(quaternion);
			translation = translation.add(centroid);
		}

		private void createVectors(final Orientation orientation,
			final double gridSize)
		{
			final int orthogonalDirection = random.nextBoolean() ? 1 : -1;
			final double halfGrid = -0.5 * gridSize;
			final double orthogonalShift = orthogonalDirection * halfGrid;
			if (orientation == XY) {
				u = new Vector3D(gridSize, 0, 0);
				v = new Vector3D(0, gridSize, 0);
				normal = new Vector3D(0, 0, orthogonalDirection);
				translation = new Vector3D(halfGrid, halfGrid, orthogonalShift);
			}
			else if (orientation == XZ) {
				u = new Vector3D(gridSize, 0, 0);
				normal = new Vector3D(0, orthogonalDirection, 0);
				v = new Vector3D(0, 0, gridSize);
				translation = new Vector3D(halfGrid, orthogonalShift, halfGrid);
			}
			else if (orientation == YZ) {
				normal = new Vector3D(orthogonalDirection, 0, 0);
				u = new Vector3D(0, gridSize, 0);
				v = new Vector3D(0, 0, gridSize);
				translation = new Vector3D(orthogonalShift, halfGrid, halfGrid);
			}
			else {
				throw new IllegalArgumentException("Unexpected orientation");
			}
		}

		private ValuePair<Vector3D, Vector3D> getSampler() {
			final double uScale = random.nextDouble();
			final double vScale = random.nextDouble();
			final double x = uScale * u.getX() + vScale * v.getX();
			final double y = uScale * u.getY() + vScale * v.getY();
			final double z = uScale * u.getZ() + vScale * v.getZ();
			final Vector3D origin = new Vector3D(x, y, z);
			origin.add(translation);
			return new ValuePair<>(origin, normal);
		}

		private static Vector3D rotate(final Vector3D v, final Quaternion q) {
			final Quaternion p = new Quaternion(0.0, v.toArray());
			final Quaternion rotated = q.multiply(p).multiply(q.getInverse());
			return new Vector3D(rotated.getVectorPart());
		}

		private void rotateVectors(final double[] quaternion) {
			final Quaternion q = new Quaternion(quaternion);
			translation = rotate(translation, q);
			u = rotate(u, q);
			v = rotate(v, q);
			normal = rotate(normal, q);
		}

		public enum Orientation {
				XY, XZ, YZ
		}
	}
	// endregion
}
