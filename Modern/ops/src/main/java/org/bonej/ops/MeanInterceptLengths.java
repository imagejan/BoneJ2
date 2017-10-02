
package org.bonej.ops;

import static org.bonej.ops.MeanInterceptLengths.SamplingPlane.Orientation.XY;
import static org.bonej.ops.MeanInterceptLengths.SamplingPlane.Orientation.XZ;
import static org.bonej.ops.MeanInterceptLengths.SamplingPlane.Orientation.YZ;

import java.util.Arrays;
import java.util.Collection;
import java.util.Objects;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.LongStream;
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
	private Long nVectors = null;

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
		final long[] dimensions = findDimensions(interval);
		fillParameters(dimensions);
		final SamplingGrid grid = new SamplingGrid(interval, rotation);
		final RandomAccess<B> access = interval.randomAccess();
		final Stream<SamplingSection> sections = grid.getSamplingSections(nVectors)
			.filter(SamplingSection::intersectsStack);
		return sections.map(section -> createMILVector(access, section, tIncrement))
			.filter(Objects::nonNull).collect(Collectors.toList());
	}

	// TODO Check RAI dimensionality
	@Override
	public boolean conforms() {
		return rotation == null || rotation.length == 4;
	}

	public static <B extends BooleanType<B>> Vector3D createMILVector(
		final RandomAccess<B> access, final SamplingSection section,
		final double increment)
	{
		final long intercepts = sampleIntercepts(access, section, increment);
		if (intercepts < 0) {
			return null;
		}
		return toMILVector(section, intercepts);
	}

	// region -- Helper methods (public static for unit testing) --
	public static <C extends NumericType<C>> long[] findDimensions(
		final RandomAccessibleInterval<C> interval)
	{
		final long[] bounds = new long[interval.numDimensions()];
		interval.dimensions(bounds);
		return bounds;
	}

	public static <B extends BooleanType<B>> long getElement(
		final RandomAccess<B> access, final Vector3D position)
	{
		access.setPosition((int) position.getX(), 0);
		access.setPosition((int) position.getY(), 1);
		access.setPosition((int) position.getZ(), 2);
		return (long) access.get().getRealDouble();
	}

	public static <B extends BooleanType<B>> long sampleIntercepts(
		final RandomAccess<B> access, final SamplingSection section,
		final double increment)
	{
		final double offset = random.nextDouble() * increment;
		final Vector3D samplingGap = section.direction.scalarMultiply(increment);
		final double tMin = section.tMin;
		final double tMax = section.tMax - offset;
		if (tMin > tMax) {
			return -1;
		}
		Vector3D point = section.direction.scalarMultiply(tMin + offset).add(
			section.origin);
		long intercepts = 0;
		long previous = getElement(access, point);
		for (double t = tMin; t <= tMax; t += increment) {
			final long current = getElement(access, point);
			intercepts += current ^ previous;
			point = point.add(samplingGap);
			previous = current;
		}
		return intercepts;
	}

	/**
	 * Creates a mean interception length (MIL) vector.
	 * <p>
	 * The MIL vector is parallel to the direction of the sampling vector. Its
	 * length equals the section of the sampling vector withing the stack
	 * (interval), and its origin is (0, 0, 0).
	 * </p>
	 *
	 * @param section
	 * @param intercepts number of foreground - background interceptions along
	 *          sampling vector.
	 * @return a new MIL vector.
	 */
	public static Vector3D toMILVector(final SamplingSection section,
		final long intercepts)
	{
		final Vector3D v = section.direction.scalarMultiply(section.length);
		if (intercepts < 2) {
			return v;
		}
		return v.scalarMultiply(1.0 / intercepts);
	}

	private void fillParameters(final long[] bounds) {
		if (nVectors == null) {
			final long n = Arrays.stream(bounds).max().orElse(0);
			nVectors = n * n;
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

	public static final class SamplingSection {

		private Vector3D origin;
		private Vector3D direction;
		private double tMin;
		private double tMax;
		private boolean intersects;
		private double length;

		public SamplingSection(final Vector3D origin, final Vector3D direction,
			final RandomAccessibleInterval<?> interval)
		{
			this.origin = new Vector3D(1.0, origin);
			this.direction = new Vector3D(1.0, direction);
			final long[] min = new long[3];
			final long[] max = new long[3];
			interval.min(min);
			interval.max(max);
			findIntervalIntersections(this.origin, this.direction, min, max);
			length = new Vector3D(tMax, direction).subtract(tMin, direction)
				.getNorm();
		}

		public boolean intersectsStack() {
			return intersects;
		}

		private void findIntervalIntersections(final Vector3D origin,
			final Vector3D direction, final long[] min, final long[] max)
		{
			intersects = false;
			// Because max = {w - 1, h - 1, d - 1} we need to add eps to the bounds to
			// get the correct intersection point. However, it can't equal 1, because
			// otherwise vector coordinates will floor to w, h or d, which will cause
			// an index out of bounds exception.
			final double eps = 1 - 1e-12;
			final double minX = direction.getX() >= 0.0 ? min[0] : max[0] + eps;
			final double maxX = direction.getX() >= 0.0 ? max[0] + eps : min[0];
			final double minY = direction.getY() >= 0.0 ? min[1] : max[1] + eps;
			final double maxY = direction.getY() >= 0.0 ? max[1] + eps : min[1];
			final double minZ = direction.getZ() >= 0.0 ? min[2] : max[2] + eps;
			final double maxZ = direction.getZ() >= 0.0 ? max[2] + eps : min[2];
			final double tX0 = (minX - origin.getX()) / direction.getX();
			final double tX1 = (maxX - origin.getX()) / direction.getX();
			final double tY0 = (minY - origin.getY()) / direction.getY();
			final double tY1 = (maxY - origin.getY()) / direction.getY();
			final double tZ0 = (minZ - origin.getZ()) / direction.getZ();
			final double tZ1 = (maxZ - origin.getZ()) / direction.getZ();
			if (tX0 > tY1 || tY0 > tX1) {
				return;
			}
			tMin = maxNan(tX0, tY0);
			tMax = minNan(tX1, tY1);
			if (tMin > tZ1 || tZ0 > tMax) {
				return;
			}
			tMin = maxNan(tZ0, tMin);
			tMax = minNan(tZ1, tMax);
			if (Double.isNaN(tMin) || Double.isNaN(tMax)) {
				return;
			}
			if (tMin > tMax) {
				final double tmp = tMin;
				tMin = tMax;
				tMax = tmp;
			}
			intersects = true;
		}

		private static double maxNan(final double a, final double b) {
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

		private static double minNan(final double a, final double b) {
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
	}

	public static final class SamplingGrid {

		private final SamplingPlane xy;
		private final SamplingPlane xz;
		private final SamplingPlane yz;
		private final RandomAccessibleInterval<?> interval;

		public SamplingGrid(final RandomAccessibleInterval<?> interval,
			final double[] quaternion)
		{
			this.interval = interval;
			final double gridSize = findGridSize(interval);
			final Vector3D centroid = findCentroid(interval);
			xy = new SamplingPlane(gridSize, XY, quaternion, centroid);
			xz = new SamplingPlane(gridSize, XZ, quaternion, centroid);
			yz = new SamplingPlane(gridSize, YZ, quaternion, centroid);
		}

		private static Vector3D findCentroid(
			final RandomAccessibleInterval<?> interval)
		{
			final double[] coordinates = IntStream.range(0, 3).mapToDouble(
				i -> (interval.max(i) + interval.min(i)) * 0.5).toArray();
			return new Vector3D(coordinates);
		}

		private static double findGridSize(
			final RandomAccessibleInterval<?> interval)
		{
			final long squaredSum = LongStream.of(interval.dimension(0), interval
				.dimension(1), interval.dimension(2)).map(d -> d * d).sum();
			return Math.sqrt(squaredSum);
		}

		private Stream<SamplingSection> getSamplingSections(final long n) {
			final Stream.Builder<ValuePair<Vector3D, Vector3D>> builder = Stream
				.builder();
			for (long i = 0; i < n; i++) {
				builder.accept(xy.getSamplingVector());
				builder.accept(xz.getSamplingVector());
				builder.accept(yz.getSamplingVector());
			}
			return builder.build().map(sampler -> new SamplingSection(sampler.a,
				sampler.b, interval));
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

		private ValuePair<Vector3D, Vector3D> getSamplingVector() {
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
			final Quaternion q = new Quaternion(quaternion[3], quaternion[0],
				quaternion[1], quaternion[2]);
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
