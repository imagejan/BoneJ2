
package org.bonej.ops;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.function.Function;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import java.util.stream.Stream;

import net.imagej.ops.Contingent;
import net.imagej.ops.Op;
import net.imagej.ops.special.function.AbstractBinaryFunctionOp;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.BooleanType;
import net.imglib2.util.ValuePair;

import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
@Plugin(type = Op.class)
public class CountInterfaces<B extends BooleanType<B>> extends
	AbstractBinaryFunctionOp<RandomAccessibleInterval<B>, Double, Long> implements
	Contingent
{

	@Parameter(required = false)
	private long vectors = 50_000;

	private static final Random random = new Random(System.currentTimeMillis());
	private List<Plane> stackPlanes;

	@Override
	public boolean conforms() {
		// TODO: add support for 2D
		return in1().numDimensions() == 3;
	}

	/**
	 * Allow setting the seed of the random generator to facilitate reproducible
	 * tests
	 */
	public static void setSeed(final long seed) {
		random.setSeed(seed);
	}

	@Override
	public Long calculate(final RandomAccessibleInterval<B> interval,
		final Double increment)
	{
		stackPlanes = createStackPlanes(interval);
		final Stream<ValuePair<Plane, Plane>> pairs = Stream.generate(
			() -> selectPlanes(stackPlanes));
		final long[] bounds = findBounds(interval);
		final Stream<Vector3d> samplePoints = createSamplePoints(pairs, increment,
			bounds);
		final double[] samples = samplePoints.mapToDouble(p -> sample(interval, p))
			.toArray();
		long intersections = 0L;
		for (int i = 0; i < samples.length - 1; i++) {
			intersections += Math.abs(samples[i] - samples[i + 1]);
		}
		return intersections;
	}

	public static long[] findBounds(final RandomAccessibleInterval interval) {
		final long[] bounds = new long[interval.numDimensions()];
		interval.dimensions(bounds);
		return bounds;
	}

	public static Stream<Vector3d> createSamplePoints(
		final Stream<ValuePair<Plane, Plane>> pairs, final double increment,
		final long[] bounds)
	{
		return pairs.parallel().map(planes -> {
			final ValuePair<Vector3d, Vector3d> segment = createSegment(planes);
			final Vector3d origin = segment.a;
			final Vector3d direction = findDirection(segment);
			final double tMax = planes.b.intersection(origin, direction);
			final int iterations = (int) Math.floor(Math.abs(tMax / increment));
			//TODO multiply by sign of tMax
			return IntStream.rangeClosed(1, iterations).mapToDouble(i -> i *
				increment).mapToObj(t -> createSamplePoint(origin, direction, t));
		}).flatMap(s -> s).filter(p -> !outOfBounds(p, bounds));
	}

	public static Vector3d createSamplePoint(final Vector3d origin,
		final Vector3d direction, final double t)
	{
		final Vector3d samplePoint = new Vector3d(direction);
		samplePoint.scale(t);
		samplePoint.add(origin);
		return samplePoint;
	}

	public static ValuePair<Plane, Plane> selectPlanes(List<Plane> planes) {
		final int n = planes.size();
		final int[] indices = IntStream.range(0, n).toArray();
		final int i = random.nextInt(n);
		swap(indices, i, n - 1);
		final int j = random.nextInt(n - 1);
		final Plane startPlane = planes.get(indices[i]);
		final Plane endPlane = planes.get(indices[j]);
		return new ValuePair<>(startPlane, endPlane);
	}

	private static void swap(final int[] indices, final int i, final int j) {
		final int tmp = indices[j];
		indices[j] = indices[i];
		indices[i] = tmp;
	}

	public static long[] toVoxelCoordinates(final Vector3d v) {
		final Function<Double, Long> rounder = d -> (long) (random.nextDouble() > 0.5  ? Math.ceil(d) : Math.floor(d));
		final long x = rounder.apply(v.x);
		final long y = rounder.apply(v.y);
		final long z = rounder.apply(v.z);
		return new long[] { x, y, z };
	}

	public static double[] getOrientation(Vector3d v) {
        final Vector3d xAxis = new Vector3d(1, 0, 0);
        final Vector3d yAxis = new Vector3d(0, 1, 0);
        final Vector3d zAxis = new Vector3d(0, 0, 1);
        return new double[]{v.dot(xAxis), v.dot(yAxis), v.dot(zAxis)};
    }

	public static <B extends BooleanType<B>> double sample(
		final RandomAccessibleInterval<B> interval, final Vector3d samplePoint)
	{
		final long[] coordinates = toVoxelCoordinates(samplePoint);
		final RandomAccess<B> access = interval.randomAccess();
		access.setPosition(coordinates);
		return access.get().getRealDouble();
	}

	public static boolean outOfBounds(final Vector3d v, final long[] bounds) {
		return (v.x < 0) || (v.x > (bounds[0] - 1)) || (v.y < 0) ||
			(v.y > (bounds[1] - 1)) || (v.z < 0) || (v.z > (bounds[2] - 1));
	}

	public static Vector3d findDirection(
		final ValuePair<Vector3d, Vector3d> segment)
	{
		final Vector3d direction = new Vector3d(segment.b);
		direction.sub(segment.a);
		direction.normalize();
		return direction;
	}

	private static ValuePair<Vector3d, Vector3d> createSegment(
		final ValuePair<Plane, Plane> planes)
	{
		final Vector3d startPoint = randomPoint(planes.a.start, planes.a.end);
		final Vector3d endPoint = randomPoint(planes.b.start, planes.b.end);
		return new ValuePair<>(startPoint, endPoint);
	}

	/** Returns a random point in the space defined by the two vectors */
	private static Vector3d randomPoint(final Vector3d start,
		final Vector3d end)
	{
		final double x = random.nextDouble() * (end.x - start.x) + start.x;
		final double y = random.nextDouble() * (end.y - start.y) + start.y;
		final double z = random.nextDouble() * (end.z - start.z) + start.z;
		return new Vector3d(x, y, z);
	}

	public static List<Plane> createStackPlanes(
		final RandomAccessibleInterval stack)
	{
		final long xMax = stack.dimension(0) - 1;
		final long yMax = stack.dimension(1) - 1;
		final long zMax = stack.dimension(2) - 1;
		final Plane bottom = new Plane(new Vector3d(xMax, 0, 0), new Vector3d(0,
			yMax, 0), new Vector3d(0, 0, 1));
		if (stack.numDimensions() == 2) {
			return Collections.singletonList(bottom);
		}
		final Plane top = new Plane(new Vector3d(xMax, 0, zMax), new Vector3d(0,
			yMax, zMax), new Vector3d(0, 0, -1));
		final Plane left = new Plane(new Vector3d(0, yMax, 0), new Vector3d(0, 0,
			zMax), new Vector3d(1, 0, 0));
		final Plane right = new Plane(new Vector3d(xMax, yMax, 0), new Vector3d(
			xMax, 0, zMax), new Vector3d(-1, 0, 0));
		final Plane front = new Plane(new Vector3d(xMax, 0, 0), new Vector3d(0, 0,
			zMax), new Vector3d(0, 1, 0));
		final Plane back = new Plane(new Vector3d(xMax, yMax, 0), new Vector3d(0,
			yMax, zMax), new Vector3d(0, -1, 0));
		return Arrays.asList(bottom, top, left, right, front, back);
	}

	public static class Plane {

		public final Vector3d u;
		public final Vector3d v;
		public final Vector3d n;
		public final Vector3d p;
		public final Vector3d start;
		public final Vector3d end;

		public Plane(final Vector3d u, final Vector3d v, final Vector3d n) {
			this.u = new Vector3d(u);
			this.v = new Vector3d(v);
			p = new Vector3d((u.x + v.x) * 0.5, (u.y + v.y) * 0.5, (u.z + v.z) * 0.5);
			this.n = new Vector3d(n);
			start = new Vector3d();
			end = new Vector3d();
			findBounds();
		}

		private void findBounds() {
			final double[] uCoords = new double[3];
			final double[] vCoords = new double[3];
			u.get(uCoords);
			v.get(vCoords);
			final double[] minCoords = new double[3];
			final double[] maxCoords = new double[3];
			for (int i = 0; i < 3; i++) {
				if (uCoords[i] <= vCoords[i]) {
					minCoords[i] = uCoords[i];
					maxCoords[i] = vCoords[i];
				}
				else {
					minCoords[i] = vCoords[i];
					maxCoords[i] = uCoords[i];
				}
			}
			start.set(minCoords);
			end.set(maxCoords);
		}

		public double intersection(Vector3d origin, Vector3d direction) {
			final Vector3d v = new Vector3d(p);
			v.sub(origin);
			return v.dot(n) / direction.dot(n);
		}
	}
}
