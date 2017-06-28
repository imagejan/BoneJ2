
package org.bonej.ops;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.function.Supplier;
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
	private long width;
	private long height;
	private long depth;
	private List<Plane> stackPlanes;

	@Override
	public boolean conforms() {
		// TODO: add support for 2D
		return in1().numDimensions() == 3;
	}

	private Supplier<ValuePair<Plane, Plane>> planePairs = this::selectPlanes;

	@Override
	public Long calculate(final RandomAccessibleInterval<B> interval,
		final Double increment)
	{
		initialize(interval);
		final Stream<Vector3d> samplePoints = Stream.generate(planePairs).parallel().map(
			planes -> {
				final ValuePair<Vector3d, Vector3d> segment = createSegment(planes);
				final Vector3d origin = segment.a;
				final Vector3d direction = findDirection(segment);
				final double tMax = planes.b.intersection(origin, direction);
				final int iterations = (int) Math.floor(Math.abs(tMax / increment));
				return IntStream.rangeClosed(1, iterations).mapToDouble(i -> i *
					increment).mapToObj(t -> createSamplePoint(origin, direction, t))
					.filter(v -> !outOfBounds(v));
			}).flatMap(s -> s).limit(vectors);
		final long[] samples = samplePoints.mapToLong(p -> sample(interval, p))
			.toArray();
		long intersections = 0L;
		for (int i = 0; i < samples.length - 1; i++) {
			intersections += Math.abs(samples[i] - samples[i + 1]);
		}
		return intersections;
	}

	public static Vector3d createSamplePoint(final Vector3d origin,
		final Vector3d direction, final double t)
	{
		final Vector3d samplePoint = new Vector3d(direction);
		samplePoint.scale(t);
		samplePoint.add(origin);
		return samplePoint;
	}

	private ValuePair<Plane, Plane> selectPlanes() {
		final int n = stackPlanes.size();
		final int[] indices = IntStream.range(0, n).toArray();
		final int i = random.nextInt(n);
		swap(indices, i, n - 1);
		final int j = random.nextInt(n - 1);
		final Plane startPlane = stackPlanes.get(indices[i]);
		final Plane endPlane = stackPlanes.get(indices[j]);
		return new ValuePair<>(startPlane, endPlane);
	}

	private void swap(final int[] indices, final int i, final int j) {
		final int tmp = indices[j];
		indices[j] = indices[i];
		indices[i] = tmp;
	}

	private long sample(final RandomAccessibleInterval<B> interval,
		final Vector3d samplePoint)
	{
		final RandomAccess<B> access = interval.randomAccess();
		final long x = Math.round(samplePoint.x);
		final long y = Math.round(samplePoint.y);
		final long z = Math.round(samplePoint.z);
		final long[] coordinates = { x, y, z };
		access.setPosition(coordinates);
		return (long) access.get().getRealDouble();
	}

	private boolean outOfBounds(final Vector3d v) {
		return (v.x < 0) || (v.x > (width - 1)) || (v.y < 0) || (v.y > (height -
			1)) || (v.z < 0) || (v.z > (depth - 1));
	}

	public static Vector3d findDirection(final ValuePair<Vector3d, Vector3d> segment) {
		final Vector3d direction = new Vector3d(segment.b);
		direction.sub(segment.a);
		direction.normalize();
		return direction;
	}

	private ValuePair<Vector3d, Vector3d> createSegment(
		final ValuePair<Plane, Plane> planes)
	{
		final Vector3d startPoint = randomPlanePoint(planes.a);
		final Vector3d endPoint = randomPlanePoint(planes.b);
		return new ValuePair<>(startPoint, endPoint);
	}

	private Vector3d randomPlanePoint(final Plane plane) {
		final double x = random.nextDouble() * (plane.end.x - plane.start.x) +
			plane.start.x;
		final double y = random.nextDouble() * (plane.end.y - plane.start.y) +
			plane.start.y;
		final double z = random.nextDouble() * (plane.end.z - plane.start.z) +
			plane.start.z;
		return new Vector3d(x, y, z);
	}

	private void initialize(final RandomAccessibleInterval<B> interval) {
		width = interval.dimension(0);
		height = interval.dimension(1);
		depth = interval.dimension(2);
		final Plane bottom = new Plane(new Vector3d(width - 1, 0, 0), new Vector3d(
			0, height - 1, 0), new Vector3d(0, 0, 1));
		final Plane top = new Plane(new Vector3d(width - 1, 0, depth - 1),
			new Vector3d(0, height - 1, depth - 1), new Vector3d(0, 0, -1));
		final Plane left = new Plane(new Vector3d(0, height - 1, 0), new Vector3d(0,
			0, depth - 1), new Vector3d(1, 0, 0));
		final Plane right = new Plane(new Vector3d(width - 1, height - 1, 0),
			new Vector3d(width - 1, 0, depth - 1), new Vector3d(-1, 0, 0));
		final Plane front = new Plane(new Vector3d(width - 1, 0, 0), new Vector3d(0,
			0, depth - 1), new Vector3d(0, 1, 0));
		final Plane back = new Plane(new Vector3d(width - 1, height - 1, 0),
			new Vector3d(0, height - 1, depth - 1), new Vector3d(0, -1, 0));
		stackPlanes = Arrays.asList(bottom, top, left, right, front, back);
	}

	private static class Plane {

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
