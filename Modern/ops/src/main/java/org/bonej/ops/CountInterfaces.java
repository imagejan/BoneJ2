
package org.bonej.ops;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

import net.imagej.ops.Contingent;
import net.imagej.ops.Op;
import net.imagej.ops.special.function.AbstractBinaryFunctionOp;
import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.BooleanType;
import net.imglib2.util.ValuePair;
import net.imglib2.view.ExtendedRandomAccessibleInterval;
import net.imglib2.view.Views;

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

	private static final Random random = new Random(System.currentTimeMillis());
	private long width;
	private long height;
	private long depth;
	private List<Plane> planes;

	@Override
	public boolean conforms() {
		// TODO: add support for 2D
		return in1().numDimensions() == 3;
	}

	@Override
	public Long calculate(final RandomAccessibleInterval<B> interval,
		final Double increment)
	{
		initialize(interval);
		final ValuePair<Plane, Plane> planes = selectPlanes();
		final ValuePair<Vector3d, Vector3d> segment = createSegment(planes);
		final Vector3d origin = segment.a;
		final Vector3d direction = findDirection(segment);
		final double tMax = planes.b.intersection(origin, direction);
		if (tMax == 0.0) {
			return 0L;
		}

		final long iterations = (long) Math.floor(tMax / increment);
		double t = increment;
		long intersections = 0;
		boolean lastPoint = sample(interval, origin, direction, t);
		for (int i = 2; i <= iterations; i++) {
			// Optional - not present when sample OoB
			final boolean point = sample(interval, origin, direction, t);
			if (lastPoint != point) {
				intersections++;
			}
			t += increment;
			lastPoint = point;

		}

		return intersections;
	}

	private ValuePair<Plane, Plane> selectPlanes() {
		int n = planes.size();
		final int[] indices = IntStream.range(0, n).toArray();
		int i = random.nextInt(n);
		final Plane startPlane = planes.get(indices[i]);
		final int tmp = indices[n - 1];
		indices[n - 1] = indices[i];
		indices[i] = tmp;
		n = n - 1;
		i = random.nextInt(n);
		final Plane endPlane = planes.get(indices[i]);
		return new ValuePair<>(startPlane, endPlane);
	}

	private boolean sample(final RandomAccessibleInterval<B> interval,
		final Vector3d origin, final Vector3d direction, final double t)
	{
		final Vector3d samplePoint = new Vector3d(direction);
		samplePoint.scale(t);
		samplePoint.add(origin);
		RandomAccess<B> access = interval.randomAccess();
		final B variable = access.get().createVariable();
		variable.set(false);
		final ExtendedRandomAccessibleInterval<B, RandomAccessibleInterval<B>> safeInterval =
			Views.extendValue(interval, variable);
		final long[] coordinates = { Math.round(samplePoint.x), Math.round(
			samplePoint.y), Math.round(samplePoint.z) };
		final String s = Arrays.toString(coordinates);
		System.out.println(s);
		access = safeInterval.randomAccess();
		access.setPosition(coordinates);
		return access.get().get();
	}

	private Vector3d findDirection(final ValuePair<Vector3d, Vector3d> segment) {
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
		final double startX = Math.min(plane.u.x, plane.v.x);
		final double endX = Math.max(plane.u.x, plane.v.x);
		final double startY = Math.min(plane.u.y, plane.v.y);
		final double endY = Math.max(plane.u.y, plane.v.y);
		final double startZ = Math.min(plane.u.z, plane.v.z);
		final double endZ = Math.max(plane.u.z, plane.v.z);
		final double x = random.nextDouble() * (endX - startX) + startX;
		final double y = random.nextDouble() * (endY - startY) + startY;
		final double z = random.nextDouble() * (endZ - startZ) + startZ;
		return new Vector3d(x, y, z);
	}

	private void initialize(final RandomAccessibleInterval<B> interval) {
		width = interval.dimension(0);
		height = interval.dimension(1);
		depth = interval.dimension(2);
		final Vector3d o = new Vector3d();
		final Plane bottom = new Plane(o, new Vector3d(0.0, height, 0.0), o,
			new Vector3d(width, 0.0, 0.0));
		final Plane top = new Plane(new Vector3d(0.0, 0.0, depth), new Vector3d(0.0,
			height, depth), new Vector3d(0.0, 0.0, depth), new Vector3d(width, 0.0,
				depth));
		final Plane left = new Plane(o, new Vector3d(0.0, height, 0.0), o,
			new Vector3d(0.0, 0.0, depth));
		final Plane right = new Plane(new Vector3d(width, 0.0, 0.0), new Vector3d(
			width, height, 0.0), new Vector3d(width, 0.0, 0.0), new Vector3d(width,
				0.0, depth));
		final Plane front = new Plane(o, new Vector3d(width, 0.0, 0.0), o,
			new Vector3d(0.0, 0.0, depth));
		final Plane back = new Plane(new Vector3d(0.0, height, 0.0), new Vector3d(
			width, height, 0.0), new Vector3d(0.0, height, 0.0), new Vector3d(0.0,
				height, depth));
		planes = Arrays.asList(bottom, top, left, right, front, back);
	}

	private static class Plane {

		public final Vector3d vMin;
		public final Vector3d vMax;
		public final Vector3d uMin;
		public final Vector3d u;
		public final Vector3d uMax;
		public final Vector3d v;
		public final Vector3d n;
		public final Vector3d p;

		public Plane(final Vector3d vMin, final Vector3d vMax, final Vector3d uMin,
			final Vector3d uMax)
		{
			this.vMin = vMin;
			this.vMax = vMax;
			this.uMin = uMin;
			this.uMax = uMax;
			u = new Vector3d(uMax);
			u.sub(uMin);
			v = new Vector3d(vMax);
			v.sub(vMin);
			p = new Vector3d((u.x + v.x) * 0.5, (u.y + v.y) * 0.5, (u.z + v.z) * 0.5);
			n = new Vector3d();
			n.cross(u, v);
			n.normalize();
		}

		public double intersection(Vector3d origin, Vector3d direction) {
			final Vector3d v = new Vector3d(p);
			v.sub(origin);
			return v.dot(n) / direction.dot(n);
		}
	}
}
