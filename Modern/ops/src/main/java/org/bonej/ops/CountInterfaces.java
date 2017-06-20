
package org.bonej.ops;

import java.util.Random;

import net.imagej.ops.Op;
import net.imagej.ops.special.function.AbstractBinaryFunctionOp;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.BooleanType;

import org.scijava.plugin.Plugin;
import org.scijava.vecmath.AxisAngle4d;
import org.scijava.vecmath.Quat4d;
import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
@Plugin(type = Op.class)
public class CountInterfaces<B extends BooleanType<B>> extends
	AbstractBinaryFunctionOp<RandomAccessibleInterval<B>, Long, Long>
{

	private static final Random random = new Random(System.currentTimeMillis());

	@Override
	public Long calculate(final RandomAccessibleInterval<B> interval,
		final Long samples)
	{
		final Vector3d centrePoint = findCentrePoint(interval);
		final double angle = random.nextDouble() * Math.PI;
		rotate(new Vector3d(1.0, 0.0, 0.0), centrePoint, angle);
		return 0L;
	}

	private static <B extends BooleanType<B>> Vector3d findCentrePoint(
		final RandomAccessibleInterval<B> interval)
	{
		final int n = interval.numDimensions();
		final long[] mins = new long[n];
		interval.min(mins);
		final long[] maxs = new long[n];
		interval.max(maxs);
		final double[] coordinates = new double[3];
		final int copy = n > 2 ? 3 : 2;
		for (int i = 0; i < copy; i++) {
			coordinates[i] = maxs[i] - mins[i];
		}
		return new Vector3d(coordinates);
	}

	public static Vector3d rotate(Vector3d v, Vector3d centrePoint,
		double angle)
	{
		final AxisAngle4d axisAngle4d = new AxisAngle4d(centrePoint, angle);
		final Quat4d q = new Quat4d();
		q.set(axisAngle4d);
		final Quat4d p = new Quat4d(v.x, v.y, v.z, 0.0);
		final Quat4d rotated = new Quat4d();
		q.mul(p, rotated);
		rotated.mulInverse(q);
		return new Vector3d(rotated.x, rotated.y, rotated.z);
	}

	public static void main(String... args) {
		final Vector3d v = new Vector3d(1, 0, 0);
		final Vector3d centre = new Vector3d(0, 0, 1);
		final Vector3d u = CountInterfaces.rotate(v, centre, Math.PI / 2.0);
		System.out.println(u.toString());
	}
}
