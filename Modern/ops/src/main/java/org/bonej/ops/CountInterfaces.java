
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

	/**
	 * Rotates the given point around the axis by theta angle
	 *
	 * @implNote Uses quaternions
	 * @param point A point in 3D space
	 * @param axis  The axis of rotation
	 * @param theta Angle of rotation in radians
	 * @return The rotated point
	 */
	public static Vector3d rotate(final Vector3d point, final Vector3d axis, final double theta) {
		final AxisAngle4d axisAngle4d = new AxisAngle4d(axis, theta);
		final Quat4d q = new Quat4d();
		q.set(axisAngle4d);
		final Quat4d p = new Quat4d();
		p.set(point.x, point.y, point.z, 0.0);
		final Quat4d qInv = new Quat4d();
		qInv.inverse(q);
		final Quat4d rotated = new Quat4d();
		rotated.mul(q, p);
		rotated.mul(qInv);
		return new Vector3d(rotated.x, rotated.y, rotated.z);
	}
}
