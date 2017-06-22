
package org.bonej.ops;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

import net.imagej.ops.Contingent;
import net.imagej.ops.Op;
import net.imagej.ops.special.function.AbstractBinaryFunctionOp;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.BooleanType;
import net.imglib2.util.ValuePair;

import org.scijava.plugin.Plugin;
import org.scijava.vecmath.AxisAngle4d;
import org.scijava.vecmath.Quat4d;
import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
@Plugin(type = Op.class)
public class CountInterfaces<B extends BooleanType<B>> extends
	AbstractBinaryFunctionOp<RandomAccessibleInterval<B>, Long, Long> implements
	Contingent
{

	private static final Random random = new Random(System.currentTimeMillis());
	private long width;
	private long height;
	private long depth;
	private long samples = 2;

	@Override
	public boolean conforms() {
		// TODO: add support for 2D
		return in1().numDimensions() == 3;
	}

	@Override
	public Long calculate(final RandomAccessibleInterval<B> interval,
		final Long samples)
	{
		findDimensions(interval);
		List<Vector3d> samplingPlane = createSamplingPlane();
		final double[] angles = { 0, 0, 0 };
		// random.doubles(3, 0.0, Math.PI).toArray();
		final Vector3d planeNormal = rotateXYZ(new Vector3d(0.0, 0.0, 1.0), angles);
		samplingPlane = samplingPlane.stream().map(v -> rotateXYZ(v, angles))
			.collect(Collectors.toList());
		final Vector3d zAdd = new Vector3d(0.0, 0.0, depth * 0.5);
		final List<List<Vector3d>> samplingLines = new ArrayList<>();
		for (final Vector3d point : samplingPlane) {
			point.add(zAdd);
			final ValuePair<Double, Double> parameters =
				findStackIntersectionParameters(point, planeNormal);
			if (parameters == null) {
				continue;
			}
			samplingLines.add(samples(planeNormal, point, parameters));

		}
		return 0L;
	}

	private List<Vector3d> samples(final Vector3d direction,
		final Vector3d origin, final ValuePair<Double, Double> parameters)
	{
		final Vector3d minVector = new Vector3d(direction);
		final Vector3d maxVector = new Vector3d(direction);
		maxVector.scale(parameters.b);
		maxVector.add(origin);
		minVector.scale(parameters.a);
		minVector.add(origin);
		final Vector3d boxSegment = new Vector3d(maxVector);
		boxSegment.sub(minVector);
		final double denominator = samples + 1.0;
		final List<Vector3d> samplePoints = new ArrayList<>((int) samples);
		for (long i = 1; i < denominator; i++) {
			final Vector3d point = new Vector3d(boxSegment);
			point.scale(i / denominator);
			point.add(minVector);
			samplePoints.add(point);
		}
		return samplePoints;
	}

	private ValuePair<Double, Double> findStackIntersectionParameters(
		final Vector3d origin, final Vector3d direction)
	{
		final double minX = 0;
		final double maxX = width;
		final double minY = 0;
		final double maxY = height;
		final double minZ = 0;
		final double maxZ = depth;
		final double tX0 = (minX - origin.x) / direction.x;
		final double tX1 = (maxX - origin.x) / direction.x;
		final double tY0 = (minY - origin.y) / direction.y;
		final double tY1 = (maxY - origin.y) / direction.y;
		final double tZ0 = (minZ - origin.z) / direction.z;
		final double tZ1 = (maxZ - origin.z) / direction.z;
		if (tX0 > tY1 || tY0 > tX1) {
			return null;
		}
		double tMin = Math.max(tX0, tY0);
		double tMax = Math.min(tX1, tY1);
		if (tMin > tZ1 || tZ0 > tMax) {
			return null;
		}
		return new ValuePair<>(Math.max(tZ0, tMin), Math.min(tZ1, tMax));
	}

	private Vector3d rotateXYZ(Vector3d v, final double[] angles) {
		final Vector3d xAxis = new Vector3d(1.0, 0.0, 0.0);
		final Vector3d yAxis = new Vector3d(0.0, 1.0, 0.0);
		final Vector3d zAxis = new Vector3d(0.0, 0.0, 1.0);
		v = rotateAboutAxis(v, xAxis, angles[0]);
		v = rotateAboutAxis(v, yAxis, angles[1]);
		v = rotateAboutAxis(v, zAxis, angles[2]);
		return v;
	}

	// Create a sampling plane centred around (0, 0, 0)
	private List<Vector3d> createSamplingPlane() {
		final int capacity = (int) (samples * samples);
		final List<Vector3d> planePoints = new ArrayList<>(capacity);
		final double denominator = samples + 1.0;
		final Vector3d p0 = new Vector3d(0.0, 0.0, 0.0);
		final Vector3d p1 = new Vector3d(width, height, 0.0);
		for (int b = 1; b <= samples; b++) {
			for (int a = 1; a <= samples; a++) {
				final Vector3d p = splitSegment(p0, p1, a, b, 0, denominator);
				planePoints.add(p);
			}
		}
		return planePoints;
	}

	private void findDimensions(final RandomAccessibleInterval<B> interval) {
		width = interval.dimension(0);
		height = interval.dimension(1);
		depth = interval.dimension(2);
	}

	/**
	 * Rotates the given point around the axis by theta angle
	 *
	 * @implNote Uses quaternions
	 * @param point A point in 3D space
	 * @param axis The axis of rotation
	 * @param theta Angle of rotation in radians
	 * @return The rotated point
	 */
	public static Vector3d rotateAboutAxis(final Vector3d point,
		final Vector3d axis, final double theta)
	{
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

	public Vector3d splitSegment(final Vector3d p0, final Vector3d p1,
		final double xNumerator, final double yNumerator, final double zNumerator,
		final double denominator)
	{
		final Vector3d split = new Vector3d();
		split.sub(p1, p0);
		split.setX(split.x * (xNumerator / denominator));
		split.setY(split.y * (yNumerator / denominator));
		split.setZ(split.z * (zNumerator / denominator));
		split.add(p0);
		return split;
	}
}
