
package org.bonej.ops;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.stream.DoubleStream;
import java.util.stream.LongStream;
import java.util.stream.Stream;

import net.imglib2.RandomAccess;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.ComplexType;
import net.imglib2.util.ValuePair;

import org.scijava.vecmath.AxisAngle4d;
import org.scijava.vecmath.Quat4d;
import org.scijava.vecmath.Vector3d;

/**
 * @author Richard Domander
 */
public class CountInterfacesGrid {

	private static Random random = new Random(0xc0ff33);
	private static final Vector3d xAxis = new Vector3d(1.0, 0.0, 0.0);
	private static final Vector3d yAxis = new Vector3d(0.0, 1.0, 0.0);
	private static final Vector3d zAxis = new Vector3d(0.0, 0.0, 1.0);

	public static <C extends ComplexType<C>> long[] findBounds(
		final RandomAccessibleInterval<C> interval)
	{
		final long[] bounds = new long[interval.numDimensions()];
		interval.dimensions(bounds);
		return bounds;
	}

	public static Stream<ValuePair<Vector3d, Vector3d>> plotSamplers(
		final long gridSize, final long sections, final long[] bounds)
	{
		final Stream.Builder<ValuePair<Vector3d, Vector3d>> samplingBuilder = Stream
			.builder();
		plotPlane(samplingBuilder, bounds, gridSize, sections, 0, 1);
		plotPlane(samplingBuilder, bounds, gridSize, sections, 0, 2);
		plotPlane(samplingBuilder, bounds, gridSize, sections, 1, 2);
		final double[] angles = random.doubles(3, 0, Math.PI).toArray();
		final Vector3d centroid = new Vector3d(bounds[0] - 1, bounds[1] - 1,
			bounds[2] - 1);
		centroid.scale(0.5);
		return samplingBuilder.build().map(s -> rotateSampler(s, angles)).map(s -> {
			s.a.add(centroid);
			return s;
		});
	}

	public static ValuePair<Vector3d, Vector3d> rotateSampler(
		final ValuePair<Vector3d, Vector3d> sampler, final double[] angles)
	{
		return new ValuePair<>(rotateXYZ(sampler.a, angles), rotateXYZ(sampler.b,
			angles));
	}

	private static Vector3d rotateXYZ(Vector3d v, final double[] angles) {
		v = rotateAboutAxis(v, xAxis, angles[0]);
		v = rotateAboutAxis(v, yAxis, angles[1]);
		v = rotateAboutAxis(v, zAxis, angles[2]);
		return v;
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

	private static Vector3d samplePoint(
		final ValuePair<Vector3d, Vector3d> sampler, final double t)
	{
		final Vector3d point = new Vector3d(sampler.b);
		point.scale(t);
		point.add(sampler.a);
		return point;
	}

	public static void plotPlane(
		final Stream.Builder<ValuePair<Vector3d, Vector3d>> samplingBuilder,
		final long[] bounds, final long gridSize, final long segments,
		final int dim0, final int dim1)
	{
		final Set<Integer> dims = new HashSet<>(Arrays.asList(0, 1, 2));
		dims.remove(dim0);
		dims.remove(dim1);
		final int[] orthogonalDims = { dims.iterator().next() };
		final int[] planeDims = { dim0, dim1 };
		final Vector3d normal = createSparseVector(orthogonalDims, 1.0);
		final Vector3d planeShift = createSparseVector(orthogonalDims, -0.5 * gridSize);
		final Vector3d p0 = createSparseVector(planeDims, -0.5 * gridSize, -0.5 *
			gridSize);
		final Vector3d p1 = createSparseVector(planeDims, 0.5 * gridSize, 0.5 *
			gridSize);
		final Vector3d gridLine = new Vector3d();
		gridLine.sub(p1, p0);
		final double[] coords = new double[3];
		gridLine.get(coords);
		// TODO Add random offset
		for (long i = 0; i <= segments; i++) {
			for (long j = 0; j <= segments; j++) {
				final Vector3d v = createSparseVector(planeDims, (1.0 * j / segments *
					coords[dim0]), (1.0 * i / segments * coords[dim1]));
				v.add(planeShift);
				v.add(p0);
				samplingBuilder.add(new ValuePair<>(v, normal));
			}
		}
	}

	public static Vector3d createSparseVector(int[] dims, double... values) {
		final double[] coordinates = new double[3];
		for (int i = 0; i < dims.length; i++) {
			coordinates[dims[i]] = values[i];
		}
		return new Vector3d(coordinates);
	}

	public static boolean outOfBounds(final long[] coordinates,
		final long[] bounds)
	{
		return (coordinates[0] < 0) || (coordinates[0] >= (bounds[0])) ||
			(coordinates[1] < 0) || (coordinates[1] >= (bounds[1])) ||
			(coordinates[2] < 0) || (coordinates[2] >= (bounds[2]));
	}

	public static long findGridSize(final long[] bounds) {
		final long sumSquared = Arrays.stream(bounds).map(i -> i * i).sum();
		return (long) Math.ceil(Math.sqrt(sumSquared));
	}

	public static long[] toVoxelCoordinates(final Vector3d v) {
		final double[] coordinates = new double[3];
		v.get(coordinates);
		final long[] voxelCoordinates = new long[3];
		for (int i = 0; i < 3; i++) {
			voxelCoordinates[i] = (long) Math.floor(coordinates[i]);
		}
		return voxelCoordinates;
	}

	public static <C extends ComplexType<C>> void sample(
		final RandomAccessibleInterval<C> interval, final long[] coordinates)
	{
		final RandomAccess<C> access = interval.randomAccess(interval);
		access.setPosition(coordinates);
		final C voxel = access.get();
		voxel.setReal(voxel.getRealDouble() + 1.0);
	}

	public static Stream<Vector3d> samplePoints(
		final Stream<ValuePair<Vector3d, Vector3d>> samplers,
		final double increment, final long gridSize)
	{
		final long iterations = (long) Math.floor(gridSize / increment) + 1;

		return samplers.flatMap(sampler -> DoubleStream.iterate(0.0, t -> t +
			increment).mapToObj(t -> CountInterfacesGrid.samplePoint(sampler, t)).limit(iterations));
	}
}
