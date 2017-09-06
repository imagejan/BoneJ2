
package org.bonej.ops;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import net.imagej.ops.Contingent;
import net.imagej.ops.Op;
import net.imagej.ops.special.function.AbstractUnaryFunctionOp;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.FastMath;
import org.scijava.plugin.Plugin;

/**
 * @author Richard Domander
 */
// TODO Contact Kaleb for licensing etc.
// TODO check calculations
@Plugin(type = Op.class)
public class FitEllipsoid extends
	AbstractUnaryFunctionOp<Collection<Vector3D>, FitEllipsoid.Solution>
	implements Contingent
{

	/**
	 * Number of terms in the surface polynomial that needs to be solved.
	 * <p>
	 * Due the number of terms, we also need at least 9 points in the input to
	 * solve the equation.
	 * </p>
	 * 
	 * @see #solveSurface(Collection)
	 */
	public static final int SURFACE_TERMS = 9;

	@Override
	public Solution calculate(final Collection<Vector3D> points) {
		final RealVector polynomial = solveSurface(points);
		final RealMatrix surface = toAlgebraicMatrix(polynomial);
		final RealVector centre = solveCentre(surface);
		final RealMatrix translatedSurface = translateToCentre(surface, centre);
		final EigenDecomposition eigenDecomposition = solveEigenvectors(
			translatedSurface);
		final double[] eigenvalues = eigenDecomposition.getRealEigenvalues();
		final double[] radii = Arrays.stream(eigenvalues).map(
			FitEllipsoid::toRadius).toArray();
		final boolean isEllipsoid = isEllipsoid(eigenvalues);
		return new Solution(surface, centre, eigenDecomposition, radii,
			isEllipsoid);
	}

	@Override
	public boolean conforms() {
		return in().size() >= SURFACE_TERMS;
	}

	// region -- Helper methods --

	/**
	 * Creates a design matrix out of a collection of points.
	 * <p>
	 * The design matrix is used in the least squares fitting of a quadric
	 * surface.
	 * </p>
	 *
	 * @see #solveSurface(Collection)
	 * @param pointCloud points in a 3D space (size >= {@link #SURFACE_TERMS}).
	 * @return a [n][9] matrix of real values.
	 */
	private static RealMatrix createDesignMatrix(
		final Collection<Vector3D> pointCloud)
	{
		final Function<Vector3D, double[]> toRow = (point) -> {
			final double x = point.getX();
			final double y = point.getY();
			final double z = point.getZ();
			// @formatter:off
			return new double[]{
					x * x, y * y, z * z,
					2 * x * y, 2 * x * z, 2 * y * z,
					2 * x, 2 * y, 2 * z
			};
			// @formatter:on
		};
		final double[][] data = pointCloud.stream().map(toRow).toArray(
			double[][]::new);
		return new Array2DRowRealMatrix(data);
	}

	private static RealMatrix createDesignMatrix2(final Collection<Vector3D> points) {
		final double[][] data = points.stream().map(p -> {
			double xx = Math.pow(p.getX(), 2);
			double yy = Math.pow(p.getY(), 2);
			double zz = Math.pow(p.getZ(), 2);
			double xy = 2 * (p.getX() * p.getY());
			double xz = 2 * (p.getX() * p.getZ());
			double yz = 2 * (p.getY() * p.getZ());
			double x = 2 * p.getX();
			double y = 2 * p.getY();
			double z = 2 * p.getZ();
			return new double[]{xx, yy, zz, xy, xz, yz, x, y, z};
		}).toArray(double[][]::new);
		return new Array2DRowRealMatrix(data);
	}

	private static RealVector solveSurface2(RealMatrix d) {
		// Multiply: d' * d
		RealMatrix dtd = d.transpose().multiply(d);

		// Create a vector of ones.
		RealVector ones = new ArrayRealVector(d.getRowDimension());
		ones.mapAddToSelf(1);

		// Multiply: d' * ones.mapAddToSelf(1)
		RealVector dtOnes = d.transpose().operate(ones);

		// Find ( d' * d )^-1
		RealMatrix dtdi = new SingularValueDecomposition(dtd)
				.getSolver().getInverse();

		// v = (( d' * d )^-1) * ( d' * ones.mapAddToSelf(1));

		return dtdi.operate(dtOnes);
	}

	private static boolean isEllipsoid(final double[] eigenvalues) {
		// the signs of the eigenvalues (diagonal elements of DD) determine the
		// type. If they are all positive, it is an ellipsoid; two positive and
		// one negative gives a hyperboloid of one sheet; one positive and two
		// negative gives a hyperboloid of two sheets. If one or more eigenvalues
		// vanish, we have a degenerate case such as a paraboloid,or a cylinder or
		// even a pair of planes.
		return Arrays.stream(eigenvalues).allMatch(x -> x > 0);
	}

	private static RealVector solveCentre(final RealMatrix quadricSurface) {
		final RealMatrix subMatrix = quadricSurface.getSubMatrix(0, 2, 0, 2).scalarMultiply(-1.0);
		final RealMatrix subInverse = MatrixUtils.inverse(subMatrix);
		// The {x,y,z} translation (from origin) part of the matrix
		final RealVector translation = quadricSurface.getRowVector(3).getSubVector(0,
			3);
		return subInverse.operate(translation);
	}

	private static EigenDecomposition solveEigenvectors(
		final RealMatrix quadricSurface)
	{
		final RealMatrix eigenMatrix = quadricSurface.getSubMatrix(0, 2, 0, 2);
		final double scalar = -1.0 / quadricSurface.getEntry(3, 3);
		eigenMatrix.scalarMultiply(scalar);
		return new EigenDecomposition(eigenMatrix);
	}

	/**
	 * Solves the quadric surface that best fits the given points.
	 * <p>
	 * The equation solved is the polynomial Ax<sup>2</sup> + By<sup>2</sup> +
	 * Cz<sup>2</sup> + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy + 2Iz, i.e. the general
	 * equation of a quadric surface.
	 * </p>
	 *
	 * @param pointCloud A collection of points in a 3D space.
	 * @return the solution vector of the surface. The solution may not form an
	 *         ellipsoid.
	 */
	private static RealVector solveSurface(
		final Collection<Vector3D> pointCloud)
	{
		final RealMatrix d = createDesignMatrix(pointCloud);
		final RealMatrix dT = d.transpose();
		// Don't use MatrixUtils.inverse - it may throw a SingularMatrixException
		final RealMatrix dTDInverse = new SingularValueDecomposition(dT.multiply(d)).getSolver().getInverse();
		final ArrayRealVector ones = new ArrayRealVector(pointCloud.size(), 1);
		return dTDInverse.operate(dT.operate(ones));
	}

	/**
	 * Creates a matrix out of a solution vector.
	 *
	 * @see #solveSurface(Collection)
	 * @return a matrix representing the polynomial solution vector in an
	 *         algebraic form.
	 */
	private static RealMatrix toAlgebraicMatrix(final RealVector v) {
		// I'm not a clever man, so I'm using named placeholder variables to better
		// follow the array assignment
		final double a = v.getEntry(0);
		final double b = v.getEntry(1);
		final double c = v.getEntry(2);
		final double d = v.getEntry(3);
		final double e = v.getEntry(4);
		final double f = v.getEntry(5);
		final double g = v.getEntry(6);
		final double h = v.getEntry(7);
		final double i = v.getEntry(8);
		// @formatter:off
		return new Array2DRowRealMatrix(new double[][] {
				{a, d, e, g},
				{d, b, f, h},
				{e, f, c, i},
				{g, h, i, -1}
		});
		// @formatter:on
	}

	private static double toRadius(final double eigenvalue) {
		return FastMath.sqrt(1.0 / eigenvalue);
	}

	private static RealMatrix translateToCentre(final RealMatrix quadricSurface,
		final RealVector centre)
	{
		final double x = centre.getEntry(0);
		final double y = centre.getEntry(1);
		final double z = centre.getEntry(2);
		// @formatter:off
		final Array2DRowRealMatrix t = new Array2DRowRealMatrix(new double[][]{
				{1, 0, 0, x},
				{0, 1, 0, y},
				{0, 0, 1, z},
				{0, 0, 0, 1}
		});
		// @formatter:on
		return t.transpose().multiply(quadricSurface).multiply(t);
	}
	// endregion

	public static void main(String... args) {
		final double a = 1;
		final double b = 2;
		final double c = 3;
		final Vector3D centre = new Vector3D(1, 1, 1);
		final List<Vector3D> ellipsoid = Stream.of(
				new Vector3D(0, 0, 0),
				new Vector3D(1 + a, 1, 1),
				new Vector3D(1 - a, 1, 1),
				new Vector3D(1, 1 + b, 1),
				new Vector3D(1, 1 - b, 1),
				new Vector3D(1, 1, 1 + c),
				new Vector3D(1, 1, 1 - c)
		).map(v -> v.add(centre)).collect(Collectors.toList());
		final RealVector surface = FitEllipsoid.solveSurface(ellipsoid);
		final RealMatrix designMatrix2 = FitEllipsoid.createDesignMatrix2(ellipsoid);
		final RealVector surface2 = FitEllipsoid.solveSurface2(designMatrix2);
		System.out.println();
	}

	// region -- Helper classes --
	public static class Solution {

		private final RealMatrix surface;
		private final RealVector centre;
		private final EigenDecomposition decomposition;
		private final double[] radii;
		private final boolean isEllipsoid;

		private Solution(final RealMatrix surface, final RealVector centre,
			final EigenDecomposition decomposition, final double[] radii,
			final boolean isEllipsoid)
		{
			this.surface = surface;
			this.centre = centre;
			this.decomposition = decomposition;
			this.radii = radii;
			this.isEllipsoid = isEllipsoid;
		}

		/**
		 * Gets the radius of the first axis of the ellipsoid.
		 *
		 * @return a real radius of the ellipsoid.
		 */
		public double getA() {
			return radii[0];
		}

		/**
		 * Gets the radius of the second axis of the ellipsoid.
		 *
		 * @return a real radius of the ellipsoid.
		 */
		public double getB() {
			return radii[1];
		}

		/**
		 * Gets the radius of the third axis of the ellipsoid.
		 *
		 * @return a real radius of the ellipsoid.
		 */
		public double getC() {
			return radii[2];
		}

		/**
		 * Gets a copy of the centre of the ellipsoid.
		 *
		 * @return A copy of the centre vector.
		 */
		public RealVector getCentre() {
			return centre.copy();
		}

		/**
		 * Gets a copy of eigenvalues of the ellipsoid surface.
		 *
		 * @return copy of the real parts of the eigenvalues.
		 */
		public double[] getEigenvalues() {
			return decomposition.getRealEigenvalues();
		}

		/**
		 * Get a copy of the ith eigenvector of the ellipsoid surface.
		 *
		 * @param i Index of the eigenvector.
		 * @return an eigenvector.
		 */
		public RealVector getEigenvector(final int i) {
			return decomposition.getEigenvector(i);
		}

		/**
		 * Gets the algebraic form of the ellipsoid surface.
		 *
		 * @return matrix of the surface terms.
		 */
		public RealMatrix getSurface() {
			return surface.copy();
		}

		/**
		 * Check whether the solved quadric surface is an ellipsoid.
		 *
		 * @return true if the surface is an ellipsoid, false otherwise.
		 */
		public boolean isEllipsoid() {
			return isEllipsoid;
		}
	}
	// endregion
}
