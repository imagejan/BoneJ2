
package org.bonej.ops;

import java.util.Arrays;
import java.util.Collection;
import java.util.Optional;
import java.util.function.Function;

import net.imagej.ops.Contingent;
import net.imagej.ops.Op;
import net.imagej.ops.special.function.AbstractUnaryFunctionOp;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.util.FastMath;
import org.scijava.plugin.Plugin;

/**
 * TODO explanation about least squares distance fitting and why it might not be
 * an ellipsoid
 *
 * @author Richard Domander
 */
// TODO Contact Kaleb for licensing etc.
// TODO check calculations
// TODO return Optional<Solution>
@Plugin(type = Op.class)
public class FitEllipsoid extends
	AbstractUnaryFunctionOp<Collection<Vector3D>, Optional<FitEllipsoid.Solution>>
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
	public Optional<Solution> calculate(final Collection<Vector3D> points) {
		final RealVector polynomial = solveSurface(points);
		final RealMatrix surface = toAlgebraicMatrix(polynomial);
		final RealVector centre = solveCentre(surface);
		final RealMatrix translatedSurface = translateToCentre(surface, centre);
		final EigenDecomposition eigenDecomposition = solveEigenvectors(
			translatedSurface);
		final double[] eigenvalues = eigenDecomposition.getRealEigenvalues();
		final double[] radii = Arrays.stream(eigenvalues).map(
			FitEllipsoid::toRadius).toArray();
        if (isEllipsoid(eigenvalues)) {
            return Optional.of(new Solution(surface, centre, eigenDecomposition, radii));
        } else {
            return Optional.empty();
        }
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

	private static boolean isEllipsoid(final double[] eigenvalues) {
		// the signs of the eigenvalues (diagonal elements of DD) determine the
		// type. If they are all positive, it is an ellipsoid; two positive and
		// one negative gives a hyperboloid of one sheet; one positive and two
		// negative gives a hyperboloid of two sheets. If one or more eigenvalues
		// vanish, we have a degenerate case such as a paraboloid,or a cylinder or
		// even a pair of planes.
		// TODO should work with tolerance x > EPS? Otherwise eigenvalues don't "vanish"
		return Arrays.stream(eigenvalues).allMatch(x -> x > 0);
	}

	private static RealVector solveCentre(final RealMatrix quadricSurface) {
		final RealMatrix subMatrix = quadricSurface.getSubMatrix(0, 2, 0, 2)
			.scalarMultiply(-1.0);
		final RealMatrix subInverse = new SingularValueDecomposition(subMatrix)
			.getSolver().getInverse();
		// The {x,y,z} translation (from origin) part of the matrix
		final RealVector translation = quadricSurface.getRowVector(3).getSubVector(
			0, 3);
		return subInverse.operate(translation);
	}

	private static EigenDecomposition solveEigenvectors(
		final RealMatrix quadricSurface)
	{
		final double scalar = -1.0 / quadricSurface.getEntry(3, 3);
		final RealMatrix eigenMatrix = quadricSurface.getSubMatrix(0, 2, 0, 2)
			.scalarMultiply(scalar);
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
		final RealMatrix dTDInverse = new SingularValueDecomposition(dT.multiply(d))
			.getSolver().getInverse();
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
	// endregion

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

	// region -- Helper classes --
	// TODO Remove isEllipsoid field
	public static class Solution {

		private final RealMatrix surface;
		private final RealVector centre;
		private final EigenDecomposition decomposition;
		private final double[] radii;

		private Solution(final RealMatrix surface, final RealVector centre,
			final EigenDecomposition decomposition, final double[] radii)
		{
			this.surface = surface;
			this.centre = centre;
			this.decomposition = decomposition;
			this.radii = radii;
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
		 * Gets a copy of the radii of the ellipsoid.
		 *
		 * @return an array of radii {a, b, c}
		 */
		public double[] getRadii() {
			return Arrays.copyOf(radii, radii.length);
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
		 * Get a copy of the eigenvectors of the ellipsoid surface.
		 *
		 * @return eigenvector array.
		 */
		public RealVector[] getEigenvectors() {
			// @formatter:off
			return new RealVector[] {
					decomposition.getEigenvector(0),
					decomposition.getEigenvector(1),
					decomposition.getEigenvector(2)
			};
			// @formatter:on
		}

		/**
		 * Gets the algebraic form of the ellipsoid surface.
		 *
		 * @return matrix of the surface terms.
		 */
		public RealMatrix getSurface() {
			return surface.copy();
		}

	}
	// endregion
}
