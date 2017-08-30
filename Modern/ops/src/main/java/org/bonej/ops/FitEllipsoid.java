
package org.bonej.ops;

import java.util.Collection;

import net.imagej.ops.Contingent;
import net.imagej.ops.special.function.AbstractUnaryFunctionOp;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

/**
 * @author Richard Domander
 */
// TODO Contact Kaleb for licensing etc.
// TODO check calculations
public class FitEllipsoid extends
	AbstractUnaryFunctionOp<Collection<Vector3D>, FitEllipsoid.Solution>
	implements Contingent
{

	private static final int SURFACE_PARAMETERS = 9;

	@Override
	public Solution calculate(final Collection<Vector3D> points) {
		final RealMatrix matrix = solveSurfacePolynomial(points);
		return new Solution(matrix);
	}

	@Override
	public boolean conforms() {
		// Need at least nine data points to do a fitting
		return in().size() >= SURFACE_PARAMETERS;
	}

	/**
	 * Solves the algebraic form matrix of the polynomial that describes the
	 * ellipsoid that best fits the point cloud
	 *
	 * @implNote Does not guarantee that the quadric quadricSurface solved is an
	 *           ellipsoid
	 * @return The algebraic form of the ellipsoid
	 */
	private RealMatrix solveSurfacePolynomial(
		final Collection<Vector3D> pointCloud)
	{
		// Solve the polynomial Ax^2 + By^2 + Cz^2 + 2Dxy + 2Exz + 2Fyz + 2Gx + 2Hy
		// + 2Iz, i.e. the general equation of a quadric quadricSurface
		final double[][] data = pointCloud.stream().map(this::toDesignRow).toArray(
			double[][]::new);
		final RealMatrix d = new Array2DRowRealMatrix(data);
		final RealMatrix dT = d.transpose();
		final RealMatrix dTDInverse = MatrixUtils.inverse(dT.multiply(d));
		final ArrayRealVector ones = new ArrayRealVector(pointCloud.size(), 1);
		final RealVector v = dTDInverse.operate(dT.operate(ones));
		// TODO return solved vector and fix javadoc
		// Convert the solution vector to matrix
		return toAlgebraicMatrix(v);
	}

	private RealMatrix toAlgebraicMatrix(final RealVector v) {
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

	private double[] toDesignRow(final Vector3D point) {
		final double x = point.getX();
		final double y = point.getY();
		final double z = point.getZ();
		// @formatter:off
		return new double[] {
				x * x, y * y, z * z,
				2 * x * y, 2 * x * z, 2 * y * z,
				2 * x, 2 * y, 2 * z
		};
		// @formatter:on
	}

	public static class Solution {

		private final boolean isEllipsoid;
		private final RealVector centre;
		private final Vector3D eigenU = new Vector3D(0.0, 0.0, 0.0);
		private final Vector3D eigenV = new Vector3D(0.0, 0.0, 0.0);
		private final Vector3D eigenW = new Vector3D(0.0, 0.0, 0.0);
		private double[] eigenvalues;
		private double a;
		private double b;
		private double c;

		public Solution(final RealMatrix quadricSurface) {
			centre = solveCentre(quadricSurface);
			final RealMatrix translatedSurface = translateToCentre(quadricSurface);
			solveEigenvectors(translatedSurface);
			solveRadii(translatedSurface);
			this.isEllipsoid = isEllipsoid(quadricSurface);
		}

		private boolean isEllipsoid(final RealMatrix quadricSurface) {
			// TODO implement
			// the signs of the eigenvalues (diagonal elements of DD) determine the
			// type. If they are all positive, it is an ellipsoid; two positive and
			// one negative gives a hyperboloid of one sheet; one positive and two
			// negative gives a hyperboloid of two sheets. If one or more eigenvalues
			// vanish, we have a degenerate case such as a paraboloid,or a cylinder or
			// even a pair of planes.
			return true;
		}

		private RealVector solveCentre(final RealMatrix quadricSurface) {
			final RealMatrix subMatrix = quadricSurface.getSubMatrix(0, 2, 0, 2);
			subMatrix.scalarMultiply(-1.0);
			final RealMatrix subInverse = MatrixUtils.inverse(subMatrix);
			final RealVector subVector = quadricSurface.getRowVector(3).getSubVector(
				0, 3);
			return subInverse.operate(subVector);
		}

		private void solveEigenvectors(final RealMatrix quadricSurface) {
			final RealMatrix eigenMatrix = quadricSurface.getSubMatrix(0, 2, 0, 2);
			final double scalar = -1.0 / quadricSurface.getEntry(3, 3);
			eigenMatrix.scalarMultiply(scalar);

		}

		private void solveRadii(final RealMatrix quadricSurface) {

		}

		private RealMatrix translateToCentre(final RealMatrix quadricSurface) {
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
	}
}
