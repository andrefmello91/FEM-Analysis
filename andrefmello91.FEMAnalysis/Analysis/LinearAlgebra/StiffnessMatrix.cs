using System;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;
using UnitsNet.Units;

namespace andrefmello91.FEMAnalysis
{

	/// <summary>
	///     Default stiffness matrix class.
	/// </summary>
	/// <remarks>
	///     Unit is <see cref="ForcePerLengthUnit" />.
	///     <para>
	///         Quantity is <see cref="ForcePerLength" />.
	///     </para>
	/// </remarks>
	public class StiffnessMatrix : StiffnessMatrix<ForcePerLength, ForcePerLengthUnit>
	{
		/// <summary>
		///		Default tolerance for stiffness matrix.
		/// </summary>
		private static ForcePerLength Tolerance { get; } = ForcePerLength.FromNewtonsPerMillimeter(1E-6);

		#region Constructors

		/// <inheritdoc cref="StiffnessMatrix(Matrix{double}, ForcePerLengthUnit)" />
		public StiffnessMatrix(double[,] value, ForcePerLengthUnit unit = ForcePerLengthUnit.NewtonPerMillimeter)
			: this(Matrix<double>.Build.DenseOfArray(value), unit)
		{
		}

		/// <summary>
		///     Create a stiffness matrix.
		/// </summary>
		/// <param name="value">The <see cref="Matrix{T}" /> or <see cref="double" /> array value.</param>
		/// <param name="unit">The unit of <paramref name="value" />'s components</param>
		public StiffnessMatrix(Matrix<double> value, ForcePerLengthUnit unit = ForcePerLengthUnit.NewtonPerMillimeter)
			: base(value, unit)
		{
		}

		#endregion

		#region Methods

# if NET5_0
		/// <inheritdoc cref="IUnitConvertible{T}.Convert" />
		public override StiffnessMatrix Convert(ForcePerLengthUnit unit) => Unit.Equals(unit)
			? Clone()
			: new StiffnessMatrix(Values.GetQuantities<ForcePerLength, ForcePerLengthUnit>(Unit).GetValues(unit), unit)
			{
				ConstraintIndex = ConstraintIndex
			};

		/// <inheritdoc cref="StiffnessMatrix{TQuantity,TUnit}.Clone"/>
		public override StiffnessMatrix Clone() => new(Values.ToMatrix(), Unit)
		{
			ConstraintIndex = ConstraintIndex
		};

		/// <inheritdoc />
		public override StiffnessMatrix Transpose() => new (Values.ToMatrix().Transpose(), Unit);
		
		/// <inheritdoc />
		public override StiffnessMatrix Transform(Matrix<double> transformationMatrix)
		{
			var value = transformationMatrix.Transpose() * Values.ToMatrix() * transformationMatrix;
			
			return
				new StiffnessMatrix(value, Unit)
				{
					ConstraintIndex = ConstraintIndex
				};
		}

#else

		/// <inheritdoc cref="IUnitConvertible{T}.Convert" />
		public override StiffnessMatrix<ForcePerLength, ForcePerLengthUnit> Convert(ForcePerLengthUnit unit) => Unit.Equals(unit)
			? Clone()
			: new StiffnessMatrix(Values.GetQuantities<ForcePerLength, ForcePerLengthUnit>(Unit).GetValues(unit), unit)
			{
				ConstraintIndex = ConstraintIndex
			};

		/// <inheritdoc cref="StiffnessMatrix{TQuantity,TUnit}.Clone" />
		public override StiffnessMatrix<ForcePerLength, ForcePerLengthUnit> Clone() => new StiffnessMatrix(Values.ToMatrix(), Unit)
		{
			ConstraintIndex = ConstraintIndex
		};

		/// <inheritdoc />
		public override StiffnessMatrix<ForcePerLength, ForcePerLengthUnit> Transpose() => new StiffnessMatrix(Values.ToMatrix().Transpose(), Unit);

		/// <inheritdoc />
		public override StiffnessMatrix<ForcePerLength, ForcePerLengthUnit> Transform(Matrix<double> transformationMatrix)
		{
			var value = transformationMatrix.Transpose() * Values.ToMatrix() * transformationMatrix;

			return
				new StiffnessMatrix(value, Unit)
				{
					ConstraintIndex = ConstraintIndex
				};
		}

#endif
		/// <inheritdoc />
		public override Matrix<double> Simplified() => Simplified(Tolerance);

		/// <summary>
		///     Solve a system d = K f.
		/// </summary>
		/// <param name="forceVector"></param>
		/// <param name="useSimplified">Use simplified matrix and vector?</param>
		/// <returns>
		///     The resulting displacement vector, with components in <see cref="LengthUnit.Millimeter" />.
		/// </returns>
		public DisplacementVector Solve(ForceVector forceVector, bool useSimplified = true)
		{
			// Convert
			var stiffnessMatrix = Unit is ForcePerLengthUnit.NewtonPerMillimeter
				? this
				: Convert(ForcePerLengthUnit.NewtonPerMillimeter);

			var forces = forceVector.Unit is ForceUnit.Newton
				? forceVector
				: forceVector.Convert(ForceUnit.Newton);

			Matrix<double> k = useSimplified
				? stiffnessMatrix.Simplified()
				: stiffnessMatrix;

			var f = useSimplified
				? forces
				: forces.Simplified();

			// Solve
			var d = k.Solve(f);

			return
				new DisplacementVector(d);
		}

		#endregion

		#region Object override

		#endregion

		/// <summary>
		///     Create a stiffness matrix with zero elements.
		/// </summary>
		/// <param name="size">The size of the matrix.</param>
		public static StiffnessMatrix Zero(int size) => new(Matrix<double>.Build.Dense(size, size));

		/// <summary>
		///     Assemble the global stiffness matrix.
		/// </summary>
		/// <param name="femInput">The <see cref="FEMInput" /></param>
		public static StiffnessMatrix Assemble(IFEMInput femInput)
		{
			var stiffness = Zero(femInput.NumberOfDoFs);
			stiffness.ConstraintIndex = femInput.ConstraintIndex;

			foreach (var element in femInput)
			{
				var dofIndex = element.DoFIndex;

				for (var i = 0; i < dofIndex.Length; i++)
				{
					// Global index
					var k = dofIndex[i];

					for (var j = 0; j < dofIndex.Length; j++)
					{
						// Global index
						var l = dofIndex[j];

						stiffness[k, l] += element.Stiffness[i, j];
					}
				}
			}

			return stiffness;
		}

		/// <summary>
		///     Get the global stiffness simplified.
		/// </summary>
		/// <param name="stiffness">The global stiffness <see cref="Matrix{T}" /> to simplify.</param>
		/// <param name="indexes">The DoF indexes to simplify matrix.</param>
		/// <param name="simplifyZeroRows">Simplify matrix at rows containing only zero elements?</param>
		public static Matrix<double> SimplifiedStiffness(Matrix<double> stiffness, IEnumerable<int> indexes, bool simplifyZeroRows = true)
		{
			var simplifiedStiffness = stiffness.Clone();

			var index = indexes.ToArray();

			// Clear the rows and columns in the stiffness matrix
			simplifiedStiffness.ClearRows(index);
			simplifiedStiffness.ClearColumns(index);

			// Set the diagonal element to 1
			foreach (var i in index)
				simplifiedStiffness[i, i] = 1;

			if (!simplifyZeroRows)
				return simplifiedStiffness;

			// Verify rows
			for (var i = 0; i < simplifiedStiffness.RowCount; i++)
			{
				// Verify what line of the matrix is composed of zeroes
				if (simplifiedStiffness.Row(i).Exists(num => !num.ApproxZero()) && simplifiedStiffness.Column(i).Exists(num => !num.ApproxZero()))
					continue;

				// The row is composed of only zeroes, so the displacement must be zero
				// Set the diagonal element to 1
				simplifiedStiffness[i, i] = 1;
			}

			return simplifiedStiffness;
		}

		/// <summary>
		///     Calculate the stiffness increment for nonlinear analysis.
		/// </summary>
		/// <param name="solver">The nonlinear solver.</param>
		/// <returns>
		///     The <see cref="StiffnessMatrix" /> increment with current unit.
		/// </returns>
		/// <inheritdoc cref="TangentIncrement" />
		public static StiffnessMatrix StiffnessIncrement(IIteration currentIteration, IIteration lastIteration, NonLinearSolver solver) =>
			solver switch
			{
				NonLinearSolver.Secant => SecantIncrement(currentIteration, lastIteration),
				_                      => TangentIncrement(currentIteration, lastIteration)
			};

		/// <summary>
		///     Calculate the secant stiffness increment.
		/// </summary>
		/// <inheritdoc cref="TangentIncrement(andrefmello91.FEMAnalysis.IIteration, andrefmello91.FEMAnalysis.IIteration)" />
		private static StiffnessMatrix SecantIncrement(IIteration currentIteration, IIteration lastIteration)
		{
			// Calculate the variation of displacements and residual as vectors
			Vector<double>
				dU = (currentIteration.Displacements - lastIteration.Displacements).Convert(LengthUnit.Millimeter),
				dR = (currentIteration.ResidualForces - lastIteration.ResidualForces).Convert(ForceUnit.Newton);

			var unit = currentIteration.Stiffness.Unit;

			Matrix<double> stiffness = unit is ForcePerLengthUnit.NewtonPerMillimeter
				? lastIteration.Stiffness
				: lastIteration.Stiffness.Convert(ForcePerLengthUnit.NewtonPerMillimeter);

			var inc = ((dR - stiffness * dU) / dU.Norm(2)).ToColumnMatrix() * dU.ToRowMatrix();

			var increment = new StiffnessMatrix(inc)
			{
				ConstraintIndex = currentIteration.Stiffness.ConstraintIndex
			};

			return
				unit is ForcePerLengthUnit.NewtonPerMillimeter
					? increment
					: (StiffnessMatrix) increment.Convert(unit);
		}

		/// <summary>
		///     Calculate the tangent stiffness increment.
		/// </summary>
		/// <param name="currentIteration">The current iteration.</param>
		/// <param name="lastIteration">The last solved iteration.</param>
		/// <returns>
		///     The <see cref="andrefmello91.FEMAnalysis.StiffnessMatrix" /> increment with current unit.
		/// </returns>
		private static StiffnessMatrix TangentIncrement(IIteration currentIteration, IIteration lastIteration)
		{
			// Get variations
			var dF = (currentIteration.InternalForces - lastIteration.InternalForces).Convert(ForceUnit.Newton);
			var dU = (currentIteration.Displacements - lastIteration.Displacements).Convert(LengthUnit.Millimeter);

			var inc = dF.ToColumnMatrix() * dU.ToRowMatrix();

			var unit = currentIteration.Stiffness.Unit;

			var increment = new StiffnessMatrix(inc)
			{
				ConstraintIndex = currentIteration.Stiffness.ConstraintIndex
			};

			return
				unit is ForcePerLengthUnit.NewtonPerMillimeter
					? increment
					: (StiffnessMatrix) increment.Convert(unit);
		}


		/// <summary>
		///     Create a force vector by multiplying the stiffness and the displacement vector.
		/// </summary>
		/// <remarks>
		///     This uses the simplified stiffness matrix and displacements.
		/// </remarks>
		/// <returns>
		///     The <see cref="ForceVector" /> with components in <see cref="ForceUnit.Newton" />.
		/// </returns>
		public static ForceVector operator *(StiffnessMatrix stiffnessMatrix, DisplacementVector displacementVector)
		{
			// Convert
			var k = stiffnessMatrix.Convert(ForcePerLengthUnit.NewtonPerMillimeter).Simplified();
			var d = displacementVector.Convert(LengthUnit.Millimeter).Simplified();

			// Multiply
			var f = k * d;

			return
				new ForceVector(f);
		}

		/// <returns>
		///     A new stiffness matrix with summed components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentOutOfRangeException">If left and right don't have the same dimensions.</exception>
		public static StiffnessMatrix operator +(StiffnessMatrix left, StiffnessMatrix right) =>
			new(left.Values.ToMatrix() + right.Convert(left.Unit), left.Unit)
			{
				ConstraintIndex = left.ConstraintIndex ?? right.ConstraintIndex
			};

		/// <returns>
		///     A new stiffness matrix with subtracted components in <paramref name="left" />'s unit.
		/// </returns>
		/// <exception cref="ArgumentOutOfRangeException">If left and right don't have the same dimensions.</exception>
		public static StiffnessMatrix operator -(StiffnessMatrix left, StiffnessMatrix right) =>
			new(left.Values.ToMatrix() - right.Convert(left.Unit), left.Unit)
			{
				ConstraintIndex = left.ConstraintIndex ?? right.ConstraintIndex
			};

		/// <returns>
		///     A new stiffness matrix with components multiplied by a value
		/// </returns>
		public static StiffnessMatrix operator *(double multiplier, StiffnessMatrix matrix) => new(multiplier * matrix.Values.ToMatrix(), matrix.Unit)
		{
			ConstraintIndex = matrix.ConstraintIndex
		};

		/// <inheritdoc cref="op_Multiply(double, StiffnessMatrix) " />
		public static StiffnessMatrix operator *(StiffnessMatrix matrix, double multiplier) => multiplier * matrix;

		/// <inheritdoc cref="Matrix{T}.op_Division(Matrix{T}, T)" />
		public static StiffnessMatrix operator /(StiffnessMatrix matrix, double divisor) => new(matrix.Values.ToMatrix() / divisor, matrix.Unit)
		{
			ConstraintIndex = matrix.ConstraintIndex
		};

		/// <inheritdoc cref="Matrix{T}.op_UnaryNegation" />
		public static StiffnessMatrix operator -(StiffnessMatrix matrix) => new(-matrix.Values.ToMatrix(), matrix.Unit)
		{
			ConstraintIndex = matrix.ConstraintIndex
		};
	}
}