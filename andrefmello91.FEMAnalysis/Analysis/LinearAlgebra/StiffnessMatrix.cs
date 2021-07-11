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
	public class StiffnessMatrix : QuantityMatrix<ForcePerLength, ForcePerLengthUnit>
	{

		#region Properties

		/// <summary>
		///     Default tolerance for stiffness matrix.
		/// </summary>
		private static ForcePerLength Tolerance { get; } = ForcePerLength.FromNewtonsPerMillimeter(1E-6);

		#endregion

		#region Constructors

		/// <inheritdoc cref="StiffnessMatrix(Matrix{double}, ForcePerLengthUnit)" />
		public StiffnessMatrix(double[,] values, ForcePerLengthUnit unit = ForcePerLengthUnit.NewtonPerMillimeter)
			: base(values, unit)
		{
		}

		/// <inheritdoc cref="StiffnessMatrix(Matrix{double}, ForcePerLengthUnit)" />
		public StiffnessMatrix(ForcePerLength[,] values)
			: base(values)
		{
		}

		/// <summary>
		///     Create a stiffness matrix.
		/// </summary>
		/// <param name="value">The <see cref="Matrix{T}" /> or array value.</param>
		/// <param name="unit">The unit of <paramref name="value" />'s components</param>
		public StiffnessMatrix(Matrix<double> value, ForcePerLengthUnit unit = ForcePerLengthUnit.NewtonPerMillimeter)
			: base(value, unit)
		{
		}

		#endregion

		#region Methods

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
		///     Create a stiffness matrix with zero elements.
		/// </summary>
		/// <param name="size">The size of the matrix.</param>
		public new static StiffnessMatrix Zero(int size) => new(new double[size, size]);

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

		/// <inheritdoc cref="QuantityMatrix{TQuantity,TUnit}.Clone" />
		public override QuantityMatrix<ForcePerLength, ForcePerLengthUnit> Clone() => new StiffnessMatrix(this, Unit)
		{
			ConstraintIndex = ConstraintIndex
		};

		/// <inheritdoc cref="IUnitConvertible{TUnit}.Convert" />
		public override QuantityMatrix<ForcePerLength, ForcePerLengthUnit> Convert(ForcePerLengthUnit unit) => Unit.Equals(unit)
			? Clone()
			: new StiffnessMatrix(this, Unit)
			{
				ConstraintIndex = ConstraintIndex
			};

		/// <inheritdoc />
		public override QuantityMatrix<ForcePerLength, ForcePerLengthUnit> Simplified() => Simplified(Tolerance);

		/// <inheritdoc />
		public override QuantityMatrix<ForcePerLength, ForcePerLengthUnit> Simplified(double? threshold)
		{
			var value = Clone();

			if (threshold.HasValue)
				value.CoerceZero(threshold.Value);

			return ConstraintIndex is not null
				? new StiffnessMatrix(SimplifiedStiffness(value, ConstraintIndex), Unit)
				: value;
		}

		/// <summary>
		///     Solve a system d = K f.
		/// </summary>
		/// <param name="forceVector">The force vector.</param>
		/// <param name="useSimplified">Use simplified matrix and vector?</param>
		/// <param name="unit">The unit to return.</param>
		/// <returns>
		///     The resulting displacement vector.
		/// </returns>
		public DisplacementVector Solve(ForceVector forceVector, bool useSimplified = true, LengthUnit unit = LengthUnit.Millimeter)
		{
			// Convert
			var stiffnessMatrix = Unit is ForcePerLengthUnit.NewtonPerMillimeter
				? this
				: Convert(ForcePerLengthUnit.NewtonPerMillimeter);

			var forces = forceVector.Unit is ForceUnit.Newton
				? forceVector
				: forceVector.Convert(ForceUnit.Newton);

			var k = useSimplified
				? stiffnessMatrix.Simplified()
				: stiffnessMatrix;

			var f = useSimplified
				? forces.Simplified(ConstraintIndex)
				: forces;

			// Solve
			var d  = k.Solve(f);
			var dv = new DisplacementVector(d);

			return
				unit is LengthUnit.Millimeter
					? dv
					: (DisplacementVector) dv.Convert(unit);
		}

		/// <summary>
		///     Solve a system d = K f.
		/// </summary>
		/// <param name="displacementVector">The displacement vector.</param>
		/// <param name="useSimplified">Use simplified matrix and vector?</param>
		/// <returns>
		///     The resulting force vector..
		/// </returns>
		public ForceVector Solve(DisplacementVector displacementVector, bool useSimplified = true, ForceUnit unit = ForceUnit.Newton)
		{
			// Convert
			var stiffnessMatrix = Unit is ForcePerLengthUnit.NewtonPerMillimeter
				? this
				: Convert(ForcePerLengthUnit.NewtonPerMillimeter);

			var displacements = displacementVector.Unit is LengthUnit.Millimeter
				? displacementVector
				: displacementVector.Convert(LengthUnit.Millimeter);

			Matrix<double> k = useSimplified
				? stiffnessMatrix.Simplified()
				: stiffnessMatrix;

			Vector<double> d = useSimplified
				? displacements.Simplified(ConstraintIndex)
				: displacements;

			// Solve
			var f  = k * d;
			var fv = new ForceVector(f);

			return
				unit is ForceUnit.Newton
					? fv
					: (ForceVector) fv.Convert(unit);
		}

		/// <inheritdoc />
		public override QuantityMatrix<ForcePerLength, ForcePerLengthUnit> Transform(Matrix<double> transformationMatrix)
		{
			var value = transformationMatrix.Transpose() * this * transformationMatrix;

			return
				new StiffnessMatrix(value, Unit)
				{
					ConstraintIndex = ConstraintIndex
				};
		}

		/// <inheritdoc />
		public override QuantityMatrix<ForcePerLength, ForcePerLengthUnit> Transpose() => new StiffnessMatrix(Build.OfStorage(Storage).Transpose(), Unit);

		#region Object override

		/// <summary>
		///     Solve the displacement vector by multiplying the inverse of stiffness and the force vector.
		/// </summary>
		/// <remarks>
		///     This uses the simplified stiffness matrix and forces.
		/// </remarks>
		/// <returns>
		///     The <see cref="DisplacementVector" /> with components in <see cref="LengthUnit.Millimeter" />.
		/// </returns>
		public static DisplacementVector operator /(StiffnessMatrix stiffnessMatrix, ForceVector forceVector) => stiffnessMatrix.Solve(forceVector);


		/// <summary>
		///     Solve the force vector by multiplying the stiffness and the displacement vector.
		/// </summary>
		/// <remarks>
		///     This uses the simplified stiffness matrix and displacements.
		/// </remarks>
		/// <returns>
		///     The <see cref="ForceVector" /> with components in <see cref="ForceUnit.Newton" />.
		/// </returns>
		public static ForceVector operator *(StiffnessMatrix stiffnessMatrix, DisplacementVector displacementVector) => stiffnessMatrix.Solve(displacementVector);

		#endregion

		#endregion

	}
}