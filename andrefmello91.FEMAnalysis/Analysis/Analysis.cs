﻿using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using UnitsNet.Units;
#nullable enable

namespace andrefmello91.FEMAnalysis
{
    /// <summary>
    ///     Linear analysis class.
    /// </summary>
    /// <typeparam name="TFiniteElement">Any type that implements <see cref="IFiniteElement"/>.</typeparam>
    public abstract class Analysis<TFiniteElement>
		where TFiniteElement : IFiniteElement
	{
		#region Properties

        /// <summary>
        ///     Get/set the displacement <see cref="Vector" />.
        /// </summary>
        /// <remarks>
        ///     Components in <see cref="LengthUnit.Millimeter" />.
        /// </remarks>
        public Vector<double>? DisplacementVector { get; protected set; }

        /// <summary>
        ///     Get the <see cref="FEMInput" />.
        /// </summary>
        public FEMInput<TFiniteElement> FemInput { get; }

		/// <inheritdoc cref="FEMInput{TFiniteElement}.ForceVector" />
		public Vector<double>? ForceVector { get; protected set; }

        /// <summary>
        ///     Get/set global stiffness <see cref="Matrix" />.
        /// </summary>
        public Matrix<double>? GlobalStiffness { get; protected set; }

		#endregion

		#region Constructors

        /// <summary>
        ///     Base analysis constructor.
        /// </summary>
        /// <param name="femInput">The <see cref="FEMInput{TFiniteElement}" /> for finite element analysis.</param>
        public Analysis(FEMInput<TFiniteElement> femInput) => FemInput = femInput;

		#endregion

		#region Methods

        /// <summary>
        ///     Assemble the global stiffness <see cref="Matrix" />.
        /// </summary>
        /// <param name="femInput">The <see cref="FEMInput{TFiniteElement}" /></param>
        public static Matrix<double> AssembleStiffness(FEMInput<TFiniteElement> femInput)
		{
			var n         = femInput.NumberOfDoFs;
			var stiffness = Matrix<double>.Build.Dense(n, n);

			stiffness.AddStiffness(femInput.Elements);

			return stiffness;
		}

        /// <summary>
        ///     Calculate the displacement <see cref="Vector" /> from an external force <see cref="Vector" /> and a global
        ///     stiffness matrix.
        /// </summary>
        /// <param name="globalStiffness">Current global stiffness <see cref="Matrix" />.</param>
        /// <param name="forceVector">Current force <see cref="Vector" />.</param>
        /// <returns>
        ///     The resultant displacement <see cref="Vector" /> or null if <paramref name="globalStiffness" /> and
        ///     <paramref name="forceVector" /> sizes don't match.
        /// </returns>
        public static Vector<double>? CalculateDisplacements(Matrix<double> globalStiffness, Vector<double> forceVector) =>
			globalStiffness.IsSquare() && globalStiffness.RowCount == forceVector.Count
				? globalStiffness.Solve(forceVector)
				: null;

        /// <summary>
        ///     Get the internal force <see cref="Vector" />.
        /// </summary>
        /// <param name="femInput">The <see cref="FemInput" />.</param>
        /// <param name="simplify">Simplify vector in constraint indexes?</param>
        public static Vector<double> InternalForces(FEMInput<TFiniteElement> femInput, bool simplify = true)
		{
			var iForces = Vector<double>.Build.Dense(femInput.NumberOfDoFs);

			iForces.AddInternalForces(femInput.Elements);

			if (!simplify)
				return iForces;

			foreach (var i in femInput.ConstraintIndex)
				iForces[i] = 0;

			return iForces;
		}

        /// <summary>
        ///     Calculate the <see cref="Vector" /> of support reactions.
        /// </summary>
        public Vector<double> GetReactions()
		{
			// Calculate forces and reactions
			var f = GlobalStiffness * DisplacementVector;

			// Subtract forces
			f -= ForceVector;

			f.CoerceZero(1E-9);

			return f;
		}

        /// <summary>
        ///     Simplify <see cref="GlobalStiffness" /> and <see cref="ForceVector" />.
        /// </summary>
        /// <param name="simplifyByConstraints">Simplify matrix and vector by constraints?</param>
        /// <param name="simplifyZeroRows">Simplify matrix and vector on rows containing only zero elements?</param>
        protected void Simplify(bool simplifyByConstraints = true, bool simplifyZeroRows = true)
		{
			if (simplifyByConstraints)
			{
				// Clear the rows and columns in the stiffness matrix
				GlobalStiffness!.ClearRows(FemInput.ConstraintIndex.ToArray());
				GlobalStiffness!.ClearColumns(FemInput.ConstraintIndex.ToArray());

				foreach (var i in FemInput.ConstraintIndex)
				{
					// Set the diagonal element to 1
					GlobalStiffness[i, i] = 1;

					// Clear the row in the force vector
					ForceVector![i] = 0;
				}
			}

			if (simplifyZeroRows)
				for (var i = 0; i < GlobalStiffness!.RowCount; i++)
				{
					// Verify what line of the matrix is composed of zeroes
					if (GlobalStiffness!.Row(i).Exists(num => !num.ApproxZero()))
						continue;

					// The row is composed of only zeroes, so the displacement must be zero
					// Set the diagonal element to 1
					GlobalStiffness[i, i] = 1;

					// Clear the row in the force vector
					ForceVector![i] = 0;
				}

			// Approximate small numbers to zero
			GlobalStiffness!.CoerceZero(1E-9);
		}

        /// <summary>
        ///     Update <see cref="GlobalStiffness" />.
        /// </summary>
        /// <param name="simplify">Simplify stiffness and force vector? (default: true)</param>
        protected virtual void UpdateStiffness(bool simplify = true)
		{
			// Initialize the global stiffness matrix
			GlobalStiffness = AssembleStiffness(FemInput);

			// Simplify stiffness matrix
			if (simplify)
				Simplify();
		}

		/// <inheritdoc />
		public override string ToString() =>
			$"{FemInput}\n" +
			"Global Stiffness:\n" +
			$"{GlobalStiffness}\n" +
			"Displacement Vector:\n" +
			$"{DisplacementVector}";

		#endregion

	}
}