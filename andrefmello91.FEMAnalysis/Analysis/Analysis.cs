using System;
using System.Collections.Generic;
using System.Linq;
using Extensions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using UnitsNet;
using UnitsNet.Units;
using static Extensions.UnitExtensions;

#nullable enable

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Linear analysis class.
	/// </summary>
	public class Analysis
	{
		#region Fields

		/// <summary>
		///     Get the <see cref="FEMAnalysis.InputData" />.
		/// </summary>
		public InputData InputData { get; }

		#endregion

		#region Properties

		/// <summary>
		///     Get/set the displacement <see cref="Vector" />.
		/// </summary>
		/// <remarks>
		///     Components in <see cref="LengthUnit.Millimeter" />.
		/// </remarks>
		public Vector<double> DisplacementVector { get; protected set; }

		/// <inheritdoc cref="InputData.ForceVector" />
		public Vector<double> ForceVector { get; protected set; }

		/// <summary>
		///     Get/set global stiffness <see cref="Matrix" />.
		/// </summary>
		public Matrix<double> GlobalStiffness { get; protected set; }

		#endregion

		#region Constructors

		/// <summary>
		///     Analysis base object.
		/// </summary>
		/// <param name="inputData">The <see cref="FEMAnalysis.InputData" /> for SPM analysis.</param>
		public Analysis(InputData inputData) => InputData = inputData;

		#endregion

		#region Methods

		/// <summary>
		///     Calculate the displacement <see cref="Vector" /> from an external force <see cref="Vector"/> and a global stiffness matrix.
		/// </summary>
		/// <param name="globalStiffness">Current global stiffness <see cref="Matrix" />.</param>
		/// <param name="forceVector">Current force <see cref="Vector" />.</param>
		/// <returns>
		///		The resultant displacement <see cref="Vector"/> or null if <paramref name="globalStiffness"/> and <paramref name="forceVector"/> sizes don't match.
		/// </returns>
		public static Vector<double>? CalculateDisplacements(Matrix<double> globalStiffness, Vector<double> forceVector) =>
			globalStiffness.IsSquare() && globalStiffness.RowCount == forceVector.Count
				? globalStiffness.Solve(forceVector)
				: null;

		/// <summary>
		///     Execute the analysis.
		/// </summary>
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis.ForceVector" /> (default: 1).</param>
		public void Execute(double loadFactor = 1)
		{
			// Set force vector
			ForceVector = InputData.ForceVector * loadFactor;

			// Assemble and simplify global stiffness and force vector
			UpdateStiffness();

			// Solve
			DisplacementVector = CalculateDisplacements(GlobalStiffness, ForceVector)!;

			// Set displacements to grips
			InputData.Grips.SetDisplacements(DisplacementVector);

			// Calculate element forces
			InputData.Elements.CalculateForces();

			// Set Reactions
			InputData.Grips.SetReactions(GetReactions());
		}

		/// <summary>
		///     Calculate the <see cref="Vector"/> of support reactions.
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
		///     Update <see cref="GlobalStiffness" />.
		/// </summary>
		/// <param name="simplify">Simplify stiffness and force vector? (default: true)</param>
		protected void UpdateStiffness(bool simplify = true)
		{
			// Initialize the global stiffness matrix
			GlobalStiffness = AssembleStiffness();

			// Simplify stiffness matrix
			if (simplify)
				Simplify();
		}

		/// <summary>
		///     Assemble the global stiffness <see cref="Matrix" />.
		/// </summary>
		protected Matrix<double> AssembleStiffness()
		{
			var n         = InputData.NumberOfDoFs;
			var stiffness = Matrix<double>.Build.Dense(n, n);

			InputData.Elements.AddToGlobalStiffness(stiffness);

			return stiffness;
		}

		/// <summary>
		///     Simplify <see cref="GlobalStiffness" /> and <see cref="ForceVector" />.
		/// </summary>
		/// <param name="simplifyByConstraints">Simplify matrix and vector by constraints?</param>
		/// <param name="simplifyZeroRows">Simplify matrix and vector on rows containing only zero elements?</param>
		protected void Simplify(bool simplifyByConstraints = true, bool simplifyZeroRows = true)
		{
			if (simplifyByConstraints)
				foreach (var i in InputData.ConstraintIndex)
				{
					// Clear the row and column [i] in the stiffness matrix (all elements will be zero)
					GlobalStiffness.ClearRow(i);
					GlobalStiffness.ClearColumn(i);

					// Set the diagonal element to 1
					GlobalStiffness[i, i] = 1;

					// Clear the row in the force vector
					ForceVector[i] = 0;

					// So ui = 0
				}

			if (simplifyZeroRows)
				foreach (var grip in InputData.Grips)
				{
					// Get DoF indexes
					var index = grip.DoFIndex;

					// Verify rows
					foreach (var i in index)
					{
						// Verify what line of the matrix is composed of zeroes
						if (GlobalStiffness.Row(i).Exists(num => !num.ApproxZero()))
							continue;

						// The row is composed of only zeroes, so the displacement must be zero
						// Set the diagonal element to 1
						GlobalStiffness[i, i] = 1;

						// Clear the row in the force vector
						ForceVector[i] = 0;
					}
				}

			// Approximate small numbers to zero
			GlobalStiffness.CoerceZero(1E-9);
		}

		public override string ToString() =>
			$"{InputData}\n" +
			"Global Stiffness:\n" +
			$"{GlobalStiffness}\n" +
			"Displacement Vector:\n" +
			$"{DisplacementVector}";

		#endregion
	}
}