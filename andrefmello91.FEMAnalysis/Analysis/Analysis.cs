using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using UnitsNet.Units;
#nullable enable

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Linear analysis class.
	/// </summary>
	/// <typeparam name="TFiniteElement">Any type that implements <see cref="IFiniteElement" />.</typeparam>
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
		public virtual Vector<double>? DisplacementVector { get; protected set; }

		/// <summary>
		///     Get the input for finite element analysis.
		/// </summary>
		public IFEMInput<TFiniteElement> FemInput { get; }

		/// <inheritdoc cref="FEMInput{TFiniteElement}.ForceVector" />
		public Vector<double>? ForceVector { get; protected set; }

		/// <summary>
		///     Get/set global stiffness <see cref="Matrix" />.
		/// </summary>
		public virtual Matrix<double>? GlobalStiffness { get; protected set; }

		#endregion

		#region Constructors

		/// <summary>
		///     Base analysis constructor.
		/// </summary>
		/// <param name="femInput">The <see cref="IFEMInput{TFiniteElement}" /> for finite element analysis.</param>
		public Analysis(IFEMInput<TFiniteElement> femInput) => FemInput = femInput;

		#endregion

		#region Methods

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
		///     Simplify the stiffness matrix and force vector at constrained DoFs.
		/// </summary>
		/// <param name="stiffness">The global stiffness <see cref="Matrix{T}" />.</param>
		/// <param name="forceVector">The global force <see cref="Vector{T}" />.</param>
		/// <param name="constraintIndex">The index of constrained DoFs.</param>
		/// <param name="simplifyZeroRows">Simplify matrix and vector at rows containing only zero elements?</param>
		internal static void Simplify(Matrix<double> stiffness, Vector<double>? forceVector, IEnumerable<int> constraintIndex, bool simplifyZeroRows = true)
		{
			var index = constraintIndex.ToArray();

			// Clear the rows and columns in the stiffness matrix
			stiffness.ClearRows(index);
			stiffness.ClearColumns(index);

			foreach (var i in index)
			{
				// Set the diagonal element to 1
				stiffness[i, i] = 1;

				// Clear the row in the force vector
				if (forceVector != null)
					forceVector[i] = 0;
			}

			if (simplifyZeroRows)
				for (var i = 0; i < stiffness.RowCount; i++)
				{
					// Verify what line of the matrix is composed of zeroes
					if (stiffness.Row(i).Exists(num => !num.ApproxZero(1E-9)))
						continue;

					// The row is composed of only zeroes, so the displacement must be zero
					// Set the diagonal element to 1
					stiffness[i, i] = 1;

					// Clear the row in the force vector
					if (forceVector != null)
						forceVector[i] = 0;
				}

			// Approximate small numbers to zero
			stiffness!.CoerceZero(1E-9);
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
		///     Update <see cref="GlobalStiffness" />.
		/// </summary>
		/// <param name="simplify">Simplify stiffness and force vector? (default: true)</param>
		protected virtual void UpdateStiffness(bool simplify = true)
		{
			// Initialize the global stiffness matrix
			GlobalStiffness = FemInput.AssembleStiffness();

			// Simplify stiffness matrix
			if (simplify)
				Simplify(GlobalStiffness, ForceVector!, FemInput.ConstraintIndex);
		}

		#region Object override

		/// <inheritdoc />
		public override string ToString() =>
			$"{FemInput}\n" +
			"Global Stiffness:\n" +
			$"{GlobalStiffness}\n" +
			"Displacement Vector:\n" +
			$"{DisplacementVector}";

		#endregion

		#endregion

	}
}