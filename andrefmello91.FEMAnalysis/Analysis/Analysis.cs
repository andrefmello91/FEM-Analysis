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
		public IFEMInput FemInput { get; }

		/// <inheritdoc cref="FEMInput{TFiniteElement}.ForceVector" />
		public Vector<double>? ForceVector { get; protected set; }

		/// <summary>
		///     Get/set global stiffness <see cref="Matrix" />.
		/// </summary>
		public virtual Matrix<double>? GlobalStiffness { get; protected set; }

		/// <summary>
		///     The DoF indexes of known displacements for displacement control.
		/// </summary>
		protected List<int> KnownDisplacementIndex { get; } = new();

		#endregion

		#region Constructors

		/// <summary>
		///     Base analysis constructor.
		/// </summary>
		/// <param name="femInput">The <see cref="IFEMInput{TFiniteElement}" /> for finite element analysis.</param>
		public Analysis(IFEMInput femInput) => FemInput = femInput;

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
		///     Get the force vector simplified by constraints.
		/// </summary>
		/// <param name="forceVector">The global force vector.</param>
		/// <param name="constraintIndexes">The constrained DoF indexes.</param>
		public static Vector<double> SimplifiedForces(Vector<double> forceVector, IEnumerable<int> constraintIndexes)
		{
			var simplifiedForces = forceVector.Clone();

			// Clear the row in the force vector
			foreach (var i in constraintIndexes)
				simplifiedForces[i] = 0;

			return simplifiedForces;
		}


		/// <summary>
		///     Get the force vector simplified by constraints and known displacements.
		/// </summary>
		/// <inheritdoc cref="SimplifiedForces(Vector{double}, IEnumerable{int})" />
		/// <param name="stiffness">The full global stiffness <see cref="Matrix{T}" /> (not simplified).</param>
		/// <param name="appliedDisplacements">The vector of known applied displacements.</param>
		public static Vector<double> SimplifiedForces(Vector<double> forceVector, IEnumerable<int> constraintIndexes, Matrix<double> stiffness, Vector<double> appliedDisplacements)
		{
			// Simplify by constraints
			var simplifiedForces = forceVector.Clone();
			var knownDispIndex   = new List<int>();

			// Clear the row in the force vector
			for (var i = 0; i < appliedDisplacements.Count; i++)
			{
				var di = appliedDisplacements[i];

				if (di.ApproxZero())
					continue;

				// Add to known displacement index
				knownDispIndex.Add(i);

				// Subtract equivalent forces from force vector
				simplifiedForces -= stiffness.Column(i) * di;
			}

			// Set displacements to force vector
			foreach (var i in knownDispIndex)
				simplifiedForces[i] = appliedDisplacements[i];


			return SimplifiedForces(simplifiedForces, constraintIndexes);
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
		protected virtual void UpdateStiffness() =>
			GlobalStiffness = FemInput.AssembleStiffness();

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