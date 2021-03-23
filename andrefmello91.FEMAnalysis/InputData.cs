using System.Collections.Generic;
using System.Linq;
using Extensions;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using UnitsNet.Units;
#nullable disable

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Input Data class.
	/// </summary>
	public class InputData
	{
		#region Properties

		/// <summary>
		///     Get the index of constrained degrees of freedom.
		/// </summary>
		public List<int> ConstraintIndex { get; }

		/// <summary>
		///     Get the elements of the finite element model.
		/// </summary>
		public IFiniteElement[] Elements { get; }

		/// <summary>
		///     Get the external force <see cref="Vector" />.
		/// </summary>
		/// <remarks>
		///     Components in <see cref="ForceUnit.Newton" />.
		/// </remarks>
		public Vector<double> ForceVector { get; }

		/// <summary>
		///     Get the grips of the finite element model.
		/// </summary>
		public IGrip[] Grips { get; }

		/// <summary>
		///     Get the number of degrees of freedom (DoFs).
		/// </summary>
		public int NumberOfDoFs { get; }

		#endregion
		#region Constructors

		/// <summary>
		///     Input Data constructor.
		/// </summary>
		/// <param name="elements">The collection containing all distinct <see cref="IFiniteElement" />'s in the model.</param>
		public InputData(IEnumerable<IFiniteElement> elements)
		{
			Elements        = elements.OrderBy(e => e.Number).ToArray();
			Grips           = elements.SelectMany(e => e.Grips).Distinct().OrderBy(g => g.Number).ToArray();
			NumberOfDoFs    = 2 * Grips.Length;
			ForceVector     = GetForceVector(Grips);
			ConstraintIndex = GetConstraintIndex(Grips);
		}

		#endregion
		#region Methods

		/// <summary>
		///     Get the indexes of constrained degrees of freedom from a collection of grips.
		/// </summary>
		/// <inheritdoc cref="GetForceVector" />
		public static List<int> GetConstraintIndex(IEnumerable<IGrip> grips)
		{
			var constraintList = new List<int>();

			foreach (var node in grips)
			{
				// Get DoF indexes
				var index = node.DoFIndex;
				int
					i = index[0],
					j = index[1];

				var constraint = node.Constraint;

				if (constraint.X)

					// There is a support in X direction
					constraintList.Add(i);

				if (constraint.Y)

					// There is a support in Y direction
					constraintList.Add(j);
			}

			return constraintList;
		}

		/// <summary>
		///     Get the force <see cref="Vector" /> from a collection of grips.
		/// </summary>
		/// <param name="grips">The collection of distinct grips of the finite element model.</param>
		/// <inheritdoc cref="ForceVector" />
		public static Vector<double> GetForceVector(IEnumerable<IGrip> grips)
		{
			// Initialize the force vector
			var f = new double[2 * grips.Count()];

			// Read the nodes data
			foreach (var grip in grips)
			{
				// Get DoF indexes
				var index = grip.DoFIndex;
				int
					i = index[0],
					j = index[1];

				// Set to force vector
				f[i] = grip.Force.X.Newtons;
				f[j] = grip.Force.Y.Newtons;
			}

			return f.ToVector();
		}

		public override string ToString() =>
			$"Number of grips: {Grips.Length}\n" +
			$"Number of elements: {Elements.Length}\n" +
			$"Force vector: \n{ForceVector}\n" +
			$"Constraint Index: {ConstraintIndex.Select(i => i.ToString()).Aggregate((i, f) => $"{i} - {f}")}";

		#endregion
	}
}