using System;
using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

#nullable enable

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///		Extensions for <see cref="IGrip"/> and <see cref="IFiniteElement"/>.
	/// </summary>
	public static class Extensions
	{
		/// <summary>
		/// Return an <see cref="INumberedElement"/> of a collection, in given <paramref name="number"/>.
		/// </summary>
		/// <param name="elements">The collection of <see cref="INumberedElement"/>'s.</param>
		/// <param name="number">The number of the element wanted.</param>
		public static INumberedElement? GetByNumber(this IEnumerable<INumberedElement> elements, int number) => elements.First(element => number == element.Number);

		/// <summary>
		/// Get global indexes of the degrees of freedom af a collection of grips.
		/// </summary>
		/// <param name="grips">A collection of <see cref="IGrip"/>'s.</param>
		public static IEnumerable<int> GlobalIndexes(IEnumerable<IGrip> grips) => grips.SelectMany(GlobalIndexes);

		/// <summary>
		/// Get global indexes of the degrees of freedom af a grip.
		/// </summary>
		/// <param name="grip">An <see cref="IGrip"/>.</param>
		public static IEnumerable<int> GlobalIndexes(IGrip grip)
		{
			var n = 2 * grip.Number;

			yield return n - 2;
			yield return n - 1;
		}
		
		/// <summary>
		///     Add the stiffness of an <see cref="IFiniteElement"/> to global stiffness <see cref="Matrix" />.
		/// </summary>
		/// <param name="element">The <see cref="IFiniteElement" /> to add to <paramref name="globalStiffness" />.</param>
		/// <param name="globalStiffness">The global stiffness <see cref="Matrix" />.</param>
		public static void AddToGlobalStiffness([NotNull] this IFiniteElement element, [NotNull] Matrix<double> globalStiffness)
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

					globalStiffness[k, l] += element.Stiffness[i, j];
				}
			}
		}

		/// <summary>
		///     Add the stiffness of each <see cref="IFiniteElement"/>'s in a collection to global stiffness <see cref="Matrix" />.
		/// </summary>
		/// <param name="elements">The <see cref="IFiniteElement" />'s to add to <paramref name="globalStiffness" />.</param>
		/// <inheritdoc cref="AddToGlobalStiffness(IFiniteElement, Matrix{double})" />
		public static void AddToGlobalStiffness([NotNull] this IEnumerable<IFiniteElement> elements, [NotNull] Matrix<double> globalStiffness)
		{
			foreach (var element in elements)
				element.AddToGlobalStiffness(globalStiffness);
		}

		/// <summary>
		///     Add the internal forces of an <see cref="IFiniteElement"/> to the global internal force <see cref="Vector" />.
		/// </summary>
		/// <param name="element">The <see cref="IFiniteElement" /> to add to <paramref name="internalForceVector" />.</param>
		/// <param name="internalForceVector">The global internal force <see cref="Vector" />.</param>
		public static void AddToInternalForces([NotNull] this IFiniteElement element, [NotNull] Vector<double> internalForceVector)
		{
			var dofIndex = element.DoFIndex;

			for (var i = 0; i < dofIndex.Length; i++)
			{
				// DoF index
				var j = dofIndex[i];

				// Add values
				internalForceVector[j] += element.Forces[i];
			}
		}

		/// <summary>
		///     Add the internal forces of a collection of <see cref="IFiniteElement"/>'s to the global internal force <see cref="Vector" />.
		/// </summary>
		/// <param name="elements">The <see cref="IFiniteElement" />'s to add to <paramref name="internalForceVector" />.</param>
		/// <inheritdoc cref="AddToInternalForces" />
		public static void AddToInternalForces([NotNull] this IEnumerable<IFiniteElement> elements, [NotNull] Vector<double> internalForceVector)
		{
			foreach (var element in elements)
				element.AddToInternalForces(internalForceVector);
		}
	}
}