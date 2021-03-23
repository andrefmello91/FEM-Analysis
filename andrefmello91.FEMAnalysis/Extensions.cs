using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using andrefmello91.OnPlaneComponents;
using andrefmello91.OnPlaneComponents.Displacement;
using andrefmello91.OnPlaneComponents.Force;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using UnitsNet.Units;
#nullable enable

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Extensions for <see cref="IGrip" /> and <see cref="IFiniteElement" />.
	/// </summary>
	public static class Extensions
	{
		#region Methods

		/// <summary>
		///     Add the stiffness of an <see cref="IFiniteElement" /> to global stiffness <see cref="Matrix" />.
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
		///     Add the stiffness of each <see cref="IFiniteElement" />'s in a collection to global stiffness <see cref="Matrix" />
		///     .
		/// </summary>
		/// <param name="elements">The <see cref="IFiniteElement" />'s to add to <paramref name="globalStiffness" />.</param>
		/// <inheritdoc cref="AddToGlobalStiffness(IFiniteElement, Matrix{double})" />
		public static void AddToGlobalStiffness([NotNull] this IEnumerable<IFiniteElement> elements, [NotNull] Matrix<double> globalStiffness)
		{
			foreach (var element in elements)
				element.AddToGlobalStiffness(globalStiffness);
		}

		/// <summary>
		///     Add the internal forces of an <see cref="IFiniteElement" /> to the global internal force <see cref="Vector" />.
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
		///     Add the internal forces of a collection of <see cref="IFiniteElement" />'s to the global internal force
		///     <see cref="Vector" />.
		/// </summary>
		/// <param name="elements">The <see cref="IFiniteElement" />'s to add to <paramref name="internalForceVector" />.</param>
		/// <inheritdoc cref="AddToInternalForces(IFiniteElement, Vector{double})" />
		public static void AddToInternalForces([NotNull] this IEnumerable<IFiniteElement> elements, [NotNull] Vector<double> internalForceVector)
		{
			foreach (var element in elements)
				element.AddToInternalForces(internalForceVector);
		}

		/// <summary>
		///     Calculate forces in each element in this collection after updating displacements in grips.
		/// </summary>
		/// <inheritdoc cref="IFiniteElement.CalculateForces" />
		public static void CalculateForces(this IEnumerable<IFiniteElement> elements)
		{
			foreach (var element in elements)
				element.CalculateForces();
		}

		/// <summary>
		///     Return an <see cref="INumberedElement" /> of a collection, in given <paramref name="number" />.
		/// </summary>
		/// <param name="elements">The collection of <see cref="INumberedElement" />'s.</param>
		/// <param name="number">The number of the element wanted.</param>
		public static INumberedElement? GetByNumber(this IEnumerable<INumberedElement> elements, int number) => elements.First(element => number == element.Number);

		/// <summary>
		///     Get global indexes of the degrees of freedom af a collection of grips.
		/// </summary>
		/// <param name="grips">A collection of <see cref="IGrip" />'s.</param>
		public static IEnumerable<int> GlobalIndexes(IEnumerable<IGrip> grips) => grips.SelectMany(GlobalIndexes);

		/// <summary>
		///     Get global indexes of the degrees of freedom af a grip.
		/// </summary>
		/// <param name="grip">An <see cref="IGrip" />.</param>
		public static IEnumerable<int> GlobalIndexes(IGrip grip)
		{
			var n = 2 * grip.Number;

			yield return n - 2;
			yield return n - 1;
		}

		/// <summary>
		///     Set displacements to an <see cref="IGrip" /> from the global displacement <see cref="Vector" />.
		/// </summary>
		/// <param name="grip">The <see cref="IGrip" /> to set displacements.</param>
		/// <param name="globalDisplacementVector">
		///     The global displacement vector, with components in
		///     <see cref="LengthUnit.Millimeter" />.
		/// </param>
		public static void SetDisplacements([NotNull] this IGrip grip, [NotNull] Vector<double> globalDisplacementVector)
		{
			var x = globalDisplacementVector[grip.DoFIndex[0]];
			var y = globalDisplacementVector[grip.DoFIndex[1]];

			grip.Displacement = new PlaneDisplacement(x, y);
		}

		/// <summary>
		///     Set displacements to a collection of <see cref="IGrip" />'s from the global displacement <see cref="Vector" />.
		/// </summary>
		/// <param name="grips">The collection of <see cref="IGrip" />'s to set displacements.</param>
		/// <inheritdoc cref="SetDisplacements(IGrip, Vector{double})" />
		public static void SetDisplacements([NotNull] this IEnumerable<IGrip> grips, [NotNull] Vector<double> globalDisplacementVector)
		{
			foreach (var grip in grips)
				grip.SetDisplacements(globalDisplacementVector);
		}

		/// <summary>
		///     Set support reactions to an <see cref="IGrip" /> from the global reaction <see cref="Vector" />.
		/// </summary>
		/// <param name="grip">The <see cref="IGrip" /> to set support reactions.</param>
		/// <param name="globalReactionVector">
		///     The global reaction <see cref="Vector" />, with components in
		///     <see cref="ForceUnit.Newton" />.
		/// </param>
		public static void SetReactions([NotNull] this IGrip grip, [NotNull] Vector<double> globalReactionVector)
		{
			// Verify if grip is free
			if (grip.Constraint.Direction is ComponentDirection.None)
				return;

			var x = globalReactionVector[grip.DoFIndex[0]];
			var y = globalReactionVector[grip.DoFIndex[1]];

			grip.Reaction = new PlaneForce(x, y);
		}

		/// <summary>
		///     Set support reactions to a collection of <see cref="IGrip" />'s from the global reaction <see cref="Vector" />.
		/// </summary>
		/// <param name="grips">The collection of <see cref="IGrip" />'s to set support reactions.</param>
		/// <inheritdoc cref="SetReactions(IGrip, Vector{double})" />
		public static void SetReactions([NotNull] this IEnumerable<IGrip> grips, [NotNull] Vector<double> globalReactionVector)
		{
			foreach (var grip in grips)
				grip.SetReactions(globalReactionVector);
		}

		#endregion
	}
}