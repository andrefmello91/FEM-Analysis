using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using andrefmello91.OnPlaneComponents;
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
		///     Calculate forces in each element in this collection after updating displacements in grips.
		/// </summary>
		/// <inheritdoc cref="IFiniteElement.CalculateForces" />
		public static void CalculateForces<TFiniteElement>(this IEnumerable<TFiniteElement> elements)
			where TFiniteElement : IFiniteElement
		{
			foreach (var element in elements)
				element.CalculateForces();
		}

		/// <summary>
		///     Return an element of a collection, in given <paramref name="number" />.
		/// </summary>
		/// <param name="elements">The collection of <see cref="INumberedElement" />'s.</param>
		/// <param name="number">The number of the element wanted.</param>
		public static TNumberedElement? GetByNumber<TNumberedElement>(this IEnumerable<TNumberedElement> elements, int number)
			where TNumberedElement : INumberedElement =>
			elements.FirstOrDefault(element => number == element.Number);

		/// <summary>
		///     Get the indexes of constrained degrees of freedom from a collection of grips.
		/// </summary>
		/// <param name="grips">The collection containing all distinct grips of the finite element model.</param>
		public static IEnumerable<int> GetConstraintIndex(this IEnumerable<IGrip> grips)
		{
			foreach (var grip in grips)
			{
				// Get DoF indexes
				var index = grip.DoFIndex;
				int
					i = index[0],
					j = index[1];

				var constraint = grip.Constraint;

				if (constraint.X)

					// There is a support in X direction
					yield return i;

				if (constraint.Y)

					// There is a support in Y direction
					yield return j;
			}
		}

		/// <summary>
		///     Get the displacement <see cref="Vector" /> from this element's grips.
		/// </summary>
		/// <returns>
		///     The displacement <see cref="Vector" />, with components in <see cref="LengthUnit.Millimeter" />.
		/// </returns>
		public static DisplacementVector GetDisplacementsFromGrips([NotNull] this IFiniteElement element) =>
			new (element.Grips
				.SelectMany(g => new[] { g.Displacement.X, g.Displacement.Y }));

		/// <summary>
		///     Get global indexes of the degrees of freedom of a collection of grips.
		/// </summary>
		/// <param name="grips">A collection of <see cref="IGrip" />'s.</param>
		public static IEnumerable<int> GlobalIndexes(params IGrip[] grips)
		{
			foreach (var grip in grips)
			{
				var n = 2 * grip.Number;

				yield return n - 2;
				yield return n - 1;
			}
		}

		/// <summary>
		///     Set displacements to an <see cref="IGrip" /> from the global displacement <see cref="Vector" />.
		/// </summary>
		/// <param name="grip">The <see cref="IGrip" /> to set displacements.</param>
		/// <param name="globalDisplacementVector">The global displacement vector.</param>
		public static void SetDisplacements([NotNull] this IGrip grip, [NotNull] DisplacementVector globalDisplacementVector)
		{
			var x = globalDisplacementVector[grip.DoFIndex[0]];
			var y = globalDisplacementVector[grip.DoFIndex[1]];

			grip.Displacement = new PlaneDisplacement(x, y);
		}

		/// <summary>
		///     Set displacements to a collection of <see cref="IGrip" />'s from the global displacement <see cref="Vector" />.
		/// </summary>
		/// <param name="grips">The collection of <see cref="IGrip" />'s to set displacements.</param>
		/// <param name="globalDisplacementVector">The global displacement vector.</param>
		public static void SetDisplacements([NotNull] this IEnumerable<IGrip> grips, [NotNull] DisplacementVector globalDisplacementVector)
		{
			foreach (var grip in grips)
				grip.SetDisplacements(globalDisplacementVector);
		}

		/// <summary>
		///     Set support reactions to an <see cref="IGrip" /> from the global reaction <see cref="Vector" />.
		/// </summary>
		/// <param name="grip">The <see cref="IGrip" /> to set support reactions.</param>
		/// <param name="globalReactionVector">The global reaction vector.</param>
		public static void SetReactions([NotNull] this IGrip grip, [NotNull] ForceVector globalReactionVector)
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
		/// <param name="globalReactionVector">The global reaction vector.</param>
		public static void SetReactions([NotNull] this IEnumerable<IGrip> grips, [NotNull] ForceVector globalReactionVector)
		{
			foreach (var grip in grips)
				grip.SetReactions(globalReactionVector);
		}

		/// <summary>
		///     Update the displacement vector for this element.
		/// </summary>
		/// <inheritdoc cref="UpdateStiffness" />
		public static void UpdateDisplacements([NotNull] this IEnumerable<IFiniteElement> elements)
		{
			foreach (var element in elements)
				element.UpdateDisplacements();
		}


		/// <summary>
		///     Update the stiffness of each element in this collection.
		/// </summary>
		/// <param name="elements">The collection of finite elements.</param>
		public static void UpdateStiffness(this IEnumerable<IFiniteElement> elements)
		{
			foreach (var element in elements)
				element.UpdateStiffness();
		}

		#endregion

	}
}