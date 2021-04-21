using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;
using MathNet.Numerics;
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
		/// <param name="globalStiffness">The global stiffness <see cref="Matrix" />.</param>
		/// <param name="element">The <see cref="IFiniteElement" /> to add to <paramref name="globalStiffness" />.</param>
		public static void AddStiffness<TFiniteElement>([NotNull] this Matrix<double> globalStiffness, [NotNull] TFiniteElement element)
			where TFiniteElement : IFiniteElement
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
		///     Add the stiffness of each <see cref="IFiniteElement" />'s in a collection to global stiffness <see cref="Matrix" />.
		/// </summary>
		/// <param name="globalStiffness">The global stiffness <see cref="Matrix" />.</param>
		/// <param name="elements">The collection <see cref="IFiniteElement" />'s to add to <paramref name="globalStiffness" />.</param>
		public static void AddStiffness<TFiniteElement>([NotNull] this Matrix<double> globalStiffness, [NotNull] IEnumerable<TFiniteElement> elements)
			where TFiniteElement : IFiniteElement
		{
			foreach (var element in elements)
				globalStiffness.AddStiffness(element);
		}

		/// <summary>
		///     Add the internal forces of an <see cref="IFiniteElement" /> to the global internal force <see cref="Vector" />.
		/// </summary>
		/// <param name="internalForceVector">The global internal force <see cref="Vector" />.</param>
		/// <param name="element">The <see cref="IFiniteElement" /> to add to <paramref name="internalForceVector" />.</param>
		public static void AddInternalForces<TFiniteElement>([NotNull] this Vector<double> internalForceVector, [NotNull] TFiniteElement element)
			where TFiniteElement : IFiniteElement
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
		/// <param name="internalForceVector">The global internal force <see cref="Vector" />.</param>
		/// <param name="elements">The <see cref="IFiniteElement" />'s to add to <paramref name="internalForceVector" />.</param>
		public static void AddInternalForces<TFiniteElement>([NotNull] this Vector<double> internalForceVector, [NotNull] IEnumerable<TFiniteElement> elements)
			where TFiniteElement : IFiniteElement
		{
			foreach (var element in elements)
				internalForceVector.AddInternalForces(element);
		}

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
			elements.First(element => number == element.Number);

		/// <summary>
		///     Update the displacement <see cref="Vector" /> from this element's grips.
		/// </summary>
		/// <returns>
		///     The displacement <see cref="Vector" />, with components in <see cref="LengthUnit.Millimeter" />.
		/// </returns>
		public static void UpdateDisplacements<TFiniteElement>([NotNull] this TFiniteElement element)
			where TFiniteElement : IFiniteElement
		{
			var displacements = element.Grips
				.SelectMany(g => new[] { g.Displacement.X.Millimeters, g.Displacement.Y.Millimeters })
				.ToVector();

			// Update nonlinear results
			if (element is INonlinearElement nonlinearElement)
			{
				nonlinearElement.LastIterationResult.Displacements    = nonlinearElement.Displacements.Clone();
				nonlinearElement.CurrentIterationResult.Displacements = displacements;
			}
			
			element.Displacements = displacements;
		}

		/// <summary>
		///     Update the displacement <see cref="Vector" /> from this element's grips.
		/// </summary>
		/// <returns>
		///     The displacement <see cref="Vector" />, with components in <see cref="LengthUnit.Millimeter" />.
		/// </returns>
		public static void UpdateDisplacements<TFiniteElement>([NotNull] this IEnumerable<TFiniteElement> elements)
			where TFiniteElement : IFiniteElement
		{
			foreach (var element in elements)
				element.UpdateDisplacements();
		}

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

		/// <summary>
		///		Update the tangent stiffness of this element for the next iteration.
		/// </summary>
		/// <param name="element">The nonlinear element.</param>
		/// <typeparam name="TNonlinearElement">Any type that implements <see cref="INonlinearElement"/>.</typeparam>
		public static void UpdateStiffness<TNonlinearElement>(this TNonlinearElement element)
			where TNonlinearElement : INonlinearElement
		{
			// Get results
			IterationResult
				last    = element.LastIterationResult,
				current = element.CurrentIterationResult;

			// Get current values
			var displacements = current.Displacements;
			var stiffness     = current.Stiffness.Clone();
			
			// Get displacement variation
			var du = displacements - last.Displacements;

			// Get stiffness variation
			var dk = stiffness - last.Stiffness;

			// Increment elements of stiffness matrix
			for (var i = 0; i < stiffness.RowCount; i++)
			for (var j = 0; j < stiffness.ColumnCount; j++)
				stiffness[i, j] += dk.Row(i) / du[j] * displacements;

			// Set new values
			last.Stiffness    = current.Stiffness;
			current.Stiffness = stiffness;
			element.Stiffness = stiffness;
		}

		/// <summary>
		///		Update the tangent stiffness of each element in this collection for the next iteration.
		/// </summary>
		/// <param name="elements">The collection of <see cref="INonlinearElement"/>'s.</param>
		/// <typeparam name="TNonlinearElement">Any type that implements <see cref="INonlinearElement"/>.</typeparam>
		public static void UpdateStiffness<TNonlinearElement>(this IEnumerable<TNonlinearElement> elements)
			where TNonlinearElement : INonlinearElement
		{
			foreach (var element in elements)
				element.UpdateStiffness();
		}
	
		#endregion

	}
}