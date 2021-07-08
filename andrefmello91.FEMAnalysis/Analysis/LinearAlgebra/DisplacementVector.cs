﻿using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;
using UnitsNet.Units;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Displacement vector class.
	/// </summary>
	/// <remarks>
	///     Unit is <see cref="LengthUnit" />.
	///     <para>
	///         Quantity is <see cref="Length" />.
	///     </para>
	/// </remarks>
	public class DisplacementVector : ComponentVector<Length, LengthUnit>
	{

		#region Constructors

		/// <inheritdoc />
		public DisplacementVector(IEnumerable<double> values, LengthUnit unit = LengthUnit.Millimeter)
			: base(values, unit)
		{
		}

		/// <inheritdoc />
		public DisplacementVector(IEnumerable<Length> values)
			: base(values)
		{
		}

		#endregion

		#region Methods

		/// <summary>
		///     Assemble the element displacement vector from it's grips.
		/// </summary>
		/// <param name="element">The finite element.</param>
		public static DisplacementVector Assemble(IFiniteElement element) =>
			new(element.Grips
				.SelectMany(g => new[] { g.Displacement.X, g.Displacement.Y }));

		/// <summary>
		///     Assemble the global displacement vector.
		/// </summary>
		/// <param name="femInput">Finite element input.</param>
		public static DisplacementVector Assemble(IFEMInput femInput)
		{
			// Initialize the force vector
			var d = Zero(femInput.NumberOfDoFs);
			d.ConstraintIndex = femInput.ConstraintIndex;

			// Read the nodes data
			foreach (var grip in femInput.Grips)
			{
				// Get DoF indexes
				var index = grip.DoFIndex;
				int
					i = index[0],
					j = index[1];

				// Set to force vector
				d[i] = grip.Displacement.X;
				d[j] = grip.Displacement.Y;
			}

			return d;
		}

		/// <summary>
		///     Create a displacement vector with zero elements.
		/// </summary>
		/// <param name="size">The size of the vector.</param>
		public static DisplacementVector Zero(int size) => new(new double[size]);

		/// <inheritdoc cref="ICloneable.Clone" />
		public new DisplacementVector Clone() => (DisplacementVector) base.Clone();

		#endregion

		#region Operators

		/// <inheritdoc cref="ComponentVector{TQuantity,TUnit}.op_Addition" />
		public static DisplacementVector operator +(DisplacementVector left, DisplacementVector right) => (DisplacementVector) ((ComponentVector<Length, LengthUnit>) left + right);

		/// <inheritdoc cref="ComponentVector{TQuantity,TUnit}.op_Subtraction" />
		public static DisplacementVector operator -(DisplacementVector left, DisplacementVector right) => (DisplacementVector) ((ComponentVector<Length, LengthUnit>) left - right);

		/// <inheritdoc cref="ComponentVector{TQuantity,TUnit}.op_Multiply(double,ComponentVector{TQuantity,TUnit}) " />
		public static DisplacementVector operator *(double multiplier, DisplacementVector vector) => (DisplacementVector) (multiplier * (ComponentVector<Length, LengthUnit>) vector);

		/// <inheritdoc cref="ComponentVector{TQuantity,TUnit}.op_Multiply(double,ComponentVector{TQuantity,TUnit}) " />
		public static DisplacementVector operator *(DisplacementVector vector, double multiplier) => multiplier * vector;

		/// <inheritdoc cref="Vector{T}.op_UnaryNegation" />
		public static DisplacementVector operator -(DisplacementVector vector) => (DisplacementVector) (-(ComponentVector<Length, LengthUnit>) vector);

		/// <inheritdoc cref="Vector{T}.op_Division(Vector{T}, T)" />
		public static DisplacementVector operator /(DisplacementVector vector, double divisor) => (DisplacementVector) ((ComponentVector<Length, LengthUnit>) vector / divisor);

		#endregion

	}
}