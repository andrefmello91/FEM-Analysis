using System;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
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

		#if NET5_0
		
		/// <inheritdoc />
		public override DisplacementVector Convert(LengthUnit unit) => new (Values.GetQuantities<Length, LengthUnit>(Unit).GetValues(unit), unit);

		/// <inheritdoc cref="ICloneable.Clone" />
		public override DisplacementVector Clone() => new (Values, Unit);
		
		#else
		
		/// <inheritdoc />
		public override ComponentVector<Length, LengthUnit> Convert(LengthUnit unit) => new DisplacementVector(Values.GetQuantities<Length, LengthUnit>(Unit).GetValues(unit), unit);

		/// <inheritdoc cref="ICloneable.Clone" />
		public override ComponentVector<Length, LengthUnit> Clone() => new DisplacementVector(Values, Unit);
		
		#endif
		
		#endregion

		#region Operators

		/// <inheritdoc cref="ForceVector.op_Addition" />
		public static DisplacementVector operator +(DisplacementVector left, DisplacementVector right) =>
			new(left.ToVector(left.Unit) + right.ToVector(left.Unit), left.Unit)
			{
				ConstraintIndex = left.ConstraintIndex ?? right.ConstraintIndex
			};
		
		/// <inheritdoc cref="ForceVector.op_Subtraction" />
		public static DisplacementVector operator -(DisplacementVector left, DisplacementVector right) => 
			new(left.ToVector(left.Unit) - right.ToVector(left.Unit), left.Unit)
			{
				ConstraintIndex = left.ConstraintIndex ?? right.ConstraintIndex
			};
		
		/// <inheritdoc cref="ForceVector.op_Multiply(double,ForceVector) " />
		public static DisplacementVector operator *(double multiplier, DisplacementVector vector) =>
			new(vector.Values.Select(v => v * multiplier), vector.Unit)
			{
				ConstraintIndex = vector.ConstraintIndex
			};

		/// <inheritdoc cref="ForceVector.op_Multiply(double,ForceVector) " />
		public static DisplacementVector operator *(DisplacementVector vector, double multiplier) => multiplier * vector;

		/// <inheritdoc cref="Vector{T}.op_UnaryNegation" />
		public static DisplacementVector operator -(DisplacementVector vector) =>
			new(vector.Select(v => -v.Value), vector.Unit)
			{
				ConstraintIndex = vector.ConstraintIndex
			};

		/// <inheritdoc cref="Vector{T}.op_Division(Vector{T}, T)" />
		public static DisplacementVector operator /(DisplacementVector vector, double divisor) =>
			new(vector.Select(v => v.Value / divisor), vector.Unit)
			{
				ConstraintIndex = vector.ConstraintIndex
			};

		#endregion

	}
}