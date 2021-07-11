using System;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;
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
	public class DisplacementVector : QuantityVector<Length, LengthUnit>
	{

		#region Properties

		/// <summary>
		///     Default tolerance for component vector.
		/// </summary>
		private static Length Tolerance { get; } = PlaneDisplacement.Tolerance;

		#endregion

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
		public new static DisplacementVector Zero(int size) => new(new double[size]);

		/// <inheritdoc cref="ICloneable.Clone" />
		public override QuantityVector<Length, LengthUnit> Clone() => new DisplacementVector(Values, Unit);

		/// <inheritdoc />
		public override QuantityVector<Length, LengthUnit> Convert(LengthUnit unit) => new DisplacementVector(Values.GetQuantities<Length, LengthUnit>(Unit).GetValues(unit), unit);

		/// <inheritdoc />
		public override QuantityVector<Length, LengthUnit> Simplified(IEnumerable<int>? indexes = null) => Simplified(Tolerance, indexes);

		/// <inheritdoc />
		public override QuantityVector<Length, LengthUnit> Simplified(double? threshold, IEnumerable<int>? indexes = null)
		{
			var simplified = Values.ToVector();

			if (indexes is not null)
				foreach (var index in indexes)
					simplified[index] = 0;

			if (threshold.HasValue)
				simplified.CoerceZero(threshold.Value);

			return
				new DisplacementVector(simplified, Unit);
		}

		#region Object override

		/// <inheritdoc cref="ForceVector.op_Addition" />
		public static DisplacementVector operator +(DisplacementVector left, DisplacementVector right)
		{
			right = right.Unit == left.Unit
				? right
				: (DisplacementVector) right.Convert(left.Unit);

			var vec = (Vector<double>) left + right;

			return
				new DisplacementVector(vec, left.Unit);
		}

		/// <inheritdoc cref="Vector{T}.op_Division(Vector{T}, T)" />
		public static DisplacementVector operator /(DisplacementVector vector, double divisor) => new(vector.Values.Select(v => v / divisor), vector.Unit);

		/// <inheritdoc cref="ForceVector.op_Multiply(double,ForceVector) " />
		public static DisplacementVector operator *(double multiplier, DisplacementVector vector) =>
			new(vector.Values.Select(v => v * multiplier), vector.Unit);

		/// <inheritdoc cref="ForceVector.op_Multiply(double,ForceVector) " />
		public static DisplacementVector operator *(DisplacementVector vector, double multiplier) => multiplier * vector;

		/// <inheritdoc cref="Matrix{T}.op_Multiply(Matrix{T}, Vector{T})" />
		public static DisplacementVector operator *(Matrix<double> left, DisplacementVector right) => new(left * (Vector<double>) right, right.Unit);

		/// <inheritdoc cref="ForceVector.op_Subtraction" />
		public static DisplacementVector operator -(DisplacementVector left, DisplacementVector right)
		{
			right = right.Unit == left.Unit
				? right
				: (DisplacementVector) right.Convert(left.Unit);

			var vec = (Vector<double>) left - right;

			return
				new DisplacementVector(vec, left.Unit);
		}

		/// <inheritdoc cref="Vector{T}.op_UnaryNegation" />
		public static DisplacementVector operator -(DisplacementVector vector) => new(vector.Values.Select(v => -v), vector.Unit);

		#endregion

		#endregion

	}
}