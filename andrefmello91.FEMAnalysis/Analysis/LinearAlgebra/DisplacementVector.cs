using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet;
using UnitsNet.Units;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///		Displacement vector class.
	/// </summary>
	/// <remarks>
	///		Unit is <see cref="LengthUnit"/>.
	///		<para>
	///		Quantity is <see cref="Length"/>.
	///		</para>
	/// </remarks>
	public class DisplacementVector : ComponentVector<Length, LengthUnit>
	{

		/// <inheritdoc />
		public DisplacementVector(IEnumerable<double> value, LengthUnit unit = LengthUnit.Millimeter)
			: base(value, unit)
		{
		}

		/// <inheritdoc />
		public DisplacementVector(IEnumerable<Length> value)
			: base(value)
		{
		}
	}

}