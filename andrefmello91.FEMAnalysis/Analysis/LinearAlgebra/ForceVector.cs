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
	///		Force vector class.
	/// </summary>
	/// <remarks>
	///		Unit is <see cref="ForceUnit"/>.
	///		<para>
	///		Quantity is <see cref="Force"/>.
	///		</para>
	/// </remarks>
	public class ForceVector : ComponentVector<Force, ForceUnit>
	{

		/// <inheritdoc />
		public ForceVector(IEnumerable<double> value, ForceUnit unit = ForceUnit.Newton)
			: base(value, unit)
		{
		}

		/// <inheritdoc />
		public ForceVector(IEnumerable<Force> value)
			: base(value)
		{
		}
	}
}