using System;
using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;
using UnitsNet;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Monitored displacement struct.
	/// </summary>
	public struct MonitoredDisplacement : IEquatable<MonitoredDisplacement>, IComparable<MonitoredDisplacement>
	{
		/// <summary>
		///     Get the displacement value, at current <see cref="LoadFactor" />.
		/// </summary>
		public Length Displacement { get; }

		/// <summary>
		///     Get the load factor associated to <see cref="Displacement" />.
		/// </summary>
		public double LoadFactor { get; }

		/// <summary>
		///     Monitored displacement constructor.
		/// </summary>
		/// <param name="displacement">The displacement value, at current <paramref name="loadFactor" />.</param>
		/// <param name="loadFactor">The load factor associated to <paramref name="displacement" />.</param>
		public MonitoredDisplacement(Length displacement, double loadFactor)
		{
			Displacement = displacement;
			LoadFactor   = loadFactor;
		}

		/// <inheritdoc />
		public bool Equals(MonitoredDisplacement other) => LoadFactor.Approx(other.LoadFactor, 1E-6) && Displacement.Approx(other.Displacement, PlaneDisplacement.Tolerance);

		/// <inheritdoc />
		public int CompareTo(MonitoredDisplacement other) =>
			LoadFactor.Approx(other.LoadFactor, 1E-6)
				? 0
				: LoadFactor > other.LoadFactor
					? 1
					: -1;

		/// <inheritdoc />
		public override string ToString() =>
			$"Displacement = {Displacement}\n" +
			$"Load Factor = {LoadFactor:0.00}";
	}
}