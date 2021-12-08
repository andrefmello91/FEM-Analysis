using System;
using andrefmello91.Extensions;

namespace andrefmello91.FEMAnalysis;

/// <summary>
///     Monitored displacement struct.
/// </summary>
public struct MonitoredValue : IMonitoredValue<double>, IEquatable<MonitoredValue>, IComparable<MonitoredValue>
{

	#region Properties

	/// <inheritdoc />
	public double LoadFactor { get; }

	/// <inheritdoc />
	public double Value { get; }

	#endregion

	#region Constructors

	/// <summary>
	///     Monitored value constructor.
	/// </summary>
	/// <param name="value">The value, at current <paramref name="loadFactor" />.</param>
	/// <param name="loadFactor">The load factor associated to <paramref name="value" />.</param>
	public MonitoredValue(double value, double loadFactor)
	{
		Value      = value;
		LoadFactor = loadFactor;
	}

	#endregion

	#region Methods

	/// <inheritdoc />
	public override string ToString() =>
		$"Value = {Value:0.00}\n" +
		$"Load Factor = {LoadFactor:0.00}";

	/// <inheritdoc />
	public int CompareTo(MonitoredValue other) =>
		LoadFactor.Approx(other.LoadFactor, 1E-6)
			? 0
			: LoadFactor > other.LoadFactor
				? 1
				: -1;

	/// <inheritdoc />
	public bool Equals(MonitoredValue other) => LoadFactor.Approx(other.LoadFactor, 1E-6) && Value.Approx(other.Value, 1E-6);

	/// <inheritdoc />
	int IComparable<IMonitoredValue<double>>.CompareTo(IMonitoredValue<double>? other) => other is MonitoredValue mv
		? CompareTo(mv)
		: 0;

	/// <inheritdoc />
	bool IEquatable<IMonitoredValue<double>>.Equals(IMonitoredValue<double>? other) => other is MonitoredValue mv && Equals(mv);

	#endregion

}