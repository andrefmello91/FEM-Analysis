using andrefmello91.Extensions;
using andrefmello91.OnPlaneComponents;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet.Units;

namespace andrefmello91.FEMAnalysis;

/// <summary>
///     Grip monitor class.
/// </summary>
public class GripMonitor : ElementMonitor
{

	#region Fields

	private readonly LengthUnit _unit;

	#endregion

	#region Constructors

	/// <inheritdoc />
	public GripMonitor(string name, LengthUnit displacementUnit = LengthUnit.Millimeter)
		: base(name, Label(displacementUnit)) =>
		_unit = displacementUnit;

	#endregion

	#region Methods

	private static string[] Label(LengthUnit displacementUnit)
	{
		var d = displacementUnit.Abbrev();

		return new[]
		{
			"Load Factor",
			$"Ux ({d})",
			$"Uy ({d})"
		};
	}

	/// <inheritdoc />
	public override void AddMonitoredValue(double loadFactor, INumberedElement element)
	{
		if (element is not IGrip grip)
			return;

		Values.Add(new MonitoredValue(loadFactor, grip.Displacement, _unit));
	}

	#endregion

	private class MonitoredValue : IVectorTransformable
	{

		#region Properties

		public double LoadFactor { get; }

		public double Ux { get; }

		public double Uy { get; }

		#endregion

		#region Constructors

		public MonitoredValue(double loadFactor, PlaneDisplacement displacement, LengthUnit unit = LengthUnit.Millimeter)
		{
			LoadFactor = loadFactor;
			Ux         = displacement.X.As(unit);
			Uy         = displacement.Y.As(unit);
		}

		#endregion

		#region Methods

		/// <inheritdoc />
		public Vector<double> AsVector() => new[]
		{
			LoadFactor,
			Ux,
			Uy
		}.ToVector();

		#endregion

	}
}