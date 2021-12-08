using System.Collections.Generic;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.Data.Text;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet.Units;

namespace andrefmello91.FEMAnalysis;

/// <summary>
///     Output data class.
/// </summary>
public class FEMOutput : List<MonitoredDisplacement>
{

	#region Fields

	private readonly List<ElementMonitor> _monitors;

	#endregion

	#region Properties

	/// <summary>
	///     Values of calculated step results.
	/// </summary>
	public List<LoadStep> LoadStepResults { get; }

	#endregion

	#region Constructors

	/// <summary>
	///     Output data constructor.
	/// </summary>
	/// <param name="loadStepResults">The values of step results.</param>
	/// <param name="elementMonitors">The element monitors.</param>
	public FEMOutput(IEnumerable<LoadStep> loadStepResults, IEnumerable<ElementMonitor> elementMonitors)
		: base(loadStepResults.Where(ls => ls.MonitoredDisplacement.HasValue).Select(ls => ls.MonitoredDisplacement!.Value))
	{
		LoadStepResults = loadStepResults
			.Where(ls => ls.Converged)
			.ToList();

		_monitors = elementMonitors.ToList();
	}

	#endregion

	#region Methods

	/// <summary>
	///     Export output data to a csv file.
	/// </summary>
	/// <param name="outputPath">The output file save location.</param>
	/// <param name="fileName">The filename, without extension.</param>
	/// <param name="unit">The required <see cref="LengthUnit" /> of displacements.</param>
	/// <param name="delimiter">The delimiter for csv file.</param>
	public void Export(string outputPath, string fileName = "FEM_Output", LengthUnit unit = LengthUnit.Millimeter, string delimiter = ";")
	{
		// Get displacements and load factors as vectors
		var disps = this
			.Select(m => m.Displacement.ToUnit(unit).Value)
			.ToVector();

		var lfs = this
			.Select(m => m.LoadFactor)
			.ToVector();

		// Create a matrix
		var result = Matrix<double>.Build.DenseOfColumnVectors(lfs, disps);

		// Create headers
		var headers = new[] { "Load Factor", $"Displacement ({unit.Abbrev()})" };

		// Set full save location
		var fullPath = $"{outputPath.TrimEnd('\u002F', '\u005C')}/{fileName}.csv";

		// Export
		DelimitedWriter.Write(fullPath, result, delimiter, headers);

		// Export monitors
		foreach (var monitor in _monitors)
			monitor.Export(outputPath, fileName, delimiter);
	}

	#endregion

}