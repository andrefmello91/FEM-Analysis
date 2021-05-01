using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.Data.Text;
using MathNet.Numerics.LinearAlgebra;
using UnitsNet.Units;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Output data class.
	/// </summary>
	public class FEMOutput
	{

		#region Properties

		/// <summary>
		///     Values of calcultated load step results.
		/// </summary>
		public List<LoadStepResult> LoadStepResults { get; }

		/// <summary>
		///     Values of monitored displacements.
		/// </summary>
		public List<MonitoredDisplacement> MonitoredDisplacements { get; }

		#endregion

		#region Constructors

		/// <summary>
		///     Output data constructor.
		/// </summary>
		/// <param name="loadStepResults">The values of load step results.</param>
		public FEMOutput([NotNull] IEnumerable<LoadStepResult> loadStepResults)
		{
			LoadStepResults = loadStepResults
				.Where(ls => ls.IsCalculated)
				.ToList();

			MonitoredDisplacements = LoadStepResults
				.Where(ls => ls.MonitoredDisplacement.HasValue)
				.Select(ls => ls.MonitoredDisplacement!.Value)
				.ToList();
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
			var disps = MonitoredDisplacements
				.Select(m => m.Displacement.ToUnit(unit).Value)
				.ToVector();

			var lfs = MonitoredDisplacements
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
		}

		#endregion

	}
}