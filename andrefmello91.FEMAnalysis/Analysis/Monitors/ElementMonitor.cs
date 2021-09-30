using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using andrefmello91.Extensions;
using MathNet.Numerics.Data.Text;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis
{

	/// <summary>
	///     Element monitor class.
	/// </summary>
	public abstract class ElementMonitor
	{

		#region Fields

		protected readonly List<IVectorTransformable> Values = new();

		#endregion

		#region Properties

		/// <summary>
		///     The name of this monitor.
		/// </summary>
		public string Name { get; }

		/// <summary>
		///     The labels for the output file.
		/// </summary>
		protected string[] Labels { get; }

		#endregion

		#region Constructors

		/// <summary>
		///     Create an element monitor.
		/// </summary>
		/// <param name="name">The name of this monitor.</param>
		/// <param name="labels">The labels for the output file.</param>
		protected ElementMonitor(string name, string[] labels)
		{
			Name   = name;
			Labels = labels;
		}

		#endregion

		#region Methods

		/// <summary>
		///     Add the monitored value for the current load factor.
		/// </summary>
		/// <param name="loadFactor">The current load factor.</param>
		/// <param name="element">The element.</param>
		public abstract void AddMonitoredValue(double loadFactor, INumberedElement element);

		/// <summary>
		///     Export data to a CSV file.
		/// </summary>
		/// <param name="outputPath">The output path.</param>
		/// <param name="fileName">The filename, without extension.</param>
		/// <param name="delimiter">CSV delimiter.</param>
		public void Export(string outputPath, string fileName = "FEM_Output", string delimiter = ";")
		{
			if (!Values.Any())
				return;

			// Create a matrix
			var result = Matrix<double>.Build.DenseOfRowVectors(Values.Select(v => v.AsVector()));

			// Set full save location
			var fullPath = $"{outputPath.TrimEnd('\u002F', '\u005C')}/{fileName}_{Name}.csv";

			// Export
			DelimitedWriter.Write(fullPath, result, delimiter, Labels, formatProvider: CultureInfo.CurrentCulture);
		}

		#endregion

	}
}