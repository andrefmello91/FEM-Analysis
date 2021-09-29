using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.Data.Text;
using MathNet.Numerics.LinearAlgebra;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Element monitor class.
	/// </summary>
	public class ElementMonitor : IEnumerable<IMonitoredValue<double>>
	{

		#region Fields

		private readonly Func<INumberedElement, double> _function;
		private readonly List<IMonitoredValue<double>> _values = new();

		#endregion

		#region Properties

		/// <summary>
		///     The name of this monitor.
		/// </summary>
		public string Name { get; }

		#endregion

		#region Constructors

		/// <summary>
		///     Create a grip monitor
		/// </summary>
		/// <param name="function">The function to get the monitored value.</param>
		/// <param name="name">The name of this monitor.</param>
		public ElementMonitor(Func<INumberedElement, double> function, string name)
		{
			_function = function;
			Name      = name;
		}

		#endregion

		#region Methods

		/// <summary>
		///     Add the monitored value for the current load factor.
		/// </summary>
		/// <param name="loadFactor">The current load factor.</param>
		/// <param name="element">The element.</param>
		public void AddMonitoredValue(double loadFactor, INumberedElement element) => _values.Add(new MonitoredValue(_function(element), loadFactor));

		/// <summary>
		///     Export data to a CSV file.
		/// </summary>
		/// <param name="outputPath">The output path.</param>
		/// <param name="delimiter">CSV delimiter.</param>
		public void Export(string outputPath, string delimiter = ";")
		{
			if (!_values.Any())
				return;

			var lfs = _values
				.Select(m => m.LoadFactor)
				.ToArray();

			var m = _values
				.Select(m => m.Value)
				.ToArray();

			// Create a matrix
			var result = Matrix<double>.Build.DenseOfColumnArrays(lfs, m);

			// Create headers
			var headers = new[] { "Load Factor", "X", "Y" };

			// Set full save location
			var fullPath = $"{outputPath.TrimEnd('\u002F', '\u005C')}/{Name}.csv";

			// Export
			DelimitedWriter.Write(fullPath, result, delimiter, headers);
		}

		/// <inheritdoc />
		public IEnumerator<IMonitoredValue<double>> GetEnumerator() => _values.GetEnumerator();

		/// <inheritdoc />
		IEnumerator IEnumerable.GetEnumerator() => GetEnumerator();

		#endregion

	}
}