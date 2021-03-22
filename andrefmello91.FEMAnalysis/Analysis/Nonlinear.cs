using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace andrefmello91.FEMAnalysis
{
	/// <summary>
	///     Nonlinear analysis class.
	/// </summary>
	public class NonlinearAnalysis : Analysis
	{
		#region Fields

		/// <summary>
		///     The displacement <see cref="Vector" /> of current iteration
		/// </summary>
		private Vector<double> _currentDisplacements;

		/// <summary>
		///     Field to store each iteration force <see cref="Vector" />
		/// </summary>
		private Vector<double> _currentForces;

		/// <summary>
		///     The residual force <see cref="Vector" /> of current iteration
		/// </summary>
		private Vector<double> _currentResidual;

		/// <summary>
		///     Field to store the DoF index for <see cref="MonitoredDisplacements" />.
		/// </summary>
		private int? _monitoredIndex;

		#endregion

		#region Properties

		/// <summary>
		///     Get/set values of displacements monitored.
		/// </summary>
		public List<double> MonitoredDisplacements { get; private set; }

		/// <summary>
		///     Get/set values of load factor associated to <see cref="MonitoredDisplacements" />.
		/// </summary>
		public List<double> MonitoredLoadFactor { get; private set; }

		/// <summary>
		///     Get/set when to stop analysis.
		/// </summary>
		public bool Stop { get; private set; }

		/// <summary>
		///     Get/set the stop message.
		/// </summary>
		public string StopMessage { get; private set; }

		#endregion

		#region Constructors

		/// <summary>
		///     Nonlinear analysis object
		/// </summary>
		/// <inheritdoc />
		public NonlinearAnalysis(InputData inputData)
			: base(inputData)
		{
		}

		#endregion

		#region  Methods

		/// <summary>
		///     Do the analysis.
		/// </summary>
		/// <param name="loadFactor">The load factor to multiply <see cref="Analysis.ForceVector" /> (default: 1).</param>
		/// <param name="monitoredIndex">The DoF index to monitor, if wanted.</param>
		/// <param name="numLoadSteps">The number of load steps to perform (default: 50).</param>
		/// <param name="tolerance">The convergence tolerance (default: 1E-3).</param>
		/// <param name="maxIterations">Maximum number of iterations for each load step (default: 1000).</param>
		public void Do(double loadFactor = 1, int? monitoredIndex = null, int numLoadSteps = 50, double tolerance = 1E-6, int maxIterations = 10000)
		{
			// Initiate lists
			Initiate(monitoredIndex);

			// Get force vector
			ForceVector = InputData.ForceVector * loadFactor;

			// Get the initial stiffness and force vector simplified
			UpdateStiffness();

			// Analysis by load steps
			StepAnalysis(numLoadSteps, tolerance, maxIterations);

			// Set displacements
			DisplacementVector = _currentDisplacements;
			NodalDisplacements(_currentDisplacements);
		}

		/// <summary>
		///     Initiate fields.
		/// </summary>
		/// <param name="monitoredIndex">The DoF index to monitor, if wanted.</param>
		private void Initiate(int? monitoredIndex)
		{
			Stop        = false;
			StopMessage = string.Empty;

			_monitoredIndex = monitoredIndex;

			MonitoredDisplacements = _monitoredIndex.HasValue ? new List<double>() : null;
			MonitoredLoadFactor    = _monitoredIndex.HasValue ? new List<double>() : null;
		}

		/// <summary>
		///     Do analysis by load steps.
		/// </summary>
		/// <param name="numLoadSteps">The number of load steps to perform (default: 50).</param>
		/// <param name="tolerance">The convergence tolerance (default: 1E-3).</param>
		/// <param name="maxIterations">Maximum number of iterations for each load step (default: 1000).</param>
		private void StepAnalysis(int numLoadSteps, double tolerance, int maxIterations)
		{
			// Solve the initial displacements
			var lf0 = (double) 1 / numLoadSteps;

			_currentDisplacements = CalculateDisplacements(GlobalStiffness, lf0 * ForceVector);

			for (var ls = 1; ls <= numLoadSteps; ls++)
			{
				// Calculate the current load factor
				var lf = (double) ls / numLoadSteps;

				// Get the force vector
				_currentForces = lf * ForceVector;

				// Iterate
				Iterate(ls, tolerance, maxIterations);

				// Verify if convergence was not reached
				if (Stop)
					break;

				// Set load step results
				SaveLoadStepResults(lf);

				// Update stiffness
				UpdateStiffness();
			}
		}

		/// <summary>
		///     Iterate to find solution.
		/// </summary>
		/// <param name="loadStep">Current load step.</param>
		/// <param name="tolerance">The convergence tolerance (default: 1E-3).</param>
		/// <param name="maxIterations">Maximum number of iterations for each load step (default: 1000).</param>
		private void Iterate(int loadStep, double tolerance, int maxIterations)
		{
			for (var it = 1; it <= maxIterations; it++)
			{
				// Calculate element displacements and forces
				ElementAnalysis(_currentDisplacements);

				// Calculate residual
				_currentResidual = ResidualForces();

				// Check convergence
				if (ConvergenceReached(tolerance, it))
					break;

				// Check if maximum number of iterations is reached
				if (it == maxIterations)
				{
					Stop = true;
					StopMessage = $"Convergence not reached at load step {loadStep}";
					break;
				}

				// Increment displacements
				_currentDisplacements += CalculateDisplacements(GlobalStiffness, _currentResidual);
			}
		}

		/// <summary>
		///     Calculate residual force <see cref="Vector" />.
		/// </summary>
		/// <returns></returns>
		private Vector<double> ResidualForces() => _currentForces - InternalForces();

		/// <summary>
		///     Calculate convergence.
		/// </summary>
		private double Convergence()
		{
			double
				num = 0,
				den = 1;

			for (var i = 0; i < _currentResidual.Count; i++)
			{
				num += _currentResidual[i] * _currentResidual[i];
				den += _currentForces[i] * _currentForces[i];
			}

			return
				num / den;
		}

		/// <summary>
		///     Returns true if achieved convergence.
		/// </summary>
		/// <param name="convergence">Calculated convergence.
		///     <para>See: <see cref="Convergence" />.</para>
		/// </param>
		/// <param name="tolerance">Given tolerance for convergence.</param>
		/// <param name="iteration">Current iteration.</param>
		/// <param name="minIterations">Minimum number of iterations to perform (default: 5).</param>
		private bool ConvergenceReached(double convergence, double tolerance, int iteration, int minIterations = 5) => convergence <= tolerance && iteration >= minIterations;

		/// <summary>
		///     Returns true if achieved convergence.
		/// </summary>
		/// <param name="tolerance">Given tolerance for convergence.</param>
		/// <param name="iteration">Current iteration.</param>
		/// <param name="minIterations">Minimum number of iterations to perform (default: 5).</param>
		private bool ConvergenceReached(double tolerance, int iteration, int minIterations = 5) => ConvergenceReached(Convergence(), tolerance, iteration, minIterations);

		/// <summary>
		///     Save load step results after achieving convergence.
		/// </summary>
		/// <param name="loadFactor">Current load factor.</param>
		private void SaveLoadStepResults(double loadFactor)
		{
			if (!_monitoredIndex.HasValue)
				return;

			MonitoredDisplacements.Add(_currentDisplacements[_monitoredIndex.Value]);
			MonitoredLoadFactor.Add(loadFactor);
		}

		#endregion
	}
}